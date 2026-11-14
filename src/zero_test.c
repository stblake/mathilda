/*
 * zero_test.c — PossibleZeroQ: hybrid symbolic-numeric zero recognition.
 *
 * Pipeline (early exit at any stage that yields a definite verdict):
 *
 *   Stage 0.5 (normalisation) — when the input carries a constant-base
 *             exponential with a NON-LINEAR symbol-dependent exponent
 *             (Power[E, f] / Exp[f], f super-linear in a free symbol), ExpandAll
 *             distributes the sums so same-base exponentials
 *             become adjacent Times factors and collapse (E^a * E^(-a) -> 1).
 *             Value-preserving, so the verdict is unchanged; it spares the numeric
 *             ladder the tiny*huge catastrophic cancellation of a same-base
 *             exponential product split across summands (a Gaussian * Erf residual
 *             of an exact-ODE integrating-factor solution otherwise hangs the
 *             ladder).
 *
 *   Stage 0 — O(1) structural shortcuts: literal Integer/Real/BigInt/MPFR
 *             zero, Complex[0, 0], List of zeros, unbound symbol, …
 *
 *   Stage 1 — Rational normalisation via Together ∘ Cancel + Expand,
 *             then is_zero_poly. Decides every identity in Q(x_1,…,x_n).
 *
 *   Stage 2 — Closed-form numeric path: numericalize at machine precision,
 *             compare |z| against the IEEE catastrophic-cancellation
 *             ambiguity threshold scale * 2^(-p/2 + 4). If ambiguous, bump
 *             precision (machine → 200 → 500 → 1000 bits) and retry.
 *             A non-zero result stays roughly constant across precisions;
 *             a true zero shrinks geometrically. Surviving the full ladder
 *             implies "True".  Every magnitude (residual, scale, threshold)
 *             is a log2, so values far outside the double range -- Pi/3^700,
 *             10^500 Pi -- are measured, not rounded to 0 or Inf.
 *
 *   Stage 3 — Schwartz–Zippel. For inputs with free symbols, substitute
 *             each free symbol with a random REAL rational of moderate
 *             magnitude, recurse into Stage 2, and require independent
 *             confirmations. Sampling is real-line only: an analytic identity
 *             holding on a real interval holds on a complex neighbourhood
 *             (identity theorem), so real points confirm it, while complex
 *             samples needlessly cross branch cuts (Log/ArcTan/Sqrt) and blow
 *             up special functions (Gamma), manufacturing false negatives.
 *             The draw stream is seeded deterministically from the input's
 *             structural hash, so the verdict is a pure function of the input
 *             (no run-to-run flakiness) and the user's RNG stream is untouched.
 */
#include "zero_test.h"
#include "zero_test_distrib.h"

#include "arithmetic.h"
#include "attr.h"
#include "eval.h"
#include "expand.h"
#include "expr.h"
#include "internal.h"
#include "numeric.h"
#include "poly/poly.h"
#include "random.h"
#include "simp.h"
#include "simp_trigexp_zero.h"
#include "sym_names.h"
#include "symtab.h"
#include "sym_intern.h"   /* intern_symbol: PossibleZeroQ opaque-atom symbols */
#include <stdio.h>

#include <gmp.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* ------------------------------------------------------------------ */
/*  Tunables                                                          */
/* ------------------------------------------------------------------ */

/* Precision ladder, in MPFR bits. The first rung is treated as machine
 * precision (53 bits, IEEE 754). Each subsequent rung roughly squares
 * the precision so a single ambiguous result cannot survive long. */
static const long PRECISION_LADDER[] = { 53, 200, 500, 1000 };
#define PRECISION_LADDER_LEN ((int)(sizeof PRECISION_LADDER / sizeof PRECISION_LADDER[0]))

/* Stage 3 runs in two phases that exploit a sharp cost asymmetry:
 *
 *   - A sample that lands on a NON-zero branch is O(1) at machine precision
 *     and is rejected at ladder rung 0 (cheap).
 *   - Confirming a genuine identity requires climbing the full MPFR ladder
 *     per sample to rule out catastrophic-cancellation false zeros (costly).
 *
 * Phase A ("screen") draws many points but evaluates each only at machine
 * precision, catching branch-dependent non-zeros cheaply. A single 4-sample
 * budget was too small here: for an expression that is identically zero on a
 * real-analytic branch but genuinely non-zero on an adjacent branch (e.g.
 * D[2 Sqrt[1-Cos[x]],x] - Sqrt[1+Cos[x]], zero on (0,pi) but not (pi,2pi)),
 * each random real point lands in a zero interval ~1/2 the time, so all-k-zero
 * — a false positive — had probability ~(1/2)^4 = 6%. With SCREEN_SAMPLES
 * points the half-and-half case drops to ~(1/2)^24, and a lopsided 3:1
 * zero:non-zero split to ~(3/4)^24 ~ 1e-3. More POINTS, spread across many
 * periods of the moderate-magnitude real range, is what defeats it.
 *
 * Phase B ("confirm") is reached only when every screen point looked zero. It
 * climbs the full precision ladder on a few fresh points to reject
 * cancellation-hidden small non-zeros (the classic Schwartz–Zippel guarantee:
 * false-positive <= (degree/|S|)^k for k laddered samples). */
#define ZT_SCREEN_SAMPLES  24
/* Complex-domain screens (Assumptions -> Element[x, Complexes]): every sign
 * orthant of the symbols is visited at least this many times, half of the
 * visits with the complex-domain symbols on the real line, up to
 * ZT_CPLX_SCREEN_MAX points (see decide_schwartz_zippel_core). */
#define ZT_CPLX_ORTHANT_VISITS 16
#define ZT_CPLX_SCREEN_MAX     2048
#define ZT_CONFIRM_SAMPLES 4

/* Region-conforming rejection sampling (coupling assumptions). When an
 * assumption couples two free symbols (e.g. x > y) the per-symbol SampleSpec
 * cannot encode it, so each drawn assignment is checked against the coupling
 * facts and RE-DRAWN if it lands outside the assumed region. A half-measure
 * region such as {x > y} conforms on ~1 in 2 draws, so this bound is reached
 * only when the feasible region is tiny / empty (then the point is skipped and,
 * if every point is unreachable, the test honestly returns UNKNOWN). 128 draws
 * put the miss probability for a half-measure region at ~2^-128. */
#define ZT_CONFORM_MAX_TRIES 128

/* Sampling ranges for Stage 3. Samples are real rationals: numerator uniform
 * on [-2^NUMERATOR_BITS, 2^NUMERATOR_BITS], denominator on
 * [1, 2^DENOMINATOR_BITS].
 *
 * The numerator magnitude is kept MODERATE on purpose. A large real argument
 * overflows Exp/Gamma to Inf (silently degrading to UNKNOWN -> True) and
 * inflates special-function magnitudes so that an identically-zero difference
 * cancels far below the noise floor — exactly the catastrophic-cancellation
 * false negatives this sampler must avoid. +-64 still spans ~10 periods of
 * 2*pi, so periodic identities are still probed across many branches, and the
 * full-granularity denominator keeps the distinct-value set large (~2^23) for
 * the Schwartz-Zippel bound. Samples are purely real: an analytic identity
 * true on a real interval is true on a complex neighbourhood (identity
 * theorem), so real points suffice to confirm it, whereas complex points
 * cross branch cuts of Log/Sqrt/ArcTan where the symbolic identity legitimately
 * fails. */
#define ZT_NUMERATOR_BITS    6
#define ZT_DENOMINATOR_BITS  16

/* Magnitude-shell ladder for the Stage-3 overflow re-draw. When a sample point
 * numericalizes out of the numeric range (E/Gamma/... overflowing to +/-Inf at
 * the moderate |v| <= 2^ZT_NUMERATOR_BITS shell), the sampler retries the SAME
 * trial from successively smaller |v| upper bounds (2^k). This rescues points
 * that only overflow at the full range — a nowhere-zero function such as
 * E^(-10 t)+E^(-100 t) then decides FALSE at a finite point instead of the whole
 * test aborting to UNKNOWN -> True. Every shell keeps |v| >= 1 (the smallest,
 * 2^0, draws whole part == 1), preserving the floor that guards against
 * polynomial Schwartz-Zippel false positives (see sample_random_value). By the
 * identity theorem an analytic identity holding on [1, 2^k) holds everywhere, so
 * a smaller shell loses no decision power. The first rung is ZT_NUMERATOR_BITS so
 * a non-overflowing input draws exactly as before (byte-for-byte stream). */
static const int ZT_OVERFLOW_SHELL_BITS[] = { ZT_NUMERATOR_BITS, 4, 2, 1, 0 };
#define ZT_OVERFLOW_SHELL_LEN \
    ((int)(sizeof ZT_OVERFLOW_SHELL_BITS / sizeof ZT_OVERFLOW_SHELL_BITS[0]))

/* For a precision-honoured ALGEBRAIC expression whose numeric residual never
 * shrinks across the ladder, a residual exceeding scale * 2^(-ZT_ALG_NONZERO_BITS)
 * is above the machine-noise floor (~2^-52) and is a genuine non-zero, not a
 * cancellation-hidden zero.  Set a few bits inside 52 so a true algebraic zero
 * pinned at machine rounding noise (residual ~ scale * 2^-52) still gets the
 * lenient zero verdict, while a resolvable small value (e.g. Sqrt[10^12+1]-10^6
 * ~ scale * 2^-42) is correctly rejected.  (A true algebraic zero with deeper
 * machine cancellation would SHRINK under MPFR and never reach this branch.) */
#define ZT_ALG_NONZERO_BITS  48

/* Safety bits in the cancellation threshold. The IEEE-cancellation rule of
 * thumb is that a result is "ambiguous" if smaller than scale * 2^(-p/2);
 * we add a couple of extra bits as a guard against fma / parsing slop. */
#define ZT_AMBIGUITY_GUARD_BITS 4

/* Deep-zero early-exit for the precision ladder. Once a residual has been
 * OBSERVED to shrink geometrically AND has fallen below scale * 2^(-N), it is a
 * genuine zero to overwhelming confidence: a real non-zero cannot shrink below
 * its own magnitude, and no algebraic/transcendental identity this system
 * produces cancels past ~N bits. Climbing the remaining (500/1000-bit) rungs
 * on a large tree is then pure cost — the dominant expense on big parametric
 * antiderivative round-trips (a 60k-leaf D[r,x]-f numericalised at 1000 bits
 * costs ~0.8 s per sample). N is set far above machine precision (52) and any
 * realistic cancellation depth, yet far below the first MPFR rung (200 bits)
 * so a true zero triggers it at that rung. */
#define ZT_DEEP_ZERO_BITS 96

/* A symbol-free constant at a huge operand scale whose residual is still
 * falling at the top of the ladder is read once more, at enough precision to
 * resolve an absolute 2^-ZT_CONST_ABS_BITS against that scale (see
 * decide_numeric_huge_scale), but never above ZT_CONST_MAX_BITS. */
#define ZT_CONST_ABS_BITS 1000
#define ZT_CONST_MAX_BITS 65536

/* ------------------------------------------------------------------ */
/*  Forward declarations                                              */
/* ------------------------------------------------------------------ */

static ZeroTestResult decide_structural(const Expr* e);
static ZeroTestResult decide_rational(const Expr* e);
static ZeroTestResult decide_numeric(const Expr* e, bool* out_overflow);
static ZeroTestResult decide_schwartz_zippel(const Expr* e);
static ZeroTestResult decide_schwartz_zippel_assuming(const Expr* e, const AssumeCtx* ctx);
static bool           has_free_symbols(const Expr* e);
static bool           is_known_constant(const char* sym_name);

/* From src/simp/simp_util.c (declared in the module-internal simp_internal.h,
 * which pulls in the whole simplifier). Matches Rule/RuleDelayed[lhs, _] whose
 * LHS is the named symbol — reused here to parse the Assumptions option. */
extern bool is_rule_with_lhs(const Expr* e, const char* lhs_symbol);

/* ------------------------------------------------------------------ */
/*  Symbol-set helper (linear set of interned pointers)               */
/* ------------------------------------------------------------------ */

typedef struct {
    const char** items;   /* each element is an interned symbol pointer  */
    size_t       count;
    size_t       cap;
} SymPtrSet;

static void sps_init(SymPtrSet* s) { s->items = NULL; s->count = 0; s->cap = 0; }
static void sps_free(SymPtrSet* s) { free(s->items); s->items = NULL; s->count = 0; s->cap = 0; }

static bool sps_contains(const SymPtrSet* s, const char* sym) {
    for (size_t i = 0; i < s->count; ++i) if (s->items[i] == sym) return true;
    return false;
}

static void sps_add(SymPtrSet* s, const char* sym) {
    if (sps_contains(s, sym)) return;
    if (s->count == s->cap) {
        s->cap = s->cap ? s->cap * 2 : 4;
        s->items = realloc(s->items, sizeof(const char*) * s->cap);
    }
    s->items[s->count++] = sym;
}

/* ------------------------------------------------------------------ */
/*  Free-symbol detection                                             */
/* ------------------------------------------------------------------ */

/* True when `sym_name` is one of the named numeric constants the
 * evaluator/numericalize already knows how to fill in. These are NOT free
 * symbols even though they're EXPR_SYMBOL. */
static bool is_known_constant(const char* sym_name) {
    return sym_name == SYM_Pi          || sym_name == SYM_E          ||
           sym_name == SYM_I           || sym_name == SYM_EulerGamma ||
           sym_name == SYM_Catalan     || sym_name == SYM_GoldenRatio||
           sym_name == SYM_Glaisher    || sym_name == SYM_Khinchin   ||
           sym_name == SYM_GoldenAngle ||
           sym_name == SYM_Degree      || sym_name == SYM_Infinity   ||
           sym_name == SYM_ComplexInfinity || sym_name == SYM_Indeterminate ||
           sym_name == SYM_True        || sym_name == SYM_False;
}

/* True when `def->own_values` is non-empty — the symbol has been
 * assigned a value and is therefore not "free". */
static bool symbol_has_own_value(const char* sym_name) {
    SymbolDef* def = symtab_get_def(sym_name);
    return def && def->own_values != NULL;
}

/* Walk `e` and collect interned-pointer identities of free symbols into
 * `out`. Symbols appearing only at *head* positions (e.g. Sin in Sin[x])
 * are intentionally ignored: substituting a function name with a number
 * would produce nonsense like 3[x]. */
static void collect_free(const Expr* e, SymPtrSet* out) {
    if (!e) return;
    if (e->type == EXPR_SYMBOL) {
        const char* name = e->data.symbol.name;
        if (is_known_constant(name)) return;
        if (symbol_has_own_value(name)) return;
        sps_add(out, name);
        return;
    }
    if (e->type == EXPR_FUNCTION) {
        /* Skip head; recurse into args only. */
        for (size_t i = 0; i < e->data.function.arg_count; ++i) {
            collect_free(e->data.function.args[i], out);
        }
    }
}

/* True for a head that is NOT holomorphic in its argument: it is defined by
 * the real/imaginary split (Conjugate, Re, Im) or by the modulus/argument
 * (Abs, Arg, Sign).  For such a head the identity theorem does not apply — a
 * relation that holds on the real line (Conjugate[x] == x, Re[x] == x,
 * Abs[x]^2 == x^2) need not hold off it — so real sample points cannot
 * certify a zero. */
static bool is_nonholomorphic_head(const char* h) {
    return h == SYM_Conjugate || h == SYM_Re  || h == SYM_Im ||
           h == SYM_Abs       || h == SYM_Arg || h == SYM_Sign;
}

/* Collect the free symbols of `e` that occur INSIDE the argument of a
 * non-holomorphic head (see is_nonholomorphic_head).  Only those symbols need
 * complex sample points: a symbol that appears solely in holomorphic positions
 * leaves the expression holomorphic in it (for every fixed value of the
 * others), so real draws for it remain decisive by the identity theorem. */
static void collect_nonholomorphic_free(const Expr* e, bool inside, SymPtrSet* out) {
    if (!e) return;
    if (e->type == EXPR_SYMBOL) {
        if (inside) collect_free(e, out);
        return;
    }
    if (e->type != EXPR_FUNCTION) return;
    const Expr* h = e->data.function.head;
    bool in = inside || (h && h->type == EXPR_SYMBOL &&
                         is_nonholomorphic_head(h->data.symbol.name));
    for (size_t i = 0; i < e->data.function.arg_count; ++i)
        collect_nonholomorphic_free(e->data.function.args[i], in, out);
}

static bool has_free_symbols(const Expr* e) {
    SymPtrSet s; sps_init(&s);
    collect_free(e, &s);
    bool any = s.count > 0;
    sps_free(&s);
    return any;
}

/* True when `e` contains an algebraic-number constant that makes the
 * Stage-1 rational normalization (Together ∘ Cancel) potentially blow up:
 * a radical or root of unity, i.e. `Sqrt[...]` or `Power[base, p/q]` with
 * a non-integer rational exponent (q != 1).  Such subexpressions live in
 * an algebraic extension Q(α); combining them via Together/Cancel — which
 * `decide_rational` invokes WITHOUT an Extension option — is
 * super-polynomial in the extension degree (cyclotomic constants such as
 * (-1)^(2/3) are the pathological case).  When the expression also has free
 * symbols, the numeric Schwartz–Zippel stage decides identities by
 * sampling and never performs symbolic combination, so it is both fast and
 * sufficient (see zero_test_decide). */
static bool expr_has_algebraic_constant(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* head = e->data.function.head;
    if (head && head->type == EXPR_SYMBOL) {
        if (head->data.symbol.name == SYM_Sqrt && e->data.function.arg_count == 1)
            return true;
        if (head->data.symbol.name == SYM_Power && e->data.function.arg_count == 2) {
            int64_t p, q;
            if (is_rational(e->data.function.args[1], &p, &q) && q != 1)
                return true;
        }
        if (expr_has_algebraic_constant(head)) return true;
    }
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (expr_has_algebraic_constant(e->data.function.args[i])) return true;
    return false;
}

/* True if `e` contains a Power with a NON-CONSTANT (free-symbol) exponent, e.g.
 * x^(1-c).  decide_rational's Stage-1 Together ∘ Cancel does not terminate on
 * such symbolic-exponent powers — the pFq second solution x^(1-c) 2F1[...] of the
 * hypergeometric equation is the motivating case.  Like the algebraic-constant
 * guard, we then skip straight to numeric Schwartz–Zippel sampling, which
 * evaluates a symbolic exponent at a numeric sample without any symbolic
 * combination.  No decision power is lost: Stage 1's only trusted verdict is
 * TRUE, which the sampler reaches too (and it additionally decides FALSE). */
static bool expr_has_symbolic_exponent(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* head = e->data.function.head;
    size_t argc = e->data.function.arg_count;
    /* A VARIABLE base with a symbolic exponent — base and exponent both carry a
     * free symbol.  This is the x^(1-c) case that blows up Together ∘ Cancel.  A
     * constant base like E^(x^2/2) or 2^x is NOT matched (E/2 have no free
     * symbol): decide_rational treats those as opaque exp atoms and does not
     * choke, and routing every exp to the sampler would be both slow and, in the
     * variable-coefficient linear-system solver, verdict-changing. */
    if (head && head->type == EXPR_SYMBOL
        && head->data.symbol.name == SYM_Power && argc == 2
        && has_free_symbols(e->data.function.args[0])
        && has_free_symbols(e->data.function.args[1]))
        return true;
    if (expr_has_symbolic_exponent(head)) return true;
    for (size_t i = 0; i < argc; i++)
        if (expr_has_symbolic_exponent(e->data.function.args[i])) return true;
    return false;
}

/* True when the exponent `f` of an exponential depends on a free symbol
 * NON-LINEARLY: it contains a Power[b, e] with a free-symbol base and e != 1
 * (x^2, Sqrt[x]), a Power with a free-symbol exponent (c^x, x^x), a Times of two
 * or more free-symbol factors (x y), or any non-arithmetic function of a free
 * symbol (Sin[x], Log[x], a nested Exp).  A purely AFFINE exponent (c x + d,
 * I x) is benign and must NOT trip the exponential-combining normalisation: at
 * the sampler's moderate range (|x| <~ 64) an affine E^(c x) stays representable
 * (e^64 ~ 1e28, well inside double), so the precision ladder resolves it without
 * climbing, whereas a super-linear E^(x^2/2) reaches e^2048 and overflows —
 * that overflow, paired with a compensating huge (Erf) in another summand, is
 * the tiny*huge cancellation that defeats the ladder.  Narrowing to non-affine
 * exponents also leaves the E^x homogeneous terms of a constant-coefficient ODE
 * untouched, so ExpandAll cannot perturb a residual that is really a trig
 * identity into a sampler-hostile form. */
static bool exp_exponent_is_nonlinear(const Expr* f) {
    if (!f || f->type != EXPR_FUNCTION) return false;   /* symbol / number: affine */
    const Expr* head = f->data.function.head;
    size_t argc = f->data.function.arg_count;
    const char* nm = (head && head->type == EXPR_SYMBOL) ? head->data.symbol.name : NULL;
    if (nm == SYM_Plus) {
        for (size_t i = 0; i < argc; i++)
            if (exp_exponent_is_nonlinear(f->data.function.args[i])) return true;
        return false;
    }
    if (nm == SYM_Times) {
        int nfree = 0;
        for (size_t i = 0; i < argc; i++) {
            const Expr* a = f->data.function.args[i];
            if (exp_exponent_is_nonlinear(a)) return true;
            if (has_free_symbols(a)) nfree++;
        }
        return nfree >= 2;   /* x*const is affine; x*y (two free factors) is not */
    }
    if (nm == SYM_Power && argc == 2) {
        const Expr* b  = f->data.function.args[0];
        const Expr* ex = f->data.function.args[1];
        if (has_free_symbols(ex)) return true;                 /* c^x, x^x */
        if (has_free_symbols(b)
            && !(ex->type == EXPR_INTEGER && ex->data.integer == 1))
            return true;                                       /* x^2, Sqrt[x] */
        return false;
    }
    /* Any other head applied to a free-symbol argument (Sin, Log, nested Exp). */
    for (size_t i = 0; i < argc; i++)
        if (has_free_symbols(f->data.function.args[i])) return true;
    return false;
}

/* True when `e` contains a CONSTANT-base exponential with a NON-LINEAR
 * symbol-dependent exponent: Power[E, f] (base is the constant E) or Exp[f],
 * with exp_exponent_is_nonlinear(f).  This is the shape that produces the
 * tiny*huge magnitude swing which defeats the numeric ladder: a Gaussian E^(-g)
 * and the E^(+g) that appears from differentiating Erf[-I x/Sqrt2] land in
 * DIFFERENT summands, so the Times like-base collector never combines them to
 * E^0 = 1, and every Schwartz-Zippel sample numericalises as a catastrophic
 * cancellation.  Deliberately distinct from
 * expr_has_symbolic_exponent, which matches a VARIABLE base (x^(1-c)) and
 * EXCLUDES the constant-base case; and gated on nonlinearity so an affine E^x
 * (representable, ladder-friendly) does not needlessly trip ExpandAll. */
static bool expr_has_symbolic_exp_kernel(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* head = e->data.function.head;
    size_t argc = e->data.function.arg_count;
    if (head && head->type == EXPR_SYMBOL) {
        const char* nm = head->data.symbol.name;
        /* Power[E, f] with a non-linear free-symbol exponent. */
        if (nm == SYM_Power && argc == 2) {
            const Expr* base = e->data.function.args[0];
            if (base && base->type == EXPR_SYMBOL && base->data.symbol.name == SYM_E
                && exp_exponent_is_nonlinear(e->data.function.args[1]))
                return true;
        }
        /* Exp[f] with a non-linear free-symbol argument. */
        if (nm == SYM_Exp && argc == 1
            && exp_exponent_is_nonlinear(e->data.function.args[0]))
            return true;
    }
    if (expr_has_symbolic_exp_kernel(head)) return true;
    for (size_t i = 0; i < argc; i++)
        if (expr_has_symbolic_exp_kernel(e->data.function.args[i])) return true;
    return false;
}

/* Exponential-combining normalisation (Stage 0.5), applied by the public entry
 * points before any staged decision runs.  When `e` carries a symbol-dependent
 * exponential (expr_has_symbolic_exp_kernel), ExpandAll distributes the sums so
 * same-base exponentials become adjacent Times factors and collapse
 * (E^a * E^(-a) -> E^0 = 1), removing the tiny*huge catastrophic cancellation
 * BEFORE the numeric ladder (Stage 2/3) can climb to 1000 bits on it.  ExpandAll
 * is a pure algebraic identity, so the verdict is provably unchanged; the gate
 * only limits the cost to inputs that can actually exhibit the cancellation.
 * expr_expand_all rebuilds every node via eval_and_free, so the Times like-base
 * collector fires during expansion and the collapse needs no separate
 * evaluate().  Returns a fresh owned tree (caller frees), or NULL when no
 * exponential kernel is present (caller then uses `e` unchanged). */
static Expr* zt_normalize_exp_kernels(const Expr* e) {
    if (!expr_has_symbolic_exp_kernel(e)) return NULL;
    return expr_expand_all((Expr*)e);
}

/* True when `e` is a PURE RATIONAL FUNCTION of its free symbols over Q: every
 * node is an exact rational coefficient (Integer / BigInt / Rational[int,int]),
 * a free symbol (NOT a named constant — GoldenRatio etc. satisfy algebraic
 * relations that a polynomial-in-the-indeterminate view would miss), or a
 * Plus / Times / Power-with-INTEGER-exponent thereof.
 *
 * For such an expression the Stage-1 normalization (Together ∘ Cancel + Expand,
 * over Q) is EXACT and COMPLETE, so a non-zero normalized numerator is a
 * RIGOROUS non-zero — the exact realization of the DeMillo–Lipton–Schwartz–
 * Zippel guarantee, with no probability of a sampling false verdict.  This lets
 * decide_rational commit a trustworthy FALSE (not just its usual TRUE), and
 * removes such inputs from the numeric sampler entirely.
 *
 * Deliberately conservative: inexact coefficients (Real/MPFR — rigor is over Q),
 * Complex/Gaussian atoms (I is algebraic), non-integer or symbolic exponents
 * (radicals / transcendental powers), and every other head (Sin, Log, Sqrt,
 * user functions) make it return false, so the input keeps its previous path. */
static bool is_pure_rational_function(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return true;
    if (e->type == EXPR_REAL) return false;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return false;
#endif
    if (e->type == EXPR_STRING) return false;
    if (e->type == EXPR_SYMBOL) return !is_known_constant(e->data.symbol.name);
    if (e->type == EXPR_FUNCTION) {
        int64_t rn = 0, rd = 1;
        if (is_rational((Expr*)e, &rn, &rd)) return true;   /* Rational[int,int] */
        Expr* re = NULL; Expr* im = NULL;
        if (is_complex((Expr*)e, &re, &im)) return false;    /* Gaussian: I is algebraic */
        const Expr* h = e->data.function.head;
        if (!h || h->type != EXPR_SYMBOL) return false;
        const char* nm = h->data.symbol.name;
        size_t argc = e->data.function.arg_count;
        if (nm == SYM_Plus || nm == SYM_Times) {
            for (size_t i = 0; i < argc; ++i)
                if (!is_pure_rational_function(e->data.function.args[i])) return false;
            return true;
        }
        if (nm == SYM_Power && argc == 2) {
            /* Exponent must be an exact integer (positive or negative) — any
             * rational/symbolic exponent introduces a radical or transcendental. */
            if (!expr_is_integer_like(e->data.function.args[1])) return false;
            return is_pure_rational_function(e->data.function.args[0]);
        }
        return false;
    }
    return false;
}

/* ------------------------------------------------------------------ */
/*  Stage 0: structural shortcuts                                     */
/* ------------------------------------------------------------------ */

static bool integer_is_zero(const Expr* e) {
    if (e->type == EXPR_INTEGER) return e->data.integer == 0;
    if (e->type == EXPR_BIGINT)  return mpz_sgn(e->data.bigint) == 0;
    return false;
}

static bool real_is_zero(const Expr* e) {
    if (e->type == EXPR_REAL) return e->data.real == 0.0;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return mpfr_zero_p(e->data.mpfr) != 0;
#endif
    return false;
}

static ZeroTestResult decide_structural(const Expr* e) {
    if (!e) return ZERO_TEST_TRUE;

    if (integer_is_zero(e) || real_is_zero(e)) return ZERO_TEST_TRUE;

    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT ||
        e->type == EXPR_REAL
#ifdef USE_MPFR
        || e->type == EXPR_MPFR
#endif
        ) {
        /* Definitively non-zero numeric. */
        return ZERO_TEST_FALSE;
    }

    if (e->type == EXPR_STRING) return ZERO_TEST_FALSE;

    if (e->type == EXPR_SYMBOL) {
        const char* name = e->data.symbol.name;
        if (name == SYM_True || name == SYM_False) return ZERO_TEST_FALSE;
        if (name == SYM_Infinity || name == SYM_ComplexInfinity)
            return ZERO_TEST_FALSE;
        /* Pi, E, I, EulerGamma, GoldenRatio, Catalan, Degree, Glaisher,
         * Khinchin — known non-zero numeric constants. */
        if (is_known_constant(name)) return ZERO_TEST_FALSE;
        /* Any other symbol is a free variable — undetermined. */
        return ZERO_TEST_UNKNOWN;
    }

    if (e->type == EXPR_FUNCTION) {
        /* Complex[re, im] — both components must vanish. */
        Expr* re = NULL; Expr* im = NULL;
        if (is_complex((Expr*)e, &re, &im)) {
            ZeroTestResult rr = decide_structural(re);
            ZeroTestResult ri = decide_structural(im);
            if (rr == ZERO_TEST_TRUE && ri == ZERO_TEST_TRUE) return ZERO_TEST_TRUE;
            if (rr == ZERO_TEST_FALSE || ri == ZERO_TEST_FALSE) return ZERO_TEST_FALSE;
            return ZERO_TEST_UNKNOWN;
        }
        /* Rational[n, d] is structurally zero iff numerator is zero. */
        int64_t n = 0, d = 1;
        if (is_rational((Expr*)e, &n, &d)) return n == 0 ? ZERO_TEST_TRUE : ZERO_TEST_FALSE;
        /* ... and of any size: a canonical Rational with a BigInt part is never
         * zero.  Left to the numeric stage, 2^2000/3^1300 (~6e-19) read as 0,
         * its BigInt parts overflowing the double operand scale. */
        if (is_rational_like(e))
            return integer_is_zero(e->data.function.args[0]) ? ZERO_TEST_TRUE : ZERO_TEST_FALSE;
    }

    return ZERO_TEST_UNKNOWN;
}

/* ------------------------------------------------------------------ */
/*  Stage 0b: structural NON-zero certificate                          */
/* ------------------------------------------------------------------ */

/* A SOUND, unconditional proof that `e` is NOT identically zero: returns true
 * only when e is finite and non-vanishing for EVERY value of its free symbols,
 * so a caller may answer ZERO_TEST_FALSE with no numeric sampling and the
 * verdict holds under any assumption.  The lever is that the exponential is
 * entire and has no zeros — E^z != 0 for every (complex) z — so any product of
 * exponentials and finite non-zero constants is non-zero.  That is exactly the
 * class the Schwartz–Zippel sampler mishandles (c * E^f, e.g. -E^(-a x)): its
 * moderate |sample| range drives a symbol-dependent exponent into IEEE
 * overflow/underflow, so every point evaluates to 0 or Inf and the residual is
 * mistaken for an identity.
 *
 * Certified non-zero:
 *   - a finite non-zero numeric literal (Integer/BigInt/Real/MPFR; Rational
 *     with non-zero numerator; Complex with a finite non-zero component);
 *   - a known non-zero constant (E, Pi, EulerGamma, …) — never Infinity, which
 *     is excluded so it is safe as a power base below;
 *   - Exp[_] and Power[E, _] — the exponential has no zeros;
 *   - Power[b, _] with b itself certified non-zero (b finite, != 0  =>  b^e =
 *     e^{e log b} != 0 for any finite e);
 *   - Times[..] with every factor certified non-zero.
 * Anything else returns false — a bare symbol, Plus (cancellation), a power
 * with a possibly-zero or possibly-infinite base, any other head.  A false is
 * always safe: it merely forgoes the shortcut and lets the later stages decide. */
static bool provably_nonzero(const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: return e->data.integer != 0;
        case EXPR_BIGINT:  return mpz_sgn(e->data.bigint) != 0;
        case EXPR_REAL:    return e->data.real != 0.0 && isfinite(e->data.real);
#ifdef USE_MPFR
        case EXPR_MPFR:    return mpfr_number_p(e->data.mpfr) != 0 &&
                                  mpfr_zero_p(e->data.mpfr) == 0;
#endif
        case EXPR_SYMBOL:  return is_known_constant(e->data.symbol.name);
        case EXPR_FUNCTION: {
            /* Complex[re, im]: non-zero iff a component is a finite non-zero. */
            Expr* re = NULL; Expr* im = NULL;
            if (is_complex((Expr*)e, &re, &im))
                return provably_nonzero(re) || provably_nonzero(im);
            int64_t rn = 0, rd = 1;
            if (is_rational(e, &rn, &rd)) return rn != 0;
            const Expr* head = e->data.function.head;
            if (!head || head->type != EXPR_SYMBOL) return false;
            const char* hn = head->data.symbol.name;
            size_t argc = e->data.function.arg_count;
            if (hn == SYM_Exp && argc == 1) return true;
            if (hn == SYM_Power && argc == 2)
                return provably_nonzero(e->data.function.args[0]);
            if (hn == SYM_Times && argc > 0) {
                for (size_t i = 0; i < argc; ++i)
                    if (!provably_nonzero(e->data.function.args[i]))
                        return false;
                return true;
            }
            return false;
        }
        default: return false;
    }
}

/* ------------------------------------------------------------------ */
/*  Stage 1: rational normalization                                   */
/* ------------------------------------------------------------------ */

static ZeroTestResult decide_rational(const Expr* e) {
    /* Attempt 1: raw Expand.  Cheap for already-polynomial inputs. */
    Expr* expanded = expr_expand((Expr*)e);
    if (expanded) {
        bool z = is_zero_poly(expanded);
        expr_free(expanded);
        if (z) return ZERO_TEST_TRUE;
    }

    /* Attempt 2: Together ∘ Cancel.  Required for sums of rational forms
     * such as 1/x + 1/y - (x + y)/(x y), which Expand alone won't reduce. */
    Expr* tg = internal_together((Expr*[]){expr_copy((Expr*)e)}, 1);
    if (!tg) return ZERO_TEST_UNKNOWN;
    Expr* canon = internal_cancel((Expr*[]){tg}, 1);
    if (!canon) return ZERO_TEST_UNKNOWN;

    bool z = is_zero_poly(canon);
    if (!z) {
        /* Numerator-only zero-test: a/b == 0 iff a == 0, regardless of b.
         * Together/Cancel may have left a Times[num, Power[den, -1]] form
         * whose numerator is the real test target. */
        Expr* num = internal_numerator((Expr*[]){expr_copy(canon)}, 1);
        if (num) {
            Expr* ne = expr_expand(num);
            if (ne) {
                z = is_zero_poly(ne);
                expr_free(ne);
            }
            expr_free(num);
        }
    }
    expr_free(canon);
    if (z) {
        /* The Together ∘ Cancel normalization returned zero.  This is a RIGOROUS
         * zero only for a pure rational function of the free symbols, where the
         * Q-normalization is exact.  When the expression carries transcendental
         * kernels, Together applies branch-unsound identities on the way — it
         * expands Log[x^2] -> 2 Log[x] and Log[a b] -> Log[a] + Log[b] (correct
         * for the Risch log-tower machinery, but NOT valid across the branch
         * cuts of the principal-value Log) — which can collapse a genuine
         * NON-identity such as Log[x^2] - 2 Log[x] to zero.  Defer those to the
         * numeric sampler, which confirms true identities and detects the branch
         * dependence of the false ones. */
        if (is_pure_rational_function(e)) return ZERO_TEST_TRUE;
        return ZERO_TEST_UNKNOWN;
    }
    /* Normalization completed but the numerator is a non-zero polynomial.  If the
     * whole input is a pure rational function of its free symbols over Q, the
     * normalization was exact and complete, so this is a RIGOROUS non-zero
     * (no sampling, no probability of error).  Otherwise stay UNKNOWN: a
     * transcendental atom (Sin[x], Log[x], …) treated as an opaque indeterminate
     * can make a genuine identity look like a non-zero polynomial, so FALSE here
     * would be unsound — defer to the numeric sampler. */
    if (is_pure_rational_function(e)) return ZERO_TEST_FALSE;
    return ZERO_TEST_UNKNOWN;
}

/* ------------------------------------------------------------------ */
/*  Stage 2: numeric precision ladder                                 */
/* ------------------------------------------------------------------ */

/* ------------------------------------------------------------------ */
/*  Magnitudes                                                         */
/* ------------------------------------------------------------------ */

/* Every magnitude in the numeric stages -- a residual |z|, an operand scale,
 * a threshold -- is carried as its base-2 LOGARITHM, never as the value
 * itself.  A double holds |v| only in [2^-1074, 2^1024): read as a double, a
 * residual of 3.4e-334 (Pi/3^700) is 0 and "zero", and an operand scale of
 * 10^500 is Inf, which turns every relative test into Inf/NaN arithmetic.  The
 * numeric values themselves are not so limited -- machine numericalization
 * promotes an out-of-range result to an extended-exponent MPFR, and the MPFR
 * rungs have a 2^(+-2^62) exponent range -- so the magnitudes must not be
 * either.  As a log2 a double spans every magnitude MPFR can represent, a
 * product of magnitudes is a sum, a ratio a difference and "x below scale *
 * 2^-k" is "x < scale - k".
 *
 * Conventions: log2(0) = -INFINITY (an exact zero residual); +INFINITY stands
 * for a value no number represents (an IEEE Inf/NaN, a pole). */

/* log2(2^a + 2^b), exact zeros (-Inf) absorbed. */
static double zt_l2_add(double a, double b) {
    if (a == -INFINITY) return b;
    if (b == -INFINITY) return a;
    double hi = a > b ? a : b, lo = a > b ? b : a;
    if (!isfinite(hi)) return hi;
    return hi + log2(1.0 + exp2(lo - hi));
}

/* log2 sqrt(2^(2a) + 2^(2b)): the modulus of re + i im from log2|re|, log2|im|. */
static double zt_l2_hypot(double a, double b) {
    if (a == -INFINITY) return b;
    if (b == -INFINITY) return a;
    double hi = a > b ? a : b, lo = a > b ? b : a;
    if (!isfinite(hi)) return hi;
    return hi + 0.5 * log2(1.0 + exp2(2.0 * (lo - hi)));
}

/* log2 |z| of an MPFR value: exact over the whole exponent range. */
#ifdef USE_MPFR
static double zt_l2_mpfr(mpfr_srcptr v) {
    if (mpfr_zero_p(v)) return -INFINITY;
    if (!mpfr_number_p(v)) return INFINITY;
    long ex = 0;
    double m = mpfr_get_d_2exp(&ex, v, MPFR_RNDN);
    return (double)ex + log2(fabs(m));
}
#endif

/* log2 |z| of a GMP integer, whatever its size. */
static double zt_l2_mpz(const mpz_t z) {
    if (mpz_sgn(z) == 0) return -INFINITY;
    long ex = 0;
    double m = mpz_get_d_2exp(&ex, z);
    return (double)ex + log2(fabs(m));
}

/* log2 |v| for a numeric value: Integer, BigInt, Real, MPFR, Rational of
 * either kind of part (also the Rational[1., 1.*^30] shape numericalize
 * emits), and Complex[re, im] (its modulus).  Returns false if the input is
 * not a recognized numeric form (e.g. still contains symbols). */
static bool expr_log2_abs(const Expr* e, double* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) {
        *out = e->data.integer == 0 ? -INFINITY
                                    : log2(fabs((double)e->data.integer));
        return true;
    }
    if (e->type == EXPR_REAL) {
        double v = e->data.real;
        *out = v == 0.0 ? -INFINITY : (isfinite(v) ? log2(fabs(v)) : INFINITY);
        return true;
    }
    if (e->type == EXPR_BIGINT) { *out = zt_l2_mpz(e->data.bigint); return true; }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) { *out = zt_l2_mpfr(e->data.mpfr); return true; }
#endif
    Expr* re = NULL; Expr* im = NULL;
    if (is_complex((Expr*)e, &re, &im)) {
        double mr = 0.0, mi = 0.0;
        if (!expr_log2_abs(re, &mr)) return false;
        if (!expr_log2_abs(im, &mi)) return false;
        *out = zt_l2_hypot(mr, mi);
        return true;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2 &&
        e->data.function.head && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Rational) {
        double mn = 0.0, md = 0.0;
        if (!expr_log2_abs(e->data.function.args[0], &mn)) return false;
        if (!expr_log2_abs(e->data.function.args[1], &md)) return false;
        if (md == -INFINITY) return false;
        *out = (mn == -INFINITY) ? -INFINITY
             : (isfinite(mn) && isfinite(md)) ? mn - md : INFINITY;
        return true;
    }
    return false;
}

/* Signed real value of a numeric leaf (Integer / Real / BigInt / MPFR /
 * Rational[int, int]).  Returns false for anything that is not a real numeric
 * scalar (symbols, Complex, unevaluated functions).  Used to read a Power's
 * exponent with its sign so magnitude_scale handles reciprocals correctly. */
static bool expr_signed_double(const Expr* e, double* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = (double)e->data.integer;      return true; }
    if (e->type == EXPR_REAL)    { *out = e->data.real;                 return true; }
    if (e->type == EXPR_BIGINT)  { *out = mpz_get_d(e->data.bigint);    return true; }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR)    { *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true; }
#endif
    int64_t n = 0, d = 1;
    if (is_rational((Expr*)e, &n, &d) && d != 0) { *out = (double)n / (double)d; return true; }
    return false;
}

/* True iff every node of `e` is an algebraic-number construction: a numeric
 * leaf, Complex[…], Rational[…], Power[base, p/q] (rational exponent), Sqrt[…],
 * Plus[…] or Times[…] of the same.  For such an expression MPFR numericalize
 * honours the requested precision exactly (no special function is silently
 * pinned to machine precision), so a residual that does NOT shrink as bits are
 * added is a genuine non-zero rather than a cancellation-hidden zero.  Bare
 * symbols (transcendental constants Pi/E/…; free variables are substituted out
 * before the numeric stage) and all other heads return false. */
static bool is_algebraic_expr(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL || e->type == EXPR_BIGINT
#ifdef USE_MPFR
        || e->type == EXPR_MPFR
#endif
        ) return true;
    if (e->type == EXPR_SYMBOL || e->type == EXPR_STRING) return false;
    if (e->type == EXPR_FUNCTION) {
        Expr* re = NULL; Expr* im = NULL;
        if (is_complex((Expr*)e, &re, &im))
            return is_algebraic_expr(re) && is_algebraic_expr(im);
        int64_t rn = 0, rd = 1;
        if (is_rational((Expr*)e, &rn, &rd)) return true;
        const Expr* h = e->data.function.head;
        if (!h || h->type != EXPR_SYMBOL) return false;
        const char* nm = h->data.symbol.name;
        size_t argc = e->data.function.arg_count;
        if (nm == SYM_Power && argc == 2) {
            int64_t pn = 0, pd = 1;
            if (!is_rational(e->data.function.args[1], &pn, &pd)) return false;
            return is_algebraic_expr(e->data.function.args[0]);
        }
        if (nm == SYM_Sqrt && argc == 1)
            return is_algebraic_expr(e->data.function.args[0]);
        if (nm == SYM_Plus || nm == SYM_Times) {
            for (size_t i = 0; i < argc; ++i)
                if (!is_algebraic_expr(e->data.function.args[i])) return false;
            return true;
        }
        return false;
    }
    return false;
}

/* Estimate the operand-magnitude scale of `e`, as log2 (see "Magnitudes"
 * above). This is the denominator in the cancellation-aware ambiguity
 * threshold: results below scale * 2^(-p/2) are indistinguishable from
 * rounding noise, so we treat them as "ambiguous" rather than zero.
 *
 * The estimate uses head-aware recursion so it captures *additive*
 * cancellation (Plus) without spuriously inflating purely multiplicative
 * forms (Times, Rational, Power), where no cancellation can occur:
 *   - Plus[a, b, …]    -> Σ |a_i|       (potential cancellation)
 *   - Times[a, b, …]   -> Π |a_i|       (exact product magnitude)
 *   - Power[a, b]      -> |a|^b         (signed exponent)
 *   - Rational[n, d]   -> |n| / |d|
 *   - Complex[re, im]  -> hypot(|re|, |im|)
 *   - generic f[args]  -> Σ |a_i|       (conservative)
 * Numeric leaves give |value| (an exact zero leaf -Inf, so Complex[0, 10^-20]
 * scores 10^-20); named constants (Pi, E, …) get their numeric value at
 * machine precision; free symbols default to 1 (log2 0).  A compound node
 * whose estimate vanishes scores 1, as a scale of 0 would void every relative
 * test; magnitude_scale, the entry point, is therefore always finite. */
static double magnitude_scale_l2(const Expr* e) {
    if (!e) return 0.0;

    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL ||
        e->type == EXPR_BIGINT
#ifdef USE_MPFR
        || e->type == EXPR_MPFR
#endif
        ) {
        double m = 0.0;
        if (expr_log2_abs(e, &m) && m != INFINITY) return m;
        return 0.0;
    }

    if (e->type == EXPR_SYMBOL) {
        if (is_known_constant(e->data.symbol.name)) {
            Expr* n = numericalize(e, numeric_machine_spec());
            double m = 0.0;
            bool ok = n && expr_log2_abs(n, &m) && isfinite(m);
            if (n) expr_free(n);
            if (ok) return m;
        }
        return 0.0;
    }

    if (e->type == EXPR_FUNCTION) {
        Expr* re = NULL; Expr* im = NULL;
        if (is_complex((Expr*)e, &re, &im))
            return zt_l2_hypot(magnitude_scale_l2(re), magnitude_scale_l2(im));
        if (is_rational_like(e)) {
            double m = 0.0;
            if (expr_log2_abs(e, &m) && m != INFINITY) return m;
            return 0.0;
        }

        const Expr* head = e->data.function.head;
        const char* head_name = (head && head->type == EXPR_SYMBOL)
                                ? head->data.symbol.name : NULL;
        size_t argc = e->data.function.arg_count;

        if (head_name == SYM_Times) {
            double prod = 0.0;
            for (size_t i = 0; i < argc; ++i) prod += magnitude_scale_l2(e->data.function.args[i]);
            return isfinite(prod) ? prod : 0.0;
        }
        if (head_name == SYM_Power && argc == 2) {
            double base = magnitude_scale_l2(e->data.function.args[0]);
            if (base == -INFINITY) base = 0.0;
            /* Use the SIGNED exponent value, not its magnitude: Power[x, -1] has
             * magnitude 1/x, not x.  Scoring a reciprocal/denominator as x^|−1|
             * grossly inflates the operand scale of any expression with a
             * denominator, which then drowns a genuine non-zero below the
             * cancellation noise floor (a Schwartz-Zippel false zero).  A
             * symbolic exponent is scored by its own magnitude. */
            double expv;
            if (!expr_signed_double(e->data.function.args[1], &expv))
                expv = exp2(magnitude_scale_l2(e->data.function.args[1]));
            double v = base * expv;
            if (!isfinite(v)) return base; /* fall back */
            return v;
        }
        if (head_name == SYM_Rational && argc == 2) {
            double n = magnitude_scale_l2(e->data.function.args[0]);
            double d = magnitude_scale_l2(e->data.function.args[1]);
            if (d == -INFINITY) return n == -INFINITY ? 0.0 : n;
            double v = n - d;
            return isfinite(v) ? v : 0.0;
        }
        if (head_name == SYM_Sqrt && argc == 1) {
            double v = magnitude_scale_l2(e->data.function.args[0]);
            return v == -INFINITY ? 0.0 : 0.5 * v;
        }
        /* Exp[a], a real number: |E^a| = 2^(a log2 E), as Power[E, a] above
         * (the generic rule would score Exp[-10^4] as 10^4, not e^-10^4). */
        if (head_name == SYM_Exp && argc == 1) {
            double a = 0.0;
            if (expr_signed_double(e->data.function.args[0], &a) && isfinite(a))
                return a * 1.4426950408889634;   /* log2 E */
        }

        /* Plus: the L1 norm of the terms.  Generic function (Sin, Cos, Log,
         * Exp, user-defined, …): treat the result magnitude as approximately
         * the sum of operand magnitudes.  Adequate for the threshold
         * heuristic. */
        double sum = -INFINITY;
        for (size_t i = 0; i < argc; ++i)
            sum = zt_l2_add(sum, magnitude_scale_l2(e->data.function.args[i]));
        return isfinite(sum) ? sum : 0.0;
    }

    return 0.0;
}

static double magnitude_scale(const Expr* e) {
    double m = magnitude_scale_l2(e);
    return isfinite(m) ? m : 0.0;
}

/* "Definitely non-zero" threshold, as log2. At precision p, anything that
 * exceeds scale * 2^(-p + ZT_NONZERO_HEADROOM) is well above the rounding-noise
 * floor and is a confident False verdict. Tuned to forgive ~30 bits of
 * cancellation so transcendental identities that cancel down to machine
 * epsilon still survive the first rung. */
#define ZT_NONZERO_HEADROOM 30
static double nonzero_threshold(double scale, long p_bits) {
    long shift = p_bits - ZT_NONZERO_HEADROOM;
    if (shift <= 0) return INFINITY;          /* never reject at trivial precision */
    return scale - (double)shift;
}

/* Rung-0 (machine) obvious-non-zero gate. A residual exceeding scale * 2^-N is
 * too large to be rounding noise (machine eps ~ 2^-52) or plausible deeper
 * cancellation, so it settles FALSE without climbing the ladder (the fast path
 * for typical non-zeros). N is set well below 52 so that a true zero whose
 * cancellation runs PAST machine precision — residual up to ~2^-12 of scale,
 * i.e. ~2^40x amplification — is sent to the precision ladder for verification
 * rather than being misreported as non-zero here. This is the surgical fix for
 * the catastrophic-cancellation false negatives (e.g. Gamma identities). */
#define ZT_OBVIOUS_NONZERO_BITS 12

/* Operand scale above which a machine-precision "obvious non-zero" is NO LONGER
 * trusted directly.  When the operands are
 * astronomically large, a (1 - near-1) subtraction inside Tanh/Cosh/... loses
 * far more than the 52-bit mantissa, so the rung-0 residual of a genuine
 * cancellation ZERO can be a large fraction of scale (ratio well above the 2^-12
 * gate above) — indistinguishable at machine precision from a true non-zero.  A
 * residual at or above this scale that trips the obvious gate is therefore sent
 * to the full MPFR ladder (decide_numeric_huge_scale) rather than fast-rejected.
 * Set well below the real cancellation cases (Weierstrass antiderivative terms
 * are 2^50-2^73; the #5 reproducers 2^375+) so they are caught with margin,
 * while ordinary-scale computations (< 2^40 ~ 1e12) keep the byte-for-byte
 * unchanged machine fast path. */
#define ZT_HUGE_SCALE_BITS 40

/* Numericalize at the given precision, returning a freshly allocated
 * Expr*. spec.bits == 53 means machine; anything else uses MPFR (if
 * available). Returns NULL on failure or when MPFR is not compiled in
 * and an MPFR rung is requested. */
static Expr* numericalize_at(const Expr* e, long bits) {
    if (bits <= 53) {
        /* the machine rung screens: a doubtful point climbs the MPFR rungs,
         * so it runs in doubles, without N's working-precision raise */
        return numericalize_machine_unplanned(e);
    }
#ifdef USE_MPFR
    NumericSpec s; s.mode = NUMERIC_MODE_MPFR; s.bits = bits;
    return numericalize(e, s);
#else
    (void)e; (void)bits;
    return NULL;
#endif
}

/* True iff `e` is a definite infinity — the symbols Infinity / ComplexInfinity
 * or a DirectedInfinity[...]. A point that numericalizes to one of these is
 * definitively NON-zero (a pole, e.g. Gamma at a non-positive integer), so the
 * sampler must settle it FALSE rather than leaving it UNKNOWN (which would
 * collapse to True). A *cancellation* of infinities such as Gamma[x] - Gamma[x]
 * numericalizes to Indeterminate, not Infinity, so it is correctly NOT caught
 * here and stays UNKNOWN -> True. */
static bool expr_is_infinity(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL)
        return e->data.symbol.name == SYM_Infinity ||
               e->data.symbol.name == SYM_ComplexInfinity;
    if (e->type == EXPR_FUNCTION && e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_DirectedInfinity)
        return true;
    return false;
}

/* True iff `e` has finished numericalizing to a definite numeric value
 * (Integer, Real, BigInt, MPFR, Complex of same). False if any symbolic
 * residue remains. */
static bool is_pure_numeric(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL ||
        e->type == EXPR_BIGINT
#ifdef USE_MPFR
        || e->type == EXPR_MPFR
#endif
        ) return true;
    Expr* re = NULL; Expr* im = NULL;
    if (is_complex((Expr*)e, &re, &im)) {
        return is_pure_numeric(re) && is_pure_numeric(im);
    }
    int64_t rn = 0, rd = 1;
    if (is_rational((Expr*)e, &rn, &rd)) return true;
    /* Rational[any-numeric, any-numeric] — handles the inexact-numeric
     * shape numericalize emits for bigint-denominator rationals. */
    if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2 &&
        e->data.function.head && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Rational) {
        return is_pure_numeric(e->data.function.args[0]) &&
               is_pure_numeric(e->data.function.args[1]);
    }
    return false;
}

/* True iff `z` is a residue the evaluator declined to reduce to a number purely
 * for MAGNITUDE / range reasons — not because any free symbol or unknown head
 * survived. That is: every node is a pure-numeric leaf, or a function whose head
 * is a NUMERIC FUNCTION (ATTR_NUMERICFUNCTION) applied to such residues. The
 * live case is N[PolyLog[2, z]] for z beyond its numeric window (|z| >~ 300):
 * the dilogarithm engine DECLINES and hands back the unevaluated head
 * PolyLog[2.0, z] (Gamma, the original
 * example, no longer declines: it returns a finite extended-range value whose
 * magnitude overflows a double and is already re-drawn by the IEEE-overflow
 * path). Such a residue is a numeric coordinate the evaluator CAN compute at a
 * smaller argument, so the Stage-3 sampler should re-draw it from a smaller
 * magnitude shell (exactly like an IEEE overflow) rather than aborting the whole
 * test to UNKNOWN -> True on the first such point. An UNDEFINED / purely symbolic
 * head (e.g. UndefinedZQHead[x]) carries NO ATTR_NUMERICFUNCTION, so it is NOT a
 * magnitude decline: it stays a genuine residue and the test still aborts to
 * UNKNOWN (preserving the trip-wire that a truly undecidable input is assumed
 * zero). */
static bool residue_is_numeric_decline(const Expr* z) {
    if (!z) return false;
    if (is_pure_numeric(z)) return true;
    if (z->type == EXPR_FUNCTION) {
        const Expr* h = z->data.function.head;
        if (!h || h->type != EXPR_SYMBOL) return false;
        /* Overflow[] / Underflow[]: arithmetic left the MPFR exponent range --
         * the same out-of-range class as an IEEE overflow, re-drawn from a
         * smaller shell.  Aborting on it made PossibleZeroQ report True for the
         * wrong antiderivative of Sinh[x]/Sqrt[Sinh[x]^2 - 49/16] (its first
         * sample, x = 48.3, overflowed through Sinh[x]^4). */
        if ((h->data.symbol.name == SYM_Overflow || h->data.symbol.name == SYM_Underflow)
            && z->data.function.arg_count == 0)
            return true;
        if (!(get_attributes(h->data.symbol.name) & ATTR_NUMERICFUNCTION))
            return false;
        for (size_t i = 0; i < z->data.function.arg_count; ++i)
            if (!residue_is_numeric_decline(z->data.function.args[i]))
                return false;
        return true;
    }
    return false;
}

/* Build a NumericSpec for the given precision in bits. */
static NumericSpec spec_at_bits(long bits) {
    if (bits <= 53) return numeric_machine_spec();
#ifdef USE_MPFR
    NumericSpec s; s.mode = NUMERIC_MODE_MPFR; s.bits = bits;
    return s;
#else
    return numeric_machine_spec();
#endif
}

/* Refined operand-magnitude scale at a given precision. For Plus at the
 * top of `e`, the scale is the L1 norm of the numericalized arguments —
 * this correctly captures cancellation depth even when the inner
 * subexpressions (e.g. Sin[Complex[20, -17i]] with sinh blowup) carry
 * intermediate values orders of magnitude larger than any naive
 * tree-walk estimate would suggest.
 *
 * For non-Plus expressions we fall back to the static, head-aware
 * estimate which is generally tight enough (e.g. Sin[Pi] has scale Pi,
 * |result| 1.2e-16, threshold ~4e-7 → ambiguous → climb). */
static double magnitude_scale_at(const Expr* e, NumericSpec spec) {
    if (!e) return 0.0;
    if (e->type == EXPR_FUNCTION && e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Plus) {
        double sum = -INFINITY;
        for (size_t i = 0; i < e->data.function.arg_count; ++i) {
            Expr* n = spec.mode == NUMERIC_MODE_MACHINE
                    ? numericalize_machine_unplanned(e->data.function.args[i])
                    : numericalize(e->data.function.args[i], spec);
            if (n) {
                double m = 0.0;
                if (expr_log2_abs(n, &m) && isfinite(m)) sum = zt_l2_add(sum, m);
                expr_free(n);
            }
        }
        if (isfinite(sum)) return sum;
        /* Fall through to the static estimate. */
    }
    return magnitude_scale(e);
}

/* True iff a numeric value carries a machine real anywhere (Complex parts
 * included).  At an MPFR rung that is the mark of a head that ignored the
 * requested precision: machine reals contaminate every sum or product they
 * enter. */
static bool zt_has_machine_real(const Expr* z) {
    if (!z) return false;
    if (z->type == EXPR_REAL) return true;
    if (z->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < z->data.function.arg_count; ++i)
        if (zt_has_machine_real(z->data.function.args[i])) return true;
    return false;
}

/* Test the expression at a single precision rung. Returns FALSE if the
 * residual decisively exceeds the noise floor; TRUE if the residual is
 * indistinguishable from zero (very small relative to scale); UNKNOWN
 * if the numericalize couldn't reduce. *out_mag receives log2 of the
 * residual magnitude when non-UNKNOWN (-Inf for an exact zero, +Inf for a
 * pole), otherwise -Inf. *out_scale (when non-NULL) receives the log2
 * operand-magnitude scale used for the threshold.  *out_pinned (when
 * non-NULL) is set when an MPFR rung's value came back as a machine real:
 * some head evaluated at machine precision whatever was asked, so this rung
 * is no more accurate than the machine one.
 *
 * `in_scale`: when not NaN, the caller supplies a precomputed operand scale and
 * we SKIP the (expensive) per-argument numericalization in magnitude_scale_at.
 * The operand magnitude is essentially precision-independent (it is a sum/
 * product of operand values, not a cancelling residual), so a scale computed
 * once at machine precision is reused for every higher rung — this removes a
 * redundant full second numericalization pass at each MPFR rung, the dominant
 * cost on large (tens-of-thousands-of-leaves) antiderivative round-trips.
 * Pass NaN to compute the scale at this rung's precision. */
static ZeroTestResult evaluate_rung_core(const Expr* e, long bits,
                                         double* out_mag, double* out_scale,
                                         double in_scale, bool* out_overflow,
                                         bool* out_pinned) {
    *out_mag = -INFINITY;
    if (out_scale) *out_scale = 0.0;
    if (out_overflow) *out_overflow = false;
    if (out_pinned) *out_pinned = false;
    Expr* z = numericalize_at(e, bits);
    if (!z) return ZERO_TEST_UNKNOWN;
    if (expr_is_infinity(z)) {
        /* A pole: definitively non-zero. Report an infinite residual so every
         * caller's "obvious non-zero" gate settles it FALSE.  (A symbolic
         * Infinity/ComplexInfinity is left decisive on purpose — it is NOT the
         * IEEE-overflow class handled below.)
         * EXCEPT a machine-precision one: a denominator that cancels to exactly
         * 0.0 in doubles (a degree-30 polynomial in a Trager Log argument at
         * x = -60.7) also numericalizes to ComplexInfinity, and read as a pole it
         * made PossibleZeroQ reject a correct antiderivative.  A genuine pole is
         * still a pole at every precision; a cancellation artifact turns finite
         * once the precision outruns the cancellation (x = -31.35 needs the
         * 1000-bit rung), so each rung but the last defers to the next. */
        expr_free(z);
        for (int i = 0; i + 1 < PRECISION_LADDER_LEN; ++i)
            if (PRECISION_LADDER[i] >= bits || (i == 0 && bits <= 53))
                return evaluate_rung_core(e, PRECISION_LADDER[i + 1], out_mag, out_scale,
                                          in_scale, out_overflow, out_pinned);
        *out_mag = INFINITY;
        if (out_scale) *out_scale = 0.0;
        return ZERO_TEST_FALSE;
    }
    if (z->type == EXPR_SYMBOL && z->data.symbol.name == SYM_Indeterminate) {
        /* Indeterminate from numbers: 0 ComplexInfinity / 0/0 formed by
         * ROUNDING -- Coth[63.9] is exactly 1.0 in doubles, so
         * (-1 + Coth[x])^(3/2) in a denominator is 1/0 and the point reads
         * Indeterminate.  Like a machine ComplexInfinity above it may resolve at
         * a higher rung; if it never does, the POINT carries no information and
         * is flagged for re-draw from a smaller shell.  It used to abort the
         * whole test to UNKNOWN, which PossibleZeroQ reports as True: the
         * wrong-signed antiderivative of Cosh[x] Coth[x]^2 passed the
         * derivative-divides verifier on its first sample.  SOUND: an
         * expression Indeterminate at every point has no informative sample
         * and still ends UNKNOWN (screen_informative == 0). */
        expr_free(z);
        for (int i = 0; i + 1 < PRECISION_LADDER_LEN; ++i)
            if (PRECISION_LADDER[i] >= bits || (i == 0 && bits <= 53))
                return evaluate_rung_core(e, PRECISION_LADDER[i + 1], out_mag, out_scale,
                                          in_scale, out_overflow, out_pinned);
        if (out_overflow) *out_overflow = true;
        return ZERO_TEST_UNKNOWN;
    }
    if (!is_pure_numeric(z)) {
        /* A residue the evaluator declined to reduce. If it is a numeric
         * function evaluated on numeric arguments that merely ran out of its
         * numeric window (N[PolyLog[2, z]] -> PolyLog[2.0, z] for large |z|), it
         * is re-drawable just like an IEEE overflow: flag it so the Stage-3 shell
         * ladder shrinks the sample to a resolvable point.
         * A genuinely symbolic residue (an
         * undefined head) is NOT flagged and aborts the test to UNKNOWN as
         * before. */
        bool redraw = residue_is_numeric_decline(z);
        expr_free(z);
        if (redraw && out_overflow) *out_overflow = true;
        return ZERO_TEST_UNKNOWN;
    }
    double mag = 0.0;
    bool ok = expr_log2_abs(z, &mag);
    if (out_pinned && bits > 53) *out_pinned = zt_has_machine_real(z);
    expr_free(z);
    if (ok && mag == INFINITY) {
        /* Numeric but OUT OF RANGE: a symbol-dependent exponent drove E/Gamma/…
         * past the IEEE/MPFR range (+/-Inf), or two independently-overflowing
         * subterms produced a NaN. This POINT carries no information: it is a
         * numeric coordinate the double/MPFR range cannot represent, not a
         * symbolic residue (an unimplemented head is caught by is_pure_numeric
         * above) and not a genuine pole (expr_is_infinity above). Flag it so the
         * Stage-3 sampler re-draws from a smaller magnitude shell (or skips it)
         * rather than aborting the whole test to UNKNOWN — which would collapse
         * to a WRONG True for a nowhere-zero function such as E^(-10 t)+E^(-100 t).
         * SOUND: an overflow is never a zero, and a tiny*huge cancellation lands
         * here as NaN (or, when it cancels symbolically, as Indeterminate on the
         * !is_pure_numeric path) — never a decisive FALSE — so re-drawing can
         * neither invent a zero nor drop a genuine one. */
        if (out_overflow) *out_overflow = true;
        return ZERO_TEST_UNKNOWN;
    }
    if (!ok) return ZERO_TEST_UNKNOWN;
    *out_mag = mag;

    double scale = !isnan(in_scale) ? in_scale
                                    : magnitude_scale_at(e, spec_at_bits(bits));
    if (out_scale) *out_scale = scale;
    double tol = nonzero_threshold(scale, bits);
    if (mag > tol) return ZERO_TEST_FALSE;
    return ZERO_TEST_TRUE;
}

/* Public rung evaluator: brackets the per-sample numeric evaluation in the
 * arithmetic-warning mute.  Every sample PossibleZeroQ draws flows through here,
 * and a sample landing on a pole of the probed expression forms 1/0 -- a
 * Power::infy that is pure sampling noise: the pole is classified decisively
 * FALSE from the Infinity value in the core above, never from the message.
 * Muting leaves the verdict bit-identical while stopping the leak from reaching
 * the user or (worse) being noted by an enclosing Check[], which would otherwise
 * read a successful probe as a failure.  Single-exit wrapper so the push/pop pair
 * holds across the core's several early returns. */
static ZeroTestResult evaluate_rung(const Expr* e, long bits,
                                    double* out_mag, double* out_scale,
                                    double in_scale, bool* out_overflow,
                                    bool* out_pinned) {
    arith_warnings_mute_push();
    ZeroTestResult r = evaluate_rung_core(e, bits, out_mag, out_scale,
                                          in_scale, out_overflow, out_pinned);
    arith_warnings_mute_pop();
    return r;
}

/* Huge-scale deep-cancellation decision.
 * Reached only when the operand scale is huge (>= 2^ZT_HUGE_SCALE_BITS) AND the
 * rung-0 residual tripped the obvious-non-zero gate — a situation machine
 * precision CANNOT adjudicate, because a genuine cancellation zero and a genuine
 * small non-zero leave the SAME large residual until the artifact clears many
 * rungs up (empirically: Cosh[14]^30 - (1-Tanh[14]^2)^-15 and that + 1 are
 * bit-identical through 500 bits; they diverge only at the 1000-bit rung).  So
 * we do NOT early-exit on the first shrink (the deep-zero exit in decide_numeric
 * would misfire here): we climb the FULL ladder and let the residual's own trend
 * decide.
 *   - On the first PLATEAU (residual stops falling by >2x), the genuine value
 *     has resolved at this rung: trust the rung's floor verdict (`rr`).  A huge
 *     non-zero with no cancellation plateaus immediately at rung 1 and rejects
 *     after a single extra MPFR evaluation.
 *   - If the residual keeps shrinking to the TOP rung, trust that rung, whose
 *     floor (~scale * 2^-(1000-30)) is finally tight enough to tell a true small
 *     value (which plateaus just below it) from zero (which has fallen far past
 *     it).
 *   - If MPFR is unavailable, keep the lenient machine verdict (TRUE) — the same
 *     fallback decide_numeric uses, and never worse than HEAD's wrong FALSE. */
/* e with every finite machine real replaced by the exact dyadic rational it
 * stands for, or NULL when e holds none.  A sample point drawn as a machine
 * real pins every MPFR rung to machine precision -- a machine number stays
 * one under N -- so climbing the ladder on it could not tell a cancellation
 * artifact from a value: the residual plateaus at its rung-0 noise and reads
 * as non-zero. */
static Expr* zt_exactify_reals(const Expr* e) {
    if (!e) return NULL;
    if (e->type == EXPR_REAL) {
        if (!isfinite(e->data.real)) return NULL;
        mpq_t q;
        mpq_init(q);
        mpq_set_d(q, e->data.real);
        Expr* r = make_rational_mpz(mpq_numref(q), mpq_denref(q));
        mpq_clear(q);
        return r;
    }
    if (e->type != EXPR_FUNCTION) return NULL;
    size_t n = e->data.function.arg_count;
    Expr* h = zt_exactify_reals(e->data.function.head);
    Expr** args = NULL;
    for (size_t i = 0; i < n; i++) {
        Expr* a = zt_exactify_reals(e->data.function.args[i]);
        if (!a) continue;
        if (!args) {
            args = malloc(sizeof(Expr*) * n);
            for (size_t j = 0; j < n; j++) args[j] = NULL;
        }
        args[i] = a;
    }
    if (!h && !args) return NULL;
    Expr** out = malloc(sizeof(Expr*) * (n ? n : 1));
    for (size_t i = 0; i < n; i++)
        out[i] = (args && args[i]) ? args[i] : expr_copy(e->data.function.args[i]);
    Expr* r = expr_new_function(h ? h : expr_copy(e->data.function.head), out, n);
    free(out);
    free(args);
    return r;
}

static ZeroTestResult decide_numeric_huge_scale(const Expr* e, double mag0,
                                                double scale, bool constant) {
    /* Climb on the exact point: the MPFR rungs then carry real precision. */
    Expr* exact = zt_exactify_reals(e);
    if (exact) e = exact;
    double prev = mag0;
    ZeroTestResult rr = ZERO_TEST_TRUE;
    /* On an exact point every rung carries real precision, so climb them all:
     * a cancellation deeper than the low rungs (1 - Tanh[x/2]^2 cubed at |x| ~
     * 60 needs ~400 bits) leaves the first rungs at the same uncancelled value,
     * which is NOT a plateau of the true value.  Decide on the top two rungs: a
     * zero still falls at least half as fast as the bits added, a genuine value
     * holds. */
    if (exact) {
        double ms[PRECISION_LADDER_LEN];
        ZeroTestResult top = ZERO_TEST_TRUE;
        int got = 0;
        ms[0] = mag0;
        for (int i = 1; i < PRECISION_LADDER_LEN; ++i) {
            double m = 0.0;
            top = evaluate_rung(e, PRECISION_LADDER[i], &m, NULL,
                                isfinite(scale) ? scale : NAN, NULL, NULL);
            if (top == ZERO_TEST_UNKNOWN) { top = ZERO_TEST_TRUE; break; }  /* MPFR gone: lenient */
            ms[i] = m; got = i;
        }
        if (got >= 2) {
            long gained = PRECISION_LADDER[got] - PRECISION_LADDER[got - 1];
            if (ms[got] <= ms[got - 1] - (double)(gained / 2)) top = ZERO_TEST_TRUE;
        }
        expr_free(exact);
        return top;
    }
    bool plateaued = false;
    for (int i = 1; i < PRECISION_LADDER_LEN; ++i) {
        double m = 0.0;
        /* Reuse the machine-precision operand scale (precision-independent). */
        rr = evaluate_rung(e, PRECISION_LADDER[i], &m, NULL,
                           isfinite(scale) ? scale : NAN, NULL, NULL);
        if (rr == ZERO_TEST_UNKNOWN) { rr = ZERO_TEST_TRUE; plateaued = true; break; }  /* MPFR gone: lenient */
        if (m >= prev - 1.0) { plateaued = true; break; }  /* value resolved, trust rung */
        prev = m;
    }
    /* A CONSTANT whose residual is still falling at the top rung has not been
     * resolved: the top rung's floor sits at scale * 2^-970, and at a scale of
     * e^10^4 that is 2^13457 -- Exp[10^4] (Sin[Pi/7]^2 + Cos[Pi/7]^2 - 1) + 1,
     * which is 1, read as 0.  A sample point is re-drawn instead; a constant
     * cannot be, so it is read once more at the precision that resolves an
     * absolute 2^-ZT_CONST_ABS_BITS against its scale, and that rung's floor
     * decides (a zero has fallen below it, a value holds above it). */
    if (constant && !plateaued && rr == ZERO_TEST_TRUE && isfinite(scale)) {
        double want = scale + (double)ZT_CONST_ABS_BITS + ZT_NONZERO_HEADROOM;
        if (want > (double)PRECISION_LADDER[PRECISION_LADDER_LEN - 1] &&
            want <= (double)ZT_CONST_MAX_BITS) {
            double m = 0.0;
            bool pinned = false;
            ZeroTestResult top = evaluate_rung(e, (long)ceil(want), &m, NULL,
                                               scale, NULL, &pinned);
            if (top != ZERO_TEST_UNKNOWN && !pinned) rr = top;
        }
    }
    return rr;   /* the plateau rung's verdict, or the top-rung floor's */
}

/* Numeric Stage 2.  Strategy:
 *   1. Machine precision. A residual that is a non-trivial fraction of the
 *      operand scale (above the ZT_OBVIOUS_NONZERO_BITS gate) cannot be
 *      rounding noise and settles FALSE immediately — the fast path for typical
 *      non-zeros (EXCEPT at huge operand scale, where machine precision is
 *      untrustworthy — see decide_numeric_huge_scale). A smaller residual is
 *      AMBIGUOUS: it may be a true zero whose
 *      cancellation runs deeper than machine precision, so we climb rather than
 *      reject it here.
 *   2. Climb the precision ladder and decide on the SHRINKAGE TREND. A true
 *      zero shrinks geometrically as bits are added; a genuine non-zero
 *      plateaus at its true value. A FALSE verdict from a high rung is only
 *      trusted once we have OBSERVED the residual shrink — proof that the
 *      requested precision is honoured downstream. Special functions that
 *      silently stay at machine precision (residual constant across rungs) must
 *      NOT be rejected by a high rung's tiny threshold; they fall back to the
 *      lenient machine verdict. This is what stops the cancellation false
 *      negatives the previous absolute-threshold loop produced. */
/* A machine-precision residual above the obvious-non-zero gate is measured
 * against the TOP-LEVEL operand scale, but cancellation can run deeper than the
 * top level: in D[F] of 1/4 Log[1 + 2x^4 + 2x^3 y + 2x^2 Sqrt[1+x^4] + 2x y^3]
 * (y = (1+x^4)^(1/4)) the Log argument cancels to ~x^0 from terms of order x^4
 * at negative x, and the machine residual of an exact zero lands far above the
 * gate -- a "decisive" FALSE that made PossibleZeroQ reject correct Trager
 * antiderivatives.  A genuine non-zero keeps its value when the precision is
 * raised; cancellation noise does not.  So the machine verdict is resolved up
 * the ladder, on the EXACT point (zt_exactify_reals): a machine-real sample
 * stays machine under N, so climbing on it made every MPFR rung a re-run of the
 * machine one, and their agreement read as a stable non-zero.  A residual that
 * collapses below the deep-zero floor (or by 2^60 from its machine value) is a
 * cancellation-hidden zero (TRUE).  Otherwise the verdict is read off the TOP
 * two rungs, as decide_numeric_huge_scale does: a zero still falls by at least
 * half the bits added, a genuine value holds.  Agreement between LOWER rungs is
 * no evidence -- at x^2 = 205, 1/2 - Tanh[x^2]/2 ~ 1e-178 rounds to 0 at both
 * 200 and 500 bits, so Sech[x^2]^2 - 4 (1/2 - Tanh/2)(1/2 + Tanh/2) sat at
 * 1e-178 on both and an identity was reported non-zero (DSolve corpus
 * 2.2.16-1585).  An unavailable MPFR path keeps the machine verdict (FALSE). */
static ZeroTestResult zt_resolve_machine_nonzero(const Expr* e, double mag, double scale,
                                                 bool* out_overflow) {
    Expr* exact = zt_exactify_reals(e);
    if (exact) e = exact;
    double ms[PRECISION_LADDER_LEN];      /* residual at each decided MPFR rung */
    long   bits[PRECISION_LADDER_LEN];    /* ... and that rung's precision      */
    int got = 0;
    ZeroTestResult verdict = ZERO_TEST_FALSE;
    for (int i = 1; i < PRECISION_LADDER_LEN; ++i) {
        double m = 0.0;
        ZeroTestResult rr = evaluate_rung(e, PRECISION_LADDER[i], &m, NULL,
                                          isfinite(scale) ? scale : NAN, NULL, NULL);
        /* an undecidable rung (Indeterminate from a still-cancelling 0 * Infinity,
         * or MPFR unavailable) says nothing: try the next one */
        if (rr == ZERO_TEST_UNKNOWN) continue;
        /* an MPFR-rung pole is a genuine pole (evaluate_rung_core already defers
         * a machine-precision 0.0 denominator up the ladder): non-zero */
        if (m == INFINITY) { got = -1; break; }
        if (isfinite(scale) && m < scale - ZT_DEEP_ZERO_BITS) {
            verdict = ZERO_TEST_TRUE;
            got = -1;                     /* decided: skip the trend tests */
            break;
        }
        ms[got] = m; bits[got] = PRECISION_LADDER[i]; ++got;
    }
    if (got >= 1) {
        /* a residual that collapsed by 2^60 or more from its machine value over
         * the ladder was cancellation noise all along */
        if (ms[got - 1] < mag - 60.0) verdict = ZERO_TEST_TRUE;
        /* the top two decided rungs: a zero falls by >= half the bits added */
        else if (got >= 2) {
            long gained = bits[got - 1] - bits[got - 2];
            if (ms[got - 1] <= ms[got - 2] - (double)(gained / 2))
                verdict = ZERO_TEST_TRUE;
        }
        /* Still at or below the top rung's precision floor, measured against the
         * STRUCTURAL operand scale (magnitude_scale reaches inside every Plus,
         * where the top-level numeric scale only sees the already-cancelled
         * summands): the cancellation outruns the ladder, and the point cannot
         * tell a zero from a value.  At u ~ -355, 1/2 + Tanh[u]/2 ~ e^-710 ~
         * 1e-308 sits unchanged through 1000 bits.  Such a point is an
         * unrepresentable coordinate, exactly like an IEEE overflow: flag it for
         * the sampler's shell re-draw instead of voting FALSE.  Sound -- a
         * re-draw neither invents a zero nor drops a genuine one: a non-zero
         * function is resolvably non-zero at the smaller-shell points. */
        if (verdict == ZERO_TEST_FALSE && out_overflow) {
            double s_struct = magnitude_scale(e);
            long top = bits[got - 1];
            if (isfinite(s_struct) && top > ZT_NONZERO_HEADROOM &&
                ms[got - 1] < s_struct - (double)(top - ZT_NONZERO_HEADROOM)) {
                *out_overflow = true;
                verdict = ZERO_TEST_UNKNOWN;
            }
        }
    }
    if (exact) expr_free(exact);
    return verdict;
}

static ZeroTestResult decide_numeric(const Expr* e, bool* out_overflow) {
    double mag = -INFINITY, scale = 0.0;
    bool ovf = false;
    ZeroTestResult r = evaluate_rung(e, PRECISION_LADDER[0], &mag, &scale, NAN,
                                     &ovf, NULL);
    if (out_overflow) *out_overflow = ovf;
    if (r == ZERO_TEST_UNKNOWN) {
        /* A constant whose machine value left the representable range (Inf - Inf
         * = NaN in 10^500 (Pi + 10^-30) - 10^500 Pi when a term stays a double)
         * cannot be re-drawn like a sample point, and an out-of-range value is
         * never evidence of a zero: read it on the MPFR rungs instead, whose
         * exponent range holds it.  A sample point (out_overflow given) goes
         * back to the sampler's shell re-draw as before. */
        if (ovf && !out_overflow) {
            double m1 = -INFINITY, s1 = 0.0;
            ZeroTestResult r1 = evaluate_rung(e, PRECISION_LADDER[1], &m1, &s1,
                                              NAN, NULL, NULL);
            if (r1 == ZERO_TEST_UNKNOWN || !isfinite(m1)) return ZERO_TEST_UNKNOWN;
            return zt_resolve_machine_nonzero(e, m1, s1, NULL);
        }
        return ZERO_TEST_UNKNOWN;
    }

    /* Huge-scale deep-cancellation guard (#5): at an astronomically large
     * operand scale, a machine residual above the obvious gate may be a
     * catastrophic-cancellation artifact rather than a genuine value, so the
     * fast reject below is NOT sound here — climb the full ladder instead. */
    if (scale >= ZT_HUGE_SCALE_BITS && mag > scale - ZT_OBVIOUS_NONZERO_BITS)
        return decide_numeric_huge_scale(e, mag, scale, out_overflow == NULL);

    if (mag > scale - ZT_OBVIOUS_NONZERO_BITS)
        return zt_resolve_machine_nonzero(e, mag, scale, out_overflow);

    double prev_mag = mag;
    bool honored = false;
    for (int i = 1; i < PRECISION_LADDER_LEN; ++i) {
        double m = -INFINITY;
        bool pinned = false;
        /* Reuse the machine-precision operand scale (precision-independent) so
         * each higher rung numericalizes the tree ONCE, not twice. */
        /* Higher rungs never overflow when rung 0 did not (MPFR's exponent range
         * exceeds the double's), so no out_overflow is threaded here. */
        ZeroTestResult rr = evaluate_rung(e, PRECISION_LADDER[i], &m, NULL,
                                          scale, NULL, &pinned);
        if (rr == ZERO_TEST_UNKNOWN) {
            /* MPFR path unavailable beyond here — accept the lenient machine
             * verdict (the rung-0 residual was below the non-zero gate). */
            return ZERO_TEST_TRUE;
        }
        if (prev_mag == -INFINITY && m != -INFINITY && !pinned) {
            /* The rung below came out an EXACT zero and this one, evaluated at
             * its full precision (no machine-pinned head in the value), did not.
             * An exact zero at the lower rung is no baseline: the residual has
             * nothing to plateau against, and the lenient "never shrank" verdict
             * below would read any value as a zero -- Exp[1 + 10^-20] - E is 0.
             * at machine precision and 2.7*10^-20 at 200 bits.  The change itself
             * proves the precision is honoured, so this rung becomes the
             * baseline and the next rung decides: a cancellation zero keeps
             * falling, a value holds and is judged against that rung's floor. */
            honored = true;
            prev_mag = m;
            continue;
        }
        if (m < prev_mag - 1.0) {
            /* Residual shrank: precision honoured, still consistent with zero. */
            honored = true;
            prev_mag = m;
            /* Deep-zero early exit: the residual has shrunk geometrically AND is
             * now far below any plausible cancellation floor, so it is a genuine
             * zero.  Stop before the costly 500/1000-bit rungs (see
             * ZT_DEEP_ZERO_BITS).  This is the surgical fix for the "correct but
             * over budget" large-antiderivative round-trips (POSSIBLE_ZEROQ_
             * FAILURES.md case B2): it does not change any verdict — a genuine
             * non-zero cannot shrink below its own magnitude — only the number of
             * rungs climbed for a confirmed zero. */
            if (m < scale - ZT_DEEP_ZERO_BITS)
                return ZERO_TEST_TRUE;
            continue;
        }
        /* Residual plateaued at this rung. */
        if (honored) {
            /* Precision is honoured (earlier shrinkage proved it) yet the
             * residual stopped falling: its true non-zero value has emerged.
             * Trust the rung's noise-floor verdict (FALSE if above the floor,
             * TRUE if still within it). */
            return rr;
        }
        /* Never shrank. For a precision-honoured ALGEBRAIC expression (MPFR is
         * exact for radicals/arithmetic), a residual that stays well ABOVE the
         * machine-noise floor is a genuine non-zero whose operand scale is merely
         * inflated by cancellation (e.g. Sqrt[10^12+1] - 10^6 + z, or a canonic
         * cyclotomic eigenprojection) — settle it FALSE rather than mistaking it
         * for a cancellation-hidden zero.  Otherwise (transcendental heads whose
         * MPFR path may silently stay at machine precision, or a residual down at
         * the noise floor) keep the lenient machine verdict (zero). */
        if (is_algebraic_expr(e) && m > scale - ZT_ALG_NONZERO_BITS)
            return ZERO_TEST_FALSE;
        return ZERO_TEST_TRUE;
    }
    return ZERO_TEST_TRUE;
}

/* Cheap machine-precision screen verdict for Stage 3. FALSE only when the point
 * is OBVIOUSLY non-zero (residual a non-trivial fraction of scale, beyond the
 * cancellation band); TRUE when zero-ish at machine precision; UNKNOWN when it
 * cannot be reduced. A single machine-precision evaluation — the deep
 * cancellation check is deferred to the confirm phase's full ladder, so a
 * borderline-cancelling true-zero point is NOT falsely rejected by the screen. */
static ZeroTestResult screen_point(const Expr* e, bool* out_overflow) {
    double mag = -INFINITY, scale = 0.0;
    ZeroTestResult r = evaluate_rung(e, PRECISION_LADDER[0], &mag, &scale, NAN,
                                     out_overflow, NULL);
    if (r == ZERO_TEST_UNKNOWN) return ZERO_TEST_UNKNOWN;
    if (mag > scale - ZT_OBVIOUS_NONZERO_BITS) {
        /* Huge-scale deep-cancellation guard (#5): at an astronomically large
         * operand scale the machine FALSE may be a catastrophic-cancellation
         * artifact.  Don't settle the screen on it — return UNKNOWN so sz_trial
         * escalates this point to the full ladder (decide_numeric ->
         * decide_numeric_huge_scale), which climbs MPFR to adjudicate. */
        if (scale >= ZT_HUGE_SCALE_BITS) return ZERO_TEST_UNKNOWN;
        /* Deep cancellation below the top-level scale: confirm up the ladder. */
        return zt_resolve_machine_nonzero(e, mag, scale, out_overflow);
    }
    return ZERO_TEST_TRUE;
}

/* ------------------------------------------------------------------ */
/*  Stage 3: Schwartz–Zippel random substitution                      */
/* ------------------------------------------------------------------ */

/* Reuse the global PRNG owned by random.c via a thin shim: construct
 * a Mathilda-level RandomInteger[{lo, hi}] call and evaluate it. This
 * adds a few allocations per sample, which is fine for k=4 trials. */
static int64_t draw_int_range(int64_t lo, int64_t hi) {
    /* random_internal_int_range, not a RandomInteger[] call: these points are
     * part of the decision procedure, and routing them through the user-facing
     * builtin tied them to whichever generator it happens to use. When
     * RandomInteger moved to xoshiro the point set changed and the integrator
     * began claiming to solve the non-elementary Integrate[E^(Log[x]^2), x].
     * This also drops a handful of Expr allocations per sample. */
    return random_internal_int_range(lo, hi);
}

/* Build a random REAL sample value of moderate magnitude. A Real leaf
 * (EXPR_REAL) is used rather than an exact Rational[n, d] because several
 * Mathilda numeric heads only take the fast numeric path when the components
 * are already Real. Sampling is real-only by design (see the ZT_NUMERATOR_BITS
 * note): complex samples cross branch cuts and inflate special functions,
 * which manufactures cancellation/branch false negatives for genuine
 * real-line identities. */
/* Sign stratification for the Stage-3 screen.  When g_zt_sign_pattern >= 0,
 * the sign of symbol number g_zt_sym_index is bit g_zt_sym_index of the pattern
 * (screen trial t uses pattern t), so the first 2^n screen points visit every
 * orthant of n symbols.  An answer wrong only where a > 0, b > 0 and
 * x^2 < Min[a, b] lives in one orthant that independent random signs reach
 * rarely -- 96 random-sign screen points still missed it.  -1 = random signs. */
static int g_zt_sign_pattern = -1;
static int g_zt_sym_index = 0;
/* Real-axis schedule for complex-domain symbols in the screen: 1 = this
 * trial draws them real, 0 = complex, -1 = a random third real (the
 * confirm phase and every non-screen draw). */
static int g_zt_cplx_real = -1;
/* Magnitude band of a real visit of the complex-domain screen (-1: none).
 * The real draws are otherwise log-uniform on [1/8, 2^b], so a failure on a
 * short real segment -- (0, 1), where both radicands of a Kummer monomial are
 * negative and sit on their cuts -- was hit by ~1 real draw in 4, and an
 * orthant's 8 real visits missed it ~12% of the time.  Band k of 4 pins the
 * magnitude to [1/8, 1/2), [1/2, 1), [1, 2) or [2, 2^b]. */
static int g_zt_cplx_band = -1;

static bool zt_negative_draw(void) {
    if (g_zt_sign_pattern >= 0 && g_zt_sym_index < 30)
        return ((g_zt_sign_pattern >> g_zt_sym_index) & 1) != 0;
    return draw_int_range(0, 1) != 0;
}

/* A magnitude log-uniform on [2^-3, 2^num_bits] (see sample_random_value for
 * why: every scale, including |v| < 1, must be sampled). */
static double zt_log_uniform_magnitude(int num_bits) {
    double lo = -3.0, hi = (double)num_bits;
    if (g_zt_cplx_band >= 0) {
        static const double band_lo[4] = { -3.0, -1.0, 0.0, 1.0 };
        static const double band_hi[3] = { -1.0, 0.0, 1.0 };
        lo = band_lo[g_zt_cplx_band & 3];
        if ((g_zt_cplx_band & 3) < 3) hi = band_hi[g_zt_cplx_band & 3];
        if (hi < lo) hi = lo;
    }
    int64_t steps = (int64_t)1 << ZT_DENOMINATOR_BITS;
    int64_t k = draw_int_range(0, steps - 1);
    return exp2(lo + (hi - lo) * ((double)k + 0.5) / (double)steps);
}

static Expr* sample_random_value(int num_bits) {
    int64_t num_bound = (int64_t)1 << num_bits;
    int64_t den_bound = (int64_t)1 << ZT_DENOMINATOR_BITS;

    /* Sample magnitude is LOG-UNIFORM on [2^-3, 2^num_bits]: every scale gets
     * the same weight, both signs equally.
     *
     * It used to be bounded below by 1 (an integer part in [1, num_bound]),
     * on the identity theorem: "an analytic identity that holds on [1, 64]
     * holds everywhere".  That is false for principal-branch expressions,
     * which are only PIECEWISE analytic: -Sqrt[x^4 - 1]/(x^3 + x^5) is an
     * antiderivative of Sqrt[(x^2+1)/(x^2-1)] (3x^4 - 2x^2 - 3)/(x^4 (x^2+1)^2)
     * for |x| > 1 and of its NEGATIVE for |x| < 1, and PossibleZeroQ said True
     * (a false positive: 19 wrong ParallelMixedTower answers of the DDI corpus
     * passed).  A uniform draw also left a failing interval like [-3.3, -0.4]
     * unsampled about half the time.
     *
     * The floor 1/8 (not 0) keeps the guard the old floor stood for: a draw
     * near 0 drives a high-degree monomial (Sqrt[2] u^7) into the rounding band
     * and misreads a genuine non-zero as an identity.  Such a non-zero is still
     * caught by the large-magnitude draws, and any single non-zero point decides
     * FALSE.  The upper bound stays moderate (<= 2^num_bits) so Exp/Gamma do not
     * overflow. */
    double val = zt_log_uniform_magnitude(num_bits);
    (void)num_bound; (void)den_bound;
    if (zt_negative_draw()) val = -val;
    return expr_new_real(val);
}

/* ------------------------------------------------------------------ */
/*  Assumption-aware sampling                                          */
/*                                                                     */
/*  PossibleZeroQ[expr, Assumptions -> a] (and an ambient Assuming[] / */
/*  $Assumptions scope) restricts the Schwartz–Zippel draw to the      */
/*  assumed region.  An identity that holds only there — Sin[n Pi] for */
/*  integer n, Sqrt[x^2]-x for x>=0, Log[Exp[z]]-z on the principal    */
/*  strip — is genuinely non-zero at a generic real point, so the      */
/*  plain sampler correctly reports FALSE; conforming the samples is   */
/*  exactly what turns those into the intended TRUE.                   */
/*                                                                     */
/*  Correctness here is SOUNDNESS-only: a spec may sample any SUBSET   */
/*  of the assumed region (over-restriction is always safe); the one   */
/*  bug to avoid is drawing a point OUTSIDE it, which would misreport  */
/*  a genuine identity as non-zero.  No Simplify is ever invoked — by  */
/*  design PossibleZeroQ is a self-contained numeric/structural test.  */
/* ------------------------------------------------------------------ */

typedef enum { SDOM_REAL, SDOM_INT, SDOM_CPLX } SampleDomain;

/* One draw channel — a real value, or the real / imaginary part of a complex
 * value.  `sign` is -1 / 0 / +1; strict vs. non-strict is irrelevant because
 * the unconstrained magnitude draw keeps |v| >= 1/8 and so never lands on 0.
 * A finite two-sided [lo, hi] range takes precedence and disables the |v|>=1
 * shell (a range such as (-1, 1) lies entirely inside it, and theta near 0
 * must be sampleable). */
typedef struct {
    bool   has_lo, has_hi;
    double lo, hi;
    int    sign;
} Channel;

typedef struct {
    SampleDomain domain;
    Channel val;      /* value channel (REAL/INT) or real-part channel (CPLX) */
    Channel im;       /* imaginary-part channel (CPLX only)                    */
    bool mentioned;   /* some assumption fact constrains this symbol; when false
                       * the REAL domain is a default, not an assumption        */
} SampleSpec;

static void channel_init(Channel* c) {
    c->has_lo = false; c->has_hi = false; c->lo = 0.0; c->hi = 0.0; c->sign = 0;
}

/* Numericalize a (possibly symbolic-constant) bound such as Pi/2 or -1 to a
 * double.  Plain numeric leaves are read directly; anything else is run
 * through machine-precision numericalize first. */
static bool bound_to_double(const Expr* e, double* out) {
    if (expr_signed_double(e, out)) return true;
    Expr* n = numericalize(e, numeric_machine_spec());
    if (!n) return false;
    bool ok = expr_signed_double(n, out);
    expr_free(n);
    return ok;
}

/* Uniform double in [0, 1) drawn from the internal integer stream. */
static double draw_unit(void) {
    int64_t gran = (int64_t)1 << ZT_DENOMINATOR_BITS;
    return (double)draw_int_range(0, gran - 1) / (double)gran;
}

/* Draw one real value honouring `c`.  `integer` rounds to the integer grid.
 * `num_bits` caps the unconstrained magnitude shell at 2^num_bits (the overflow
 * re-draw shrinks it); it is IGNORED by the two-sided finite-range branch, whose
 * magnitude is pinned by the assumed interval and cannot shrink without leaving
 * the assumed region. */
static double sample_channel(const Channel* c, bool integer, int num_bits) {
    if (c->has_lo && c->has_hi) {
        /* Two-sided finite range: sample inside it (shell dropped). */
        double lo = c->lo, hi = c->hi;
        if (hi < lo) { double t = lo; lo = hi; hi = t; }
        if (integer) {
            int64_t ilo = (int64_t)ceil(lo);
            int64_t ihi = (int64_t)floor(hi);
            if (ihi < ilo) return (double)ilo;         /* empty grid: clamp */
            return (double)draw_int_range(ilo, ihi);
        }
        return lo + (hi - lo) * draw_unit();
    }
    /* Unconstrained magnitude: moderate, sign-honoured; integers on [1, 2^b],
     * reals log-uniform on [1/8, 2^b] (as sample_random_value). */
    int64_t num_bound = (int64_t)1 << num_bits;
    double v;
    if (integer) {
        v = (double)draw_int_range(1, num_bound);
    } else {
        v = zt_log_uniform_magnitude(num_bits);
    }
    if (c->sign < 0) v = -v;
    else if (c->sign == 0 && zt_negative_draw()) v = -v;
    return v;
}

/* Draw a conforming sample value for one free symbol.  spec == NULL restores
 * the legacy real-only draw byte-for-byte (see sample_random_value). */
static Expr* sample_random_value_spec(const SampleSpec* spec, int num_bits) {
    if (!spec) return sample_random_value(num_bits);
    if (spec->domain == SDOM_INT)
        return expr_new_integer((int64_t)sample_channel(&spec->val, true, num_bits));
    if (spec->domain == SDOM_CPLX) {
        bool real_visit = !spec->im.has_lo && !spec->im.has_hi && spec->im.sign == 0 &&
            (g_zt_cplx_real >= 0 ? g_zt_cplx_real == 1 : draw_int_range(0, 2) == 0);
        /* a real visit of the screen draws its magnitude from its band */
        int saved_band = g_zt_cplx_band;
        if (!real_visit) g_zt_cplx_band = -1;
        double re = sample_channel(&spec->val, false, num_bits);
        g_zt_cplx_band = -1;
        double im = sample_channel(&spec->im, false, num_bits);
        g_zt_cplx_band = saved_band;
        /* The real line is part of C, and a measure-zero part that random
         * complex draws never reach: an antiderivative that fails only on a
         * real-axis branch cut (Sqrt[(1-x)/(1+x)] on x < -1, where the integrand
         * is defined) must be caught too.  So a third of the draws are real,
         * unless Im is constrained (an explicit Im fact keeps its channel). */
        /* (a plain Real, not Complex[re, 0.]: a signed zero would put the point
         *  just beside the cut, where the identity still holds) */
        if (real_visit)
            return expr_new_real(re);
        Expr* parts[2] = { expr_new_real(re), expr_new_real(im) };
        return expr_new_function(expr_new_symbol(SYM_Complex), parts, 2);
    }
    return expr_new_real(sample_channel(&spec->val, false, num_bits));
}

/* Which channel a fact's variable operand refers to, for free symbol `sym`:
 * 0 = the bare symbol, 1 = Re[sym], 2 = Im[sym], -1 = not this symbol. */
static int channel_of_operand(const Expr* e, const char* sym) {
    if (!e) return -1;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == sym ? 0 : -1;
    if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 1 &&
        e->data.function.head && e->data.function.head->type == EXPR_SYMBOL) {
        const char* h = e->data.function.head->data.symbol.name;
        const Expr* a = e->data.function.args[0];
        if (a && a->type == EXPR_SYMBOL && a->data.symbol.name == sym) {
            if (h == SYM_Re) return 1;
            if (h == SYM_Im) return 2;
        }
    }
    return -1;
}

/* Intersect two channel constraints (b tightens a). */
static Channel merge_channel(Channel a, Channel b) {
    Channel r = a;
    if (b.has_lo && (!r.has_lo || b.lo > r.lo)) { r.has_lo = true; r.lo = b.lo; }
    if (b.has_hi && (!r.has_hi || b.hi < r.hi)) { r.has_hi = true; r.hi = b.hi; }
    if (r.sign == 0) r.sign = b.sign;
    return r;
}

/* Resolve a one-sided bound into either the signed magnitude shell (bound at
 * 0 — keeps |v| >= 1, avoids reintroducing the small-magnitude shell that
 * manufactures polynomial false positives) or an explicit range (nonzero
 * bound). A two-sided range is left as-is. */
static void finalize_channel(Channel* c) {
    double span = (double)((int64_t)2 << ZT_NUMERATOR_BITS);
    if (c->has_lo && c->has_hi) return;
    if (c->has_lo && !c->has_hi) {
        if (c->lo == 0.0) { c->sign = 1; c->has_lo = false; }
        else { c->has_hi = true; c->hi = c->lo + span; }
    } else if (c->has_hi && !c->has_lo) {
        if (c->hi == 0.0) { c->sign = -1; c->has_hi = false; }
        else { c->has_lo = true; c->lo = c->hi - span; }
    }
}

/* Build a SampleSpec for `sym` by scanning the assumption fact set.  Handles
 * Element[_, Domain], the binary relations (Greater/GreaterEqual/Less/
 * LessEqual/Equal), and Inequality[a, op, MID, op, b], where the constrained
 * operand may be the bare symbol, Re[sym], or Im[sym].  Facts it does not
 * recognise (Unequal, Or, ...) are ignored — a generic sample satisfies them
 * with probability 1. */
static SampleSpec extract_spec(const AssumeCtx* ctx, const char* sym) {
    Channel ch_bare, ch_re, ch_im;
    channel_init(&ch_bare); channel_init(&ch_re); channel_init(&ch_im);
    bool seen_int = false, seen_cplx = false;
    bool seen_re = false, seen_im = false, seen_im_zero = false;
    bool mentioned = false;

    for (size_t i = 0; ctx && i < ctx->count; ++i) {
        const Expr* f = ctx->facts[i];
        if (!f || f->type != EXPR_FUNCTION || !f->data.function.head ||
            f->data.function.head->type != EXPR_SYMBOL) continue;
        const char* h = f->data.function.head->data.symbol.name;
        size_t argc = f->data.function.arg_count;

        if (h == SYM_Element && argc == 2) {
            if (channel_of_operand(f->data.function.args[0], sym) != 0) continue;
            mentioned = true;
            const Expr* dom = f->data.function.args[1];
            if (!dom || dom->type != EXPR_SYMBOL) continue;
            const char* dn = dom->data.symbol.name;
            if      (dn == SYM_Integers)            seen_int = true;
            else if (strcmp(dn, "PositiveIntegers") == 0)    { seen_int = true; ch_bare.has_lo = true; ch_bare.lo = 0.0; }
            else if (strcmp(dn, "NonnegativeIntegers") == 0) { seen_int = true; ch_bare.has_lo = true; ch_bare.lo = 0.0; }
            else if (strcmp(dn, "NegativeIntegers") == 0)    { seen_int = true; ch_bare.has_hi = true; ch_bare.hi = 0.0; }
            else if (strcmp(dn, "NonpositiveIntegers") == 0) { seen_int = true; ch_bare.has_hi = true; ch_bare.hi = 0.0; }
            else if (dn == SYM_Complexes)           seen_cplx = true;
            /* Reals / Rationals / Algebraics confirm the REAL default. */
            continue;
        }

        if ((h == SYM_Greater || h == SYM_GreaterEqual ||
             h == SYM_Less    || h == SYM_LessEqual   || h == SYM_Equal) && argc == 2) {
            int ch; const Expr* cst; bool var_left;
            ch = channel_of_operand(f->data.function.args[0], sym);
            if (ch >= 0) { cst = f->data.function.args[1]; var_left = true; }
            else { ch = channel_of_operand(f->data.function.args[1], sym); cst = f->data.function.args[0]; var_left = false; }
            if (ch < 0) continue;
            mentioned = true;
            double c;
            if (!bound_to_double(cst, &c)) continue;
            if (ch == 1) seen_re = true;
            if (ch == 2) seen_im = true;
            Channel* tgt = (ch == 1) ? &ch_re : (ch == 2) ? &ch_im : &ch_bare;
            if (h == SYM_Equal) {
                if (ch == 2 && c == 0.0) seen_im_zero = true;
                tgt->has_lo = true; tgt->lo = c; tgt->has_hi = true; tgt->hi = c;
                continue;
            }
            bool ge = (h == SYM_Greater || h == SYM_GreaterEqual);
            bool lower = var_left ? ge : !ge;   /* does this lower-bound the var? */
            if (lower) { tgt->has_lo = true; tgt->lo = c; }
            else       { tgt->has_hi = true; tgt->hi = c; }
            continue;
        }

        if (h == SYM_Inequality && argc == 5) {
            int ch = channel_of_operand(f->data.function.args[2], sym);
            if (ch < 0) continue;
            mentioned = true;
            if (ch == 1) seen_re = true;
            if (ch == 2) seen_im = true;
            Channel* tgt = (ch == 1) ? &ch_re : (ch == 2) ? &ch_im : &ch_bare;
            double a, b;
            if (bound_to_double(f->data.function.args[0], &a)) { tgt->has_lo = true; tgt->lo = a; }
            if (bound_to_double(f->data.function.args[4], &b)) { tgt->has_hi = true; tgt->hi = b; }
            continue;
        }
        /* Unequal / Or / everything else: generic samples satisfy them. */
    }

    SampleSpec spec;
    spec.mentioned = mentioned;
    if (seen_int)                                    spec.domain = SDOM_INT;
    else if (seen_im_zero)                           spec.domain = SDOM_REAL;
    else if (seen_cplx || seen_re || seen_im)        spec.domain = SDOM_CPLX;
    else                                             spec.domain = SDOM_REAL;

    channel_init(&spec.val); channel_init(&spec.im);
    /* For REAL/INT (incl. Im[s]==0) the real value IS the real part, so fold
     * bare + Re[s] constraints together. For CPLX, Re[s] drives the real part
     * and Im[s] the imaginary part. */
    spec.val = merge_channel(ch_bare, ch_re);
    if (spec.domain == SDOM_CPLX) spec.im = ch_im;
    finalize_channel(&spec.val);
    if (spec.domain == SDOM_CPLX) finalize_channel(&spec.im);
    return spec;
}

/* Substitute every free symbol in `e` with an entry from `(syms, vals)`,
 * returning a fresh expression. Walks the tree directly so we don't pay
 * the pattern-matcher overhead of building rules. */
static Expr* substitute_symbols(const Expr* e, const char** syms, Expr** vals, size_t n) {
    if (!e) return NULL;

    if (e->type == EXPR_SYMBOL) {
        for (size_t i = 0; i < n; ++i) {
            if (e->data.symbol.name == syms[i]) return expr_copy(vals[i]);
        }
        return expr_copy((Expr*)e);
    }

    if (e->type == EXPR_FUNCTION) {
        Expr* new_head = substitute_symbols(e->data.function.head, syms, vals, n);
        size_t argc = e->data.function.arg_count;
        Expr** new_args = malloc(sizeof(Expr*) * (argc > 0 ? argc : 1));
        for (size_t i = 0; i < argc; ++i) {
            new_args[i] = substitute_symbols(e->data.function.args[i], syms, vals, n);
        }
        Expr* fn = expr_new_function(new_head, new_args, argc);
        free(new_args);
        return fn;
    }

    return expr_copy((Expr*)e);
}

/* ------------------------------------------------------------------ */
/*  Region-conforming rejection sampling for coupling assumptions      */
/* ------------------------------------------------------------------ */

/* True if `h` names a DENSE (positive-measure) relation. These carve a
 * half-space / interval, so a point can be rejection-sampled into the feasible
 * region in O(1) expected draws. Equal is deliberately absent: it pins a
 * MEASURE-ZERO variety the sampler can never hit, so a coupling equality keeps
 * the fact_keeps_false_sound FALSE->UNKNOWN downgrade instead. */
static bool is_dense_relation_head(const char* h) {
    return h == SYM_Greater || h == SYM_GreaterEqual ||
           h == SYM_Less    || h == SYM_LessEqual    ||
           h == SYM_Inequality || h == SYM_Unequal;
}

/* Evaluate a 2-operand comparison on concrete doubles. An unrecognised
 * operator (should not occur for the heads routed here) reports "does not
 * hold" — a conservative reject. */
static bool relation_true_binary(const char* h, double a, double b) {
    if (h == SYM_Greater)      return a >  b;
    if (h == SYM_GreaterEqual) return a >= b;
    if (h == SYM_Less)         return a <  b;
    if (h == SYM_LessEqual)    return a <= b;
    if (h == SYM_Equal)        return a == b;
    if (h == SYM_Unequal)      return a != b;
    return false;
}

/* Collect the variables appearing in a relation fact's OPERANDS (skipping the
 * operator symbols of an Inequality[a, op, b, op, c, ...], which sit at the
 * odd indices and would otherwise be miscounted as free variables). */
static void collect_relation_vars(const Expr* f, SymPtrSet* out) {
    if (!f || f->type != EXPR_FUNCTION || !f->data.function.head ||
        f->data.function.head->type != EXPR_SYMBOL) return;
    const char* h = f->data.function.head->data.symbol.name;
    size_t argc = f->data.function.arg_count;
    if (h == SYM_Inequality) {
        for (size_t i = 0; i < argc; i += 2)        /* operands at even indices */
            collect_free(f->data.function.args[i], out);
    } else {
        for (size_t i = 0; i < argc; ++i)
            collect_free(f->data.function.args[i], out);
    }
}

/* True iff the fully-substituted (concrete) relation `fs` holds numerically.
 * Operands are numericalised to doubles via bound_to_double; any operand that
 * fails to numericalise (e.g. a symbol the draw did not ground), or an
 * unrecognised shape, counts as NOT holding — a conservative reject (sampling
 * a few points fewer is always sound; sampling outside the region is the bug
 * this guards against). No Simplify, no evaluator re-entry. */
static bool concrete_fact_holds(const Expr* fs) {
    if (!fs || fs->type != EXPR_FUNCTION || !fs->data.function.head ||
        fs->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = fs->data.function.head->data.symbol.name;
    size_t argc = fs->data.function.arg_count;

    if (h == SYM_Inequality) {
        /* Inequality[x0, op1, x1, op2, x2, ...]: every (xi opi x{i+1}) must hold. */
        if (argc < 3 || (argc % 2) == 0) return false;
        double prev;
        if (!bound_to_double(fs->data.function.args[0], &prev)) return false;
        for (size_t i = 1; i + 1 < argc; i += 2) {
            const Expr* op = fs->data.function.args[i];
            if (!op || op->type != EXPR_SYMBOL) return false;
            double cur;
            if (!bound_to_double(fs->data.function.args[i + 1], &cur)) return false;
            if (!relation_true_binary(op->data.symbol.name, prev, cur)) return false;
            prev = cur;
        }
        return true;
    }

    if (argc == 2 && is_dense_relation_head(h)) {
        double a, b;
        if (!bound_to_double(fs->data.function.args[0], &a)) return false;
        if (!bound_to_double(fs->data.function.args[1], &b)) return false;
        return relation_true_binary(h, a, b);
    }
    return false;
}

/* True iff the drawn assignment (syms -> vals) satisfies EVERY coupling fact.
 * Each fact is grounded with the drawn values and checked numerically. */
static bool assignment_conforms(const Expr** facts, size_t nfacts,
                                const char** syms, Expr** vals, size_t nsyms) {
    for (size_t i = 0; i < nfacts; ++i) {
        Expr* fs = substitute_symbols(facts[i], syms, vals, nsyms);
        bool ok = concrete_fact_holds(fs);
        expr_free(fs);
        if (!ok) return false;
    }
    return true;
}

/* The verdict of one substituted point (screen or full ladder). */
static ZeroTestResult zt_point_verdict_raw(const Expr* sub, bool screen, bool* out_overflow);
#ifdef USE_MPFR
/* Is the sample point `sub` (a substituted, unevaluated tree) ON the branch cut
 * of one of its radicals or logarithms in a way no precision can resolve?  That
 * is a radicand / Log argument b that is mathematically real and negative but
 * reaches the cut only through complex arithmetic (TrigToExp forms: (I (E^(-I x)
 * - E^(I x)))/(E^(-I x) + E^(I x)) is Tan[x], real at real x), so its computed
 * imaginary part is rounding noise of either sign at every rung, and Sqrt or Log
 * of it lands on either side of the cut: the point's value is undetermined.
 * Read at 200 bits, such a b has Re b < 0 and 0 < |Im b| < 2^-150 |Re b| (a
 * true zero leaves ~2^-190 relative).  A b computed EXACTLY real (Im == 0, as
 * for (1 - x)/(1 + x) at a real x) is not ambiguous: the principal value is
 * defined, and a failure on that cut is a genuine one the sampler must see. */
static bool zt_on_ambiguous_cut(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    size_t n = e->data.function.arg_count;
    const Expr* b = NULL;
    if (h && h->type == EXPR_SYMBOL) {
        const char* hn = h->data.symbol.name;
        if (hn == SYM_Log && n == 1) b = e->data.function.args[0];
        else if (hn == SYM_Power && n == 2) {
            const Expr* p = e->data.function.args[1];
            if (!(p->type == EXPR_INTEGER || p->type == EXPR_BIGINT)) b = e->data.function.args[0];
        }
    }
    if (b && b->type == EXPR_FUNCTION) {
        Expr* v = numericalize_at(b, 200);
        Expr* re = NULL; Expr* im = NULL;
        bool amb = false;
        double r = 0.0, l2r = 0.0, l2i = 0.0;
        if (v && is_complex(v, &re, &im) && expr_signed_double(re, &r) && r < 0.0 &&
            expr_log2_abs(re, &l2r) && expr_log2_abs(im, &l2i) &&
            l2i != -INFINITY && l2i < l2r - 150.0)
            amb = true;
        if (v) expr_free(v);
        if (amb) return true;
    }
    for (size_t i = 0; i < n; i++)
        if (zt_on_ambiguous_cut(e->data.function.args[i])) return true;
    return false;
}
#endif

static ZeroTestResult zt_point_verdict(const Expr* sub, bool screen, bool* out_overflow) {
    ZeroTestResult r = zt_point_verdict_raw(sub, screen, out_overflow);
#ifdef USE_MPFR
    /* A non-zero reading at a point on an unresolvable cut (see
     * zt_on_ambiguous_cut) is no evidence: send the point back for a re-draw,
     * like an out-of-range coordinate.  Checked only on a FALSE, so the zero
     * path pays nothing. */
    if (r == ZERO_TEST_FALSE && out_overflow && zt_on_ambiguous_cut(sub)) {
        *out_overflow = true;
        return ZERO_TEST_UNKNOWN;
    }
#endif
    return r;
}

static ZeroTestResult zt_point_verdict_raw(const Expr* sub, bool screen, bool* out_overflow) {
    ZeroTestResult r = decide_structural(sub);
    if (r == ZERO_TEST_UNKNOWN) {
        if (screen) {
            r = screen_point(sub, out_overflow);
            /* a point flagged for re-draw (overflow / unresolvable cancellation)
             * goes back to the sampler; only an unreduced one climbs the ladder */
            if (r == ZERO_TEST_UNKNOWN && !(out_overflow && *out_overflow))
                r = decide_numeric(sub, out_overflow);
        } else {
            r = decide_numeric(sub, out_overflow);
        }
    }
    return r;
}

/* One Stage-3 trial: draw a fresh value for every free symbol, substitute it
 * into `e`, and return that point's verdict.
 *
 * Substitution is structural (no evaluation): the Mathilda evaluator would
 * eagerly numericalize, collapsing the top-level Plus into a single Complex
 * residue and discarding the operand magnitudes decide_numeric relies on for
 * its cancellation threshold. The structural shape (Plus / Times / Sin / ...)
 * is preserved through substitute_symbols and numericalized per rung.
 *
 *   - screen == true : machine precision only (screen_point). Cheap; rejects
 *     points that are OBVIOUSLY non-zero at this point. If machine precision
 *     cannot reduce the point at all (UNKNOWN), fall back to the full ladder so
 *     the point is still classified rather than silently passing the screen.
 *   - screen == false: climb the full precision ladder (decide_numeric) to
 *     reject cancellation-hidden small non-zeros. */
static ZeroTestResult sz_trial(const Expr* e, const char** syms,
                               const SampleSpec* specs, size_t nsyms,
                               bool screen, int num_bits, bool* out_overflow,
                               const Expr** conform_facts, size_t nconform,
                               bool* out_skip) {
    if (out_overflow) *out_overflow = false;
    if (out_skip) *out_skip = false;
    Expr** vals = malloc(sizeof(Expr*) * nsyms);

    /* Draw an assignment. When coupling facts are present (nconform > 0),
     * REJECT and re-draw any assignment outside the assumed region so only
     * conforming points are ever tested — this is what lets the per-symbol
     * sampler honour a cross-variable constraint such as x > y. When
     * nconform == 0 this draws exactly once, consuming the RNG stream
     * byte-for-byte as before (no behavioural change for non-coupling inputs). */
    int tries = 0;
    for (;;) {
        for (size_t i = 0; i < nsyms; ++i) {
            g_zt_sym_index = (int)i;
            vals[i] = sample_random_value_spec(specs ? &specs[i] : NULL, num_bits);
        }
        if (nconform == 0 ||
            assignment_conforms(conform_facts, nconform, syms, vals, nsyms))
            break;
        for (size_t i = 0; i < nsyms; ++i) expr_free(vals[i]);
        if (++tries >= ZT_CONFORM_MAX_TRIES) {
            free(vals);
            if (out_skip) *out_skip = true;   /* region unreachable: skip point */
            return ZERO_TEST_UNKNOWN;
        }
    }

    Expr* sub = substitute_symbols(e, syms, vals, nsyms);

    ZeroTestResult r = ZERO_TEST_UNKNOWN;
    if (sub) {
        r = zt_point_verdict(sub, screen, out_overflow);
        expr_free(sub);
    }

    for (size_t i = 0; i < nsyms; ++i) expr_free(vals[i]);
    free(vals);
    return r;
}

/* Salt mixed into the structural-hash seed so the sample distribution can be
 * re-tuned later (bump the salt) without colliding with any cached behaviour. */
#define ZT_SEED_SALT 0x5A3D9E1Bull

/* True if `op` denotes a single symbol's real value or its Re/Im channel: a
 * bare symbol, or Re[sym] / Im[sym] with sym a bare symbol. These are exactly
 * the operands extract_spec can fold into a per-symbol SampleSpec. */
static bool is_single_symbol_operand(const Expr* op) {
    if (!op) return false;
    if (op->type == EXPR_SYMBOL) return true;
    if (op->type == EXPR_FUNCTION && op->data.function.head &&
        op->data.function.head->type == EXPR_SYMBOL &&
        (op->data.function.head->data.symbol.name == SYM_Re ||
         op->data.function.head->data.symbol.name == SYM_Im) &&
        op->data.function.arg_count == 1 &&
        op->data.function.args[0] &&
        op->data.function.args[0]->type == EXPR_SYMBOL)
        return true;
    return false;
}

/* A ctx fact under which a sampler FALSE stays SOUND. The per-symbol sampler
 * (extract_spec) can only draw from single-variable regions, so it silently
 * ignores any fact that couples variables or constrains a compound. When such
 * a fact is present, an unconstrained draw may land at a point the assumptions
 * EXCLUDE, and a "found a non-zero point" FALSE is then unsound (this is the
 * a==b / a-b==0 bug: the sampler drew independent a,b, saw a-b!=0, and wrongly
 * returned FALSE). A fact keeps FALSE sound iff it is:
 *   - True / an opaque symbol            (no restriction),
 *   - Unequal[...]                       (excludes only a measure-zero set),
 *   - Element[single-symbol, Dom], or Element[_, Reals/Rationals/Algebraics/
 *     Complexes]                         (folded, or never excludes a draw),
 *   - a binary relation (>,>=,<,<=,==) with a single-symbol channel on one
 *     side and a numeric literal on the other  (folded into the channel), or
 *   - Inequality[num, op, single-symbol, op, num]  (folded).
 * Everything else — a coupling equality, an Equal/relation on a compound, a
 * discrete-domain membership of a compound, Or/Implies/Not/Xor — restricts the
 * space in a way the sampler cannot honour. TRUE is unaffected: an expression
 * that samples zero everywhere is zero on any sub-region too. */
static bool fact_keeps_false_sound(const Expr* f) {
    if (!f) return true;
    if (f->type != EXPR_FUNCTION || !f->data.function.head ||
        f->data.function.head->type != EXPR_SYMBOL) return true;
    const char* h = f->data.function.head->data.symbol.name;
    size_t argc = f->data.function.arg_count;

    if (h == SYM_Unequal) return true;

    if (h == SYM_Element && argc == 2) {
        const Expr* x = f->data.function.args[0];
        const Expr* d = f->data.function.args[1];
        if (is_single_symbol_operand(x)) return true;
        if (d && d->type == EXPR_SYMBOL) {
            const char* dn = d->data.symbol.name;
            if (dn == SYM_Reals || dn == SYM_Rationals ||
                dn == SYM_Algebraics || dn == SYM_Complexes)
                return true;   /* dense in the sampled space: never excludes */
        }
        return false;          /* Element[compound, discrete] restricts */
    }

    /* A single-symbol inequality (x > 0, a < s < b) is folded into the per-symbol
     * SampleSpec, so FALSE is sound. A COUPLING inequality (a > b) is NOT folded,
     * but it is not downgraded here either: decide_schwartz_zippel_core
     * rejection-samples it into its feasible region (conform_facts), so the draws
     * lie INSIDE the region and FALSE stays sound. (The earlier rationale — "an
     * inequality is full-measure, so identically-zero = zero-everywhere" — was
     * WRONG for a piecewise head like Max[x,y], which is a different function off
     * the region; the sampler must actually be confined to the region, which it
     * now is.) Only EQUALITY drops to a measure-zero variety the sampler cannot
     * hit, so a coupling equality is downgraded below. */
    if ((h == SYM_Greater || h == SYM_GreaterEqual || h == SYM_Less ||
         h == SYM_LessEqual) && argc == 2) return true;
    if (h == SYM_Inequality) return true;

    /* Equality restricts to a lower-dimensional variety where an expression can
     * be identically zero yet generically non-zero (a-b on a==b). Trust FALSE
     * only when extract_spec pins and samples the exact point: a single-symbol
     * channel == numeric literal. A coupling/compound equality is downgraded. */
    if (h == SYM_Equal && argc == 2) {
        const Expr* a = f->data.function.args[0];
        const Expr* b = f->data.function.args[1];
        double d;
        if (is_single_symbol_operand(a) && bound_to_double(b, &d)) return true;
        if (is_single_symbol_operand(b) && bound_to_double(a, &d)) return true;
        return false;          /* a==b, a-b==0, f[..]==g[..] restrict to a variety */
    }

    return false;              /* Or / Implies / Not / Xor / ... downgrade */
}

/* True if `e` contains any free symbol present in `set`. */
static bool expr_shares_symbol(const Expr* e, const SymPtrSet* set) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return sps_contains(set, e->data.symbol.name);
    if (e->type == EXPR_FUNCTION) {
        if (expr_shares_symbol(e->data.function.head, set)) return true;
        for (size_t i = 0; i < e->data.function.arg_count; ++i)
            if (expr_shares_symbol(e->data.function.args[i], set)) return true;
    }
    return false;
}

/* Run one Stage-3 trial through the overflow-shell ladder. If the drawn point
 * numericalizes OUT OF RANGE (a symbol-dependent exponent drove E/Gamma/… past
 * the IEEE/MPFR range to +/-Inf, or two independently-overflowing subterms gave a
 * NaN — evaluate_rung sets *overflow), re-draw the SAME trial from a successively
 * smaller |value| shell so the point falls back into range and becomes
 * informative. Returns the first DECISIVE (FALSE/TRUE) or non-overflow UNKNOWN
 * (genuine symbolic residue / Indeterminate — left to the caller to abort on,
 * exactly as before this fix) verdict; if EVERY shell was re-drawable, returns
 * UNKNOWN with *overflow == true so the caller can SKIP the uninformative point
 * instead of aborting the whole test. Re-drawn points are the IEEE-overflow class
 * AND a numeric-function-head magnitude decline (residue_is_numeric_decline, e.g.
 * N[PolyLog[2, z]] for large |z| — #4); a genuinely symbolic residue (undefined head) is NOT
 * re-drawn, so the sampler's behaviour (and its draw stream) on inputs that are
 * neither out-of-range nor declined is byte-for-byte unchanged. */
static ZeroTestResult sz_trial_shelled(const Expr* e, const char** syms,
                                       const SampleSpec* specs, size_t nsyms,
                                       bool screen, bool* overflow,
                                       const Expr** conform_facts, size_t nconform,
                                       bool* skip) {
    ZeroTestResult r = ZERO_TEST_UNKNOWN;
    for (int s = 0; s < ZT_OVERFLOW_SHELL_LEN; ++s) {
        *overflow = false;
        if (skip) *skip = false;
        r = sz_trial(e, syms, specs, nsyms, screen,
                     ZT_OVERFLOW_SHELL_BITS[s], overflow,
                     conform_facts, nconform, skip);
        if (skip && *skip) return r;   /* region unreachable: smaller shells can't help */
        if (!(r == ZERO_TEST_UNKNOWN && *overflow)) return r;
    }
    return r;   /* all shells overflowed: UNKNOWN, *overflow == true */
}

/* Shared engine for the plain and assumption-aware samplers. When `ctx` is
 * non-NULL a per-symbol SampleSpec restricts each draw to the assumed region;
 * ctx == NULL reproduces the legacy unconstrained real-only draw exactly. */
static ZeroTestResult decide_schwartz_zippel_core(const Expr* e, const AssumeCtx* ctx) {
    SymPtrSet syms; sps_init(&syms);
    collect_free(e, &syms);
    if (syms.count == 0) {
        sps_free(&syms);
        return ZERO_TEST_UNKNOWN;
    }

    /* Symbols under a non-holomorphic head (Conjugate/Re/Im/Abs/Arg/Sign) must
     * be drawn off the real line: Conjugate[x] - x vanishes at every real x but
     * is not zero.  A symbol whose REAL domain is only the default (no
     * assumption fact constrains it) is therefore sampled at a generic complex
     * point; an explicit Element[x, Reals] / x > 0 / ... keeps it real.  With
     * no such symbol the draw stream is byte-for-byte the legacy one. */
    SymPtrSet nonhol; sps_init(&nonhol);
    collect_nonholomorphic_free(e, false, &nonhol);

    SampleSpec* specs = NULL;
    if (ctx || nonhol.count) {
        specs = malloc(sizeof(SampleSpec) * syms.count);
        if (specs)
            for (size_t i = 0; i < syms.count; ++i) {
                /* extract_spec(NULL, ..) is the unconstrained REAL spec, whose
                 * draw equals sample_random_value's exactly. */
                specs[i] = extract_spec(ctx, syms.items[i]);
                if (!specs[i].mentioned && sps_contains(&nonhol, syms.items[i])) {
                    specs[i].domain = SDOM_CPLX;
                    channel_init(&specs[i].im);
                }
            }
    }
    sps_free(&nonhol);

    /* Classify every ctx fact touching e's symbols into one of three handlings,
     * because the per-symbol SampleSpec can only express single-variable regions:
     *
     *  (a) A fact that restricts to a MEASURE-ZERO set the sampler cannot hit
     *      (a==b, a-b==0, Element[a+b,Integers], Or[...]) — !fact_keeps_false_sound.
     *      A "found a non-zero point" FALSE is then unsound, so we downgrade
     *      FALSE -> UNKNOWN (false_untrusted). TRUE stays valid.
     *
     *  (b) A DENSE relation that COUPLES >= 2 of e's free symbols (x > y, the
     *      half-space the per-symbol spec drops) and is fully groundable by the
     *      draw — collected into conform_facts[]. The draw is rejection-sampled
     *      into {fact holds} (sz_trial), so the points actually lie INSIDE the
     *      assumed region and a FALSE is sound again. This is the fix for the
     *      Max[x,y]-x under x>y class, where the old "inequalities keep FALSE
     *      sound" shortcut was WRONG: a piecewise head is a different function
     *      off the region, so sampling the whole space mis-reported it non-zero.
     *
     *  (c) A coupling dense relation that references a symbol NOT in e (cannot be
     *      grounded by the draw, e.g. x > w with w absent from e): we cannot
     *      rejection-sample it, so fall back to the sound downgrade (a).
     *
     * A single-symbol-vs-literal inequality (x > 0, a < s < b) shares < 2 of e's
     * symbols, is already folded into the SampleSpec by extract_spec, and is left
     * untouched here — so conform_facts stays empty and the draw stream is
     * byte-for-byte unchanged for every non-coupling input. */
    bool false_untrusted = false;
    const Expr** conform_facts = NULL;
    size_t nconform = 0;
    if (ctx && ctx->count) {
        conform_facts = malloc(sizeof(Expr*) * ctx->count);
        for (size_t i = 0; i < ctx->count; ++i) {
            const Expr* f = ctx->facts[i];
            if (!expr_shares_symbol(f, &syms)) continue;
            if (!fact_keeps_false_sound(f)) { false_untrusted = true; continue; }  /* (a) */
            if (f->type == EXPR_FUNCTION && f->data.function.head &&
                f->data.function.head->type == EXPR_SYMBOL &&
                is_dense_relation_head(f->data.function.head->data.symbol.name)) {
                SymPtrSet fv; sps_init(&fv);
                collect_relation_vars(f, &fv);
                size_t shared = 0; bool grounded = true;
                for (size_t k = 0; k < fv.count; ++k) {
                    if (sps_contains(&syms, fv.items[k])) shared++;
                    else grounded = false;
                }
                sps_free(&fv);
                if (shared >= 2) {
                    if (grounded && conform_facts) conform_facts[nconform++] = f;  /* (b) */
                    else false_untrusted = true;                                   /* (c) */
                }
            }
        }
    }

    /* Seed the draw stream deterministically from the input's structural hash
     * so the verdict is a pure function of the input (no run-to-run flakiness).
     * The push/pop pair leaves the user's RandomInteger/SeedRandom stream
     * exactly as it was. */
    random_push_seed(expr_hash(e) ^ ZT_SEED_SALT);

    ZeroTestResult verdict = ZERO_TEST_TRUE;

    /* Phase A — screen: many cheap machine-precision points catch
     * branch-dependent non-zeros. A single decisively non-zero point settles
     * the whole test. A RE-DRAWABLE point (an out-of-range numeric coordinate, or
     * a numeric-function-head magnitude decline such as N[PolyLog[2, z]], |z| large) carries
     * no information: it is re-drawn from a smaller shell and, failing that,
     * SKIPPED — never allowed to abort the test (which would collapse to a wrong
     * True). Only a genuinely symbolic residue / Indeterminate UNKNOWN aborts. */
    int screen_informative = 0;
    /* More free symbols, more conditions a failing region can hide behind: an
     * antiderivative wrong only where a > 0, b > 0 and x^2 < Min[a, b] is hit by
     * ~2.5% of independent draws, and 24 points missed it about half the time.
     * The screen is machine precision only, so it scales with the symbol count
     * (24 per symbol, at most 96 points). */
    int screen_n = ZT_SCREEN_SAMPLES * (syms.count > 1 ? (int)syms.count : 1);
    if (screen_n > 4 * ZT_SCREEN_SAMPLES) screen_n = 4 * ZT_SCREEN_SAMPLES;
    /* Complex domain (the antiderivative verifiers ask for validity for all
     * complex x): a failure confined to a real-axis cut lives on the real slice
     * of ONE orthant of the symbols -- an answer of ParallelMixedTower for
     * (x^4 - a b)/((d x^4 - (a d + b d - c) x^2 + a b d) Sqrt[x^2 - a] Sqrt[x^2 - b])
     * is -f on parts of x^2 < a when a, b, c, d > 0: 2 orthants of 32, ~15% of
     * their real points, so ~1 real draw per orthant missed it 3 times in 4.
     * Each orthant is visited ZT_CPLX_ORTHANT_VISITS times, every other visit
     * with the complex-domain symbols real.  The screen is machine precision
     * (~0.5 ms a point), and only complex-domain tests pay for it. */
    int cplx_bits = -1;
    if (specs)
        for (size_t i = 0; i < syms.count; ++i)
            if (specs[i].domain == SDOM_CPLX) { cplx_bits = syms.count < 7 ? (int)syms.count : 7; break; }
    if (cplx_bits >= 0) {
        int n_cplx = ZT_CPLX_ORTHANT_VISITS << cplx_bits;
        if (n_cplx > ZT_CPLX_SCREEN_MAX) n_cplx = ZT_CPLX_SCREEN_MAX;
        if (n_cplx > screen_n) screen_n = n_cplx;
    }
    int saved_pattern = g_zt_sign_pattern;
    int saved_cplx_real = g_zt_cplx_real;
    for (int trial = 0; trial < screen_n; ++trial) {
        bool overflow = false, skip = false;
        g_zt_sign_pattern = trial;      /* orthant stratification (see zt_negative_draw) */
        if (cplx_bits >= 0) {
            g_zt_cplx_real = ((trial >> cplx_bits) & 1) == 0;
            g_zt_cplx_band = g_zt_cplx_real ? (trial >> (cplx_bits + 1)) & 3 : -1;
        }
        ZeroTestResult r = sz_trial_shelled(e, syms.items, specs, syms.count, true,
                                            &overflow, conform_facts, nconform, &skip);
        if (skip) continue;                               /* off-region unreachable: skip */
        if (r == ZERO_TEST_FALSE)   { verdict = false_untrusted ? ZERO_TEST_UNKNOWN : ZERO_TEST_FALSE; goto done; }
        if (r == ZERO_TEST_UNKNOWN) {
            if (overflow) continue;                       /* uninformative: skip */
            verdict = ZERO_TEST_UNKNOWN; goto done;       /* residue: cannot decide */
        }
        ++screen_informative;                             /* TRUE: informative zero-ish point */
    }
    g_zt_sign_pattern = saved_pattern;   /* confirm phase: random signs */
    g_zt_cplx_real = saved_cplx_real;
    g_zt_cplx_band = -1;
    /* Every screen point overflowed out of range — no informative evidence at
     * all — so the identity cannot be affirmed: return UNKNOWN honestly. */
    if (screen_informative == 0) { verdict = ZERO_TEST_UNKNOWN; goto done; }

    /* Phase B — confirm: every screened point looked zero. Climb the full
     * ladder on a few fresh points to reject cancellation-hidden small
     * non-zeros before declaring a genuine identity. Overflow points are again
     * re-drawn / skipped; if all confirm points overflow, the screen phase has
     * already established the finite points are zero-ish, so verdict stays TRUE
     * (mirrors the MPFR-unavailable fallback in decide_numeric). */
    for (int trial = 0; trial < ZT_CONFIRM_SAMPLES; ++trial) {
        bool overflow = false, skip = false;
        ZeroTestResult r = sz_trial_shelled(e, syms.items, specs, syms.count, false,
                                            &overflow, conform_facts, nconform, &skip);
        if (skip) continue;                               /* off-region unreachable: skip */
        if (r == ZERO_TEST_FALSE)   { verdict = ZERO_TEST_FALSE;   goto done; }
        if (r == ZERO_TEST_UNKNOWN) {
            if (overflow) continue;                       /* uninformative: skip */
            verdict = ZERO_TEST_UNKNOWN; goto done;       /* residue: cannot decide */
        }
    }

done:
    g_zt_sign_pattern = saved_pattern;
    g_zt_cplx_real = saved_cplx_real;
    g_zt_cplx_band = -1;
    random_pop_seed();
    free(specs);
    free((void*)conform_facts);
    sps_free(&syms);
    return verdict;
}

static ZeroTestResult decide_schwartz_zippel(const Expr* e) {
    return decide_schwartz_zippel_core(e, NULL);
}

static ZeroTestResult decide_schwartz_zippel_assuming(const Expr* e, const AssumeCtx* ctx) {
    return decide_schwartz_zippel_core(e, ctx);
}

/* ------------------------------------------------------------------ */
/*  Public entry points                                               */
/* ------------------------------------------------------------------ */

/* The staged decision procedure proper.  Operates on an expression already
 * put through zt_normalize_exp_kernels by the public entry point below. */
static ZeroTestResult zt_decide_core(const Expr* e) {
    ZeroTestResult r;

    r = decide_structural(e);
    if (r != ZERO_TEST_UNKNOWN) return r;

    /* Stage 0a: a generalized function (DiracDelta / its derivatives /
     * HeavisideTheta).  Pointwise sampling cannot decide one -- a delta vanishes
     * at every sample, and Theta does not evaluate at a complex one -- so it is
     * decided as a distribution: sifted delta coefficients and the regular piece
     * on each interval between breakpoints, each by this zero test
     * (zero_test_distrib.c).  An unsupported shape falls through. */
    if (zt_has_distribution(e)) {
        r = zt_decide_distribution(e, zero_test_decide);
        if (r != ZERO_TEST_UNKNOWN) return r;
    }

    /* Stage 0b: a structural non-zero certificate (c * E^f and the like).  Sound
     * and unconditional, and it must run BEFORE the algebraic/symbolic-exponent
     * guards below, which route straight to Schwartz–Zippel — the sampler cannot
     * see these, since a symbol-dependent exponent overflows/underflows E at its
     * moderate range so every point reads as 0/Inf (e.g. -E^(-a x)). */
    if (provably_nonzero(e)) return ZERO_TEST_FALSE;

    /* Phase 2 (SIMPLIFY_IMPROVEMENT_PLAN): when the expression mixes free
     * symbols with an algebraic-number constant (radical / root of unity),
     * skip the Stage-1 Together ∘ Cancel — over an extension Q(α) it blows
     * up super-polynomially (cyclotomic constants are the worst case).  Go
     * straight to numeric Schwartz–Zippel sampling, which decides true
     * identities without any symbolic combination.  This loses no decision
     * power: Stage 1's only trustworthy verdict is TRUE, and the sampler
     * reaches the same TRUE for genuine identities (and FALSE for genuine
     * non-zeros), so the public PossibleZeroQ result is unchanged — only
     * the non-terminating symbolic path is avoided. */
    if (has_free_symbols(e) && expr_has_algebraic_constant(e))
        return decide_schwartz_zippel(e);
    /* Symbolic-exponent guard (mirrors the algebraic-constant one): a x^(1-c)
     * power sends Stage-1 Together ∘ Cancel non-terminating; sample instead. */
    if (has_free_symbols(e) && expr_has_symbolic_exponent(e))
        return decide_schwartz_zippel(e);

    r = decide_rational(e);
    /* Trust TRUE always; trust FALSE only for the rigorous pure-rational path
     * (decide_rational returns FALSE exclusively when is_pure_rational_function
     * holds, where the Q-normalization is exact and complete). */
    if (r != ZERO_TEST_UNKNOWN) return r;

    /* Exact trig/exp kernel zero-test.  A rational function of a single kernel
     * t = E^(i x) (Sin[k x], Cos[k x], Sec[k x], …, with opaque Log/Sqrt terms
     * and parameters as independent generators) is decided by EXACT rational
     * point-evaluation on a Nullstellensatz grid — a symbolic decision, no
     * numeric sampling.  decide_rational leaves these UNKNOWN because it treats
     * every transcendental as an opaque indeterminate; this stage closes the
     * Sec^n/Csc^n and symbolic-parameter diff-back identities exactly and fast.
     * Declines (UNKNOWN → fall through) on non-single-kernel forms. */
    {
        TrigExpZeroResult tz = trigexp_rational_is_zero(e);
        if (tz == TRIGEXP_ZERO_TRUE)  return ZERO_TEST_TRUE;
        if (tz == TRIGEXP_ZERO_FALSE) return ZERO_TEST_FALSE;
    }

    if (!has_free_symbols(e)) {
        /* Closed-form: no re-draw is possible, so overflow is not tracked here.
         * A symbol-free overflow (e.g. E^10000) is a genuine huge value that the
         * pole path / lenient ladder already handle. */
        r = decide_numeric(e, NULL);
        if (r != ZERO_TEST_UNKNOWN) return r;
        return ZERO_TEST_UNKNOWN;
    }

    return decide_schwartz_zippel(e);
}

/* Public entry.  Apply the exponential-combining normalisation (Stage 0.5),
 * then run the staged decision procedure on the normalised expression.  The
 * normalisation is value-preserving, so the verdict is unchanged; it only spares
 * the numeric ladder the tiny*huge cancellation of a same-base exponential
 * product split across summands. */
ZeroTestResult zero_test_decide(const Expr* e) {
    Expr* norm = zt_normalize_exp_kernels(e);
    ZeroTestResult r = zt_decide_core(norm ? norm : e);
    if (norm) expr_free(norm);
    return r;
}

/* The assumption-aware staged procedure proper (ctx is guaranteed non-empty and
 * the input already normalised by the public entry below). */
static ZeroTestResult zt_decide_assuming_core(const Expr* e, const struct AssumeCtx* ctx) {
    /* A structurally decided literal (0, a non-zero constant, Complex[0,0], …)
     * is unconditional, so its verdict holds under any assumption. */
    ZeroTestResult r = decide_structural(e);
    if (r != ZERO_TEST_UNKNOWN) return r;

    /* Stage 0b: unconditional non-zero certificate.  E^z has no zeros for any z,
     * so c * E^f is non-zero for every value the symbols may take — the FALSE is
     * sound even with assumptions in scope, and it must precede the guards that
     * route to constrained sampling (which cannot see these). */
    if (provably_nonzero(e)) return ZERO_TEST_FALSE;

    /* Algebraic-constant guard (mirrors zero_test_decide): a Sqrt / rational
     * power over free symbols must skip Together ∘ Cancel and go straight to
     * constrained sampling. */
    if (has_free_symbols(e) && expr_has_algebraic_constant(e))
        return decide_schwartz_zippel_assuming(e, ctx);
    if (has_free_symbols(e) && expr_has_symbolic_exponent(e))
        return decide_schwartz_zippel_assuming(e, ctx);

    /* Trust ONLY the TRUE verdict from the unconditional exact stages. An
     * unconditional FALSE can be wrong under an assumption — trigexp proves
     * Exp[2 Pi I k] - 1 nonzero as a function of a continuous k, yet it is
     * identically zero for integer k. So the exact stages may fast-path a
     * genuine identity to TRUE, but constrained sampling is the SOLE arbiter
     * of FALSE (and UNKNOWN) once assumptions are in play. */
    if (decide_rational(e) == ZERO_TEST_TRUE) return ZERO_TEST_TRUE;
    if (trigexp_rational_is_zero(e) == TRIGEXP_ZERO_TRUE) return ZERO_TEST_TRUE;

    if (has_free_symbols(e))
        return decide_schwartz_zippel_assuming(e, ctx);

    /* Closed-form constant: assumptions are irrelevant to a definite number.
     * No re-draw is possible for a symbol-free input, so overflow is not tracked. */
    return decide_numeric(e, NULL);
}

/* Public entry (assumption-aware).  Same exponential-combining normalisation as
 * zero_test_decide, then the assumption-aware staged procedure. */
ZeroTestResult zero_test_decide_assuming(const Expr* e, const struct AssumeCtx* ctx) {
    /* No usable assumptions → the legacy path, byte-for-byte. */
    if (!ctx || ctx->count == 0) return zero_test_decide(e);
    Expr* norm = zt_normalize_exp_kernels(e);
    ZeroTestResult r = zt_decide_assuming_core(norm ? norm : e, ctx);
    if (norm) expr_free(norm);
    return r;
}

/* Read the $Assumptions OwnValue WITHOUT evaluating it (mirrors the reader in
 * simp_builtins.c): evaluating would make Element recurse through the very
 * OwnValue we are reading. Returns an owned copy, or True if unset. */
static Expr* pzq_read_dollar_assumptions(void) {
    Rule* r = symtab_get_own_values("$Assumptions");
    if (!r || !r->replacement) return expr_new_symbol(SYM_True);
    return expr_copy(r->replacement);
}

/* And-combine two owned assumption expressions and evaluate (canonicalises
 * And[True, x] -> x, flattens nested And). Both inputs are consumed. */
static Expr* pzq_combine_with_and(Expr* a, Expr* b) {
    Expr* args[2] = { a, b };
    Expr* call = expr_new_function(expr_new_symbol(SYM_And), args, 2);
    Expr* out = evaluate(call);
    expr_free(call);
    return out;
}

/* True when the effective assumption carries no usable facts (literal True or
 * an empty context). */
static bool assumption_is_trivial(const Expr* a) {
    return !a || (a->type == EXPR_SYMBOL && a->data.symbol.name == SYM_True);
}

/* Collapse a ZeroTestResult to the PossibleZeroQ public boolean: only a
 * proven/strongly-believed non-zero is False; UNKNOWN collapses to True, the
 * documented "assume zero when uncertain" behaviour. */
static Expr* pzq_bool(ZeroTestResult r) {
    return expr_new_symbol(r == ZERO_TEST_FALSE ? SYM_False : SYM_True);
}

/* ------------------------------------------------------------------ *
 *  Opaque atoms for the PossibleZeroQ builtin.                        *
 *                                                                     *
 *  The sampler substitutes numbers for symbols only, so an inert call  *
 *  -- f[a], C[1], a package's Pair[Momentum[p], Momentum[q]] -- stays  *
 *  symbolic at every point, the test ends UNKNOWN, and the builtin     *
 *  reported True: PossibleZeroQ[f[a]] was True and FeynCalc accepted a *
 *  momentum shift that changes a propagator.  The reference treats     *
 *  such a call as an independent quantity, so the builtin abstracts    *
 *  each distinct maximal one to a fresh symbol (equal calls share it,  *
 *  so f[a] - f[a] stays True) and samples that.  A call is opaque when *
 *  its head is not a symbol (Derivative[1][f][x]) or names neither a   *
 *  builtin nor a NumericFunction.  Only the builtin does this:         *
 *  zero_test_decide's internal callers read UNKNOWN as UNKNOWN.        *
 * ------------------------------------------------------------------ */

typedef struct { Expr** from; Expr** to; size_t n, cap; } PzqAtoms;

static bool pzq_opaque_call(const Expr* e) {
    if (e->type != EXPR_FUNCTION) return false;
    const Expr* h = e->data.function.head;
    if (h->type != EXPR_SYMBOL) return true;
    const char* nm = h->data.symbol.name;
    if (nm == SYM_Complex || nm == SYM_Rational) return false;
    SymbolDef* d = symtab_lookup(nm);
    if (!d) return true;
    return d->builtin_func == NULL && !(d->attributes & ATTR_NUMERICFUNCTION);
}

static ZeroTestResult pzq_decide_abstracted(const Expr* e);

/* Do two opaque calls denote the SAME quantity?  An inert f is still a
 * function: f[a] == f[b] whenever a == b, so the atom key is the head plus the
 * VALUES of the arguments, not their spelling.  Keyed on spelling alone,
 * PossibleZeroQ[f[2 (x + 1)] - f[2 x + 2]] was False, and so was the residual
 * of every characteristic invariant whose Leibniz derivative re-spells the
 * argument (f[c + a x + b y] against f[c + b (a x/b + y)], the PDE corpus's
 * pde-804).  Same head, same arity, and each argument pair proved equal by the
 * same abstracted decision procedure; structural equality is the fast path. */
static bool pzq_same_atom(const Expr* a, const Expr* b) {
    if (expr_eq((Expr*)a, (Expr*)b)) return true;
    if (a->type != EXPR_FUNCTION || b->type != EXPR_FUNCTION) return false;
    size_t n = a->data.function.arg_count;
    if (n == 0 || n != b->data.function.arg_count) return false;
    if (!expr_eq(a->data.function.head, b->data.function.head)) return false;
    for (size_t i = 0; i < n; i++) {
        const Expr* ai = a->data.function.args[i];
        const Expr* bi = b->data.function.args[i];
        if (expr_eq((Expr*)ai, (Expr*)bi)) continue;
        Expr* diff = eval_and_free(expr_new_function(expr_new_symbol(SYM_Subtract),
                         (Expr*[]){ expr_copy((Expr*)ai), expr_copy((Expr*)bi) }, 2));
        ZeroTestResult r = pzq_decide_abstracted(diff);
        expr_free(diff);
        if (r != ZERO_TEST_TRUE) return false;
    }
    return true;
}

/* Owned copy of `e` with every opaque call replaced; NULL if none. */
static Expr* pzq_abstract(const Expr* e, PzqAtoms* at) {
    if (e->type != EXPR_FUNCTION) return NULL;
    /* A pure function is code, not an independent quantity: Function has no
     * builtin (the evaluator applies it), so it read as opaque and
     * RootSum[1/256 + #^4 &, Log[1 + 64 v^(1/4) #] # &] lost its v.  Kept
     * whole, body included (its slots are not variables either). */
    if (e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Function)
        return NULL;
    if (pzq_opaque_call(e)) {
        for (size_t i = 0; i < at->n; i++)
            if (pzq_same_atom(at->from[i], e)) return expr_copy(at->to[i]);
        if (at->n == at->cap) {
            at->cap = at->cap ? at->cap * 2 : 8;
            at->from = (Expr**)realloc(at->from, at->cap * sizeof(Expr*));
            at->to = (Expr**)realloc(at->to, at->cap * sizeof(Expr*));
        }
        char name[64];
        snprintf(name, sizeof name, "PossibleZeroQ`Private`atom%zu", at->n + 1);
        at->from[at->n] = expr_copy((Expr*)e);
        at->to[at->n] = expr_new_symbol(intern_symbol(name));
        at->n++;
        return expr_copy(at->to[at->n - 1]);
    }
    size_t n = e->data.function.arg_count;
    Expr** args = NULL;
    for (size_t i = 0; i < n; i++) {
        Expr* r = pzq_abstract(e->data.function.args[i], at);
        if (r && !args) {
            args = (Expr**)malloc(sizeof(Expr*) * n);
            for (size_t j = 0; j < i; j++) args[j] = expr_copy(e->data.function.args[j]);
        }
        if (args) args[i] = r ? r : expr_copy(e->data.function.args[i]);
    }
    if (!args) return NULL;
    Expr* out = expr_new_function(expr_copy(e->data.function.head), args, n);
    free(args);
    return out;
}

static void pzq_atoms_free(PzqAtoms* at) {
    for (size_t i = 0; i < at->n; i++) { expr_free(at->from[i]); expr_free(at->to[i]); }
    free(at->from);
    free(at->to);
}

/* zero_test_decide on `e` with its opaque calls abstracted to fresh symbols,
 * exactly as a top-level unconditional PossibleZeroQ call does. */
static ZeroTestResult pzq_decide_abstracted(const Expr* e) {
    PzqAtoms atoms = { NULL, NULL, 0, 0 };
    Expr* abstracted = pzq_abstract(e, &atoms);
    ZeroTestResult r = zero_test_decide(abstracted ? abstracted : e);
    if (abstracted) expr_free(abstracted);
    pzq_atoms_free(&atoms);
    return r;
}

Expr* builtin_possible_zero_q(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return NULL;

    Expr* arg = res->data.function.args[0];

    /* Manual threading over the FIRST argument only (PossibleZeroQ is no longer
     * ATTR_LISTABLE, so the assumption argument is never mis-threaded): a list
     * first argument maps element-wise while every option/assumption argument
     * is broadcast unchanged. */
    if (arg->type == EXPR_FUNCTION && arg->data.function.head &&
        arg->data.function.head->type == EXPR_SYMBOL &&
        arg->data.function.head->data.symbol.name == SYM_List) {
        size_t n = arg->data.function.arg_count;
        Expr** out = malloc(sizeof(Expr*) * (n > 0 ? n : 1));
        for (size_t i = 0; i < n; ++i) {
            /* Rebuild PossibleZeroQ[elem, <rest>] and evaluate so each element
             * re-enters this builtin with the same options. */
            Expr** call_args = malloc(sizeof(Expr*) * argc);
            call_args[0] = expr_copy(arg->data.function.args[i]);
            for (size_t j = 1; j < argc; ++j)
                call_args[j] = expr_copy(res->data.function.args[j]);
            Expr* call = expr_new_function(expr_new_symbol(SYM_PossibleZeroQ),
                                           call_args, argc);
            free(call_args);
            out[i] = evaluate(call);
            expr_free(call);
        }
        Expr* list = expr_new_function(expr_new_symbol(SYM_List), out, n);
        free(out);
        return list;
    }

    /* Parse the trailing arguments: a Rule[Assumptions, X] option (which
     * overrides $Assumptions) and/or a single positional assumption (which is
     * And-combined with $Assumptions). */
    Expr* opt_assumptions = NULL;
    Expr* positional_assum = NULL;
    for (size_t i = 1; i < argc; ++i) {
        Expr* a = res->data.function.args[i];
        if (is_rule_with_lhs(a, "Assumptions")) {
            opt_assumptions = a->data.function.args[1];
        } else if (positional_assum == NULL) {
            positional_assum = a;
        }
    }

    /* Effective assumption: option overrides; else positional && $Assumptions;
     * else $Assumptions. A list assumption is normalised to a conjunction by
     * assume_ctx_from_expr / And evaluation. */
    Expr* effective;
    if (opt_assumptions) {
        if (positional_assum) {
            effective = pzq_combine_with_and(expr_copy(positional_assum),
                                             expr_copy(opt_assumptions));
        } else {
            /* evaluate() borrows its argument (does an internal copy and never
             * frees it), so the copy must be freed here or it leaks. */
            Expr* oc = expr_copy(opt_assumptions);
            effective = evaluate(oc);
            expr_free(oc);
        }
    } else {
        Expr* dollar = pzq_read_dollar_assumptions();
        if (positional_assum)
            effective = pzq_combine_with_and(expr_copy(positional_assum), dollar);
        else
            effective = dollar;
    }

    /* A generalized function is decided as a distribution BEFORE the opaque-
     * call abstraction below, which would otherwise turn DiracDelta[t - a] and
     * HeavisideTheta[t - a] into free symbols and sample them as numbers.  Its
     * scalar sub-questions (sifted coefficients, interval pieces) go through the
     * same abstract-then-decide path as a top-level call.  Unconditional only:
     * under assumptions the constrained sampler stays the arbiter. */
    if (assumption_is_trivial(effective) && zt_has_distribution(arg)) {
        ZeroTestResult rd = zt_decide_distribution(arg, pzq_decide_abstracted);
        if (rd != ZERO_TEST_UNKNOWN) {
            expr_free(effective);
            return pzq_bool(rd);
        }
    }

    /* Opaque calls become fresh symbols, in the expression and in the
     * assumptions alike (f[a] > 0 must still constrain the atom). */
    PzqAtoms atoms = { NULL, NULL, 0, 0 };
    Expr* abstracted = pzq_abstract(arg, &atoms);
    Expr* target = abstracted ? abstracted : arg;
    if (abstracted && effective) {
        Expr* ea = pzq_abstract(effective, &atoms);
        if (ea) { expr_free(effective); effective = ea; }
    }

    ZeroTestResult r;
    if (assumption_is_trivial(effective)) {
        r = zero_test_decide(target);              /* legacy path, no ctx */
    } else {
        AssumeCtx* ctx = assume_ctx_from_expr(effective);
        r = zero_test_decide_assuming(target, ctx);
        assume_ctx_free(ctx);                      /* safe no-op if ctx == NULL */
    }
    expr_free(effective);
    if (abstracted) expr_free(abstracted);
    pzq_atoms_free(&atoms);
    return pzq_bool(r);
}

void zero_test_init(void) {
    symtab_add_builtin("PossibleZeroQ", builtin_possible_zero_q);
    SymbolDef* def = symtab_get_def("PossibleZeroQ");
    /* NOT ATTR_LISTABLE: the first argument is threaded manually inside the
     * builtin so a positional list assumption is treated as one conjunction
     * rather than mis-threaded against the expression list. */
    if (def) def->attributes |= ATTR_PROTECTED;
}
