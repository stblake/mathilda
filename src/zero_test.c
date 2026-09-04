/*
 * zero_test.c — PossibleZeroQ: hybrid symbolic-numeric zero recognition.
 *
 * Pipeline (early exit at any stage that yields a definite verdict):
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
 *             implies "True".
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
 *
 * See ZERO_RECOGNISE_PLAN.md for design notes and references.
 */
#include "zero_test.h"

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
#define ZT_CONFIRM_SAMPLES 4

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

/* ------------------------------------------------------------------ */
/*  Forward declarations                                              */
/* ------------------------------------------------------------------ */

static ZeroTestResult decide_structural(const Expr* e);
static ZeroTestResult decide_rational(const Expr* e);
static ZeroTestResult decide_numeric(const Expr* e);
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
    }

    return ZERO_TEST_UNKNOWN;
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

/* Extract |v| as a double. Handles Integer, BigInt, Real, MPFR, and
 * Complex[re, im] (returns hypot of magnitudes). Returns false if the
 * input is not a recognized numeric form (e.g. still contains symbols). */
static bool expr_abs_double(const Expr* e, double* out) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) { *out = fabs((double)e->data.integer); return true; }
    if (e->type == EXPR_REAL)    { *out = fabs(e->data.real);             return true; }
    if (e->type == EXPR_BIGINT)  { *out = fabs(mpz_get_d(e->data.bigint));return true; }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) {
        *out = fabs(mpfr_get_d(e->data.mpfr, MPFR_RNDN));
        /* mpfr_get_d underflows extremely small values to 0.0; this is
         * the behaviour we want when checking "is this approximately
         * zero?", so no special-case is needed. */
        return true;
    }
#endif
    Expr* re = NULL; Expr* im = NULL;
    if (is_complex((Expr*)e, &re, &im)) {
        double mr = 0.0, mi = 0.0;
        if (!expr_abs_double(re, &mr)) return false;
        if (!expr_abs_double(im, &mi)) return false;
        *out = hypot(mr, mi);
        return true;
    }
    int64_t rn = 0, rd = 1;
    if (is_rational((Expr*)e, &rn, &rd)) {
        if (rd == 0) return false;
        *out = fabs((double)rn / (double)rd);
        return true;
    }
    /* Rational[n, d] with non-int64 args (e.g. Rational[1.0, 1e+30] that
     * numericalize emits for bigint-denominator rationals at machine
     * precision). Recurse on components. */
    if (e->type == EXPR_FUNCTION && e->data.function.arg_count == 2 &&
        e->data.function.head && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Rational) {
        double mn = 0.0, md = 0.0;
        if (!expr_abs_double(e->data.function.args[0], &mn)) return false;
        if (!expr_abs_double(e->data.function.args[1], &md)) return false;
        if (md == 0.0) return false;
        *out = mn / md;
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

/* Estimate the operand-magnitude scale of `e` at machine precision. This
 * is the denominator in the cancellation-aware ambiguity threshold:
 * results below scale * 2^(-p/2) are indistinguishable from rounding
 * noise, so we treat them as "ambiguous" rather than zero.
 *
 * The estimate uses head-aware recursion so it captures *additive*
 * cancellation (Plus) without spuriously inflating purely multiplicative
 * forms (Times, Rational, Power), where no cancellation can occur:
 *   - Plus[a, b, …]    -> Σ |a_i|       (potential cancellation)
 *   - Times[a, b, …]   -> Π |a_i|       (exact product magnitude)
 *   - Power[a, b]      -> |a|^|b|       (treating exponent magnitude)
 *   - Rational[n, d]   -> |n| / |d|
 *   - Complex[re, im]  -> hypot(|re|, |im|)
 *   - generic f[args]  -> Σ |a_i|       (conservative)
 * Numeric leaves give |value|; named constants (Pi, E, …) get their
 * numeric value at machine precision; free symbols default to 1. */
static double magnitude_scale(const Expr* e) {
    if (!e) return 1.0;

    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL ||
        e->type == EXPR_BIGINT
#ifdef USE_MPFR
        || e->type == EXPR_MPFR
#endif
        ) {
        double m = 0.0;
        if (expr_abs_double(e, &m)) return m;
        return 1.0;
    }

    if (e->type == EXPR_SYMBOL) {
        if (is_known_constant(e->data.symbol.name)) {
            Expr* n = numericalize(e, numeric_machine_spec());
            double m = 0.0;
            if (n && expr_abs_double(n, &m)) { expr_free(n); return m > 0.0 ? m : 1.0; }
            if (n) expr_free(n);
        }
        return 1.0;
    }

    if (e->type == EXPR_FUNCTION) {
        Expr* re = NULL; Expr* im = NULL;
        if (is_complex((Expr*)e, &re, &im)) {
            double mr = magnitude_scale(re);
            double mi = magnitude_scale(im);
            return hypot(mr, mi);
        }
        int64_t rn = 0, rd = 1;
        if (is_rational((Expr*)e, &rn, &rd) && rd != 0) {
            return fabs((double)rn / (double)rd);
        }

        const Expr* head = e->data.function.head;
        const char* head_name = (head && head->type == EXPR_SYMBOL)
                                ? head->data.symbol.name : NULL;
        size_t argc = e->data.function.arg_count;

        if (head_name == SYM_Plus) {
            double sum = 0.0;
            for (size_t i = 0; i < argc; ++i) sum += magnitude_scale(e->data.function.args[i]);
            return sum > 0.0 ? sum : 1.0;
        }
        if (head_name == SYM_Times) {
            double prod = 1.0;
            for (size_t i = 0; i < argc; ++i) prod *= magnitude_scale(e->data.function.args[i]);
            if (!isfinite(prod) || prod <= 0.0) return 1.0;
            return prod;
        }
        if (head_name == SYM_Power && argc == 2) {
            double base = magnitude_scale(e->data.function.args[0]);
            if (base <= 0.0) base = 1.0;
            /* Use the SIGNED exponent value, not its magnitude: Power[x, -1] has
             * magnitude 1/x, not x.  Scoring a reciprocal/denominator as x^|−1|
             * grossly inflates the operand scale of any expression with a
             * denominator, which then drowns a genuine non-zero below the
             * cancellation noise floor (a Schwartz-Zippel false zero). */
            double expv;
            if (!expr_signed_double(e->data.function.args[1], &expv))
                expv = magnitude_scale(e->data.function.args[1]); /* symbolic exp */
            double v = pow(base, expv);
            if (!isfinite(v) || v <= 0.0) return base; /* fall back */
            return v;
        }
        if (head_name == SYM_Rational && argc == 2) {
            double n = magnitude_scale(e->data.function.args[0]);
            double d = magnitude_scale(e->data.function.args[1]);
            if (d <= 0.0) return n > 0.0 ? n : 1.0;
            double v = n / d;
            if (!isfinite(v)) return 1.0;
            return v > 0.0 ? v : 1.0;
        }
        if (head_name == SYM_Sqrt && argc == 1) {
            double v = magnitude_scale(e->data.function.args[0]);
            return v > 0.0 ? sqrt(v) : 1.0;
        }

        /* Generic function (Sin, Cos, Log, Exp, user-defined, …): treat
         * the result magnitude as approximately the sum of operand
         * magnitudes.  Adequate for the threshold heuristic. */
        double sum = 0.0;
        for (size_t i = 0; i < argc; ++i) sum += magnitude_scale(e->data.function.args[i]);
        return sum > 0.0 ? sum : 1.0;
    }

    return 1.0;
}

/* "Definitely non-zero" threshold. At precision p, anything that exceeds
 * scale * 2^(-p + ZT_NONZERO_HEADROOM) is well above the rounding-noise
 * floor and is a confident False verdict. Tuned to forgive ~30 bits of
 * cancellation so transcendental identities that cancel down to machine
 * epsilon still survive the first rung. */
#define ZT_NONZERO_HEADROOM 30
static double nonzero_threshold(double scale, long p_bits) {
    long shift = p_bits - ZT_NONZERO_HEADROOM;
    if (shift <= 0) return INFINITY;          /* never reject at trivial precision */
    return scale * ldexp(1.0, -(int)shift);
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

/* Numericalize at the given precision, returning a freshly allocated
 * Expr*. spec.bits == 53 means machine; anything else uses MPFR (if
 * available). Returns NULL on failure or when MPFR is not compiled in
 * and an MPFR rung is requested. */
static Expr* numericalize_at(const Expr* e, long bits) {
    if (bits <= 53) {
        return numericalize(e, numeric_machine_spec());
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
    if (!e) return 1.0;
    if (e->type == EXPR_FUNCTION && e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Plus) {
        double sum = 0.0;
        for (size_t i = 0; i < e->data.function.arg_count; ++i) {
            Expr* n = numericalize(e->data.function.args[i], spec);
            if (n) {
                double m = 0.0;
                if (expr_abs_double(n, &m) && isfinite(m)) sum += m;
                expr_free(n);
            }
        }
        if (sum > 0.0 && isfinite(sum)) return sum;
        /* Fall through to the static estimate. */
    }
    double m = magnitude_scale(e);
    if (!isfinite(m) || m <= 0.0) return 1.0;
    return m;
}

/* Test the expression at a single precision rung. Returns FALSE if the
 * residual decisively exceeds the noise floor; TRUE if the residual is
 * indistinguishable from zero (very small relative to scale); UNKNOWN
 * if the numericalize couldn't reduce. *out_mag receives the residual
 * magnitude when non-UNKNOWN, otherwise 0.0. *out_scale (when non-NULL)
 * receives the operand-magnitude scale used for the threshold.
 *
 * `in_scale`: when >= 0, the caller supplies a precomputed operand scale and
 * we SKIP the (expensive) per-argument numericalization in magnitude_scale_at.
 * The operand magnitude is essentially precision-independent (it is a sum/
 * product of operand values, not a cancelling residual), so a scale computed
 * once at machine precision is reused for every higher rung — this removes a
 * redundant full second numericalization pass at each MPFR rung, the dominant
 * cost on large (tens-of-thousands-of-leaves) antiderivative round-trips.
 * Pass a negative value to compute the scale at this rung's precision. */
static ZeroTestResult evaluate_rung(const Expr* e, long bits,
                                    double* out_mag, double* out_scale,
                                    double in_scale) {
    *out_mag = 0.0;
    if (out_scale) *out_scale = 1.0;
    Expr* z = numericalize_at(e, bits);
    if (!z) return ZERO_TEST_UNKNOWN;
    if (expr_is_infinity(z)) {
        /* A pole: definitively non-zero. Report an infinite residual so every
         * caller's "obvious non-zero" gate settles it FALSE. */
        expr_free(z);
        *out_mag = INFINITY;
        if (out_scale) *out_scale = 1.0;
        return ZERO_TEST_FALSE;
    }
    if (!is_pure_numeric(z)) { expr_free(z); return ZERO_TEST_UNKNOWN; }
    double mag = 0.0;
    bool ok = expr_abs_double(z, &mag);
    expr_free(z);
    if (!ok || !isfinite(mag)) return ZERO_TEST_UNKNOWN;
    *out_mag = mag;

    double scale = (in_scale >= 0.0) ? in_scale
                                     : magnitude_scale_at(e, spec_at_bits(bits));
    if (out_scale) *out_scale = scale;
    double tol = nonzero_threshold(scale, bits);
    if (mag > tol) return ZERO_TEST_FALSE;
    return ZERO_TEST_TRUE;
}

/* Numeric Stage 2.  Strategy:
 *   1. Machine precision. A residual that is a non-trivial fraction of the
 *      operand scale (above the ZT_OBVIOUS_NONZERO_BITS gate) cannot be
 *      rounding noise and settles FALSE immediately — the fast path for typical
 *      non-zeros. A smaller residual is AMBIGUOUS: it may be a true zero whose
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
static ZeroTestResult decide_numeric(const Expr* e) {
    double mag = 0.0, scale = 1.0;
    ZeroTestResult r = evaluate_rung(e, PRECISION_LADDER[0], &mag, &scale, -1.0);
    if (r == ZERO_TEST_UNKNOWN) return ZERO_TEST_UNKNOWN;

    if (scale > 0.0 && mag > scale * ldexp(1.0, -ZT_OBVIOUS_NONZERO_BITS))
        return ZERO_TEST_FALSE;

    double prev_mag = mag;
    bool honored = false;
    for (int i = 1; i < PRECISION_LADDER_LEN; ++i) {
        double m = 0.0;
        /* Reuse the machine-precision operand scale (precision-independent) so
         * each higher rung numericalizes the tree ONCE, not twice. */
        ZeroTestResult rr = evaluate_rung(e, PRECISION_LADDER[i], &m, NULL,
                                          scale > 0.0 ? scale : -1.0);
        if (rr == ZERO_TEST_UNKNOWN) {
            /* MPFR path unavailable beyond here — accept the lenient machine
             * verdict (the rung-0 residual was below the non-zero gate). */
            return ZERO_TEST_TRUE;
        }
        if (m < prev_mag * 0.5) {
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
            if (scale > 0.0 && m < scale * ldexp(1.0, -ZT_DEEP_ZERO_BITS))
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
        if (is_algebraic_expr(e) && scale > 0.0 &&
            m > scale * ldexp(1.0, -ZT_ALG_NONZERO_BITS))
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
static ZeroTestResult screen_point(const Expr* e) {
    double mag = 0.0, scale = 1.0;
    ZeroTestResult r = evaluate_rung(e, PRECISION_LADDER[0], &mag, &scale, -1.0);
    if (r == ZERO_TEST_UNKNOWN) return ZERO_TEST_UNKNOWN;
    if (scale > 0.0 && mag > scale * ldexp(1.0, -ZT_OBVIOUS_NONZERO_BITS))
        return ZERO_TEST_FALSE;
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
static Expr* sample_random_value(void) {
    int64_t num_bound = (int64_t)1 << ZT_NUMERATOR_BITS;
    int64_t den_bound = (int64_t)1 << ZT_DENOMINATOR_BITS;

    /* Sample magnitude is bounded BELOW by 1 (an integer part in [1, num_bound])
     * plus a full-granularity fractional part for a rich distinct-value set.
     * Keeping |value| >= 1 is essential for polynomial/algebraic identity
     * testing: a sample drawn from (-1, 1) drives a high-degree monomial
     * (u^2, u^3, ...) far below the operand-magnitude scale, so a GENUINE
     * non-zero such as Sqrt[2] u^3 collapses into the rounding-noise band and is
     * misread as an identity (Schwartz-Zippel false positive).  An analytic
     * identity that holds on [1, num_bound] holds everywhere it is analytic
     * (identity theorem), so excluding the small-magnitude shell loses no
     * decision power while removing the false-zero failure mode.  The numerator
     * magnitude stays moderate (<= ~num_bound) to avoid overflowing Exp/Gamma. */
    int64_t whole  = draw_int_range(1, num_bound);
    int64_t frac_n = draw_int_range(0, den_bound - 1);
    double  val    = (double)whole + (double)frac_n / (double)den_bound;
    if (draw_int_range(0, 1)) val = -val;
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
 * the unconstrained magnitude draw keeps |v| >= 1 and so never lands on 0.
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

/* Draw one real value honouring `c`.  `integer` rounds to the integer grid. */
static double sample_channel(const Channel* c, bool integer) {
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
    /* Unconstrained magnitude: |v| >= 1, moderate, sign-honoured. */
    int64_t num_bound = (int64_t)1 << ZT_NUMERATOR_BITS;
    double v;
    if (integer) {
        v = (double)draw_int_range(1, num_bound);
    } else {
        int64_t den_bound = (int64_t)1 << ZT_DENOMINATOR_BITS;
        int64_t whole  = draw_int_range(1, num_bound);
        int64_t frac_n = draw_int_range(0, den_bound - 1);
        v = (double)whole + (double)frac_n / (double)den_bound;
    }
    if (c->sign < 0) v = -v;
    else if (c->sign == 0 && draw_int_range(0, 1)) v = -v;
    return v;
}

/* Draw a conforming sample value for one free symbol.  spec == NULL restores
 * the legacy real-only draw byte-for-byte (see sample_random_value). */
static Expr* sample_random_value_spec(const SampleSpec* spec) {
    if (!spec) return sample_random_value();
    if (spec->domain == SDOM_INT)
        return expr_new_integer((int64_t)sample_channel(&spec->val, true));
    if (spec->domain == SDOM_CPLX) {
        double re = sample_channel(&spec->val, false);
        double im = sample_channel(&spec->im, false);
        Expr* parts[2] = { expr_new_real(re), expr_new_real(im) };
        return expr_new_function(expr_new_symbol(SYM_Complex), parts, 2);
    }
    return expr_new_real(sample_channel(&spec->val, false));
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

    for (size_t i = 0; ctx && i < ctx->count; ++i) {
        const Expr* f = ctx->facts[i];
        if (!f || f->type != EXPR_FUNCTION || !f->data.function.head ||
            f->data.function.head->type != EXPR_SYMBOL) continue;
        const char* h = f->data.function.head->data.symbol.name;
        size_t argc = f->data.function.arg_count;

        if (h == SYM_Element && argc == 2) {
            if (channel_of_operand(f->data.function.args[0], sym) != 0) continue;
            const Expr* dom = f->data.function.args[1];
            if (!dom || dom->type != EXPR_SYMBOL) continue;
            const char* dn = dom->data.symbol.name;
            if      (strcmp(dn, "Integers") == 0)            seen_int = true;
            else if (strcmp(dn, "PositiveIntegers") == 0)    { seen_int = true; ch_bare.has_lo = true; ch_bare.lo = 0.0; }
            else if (strcmp(dn, "NonnegativeIntegers") == 0) { seen_int = true; ch_bare.has_lo = true; ch_bare.lo = 0.0; }
            else if (strcmp(dn, "NegativeIntegers") == 0)    { seen_int = true; ch_bare.has_hi = true; ch_bare.hi = 0.0; }
            else if (strcmp(dn, "NonpositiveIntegers") == 0) { seen_int = true; ch_bare.has_hi = true; ch_bare.hi = 0.0; }
            else if (strcmp(dn, "Complexes") == 0)           seen_cplx = true;
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
                               bool screen) {
    Expr** vals = malloc(sizeof(Expr*) * nsyms);
    for (size_t i = 0; i < nsyms; ++i)
        vals[i] = sample_random_value_spec(specs ? &specs[i] : NULL);

    Expr* sub = substitute_symbols(e, syms, vals, nsyms);

    ZeroTestResult r = ZERO_TEST_UNKNOWN;
    if (sub) {
        r = decide_structural(sub);
        if (r == ZERO_TEST_UNKNOWN) {
            if (screen) {
                r = screen_point(sub);
                if (r == ZERO_TEST_UNKNOWN) r = decide_numeric(sub);
            } else {
                r = decide_numeric(sub);
            }
        }
        expr_free(sub);
    }

    for (size_t i = 0; i < nsyms; ++i) expr_free(vals[i]);
    free(vals);
    return r;
}

/* Salt mixed into the structural-hash seed so the sample distribution can be
 * re-tuned later (bump the salt) without colliding with any cached behaviour. */
#define ZT_SEED_SALT 0x5A3D9E1Bull

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

    SampleSpec* specs = NULL;
    if (ctx) {
        specs = malloc(sizeof(SampleSpec) * syms.count);
        if (specs)
            for (size_t i = 0; i < syms.count; ++i)
                specs[i] = extract_spec(ctx, syms.items[i]);
    }

    /* Seed the draw stream deterministically from the input's structural hash
     * so the verdict is a pure function of the input (no run-to-run flakiness).
     * The push/pop pair leaves the user's RandomInteger/SeedRandom stream
     * exactly as it was. */
    random_push_seed(expr_hash(e) ^ ZT_SEED_SALT);

    ZeroTestResult verdict = ZERO_TEST_TRUE;

    /* Phase A — screen: many cheap machine-precision points catch
     * branch-dependent non-zeros. A single decisively non-zero point settles
     * the whole test. */
    for (int trial = 0; trial < ZT_SCREEN_SAMPLES; ++trial) {
        ZeroTestResult r = sz_trial(e, syms.items, specs, syms.count, true);
        if (r == ZERO_TEST_FALSE)   { verdict = ZERO_TEST_FALSE;   goto done; }
        if (r == ZERO_TEST_UNKNOWN) { verdict = ZERO_TEST_UNKNOWN; goto done; }
        /* r == TRUE: zero-ish at machine precision, keep screening. */
    }

    /* Phase B — confirm: every screened point looked zero. Climb the full
     * ladder on a few fresh points to reject cancellation-hidden small
     * non-zeros before declaring a genuine identity. */
    for (int trial = 0; trial < ZT_CONFIRM_SAMPLES; ++trial) {
        ZeroTestResult r = sz_trial(e, syms.items, specs, syms.count, false);
        if (r == ZERO_TEST_FALSE)   { verdict = ZERO_TEST_FALSE;   goto done; }
        if (r == ZERO_TEST_UNKNOWN) { verdict = ZERO_TEST_UNKNOWN; goto done; }
    }

done:
    random_pop_seed();
    free(specs);
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

ZeroTestResult zero_test_decide(const Expr* e) {
    ZeroTestResult r;

    r = decide_structural(e);
    if (r != ZERO_TEST_UNKNOWN) return r;

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
        r = decide_numeric(e);
        if (r != ZERO_TEST_UNKNOWN) return r;
        return ZERO_TEST_UNKNOWN;
    }

    return decide_schwartz_zippel(e);
}

ZeroTestResult zero_test_decide_assuming(const Expr* e, const struct AssumeCtx* ctx) {
    /* No usable assumptions → the legacy path, byte-for-byte. */
    if (!ctx || ctx->count == 0) return zero_test_decide(e);

    /* A structurally decided literal (0, a non-zero constant, Complex[0,0], …)
     * is unconditional, so its verdict holds under any assumption. */
    ZeroTestResult r = decide_structural(e);
    if (r != ZERO_TEST_UNKNOWN) return r;

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

    /* Closed-form constant: assumptions are irrelevant to a definite number. */
    return decide_numeric(e);
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

    if (assumption_is_trivial(effective)) {
        expr_free(effective);
        return pzq_bool(zero_test_decide(arg));   /* legacy path, no ctx */
    }

    AssumeCtx* ctx = assume_ctx_from_expr(effective);
    expr_free(effective);
    ZeroTestResult r = zero_test_decide_assuming(arg, ctx);
    assume_ctx_free(ctx);                          /* safe no-op if ctx == NULL */
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
