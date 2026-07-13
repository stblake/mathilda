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
static bool           has_free_symbols(const Expr* e);
static bool           is_known_constant(const char* sym_name);

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
    if (z) return ZERO_TEST_TRUE;
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
    Expr* low  = expr_new_integer(lo);
    Expr* high = expr_new_integer(hi);
    Expr* list_args[] = {low, high};
    Expr* range = expr_new_function(expr_new_symbol(SYM_List), list_args, 2);
    Expr* call_args[] = {range};
    Expr* call = expr_new_function(expr_new_symbol(SYM_RandomInteger), call_args, 1);
    /* eval_and_free consumes `call`; plain evaluate() would leak the
     * RandomInteger[List[...]] wrapper, which it does not take ownership of. */
    Expr* r = eval_and_free(call);
    int64_t out = 0;
    if (r) {
        if (r->type == EXPR_INTEGER) out = r->data.integer;
        else if (r->type == EXPR_BIGINT) out = (int64_t)mpz_get_si(r->data.bigint);
        expr_free(r);
    }
    return out;
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
static ZeroTestResult sz_trial(const Expr* e, const char** syms, size_t nsyms,
                               bool screen) {
    Expr** vals = malloc(sizeof(Expr*) * nsyms);
    for (size_t i = 0; i < nsyms; ++i) vals[i] = sample_random_value();

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

static ZeroTestResult decide_schwartz_zippel(const Expr* e) {
    SymPtrSet syms; sps_init(&syms);
    collect_free(e, &syms);
    if (syms.count == 0) {
        sps_free(&syms);
        return ZERO_TEST_UNKNOWN;
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
        ZeroTestResult r = sz_trial(e, syms.items, syms.count, true);
        if (r == ZERO_TEST_FALSE)   { verdict = ZERO_TEST_FALSE;   goto done; }
        if (r == ZERO_TEST_UNKNOWN) { verdict = ZERO_TEST_UNKNOWN; goto done; }
        /* r == TRUE: zero-ish at machine precision, keep screening. */
    }

    /* Phase B — confirm: every screened point looked zero. Climb the full
     * ladder on a few fresh points to reject cancellation-hidden small
     * non-zeros before declaring a genuine identity. */
    for (int trial = 0; trial < ZT_CONFIRM_SAMPLES; ++trial) {
        ZeroTestResult r = sz_trial(e, syms.items, syms.count, false);
        if (r == ZERO_TEST_FALSE)   { verdict = ZERO_TEST_FALSE;   goto done; }
        if (r == ZERO_TEST_UNKNOWN) { verdict = ZERO_TEST_UNKNOWN; goto done; }
    }

done:
    random_pop_seed();
    sps_free(&syms);
    return verdict;
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

    r = decide_rational(e);
    /* Trust TRUE always; trust FALSE only for the rigorous pure-rational path
     * (decide_rational returns FALSE exclusively when is_pure_rational_function
     * holds, where the Q-normalization is exact and complete). */
    if (r != ZERO_TEST_UNKNOWN) return r;

    if (!has_free_symbols(e)) {
        r = decide_numeric(e);
        if (r != ZERO_TEST_UNKNOWN) return r;
        return ZERO_TEST_UNKNOWN;
    }

    return decide_schwartz_zippel(e);
}

Expr* builtin_possible_zero_q(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1) return NULL;

    Expr* arg = res->data.function.args[0];
    ZeroTestResult r = zero_test_decide(arg);

    /* Mathematica's documented behaviour: when uncertain, *assume* zero
     * and emit a PossibleZeroQ::ztest1 message. Mathilda currently has
     * no general message channel; we silently collapse UNKNOWN to True
     * to preserve the documented public-facing return value. */
    if (r == ZERO_TEST_FALSE) return expr_new_symbol(SYM_False);
    return expr_new_symbol(SYM_True);
}

void zero_test_init(void) {
    symtab_add_builtin("PossibleZeroQ", builtin_possible_zero_q);
    SymbolDef* def = symtab_get_def("PossibleZeroQ");
    if (def) def->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
}
