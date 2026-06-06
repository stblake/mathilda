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
 *             each free symbol with a random rational (drawn from Q[i] so
 *             branch cuts are tested too), recurse into Stage 2, and
 *             require ZT_NUM_SAMPLES independent confirmations.
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

/* Number of independent random-substitution trials for Stage 3. Each
 * trial reduces the false-positive probability geometrically (see
 * Schwartz–Zippel: false-positive ≤ (degree/|S|)^k for k samples). */
#define ZT_NUM_SAMPLES 4

/* Sampling ranges for Stage 3. Numerators are uniform on
 * [-2^NUMERATOR_BITS, 2^NUMERATOR_BITS]; denominators on
 * [1, 2^DENOMINATOR_BITS]. Imaginary part draws a coin and, if set,
 * picks a second rational of the same form. */
#define ZT_NUMERATOR_BITS    20
#define ZT_DENOMINATOR_BITS  16

/* Safety bits in the cancellation threshold. The IEEE-cancellation rule of
 * thumb is that a result is "ambiguous" if smaller than scale * 2^(-p/2);
 * we add a couple of extra bits as a guard against fma / parsing slop. */
#define ZT_AMBIGUITY_GUARD_BITS 4

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
        const char* name = e->data.symbol;
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
        const char* name = e->data.symbol;
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
    return z ? ZERO_TEST_TRUE : ZERO_TEST_UNKNOWN;
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
        e->data.function.head->data.symbol == SYM_Rational) {
        double mn = 0.0, md = 0.0;
        if (!expr_abs_double(e->data.function.args[0], &mn)) return false;
        if (!expr_abs_double(e->data.function.args[1], &md)) return false;
        if (md == 0.0) return false;
        *out = mn / md;
        return true;
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
        if (is_known_constant(e->data.symbol)) {
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
                                ? head->data.symbol : NULL;
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
            double exp  = magnitude_scale(e->data.function.args[1]);
            if (base <= 0.0) base = 1.0;
            double v = pow(base, exp);
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
        e->data.function.head->data.symbol == SYM_Rational) {
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
        e->data.function.head->data.symbol == SYM_Plus) {
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
 * magnitude when non-UNKNOWN, otherwise 0.0. */
static ZeroTestResult evaluate_rung(const Expr* e, long bits, double* out_mag) {
    *out_mag = 0.0;
    Expr* z = numericalize_at(e, bits);
    if (!z) return ZERO_TEST_UNKNOWN;
    if (!is_pure_numeric(z)) { expr_free(z); return ZERO_TEST_UNKNOWN; }
    double mag = 0.0;
    bool ok = expr_abs_double(z, &mag);
    expr_free(z);
    if (!ok || !isfinite(mag)) return ZERO_TEST_UNKNOWN;
    *out_mag = mag;

    double scale = magnitude_scale_at(e, spec_at_bits(bits));
    double tol = nonzero_threshold(scale, bits);
    if (mag > tol) return ZERO_TEST_FALSE;
    return ZERO_TEST_TRUE;
}

/* Numeric Stage 2.  Strategy:
 *   1. Evaluate at machine precision. Verdict FALSE if mag exceeds the
 *      lenient nonzero threshold; otherwise TRUE-so-far.
 *   2. Climb the precision ladder. At each rung either confirm shrinkage
 *      (mag dropped substantially → still consistent with zero, continue)
 *      or detect non-shrinkage (the higher-precision call returned ~the
 *      same value, meaning the underlying numericalize path doesn't
 *      actually use the requested precision — e.g. Sin/Cos of
 *      Complex[Real, Real] which Mathilda evaluates with double-only).
 *      In the latter case fall back to the previous verdict. */
static ZeroTestResult decide_numeric(const Expr* e) {
    double mag = 0.0;
    ZeroTestResult r = evaluate_rung(e, PRECISION_LADDER[0], &mag);
    if (r == ZERO_TEST_UNKNOWN) return ZERO_TEST_UNKNOWN;
    if (r == ZERO_TEST_FALSE)   return ZERO_TEST_FALSE;

    double prev_mag = mag;
    for (int i = 1; i < PRECISION_LADDER_LEN; ++i) {
        long bits = PRECISION_LADDER[i];
        double m = 0.0;
        ZeroTestResult rr = evaluate_rung(e, bits, &m);
        if (rr == ZERO_TEST_UNKNOWN) {
            /* MPFR path unavailable for this expression — accept the
             * machine-precision verdict. */
            return ZERO_TEST_TRUE;
        }
        if (rr == ZERO_TEST_FALSE) return ZERO_TEST_FALSE;
        /* A genuine zero shrinks geometrically with precision; if the
         * residual is essentially unchanged, the requested precision is
         * not being honoured downstream. Accept the prior verdict. */
        if (m > prev_mag * 0.5) return ZERO_TEST_TRUE;
        prev_mag = m;
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
    Expr* low  = expr_new_integer(lo);
    Expr* high = expr_new_integer(hi);
    Expr* list_args[] = {low, high};
    Expr* range = expr_new_function(expr_new_symbol("List"), list_args, 2);
    Expr* call_args[] = {range};
    Expr* call = expr_new_function(expr_new_symbol("RandomInteger"), call_args, 1);
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

/* Build a random (possibly complex) sample value. Real-valued leaves
 * (EXPR_REAL) are used rather than exact Rational[n, d] because several
 * Mathilda numeric heads (notably Log of Complex[Rational, …]) only take
 * the fast numeric path when the components are already Real. Imaginary
 * part is drawn with probability 1/2 so branch cuts get exercised. */
static Expr* sample_random_value(void) {
    int64_t num_bound = (int64_t)1 << ZT_NUMERATOR_BITS;
    int64_t den_bound = (int64_t)1 << ZT_DENOMINATOR_BITS;

    int64_t n_re = draw_int_range(-num_bound, num_bound);
    int64_t d_re = draw_int_range(1, den_bound);
    if (d_re == 0) d_re = 1;
    Expr* re = expr_new_real((double)n_re / (double)d_re);

    int64_t coin = draw_int_range(0, 1);
    if (coin == 0) return re;

    int64_t n_im = draw_int_range(-num_bound, num_bound);
    if (n_im == 0) n_im = 1;  /* keep imaginary part non-trivial */
    int64_t d_im = draw_int_range(1, den_bound);
    if (d_im == 0) d_im = 1;
    Expr* im = expr_new_real((double)n_im / (double)d_im);

    return make_complex(re, im);
}

/* Substitute every free symbol in `e` with an entry from `(syms, vals)`,
 * returning a fresh expression. Walks the tree directly so we don't pay
 * the pattern-matcher overhead of building rules. */
static Expr* substitute_symbols(const Expr* e, const char** syms, Expr** vals, size_t n) {
    if (!e) return NULL;

    if (e->type == EXPR_SYMBOL) {
        for (size_t i = 0; i < n; ++i) {
            if (e->data.symbol == syms[i]) return expr_copy(vals[i]);
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

static ZeroTestResult decide_schwartz_zippel(const Expr* e) {
    SymPtrSet syms; sps_init(&syms);
    collect_free(e, &syms);
    if (syms.count == 0) {
        sps_free(&syms);
        return ZERO_TEST_UNKNOWN;
    }

    ZeroTestResult verdict = ZERO_TEST_TRUE;
    for (int trial = 0; trial < ZT_NUM_SAMPLES; ++trial) {
        Expr** vals = malloc(sizeof(Expr*) * syms.count);
        for (size_t i = 0; i < syms.count; ++i) vals[i] = sample_random_value();

        /* Substitute without evaluating: the Mathilda evaluator would
         * eagerly numericalize, collapsing the top-level Plus into a
         * single Complex residue and discarding the operand magnitudes
         * decide_numeric relies on for its threshold. The structural
         * shape (Plus / Times / Sin / Cos / ...) is preserved through
         * substitute_symbols, and decide_numeric will numericalize the
         * whole expression at each ladder rung. */
        Expr* sub = substitute_symbols(e, syms.items, vals, syms.count);

        ZeroTestResult r = ZERO_TEST_UNKNOWN;
        if (sub) {
            r = decide_structural(sub);
            if (r == ZERO_TEST_UNKNOWN) r = decide_numeric(sub);
            expr_free(sub);
        }

        for (size_t i = 0; i < syms.count; ++i) expr_free(vals[i]);
        free(vals);

        if (r == ZERO_TEST_FALSE) { verdict = ZERO_TEST_FALSE; break; }
        if (r == ZERO_TEST_UNKNOWN) { verdict = ZERO_TEST_UNKNOWN; break; }
        /* r == TRUE: continue, need k consecutive hits */
    }

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

    r = decide_rational(e);
    if (r == ZERO_TEST_TRUE) return r;   /* never trust False from Stage 1 alone */

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
    if (r == ZERO_TEST_FALSE) return expr_new_symbol("False");
    return expr_new_symbol("True");
}

void zero_test_init(void) {
    symtab_add_builtin("PossibleZeroQ", builtin_possible_zero_q);
    SymbolDef* def = symtab_get_def("PossibleZeroQ");
    if (def) def->attributes |= ATTR_LISTABLE | ATTR_PROTECTED;
}
