/* Mathilda — numeric evaluation implementation.
 *
 * See numeric.h for the module-level overview and extensibility notes.
 *
 * This file implements `N[expr]` / `N[expr, prec]`. Phase 1 targets
 * machine-precision IEEE doubles; Phase 2 (gated behind USE_MPFR) adds
 * MPFR arbitrary precision. Phase-2 extension points are marked with
 * an inline "Phase 2" marker so the eventual additions are obvious.
 */

#include "numeric.h"

#include "arithmetic.h"
#include "attr.h"
#include "eval.h"
#include "symtab.h"
#include "sym_names.h"
#include "root_numeric.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>

/* M_PI and M_E are POSIX/GNU extensions, not part of C99. glibc hides them
 * under -std=c99, so provide portable fallbacks here. */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E
#define M_E  2.71828182845904523536
#endif

/* ------------------------------------------------------------------------
 *  Precision unit conversion
 *
 *  Mathematica talks about decimal digits of precision; MPFR talks about
 *  binary bits. The ratio is log2(10) ≈ 3.3219280948873626. We round up on
 *  conversion so the stored value always carries at least the requested
 *  information content.
 * ---------------------------------------------------------------------- */

static const double LOG2_10 = 3.3219280948873626;

long numeric_digits_to_bits(double digits) {
    if (digits <= 0.0) return 1;
    return (long)ceil(digits * LOG2_10);
}

double numeric_bits_to_digits(long bits) {
    if (bits <= 0) return 0.0;
    return (double)bits / LOG2_10;
}

/* ------------------------------------------------------------------------
 *  Constants registry
 *
 *  To add a new named constant: append a row to `kConstants` below. The
 *  machine value is the double you want `N[Name]` to produce. The optional
 *  `mpfr_fill` callback (Phase 2) writes the constant at arbitrary
 *  precision; leave NULL if no MPFR implementation exists yet.
 *
 *  Passing through unevaluated symbols (Infinity, True, etc.) is handled
 *  by the default path in `numericalize` — they are simply not listed here.
 * ---------------------------------------------------------------------- */

#ifdef USE_MPFR
#include <mpfr.h>
typedef void (*NumericMpfrFill)(mpfr_t out, mpfr_prec_t bits);
#else
typedef void (*NumericMpfrFill)(void);  /* placeholder for ABI stability */
#endif

typedef struct {
    const char*      name;
    double           machine_value;
    NumericMpfrFill  mpfr_fill;     /* NULL if not implemented */
} NumericConstant;

/* Forward declarations for future MPFR fillers. Not implemented in Phase 1. */
/* Phase 2 */
#ifdef USE_MPFR
static void fill_mpfr_pi(mpfr_t out, mpfr_prec_t bits);
static void fill_mpfr_e (mpfr_t out, mpfr_prec_t bits);
static void fill_mpfr_eulergamma(mpfr_t out, mpfr_prec_t bits);
static void fill_mpfr_catalan(mpfr_t out, mpfr_prec_t bits);
static void fill_mpfr_goldenratio(mpfr_t out, mpfr_prec_t bits);
static void fill_mpfr_degree(mpfr_t out, mpfr_prec_t bits);
static void fill_mpfr_goldenangle(mpfr_t out, mpfr_prec_t bits);
static void fill_mpfr_glaisher(mpfr_t out, mpfr_prec_t bits);
static void fill_mpfr_khinchin(mpfr_t out, mpfr_prec_t bits);
#define MPFR_FILL(fn) fn
#else
#define MPFR_FILL(fn) NULL
#endif

static const NumericConstant kConstants[] = {
    { "Pi",          M_PI,                                       MPFR_FILL(fill_mpfr_pi)          },
    { "E",           M_E,                                        MPFR_FILL(fill_mpfr_e)           },
    { "EulerGamma",  0.5772156649015328606065120900824024310421, MPFR_FILL(fill_mpfr_eulergamma)  },
    { "Catalan",     0.9159655941772190150546035149323841107741, MPFR_FILL(fill_mpfr_catalan)     },
    { "GoldenRatio", 1.6180339887498948482045868343656381177203, MPFR_FILL(fill_mpfr_goldenratio) },
    { "Degree",      M_PI / 180.0,                               MPFR_FILL(fill_mpfr_degree)      },
    { "GoldenAngle", 2.3999632297286533222315555066336138531249, MPFR_FILL(fill_mpfr_goldenangle) },
    { "Glaisher",    1.2824271291006226368753425688697917277676, MPFR_FILL(fill_mpfr_glaisher)    },
    { "Khinchin",    2.6854520010653064453097148354817956938203, MPFR_FILL(fill_mpfr_khinchin)    },
};
static const size_t kConstantCount = sizeof(kConstants) / sizeof(kConstants[0]);

static const NumericConstant* find_constant(const char* name) {
    if (!name) return NULL;
    for (size_t i = 0; i < kConstantCount; ++i) {
        if (strcmp(kConstants[i].name, name) == 0) return &kConstants[i];
    }
    return NULL;
}

bool numeric_constant_machine_value(const char* name, double* out) {
    const NumericConstant* c = find_constant(name);
    if (!c) return false;
    if (out) *out = c->machine_value;
    return true;
}

#ifdef USE_MPFR
/* ------------------------------------------------------------------------
 *  MPFR propagation helpers used by Plus/Times/Power/trig/hyperbolic/log.
 *
 *  These keep each builtin's per-function code focused on the math; the
 *  precision arithmetic and type juggling is centralized here.
 * ---------------------------------------------------------------------- */

/* Max precision of the MPFR values contained in `e`, or 0 if `e` carries
 * none. Descends through Complex[...] so that e.g. Complex[3.14`50, 0]
 * contributes 50 digits. */
static long expr_max_mpfr_prec(const Expr* e) {
    if (!e) return 0;
    if (e->type == EXPR_MPFR) return (long)mpfr_get_prec(e->data.mpfr);
    if (e->type == EXPR_FUNCTION) {
        Expr *re, *im;
        if (is_complex((Expr*)e, &re, &im)) {
            long a = expr_max_mpfr_prec(re);
            long b = expr_max_mpfr_prec(im);
            return a > b ? a : b;
        }
    }
    return 0;
}

long numeric_combined_bits(const Expr* a, const Expr* b, long default_bits) {
    long pa = expr_max_mpfr_prec(a);
    long pb = expr_max_mpfr_prec(b);
    long m  = pa > pb ? pa : pb;
    if (m <= 0) m = default_bits;
    if (m <= 0) m = 53;  /* safety floor: IEEE double precision */
    return m;
}

bool numeric_any_mpfr(const Expr* a, const Expr* b) {
    return expr_max_mpfr_prec(a) > 0 || expr_max_mpfr_prec(b) > 0;
}

/* Minimum precision among the inexact (Real/MPFR) leaves of `e`, descending
 * into Complex[...]; 0 if there is no inexact leaf. A machine Real counts as
 * 53 bits, an MPFR as its own precision; exact atoms impose no constraint. */
long numeric_min_inexact_bits(const Expr* e) {
    if (!e) return 0;
    if (e->type == EXPR_REAL) return 53;
    if (e->type == EXPR_MPFR) return (long)mpfr_get_prec(e->data.mpfr);
    if (e->type == EXPR_FUNCTION) {
        Expr *re, *im;
        if (is_complex((Expr*)e, &re, &im)) {
            long a = numeric_min_inexact_bits(re);
            long b = numeric_min_inexact_bits(im);
            if (a == 0) return b;          /* re exact: only im constrains */
            if (b == 0) return a;          /* im exact: only re constrains */
            return a < b ? a : b;          /* both inexact: tighter wins */
        }
    }
    return 0;                              /* exact atom: no constraint */
}

typedef int (*MpfrBinaryOp)(mpfr_t, const mpfr_t, const mpfr_t, mpfr_rnd_t);

/* Generic MPFR binary dispatcher. Handles the common shape: extract real
 * MPFR operands from `a` and `b`, apply `op`, return a fresh EXPR_MPFR.
 * Returns NULL if either operand isn't a pure real (Complex → caller
 * handles separately). */
static Expr* mpfr_binop(const Expr* a, const Expr* b, long default_bits,
                        MpfrBinaryOp op) {
    long bits = numeric_combined_bits(a, b, default_bits);
    mpfr_t ra, ia, rb, ib, out;
    mpfr_init2(ra, bits);  mpfr_init2(ia, bits);
    mpfr_init2(rb, bits);  mpfr_init2(ib, bits);
    bool ok_a = get_approx_mpfr(a, ra, ia, NULL);
    bool ok_b = get_approx_mpfr(b, rb, ib, NULL);
    if (!ok_a || !ok_b || !mpfr_zero_p(ia) || !mpfr_zero_p(ib)) {
        mpfr_clear(ra); mpfr_clear(ia);
        mpfr_clear(rb); mpfr_clear(ib);
        return NULL;
    }
    mpfr_init2(out, bits);
    op(out, ra, rb, MPFR_RNDN);
    mpfr_clear(ra); mpfr_clear(ia);
    mpfr_clear(rb); mpfr_clear(ib);
    return expr_new_mpfr_move(out);
}

Expr* numeric_mpfr_add(const Expr* a, const Expr* b, long default_bits) {
    return mpfr_binop(a, b, default_bits, mpfr_add);
}
Expr* numeric_mpfr_sub(const Expr* a, const Expr* b, long default_bits) {
    return mpfr_binop(a, b, default_bits, mpfr_sub);
}
Expr* numeric_mpfr_mul(const Expr* a, const Expr* b, long default_bits) {
    return mpfr_binop(a, b, default_bits, mpfr_mul);
}
Expr* numeric_mpfr_div(const Expr* a, const Expr* b, long default_bits) {
    return mpfr_binop(a, b, default_bits, mpfr_div);
}
Expr* numeric_mpfr_pow(const Expr* a, const Expr* b, long default_bits) {
    return mpfr_binop(a, b, default_bits, mpfr_pow);
}

bool numeric_expr_is_mpfr(const Expr* e) {
    return expr_max_mpfr_prec(e) > 0;
}

Expr* numeric_mpfr_apply_unary(const Expr* e, long default_bits, MpfrUnaryOp op) {
    long bits = numeric_combined_bits(e, NULL, default_bits);
    mpfr_t r, i, out;
    mpfr_init2(r, bits);
    mpfr_init2(i, bits);
    bool ok = get_approx_mpfr(e, r, i, NULL);
    if (!ok || !mpfr_zero_p(i)) {
        mpfr_clear(r); mpfr_clear(i);
        return NULL;
    }
    mpfr_init2(out, bits);
    op(out, r, MPFR_RNDN);
    mpfr_clear(r); mpfr_clear(i);
    /* Domain failure (e.g. asin(2), acosh(0.5), log(-1)) — return NULL so
     * the caller can route through a complex-aware fallback rather than
     * propagating NaN into the symbolic result. */
    if (mpfr_nan_p(out)) {
        mpfr_clear(out);
        return NULL;
    }
    return expr_new_mpfr_move(out);
}

/* Complex power at MPFR precision via polar form.
 *
 * Writing base = a + b*I and exp = c + d*I, set r = |base| and theta = arg(base):
 *     base^exp = r^c * exp(-d*theta) * (cos(c*theta + d*ln(r))
 *                                       + I*sin(c*theta + d*ln(r)))
 *
 * This is the analytic continuation MPFR needs because there's no MPC linkage
 * in the build — without it, the cpow fallback in power.c silently coerces
 * MPFR operands to zero and produces NaN. See power.c (negative-base/complex-
 * exponent path) for the caller.
 *
 * Returns Complex[MPFR, MPFR], a pure MPFR when the imaginary part rounds to
 * zero, or NULL if either operand isn't an MPFR-promotable numeric value or
 * the base is zero (handled by the dedicated 0^x path elsewhere). */
Expr* numeric_mpfr_complex_pow(const Expr* base, const Expr* exp,
                               long default_bits) {
    long bits = numeric_combined_bits(base, exp, default_bits);
    if (bits <= 0) bits = 53;

    mpfr_t a, b, c, d;
    mpfr_init2(a, bits); mpfr_init2(b, bits);
    mpfr_init2(c, bits); mpfr_init2(d, bits);
    bool ok = get_approx_mpfr(base, a, b, NULL)
           && get_approx_mpfr(exp,  c, d, NULL);
    if (!ok) {
        mpfr_clear(a); mpfr_clear(b);
        mpfr_clear(c); mpfr_clear(d);
        return NULL;
    }

    mpfr_t r, theta;
    mpfr_init2(r, bits);
    mpfr_init2(theta, bits);
    mpfr_hypot(r, a, b, MPFR_RNDN);

    if (mpfr_zero_p(r)) {
        mpfr_clear(a); mpfr_clear(b); mpfr_clear(c); mpfr_clear(d);
        mpfr_clear(r); mpfr_clear(theta);
        return NULL;
    }

    mpfr_atan2(theta, b, a, MPFR_RNDN);  /* (-pi, pi] */

    mpfr_t ln_r, mag, angle, t1, t2;
    mpfr_init2(ln_r, bits); mpfr_init2(mag, bits); mpfr_init2(angle, bits);
    mpfr_init2(t1, bits);   mpfr_init2(t2, bits);

    mpfr_log(ln_r, r, MPFR_RNDN);

    /* mag = r^c * exp(-d*theta). For d == 0 the exp factor is exactly 1; we
     * still compute it generically — MPFR returns 1 with no precision loss. */
    mpfr_pow(t1, r, c, MPFR_RNDN);
    mpfr_mul(t2, d, theta, MPFR_RNDN);
    mpfr_neg(t2, t2, MPFR_RNDN);
    mpfr_exp(t2, t2, MPFR_RNDN);
    mpfr_mul(mag, t1, t2, MPFR_RNDN);

    /* angle = c*theta + d*ln(r). */
    mpfr_mul(t1, c, theta, MPFR_RNDN);
    mpfr_mul(t2, d, ln_r, MPFR_RNDN);
    mpfr_add(angle, t1, t2, MPFR_RNDN);

    Expr* re_expr = expr_new_mpfr_bits(bits);
    Expr* im_expr = expr_new_mpfr_bits(bits);
    mpfr_cos(t1, angle, MPFR_RNDN);
    mpfr_sin(t2, angle, MPFR_RNDN);
    mpfr_mul(re_expr->data.mpfr, mag, t1, MPFR_RNDN);
    mpfr_mul(im_expr->data.mpfr, mag, t2, MPFR_RNDN);

    bool im_is_zero = mpfr_zero_p(im_expr->data.mpfr);

    mpfr_clear(a); mpfr_clear(b); mpfr_clear(c); mpfr_clear(d);
    mpfr_clear(r); mpfr_clear(theta);
    mpfr_clear(ln_r); mpfr_clear(mag); mpfr_clear(angle);
    mpfr_clear(t1); mpfr_clear(t2);

    if (im_is_zero) {
        expr_free(im_expr);
        return re_expr;
    }
    return make_complex(re_expr, im_expr);
}

bool numeric_constant_mpfr_value(const char* name, mpfr_t out, long bits) {
    const NumericConstant* c = find_constant(name);
    if (!c || !c->mpfr_fill) return false;
    c->mpfr_fill(out, bits);
    return true;
}

/* ------------------------------------------------------------------------
 *  get_approx_mpfr — shared helper for MPFR-aware per-function builtins.
 *
 *  Models the existing `get_approx` pattern (see e.g. trig.c:158) but
 *  fills MPFR reals instead of a double complex. Precision is inherited
 *  from the caller's already-initialized `re` / `im`.
 * ---------------------------------------------------------------------- */
bool get_approx_mpfr(const Expr* e, mpfr_t re, mpfr_t im, bool* is_inexact) {
    if (!e) return false;
    bool inexact_local = false;
    if (!is_inexact) is_inexact = &inexact_local;

    switch (e->type) {
        case EXPR_INTEGER:
            mpfr_set_si(re, (long)e->data.integer, MPFR_RNDN);
            mpfr_set_zero(im, +1);
            return true;
        case EXPR_BIGINT:
            mpfr_set_z(re, e->data.bigint, MPFR_RNDN);
            mpfr_set_zero(im, +1);
            return true;
        case EXPR_REAL:
            mpfr_set_d(re, e->data.real, MPFR_RNDN);
            mpfr_set_zero(im, +1);
            *is_inexact = true;
            return true;
        case EXPR_MPFR:
            mpfr_set(re, e->data.mpfr, MPFR_RNDN);
            mpfr_set_zero(im, +1);
            *is_inexact = true;
            return true;
        case EXPR_FUNCTION: {
            int64_t n, d;
            if (is_rational((Expr*)e, &n, &d)) {
                mpfr_set_si(re, (long)n, MPFR_RNDN);
                mpfr_div_si(re, re, (long)d, MPFR_RNDN);
                mpfr_set_zero(im, +1);
                return true;
            }
            Expr *rpart, *ipart;
            if (is_complex((Expr*)e, &rpart, &ipart)) {
                /* Allocate scratch MPFRs at the requested precision. */
                mpfr_prec_t p = mpfr_get_prec(re);
                mpfr_t scratch_r, scratch_i;
                mpfr_init2(scratch_r, p);
                mpfr_init2(scratch_i, p);
                bool ok_r = get_approx_mpfr(rpart, re, scratch_i, is_inexact);
                if (!ok_r || mpfr_zero_p(scratch_i) == 0) {
                    /* A nested Complex in the real component is rare; we
                     * don't try to simplify it here and return failure
                     * rather than producing a wrong answer. */
                    mpfr_clear(scratch_r);
                    mpfr_clear(scratch_i);
                    return false;
                }
                bool ok_i = get_approx_mpfr(ipart, im, scratch_r, is_inexact);
                mpfr_clear(scratch_r);
                mpfr_clear(scratch_i);
                return ok_i;
            }
            return false;
        }
        default:
            return false;
    }
}
#endif

/* ------------------------------------------------------------------------
 *  Hold-form detection
 *
 *  N preserves Hold / HoldForm / Unevaluated wrappers — their job is to
 *  block evaluation, and that applies to N just as it does to ordinary
 *  rewriting.
 * ---------------------------------------------------------------------- */
static bool is_hold_head(const Expr* head) {
    if (!head || head->type != EXPR_SYMBOL) return false;
    const char* s = head->data.symbol.name;
    return s == SYM_HoldForm
        || s == SYM_Hold
        || s == SYM_HoldComplete
        || s == SYM_HoldPattern
        || s == SYM_Unevaluated;
}

/* ------------------------------------------------------------------------
 *  Leaf conversion
 *
 *  Turn a concrete numeric leaf into a `NumericMode`-appropriate value.
 *  Phase 1 always returns EXPR_REAL. Phase 2 will produce EXPR_MPFR when
 *  spec.mode == NUMERIC_MODE_MPFR.
 * ---------------------------------------------------------------------- */

static Expr* leaf_from_double(double v, NumericSpec spec) {
#ifdef USE_MPFR
    if (numeric_spec_is_mpfr(spec)) {
        return expr_new_mpfr_from_d(v, spec.bits);
    }
#else
    (void)spec;
#endif
    return expr_new_real(v);
}

static Expr* leaf_from_integer(int64_t v, NumericSpec spec) {
#ifdef USE_MPFR
    if (numeric_spec_is_mpfr(spec)) {
        /* long may be 32 bits on some platforms; use a path that carries
         * the full int64 range through the mpz_t bridge. */
        if (v >= (int64_t)LONG_MIN && v <= (int64_t)LONG_MAX) {
            return expr_new_mpfr_from_si((long)v, spec.bits);
        }
        mpz_t tmp;
        mpz_init_set_si(tmp, (long)(v >> 32));
        mpz_mul_2exp(tmp, tmp, 32);
        mpz_add_ui(tmp, tmp, (unsigned long)(v & 0xFFFFFFFFu));
        Expr* e = expr_new_mpfr_from_mpz(tmp, spec.bits);
        mpz_clear(tmp);
        return e;
    }
#endif
    return leaf_from_double((double)v, spec);
}

static Expr* leaf_from_bigint(const mpz_t v, NumericSpec spec) {
#ifdef USE_MPFR
    if (numeric_spec_is_mpfr(spec)) {
        return expr_new_mpfr_from_mpz(v, spec.bits);
    }
    /* Machine mode: a plain (double) conversion overflows to +/-inf once the
     * magnitude exceeds DBL_MAX (~1.8e308), e.g. N[1001!]. Mathematica's
     * machine-precision numbers carry an arbitrary exponent, so fall back to a
     * machine-precision (DBL_MANT_DIG-bit) MPFR real, which is finite for any
     * magnitude. We keep the IEEE double for in-range values so ordinary
     * machine arithmetic is unaffected. */
    {
        double d = mpz_get_d(v);
        if (isinf(d)) {
            return expr_new_mpfr_from_mpz(v, DBL_MANT_DIG);
        }
        return leaf_from_double(d, spec);
    }
#else
    return leaf_from_double(mpz_get_d(v), spec);
#endif
}

/* ------------------------------------------------------------------------
 *  Symbol numericalization
 *
 *  A recognized constant → its numeric value.
 *  Anything else (Infinity, True, user symbols) → a fresh copy.
 * ---------------------------------------------------------------------- */
static Expr* numericalize_symbol(const Expr* e, NumericSpec spec) {
    const NumericConstant* c = find_constant(e->data.symbol.name);
    if (c) {
#ifdef USE_MPFR
        if (numeric_spec_is_mpfr(spec) && c->mpfr_fill) {
            mpfr_t tmp;
            mpfr_init2(tmp, spec.bits);
            c->mpfr_fill(tmp, spec.bits);
            return expr_new_mpfr_move(tmp);
        }
#endif
        return leaf_from_double(c->machine_value, spec);
    }
    return expr_copy((Expr*)e);
}

/* ------------------------------------------------------------------------
 *  Function numericalization
 *
 *  - Rational[n, d] and Complex[re, im] are handled directly for speed and
 *    to keep Complex[N[a], N[b]] structurally intact.
 *  - Hold-like heads pass through unchanged.
 *  - Everything else: numericalize each argument, rebuild, re-evaluate.
 *
 *  The re-evaluation step is critical — it is how Sin[1.0] actually
 *  invokes `csin`, how Plus[2.0, 3.0] collapses to 5.0, and how any future
 *  MPFR/GSL-aware per-function code gets reached.
 * ---------------------------------------------------------------------- */
/* Forward declaration — the non-static definition follows below, since
 * `numericalize` is part of the public module API (declared in numeric.h). */

static Expr* numericalize_function(const Expr* e, NumericSpec spec) {
    /* Root[Function[p_in_slot1], k] → companion-matrix all-roots backend,
     * canonical sort, Newton refinement.  Runs first because Root is
     * HoldAll and we must NOT numericalize its arguments (Slot[1] and the
     * symbolic polynomial body would lose structure).  On NULL we fall
     * through to the generic rebuild, which preserves the original Root
     * call unchanged (the Root builtin itself returns NULL). */
    if (e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Root) {
        Expr* out = root_numericalize(e, spec);
        if (out) return out;
        /* Failed (e.g. non-integer coefficients, out-of-range k).  Don't
         * descend into the held body — return the Root call verbatim. */
        return expr_copy((Expr*)e);
    }

    /* Rational[n, d] → direct numeric quotient. Compute at full target
     * precision; a plain (double)n/d loses bits beyond 1e-15 and would
     * destroy the request for e.g. N[1/3, 40]. */
    int64_t rn, rd;
    if (is_rational((Expr*)e, &rn, &rd)) {
#ifdef USE_MPFR
        if (numeric_spec_is_mpfr(spec)) {
            Expr* r = expr_new_mpfr_from_si((long)rn, spec.bits);
            if (r) mpfr_div_si(r->data.mpfr, r->data.mpfr, (long)rd, MPFR_RNDN);
            return r;
        }
#endif
        return leaf_from_double((double)rn / (double)rd, spec);
    }

    /* Bigint-aware Rational[n, d] → direct numeric quotient. The int64 path
     * above only fires when both components fit int64; once a numerator or
     * denominator overflows into a BigInt, is_rational() reports false and the
     * generic rebuild below would numericalize each component independently and
     * leave a frozen Rational[Real, Real] (e.g. N[1/10^30] → 1.0/1.0e+30).
     * Compute the true quotient here instead, exactly like the int64 path. */
    if (e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Rational
        && e->data.function.arg_count == 2
        && expr_is_integer_like(e->data.function.args[0])
        && expr_is_integer_like(e->data.function.args[1])) {
        mpq_t q;
        mpq_init(q);
        expr_to_mpz(e->data.function.args[0], mpq_numref(q));
        expr_to_mpz(e->data.function.args[1], mpq_denref(q));
        mpq_canonicalize(q);
#ifdef USE_MPFR
        if (numeric_spec_is_mpfr(spec)) {
            mpfr_t r;
            mpfr_init2(r, spec.bits);
            mpfr_set_q(r, q, MPFR_RNDN);
            mpq_clear(q);
            return expr_new_mpfr_move(r);
        }
        /* Machine mode: a plain double quotient overflows to +/-inf or
         * underflows to 0 once either component exceeds DBL_MAX (~1.8e308) —
         * e.g. N[10^400/3] or N[1/10^400]. Mathematica's machine-precision
         * numbers carry an arbitrary exponent, so fall back to a
         * machine-precision (DBL_MANT_DIG-bit) MPFR, which is finite and
         * nonzero for any in-range rational. In-range values keep the IEEE
         * double so ordinary machine arithmetic is unaffected. */
        double dq = mpq_get_d(q);
        if (!isfinite(dq) || (dq == 0.0 && mpq_sgn(q) != 0)) {
            mpfr_t r;
            mpfr_init2(r, DBL_MANT_DIG);
            mpfr_set_q(r, q, MPFR_RNDN);
            mpq_clear(q);
            return expr_new_mpfr_move(r);
        }
        mpq_clear(q);
        return leaf_from_double(dq, spec);
#else
        double dq = mpq_get_d(q);
        mpq_clear(q);
        return leaf_from_double(dq, spec);
#endif
    }

    /* Complex[re, im] → numericalize components, rebuild via make_complex,
     * which normalizes (im == 0 → return the real component). */
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im)) {
        Expr* nre = numericalize(re, spec);
        Expr* nim = numericalize(im, spec);
        return make_complex(nre, nim);
    }

    /* Hold-forms: return a plain deep copy. */
    if (is_hold_head(e->data.function.head)) {
        return expr_copy((Expr*)e);
    }

    /* General function f[a, b, ...] — rebuild with numericalized args, then
     * feed the result to the evaluator to trigger any numeric fast paths in
     * the per-function builtins (e.g. trig.c's csin, power.c's cpow). */
    const size_t n = e->data.function.arg_count;
    Expr* new_head = numericalize(e->data.function.head, spec);

    /* Power[base, exp] needs special handling for two cases:
     *   1. Integer/bigint exponent: preserve verbatim. Matches Mathematica —
     *      N[x^2] is x^2, not x^2.0 — and keeps polynomial structure intact.
     *   2. Symbolic base + rational exponent: preserve the rational exponent
     *      so Sqrt[x] stays Sqrt[x] under contagion (1.0 * Sqrt[x] should not
     *      become x^0.5). Only when the base numericalizes to an actual number
     *      (e.g. Sqrt[2], Sqrt[Pi]) do we want the exponent demoted so the
     *      Power can collapse to a numeric value. */
    const bool is_power =
        (n == 2)
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Power;
    const bool exp_is_int =
        is_power
        && (e->data.function.args[1]->type == EXPR_INTEGER
            || e->data.function.args[1]->type == EXPR_BIGINT);

    Expr** new_args = (Expr**)malloc(n * sizeof(Expr*));
    if (!new_args) {
        /* Out-of-memory: return an uncomputed copy rather than crash. */
        expr_free(new_head);
        return expr_copy((Expr*)e);
    }
    /* Numericalize args left-to-right. For Power, decide whether to preserve
     * the exponent based on the (already-computed) numericalized base. */
    bool preserve_power_exp = exp_is_int;
    for (size_t i = 0; i < n; ++i) {
        if (is_power && i == 1) {
            /* By now new_args[0] holds the numericalized base. If the base
             * is still non-numeric, keep the exponent symbolic so e.g.
             * Power[x, 1/2] does not become Power[x, 0.5]. */
            if (!preserve_power_exp && !expr_is_numeric_like(new_args[0])) {
                preserve_power_exp = true;
            }
            if (preserve_power_exp) {
                new_args[i] = expr_copy(e->data.function.args[i]);
                continue;
            }
        }
        new_args[i] = numericalize(e->data.function.args[i], spec);
    }
    Expr* rebuilt = expr_new_function(new_head, new_args, n);
    free(new_args);  /* expr_new_function copies the pointer list. */

    /* eval_and_free consumes `rebuilt` and returns a fresh evaluated Expr. */
    return eval_and_free(rebuilt);
}

/* ------------------------------------------------------------------------
 *  Top-level dispatch
 * ---------------------------------------------------------------------- */
Expr* numericalize(const Expr* e, NumericSpec spec) {
    if (!e) return NULL;
    switch (e->type) {
        case EXPR_INTEGER: return leaf_from_integer(e->data.integer, spec);
        case EXPR_BIGINT:  return leaf_from_bigint(e->data.bigint, spec);
        case EXPR_REAL:
#ifdef USE_MPFR
            if (spec.mode == NUMERIC_MODE_MPFR) {
                /* Promote double to MPFR at the requested precision.
                 * The value is exact to 53 bits; bits beyond are zero.
                 * This matches Mathematica's "pad zero" promotion, which is
                 * what SetPrecision wants. N[x, p] uses NUMERIC_MODE_MPFR_CAP
                 * and deliberately falls through: a machine Real carries only
                 * ~15.95 digits and N must never manufacture precision it
                 * doesn't have, so it stays a machine Real (MachinePrecision). */
                return expr_new_mpfr_from_d(e->data.real, spec.bits);
            }
#endif
            return expr_new_real(e->data.real);
        case EXPR_STRING:  return expr_copy((Expr*)e);
        case EXPR_NDARRAY:  return expr_copy((Expr*)e); /* already machine-precision */
        case EXPR_SYMBOL:  return numericalize_symbol(e, spec);
        case EXPR_FUNCTION: return numericalize_function(e, spec);
#ifdef USE_MPFR
        case EXPR_MPFR:
            if (spec.mode == NUMERIC_MODE_MACHINE) {
                /* Bare N[expr] leaves already-approximate numbers at their
                 * existing precision — only exact quantities are numericalized
                 * to machine. Without this, N[N[Pi, 100]] (i.e. N[Pi,100]//N)
                 * would collapse 100 digits to machine precision. The flag is
                 * set solely by the one-argument N builtin; contagion and the
                 * explicit two-argument N[..., MachinePrecision] form leave it
                 * clear, so 1. + N[Pi,100] still collapses to machine. */
                if (spec.preserve_inexact) {
                    return expr_new_mpfr_copy(e->data.mpfr);
                }
                /* Down-convert to machine precision. A finite MPFR value can
                 * still exceed DBL_MAX (~1.8e308) — e.g. N[1.5 + 1001!], whose
                 * argument is already an MPFR ~4e2570. Plain mpfr_get_d would
                 * overflow to +/-inf, so keep such values as a machine-precision
                 * (DBL_MANT_DIG-bit) MPFR, which has an arbitrary exponent. */
                double d = mpfr_get_d(e->data.mpfr, MPFR_RNDN);
                if (isinf(d) && mpfr_number_p(e->data.mpfr)) {
                    Expr* r = expr_new_mpfr_bits(DBL_MANT_DIG);
                    if (r) mpfr_set(r->data.mpfr, e->data.mpfr, MPFR_RNDN);
                    return r;
                }
                return expr_new_real(d);
            }
            /* MPFR → MPFR. NUMERIC_MODE_MPFR sets the precision to exactly
             * spec.bits (SetPrecision pads up when the request exceeds the
             * value's current precision). NUMERIC_MODE_MPFR_CAP (the two-arg
             * N[x, p]) must never *increase* precision: cap the target at
             * min(existing, requested) so a 30-digit value stays 30 digits
             * under N[.., 50] but a 50-digit value still reduces under
             * N[.., 30]. */
            {
                long cur    = (long)mpfr_get_prec(e->data.mpfr);
                long target = spec.bits;
                if (spec.mode == NUMERIC_MODE_MPFR_CAP && cur < target) {
                    target = cur;
                }
                if (cur == target) {
                    return expr_new_mpfr_copy(e->data.mpfr);
                }
                Expr* r = expr_new_mpfr_bits(target);
                if (r) mpfr_set(r->data.mpfr, e->data.mpfr, MPFR_RNDN);
                return r;
            }
#endif
    }
    return expr_copy((Expr*)e);
}

/* ------------------------------------------------------------------------
 *  Inexact contagion
 *
 *  Mathematica's Plus / Times numericalize exact parts (Pi, E, Sqrt[2], …)
 *  whenever any other summand / factor is inexact. Without this rule the
 *  user sees `1. Pi` frozen as `1. Pi` instead of `3.14159`. The helper
 *  below is shared so Plus and Times stay in sync and future numeric
 *  heads can opt in with one call.
 * ---------------------------------------------------------------------- */

static bool arg_is_inexact(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
    if (expr_max_mpfr_prec(e) > 0) return true;
#endif
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im)) {
        if (re && re->type == EXPR_REAL) return true;
        if (im && im->type == EXPR_REAL) return true;
    }
    return false;
}

/* True iff the inexact content of `e` is machine precision — i.e. there
 * is an EXPR_REAL somewhere (possibly inside Complex[...]) and no MPFR
 * value alongside. MachinePrecision wins contagion: mixing `1.` (double)
 * with `1.0`50` (MPFR) must collapse to machine, not preserve 50 digits. */
static bool arg_has_machine_real(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im)) {
        if (re && re->type == EXPR_REAL) return true;
        if (im && im->type == EXPR_REAL) return true;
    }
    return false;
}

bool numeric_contagion_args(Expr* const* args, size_t n, Expr** out) {
    if (!args || n == 0) return false;

    bool any_inexact = false;
    bool any_machine = false;
    bool all_plain_real = true;   /* every arg already a bare machine double? */
#ifdef USE_MPFR
    long min_mpfr_prec = 0;  /* 0 = not yet observed */
#endif
    for (size_t i = 0; i < n; i++) {
        if (args[i]->type != EXPR_REAL) all_plain_real = false;
        if (arg_is_inexact(args[i])) {
            any_inexact = true;
            if (arg_has_machine_real(args[i])) any_machine = true;
#ifdef USE_MPFR
            long p = expr_max_mpfr_prec(args[i]);
            if (p > 0 && (min_mpfr_prec == 0 || p < min_mpfr_prec)) {
                min_mpfr_prec = p;
            }
#endif
        }
    }
    if (!any_inexact) return false;
    /* Fast path for the dominant case in tight float loops: every operand is
     * already a bare machine double, so contagion has nothing exact to
     * numericalize and would only clone each real into an identical one. Skip
     * the output array and the per-arg allocations entirely. */
    if (all_plain_real) return false;

    /* Precision contagion follows Mathematica: the *lowest* precision
     * among inexact operands wins. MachinePrecision is the floor — any
     * EXPR_REAL collapses the result to machine even alongside MPFR
     * values carrying more digits. With only MPFR operands, pick the
     * minimum MPFR precision so e.g. 1.0`50 + 1.0`20 lands at 20 digits
     * rather than preserving the 50-digit operand. */
    NumericSpec spec = numeric_machine_spec();
#ifdef USE_MPFR
    if (!any_machine && min_mpfr_prec > 0) {
        spec.mode = NUMERIC_MODE_MPFR;
        spec.bits = min_mpfr_prec;
    }
#endif

    for (size_t i = 0; i < n; i++) {
        out[i] = numericalize(args[i], spec);
    }
    return true;
}

/* ------------------------------------------------------------------------
 *  Precision argument parsing
 *
 *  N[expr]                   → machine precision
 *  N[expr, MachinePrecision] → machine precision
 *  N[expr, n] with Integer n → MPFR at n decimal digits (Phase 2)
 *
 *  Returns true on success and fills *out_spec. If the precision value is
 *  symbolic or non-numeric, returns false and the caller should keep the
 *  expression unevaluated (builtin returns NULL).
 * ---------------------------------------------------------------------- */
static bool parse_precision_arg(const Expr* prec, NumericSpec* out_spec) {
    /* MachinePrecision symbol → machine mode. */
    if (prec->type == EXPR_SYMBOL
        && prec->data.symbol.name == SYM_MachinePrecision) {
        *out_spec = numeric_machine_spec();
        return true;
    }

    /* Numeric prec → MPFR when available, otherwise warn + fall back. */
    int64_t rn, rd;
    double digits = 0.0;
    if (prec->type == EXPR_INTEGER) {
        digits = (double)prec->data.integer;
    } else if (prec->type == EXPR_REAL) {
        digits = prec->data.real;
    } else if (is_rational((Expr*)prec, &rn, &rd)) {
        digits = (double)rn / (double)rd;
    } else {
        return false;  /* symbolic / unknown — keep unevaluated */
    }

    if (digits <= 0.0) return false;

#ifdef USE_MPFR
    /* N[expr, p] caps inexact leaves at their existing precision — it never
     * manufactures digits. NUMERIC_MODE_MPFR_CAP encodes that; SetPrecision /
     * SetAccuracy (their own parsers) use plain NUMERIC_MODE_MPFR and pad up. */
    out_spec->mode = NUMERIC_MODE_MPFR_CAP;
    out_spec->bits = numeric_digits_to_bits(digits);
    out_spec->preserve_inexact = false;
    return true;
#else
    /* Phase 1 fallback: emit a one-shot warning, then use machine. */
    static bool warned = false;
    if (!warned) {
        fprintf(stderr,
                "N::prec: arbitrary precision unavailable in this build "
                "(USE_MPFR=0); using machine precision.\n");
        warned = true;
    }
    *out_spec = numeric_machine_spec();
    return true;
#endif
}

/* ------------------------------------------------------------------------
 *  Builtin: N
 * ---------------------------------------------------------------------- */
Expr* builtin_n(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    const size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return NULL;

    NumericSpec spec = numeric_machine_spec();
    if (argc == 2) {
        if (!parse_precision_arg(res->data.function.args[1], &spec)) {
            return NULL;  /* non-numeric precision → remain unevaluated */
        }
    } else {
        /* One-argument N[expr]: Mathematica's N numericalizes only the exact
         * parts of expr and preserves the precision of numbers that are
         * already approximate. So N[N[Pi, 100]] stays at 100 digits rather
         * than collapsing to machine precision (see numericalize's EXPR_MPFR
         * branch). The two-argument form is an explicit precision request and
         * deliberately does not set this. */
        spec.preserve_inexact = true;
    }

    /* Ownership: Mathilda evaluator frees `res` after a non-NULL return
     * (see eval.c around the builtin dispatch site), so we must NOT free
     * it ourselves. Returning NULL means "keep `res` unevaluated" — the
     * evaluator retains ownership in that path. */
    return numericalize(res->data.function.args[0], spec);
}

/* ------------------------------------------------------------------------
 *  Phase-2 MPFR fillers (stubs — implemented alongside the EXPR_MPFR work).
 * ---------------------------------------------------------------------- */
#ifdef USE_MPFR
static void fill_mpfr_pi(mpfr_t out, mpfr_prec_t bits) {
    mpfr_set_prec(out, bits);
    mpfr_const_pi(out, MPFR_RNDN);
}
static void fill_mpfr_e(mpfr_t out, mpfr_prec_t bits) {
    mpfr_set_prec(out, bits);
    mpfr_t one; mpfr_init2(one, bits); mpfr_set_ui(one, 1, MPFR_RNDN);
    mpfr_exp(out, one, MPFR_RNDN);
    mpfr_clear(one);
}
static void fill_mpfr_eulergamma(mpfr_t out, mpfr_prec_t bits) {
    mpfr_set_prec(out, bits);
    mpfr_const_euler(out, MPFR_RNDN);
}
static void fill_mpfr_catalan(mpfr_t out, mpfr_prec_t bits) {
    mpfr_set_prec(out, bits);
    mpfr_const_catalan(out, MPFR_RNDN);
}
static void fill_mpfr_goldenratio(mpfr_t out, mpfr_prec_t bits) {
    /* phi = (1 + sqrt(5)) / 2. Compute at bits + 20 guard then round. */
    mpfr_prec_t guard = bits + 20;
    mpfr_t tmp; mpfr_init2(tmp, guard);
    mpfr_set_ui(tmp, 5, MPFR_RNDN);
    mpfr_sqrt(tmp, tmp, MPFR_RNDN);
    mpfr_add_ui(tmp, tmp, 1, MPFR_RNDN);
    mpfr_div_ui(tmp, tmp, 2, MPFR_RNDN);
    mpfr_set_prec(out, bits);
    mpfr_set(out, tmp, MPFR_RNDN);
    mpfr_clear(tmp);
}
static void fill_mpfr_degree(mpfr_t out, mpfr_prec_t bits) {
    mpfr_prec_t guard = bits + 10;
    mpfr_t pi; mpfr_init2(pi, guard);
    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_set_prec(out, bits);
    mpfr_div_ui(out, pi, 180, MPFR_RNDN);
    mpfr_clear(pi);
}

/* GoldenAngle = (3 - sqrt(5)) pi = 2 pi / GoldenRatio^2, the closed form the
 * Mathematica constant carries. Computed at guard precision, then rounded. */
static void fill_mpfr_goldenangle(mpfr_t out, mpfr_prec_t bits) {
    mpfr_prec_t guard = bits + 20;
    mpfr_t s, pi; mpfr_init2(s, guard); mpfr_init2(pi, guard);
    mpfr_set_ui(s, 5, MPFR_RNDN);
    mpfr_sqrt(s, s, MPFR_RNDN);          /* sqrt(5)        */
    mpfr_ui_sub(s, 3, s, MPFR_RNDN);     /* 3 - sqrt(5)    */
    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_set_prec(out, bits);
    mpfr_mul(out, s, pi, MPFR_RNDN);     /* (3-sqrt(5)) pi */
    mpfr_clear(s); mpfr_clear(pi);
}

/* Glaisher-Kinkelin constant A, via
 *
 *     ln A = (gamma + ln(2 pi)) / 12  -  zeta'(2) / (2 pi^2).
 *
 * MPFR has no zeta-derivative primitive, so zeta'(2) = -sum_{n>=1} ln(n)/n^2
 * is evaluated by Euler-Maclaurin summation. With f(x) = ln(x) x^-2,
 *
 *     sum_{n>=1} f(n) = sum_{n=1}^{N-1} f(n) + int_N^inf f + f(N)/2
 *                       - sum_{k>=1} B_{2k}/(2k)! f^{(2k-1)}(N),
 *
 *     int_N^inf f = (ln N + 1)/N,
 *     f^{(2k-1)}(N) = (2k-1)! N^-(2k+1) [2k (H_{2k-1} - ln N) - (2k-1)].
 *
 * The required even Bernoulli numbers come from MPFR's zeta:
 *     B_{2k} = (-1)^(k+1) 2 (2k)! zeta(2k) / (2 pi)^(2k),
 * and (2k)! cancels against the 1/(2k)! prefactor, so the correction term is
 *     B_{2k}/(2k)! f^{(2k-1)}(N)
 *       = (-1)^(k+1) 2 zeta(2k) (2k-1)! / ((2 pi N)^(2k) N)
 *         * [2k (H_{2k-1} - ln N) - (2k-1)].
 * The asymptotic correction series is divergent, so we stop once a term stops
 * shrinking (smallest-term truncation). N ~ bits/6 makes that floor lie below
 * the target precision. */
static void fill_mpfr_glaisher(mpfr_t out, mpfr_prec_t bits) {
    mpfr_prec_t guard = bits + 64;
    unsigned long N = (unsigned long)(bits / 6);
    if (N < 40) N = 40;

    mpfr_t S, t, lnN, Hk, twopiN, num, prevmag, mag;
    mpfr_init2(S, guard);   mpfr_init2(t, guard);    mpfr_init2(lnN, guard);
    mpfr_init2(Hk, guard);  mpfr_init2(twopiN, guard);
    mpfr_init2(num, guard); mpfr_init2(prevmag, guard); mpfr_init2(mag, guard);

    /* head: sum_{n=2}^{N-1} ln(n)/n^2 (n=1 term is 0) */
    mpfr_set_zero(S, +1);
    for (unsigned long n = 2; n < N; ++n) {
        mpfr_set_ui(t, n, MPFR_RNDN);
        mpfr_log(t, t, MPFR_RNDN);          /* ln n          */
        mpfr_div_ui(t, t, n, MPFR_RNDN);
        mpfr_div_ui(t, t, n, MPFR_RNDN);    /* ln n / n^2    */
        mpfr_add(S, S, t, MPFR_RNDN);
    }
    /* ln N */
    mpfr_set_ui(lnN, N, MPFR_RNDN);
    mpfr_log(lnN, lnN, MPFR_RNDN);
    /* integral tail (ln N + 1)/N */
    mpfr_add_ui(t, lnN, 1, MPFR_RNDN);
    mpfr_div_ui(t, t, N, MPFR_RNDN);
    mpfr_add(S, S, t, MPFR_RNDN);
    /* f(N)/2 = ln N / (2 N^2) */
    mpfr_div_ui(t, lnN, N, MPFR_RNDN);
    mpfr_div_ui(t, t, N, MPFR_RNDN);
    mpfr_div_ui(t, t, 2, MPFR_RNDN);
    mpfr_add(S, S, t, MPFR_RNDN);

    /* (2 pi N) for the (2 pi N)^(2k) denominator */
    mpfr_const_pi(twopiN, MPFR_RNDN);
    mpfr_mul_ui(twopiN, twopiN, 2, MPFR_RNDN);
    mpfr_mul_ui(twopiN, twopiN, N, MPFR_RNDN);

    mpfr_set_inf(prevmag, +1);
    mpfr_set_zero(Hk, +1);                  /* running H_{2k-1} */
    for (unsigned long k = 1; k < 4 * N; ++k) {
        unsigned long two_k = 2 * k;
        /* H_{2k-1} = H_{2k-3} + 1/(2k-2) + 1/(2k-1) */
        if (k > 1) {
            mpfr_set_ui(t, 1, MPFR_RNDN);
            mpfr_div_ui(t, t, two_k - 2, MPFR_RNDN);
            mpfr_add(Hk, Hk, t, MPFR_RNDN);
        }
        mpfr_set_ui(t, 1, MPFR_RNDN);
        mpfr_div_ui(t, t, two_k - 1, MPFR_RNDN);
        mpfr_add(Hk, Hk, t, MPFR_RNDN);

        /* bracket = 2k (H_{2k-1} - ln N) - (2k-1) */
        mpfr_sub(num, Hk, lnN, MPFR_RNDN);
        mpfr_mul_ui(num, num, two_k, MPFR_RNDN);
        mpfr_sub_ui(num, num, two_k - 1, MPFR_RNDN);

        /* zeta(2k) */
        mpfr_zeta_ui(t, two_k, MPFR_RNDN);
        mpfr_mul(num, num, t, MPFR_RNDN);
        /* * 2 (2k-1)! */
        mpfr_mul_ui(num, num, 2, MPFR_RNDN);
        mpfr_set_ui(t, 1, MPFR_RNDN);
        for (unsigned long j = 1; j <= two_k - 1; ++j)
            mpfr_mul_ui(t, t, j, MPFR_RNDN);   /* (2k-1)! */
        mpfr_mul(num, num, t, MPFR_RNDN);
        /* / (2pi N)^(2k) */
        mpfr_pow_ui(t, twopiN, two_k, MPFR_RNDN);
        mpfr_div(num, num, t, MPFR_RNDN);
        /* / N */
        mpfr_div_ui(num, num, N, MPFR_RNDN);
        /* sign (-1)^(k+1) */
        if ((k % 2) == 0) mpfr_neg(num, num, MPFR_RNDN);

        /* smallest-term truncation: stop when the asymptotic term grows */
        mpfr_abs(mag, num, MPFR_RNDN);
        if (mpfr_cmp(mag, prevmag) > 0) break;
        mpfr_set(prevmag, mag, MPFR_RNDN);

        /* S -= B_{2k}/(2k)! f^{(2k-1)}(N) == num */
        mpfr_sub(S, S, num, MPFR_RNDN);
    }
    /* S now == sum_{n>=1} ln(n)/n^2 == -zeta'(2). So zeta'(2) = -S. */

    /* ln A = (gamma + ln(2 pi))/12 - zeta'(2)/(2 pi^2)
     *      = (gamma + ln(2 pi))/12 + S/(2 pi^2). */
    mpfr_t pi, lnA, term;
    mpfr_init2(pi, guard); mpfr_init2(lnA, guard); mpfr_init2(term, guard);
    mpfr_const_pi(pi, MPFR_RNDN);
    mpfr_const_euler(lnA, MPFR_RNDN);       /* gamma */
    mpfr_mul_ui(term, pi, 2, MPFR_RNDN);
    mpfr_log(term, term, MPFR_RNDN);        /* ln(2 pi) */
    mpfr_add(lnA, lnA, term, MPFR_RNDN);
    mpfr_div_ui(lnA, lnA, 12, MPFR_RNDN);   /* (gamma + ln(2 pi))/12 */
    mpfr_mul(term, pi, pi, MPFR_RNDN);
    mpfr_mul_ui(term, term, 2, MPFR_RNDN);  /* 2 pi^2 */
    mpfr_div(term, S, term, MPFR_RNDN);     /* S/(2 pi^2) */
    mpfr_add(lnA, lnA, term, MPFR_RNDN);

    mpfr_set_prec(out, bits);
    mpfr_exp(out, lnA, MPFR_RNDN);

    mpfr_clear(S); mpfr_clear(t); mpfr_clear(lnN); mpfr_clear(Hk);
    mpfr_clear(twopiN); mpfr_clear(num); mpfr_clear(prevmag); mpfr_clear(mag);
    mpfr_clear(pi); mpfr_clear(lnA); mpfr_clear(term);
}

/* Khinchin's constant K, via the geometrically convergent zeta series
 * (Bailey-Borwein-Crandall):
 *
 *     ln K * ln 2 = sum_{n>=1} (zeta(2n) - 1)/n * sum_{k=1}^{2n-1} (-1)^(k+1)/k.
 *
 * zeta(2n) - 1 ~ 4^-n, so ~bits/2 terms reach the target precision; we stop
 * once a term drops below the working epsilon. */
static void fill_mpfr_khinchin(mpfr_t out, mpfr_prec_t bits) {
    mpfr_prec_t guard = bits + 64;
    mpfr_t S, z, inner, h, term, mag, eps;
    mpfr_init2(S, guard);    mpfr_init2(z, guard);   mpfr_init2(inner, guard);
    mpfr_init2(h, guard);    mpfr_init2(term, guard);
    mpfr_init2(mag, guard);  mpfr_init2(eps, guard);

    /* eps = 2^-(bits+32): term cutoff */
    mpfr_set_ui(eps, 1, MPFR_RNDN);
    mpfr_div_2ui(eps, eps, bits + 32, MPFR_RNDN);

    mpfr_set_zero(S, +1);
    mpfr_set_zero(inner, +1);               /* running alternating harmonic */
    unsigned long last_k = 0;               /* last index summed into inner  */
    for (unsigned long n = 1; n < 8 * (unsigned long)bits + 16; ++n) {
        unsigned long upper = 2 * n - 1;
        /* extend inner sum_{k=1}^{2n-1} (-1)^(k+1)/k */
        for (unsigned long k = last_k + 1; k <= upper; ++k) {
            mpfr_set_ui(h, 1, MPFR_RNDN);
            mpfr_div_ui(h, h, k, MPFR_RNDN);
            if ((k % 2) == 0) mpfr_neg(h, h, MPFR_RNDN);
            mpfr_add(inner, inner, h, MPFR_RNDN);
        }
        last_k = upper;

        mpfr_zeta_ui(z, 2 * n, MPFR_RNDN);
        mpfr_sub_ui(z, z, 1, MPFR_RNDN);    /* zeta(2n) - 1 */
        mpfr_div_ui(z, z, n, MPFR_RNDN);    /* /n           */
        mpfr_mul(term, z, inner, MPFR_RNDN);
        mpfr_add(S, S, term, MPFR_RNDN);

        mpfr_abs(mag, term, MPFR_RNDN);
        if (n > 3 && mpfr_cmp(mag, eps) < 0) break;
    }

    /* K = exp(S / ln 2) */
    mpfr_const_log2(h, MPFR_RNDN);
    mpfr_div(S, S, h, MPFR_RNDN);
    mpfr_set_prec(out, bits);
    mpfr_exp(out, S, MPFR_RNDN);

    mpfr_clear(S); mpfr_clear(z); mpfr_clear(inner); mpfr_clear(h);
    mpfr_clear(term); mpfr_clear(mag); mpfr_clear(eps);
}
#endif /* USE_MPFR */

/* ------------------------------------------------------------------------
 *  Builtin: MachineNumberQ
 *
 *  MachineNumberQ[expr] returns True iff expr is a machine-precision real
 *  or complex number. The definition mirrors Mathematica:
 *
 *   - A finite EXPR_REAL is a machine number; +/-inf and NaN are not, even
 *     though they share the C `double` representation. This matches WL,
 *     where MachineNumberQ[Exp[1000.]] is False because Exp[1000.] overflows
 *     IEEE range (Mathilda's cexp returns +inf, hence the isfinite check).
 *   - EXPR_INTEGER, EXPR_BIGINT, Rational[_,_] are *exact* numbers — not
 *     machine numbers — so they return False.
 *   - EXPR_MPFR is arbitrary-precision (any precision other than machine),
 *     so it returns False even when the precision happens to be 53 bits.
 *   - Complex[re, im] is a machine number iff BOTH re and im are finite
 *     machine reals. Complex with one or both exact integer parts (e.g.
 *     Complex[1, 2], a Gaussian integer) is False — consistent with WL.
 *
 *  Any non-numeric input (symbol, head other than Complex, etc.) → False.
 *  Arity != 1 keeps the call unevaluated.
 * ---------------------------------------------------------------------- */

static bool is_machine_real_leaf(const Expr* e) {
    return e
        && e->type == EXPR_REAL
        && isfinite(e->data.real);
}

Expr* builtin_machinenumberq(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION
        || res->data.function.arg_count != 1) {
        return NULL;
    }
    Expr* arg = res->data.function.args[0];

    bool is_machine = false;
    if (is_machine_real_leaf(arg)) {
        is_machine = true;
    } else {
        Expr *re = NULL, *im = NULL;
        if (is_complex(arg, &re, &im)
            && is_machine_real_leaf(re)
            && is_machine_real_leaf(im)) {
            is_machine = true;
        }
    }
    return expr_new_symbol(is_machine ? "True" : "False");
}

/* ------------------------------------------------------------------------
 *  Registration
 * ---------------------------------------------------------------------- */
void numeric_init(void) {
    symtab_add_builtin("N", builtin_n);
    symtab_get_def("N")->attributes |= ATTR_PROTECTED | ATTR_LISTABLE;
    symtab_set_docstring("N",
        "N[expr]\n"
        "\tGives a machine-precision numerical approximation of expr.\n"
        "N[expr, n]\n"
        "\tGives a numerical approximation to n decimal digits. Requires\n"
        "\ta USE_MPFR build; without it, a warning is emitted and machine\n"
        "\tprecision is used.\n");

    symtab_add_builtin("MachineNumberQ", builtin_machinenumberq);
    symtab_get_def("MachineNumberQ")->attributes |= ATTR_PROTECTED;

    /* MachinePrecision is a reserved name — mark it protected so users
     * can't accidentally overwrite it. */
    symtab_get_def("MachinePrecision")->attributes |= ATTR_PROTECTED;

    /* Mathematical constants (Pi, E, EulerGamma, Catalan, GoldenRatio,
     * Degree) are first-class constant symbols. Stamp each Constant (so D
     * treats them as constants and they read as genuine constants) and
     * Protected (so they cannot be reassigned). Their numeric values live in
     * the kConstants table above, the single source of truth for the set;
     * EulerGamma's identity is additionally stamped in src/eulergamma.c. */
    for (size_t i = 0; i < kConstantCount; ++i) {
        SymbolDef* cdef = symtab_get_def(kConstants[i].name);
        if (cdef) cdef->attributes |= (ATTR_CONSTANT | ATTR_PROTECTED);
    }
}
