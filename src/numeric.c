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
#include "match.h"
#include "sym_names.h"
#include "root_numeric.h"
#include "root.h"        /* rootsum_expand_over_roots: N over a RootSum */
#include "sparse/sparse.h"  /* N over a SparseArray */
#include "ndarray.h"    /* ndt_get_i / ndarray_to_nested_list: N over a packed
                         * integer buffer widens it to float64 in place. */

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include "message.h"   /* mth_message: Quiet/Check funnel */
#include "print.h"     /* expr_to_string: N::meprec text */
#include "stieltjesgamma.h" /* numeric_stieltjes_n: N of an NHoldFirst constant */
#if defined(USE_FLINT) && defined(USE_MPFR)
#include <flint/arb.h>   /* arb twins of the MPFR elementary functions */
#include <flint/arf.h>
#endif

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

/* See numeric.h: a w-bit mantissa carries w - NUMERIC_PREC_GUARD_BITS bits of
 * precision, except the reserved machine width DBL_MANT_DIG (no guard) and the
 * width one above it, which stands in for the precision that would land on 53. */
long numeric_precision_bits_width(long prec_bits) {
    if (prec_bits < 1) prec_bits = 1;
    long w = prec_bits + NUMERIC_PREC_GUARD_BITS;
    if (w == (long)DBL_MANT_DIG) w++;
    return w;
}

long numeric_width_precision_bits(long width) {
    if (width == (long)DBL_MANT_DIG) return (long)DBL_MANT_DIG;
    if (width == (long)DBL_MANT_DIG + 1) return (long)DBL_MANT_DIG - NUMERIC_PREC_GUARD_BITS;
    return width > NUMERIC_PREC_GUARD_BITS ? width - NUMERIC_PREC_GUARD_BITS : 1;
}

long numeric_digits_to_bits(double digits) {
    if (!(digits > 0.0)) return numeric_precision_bits_width(1);
    if (digits > 1e15) digits = 1e15;
    return numeric_precision_bits_width((long)ceil(digits * LOG2_10));
}

double numeric_bits_to_digits(long bits) {
    if (bits <= 0) return 0.0;
    double p = (double)numeric_width_precision_bits(bits) / LOG2_10;
    if (bits != (long)DBL_MANT_DIG) {
        /* the width N[x, d] gives for a whole d reads back as exactly d */
        long d = (long)floor(p + 1e-9);
        if (d >= 1 && numeric_digits_to_bits((double)d) == bits) return (double)d;
    }
    return p;
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

/* An approximate zero carries only accuracy: the precision of its width IS
 * its accuracy in bits (Accuracy[0``a] = a; the width carries the guard bits
 * on top, numeric.h). When x + y
 * cancels exactly, the zero knows the sum only as well as the operands knew
 * their values, p - log2 max(|x|, |y|) bits below the units place -- not the
 * bare p bits it was computed with, which would make N[1/34, 30] -
 * N[1/34, 30] ~ 2^-100 uncertain where it is ~ 2^-105, and N[10^10, 30] -
 * N[10^10, 30] ~ 2^-100 where it is ~ 2^-67. Rounded down, so the accuracy
 * is never overstated. */
/* The precision in bits of a width, as a double for the accuracy formulas. */
static double width_prec(mpfr_prec_t w) {
    return (double)numeric_width_precision_bits((long)w);
}

/* Make `z` a zero of accuracy `acc` bits (stored at the width that carries
 * it), unless the accuracy is absurd. */
static void set_zero_accuracy(mpfr_t z, double acc) {
    if (acc < 1.0) acc = 1.0;
    if (acc > (double)MPFR_PREC_MAX / 4) return;
    mpfr_set_prec(z, (mpfr_prec_t)numeric_precision_bits_width((long)acc));
    mpfr_set_zero(z, +1);
}

static void cancelled_zero_accuracy(mpfr_t s, const mpfr_t x, const mpfr_t y) {
    if (!mpfr_zero_p(s) || (!mpfr_regular_p(x) && !mpfr_regular_p(y))) return;
    const mpfr_srcptr big = !mpfr_regular_p(y) ? x
                          : !mpfr_regular_p(x) ? y
                          : (mpfr_cmpabs(x, y) >= 0 ? x : y);
    long e;
    double m = mpfr_get_d_2exp(&e, big, MPFR_RNDN);
    set_zero_accuracy(s, floor(width_prec(mpfr_get_prec(s)) - ((double)e + log2(fabs(m)))));
}

/* The scaled accuracy of an approximate zero in a product or quotient. A
 * zero known to within 2^-a (a bits of accuracy, its bit precision) times a
 * nonzero y is known to within |y| 2^-a, so the product keeps a - log2|y|
 * bits; 0/y keeps a + log2|y|; 0 times 0 adds the two accuracies. The bits
 * the operation was computed with say nothing about a zero's uncertainty:
 * without this, 0``80 N[10^50, 100] stayed 0``80 (it is 0``30) and
 * 0``80 / 10^35 stayed 0``80 (it is 0``115). Rounded down, so never
 * overstated; `za`/`zb` are the operands as given (a zero's own bits). */
static void scaled_zero_accuracy(mpfr_t out, const Expr* za, const Expr* zb,
                                 const mpfr_t x, const mpfr_t y, bool is_div) {
    if (!mpfr_zero_p(out)) return;
    const bool x0 = za->type == EXPR_MPFR && mpfr_zero_p(za->data.mpfr);
    const bool y0 = zb->type == EXPR_MPFR && mpfr_zero_p(zb->data.mpfr);
    double acc;
    if (x0 && y0) {
        if (is_div) return;
        acc = width_prec(mpfr_get_prec(za->data.mpfr)) + width_prec(mpfr_get_prec(zb->data.mpfr));
    } else if (x0 && mpfr_regular_p(y)) {
        long e;
        double m = mpfr_get_d_2exp(&e, y, MPFR_RNDN);
        double l = (double)e + log2(fabs(m));
        acc = floor(width_prec(mpfr_get_prec(za->data.mpfr)) + (is_div ? l : -l));
    } else if (y0 && !is_div && mpfr_regular_p(x)) {
        long e;
        double m = mpfr_get_d_2exp(&e, x, MPFR_RNDN);
        acc = floor(width_prec(mpfr_get_prec(zb->data.mpfr)) - ((double)e + log2(fabs(m))));
    } else {
        return;
    }
    set_zero_accuracy(out, acc);
}

/* Generic MPFR binary dispatcher. Handles the common shape: extract real
 * MPFR operands from `a` and `b`, apply `op`, return a fresh EXPR_MPFR.
 * Returns NULL if either operand isn't a pure real (Complex → caller
 * handles separately). */
/* Worst cancellation seen by an MPFR addition / subtraction since the last
 * reset: the leading bits a sum lost (largest operand exponent minus the
 * result's; the whole precision for an exact zero from nonzero operands).
 * N[expr, p]'s adaptive precision re-evaluates only when this exceeded the
 * guard headroom -- the only way a guarded result can fall short of p digits
 * through arithmetic -- so an expression that does not cancel pays nothing. */
static long g_numeric_cancel_bits = 0;
/* The fewest bits any cancelling addition since the last reset kept: its
 * result's precision less the bits it lost, 0 for a cancellation to an exact
 * zero (which keeps nothing). N[expr, p] certifies a pass only when this is
 * at least p: an exact 0. from a cancellation agrees with every other one. */
static long g_numeric_cancel_left = LONG_MAX;
/* The deepest operand an MPFR addition since the last reset absorbed whole:
 * the exponent gap, in bits, between a sum's operands when the smaller lies
 * entirely below the result's precision (1/1000 - E^-1000 at 200 bits drops
 * an operand 1433 bits down).  When a pass then cancels to zero, the value
 * that vanished may be exactly such an absorbed part, and this is how much
 * precision it takes to see it -- N's next step uses it instead of only
 * doubling (see numericalize_top). */
static long g_numeric_absorb_bits = 0;

static void note_cancellation(const mpfr_t x, const mpfr_t y, const mpfr_t s) {
    if (!mpfr_regular_p(x) && !mpfr_regular_p(y)) return;
    if (mpfr_regular_p(x) && mpfr_regular_p(y)) {
        long ex = (long)mpfr_get_exp(x), ey = (long)mpfr_get_exp(y);
        long gap = ex > ey ? ex - ey : ey - ex;
        if (gap > (long)mpfr_get_prec(s) && gap > g_numeric_absorb_bits)
            g_numeric_absorb_bits = gap;
    }
    long top = LONG_MIN;
    if (mpfr_regular_p(x)) top = (long)mpfr_get_exp(x);
    if (mpfr_regular_p(y) && (long)mpfr_get_exp(y) > top) top = (long)mpfr_get_exp(y);
    long lost = mpfr_zero_p(s) ? (long)mpfr_get_prec(s)
              : mpfr_regular_p(s) ? top - (long)mpfr_get_exp(s) : 0;
    if (lost > g_numeric_cancel_bits) g_numeric_cancel_bits = lost;
    if (lost > 0) {
        long left = mpfr_zero_p(s) ? 0 : (long)mpfr_get_prec(s) - lost;
        if (left < g_numeric_cancel_left) g_numeric_cancel_left = left;
    }
}

void numeric_note_cancellation_bits(long lost, long left) {
    if (lost > g_numeric_cancel_bits) g_numeric_cancel_bits = lost;
    if (lost > 0) {
        if (left < 0) left = 0;
        if (left < g_numeric_cancel_left) g_numeric_cancel_left = left;
    }
}

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
    if (op == mpfr_add || op == mpfr_sub) {
        note_cancellation(ra, rb, out);
        cancelled_zero_accuracy(out, ra, rb);
    } else if (op == mpfr_mul || op == mpfr_div) {
        scaled_zero_accuracy(out, a, b, ra, rb, op == mpfr_div);
    }
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
/* x^y for real MPFR operands, the kernel behind real MPFR Power (interpreted
 * and compiled).  mpfr_pow is correctly rounded, and pays for it: before the
 * exp(y log x) it runs exactness and overflow probes, and its Ziv loop starts
 * from a precision it then often has to raise -- 2.6x the cost of the log and
 * exp it amounts to (102 us against 39 us for x^x^x at 397 bits).  For a
 * positive base and a non-integer exponent, exp(y log x) evaluated with
 * log2|y log x| + 24 guard bits is faithful to the target precision: the
 * relative error of t = y log x is ~2^-p' and becomes an absolute error |t|
 * 2^-p' in the exponent, i.e. the same relative error in exp(t).  Integer
 * exponents, non-positive bases, and the special values keep mpfr_pow (exact
 * cases, signs, and 0^y / Inf rules live there). */
int numeric_mpfr_pow_fast(mpfr_ptr rop, mpfr_srcptr x, mpfr_srcptr y,
                          mpfr_rnd_t rnd) {
    if (!mpfr_regular_p(x) || !mpfr_regular_p(y) || mpfr_sgn(x) <= 0
        || mpfr_integer_p(y) || mpfr_cmp_ui(x, 1) == 0)
        return mpfr_pow(rop, x, y, rnd);
    /* |log x| <= (|exp(x)| + 1) ln 2, so log2|y log x| <= exp(y) + log2(|exp(x)|+1) */
    long ex = (long)mpfr_get_exp(x);
    long lx = 1;
    for (long v = (ex < 0 ? -ex : ex) + 1; v > 1; v >>= 1) lx++;
    long ey = (long)mpfr_get_exp(y);
    long guard = 24 + (ey > 0 ? ey : 0) + lx;
    mpfr_prec_t wp = mpfr_get_prec(rop) + (mpfr_prec_t)guard;
    mpfr_t t;
    mpfr_init2(t, wp);
    mpfr_log(t, x, MPFR_RNDN);
    mpfr_mul(t, t, y, MPFR_RNDN);
    int r = mpfr_exp(rop, t, rnd);
    mpfr_clear(t);
    return r;
}

Expr* numeric_mpfr_pow(const Expr* a, const Expr* b, long default_bits) {
    return mpfr_binop(a, b, default_bits, numeric_mpfr_pow_fast);
}

bool numeric_expr_is_mpfr(const Expr* e) {
    return expr_max_mpfr_prec(e) > 0;
}

#ifdef USE_FLINT
/* ---- elementary functions through arb, correctly rounded -------------------
 *
 * MPFR's elementary functions are correctly rounded but slow: at 1000 bits
 * mpfr_atan takes 27 us and mpfr_asin 33 us, where FLINT's arb_atan / arb_asin
 * take 3.8 / 4.9 us; sin, tan, sinh, tanh are 2-3x faster in arb from 200 to
 * 6000 bits.  N of a composition of elementary functions deep in a
 * cancellation spends essentially all its time in these calls.  arb returns a
 * ball [m - r, m + r] that contains the true value, so the value is kept only
 * when both ends of the ball round (to nearest, at the target precision) to the
 * same number -- which is then the correctly rounded value, bit for bit what
 * the MPFR function returns.  Otherwise (a value within the ball's width of a
 * rounding boundary: probability ~2^-30; a domain boundary; a result outside
 * the MPFR exponent range) the MPFR function decides. */
typedef void (*ArbUnaryFn)(arb_t, const arb_t, slong);

static ArbUnaryFn arb_twin_of(MpfrUnaryOp op) {
    if (op == (MpfrUnaryOp)mpfr_sin)   return arb_sin;
    if (op == (MpfrUnaryOp)mpfr_cos)   return arb_cos;
    if (op == (MpfrUnaryOp)mpfr_tan)   return arb_tan;
    if (op == (MpfrUnaryOp)mpfr_sinh)  return arb_sinh;
    if (op == (MpfrUnaryOp)mpfr_cosh)  return arb_cosh;
    if (op == (MpfrUnaryOp)mpfr_tanh)  return arb_tanh;
    if (op == (MpfrUnaryOp)mpfr_atan)  return arb_atan;
    if (op == (MpfrUnaryOp)mpfr_asin)  return arb_asin;
    if (op == (MpfrUnaryOp)mpfr_acos)  return arb_acos;
    if (op == (MpfrUnaryOp)mpfr_atanh) return arb_atanh;
    if (op == (MpfrUnaryOp)mpfr_asinh) return arb_asinh;
    if (op == (MpfrUnaryOp)mpfr_acosh) return arb_acosh;
    if (op == (MpfrUnaryOp)mpfr_cot)   return arb_cot;
    if (op == (MpfrUnaryOp)mpfr_sec)   return arb_sec;
    if (op == (MpfrUnaryOp)mpfr_csc)   return arb_csc;
    if (op == (MpfrUnaryOp)mpfr_coth)  return arb_coth;
    if (op == (MpfrUnaryOp)mpfr_sech)  return arb_sech;
    if (op == (MpfrUnaryOp)mpfr_csch)  return arb_csch;
    if (op == (MpfrUnaryOp)mpfr_exp)   return arb_exp;
    if (op == (MpfrUnaryOp)mpfr_log)   return arb_log;
    return NULL;
}

/* True when the finite arf `x` rounds into MPFR's current exponent range. */
static bool arf_fits_mpfr(const arf_t x) {
    if (arf_is_zero(x)) return true;
    if (!arf_is_finite(x)) return false;
    return arf_cmpabs_2exp_si(x, (slong)mpfr_get_emax() - 2) < 0
        && arf_cmpabs_2exp_si(x, (slong)mpfr_get_emin() + 2) > 0;
}

/* out = f(x) correctly rounded to nearest at out's precision, or false. */
static bool arb_unary_round(mpfr_t out, const mpfr_t x, ArbUnaryFn f) {
    if (!mpfr_number_p(x)) return false;
    mpfr_prec_t p = mpfr_get_prec(out);
    arb_t a, b;
    arf_t lo, hi;
    arb_init(a); arb_init(b);
    arf_init(lo); arf_init(hi);
    arf_set_mpfr(arb_midref(a), x);
    mag_zero(arb_radref(a));
    f(b, a, (slong)p + 32);
    bool ok = false;
    if (arb_is_finite(b)) {
        arf_set_mag(hi, arb_radref(b));
        arf_sub(lo, arb_midref(b), hi, ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_add(hi, arb_midref(b), hi, ARF_PREC_EXACT, ARF_RND_DOWN);
        if (arf_fits_mpfr(lo) && arf_fits_mpfr(hi)) {
            mpfr_t rh;
            mpfr_init2(rh, p);
            arf_get_mpfr(out, lo, MPFR_RNDN);
            arf_get_mpfr(rh, hi, MPFR_RNDN);
            ok = mpfr_equal_p(out, rh) != 0;
            mpfr_clear(rh);
        }
    }
    arf_clear(lo); arf_clear(hi);
    arb_clear(a); arb_clear(b);
    return ok;
}
#endif

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
#ifdef USE_FLINT
    ArbUnaryFn twin = bits > 64 ? arb_twin_of(op) : NULL;
    if (!twin || !arb_unary_round(out, r, twin))
        op(out, r, MPFR_RNDN);
#else
    op(out, r, MPFR_RNDN);
#endif
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
/* True for an exact zero: Integer 0 (the exact part of Complex[0, y]). */
static bool is_exact_zero(const Expr* e) {
    return e && e->type == EXPR_INTEGER && e->data.integer == 0;
}

/* base^exp for a base on an axis -- a negative MPFR real, or Complex[0, y]
 * with an exact zero real part -- and an exact exponent n/2 that keeps the
 * result on an axis: base = |m| I^q exactly (q = 2 for a negative real, 1 or
 * 3 for +-I y), so base^(n/2) = |m|^(n/2) I^(q n/2) when q n is even (the
 * principal value: arg base = q Pi/2 lies in (-Pi, Pi]). The zero part is
 * exact, as it is in the base: the polar form computed it as cos(Pi/2)-size
 * noise (Sqrt[N[-2, 20]] had a real part 2.9*10^-21, Complex[0, 2.`20]^2 an
 * imaginary part 1.6*10^-20). NULL for any other shape. */
static Expr* mpfr_axial_pow(const Expr* base, const Expr* exp, long bits) {
    int64_t n2;                      /* the exponent times 2 */
    if (exp->type == EXPR_INTEGER) {
        if (exp->data.integer > INT64_MAX / 2 || exp->data.integer < INT64_MIN / 2) return NULL;
        n2 = 2 * exp->data.integer;
    } else {
        int64_t en, ed;
        if (!is_rational((Expr*)exp, &en, &ed) || ed != 2) return NULL;
        n2 = en;
    }
    const Expr* mag;                 /* the MPFR magnitude, signed */
    int q;
    Expr *re, *im;
    if (base->type == EXPR_MPFR && mpfr_sgn(base->data.mpfr) < 0) {
        mag = base; q = 2;
    } else if (is_complex((Expr*)base, &re, &im) && is_exact_zero(re)
               && im->type == EXPR_MPFR && mpfr_regular_p(im->data.mpfr)) {
        mag = im; q = mpfr_sgn(im->data.mpfr) > 0 ? 1 : 3;
    } else {
        return NULL;
    }
    if (((int64_t)q * n2) % 2 != 0) return NULL;      /* off the axes */
    /* k = (q n2 / 2) mod 4, from q n2 mod 8 (no overflow; q n2 is even) */
    int64_t k = (((int64_t)q * (n2 % 8)) % 8 + 8) % 8 / 2;
    long pb = (long)mpfr_get_prec(mag->data.mpfr);
    if (pb < bits) bits = pb;
    mpfr_t m, ex, v;
    mpfr_init2(m, bits); mpfr_init2(ex, 64); mpfr_init2(v, bits);
    mpfr_abs(m, mag->data.mpfr, MPFR_RNDN);
    mpfr_set_si(ex, (long)n2, MPFR_RNDN);
    mpfr_div_2ui(ex, ex, 1, MPFR_RNDN);              /* n2 / 2, exact */
    mpfr_pow(v, m, ex, MPFR_RNDN);
    mpfr_clear(m); mpfr_clear(ex);
    if (!mpfr_number_p(v)) { mpfr_clear(v); return NULL; }
    if (k == 2 || k == 3) mpfr_neg(v, v, MPFR_RNDN);
    Expr* val = expr_new_mpfr_move(v);
    if (k == 0 || k == 2) return val;
    return make_complex(expr_new_integer(0), val);
}

Expr* numeric_mpfr_imag_axis(const Expr* arg, const char* partner) {
    Expr *re, *im;
    if (!arg || !is_complex((Expr*)arg, &re, &im) || !is_exact_zero(re)
        || im->type != EXPR_MPFR || !mpfr_number_p(im->data.mpfr)) return NULL;
    Expr* g = eval_and_free(expr_new_function(expr_new_symbol(partner),
                                              (Expr*[]){ expr_copy(im) }, 1));
    if (!g) return NULL;
    if (g->type != EXPR_MPFR || !mpfr_number_p(g->data.mpfr)) { expr_free(g); return NULL; }
    return make_complex(expr_new_integer(0), g);
}

Expr* numeric_mpfr_complex_pow(const Expr* base, const Expr* exp,
                               long default_bits) {
    long bits = numeric_combined_bits(base, exp, default_bits);
    if (bits <= 0) bits = 53;
    {
        Expr* ax = mpfr_axial_pow(base, exp, bits);
        if (ax) return ax;
    }

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

#ifdef USE_MPFR
/* ------------------------------------------------------------------------
 *  Machine-precision rendering of an MPFR value
 *
 *  A machine number in Mathematica carries a 53-bit mantissa but an
 *  *arbitrary* exponent, so N[1001!] is 4.0279e2570 rather than Infinity.
 *  An IEEE double cannot do that, so whenever the value would overflow to
 *  +/-inf or flush a nonzero value to zero we keep a DBL_MANT_DIG-bit MPFR
 *  instead: same mantissa, unlimited exponent range.
 *
 *  Shared by the EXPR_MPFR down-convert, the bigint / rational leaf escapes
 *  and the working-precision round-back so all four agree by construction.
 * ---------------------------------------------------------------------- */
static Expr* mpfr_to_machine(const mpfr_t v) {
    double d = mpfr_get_d(v, MPFR_RNDN);
    if (mpfr_number_p(v)
        && ((isinf(d) && !mpfr_inf_p(v)) || (d == 0.0 && !mpfr_zero_p(v)))) {
        Expr* r = expr_new_mpfr_bits(DBL_MANT_DIG);
        if (r) {
            mpfr_set(r->data.mpfr, v, MPFR_RNDN);
            return r;
        }
    }
    return expr_new_real(d);
}
#endif

/* ------------------------------------------------------------------------
 *  Leaf conversion
 *
 *  Turn a concrete numeric leaf into a `NumericMode`-appropriate value.
 *  Phase 1 always returns EXPR_REAL. Phase 2 will produce EXPR_MPFR when
 *  spec.mode == NUMERIC_MODE_MPFR.
 * ---------------------------------------------------------------------- */

/* GMP's mpz_get_d / mpq_get_d truncate toward zero; IEEE conversion — and
 * therefore the decimal literal the user could have typed instead — rounds
 * to nearest. Left alone, N[10^25, MachinePrecision] and 1.0*^25 name two
 * different doubles, 2^31 apart. Route the inexact cases through MPFR,
 * which rounds; values that fit exactly need no help. Without MPFR there is
 * no rounding conversion available and GMP's truncation stands. */
static double mpz_to_double_rn(const mpz_t v) {
#ifdef USE_MPFR
    if (mpz_sizeinbase(v, 2) > (size_t)DBL_MANT_DIG) {
        mpfr_t t;
        mpfr_init2(t, DBL_MANT_DIG);
        mpfr_set_z(t, v, MPFR_RNDN);
        double d = mpfr_get_d(t, MPFR_RNDN);
        mpfr_clear(t);
        return d;
    }
#endif
    return mpz_get_d(v);
}

/* See numeric.h.  |q| = n/d is scaled by 2^-s so that the integer quotient
 * holds the double's significand: s = e - (DBL_MANT_DIG - 1) for
 * 2^e <= |q| < 2^(e+1), but never below the subnormal spacing exponent
 * DBL_MIN_EXP - DBL_MANT_DIG (-1074), so a subnormal keeps only the bits it
 * has.  The remainder decides the rounding (half to even), and
 * m 2^s is then exact in a double (m <= 2^53), so ldexp does not round. */
double numeric_mpq_to_double(const void* qv) {
    mpq_srcptr q = (mpq_srcptr)qv;
    int sgn = mpq_sgn(q);
    if (sgn == 0) return 0.0;
    mpz_t n, d, m, r;
    mpz_inits(n, d, m, r, NULL);
    mpz_abs(n, mpq_numref(q));
    mpz_set(d, mpq_denref(q));
    /* e = floor(log2(n/d)): n/d lies in [2^(e0-1), 2^(e0+1)) for
     * e0 = bits(n) - bits(d); one comparison picks the side. */
    long e = (long)mpz_sizeinbase(n, 2) - (long)mpz_sizeinbase(d, 2);
    mpz_set(m, d);
    if (e >= 0) mpz_mul_2exp(m, m, (mp_bitcnt_t)e);
    else        mpz_mul_2exp(r, n, (mp_bitcnt_t)(-e));
    if (e >= 0 ? (mpz_cmp(n, m) < 0) : (mpz_cmp(r, d) < 0)) e--;
    double out;
    if (e >= DBL_MAX_EXP) {
        out = HUGE_VAL;                               /* past DBL_MAX */
    } else {
        long s = e - (DBL_MANT_DIG - 1);
        const long smin = (long)DBL_MIN_EXP - DBL_MANT_DIG;   /* -1074 */
        if (s < smin) s = smin;
        /* m = round(n / (d 2^s)), half to even */
        if (s >= 0) { mpz_mul_2exp(d, d, (mp_bitcnt_t)s); }
        else        { mpz_mul_2exp(n, n, (mp_bitcnt_t)(-s)); }
        mpz_fdiv_qr(m, r, n, d);
        mpz_mul_2exp(r, r, 1);
        int c = mpz_cmp(r, d);
        if (c > 0 || (c == 0 && mpz_odd_p(m))) mpz_add_ui(m, m, 1);
        out = ldexp(mpz_get_d(m), (int)s);            /* exact, or inf */
    }
    mpz_clears(n, d, m, r, NULL);
    return sgn < 0 ? -out : out;
}

#ifdef USE_MPFR
static double mpq_to_double_rn(const mpq_t q) {
    return numeric_mpq_to_double(q);
}
#endif

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

static bool g_n_contagion;   /* defined below, with its contract */

static Expr* leaf_from_integer(int64_t v, NumericSpec spec) {
#ifdef USE_MPFR
    /* N[0, p] is the exact 0: a zero has no digits to approximate, so it
     * carries no precision, and an exact-zero part stays exact (N[{0, 1/3},
     * 20] is {0, 0.33...}, N[2 I, 20] is Complex[0, 2.`20]). Only N's own
     * capping mode: SetPrecision (plain MPFR mode) gives a zero a precision,
     * and machine N[0] is 0. as before. The exact operands of an
     * arbitrary-precision Plus/Times (inexact contagion, plain MPFR mode)
     * follow the same rule: 0 N[Pi, 20] is 0, Complex[0, 1] N[Pi, 20] is
     * Complex[0, 3.14...`20]. */
    if (v == 0 && (spec.mode == NUMERIC_MODE_MPFR_CAP
                   || (spec.mode == NUMERIC_MODE_MPFR && g_n_contagion)))
        return expr_new_integer(0);
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
        double d = mpz_to_double_rn(v);
        if (isinf(d)) {
            return expr_new_mpfr_from_mpz(v, DBL_MANT_DIG);
        }
        return leaf_from_double(d, spec);
    }
#else
    return leaf_from_double(mpz_to_double_rn(v), spec);
#endif
}

/* ------------------------------------------------------------------------
 *  Symbol numericalization
 *
 *  A recognized constant → its numeric value.
 *  Anything else (Infinity, True, user symbols) → a fresh copy.
 * ---------------------------------------------------------------------- */
static Expr* numericalize_rec(const Expr* e, NumericSpec spec);

/* Depth of N-value rewrites in progress (n_value below): a rule whose value
 * mentions its own tag again must not recurse without bound. */
static int g_n_value_depth = 0;

/* The value an N rule gives the symbol or call `e` whose tag is `tag`, at the
 * working precision -- N[e] at machine precision, N[e, digits] otherwise -- then
 * numericalized itself (z2 /: N[z2] = 3/2 gives 1.5); NULL when no N rule of
 * tag applies.  This is how N[2 + z2] reaches z2's rule inside a sum, as in the
 * reference.  A rule written N[z2] = ... has no precision argument, so like the
 * reference's N-value for {MachinePrecision, MachinePrecision} it serves only
 * machine N. */
static Expr* n_value(const Expr* e, const char* tag, NumericSpec spec) {
    if (symtab_n_value_count == 0 || g_n_value_depth > 64) return NULL;
    Rule* rules = symtab_get_n_values(tag);
    if (!rules) return NULL;
    /* Machine N asks for N[e, {MachinePrecision, MachinePrecision}], the key
     * an N[e] = v rule is stored under; arbitrary precision for N[e, digits]. */
    Expr* call;
    if (spec.mode == NUMERIC_MODE_MACHINE) {
        Expr* mp = expr_new_function(expr_new_symbol(SYM_List),
                       (Expr*[]){ expr_new_symbol(SYM_MachinePrecision),
                                  expr_new_symbol(SYM_MachinePrecision) }, 2);
        Expr* a[2] = { expr_copy((Expr*)e), mp };
        call = expr_new_function(expr_new_symbol(SYM_N), a, 2);
    } else {
        Expr* a[2] = { expr_copy((Expr*)e), expr_new_integer((int64_t)floor(numeric_bits_to_digits(spec.bits) + 0.5)) };  /* the precision, not the guarded width */
        call = expr_new_function(expr_new_symbol(SYM_N), a, 2);
    }
    Expr* v = NULL;
    for (Rule* r = rules; r && !v; r = r->next) {
        MatchEnv* env = env_new();
        if (match(call, r->pattern, env)) v = replace_bindings(r->replacement, env);
        env_free(env);
    }
    expr_free(call);
    if (!v) return NULL;
    g_n_value_depth++;
    Expr* ev = evaluate(v);
    expr_free(v);
    Expr* out = numericalize_rec(ev, spec);
    expr_free(ev);
    g_n_value_depth--;
    return out;
}

static Expr* numericalize_symbol(const Expr* e, NumericSpec spec) {
    {
        Expr* nv = n_value(e, e->data.symbol.name, spec);
        if (nv) return nv;
    }
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
/* Forward declaration — the recursion step is defined below. The public
 * `numericalize` (declared in numeric.h) is a thin wrapper around it that
 * picks the working precision and rounds the result back. */
static Expr* numericalize_rec(const Expr* e, NumericSpec spec);

/* > 0 while numericalizing inside an argument the enclosing head HOLDS (If's
 * branches, a Function body): the structure is numericalized but never
 * evaluated, so N[If[c, Print[a], b]] is If[c., Print[a], b.] and nothing
 * prints.  A counter, not a NumericSpec field: callers across the tree build
 * specs field by field, and a new field would be stack garbage there.  The
 * public numericalize() saves, zeroes and restores it, so a nested N and any
 * state an unwind left behind never leak into an unrelated call. */
static int g_n_held_depth = 0;

/* Inexact-contagion intent, set only while numeric_contagion_args
 * numericalizes the exact operands of a Plus / Times that has an inexact one:
 * then numericalize must NOT thread into the arguments of a non-numeric head
 * (`1.0 x[1]` stays `1.0 x[1]`, never `1.0 x[1.0]`); numeric heads (Plus,
 * Times, Power, Sin, ...) still recurse, so `1.0 + 2 x` -> `1. + 2. x`. Every
 * other entry -- explicit N[], and each internal caller -- runs with it off
 * (N[x[1]] == x[1.]), because the public numericalize() clears it. It used to
 * be a NumericSpec field, which callers that build their spec by hand left
 * uninitialised. */
static bool g_n_contagion = false;

/* ---- common subexpressions within one numericalization pass --------------
 *
 * N of a composition substitutes whole expressions for a variable, so the
 * same exact subtree recurs: f[g[x]] with a twelve-term g puts g's tree in
 * every one of f's twelve leaf positions, and each pass of N's adaptive
 * precision (at up to thousands of bits) re-evaluated every copy -- 468
 * elementary-function calls at the precision where 72 distinct ones are the
 * whole job.  A pass therefore memoizes the value of each call of a
 * NumericFunction head (Plus, Times, Power, Sin, ArcTan, ...: pure functions
 * of their numeric arguments, so evaluating a second identical copy can only
 * give the same value again) keyed by the subtree's structure, the spec and
 * the contagion mode.  The memo lives for one pass of numericalize_top --
 * a nested numericalize_top, which reads the pass's cancellation statistics,
 * gets its own -- so the statistics a hit skips were recorded by the first
 * evaluation in the same pass. */
typedef struct {
    Expr* key;          /* owned copy (a refcount) of the exact subtree */
    Expr* val;          /* owned numericalized value */
    uint64_t h;
    long bits;
    int mode;
    bool contagion;
} NMemoEntry;

typedef struct {
    NMemoEntry* tab;    /* open addressing, cap a power of two */
    size_t cap, count;
} NMemo;

static NMemo* g_nmemo = NULL;

static bool nmemo_eligible(const Expr* e) {
    if (!g_nmemo || g_n_held_depth > 0 || e->data.function.arg_count == 0) return false;
    const Expr* h = e->data.function.head;
    return h && h->type == EXPR_SYMBOL
        && (get_attributes(h->data.symbol.name) & ATTR_NUMERICFUNCTION);
}

static Expr* nmemo_find(const Expr* e, uint64_t h, NumericSpec spec) {
    NMemo* m = g_nmemo;
    if (!m->cap) return NULL;
    for (size_t i = (size_t)h & (m->cap - 1);; i = (i + 1) & (m->cap - 1)) {
        NMemoEntry* t = &m->tab[i];
        if (!t->key) return NULL;
        if (t->h == h && t->bits == spec.bits && t->mode == (int)spec.mode
            && t->contagion == g_n_contagion && expr_eq(t->key, (Expr*)e))
            return expr_copy(t->val);
    }
}

static void nmemo_put(const Expr* e, uint64_t h, NumericSpec spec, const Expr* val) {
    NMemo* m = g_nmemo;
    if (!val) return;
    if (2 * (m->count + 1) > m->cap) {
        size_t nc = m->cap ? 2 * m->cap : 64;
        NMemoEntry* nt = (NMemoEntry*)calloc(nc, sizeof(NMemoEntry));
        if (!nt) return;
        for (size_t i = 0; i < m->cap; i++) {
            if (!m->tab[i].key) continue;
            size_t j = (size_t)m->tab[i].h & (nc - 1);
            while (nt[j].key) j = (j + 1) & (nc - 1);
            nt[j] = m->tab[i];
        }
        free(m->tab);
        m->tab = nt;
        m->cap = nc;
    }
    size_t j = (size_t)h & (m->cap - 1);
    while (m->tab[j].key) j = (j + 1) & (m->cap - 1);
    NMemoEntry* t = &m->tab[j];
    t->key = expr_copy((Expr*)e);
    t->val = expr_copy((Expr*)val);
    t->h = h;
    t->bits = spec.bits;
    t->mode = (int)spec.mode;
    t->contagion = g_n_contagion;
    m->count++;
}

static void nmemo_clear(NMemo* m) {
    for (size_t i = 0; i < m->cap; i++)
        if (m->tab[i].key) { expr_free(m->tab[i].key); expr_free(m->tab[i].val); }
    free(m->tab);
    m->tab = NULL;
    m->cap = m->count = 0;
}

static Expr* numericalize_function_uncached(const Expr* e, NumericSpec spec);
static Expr* numericalize_pass(const Expr* e, NumericSpec spec);

static Expr* numericalize_function(const Expr* e, NumericSpec spec) {
    if (!nmemo_eligible(e)) return numericalize_function_uncached(e, spec);
    uint64_t h = expr_hash(e);
    Expr* hit = nmemo_find(e, h, spec);
    if (hit) return hit;
    Expr* r = numericalize_function_uncached(e, spec);
    nmemo_put(e, h, spec, r);
    return r;
}

static Expr* numericalize_function_uncached(const Expr* e, NumericSpec spec) {
    /* An N rule on the head (h /: N[h[x_]] := ...) gives the whole call's value. */
    if (e->data.function.head && e->data.function.head->type == EXPR_SYMBOL) {
        Expr* nv = n_value(e, e->data.function.head->data.symbol.name, spec);
        if (nv) return nv;
    }
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

    /* RootSum[Function[p], Function[body]] → numericalise the explicit sum of
     * body over the roots of p.  RootSum is HoldAll and often has no closed
     * rational form (a Log logand, an (x + t) denominator), so N sums it
     * directly: expand to Plus[body[Root[p, k]], ...] and numericalise that,
     * reusing the Root numeric backend.  (docs/notes/MATHILDA_DIVERGENCES.md A12.) */
    if (e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_RootSum) {
        Expr* expanded = rootsum_expand_over_roots(e);
        if (expanded) {
            Expr* out = numericalize_rec(expanded, spec);
            expr_free(expanded);
            return out;
        }
    }

    /* AlgebraicNumber[theta, {c0..cm}] → numericalize sum ci theta^i.  Runs
     * before the generic rebuild because AlgebraicNumber is NHoldAll and
     * numericalize does NOT honour NHold: the generic path would numericalize
     * the generator and coefficient List in place and freeze the object (the
     * AlgebraicNumber builtin declines on a machine-Real generator). Building
     * the power sum and recursing reuses the Root/radical numeric backend, so
     * this carries N to any requested precision. */
    if (e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_AlgebraicNumber
        && e->data.function.arg_count == 2
        && e->data.function.args[1]->type == EXPR_FUNCTION
        && e->data.function.args[1]->data.function.head->type == EXPR_SYMBOL
        && e->data.function.args[1]->data.function.head->data.symbol.name == SYM_List) {
        const Expr* theta = e->data.function.args[0];
        const Expr* cl    = e->data.function.args[1];
        size_t m = cl->data.function.arg_count;
        Expr** terms = (Expr**)malloc((m ? m : 1) * sizeof(Expr*));
        size_t nt = 0;
        for (size_t i = 0; i < m; i++) {
            Expr* ci = expr_copy(cl->data.function.args[i]);
            Expr* term;
            if (i == 0) {
                term = ci;                       /* constant term */
            } else {
                Expr* pw = (i == 1)
                    ? expr_copy((Expr*)theta)
                    : expr_new_function(expr_new_symbol(SYM_Power),
                          (Expr*[]){ expr_copy((Expr*)theta),
                                     expr_new_integer((int64_t)i) }, 2);
                term = expr_new_function(expr_new_symbol(SYM_Times),
                          (Expr*[]){ ci, pw }, 2);
            }
            terms[nt++] = term;
        }
        Expr* poly = (nt == 0) ? expr_new_integer(0)
                   : (nt == 1) ? terms[0]
                   : expr_new_function(expr_new_symbol(SYM_Plus), terms, nt);
        free(terms);
        Expr* out = numericalize_rec(poly, spec);
        expr_free(poly);
        return out;
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
        /* One IEEE division is one correct rounding when both parts are
         * exact doubles (|.| <= 2^53); otherwise converting them first would
         * round twice. */
        const int64_t exact_max = INT64_C(1) << DBL_MANT_DIG;
        if (rn <= exact_max && rn >= -exact_max && rd <= exact_max && rd >= -exact_max)
            return leaf_from_double((double)rn / (double)rd, spec);
        mpq_t q;
        mpq_init(q);
        mpz_set_si(mpq_numref(q), 0);
        {
            /* int64 -> mpz without assuming long is 64 bits (LLP64) */
            mpz_t hi;
            mpz_init_set_si(hi, (long)(rn >> 32));
            mpz_mul_2exp(hi, hi, 32);
            mpz_add_ui(hi, hi, (unsigned long)(rn & 0xFFFFFFFFu));
            mpz_set(mpq_numref(q), hi);
            mpz_set_si(hi, (long)(rd >> 32));
            mpz_mul_2exp(hi, hi, 32);
            mpz_add_ui(hi, hi, (unsigned long)(rd & 0xFFFFFFFFu));
            mpz_set(mpq_denref(q), hi);
            mpz_clear(hi);
        }
        mpq_canonicalize(q);
        double dq = numeric_mpq_to_double(q);
        mpq_clear(q);
        return leaf_from_double(dq, spec);
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
        /* No mpq_init: expr_to_mpz initialises its target itself, so
         * pre-initialising would allocate limbs and immediately leak them. */
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
        double dq = mpq_to_double_rn(q);
        if (!isfinite(dq) || (dq == 0.0 && mpq_sgn(q) != 0)) {
            mpfr_t r;
            mpfr_init2(r, DBL_MANT_DIG);
            mpfr_set_q(r, q, MPFR_RNDN);
            mpq_clear(q);
            Expr* out = mpfr_to_machine(r);
            mpfr_clear(r);
            return out;
        }
        mpq_clear(q);
        return leaf_from_double(dq, spec);
#else
        double dq = numeric_mpq_to_double(q);
        mpq_clear(q);
        return leaf_from_double(dq, spec);
#endif
    }

    /* Complex[re, im] → numericalize components, rebuild via make_complex,
     * which normalizes (im == 0 → return the real component). */
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im)) {
        Expr* nre = numericalize_rec(re, spec);
        Expr* nim = numericalize_rec(im, spec);
        return make_complex(nre, nim);
    }

    /* Hold-forms: return a plain deep copy. */
    if (is_hold_head(e->data.function.head)) {
        return expr_copy((Expr*)e);
    }

    /* Inexact contagion must not thread into the arguments of a non-numeric
     * head: `1.0 x[1]` stays `1.0 x[1]`, never `1.0 x[1.0]` (which would leave
     * an unbound x[1.0] and silently break, e.g., NMinimize's indexed-variable
     * substitution). A numeric head (Plus, Times, Power, Sin, ...) still
     * recurses below, so `1.0 + 2 x` → `1. + 2. x` and `1.0 + Sin[2]` numericalizes.
     * Explicit N[] runs with g_n_contagion off and keeps threading everywhere,
     * matching Mathematica's N[x[1]] == x[1.]. */
    if (g_n_contagion) {
        const Expr* h = e->data.function.head;
        bool numeric_head = h && h->type == EXPR_SYMBOL
            && (get_attributes(h->data.symbol.name) & ATTR_NUMERICFUNCTION);
        if (!numeric_head) return expr_copy((Expr*)e);
    }

    /* StieltjesGamma[n] has nothing N could numericalize (its order is
     * NHoldFirst), so its value is asked for directly. */
    if (e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_StieltjesGamma) {
        Expr* r = numeric_stieltjes_n(e, spec);
        if (r) return r;
    }

#ifdef USE_MPFR
    /* E^y at arbitrary precision is exp(y), one correctly rounded call.  The
     * generic rebuild numericalized E itself (an exp(1)) and then evaluated
     * Power[e, y] as exp(y log e) -- three transcendental evaluations for
     * every Exp[..] N meets, at every precision of its adaptive climb.  An
     * approximate machine exponent, a complex one, or an exp outside the
     * MPFR range takes the generic path. */
    if (numeric_spec_is_mpfr(spec) && g_n_held_depth == 0
        && e->data.function.arg_count == 2
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Power
        && e->data.function.args[0]->type == EXPR_SYMBOL
        && e->data.function.args[0]->data.symbol.name == SYM_E) {
        Expr* y = numericalize_rec(e->data.function.args[1], spec);
        if (y && y->type == EXPR_MPFR && mpfr_number_p(y->data.mpfr)) {
            Expr* r = numeric_mpfr_apply_unary(y, 0, mpfr_exp);
            if (r && r->type == EXPR_MPFR && mpfr_number_p(r->data.mpfr)) {
                expr_free(y);
                return r;
            }
            if (r) expr_free(r);
        }
        if (y) expr_free(y);
    }
#endif

    /* General function f[a, b, ...] — rebuild with numericalized args, then
     * feed the result to the evaluator to trigger any numeric fast paths in
     * the per-function builtins (e.g. trig.c's csin, power.c's cpow). */
    const size_t n = e->data.function.arg_count;
    Expr* new_head = numericalize_rec(e->data.function.head, spec);

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
    /* Integer-label arguments that N must leave exact. DirichletL[k, j, s]
     * labels a character by (k, j); Mathematica keeps them exact (NHoldAll
     * plus an N rule for s). Numericalizing them would give DirichletL[3., 2.,
     * 0.5`50], which is not a character label and stays unevaluated. */
    size_t n_label_held = 0;
    bool dirichlet_l = n == 3 && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_DirichletL;
    if (dirichlet_l)
        n_label_held = 2;
    /* NHoldFirst / NHoldRest: positions N leaves exact (StieltjesGamma's order,
     * Part's indices). DirichletL is NHoldAll yet numericalizes s, so its
     * labels are the explicit rule above. */
    uint32_t nhold = 0;
    if (!dirichlet_l && e->data.function.head->type == EXPR_SYMBOL)
        nhold = get_attributes(e->data.function.head->data.symbol.name) &
                (ATTR_NHOLDFIRST | ATTR_NHOLDREST);
    /* Evaluation-hold positions (If's branches, Function / Module bodies):
     * numericalized, but under no_eval, so a held Print never runs. */
    uint32_t ehold = 0;
    if (e->data.function.head->type == EXPR_SYMBOL)
        ehold = get_attributes(e->data.function.head->data.symbol.name) &
                (ATTR_HOLDFIRST | ATTR_HOLDREST | ATTR_HOLDALLCOMPLETE);

    /* Numericalize args left-to-right. For Power, decide whether to preserve
     * the exponent based on the (already-computed) numericalized base. */
    bool preserve_power_exp = exp_is_int;
    for (size_t i = 0; i < n; ++i) {
        if (i < n_label_held || (i == 0 && (nhold & ATTR_NHOLDFIRST)) ||
            (i > 0 && (nhold & ATTR_NHOLDREST))) {
            new_args[i] = expr_copy(e->data.function.args[i]);
            continue;
        }
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
        bool held = (ehold & ATTR_HOLDALLCOMPLETE) || (i == 0 && (ehold & ATTR_HOLDFIRST)) ||
                    (i > 0 && (ehold & ATTR_HOLDREST));
        if (held) g_n_held_depth++;
        new_args[i] = numericalize_rec(e->data.function.args[i], spec);
        if (held) g_n_held_depth--;
    }
    Expr* rebuilt = expr_new_function(new_head, new_args, n);
    free(new_args);  /* expr_new_function copies the pointer list. */
    if (g_n_held_depth > 0) return rebuilt;

    /* eval_and_free consumes `rebuilt` and returns a fresh evaluated Expr. */
    return eval_and_free(rebuilt);
}

/* ------------------------------------------------------------------------
 *  Recursive dispatch
 *
 *  This is the workhorse. It runs at the *working* precision chosen by
 *  numericalize() below, which is not necessarily the precision the caller
 *  asked for — see the "Working precision" section further down.
 * ---------------------------------------------------------------------- */
/* N over a SparseArray: every explicit value and the background, sharing the
 * pattern. At machine precision a typed value buffer converts directly. */
static Expr* numericalize_sparse_cb(const Expr* v, void* ctx) {
    return numericalize_rec(v, *(const NumericSpec*)ctx);
}

static Expr* numericalize_sparse(const Expr* e, NumericSpec spec) {
    /* Inexact contagion never reaches into a raw object. */
    if (g_n_contagion) return expr_copy((Expr*)e);
    Expr* r = sp_map_values(e, numericalize_sparse_cb, &spec, spec.mode == NUMERIC_MODE_MACHINE);
    return r ? r : expr_copy((Expr*)e);
}

static Expr* numericalize_rec(const Expr* e, NumericSpec spec) {
    if (!e) return NULL;
    switch (e->type) {
        case EXPR_SPARSE: return numericalize_sparse(e, spec);
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
        case EXPR_NDARRAY:
            /* A machine (real/complex) buffer is already numeric and passes
             * through. An INTEGER buffer still has to be numericalized: N of a
             * packed integer array must give a packed float64 array -- the same
             * machine reals the interpreter's N gives element by element -- not
             * the int64 buffer copied verbatim, and not an unpacked list of
             * boxed reals. Widen the whole buffer in one pass, inheriting the
             * source's presentation (a packed List stays a packed List; a
             * visible NDArray stays visible). This is what lets N claim
             * packed_int64_ok (src/pack.c): the transparency gate no longer
             * materialises N's integer argument, so control reaches here rather
             * than the element-by-element rebuild in numericalize_function. */
            if (e->data.ndarray.dtype == NDT_INT64) {
                if (spec.mode == NUMERIC_MODE_MACHINE) {
                    const int rank = e->data.ndarray.rank;
                    const int64_t* dims = e->data.ndarray.dims;
                    size_t sz = 1;
                    for (int i = 0; i < rank; i++) sz *= (size_t)dims[i];
                    double* out = (double*)malloc(sz * sizeof(double));
                    if (out) {
                        const void* in = e->data.ndarray.data;
                        /* Exact for |v| < 2^53; larger int64 values round to the
                         * nearest double -- which is exactly what N does to the
                         * same integer as a scalar, so the two surfaces agree. */
                        for (size_t i = 0; i < sz; i++)
                            out[i] = (double)ndt_get_i(in, i, NDT_INT64);
                        return expr_new_ndarray_like(e, rank, dims, out,
                                                     NDT_FLOAT64);
                    }
                    /* Out of memory: fall through to the verbatim copy. */
                } else {
                    /* Precision request N[arr, p]: MPFR is not a machine dtype,
                     * so there is no packed answer. Materialise and numericalize
                     * element by element -- exactly the path the gate used to
                     * force by materialising an int64 argument for a head
                     * without the int64 claim, so the answer is unchanged. */
                    Expr* nested = ndarray_to_nested_list(e);
                    Expr* r = numericalize_rec(nested, spec);
                    expr_free(nested);
                    return r;
                }
            }
            return expr_copy((Expr*)e); /* already machine-precision */
        case EXPR_COMPILED: return expr_copy((Expr*)e); /* opaque object; N[] is a no-op */
        case EXPR_SYMBOL:  return numericalize_symbol(e, spec);
        case EXPR_FUNCTION: return numericalize_function(e, spec);
#ifdef USE_MPFR
        case EXPR_MPFR:
            if (spec.mode == NUMERIC_MODE_MACHINE) {
                /* Down-convert an already-approximate arbitrary-precision value
                 * to a machine double. Mathematica's one-argument N[expr]
                 * targets machine precision even when expr is already inexact,
                 * so N[N[Pi, 100]] (i.e. N[Pi,100]//N) comes back at machine
                 * precision, not 100 digits. A specific precision is requested
                 * with the two-argument N[expr, p], which takes the MPFR branch
                 * below instead. A finite MPFR value can
                 * still exceed DBL_MAX (~1.8e308) — e.g. N[1.5 + 1001!], whose
                 * argument is already an MPFR ~4e2570 — or sit below DBL_MIN.
                 * mpfr_to_machine keeps those as a machine-precision
                 * (DBL_MANT_DIG-bit) MPFR, which has an arbitrary exponent. */
                return mpfr_to_machine(e->data.mpfr);
            }
            /* MPFR → MPFR. NUMERIC_MODE_MPFR sets the precision to exactly
             * spec.bits (SetPrecision pads up when the request exceeds the
             * value's current precision). NUMERIC_MODE_MPFR_CAP (the two-arg
             * N[x, p]) must never *increase* precision: cap the target at
             * min(existing, requested) so a 30-digit value stays 30 digits
             * under N[.., 50] but a 50-digit value still reduces under
             * N[.., 30]. */
            {
                /* An approximate zero has no precision, only accuracy (its
                 * bits), which N[x, p] does not change: N[0``80, 20] is
                 * 0``80, not a zero known only to 20 digits. */
                if (spec.mode == NUMERIC_MODE_MPFR_CAP && mpfr_zero_p(e->data.mpfr))
                    return expr_new_mpfr_copy(e->data.mpfr);
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
 *  Working precision
 *
 *  `N` replaces leaves and re-evaluates, so the leaf conversion is the one
 *  step that must never be the thing that loses information. Rounding an
 *  exact leaf to a double first costs a relative 2^-53, and a function
 *  amplifies that by its condition number — for Sin/Cos/Tan that factor is
 *  |x|, so N[Sin[3141592653589793238]] used to answer -0.641653 instead of
 *  -0.446315: the double nearest that integer is 3141592653589793280, and
 *  the 42 discarded bits are exactly the ones the argument reduction needs.
 *  Sin[3141592653589793238.] was right only because the parser builds that
 *  literal as a 62-bit MPFR straight from its text.
 *
 *  So: scan the input, work out the precision at which no exact leaf is
 *  lossy, evaluate there, and round the result *once* at the end. Both
 *  quantities below are derived, not tuned:
 *
 *    exact integer (and dyadic rational)  needs its mantissa bit count, at
 *      which it is represented exactly and contributes no error at all;
 *
 *    every other exact rational  needs out_bits + ceil(log2|x|) + guard,
 *      since a relative 2^-W on an argument of size |x| shows up in the
 *      answer amplified by |x|.
 *
 *  The guard is headroom for cancellation *within* the expression, which is
 *  the general significance-arithmetic problem and not what this is solving;
 *  64 bits matches the existing convention in bessel.c and contfrac.c.
 * ---------------------------------------------------------------------- */

/* The whole working-precision scan below — ExactScan and its helpers
 * (i64_lossy_bits, mpz_lossy_bits, scan_rational, scan_exact_leaves) plus
 * numeric_plan_working_spec / numeric_round_result — exists only to drive the
 * arbitrary-precision retry path, which is compiled only when MPFR is present.
 * Without MPFR every leaf is a machine double and there is nothing to scan, so
 * the entire cluster is guarded; otherwise it is dead code that trips
 * -Werror=unused-function on the no-MPFR build. */
#ifdef USE_MPFR
#define NUMERIC_GUARD_BITS 64

typedef struct {
    /* Max mantissa width over exact leaves that ARE finitely representable
     * (integers, dyadic rationals). Working precision must reach this for
     * them to survive the conversion intact. */
    long finite_bits;
    /* Max ceil(log2|x|) over exact leaves that are NOT finitely
     * representable (ordinary rationals). Their conversion error is
     * amplified by |x|, so they need that many bits *on top of* the output
     * precision. Zero for |x| < 2, where the amplification is nil. */
    long amplify_bits;
    /* Max log2|n| over exact exponents of Power[b, n]: the base's relative
     * error times |n|. Applied to arbitrary-precision requests only, so a
     * machine N[b^n] keeps its machine evaluation. */
    long power_amplify_bits;
    bool has_call;         /* a real function application, not just a number */
    bool has_infinity;     /* Infinity / ComplexInfinity / DirectedInfinity */
    bool has_zero;         /* an exact or inexact zero leaf */
} ExactScan;

/* Mantissa bits an int64 needs to be an exact double, or 0 when it already
 * is one. Values below 2^DBL_MANT_DIG always are; above it, only those whose
 * low bits are zero (2^53 yes, 2^53+1 no — and today both answer with the
 * same Sin, which is the bug in miniature). */
static long i64_lossy_bits(int64_t v) {
    if (v == 0) return 0;
    uint64_t a = (v < 0) ? (uint64_t)(-(v + 1)) + 1u : (uint64_t)v;
    long bits = 0;
    for (uint64_t t = a; t; t >>= 1) bits++;
    if (bits <= DBL_MANT_DIG) return 0;
    if ((a & ((((uint64_t)1) << (bits - DBL_MANT_DIG)) - 1u)) == 0) return 0;
    return bits;
}

/* Same question for an arbitrary-size integer. mpz_get_d truncates toward
 * zero, so the round-trip through a double is an exact representability
 * test; an out-of-range magnitude short-circuits it. */
static long mpz_lossy_bits(const mpz_t v) {
    if (mpz_sgn(v) == 0) return 0;
    size_t bits = mpz_sizeinbase(v, 2);
    if (bits <= (size_t)DBL_MANT_DIG) return 0;
    double d = mpz_get_d(v);
    if (isfinite(d)) {
        mpz_t back;
        mpz_init(back);
        mpz_set_d(back, d);
        int exact = (mpz_cmp(back, v) == 0);
        mpz_clear(back);
        if (exact) return 0;
    }
    return (long)bits;
}

static void scan_exact_leaves(const Expr* e, ExactScan* s);

/* Rational[n, d] with integer-like components, of any size. Dyadic
 * denominators are finitely representable and go in the `finite_bits`
 * bucket; everything else lands in `amplify_bits`. */
static long u64_bitlen(uint64_t a) {
    long bits = 0;
    for (; a; a >>= 1) bits++;
    return bits;
}

static void scan_rational(const Expr* num, const Expr* den, ExactScan* s) {
    /* Machine-integer parts (every Rational of a Table of 1/(i + j - 1),
     * ((i - 1)/(n - 1))^k, ...): the same three answers without a single
     * GMP allocation -- the mpq path below cost ~1 us a leaf, more than
     * the conversion N then does. */
    if (num->type == EXPR_INTEGER && den->type == EXPR_INTEGER &&
        num->data.integer != INT64_MIN && den->data.integer > 0) {
        int64_t n = num->data.integer;
        if (n == 0) { s->has_zero = true; return; }
        uint64_t a = n < 0 ? (uint64_t)(-n) : (uint64_t)n;
        uint64_t b = (uint64_t)den->data.integer;
        uint64_t x = a, y = b;
        while (y) { uint64_t t = x % y; x = y; y = t; }
        a /= x; b /= x;                                   /* canonical */
        if (b == 1) {
            long l = i64_lossy_bits((int64_t)a);
            if (l > s->finite_bits) s->finite_bits = l;
            return;
        }
        /* a / 2^k is its own double exactly when a is */
        if ((b & (b - 1)) == 0 && i64_lossy_bits((int64_t)a) == 0) return;
        long mag = u64_bitlen(a) - u64_bitlen(b);
        if (mag > s->amplify_bits) s->amplify_bits = mag;
        return;
    }
    mpq_t q;
    /* expr_to_mpz *initialises* its target (mpz_init_set), so it must be
     * handed raw slots. An mpq_init first would allocate limbs that the
     * re-initialisation then leaks. */
    expr_to_mpz(num, mpq_numref(q));
    expr_to_mpz(den, mpq_denref(q));
    if (mpz_sgn(mpq_denref(q)) == 0) {   /* malformed Rational[n, 0] */
        mpq_clear(q);
        return;
    }
    if (mpq_sgn(q) == 0) {
        s->has_zero = true;
        mpq_clear(q);
        return;
    }
    mpq_canonicalize(q);

    if (mpz_cmp_ui(mpq_denref(q), 1) == 0) {
        long l = mpz_lossy_bits(mpq_numref(q));
        if (l > s->finite_bits) s->finite_bits = l;
        mpq_clear(q);
        return;
    }

    double d = mpq_get_d(q);
    if (isfinite(d) && d != 0.0) {
        mpq_t back;
        mpq_init(back);
        mpq_set_d(back, d);
        int exact = (mpq_equal(back, q) != 0);
        mpq_clear(back);
        if (exact) {                       /* dyadic: the double is the value */
            mpq_clear(q);
            return;
        }
    }
    long mag = (long)mpz_sizeinbase(mpq_numref(q), 2)
             - (long)mpz_sizeinbase(mpq_denref(q), 2);
    if (mag > s->amplify_bits) s->amplify_bits = mag;
    mpq_clear(q);
}

static void scan_exact_leaves(const Expr* e, ExactScan* s) {
    if (!e) return;
    switch (e->type) {
        case EXPR_INTEGER: {
            if (e->data.integer == 0) { s->has_zero = true; break; }
            long l = i64_lossy_bits(e->data.integer);
            if (l > s->finite_bits) s->finite_bits = l;
            break;
        }
        case EXPR_BIGINT: {
            long l = mpz_lossy_bits(e->data.bigint);
            if (l > s->finite_bits) s->finite_bits = l;
            break;
        }
        case EXPR_REAL:
            /* A machine Real *is* its own exact binary value; there is
             * nothing to recover. Sin[1.*^25] is libm's answer for that
             * double, and that is the correct, Mathematica-matching result. */
            if (e->data.real == 0.0) s->has_zero = true;
            break;
        case EXPR_SYMBOL:
            if (e->data.symbol.name == SYM_Infinity
                || e->data.symbol.name == SYM_ComplexInfinity
                || e->data.symbol.name == SYM_Indeterminate) {
                s->has_infinity = true;
            }
            break;
        case EXPR_FUNCTION: {
            const Expr* head = e->data.function.head;
            /* Hold-forms are returned verbatim; never look inside. */
            if (is_hold_head(head)) return;

            if (head && head->type == EXPR_SYMBOL) {
                const char* h = head->data.symbol.name;
                if (h == SYM_DirectedInfinity) { s->has_infinity = true; return; }
                if (h == SYM_Rational && e->data.function.arg_count == 2
                    && expr_is_integer_like(e->data.function.args[0])
                    && expr_is_integer_like(e->data.function.args[1])) {
                    scan_rational(e->data.function.args[0],
                                  e->data.function.args[1], s);
                    return;
                }
                if (h == SYM_Complex) {
                    /* A structural wrapper, like Rational — descend without
                     * claiming this is a function that amplifies error. */
                    for (size_t i = 0; i < e->data.function.arg_count; ++i)
                        scan_exact_leaves(e->data.function.args[i], s);
                    return;
                }
                /* b^n with an exact exponent multiplies the relative error of
                 * b by |n|: (1 + 10^-40)^(10^40) at 98 working bits, where
                 * the base is rounded at 2^-98 and the exponent is 2^133,
                 * came back 6.27 instead of E. Such a power is left
                 * unevaluated only when n is huge, so its size is the
                 * amplification. */
                if (h == SYM_Power && e->data.function.arg_count == 2) {
                    const Expr* n = e->data.function.args[1];
                    long mag = 0;
                    if (expr_is_integer_like(n)) {
                        mpz_t z;
                        expr_to_mpz(n, z);
                        if (mpz_sgn(z) != 0) mag = (long)mpz_sizeinbase(z, 2);
                        mpz_clear(z);
                    } else if (n->type == EXPR_FUNCTION && n->data.function.head->type == EXPR_SYMBOL
                               && n->data.function.head->data.symbol.name == SYM_Rational
                               && n->data.function.arg_count == 2
                               && expr_is_integer_like(n->data.function.args[0])
                               && expr_is_integer_like(n->data.function.args[1])) {
                        mpz_t a, d;
                        expr_to_mpz(n->data.function.args[0], a);
                        expr_to_mpz(n->data.function.args[1], d);
                        if (mpz_sgn(a) != 0)
                            mag = (long)mpz_sizeinbase(a, 2) - (long)mpz_sizeinbase(d, 2);
                        mpz_clear(a);
                        mpz_clear(d);
                    }
                    if (mag > s->power_amplify_bits) s->power_amplify_bits = mag;
                }
            }
            s->has_call = true;
            scan_exact_leaves(head, s);
            for (size_t i = 0; i < e->data.function.arg_count; ++i)
                scan_exact_leaves(e->data.function.args[i], s);
            break;
        }
        default:
#ifdef USE_MPFR
            if (e->type == EXPR_MPFR) {
                if (mpfr_zero_p(e->data.mpfr)) s->has_zero = true;
            }
#endif
            break;
    }
}

/* Decide the working precision. Returns true (and fills *work) when the
 * evaluation must be raised, meaning the caller has to round the result
 * back down afterwards. */
static bool numeric_plan_working_spec(const ExactScan* s, NumericSpec spec,
                                      NumericSpec* work) {
    *work = spec;

    /* A bare number is already correctly rounded by the leaf conversion —
     * there is no function to amplify anything, so leave the fast path
     * untouched. */
    if (!s->has_call) return false;

    /* The precision the caller wants back, which is also the precision the
     * unraised evaluation would run at. */
    const long out_bits = numeric_spec_is_mpfr(spec) && spec.bits > 0
                        ? spec.bits : (long)DBL_MANT_DIG;

    long required = s->finite_bits;
    long amplify = s->amplify_bits;
    if (numeric_spec_is_mpfr(spec) && s->power_amplify_bits > amplify)
        amplify = s->power_amplify_bits;
    if (amplify > 0) {
        long amplified = out_bits + amplify + NUMERIC_GUARD_BITS;
        if (amplified > required) required = amplified;
    }

    /* An explicit N[expr, d] asks for d GOOD digits, and the amount a function
     * loses is invisible to the leaf scan above: it depends on the function's
     * CONDITIONING at the point, not on how the input is spelled. Measured, with
     * an exact rational input so no leaf was lossy and nothing was raised:
     *
     *   N[EllipticK[99999999999999999/10^17], 20]  ->  20.9581794072978496314
     *   true                                           20.958267651569278982883...
     *
     * -- six good digits of twenty, silently, because dK/dm ~ 1/(2(1-m)) ~ 2^56
     * there and the whole evaluation ran at exactly the 67 bits that 20 digits
     * asks for. Zeta at 1 + 10^-20 came back as 2^66 + 1 the same way. So give
     * every arbitrary-precision request the guard, not just the ones whose
     * LEAVES needed raising; `floor_bits` below is the headroom this scan
     * already intended for "cancellation within the expression", and a
     * conditioning loss is that same thing reached from the other side.
     *
     * Machine mode is deliberately excluded: 53 bits is all a double has, there
     * is no larger request to satisfy, and raising it would push every plain
     * N[Sin[0.5]] onto MPFR. The guard is 64 bits, so this fixes a loss of up to
     * ~19 digits and not an unbounded one -- a truly ill-conditioned point still
     * needs the caller to ask for more, exactly as $MaxExtraPrecision bounds the
     * same ladder elsewhere. */
    if (!numeric_spec_is_mpfr(spec) && required <= out_bits) return false;
    if (numeric_spec_is_mpfr(spec) && required <= out_bits
        && s->finite_bits <= 0 && amplify <= 0) {
        required = out_bits + NUMERIC_GUARD_BITS;
    }
    if (required <= out_bits) return false;   /* already precise enough */

    long floor_bits = out_bits + NUMERIC_GUARD_BITS;
    work->bits = required > floor_bits ? required : floor_bits;
    if (spec.mode == NUMERIC_MODE_MACHINE) {
        /* Plain MPFR, not MPFR_CAP: exact leaves must be produced at the
         * full working precision. A machine Real leaf is padded with zeros,
         * which is exactly right — its bit pattern is the value. */
        work->mode = NUMERIC_MODE_MPFR;
    }
    return true;
}

/* True when some EXPR_MPFR node would be changed by rounding to `out`. */
static bool numeric_needs_round(const Expr* e, NumericSpec out) {
    if (!e) return false;
    if (e->type == EXPR_MPFR) {
        if (out.mode == NUMERIC_MODE_MACHINE) return true;
        /* a zero's bits are its accuracy, which rounding to p digits of
         * precision does not touch (see numericalize's EXPR_MPFR branch) */
        if (mpfr_zero_p(e->data.mpfr)) return false;
        return out.bits > 0 && (long)mpfr_get_prec(e->data.mpfr) > out.bits;
    }
    if (e->type != EXPR_FUNCTION) return false;
    if (numeric_needs_round(e->data.function.head, out)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; ++i) {
        if (numeric_needs_round(e->data.function.args[i], out)) return true;
    }
    return false;
}

/* Round a working-precision result down to what the caller asked for.
 * Rewrites EXPR_MPFR nodes only — structure, symbols and exact numbers pass
 * through, and nothing is re-evaluated. Consumes `e` and returns a value the
 * caller owns. Nodes are rebuilt rather than written through: expr_copy is a
 * refcount bump, so the tree may be shared. */
static Expr* numeric_round_result(Expr* e, NumericSpec out) {
    if (!e || !numeric_needs_round(e, out)) return e;

    if (e->type == EXPR_MPFR) {
        Expr* r;
        if (out.mode == NUMERIC_MODE_MACHINE) {
            r = mpfr_to_machine(e->data.mpfr);
        } else {
            r = expr_new_mpfr_bits(out.bits);
            if (r) mpfr_set(r->data.mpfr, e->data.mpfr, MPFR_RNDN);
        }
        if (!r) return e;
        expr_free(e);
        return r;
    }

    const size_t n = e->data.function.arg_count;
    Expr** args = n ? (Expr**)malloc(n * sizeof(Expr*)) : NULL;
    if (n && !args) return e;
    Expr* head = numeric_round_result(expr_copy(e->data.function.head), out);
    for (size_t i = 0; i < n; ++i) {
        args[i] = numeric_round_result(expr_copy(e->data.function.args[i]), out);
    }
    /* A machine complex number has two machine parts: the exact zero part an
     * arbitrary-precision kernel keeps (ArcCosh[x] = I ArcCos[x], |x| < 1)
     * is 0. once the other part is a machine Real, as machine N always gave. */
    if (out.mode == NUMERIC_MODE_MACHINE && n == 2 && head->type == EXPR_SYMBOL
        && head->data.symbol.name == SYM_Complex) {
        for (size_t i = 0; i < 2; ++i) {
            if (args[i]->type == EXPR_INTEGER && args[i]->data.integer == 0
                && args[1 - i]->type == EXPR_REAL) {
                expr_free(args[i]);
                args[i] = expr_new_real(0.0);
            }
        }
    }
    Expr* rebuilt = expr_new_function(head, args, n);
    free(args);
    if (!rebuilt) return e;
    expr_free(e);
    return rebuilt;
}

/* True when some EXPR_REAL in the tree is Inf or NaN. */
static bool expr_has_nonfinite_real(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return !isfinite(e->data.real);
    if (e->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < e->data.function.arg_count; ++i) {
        if (expr_has_nonfinite_real(e->data.function.args[i])) return true;
    }
    return false;
}

/* True when some EXPR_REAL in the tree is an exact zero. */
static bool expr_has_zero_real(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return e->data.real == 0.0;
    if (e->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < e->data.function.arg_count; ++i) {
        if (expr_has_zero_real(e->data.function.args[i])) return true;
    }
    return false;
}

/* True when the tree carries an Overflow[] or Indeterminate. The evaluator's
 * inf/nan net converts a leaked machine overflow into one of these *during* the
 * machine-mode inner evaluation, so N[]'s own MPFR retry (which used to detect
 * the raw Inf/NaN) must recognise them as the same degeneracy signal. */
static bool expr_has_overflow_or_indeterminate(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_SYMBOL) return e->data.symbol.name == SYM_Indeterminate;
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head && e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Overflow) return true;
    for (size_t i = 0; i < e->data.function.arg_count; ++i) {
        if (expr_has_overflow_or_indeterminate(e->data.function.args[i])) return true;
    }
    return false;
}

/* True when some MPFR leaf is finite and nonzero yet its magnitude is outside
 * the machine-double range — a value only an extended-exponent (promoted)
 * representation can hold (Exp[1000.] -> 1.97*10^434, Exp[-1000.] -> 1.9*10^-434).
 * This distinguishes a promoted overflow/underflow, which a machine-mode N[]
 * should recompute from the exact form for full accuracy, from an ordinary
 * in-range MPFR result (which it should leave alone). */
static bool expr_has_extended_range_real(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_MPFR) {
        if (mpfr_zero_p(e->data.mpfr) || mpfr_nan_p(e->data.mpfr) || mpfr_inf_p(e->data.mpfr))
            return false;
        double d = mpfr_get_d(e->data.mpfr, MPFR_RNDN);
        return !isfinite(d) || d == 0.0;   /* overflowed or underflowed a double */
    }
    if (e->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < e->data.function.arg_count; ++i)
        if (expr_has_extended_range_real(e->data.function.args[i])) return true;
    return false;
}
#endif  /* USE_MPFR */

/* ------------------------------------------------------------------------
 *  inf / nan never escape a numeric result
 *
 *  A raw IEEE infinity or NaN is an artifact of the C math library, not a
 *  Mathilda value: Exp[800.] overflowing a double must never surface as
 *  `inf.0`, nor Binomial[2000., 1000.] as `nan`. The evaluator routes every
 *  NumericFunction result through numeric_sanitize_nonfinite() as a
 *  last-resort backstop, rewriting a leaked non-finite machine real into a
 *  genuine construct:
 *      +Inf -> Overflow[]      -Inf -> -Overflow[]      NaN -> Indeterminate
 *  The representable-overflow cases (Exp[800.] -> 2.7*10^347, ...) are promoted
 *  to an extended-exponent MPFR real upstream in each head's machine-real path,
 *  so this net only fires for magnitudes past even MPFR's exponent range, where
 *  Overflow[] is the correct answer (as it is in Mathematica). Always compiled:
 *  the "no inf/nan" invariant must hold with or without MPFR.
 * ---------------------------------------------------------------------- */

static Expr* make_overflow_expr(void) {
    return expr_new_function(expr_new_symbol(SYM_Overflow), NULL, 0);
}

/* Cheap top-level test used to gate the rewrite: is `e` a machine real, or a
 * Complex of machine reals, carrying an Inf or NaN? The common finite path
 * pays only this type check. */
/* A single leaf carrying an Inf or NaN (machine or MPFR). */
static bool leaf_is_nonfinite(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return !isfinite(e->data.real);
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return mpfr_nan_p(e->data.mpfr) || mpfr_inf_p(e->data.mpfr);
#endif
    return false;
}

bool numeric_result_has_nonfinite(const Expr* e) {
    if (!e) return false;
    if (leaf_is_nonfinite(e)) return true;
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im)) {
        if (leaf_is_nonfinite(re) || leaf_is_nonfinite(im)) return true;
    }
    return false;
}

/* NaN? / (positive) Inf? for a single machine-or-MPFR leaf. */
static bool leaf_is_nan(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return isnan(e->data.real);
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return mpfr_nan_p(e->data.mpfr) != 0;
#endif
    return false;
}
static bool leaf_is_neg_inf(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return isinf(e->data.real) && e->data.real < 0.0;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return mpfr_inf_p(e->data.mpfr) && mpfr_sgn(e->data.mpfr) < 0;
#endif
    return false;
}

/* Rewrite a leaked non-finite machine real into a Mathilda construct. Takes
 * ownership of `e` and returns the replacement (freeing `e`), or `e` unchanged
 * when it holds nothing non-finite — so it is safe to call unconditionally,
 * though the evaluator gates it on numeric_result_has_nonfinite() to keep the
 * hot path free. */
/* Build the construct a non-finite leaf collapses to. `nan` wins over sign. */
static Expr* nonfinite_construct(bool nan, bool neg) {
    if (nan) return expr_new_symbol(SYM_Indeterminate);
    if (neg) {
        Expr* a[2] = { expr_new_integer(-1), make_overflow_expr() };
        return expr_new_function(expr_new_symbol(SYM_Times), a, 2);   /* -Overflow[] */
    }
    return make_overflow_expr();
}

Expr* numeric_sanitize_nonfinite(Expr* e) {
    if (!e) return e;
    if (leaf_is_nonfinite(e)) {
        Expr* r = nonfinite_construct(leaf_is_nan(e), leaf_is_neg_inf(e));
        expr_free(e);
        return r;
    }
    Expr *re, *im;
    if (is_complex(e, &re, &im)) {
        if (leaf_is_nonfinite(re) || leaf_is_nonfinite(im)) {
            /* A complex overflow has no meaningful sign; NaN anywhere -> Indeterminate. */
            Expr* r = nonfinite_construct(leaf_is_nan(re) || leaf_is_nan(im), false);
            expr_free(e);
            return r;
        }
    }
    return e;
}

#ifdef USE_MPFR
/* ------------------------------------------------------------------------
 *  L2: promote a degenerate machine-real call to an extended-exponent real
 *
 *  Every head with a machine-real fast path also has an MPFR path (N[f[exact]]
 *  promotes). When the machine path yields a degenerate double — Inf/NaN from
 *  overflow, or 0 from underflow — the fix is to recompute the SAME call with
 *  its machine-real arguments promoted to MPFR (whose exponent range is
 *  effectively unbounded) and present the result at machine precision. This is
 *  the single mechanism behind Gamma[256.] -> 3.35*10^504, Exp[800.], Sinh,
 *  Binomial, BernoulliB, ... : each head calls numeric_promote_machine_call()
 *  from its own (cold) degeneracy branch, where the arguments are still live.
 * ---------------------------------------------------------------------- */

/* Promote a machine-real leaf (EXPR_REAL, or a Complex of them) to an MPFR at
 * `bits`; copy any other atom unchanged. Sets *saw when a machine real was
 * actually promoted, so a call with nothing to promote can bail cheaply. */
static Expr* promote_leaf_to_mpfr(const Expr* e, mpfr_prec_t bits, bool* saw) {
    if (!e) return NULL;
    if (e->type == EXPR_REAL) {
        *saw = true;
        Expr* m = expr_new_mpfr_bits(bits);
        mpfr_set_d(m->data.mpfr, e->data.real, MPFR_RNDN);
        return m;
    }
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im) &&
        (re->type == EXPR_REAL || im->type == EXPR_REAL)) {
        return make_complex(promote_leaf_to_mpfr(re, bits, saw),
                            promote_leaf_to_mpfr(im, bits, saw));
    }
    /* Recurse into function arguments so machine reals inside List parameters
     * are promoted too (HypergeometricPFQ[{1.}, {2.}, 900.] -> all-MPFR). The
     * call is a numeric computation, so every real leaf is fair game. */
    if (e->type == EXPR_FUNCTION) {
        size_t n = e->data.function.arg_count;
        Expr** a = (Expr**)malloc(n * sizeof(Expr*));
        for (size_t i = 0; i < n; ++i)
            a[i] = promote_leaf_to_mpfr(e->data.function.args[i], bits, saw);
        Expr* out = expr_new_function(expr_copy(e->data.function.head), a, n);
        free(a);
        return out;
    }
    return expr_copy((Expr*)e);
}

/* Classify a numeric leaf: finite? and, when finite, a (signed) zero?
 * A symbolic / non-numeric result is reported as not finite. */
static bool numeric_leaf_is_finite(const Expr* e, bool* is_zero) {
    *is_zero = false;
    if (!e) return false;
    if (e->type == EXPR_REAL) { *is_zero = (e->data.real == 0.0); return isfinite(e->data.real); }
    if (e->type == EXPR_MPFR) {
        if (mpfr_nan_p(e->data.mpfr) || mpfr_inf_p(e->data.mpfr)) return false;
        *is_zero = (mpfr_zero_p(e->data.mpfr) != 0);
        return true;
    }
    if (e->type == EXPR_INTEGER) { *is_zero = (e->data.integer == 0); return true; }
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im)) {
        bool zr, zi;
        bool fr = numeric_leaf_is_finite(re, &zr);
        bool fi = numeric_leaf_is_finite(im, &zi);
        *is_zero = zr && zi;
        return fr && fi;
    }
    return false;
}

/* Round an MPFR result (or the parts of a Complex) down to a machine 53-bit
 * mantissa in place: a value computed at a higher working precision is then
 * presented as a machine number (Precision ~ 15.95), matching N[f[exact]]. The
 * exponent is untouched, so an overflowing magnitude survives the round. */
static void round_leaf_to_machine(Expr* e) {
    if (!e) return;
    if (e->type == EXPR_MPFR) { mpfr_prec_round(e->data.mpfr, 53, MPFR_RNDN); return; }
    Expr *re, *im;
    if (is_complex(e, &re, &im)) { round_leaf_to_machine(re); round_leaf_to_machine(im); }
}

/* If a promoted pure-real MPFR value sits inside the normal machine-double range
 * (Binomial[1000., 500.] -> 2.7*10^299, whose tgamma ratio only overflowed as an
 * intermediate), hand it back as a machine Real so MachineNumberQ agrees with
 * Mathematica. A magnitude past DBL_MAX (Gamma[256.]) stays MPFR; one below
 * DBL_MIN (Exp[-800.] -> 3.7*10^-348) stays MPFR too, since a double would
 * re-underflow it to zero. Complex results are left as MPFR. */
static Expr* demote_to_machine_if_fits(Expr* r) {
    if (r && r->type == EXPR_MPFR && !mpfr_zero_p(r->data.mpfr)) {
        double d = mpfr_get_d(r->data.mpfr, MPFR_RNDN);
        if (isfinite(d) && d != 0.0) {   /* normal double: exact for a 53-bit source */
            expr_free(r);
            return expr_new_real(d);
        }
    }
    return r;
}

/* Public: lower a (consumed) MPFR or Complex[MPFR|Real, MPFR|Real] value to
 * machine Reals when every part fits a normal double (a part that is exactly
 * zero is fine); otherwise return it unchanged -- the extended-exponent
 * convention for a machine call whose value lies outside the double range.
 * For heads that compute a machine-precision call through arb (53-bit MPFR
 * out). */
Expr* numeric_lower_mpfr_to_machine(Expr* r) {
    if (!r) return r;
    Expr *re, *im;
    bool cx = is_complex(r, &re, &im);
    if (!cx) { re = r; im = NULL; }
    Expr* parts[2] = { re, im };
    double d[2] = { 0.0, 0.0 };
    for (int i = 0; i < 2; i++) {
        Expr* p = parts[i];
        if (!p) continue;
        if (p->type == EXPR_REAL) { d[i] = p->data.real; continue; }
        if (p->type == EXPR_INTEGER && p->data.integer == 0) continue;
        if (p->type != EXPR_MPFR) return r;
        if (mpfr_zero_p(p->data.mpfr)) continue;
        d[i] = mpfr_get_d(p->data.mpfr, MPFR_RNDN);
        if (!isfinite(d[i]) || fabs(d[i]) < DBL_MIN) return r;
    }
    expr_free(r);
    if (!cx) return expr_new_real(d[0]);
    return make_complex(expr_new_real(d[0]), expr_new_real(d[1]));
}

Expr* numeric_promote_machine_call(const Expr* call) {
    if (!call || call->type != EXPR_FUNCTION) return NULL;
    size_t n = call->data.function.arg_count;
    if (n == 0) return NULL;

    /* Machine precision first (so the output matches N[f[exact]]'s ~15.95
     * digits); a high-precision retry only if that underflowed to zero — some
     * MPFR kernels (BesselK) drop the exponent at 53 bits though the value is
     * representable. Rounding the retry back to 53 bits keeps machine output. */
    const mpfr_prec_t passes[2] = { 53, (mpfr_prec_t)(53 + NUMERIC_GUARD_BITS + 64) };
    for (int p = 0; p < 2; ++p) {
        mpfr_prec_t bits = passes[p];
        bool saw = false;
        Expr** args = (Expr**)malloc(n * sizeof(Expr*));
        if (!args) return NULL;
        for (size_t i = 0; i < n; ++i)
            args[i] = promote_leaf_to_mpfr(call->data.function.args[i], bits, &saw);
        if (!saw) {                       /* no machine real to promote */
            for (size_t i = 0; i < n; ++i) expr_free(args[i]);
            free(args);
            return NULL;
        }
        Expr* promoted = expr_new_function(expr_copy(call->data.function.head), args, n);
        free(args);
        Expr* r = eval_and_free(promoted);

        bool is_zero = false;
        bool finite = numeric_leaf_is_finite(r, &is_zero);

        if (finite && !is_zero) {
            if (bits > 53) round_leaf_to_machine(r);
            return demote_to_machine_if_fits(r);   /* machine when representable */
        }
        /* Degenerate at this precision. The machine-precision (53-bit) pass can
         * be misleading: a kernel may flush to zero (BesselK), or take an
         * internal machine path that overflows though the value is representable
         * (Erfi, ExpIntegralEi). Retry once at a high working precision before
         * concluding anything. */
        if (p == 0) { expr_free(r); continue; }

        /* Second pass still degenerate: this is the real shape of the value. */
        if (finite && is_zero) { expr_free(r); return NULL; }  /* true zero: keep machine 0 */
        if (r && (r->type == EXPR_REAL || r->type == EXPR_MPFR)) {
            /* Non-finite even in MPFR: past the exponent range -> Overflow[] etc. */
            return numeric_sanitize_nonfinite(r);
        }
        expr_free(r);                     /* symbolic decline: keep caller's behaviour */
        return NULL;
    }
    return NULL;
}
#else  /* !USE_MPFR: no arbitrary-precision path to promote into */
Expr* numeric_promote_machine_call(const Expr* call) { (void)call; return NULL; }
Expr* numeric_lower_mpfr_to_machine(Expr* r) { return r; }
#endif /* USE_MPFR */

/* Wrap the pure-real machine result `v` of a unary numeric `call` (head[x]),
 * where `xin` is the real input. A degenerate `v` — Inf/NaN from overflow, or a
 * flush-to-zero from a nonzero argument — is recomputed in MPFR via
 * numeric_promote_machine_call so Exp[800.], Sinh[800.], Sech[800.], ... return
 * the true extended-exponent value instead of leaking a libm artifact. In a
 * build without MPFR this is exactly expr_new_real(v) (and the evaluator's net
 * still stops any Inf/NaN from escaping). */
Expr* numeric_machine_real_or_promote(double v, double xin, const Expr* call) {
#ifdef USE_MPFR
    if (!isfinite(v) || (v == 0.0 && xin != 0.0)) {
        Expr* p = numeric_promote_machine_call(call);
        if (p) return p;
    }
#else
    (void)xin; (void)call;
#endif
    return expr_new_real(v);
}

#ifdef USE_MPFR
/* True when `call` has at least one non-zero machine-real argument (the usual
 * signature of an argument big/small enough to have driven an over/underflow).
 * Only consulted on the MPFR promotion path. */
static bool call_has_nonzero_real_arg(const Expr* call) {
    if (!call || call->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < call->data.function.arg_count; ++i) {
        Expr* a = call->data.function.args[i];
        if (a && a->type == EXPR_REAL && a->data.real != 0.0) return true;
        Expr *re, *im;
        if (a && is_complex(a, &re, &im)) {
            if (re->type == EXPR_REAL && re->data.real != 0.0) return true;
            if (im->type == EXPR_REAL && im->data.real != 0.0) return true;
        }
    }
    return false;
}
#endif

/* Wrap the *result* `r` of a numeric `call` computed via a machine fast path.
 * If `r` is degenerate — Inf/NaN, or a flush-to-zero while the call has a
 * non-zero real argument — recompute the call in MPFR (numeric_promote_machine_call)
 * and return that, freeing `r`. Otherwise return `r` unchanged. For heads whose
 * numeric core is a helper that only sees the argument (Erfc/Erfi/ExpIntegralEi,
 * ...), so the argument-live promotion is applied here at the builtin entry.
 * Always compiled; a no-op without MPFR (the evaluator's net still catches Inf/
 * NaN). NULL `r` (declined) is returned unchanged. */
Expr* numeric_promote_result_if_degenerate(Expr* r, const Expr* call) {
#ifdef USE_MPFR
    if (!r) return r;
    bool degen = false;
    if (r->type == EXPR_REAL) {
        double v = r->data.real;
        degen = !isfinite(v) || (v == 0.0 && call_has_nonzero_real_arg(call));
    } else {
        Expr *re, *im;
        if (is_complex(r, &re, &im))
            degen = leaf_is_nonfinite(re) || leaf_is_nonfinite(im);
    }
    if (degen) {
        Expr* p = numeric_promote_machine_call(call);
        if (p) { expr_free(r); return p; }
    }
#else
    (void)call;
#endif
    return r;
}

/* ------------------------------------------------------------------------
 *  Top-level dispatch
 * ---------------------------------------------------------------------- */
static Expr* numericalize_top(const Expr* e, NumericSpec spec);

/* Set by the N builtin for the numericalization it is about to run and
 * consumed by the first numericalize_top that sees it, so a nested one does
 * not. A part that vanishes at the precision limit is reported (N::meprec)
 * only to N's own caller: it cannot be told from a nonzero value below the
 * noise there, and N[expr, p] promised p digits of it. Contagion inside
 * Plus/Times (N[1, 30] + (Sin[Pi/7]^2 + Cos[Pi/7]^2 - 1)), which only needs
 * the part to the precision of its partner, keeps it as a silent 0. */
static bool g_numeric_from_n = false;

/* Count of adaptive re-evaluations (numeric.h). */
static long g_numeric_reevaluations = 0;
long numeric_adaptive_reevaluations(void) { return g_numeric_reevaluations; }


#ifdef USE_MPFR
/* ---- adaptive precision for N[expr, p] ------------------------------------
 *
 * N[expr, p] promises p good digits. One evaluation at a fixed working
 * precision cannot keep that promise once the expression cancels: the result
 * then has as many good bits as the working precision exceeds the
 * cancellation, and they are reported as p digits regardless. Measured at
 * p = 32 (170 working bits): N[Cos[10^-20] - 1] had 11 good digits of 32 and
 * at p = 20 came back as an exact 0.0 (true value -5*10^-41); the exact zeros
 * E^(I Pi/4) - (-1)^(1/4) and Sin[Pi/7]^2 + Cos[Pi/7]^2 - 1 came back as
 * 10^-52-size "values".
 *
 * So N re-evaluates at a higher working precision W until the result is
 * certified, up to p + $MaxExtraPrecision digits (plus the guard). A result at
 * W is certified when
 *   (a) it agrees with an evaluation at least 32 bits lower to p bits, either
 *       directly or by extrapolation: a loss of L bits leaves W - L good bits
 *       at every W, so a lower result that agrees on A bits has lost
 *       W_lo - A and the upper one has A + (W - W_lo) good bits (less a
 *       16-bit margin, and only from A >= 32 so noise cannot pass); and
 *   (b) every cancelling addition of the pass itself kept at least p bits
 *       (g_numeric_cancel_left): one that cancelled to an exact zero kept
 *       none, so two zeros never certify each other.
 * (b) is what makes an exact zero honest. E^1000 (E^(1/1000 - E^-1000) -
 * E^(1/1000)) is an exact 0.0 at every W below ~1450 bits (E^-1000 is absorbed
 * by 1/1000, the two exponentials are then equal) and two such zeros "agreed",
 * so N answered 0.0 for -1.001...; (a)'s 32-bit gap is what keeps a final
 * short step to the cap (330 -> 331 bits) from agreeing on rounding noise,
 * which certified 18 wrong digits of N[(Log[Cos[x]] + x^2/2 + x^4/12 +
 * x^6/45)/x^8 /. x -> 10^-10, 30].
 *
 * Each step goes to at least W + 64 and, when the pass cancelled L < W bits,
 * to p + L + 32, the precision that cancellation needs; after a cancellation
 * to zero (L unknown) the extra precision doubles.
 *
 * At the cap a part is kept if it is certified by (a). A part that is zero, or
 * keeps shrinking as W rises, is 0 to the precision reached: it is returned as
 * 0. and, to the N builtin's caller, N::meprec says it is not certified -- a
 * value below the noise at the limit (E^1000 (E^(1/1000 - E^-1000) -
 * E^(1/1000)), or (1 + 10^-40)^(10^40) - E at 10 digits) shrinks with W
 * exactly as the residue of an algebraic zero does. Any other unsettled part
 * is kept at the bits it can claim, with N::meprec. */

/* Leading bits on which two MPFR values agree (a large number when equal). */
static long mpfr_agree_bits(const mpfr_t a, const mpfr_t b) {
    if (mpfr_zero_p(a) && mpfr_zero_p(b)) return 1L << 30;
    if (mpfr_zero_p(a) || mpfr_zero_p(b) || !mpfr_number_p(a) || !mpfr_number_p(b)) return 0;
    mpfr_t d;
    mpfr_init2(d, mpfr_get_prec(a) > mpfr_get_prec(b) ? mpfr_get_prec(a) : mpfr_get_prec(b));
    mpfr_sub(d, a, b, MPFR_RNDN);
    long r = mpfr_zero_p(d) ? (1L << 30) : (long)(mpfr_get_exp(a) - mpfr_get_exp(d));
    mpfr_clear(d);
    return r;
}

/* Minimum agreement over corresponding MPFR leaves of two results; -1 when
 * their shapes differ (a part changed kind between the two precisions). */
static long expr_agree_bits(const Expr* a, const Expr* b) {
    if (!a || !b) return -1;
    if (a->type == EXPR_MPFR || b->type == EXPR_MPFR) {
        if (a->type != EXPR_MPFR || b->type != EXPR_MPFR) return -1;
        return mpfr_agree_bits(a->data.mpfr, b->data.mpfr);
    }
    if (a->type != EXPR_FUNCTION || b->type != EXPR_FUNCTION)
        return expr_eq((Expr*)a, (Expr*)b) ? (1L << 30) : -1;
    if (a->data.function.arg_count != b->data.function.arg_count
        || !expr_eq(a->data.function.head, b->data.function.head)) return -1;
    long m = 1L << 30;
    for (size_t i = 0; i < a->data.function.arg_count; i++) {
        long k = expr_agree_bits(a->data.function.args[i], b->data.function.args[i]);
        if (k < 0) return -1;
        if (k < m) m = k;
    }
    return m;
}

/* Bits of `bits` an evaluation at W can claim from agreeing on `agree`
 * leading bits with one `gap` bits lower (see the block comment above):
 * `agree` itself, or the extrapolated agree + gap - 16 when agree >= 32. */
static long adaptive_good_bits(long agree, long gap) {
    if (agree >= (1L << 29)) return agree;          /* identical */
    if (agree >= 32 && gap >= 32) {
        long g = agree + gap - 16;
        if (g > agree) return g;
    }
    return agree;
}

/* The verdict at the precision cap, part by part: `cur` (at the cap) against
 * `ref` (at least 32 bits lower) and `first` (the first evaluation). A fresh
 * tree, never mutating. A certified MPFR leaf is kept. A zero or shrinking
 * one -- or one whose `ref` is zero while some addition of the pass at the cap
 * kept at most 16 bits (`total_cancel`: noise in the last bits) -- becomes a
 * zero of `bits` precision and sets *vanished. Any other unsettled leaf is
 * kept at the bits it can claim (at least 1) and sets *unsettled, as does a
 * part whose shape changed. */
static Expr* expr_cap_verdict(const Expr* cur, const Expr* ref, const Expr* first,
                              long gap, long bits, bool total_cancel,
                              bool* vanished, bool* unsettled) {
    if (cur->type == EXPR_MPFR && ref && ref->type == EXPR_MPFR) {
        const mpfr_srcptr c = cur->data.mpfr, r = ref->data.mpfr;
        if (!mpfr_zero_p(c)) {
            long good = adaptive_good_bits(mpfr_agree_bits(c, r), gap);
            if (good >= bits) return expr_copy((Expr*)cur);
            /* shrinking: smaller than the lower evaluation by half the
             * precision gained (noise of size 2^-W shrinks by all of it) and
             * smaller than the first evaluation too */
            long drop = gap / 2 > 16 ? gap / 2 : 16;
            bool shrink;
            if (mpfr_zero_p(r)) {
                shrink = total_cancel;
            } else {
                shrink = mpfr_regular_p(r)
                    && mpfr_get_exp(c) <= mpfr_get_exp(r) - drop
                    && (!first || first->type != EXPR_MPFR
                        || !mpfr_regular_p(first->data.mpfr)
                        || mpfr_get_exp(c) <= mpfr_get_exp(first->data.mpfr) - 16);
            }
            if (!shrink) {
                *unsettled = true;
                if (good < 1) good = 1;
                /* `good` bits are known: the precision, stored with its guard */
                long kw = numeric_precision_bits_width(good);
                Expr* k = expr_new_mpfr_bits((mpfr_prec_t)(kw < bits ? kw : bits));
                if (!k) return expr_copy((Expr*)cur);
                mpfr_set(k->data.mpfr, c, MPFR_RNDN);
                return k;
            }
        }
        *vanished = true;
        Expr* z = expr_new_mpfr_bits((mpfr_prec_t)bits);
        if (z) mpfr_set_zero(z->data.mpfr, 1);
        return z ? z : expr_copy((Expr*)cur);
    }
    if (cur->type != EXPR_FUNCTION || !ref || ref->type != EXPR_FUNCTION
        || cur->data.function.arg_count != ref->data.function.arg_count) {
        if (expr_agree_bits(cur, ref) < bits) *unsettled = true;
        return expr_copy((Expr*)cur);
    }
    bool first_ok = first && first->type == EXPR_FUNCTION
        && first->data.function.arg_count == cur->data.function.arg_count;
    size_t n = cur->data.function.arg_count;
    Expr** args = (Expr**)malloc((n ? n : 1) * sizeof(Expr*));
    if (!args) return expr_copy((Expr*)cur);
    for (size_t i = 0; i < n; i++)
        args[i] = expr_cap_verdict(cur->data.function.args[i], ref->data.function.args[i],
                                   first_ok ? first->data.function.args[i] : NULL,
                                   gap, bits, total_cancel, vanished, unsettled);
    Expr* out = expr_new_function(expr_copy(cur->data.function.head), args, n);
    free(args);
    return out;
}

/* True when `e` contains a Plus: the only way an evaluation can cancel. */
static bool expr_has_plus(const Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Plus) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (expr_has_plus(e->data.function.args[i])) return true;
    return false;
}

static long numeric_max_extra_bits(double* digits_out) {
    Rule* r = symtab_get_own_values("$MaxExtraPrecision");
    double digits = 50.0;
    if (r && r->replacement) {
        if (r->replacement->type == EXPR_REAL) digits = r->replacement->data.real;
        else if (r->replacement->type == EXPR_INTEGER) digits = (double)r->replacement->data.integer;
    }
    if (!(digits >= 0.0)) digits = 0.0;
    if (digits > 100000.0) digits = 100000.0;
    if (digits_out) *digits_out = digits;
    return (long)ceil(digits * 3.3219280948873622);   /* log2(10) */
}
#endif

Expr* numericalize(const Expr* e, NumericSpec spec) {
    int saved = g_n_held_depth;
    bool saved_contagion = g_n_contagion;
    g_n_held_depth = 0;
    g_n_contagion = false;
    Expr* r = numericalize_top(e, spec);
    g_n_held_depth = saved;
    g_n_contagion = saved_contagion;
    return r;
}

/* numericalize_machine_unplanned: machine-precision numericalization WITHOUT
 * the working-precision plan.  numericalize raises a machine request to MPFR
 * when an exact leaf could lose accuracy in conversion (any non-dyadic rational
 * |q| >= 2: 10/3 lifts the whole evaluation to ~119 bits), which is what N
 * promises.  A caller that only screens -- PossibleZeroQ's machine rung, whose
 * threshold is relative (2^-12) and which escalates any doubtful point up its
 * own MPFR ladder -- pays ~10x for nothing: on a 10k-leaf antiderivative
 * residual, 85 ms a point against a few. */
Expr* numericalize_machine_unplanned(const Expr* e) {
    int saved = g_n_held_depth;
    bool saved_contagion = g_n_contagion;
    g_n_held_depth = 0;
    g_n_contagion = false;
    Expr* r = numericalize_pass(e, numeric_machine_spec());
    g_n_held_depth = saved;
    g_n_contagion = saved_contagion;
    return r;
}

/* numericalize for an exact operand of an inexact Plus / Times: as above,
 * but with g_n_contagion on for the duration. */
static Expr* numericalize_contagion(const Expr* e, NumericSpec spec) {
    int saved = g_n_held_depth;
    bool saved_contagion = g_n_contagion;
    g_n_held_depth = 0;
    g_n_contagion = true;
    Expr* r = numericalize_top(e, spec);
    g_n_held_depth = saved;
    g_n_contagion = saved_contagion;
    return r;
}

/* One evaluation of `e` at `spec` with its own subexpression memo (see
 * NMemo), restoring the enclosing pass's memo afterwards. */
static Expr* numericalize_pass(const Expr* e, NumericSpec spec) {
    NMemo memo = { NULL, 0, 0 };
    NMemo* saved = g_nmemo;
    /* Machine-precision passes run without it: a machine call costs about
     * what the lookup does, so only an arbitrary-precision pass gains. */
#ifdef USE_MPFR
    g_nmemo = numeric_spec_is_mpfr(spec) ? &memo : NULL;
#else
    g_nmemo = NULL;
#endif
    Expr* r = numericalize_rec(e, spec);
    g_nmemo = saved;
    nmemo_clear(&memo);
    return r;
}

static Expr* numericalize_top(const Expr* e, NumericSpec spec) {
    if (!e) return NULL;

#ifndef USE_MPFR
    return numericalize_pass(e, spec);
#else
    const bool from_n = g_numeric_from_n;
    g_numeric_from_n = false;
    ExactScan s;
    memset(&s, 0, sizeof(s));
    scan_exact_leaves(e, &s);

    NumericSpec work;
    bool raised = numeric_plan_working_spec(&s, spec, &work);
    long saved_cancel = g_numeric_cancel_bits;
    long saved_left = g_numeric_cancel_left;
    long saved_absorb = g_numeric_absorb_bits;
    g_numeric_cancel_bits = 0;
    g_numeric_absorb_bits = 0;
    Expr* r = numericalize_pass(e, work);
    long cancelled = g_numeric_cancel_bits;
    long absorbed = g_numeric_absorb_bits;
    /* Adaptive precision is the promise of the N builtin's N[expr, p] (the
     * capping mode is produced solely by N); internal working-precision
     * numericalizations keep a single evaluation. */
    /* The headroom a cancellation may use is what the working precision has
     * above p LESS the bits planned for amplification by a large exact leaf
     * or exponent: those are spent by the leaf, not available to the sum.
     * N[(1 + 10^-40)^(10^40) - E, 10] ran at 231 bits, 133 of them taken by
     * the power, and its 133-bit cancellation passed as "within the guard". */
    long amp = s.amplify_bits > s.power_amplify_bits ? s.amplify_bits : s.power_amplify_bits;
    if (amp < 0 || !raised) amp = 0;
    long headroom = (raised ? work.bits : spec.bits) - spec.bits - amp;
    if (headroom < 0) headroom = 0;
    bool adaptive = spec.mode == NUMERIC_MODE_MPFR_CAP && spec.bits > 0 && r && expr_has_plus(e)
        && cancelled > headroom;   /* lost past the guard */
    Expr* r_raw = adaptive ? expr_copy(r) : NULL;   /* the result at W, unrounded */
    if (raised) r = numeric_round_result(r, spec);

    /* Adaptive precision (see the block comment above mpfr_agree_bits): only
     * an N[expr, p] whose expression cancelled pays for a second evaluation. */
    if (adaptive) {
        const long p = spec.bits;
        double mep_digits = 50.0;
        const long Wmax = p + numeric_max_extra_bits(&mep_digits) + NUMERIC_GUARD_BITS;
        /* The two most recent results, newest first, with their precisions
         * and pass cancellations; `first` is the shrink reference at the cap. */
        Expr* prev = r_raw;  long Wp = raised ? work.bits : p;  long Lp = cancelled;
        long Ap = absorbed;
        Expr* prev2 = NULL;  long Wp2 = 0;
        Expr* first = expr_copy(r_raw);
        while (prev) {
            long W2 = Wp + 64;
            if (Lp < Wp) {                       /* a measured loss: go where it is covered */
                if (p + amp + Lp + 32 > W2) W2 = p + amp + Lp + 32;
            } else {
                /* cancelled to zero: double the extra precision -- or, when
                 * an addition of the pass absorbed an operand whole, go where
                 * that operand is seen (Exp[x] (Exp[1/x - Exp[-x]] -
                 * Exp[1/x]) at x = 1000: the E^-1000 that 1/x absorbed is the
                 * whole answer; one step to ~1500 bits instead of five
                 * doublings) */
                if (Wp - p > 64) W2 = Wp + (Wp - p);
                if (Ap > 0 && p + amp + Ap + 64 > W2) W2 = p + amp + Ap + 64;
            }
            /* the pass at the cap is compared with one 64 bits below it */
            if (W2 > Wmax - 64 && Wp < Wmax - 64) W2 = Wmax - 64;
            if (W2 > Wmax) W2 = Wmax;
            if (W2 <= Wp) break;
            NumericSpec w = work; w.bits = W2;   /* work's mode: CAP keeps inexact leaves as they are */
            g_numeric_cancel_bits = 0;
            g_numeric_cancel_left = LONG_MAX;
            g_numeric_absorb_bits = 0;
            g_numeric_reevaluations++;
            Expr* cur = numericalize_pass(e, w);
            long Lc = g_numeric_cancel_bits;
            long Kc = g_numeric_cancel_left;
            long Ac = g_numeric_absorb_bits;
            if (Lc > cancelled) cancelled = Lc;
            if (Ac > absorbed) absorbed = Ac;
            if (!cur) break;
            /* A higher precision that comes back WORSE (a kernel that fails
             * there: N[BarnesG[37/10], 80] was Overflow[]) ends the climb on
             * the last good result instead of chasing it to the cap. */
            if ((expr_has_nonfinite_real(cur) || expr_has_overflow_or_indeterminate(cur))
                && !(expr_has_nonfinite_real(prev) || expr_has_overflow_or_indeterminate(prev))) {
                expr_free(cur);
                break;
            }
            /* compare with the newest result at least 32 bits lower */
            const Expr* ref = prev;  long Wr = Wp;
            if (W2 - Wp < 32 && prev2) { ref = prev2; Wr = Wp2; }
            long gap = W2 - Wr;
            long good = adaptive_good_bits(expr_agree_bits(cur, ref), gap);
            if (good >= p && Kc >= p + amp) {     /* certified: p good digits */
                expr_free(prev);
                prev = cur;
                break;
            }
            if (W2 >= Wmax) {                    /* at the cap */
                bool vanished = false, unsettled = false;
                Expr* z = expr_cap_verdict(cur, ref, first, gap, p, Kc <= 16,
                                           &vanished, &unsettled);
                /* Overflow[] / Indeterminate at the cap is a division by a
                 * difference that cancelled to zero at every precision tried
                 * (f[g[x]]/g[f[x]] at x deep in the cancellation): it is not
                 * an answer either. */
                if (expr_has_nonfinite_real(cur) || expr_has_overflow_or_indeterminate(cur))
                    unsettled = true;
                expr_free(cur);
                expr_free(prev);
                prev = z;
                if (unsettled || (vanished && from_n)) {
                    char* str = expr_to_string((Expr*)e);
                    mth_message("N", "meprec",
                        "Internal precision limit $MaxExtraPrecision = %g reached while evaluating %s.",
                        mep_digits, str ? str : "");
                    free(str);
                }
                break;
            }
            expr_free(prev2);
            prev2 = prev;  Wp2 = Wp;
            prev = cur;    Wp = W2;  Lp = Lc;  Ap = Ac;
        }
        expr_free(prev2);
        expr_free(first);
        if (prev) {
            expr_free(r);
            r = numeric_round_result(prev, spec);
        }
    }
    /* an enclosing numericalization sees the worst cancellation of any pass */
    g_numeric_cancel_bits = saved_cancel > cancelled ? saved_cancel : cancelled;
    if (saved_left < g_numeric_cancel_left) g_numeric_cancel_left = saved_left;
    g_numeric_absorb_bits = saved_absorb > absorbed ? saved_absorb : absorbed;

    /* A machine double has a 53-bit mantissa but only a ~1e308 exponent
     * range, so Exp[1000] flushes to inf and Exp[-1000] to zero even though
     * a machine *number* (mantissa + arbitrary exponent, as N[1001!]
     * already returns) represents both fine. When that happens, and the
     * input itself held no infinity or zero to explain it, recompute once
     * in MPFR and keep the answer only if it is genuinely better. Fires
     * only on a degenerate result, so the normal path pays nothing. */
    if (spec.mode == NUMERIC_MODE_MACHINE && s.has_call && !s.has_infinity) {
        /* A raw Inf/NaN, or the Overflow[]/Indeterminate the net rewrote one
         * into, both signal a machine overflow that MPFR may resolve. An MPFR
         * leaf in a machine-mode result likewise means a head already promoted
         * an overflow internally (Exp[1000.] -> 1.97*10^434): recompute from the
         * exact form at working precision, since a head that promoted a *rounded*
         * machine argument (E^1000 as (2.718...53-bit)^1000) is less accurate
         * than deriving E to full precision here. */
        bool inf  = expr_has_nonfinite_real(r) || expr_has_overflow_or_indeterminate(r)
                 || expr_has_extended_range_real(r);
        bool zero = !inf && !s.has_zero && expr_has_zero_real(r);
        if (inf || zero) {
            NumericSpec retry = spec;
            retry.mode = NUMERIC_MODE_MPFR;
            retry.bits = (raised && work.bits > (long)DBL_MANT_DIG + NUMERIC_GUARD_BITS)
                       ? work.bits : (long)DBL_MANT_DIG + NUMERIC_GUARD_BITS;
            Expr* alt = numeric_round_result(numericalize_rec(e, retry), spec);
            bool better = alt && !expr_has_nonfinite_real(alt)
                       && !expr_has_overflow_or_indeterminate(alt)
                       && (inf || !expr_has_zero_real(alt));
            if (better) {
                expr_free(r);
                return alt;
            }
            expr_free(alt);
        }
    }
    return r;
#endif
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

bool arg_is_inexact(const Expr* e) {
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
    NumericSpec spec = numeric_machine_spec();   /* numericalize_contagion below */
#ifdef USE_MPFR
    if (!any_machine && min_mpfr_prec > 0) {
        spec.mode = NUMERIC_MODE_MPFR;
        spec.bits = min_mpfr_prec;
    }
#endif

    /* Only the numeric arguments take the inexact precision, as in the
     * reference: 1.5 + Pi + a is 4.64159 + a, but 1.5 - a stays 1.5 - a (not
     * 1.5 - 1.0 a), 1.5 + a - a is 1.5, and 1.5 (2 + a) and 1.5 + Sqrt[2] a
     * keep their exact parts. */
    for (size_t i = 0; i < n; i++) {
        out[i] = expr_is_numeric_quantity(args[i]) ? numericalize_contagion(args[i], spec)
                                                   : expr_copy(args[i]);
    }
    return true;
}

/* ------------------------------------------------------------------------
 *  Mixed exact/inexact arguments of a numeric builtin
 *
 *  Zeta[2., Sqrt[2]] is a machine number in the reference, exactly as
 *  Sin[2.] is: a numeric builtin given an inexact argument numericalizes its
 *  exact numeric arguments to the inexact one's precision. Plus / Times do
 *  that in their own contagion pass; every other head used to stay
 *  unevaluated, because each one's numeric path tests "all arguments are
 *  numbers". numeric_mixed_args_call() supplies the rule once, for every
 *  NumericFunction builtin, from the evaluator (step 5 of evaluate_step).
 * ---------------------------------------------------------------------- */

/* Rational[p, q], or Complex[a, b] with integer / rational parts. */
static bool mixed_is_rational_atom(const Expr* e) {
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return true;
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL) return false;
    const char* h = e->data.function.head->data.symbol.name;
    if (h == SYM_Rational) return true;
    if (h != SYM_Complex || e->data.function.arg_count != 2) return false;
    return mixed_is_rational_atom(e->data.function.args[0]) &&
           mixed_is_rational_atom(e->data.function.args[1]);
}

/* Classify the leaves of one argument: a number (inexact or exact), a numeric
 * quantity (Sqrt[2], Pi, Gamma[1/3]: exact), or a List of those. Returns false
 * on anything else (a symbol, a string, a rule), which ends the attempt. */
static bool mixed_arg_scan(const Expr* e, bool* any_inexact, bool* any_exact,
                           bool* any_machine, long* min_bits, int depth) {
    if (e->type == EXPR_REAL) { *any_inexact = true; *any_machine = true; return true; }
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) {
        long p = (long)mpfr_get_prec(e->data.mpfr);
        if (*min_bits == 0 || p < *min_bits) *min_bits = p;
        *any_inexact = true;
        return true;
    }
#endif
    /* Integers, rationals and exact Gaussian rationals are already what every
     * numeric path accepts beside a real (and often a structural index:
     * BesselJ's order, ProductLog's branch), so they neither trigger the
     * rewrite nor change under it. Only an irrational numeric quantity does. */
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT) return true;
    if (e->type != EXPR_FUNCTION) {
        if (e->type == EXPR_SYMBOL && expr_is_numeric_quantity(e)) { *any_exact = true; return true; }
        return false;
    }
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im) || e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.head->data.symbol.name != SYM_List) {
        if (arg_is_inexact(e)) {
            if (arg_has_machine_real(e)) *any_machine = true;
#ifdef USE_MPFR
            long p = numeric_min_inexact_bits(e);
            if (p > 0 && !arg_has_machine_real(e) && (*min_bits == 0 || p < *min_bits)) *min_bits = p;
#endif
            /* A compound holding an inexact leaf (Sqrt[2.]-style residue an
             * earlier pass left) is not a number; leave it alone. */
            if (!is_complex((Expr*)e, &re, &im)) return false;
            *any_inexact = true;
            return true;
        }
        if (mixed_is_rational_atom(e)) return true;
        if (!expr_is_numeric_quantity(e)) return false;
        *any_exact = true;
        return true;
    }
    /* A List argument (HypergeometricPFQ's parameter lists): its leaves count. */
    if (depth > 4) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!mixed_arg_scan(e->data.function.args[i], any_inexact, any_exact,
                            any_machine, min_bits, depth + 1))
            return false;
    return true;
}

static bool mixed_is_list(const Expr* e) {
    return e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL &&
           e->data.function.head->data.symbol.name == SYM_List;
}

/* True iff every leaf of `e` is an inexact number or a rational atom: nothing
 * left for the rewrite to convert (and the retry's postcondition). */
static bool mixed_arg_settled(const Expr* e) {
    if (e->type == EXPR_REAL || mixed_is_rational_atom(e)) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    if (e->type != EXPR_FUNCTION) return false;
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im)) return arg_is_inexact(e);
    if (!mixed_is_list(e)) return false;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (!mixed_arg_settled(e->data.function.args[i])) return false;
    return true;
}

/* N of one argument, leaf by leaf through Lists so a rational parameter
 * beside the irrational one keeps its exact form. */
static Expr* mixed_numericalize(const Expr* e, NumericSpec spec) {
    if (mixed_arg_settled(e)) return expr_copy((Expr*)e);
    if (!mixed_is_list(e)) return numericalize(e, spec);
    size_t n = e->data.function.arg_count;
    Expr** a = (Expr**)malloc((n ? n : 1) * sizeof(Expr*));
    if (!a) return NULL;
    for (size_t i = 0; i < n; i++) {
        a[i] = mixed_numericalize(e->data.function.args[i], spec);
        if (!a[i]) {
            for (size_t j = 0; j < i; j++) expr_free(a[j]);
            free(a);
            return NULL;
        }
    }
    Expr* out = expr_new_function(expr_copy(e->data.function.head), a, n);
    free(a);
    return out;
}

Expr* numeric_mixed_args_call(const Expr* call, uint32_t attrs) {
    if (!call || call->type != EXPR_FUNCTION) return NULL;
    size_t n = call->data.function.arg_count;
    if (n < 2) return NULL;
    /* Cheap gate first: this runs on every declined NumericFunction call
     * (Power[x, 2], ...), and an inexact argument can only be a top-level
     * Real / MPFR, a Complex, or sit inside a List. */
    bool candidate = false;
    for (size_t i = 0; i < n && !candidate; i++) {
        const Expr* a = call->data.function.args[i];
        if (a->type == EXPR_REAL) candidate = true;
#ifdef USE_MPFR
        else if (a->type == EXPR_MPFR) candidate = true;
#endif
        else if (a->type == EXPR_FUNCTION && a->data.function.head->type == EXPR_SYMBOL &&
                 (a->data.function.head->data.symbol.name == SYM_Complex ||
                  a->data.function.head->data.symbol.name == SYM_List))
            candidate = true;
    }
    if (!candidate) return NULL;
    bool any_inexact = false, any_exact = false, any_machine = false;
    long min_bits = 0;
    for (size_t i = 0; i < n; i++) {
        if ((i == 0 && (attrs & ATTR_NHOLDFIRST)) || (i > 0 && (attrs & ATTR_NHOLDREST)))
            continue;
        if (!mixed_arg_scan(call->data.function.args[i], &any_inexact, &any_exact,
                            &any_machine, &min_bits, 0))
            return NULL;
    }
    if (!any_inexact || !any_exact) return NULL;

    /* The lowest precision governs, and machine precision is the floor. */
    NumericSpec spec = numeric_machine_spec();
#ifdef USE_MPFR
    if (!any_machine && min_bits > 0) {
        spec.mode = NUMERIC_MODE_MPFR_CAP;
        spec.bits = min_bits;
    }
#else
    (void)min_bits;
#endif
    Expr** args = (Expr**)malloc(n * sizeof(Expr*));
    if (!args) return NULL;
    bool ok = true;
    for (size_t i = 0; i < n; i++) {
        const Expr* a = call->data.function.args[i];
        bool held = (i == 0 && (attrs & ATTR_NHOLDFIRST)) || (i > 0 && (attrs & ATTR_NHOLDREST));
        if (held || mixed_arg_settled(a)) { args[i] = expr_copy((Expr*)a); continue; }
        args[i] = mixed_numericalize(a, spec);
        if (!args[i] || !mixed_arg_settled(args[i])) {
            if (!args[i]) args[i] = expr_copy((Expr*)a);
            ok = false;
        }
    }
    if (!ok) {
        for (size_t i = 0; i < n; i++) expr_free(args[i]);
        free(args);
        return NULL;
    }
    Expr* out = expr_new_function(expr_copy(call->data.function.head), args, n);
    free(args);
    return out;
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
    return true;
#else
    /* Phase 1 fallback: emit a one-shot warning, then use machine. */
    static bool warned = false;
    if (!warned) {
        mth_message("N", "prec", "arbitrary precision unavailable in this build (USE_MPFR=0); using machine precision.");
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
        const Expr* p = res->data.function.args[1];
        /* N[expr, {p, a}]: a precision goal and an accuracy goal. The
         * precision goal decides (N[Pi, {MachinePrecision, MachinePrecision}]
         * is machine Pi, N[Pi, {20, 20}] has 20 digits); an Infinity
         * precision goal leaves the accuracy goal, used as the digit count. */
        if (p->type == EXPR_FUNCTION && p->data.function.head->type == EXPR_SYMBOL &&
            p->data.function.head->data.symbol.name == SYM_List && p->data.function.arg_count == 2) {
            const Expr* pg = p->data.function.args[0];
            if (pg->type == EXPR_SYMBOL && pg->data.symbol.name == SYM_Infinity)
                pg = p->data.function.args[1];
            p = pg;
        }
        if (!parse_precision_arg(p, &spec)) {
            return NULL;  /* non-numeric precision → remain unevaluated */
        }
    }
    /* One-argument N[expr] keeps the machine spec: Mathematica's N[expr]
     * targets machine precision even when expr is already an approximate
     * arbitrary-precision number, so N[N[Pi, 100]] comes back at machine
     * precision (see numericalize's EXPR_MPFR branch). Use the two-argument
     * N[expr, p] to request a specific number of digits. */

    /* Ownership: Mathilda evaluator frees `res` after a non-NULL return
     * (see eval.c around the builtin dispatch site), so we must NOT free
     * it ourselves. Returning NULL means "keep `res` unevaluated" — the
     * evaluator retains ownership in that path. */
    g_numeric_from_n = true;
    Expr* out = numericalize(res->data.function.args[0], spec);
    g_numeric_from_n = false;
    return out;
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
    /* Not Listable (the reference's N is only Protected): numericalize
     * descends into lists itself, and N[x, {p, a}] is a precision/accuracy
     * pair, not two precisions. */
    symtab_get_def("N")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("N",
        "N[expr]\n"
        "\tGives a machine-precision numerical approximation of expr.\n"
        "N[expr, n]\n"
        "\tGives a numerical approximation to n decimal digits. Requires\n"
        "\ta USE_MPFR build; without it, a warning is emitted and machine\n"
        "\tprecision is used. When the evaluation cancels, the working\n"
        "\tprecision rises until n digits are certified, up to n +\n"
        "\t$MaxExtraPrecision digits (N::meprec beyond that).\n"
        "N[expr, {p, a}]\n"
        "\tPrecision goal p and accuracy goal a; p is used (a when p is Infinity).\n");

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
