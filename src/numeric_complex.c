/* Mathilda — MPFR-complex transcendental helpers.
 *
 * See numeric_complex.h for the contract. This file implements the
 * concrete ops, the small dispatch helper used by per-builtin wirings,
 * and a single private MPFR-complex multiply / sqrt used by the
 * inverse-trig / inverse-hyperbolic compositions.
 *
 * Memory contract: every mpfr_inits2/mpfr_init2 is paired with a
 * mpfr_clear / mpfr_clears at the bottom of the function. No allocation
 * survives a function return.
 */

#ifdef USE_MPFR

#include "numeric_complex.h"
#include "numeric.h"
#include "arithmetic.h"   /* make_complex */
#include "expr.h"
#include <mpfr.h>
#include <stdbool.h>

/* --------------------------------------------------------------------
 *  Construction
 * ------------------------------------------------------------------ */

Expr* numeric_mpfr_make_complex(const mpfr_t out_re, const mpfr_t out_im) {
    if (mpfr_zero_p(out_im)) {
        return expr_new_mpfr_copy(out_re);
    }
    Expr* re = expr_new_mpfr_copy(out_re);
    Expr* im = expr_new_mpfr_copy(out_im);
    return make_complex(re, im);
}

/* --------------------------------------------------------------------
 *  Private MPFR-complex algebra used by the unary ops below
 * ------------------------------------------------------------------ */

/* Complex sqrt via the cancellation-free identity. Branch: principal,
 * matching Mathematica (Re >= 0, or Re == 0 and Im >= 0 on the cut). */
static void mc_sqrt(mpfr_t out_re, mpfr_t out_im,
                    const mpfr_t a_re, const mpfr_t a_im) {
    mpfr_prec_t p = mpfr_get_prec(out_re);
    if (mpfr_zero_p(a_im)) {
        if (mpfr_sgn(a_re) >= 0) {
            mpfr_sqrt(out_re, a_re, MPFR_RNDN);
            mpfr_set_zero(out_im, +1);
        } else {
            mpfr_t neg;
            mpfr_init2(neg, p);
            mpfr_neg(neg, a_re, MPFR_RNDN);
            mpfr_sqrt(out_im, neg, MPFR_RNDN);
            mpfr_set_zero(out_re, +1);
            mpfr_clear(neg);
        }
        return;
    }
    mpfr_t mag, abs_a, half_sum, w, two_w;
    mpfr_inits2(p, mag, abs_a, half_sum, w, two_w, (mpfr_ptr)NULL);
    mpfr_hypot(mag, a_re, a_im, MPFR_RNDN);
    mpfr_abs(abs_a, a_re, MPFR_RNDN);
    mpfr_add(half_sum, mag, abs_a, MPFR_RNDN);
    mpfr_div_2ui(half_sum, half_sum, 1, MPFR_RNDN);  /* (|z| + |a|) / 2 */
    mpfr_sqrt(w, half_sum, MPFR_RNDN);
    if (mpfr_sgn(a_re) >= 0) {
        mpfr_set(out_re, w, MPFR_RNDN);
        mpfr_mul_2ui(two_w, w, 1, MPFR_RNDN);
        mpfr_div(out_im, a_im, two_w, MPFR_RNDN);
    } else {
        if (mpfr_sgn(a_im) >= 0) mpfr_set(out_im, w, MPFR_RNDN);
        else                     mpfr_neg(out_im, w, MPFR_RNDN);
        mpfr_mul_2ui(two_w, out_im, 1, MPFR_RNDN);
        mpfr_div(out_re, a_im, two_w, MPFR_RNDN);
    }
    mpfr_clears(mag, abs_a, half_sum, w, two_w, (mpfr_ptr)NULL);
}

void mpfr_complex_div(mpfr_t out_re, mpfr_t out_im,
                      const mpfr_t a_re, const mpfr_t a_im,
                      const mpfr_t b_re, const mpfr_t b_im) {
    mpfr_prec_t p = mpfr_get_prec(out_re);
    mpfr_t abs_br, abs_bi, r, den, t1, t2, re_scratch;
    mpfr_inits2(p, abs_br, abs_bi, r, den, t1, t2, re_scratch, (mpfr_ptr)NULL);
    mpfr_abs(abs_br, b_re, MPFR_RNDN);
    mpfr_abs(abs_bi, b_im, MPFR_RNDN);
    /* Smith's algorithm: divide by the larger magnitude axis to avoid
     * over/underflow when |b_re| and |b_im| differ wildly. */
    if (mpfr_cmp(abs_br, abs_bi) >= 0) {
        mpfr_div(r, b_im, b_re, MPFR_RNDN);
        mpfr_mul(t1, r, b_im, MPFR_RNDN);
        mpfr_add(den, b_re, t1, MPFR_RNDN);
        mpfr_mul(t1, r, a_im, MPFR_RNDN);
        mpfr_add(t2, a_re, t1, MPFR_RNDN);
        mpfr_div(re_scratch, t2, den, MPFR_RNDN);
        mpfr_mul(t1, r, a_re, MPFR_RNDN);
        mpfr_sub(t2, a_im, t1, MPFR_RNDN);
        mpfr_div(out_im, t2, den, MPFR_RNDN);
    } else {
        mpfr_div(r, b_re, b_im, MPFR_RNDN);
        mpfr_mul(t1, r, b_re, MPFR_RNDN);
        mpfr_add(den, t1, b_im, MPFR_RNDN);
        mpfr_mul(t1, r, a_re, MPFR_RNDN);
        mpfr_add(t2, t1, a_im, MPFR_RNDN);
        mpfr_div(re_scratch, t2, den, MPFR_RNDN);
        mpfr_mul(t1, r, a_im, MPFR_RNDN);
        mpfr_sub(t2, t1, a_re, MPFR_RNDN);
        mpfr_div(out_im, t2, den, MPFR_RNDN);
    }
    mpfr_set(out_re, re_scratch, MPFR_RNDN);
    mpfr_clears(abs_br, abs_bi, r, den, t1, t2, re_scratch, (mpfr_ptr)NULL);
}

/* --------------------------------------------------------------------
 *  Forward transcendentals
 *
 * Every function here is cancellation-free in floating point -- each sum
 * it forms has terms of one sign -- so the result is accurate to a few
 * ulps at any precision, near zeros, poles and branch points included.  The
 * naive log-form inverses this file used to carry (asin = -i log(iz +
 * sqrt(1 - z^2)) and friends) cancelled catastrophically on half the plane:
 * ArcSin[SetPrecision[-10^30, 40]] was Overflow[] (iz + sqrt(1 - z^2)
 * rounded to exactly 0), ArcSin[SetPrecision[-10^30, 60]] six digits,
 * ArcTan[SetPrecision[10^-30 I, 30]] six digits, Log[1 + 10^-20 I] at 30
 * digits had Re = 0.
 * ------------------------------------------------------------------ */

/* Guard bits for the internal working precision of the compositions below;
 * the result is rounded once, into the caller's output precision. */
#define NC_GUARD 32

/* r = x * y, where a zero factor gives an exact zero even when the other
 * factor is infinite.  sin(0 + i 10^30) = sin(0) cosh(10^30) + ... has an
 * exactly-zero real part; 0 * inf = NaN made it Indeterminate. */
static void nc_mul0(mpfr_t r, const mpfr_t x, const mpfr_t y) {
    if (mpfr_zero_p(x) || mpfr_zero_p(y)) {
        bool neg = (mpfr_signbit(x) != 0) != (mpfr_signbit(y) != 0);
        mpfr_set_zero(r, neg ? -1 : 1);
    } else {
        mpfr_mul(r, x, y, MPFR_RNDN);
    }
}

void mpfr_complex_exp(mpfr_t out_re, mpfr_t out_im,
                      const mpfr_t a, const mpfr_t b) {
    mpfr_prec_t p = mpfr_get_prec(out_re);
    mpfr_t ea, cb, sb;
    mpfr_inits2(p, ea, cb, sb, (mpfr_ptr)NULL);
    mpfr_exp(ea, a, MPFR_RNDN);
    mpfr_sin_cos(sb, cb, b, MPFR_RNDN);
    nc_mul0(out_re, ea, cb);
    nc_mul0(out_im, ea, sb);
    mpfr_clears(ea, cb, sb, (mpfr_ptr)NULL);
}

/* log(a + bi) = log|z| + i atan2(b, a).  Near the unit circle log|z| is the
 * difference of two nearly equal numbers, so there it is log1p(|z|^2 - 1)/2
 * with |z|^2 - 1 = (a - 1)(a + 1) + b^2 formed exactly and rounded once
 * (a - 1 is exact for a in [1/2, 2], the products are exact at 2p bits).
 * The imaginary part keeps atan2's IEEE signed-zero reading. */
void mpfr_complex_log(mpfr_t out_re, mpfr_t out_im,
                      const mpfr_t a, const mpfr_t b) {
    mpfr_prec_t p = mpfr_get_prec(out_re);
    mpfr_prec_t pa = mpfr_get_prec(a), pb = mpfr_get_prec(b);
    mpfr_prec_t pin = pa > pb ? pa : pb;
    if (pin < p) pin = p;
    mpfr_t r, im;
    mpfr_init2(r, p + NC_GUARD);
    mpfr_init2(im, p);
    mpfr_atan2(im, b, a, MPFR_RNDN);      /* before out_re may alias a */
    mpfr_hypot(r, a, b, MPFR_RNDN);
    if (mpfr_number_p(r) && mpfr_cmp_d(r, 0.5) >= 0 && mpfr_cmp_ui(r, 2) <= 0) {
        mpfr_t am1, ap1, prod, b2;
        mpfr_inits2(2 * pin + 8, am1, ap1, prod, b2, (mpfr_ptr)NULL);
        mpfr_sub_ui(am1, a, 1, MPFR_RNDN);
        mpfr_add_ui(ap1, a, 1, MPFR_RNDN);
        mpfr_mul(prod, am1, ap1, MPFR_RNDN);
        mpfr_sqr(b2, b, MPFR_RNDN);
        mpfr_add(r, prod, b2, MPFR_RNDN);   /* |z|^2 - 1, one rounding */
        mpfr_log1p(r, r, MPFR_RNDN);
        mpfr_div_2ui(out_re, r, 1, MPFR_RNDN);
        mpfr_clears(am1, ap1, prod, b2, (mpfr_ptr)NULL);
    } else {
        mpfr_log(out_re, r, MPFR_RNDN);
    }
    mpfr_set(out_im, im, MPFR_RNDN);
    mpfr_clear(r);
    mpfr_clear(im);
}

void mpfr_complex_sin(mpfr_t out_re, mpfr_t out_im,
                      const mpfr_t a, const mpfr_t b) {
    mpfr_prec_t p = mpfr_get_prec(out_re);
    mpfr_t sa, ca, sb, cb;
    mpfr_inits2(p, sa, ca, sb, cb, (mpfr_ptr)NULL);
    mpfr_sin_cos(sa, ca, a, MPFR_RNDN);
    mpfr_sinh_cosh(sb, cb, b, MPFR_RNDN);
    nc_mul0(out_re, sa, cb);
    nc_mul0(out_im, ca, sb);
    mpfr_clears(sa, ca, sb, cb, (mpfr_ptr)NULL);
}

void mpfr_complex_cos(mpfr_t out_re, mpfr_t out_im,
                      const mpfr_t a, const mpfr_t b) {
    mpfr_prec_t p = mpfr_get_prec(out_re);
    mpfr_t sa, ca, sb, cb;
    mpfr_inits2(p, sa, ca, sb, cb, (mpfr_ptr)NULL);
    mpfr_sin_cos(sa, ca, a, MPFR_RNDN);
    mpfr_sinh_cosh(sb, cb, b, MPFR_RNDN);
    nc_mul0(out_re, ca, cb);
    nc_mul0(out_im, sa, sb);
    mpfr_neg(out_im, out_im, MPFR_RNDN);
    mpfr_clears(sa, ca, sb, cb, (mpfr_ptr)NULL);
}

/* tan, cot, tanh, coth, sec, csc, sech, csch of a + bi, from one
 * cancellation-free, overflow-free form.  With c the circular argument and h
 * the hyperbolic one (c = a, h = b for the circular functions; c = b, h = a
 * for the hyperbolic ones), s = sech h, t = tanh h:
 *
 *   tan  = (sin c cos c s^2 + i t) / (cos^2 c s^2 + t^2)
 *   cot  = (sin c cos c s^2 - i t) / (sin^2 c s^2 + t^2)
 *   tanh = (t + i sin c cos c s^2) / (cos^2 c s^2 + t^2)
 *   coth = (t - i sin c cos c s^2) / (sin^2 c s^2 + t^2)
 *   sec  = (cos c s + i sin c t s) / (cos^2 c s^2 + t^2)
 *   csc  = (sin c s - i cos c t s) / (sin^2 c s^2 + t^2)
 *   sech = (cos c s - i sin c t s) / (cos^2 c s^2 + t^2)
 *   csch = (cos c t s - i sin c s) / (sin^2 c s^2 + t^2)
 *
 * They are the textbook (sin a cos a + i sinh b cosh b) / (cos^2 a + sinh^2 b)
 * and its relatives divided through by cosh^2 h.  The denominator is a sum of
 * two squares, so it never cancels next to a pole (the double-angle spelling
 * cos 2a + cosh 2b cancels to nothing at a = -Pi/2 + d, b tiny), and s, t
 * stay in [0, 1]: Tan[SetPrecision[10^30 I, 30]] is I, where cosh^2 b
 * overflowed and inf/inf made it Indeterminate. */
typedef enum { NC_TAN, NC_COT, NC_SEC, NC_CSC } NcTanKind;

static void complex_tan_family(mpfr_t out_re, mpfr_t out_im,
                               const mpfr_t a, const mpfr_t b,
                               bool hyp, NcTanKind kind) {
    mpfr_prec_t p = mpfr_get_prec(out_re) + 8;
    mpfr_t sc, cc, s, t, s2, re, im, d;
    mpfr_inits2(p, sc, cc, s, t, s2, re, im, d, (mpfr_ptr)NULL);
    mpfr_sin_cos(sc, cc, hyp ? b : a, MPFR_RNDN);
    mpfr_sech(s, hyp ? a : b, MPFR_RNDN);
    mpfr_tanh(t, hyp ? a : b, MPFR_RNDN);
    mpfr_sqr(s2, s, MPFR_RNDN);
    /* denominator: (cos^2 c or sin^2 c) s^2 + t^2 */
    bool use_sin = (kind == NC_COT || kind == NC_CSC);
    mpfr_sqr(d, use_sin ? sc : cc, MPFR_RNDN);
    nc_mul0(d, d, s2);
    mpfr_sqr(re, t, MPFR_RNDN);
    mpfr_add(d, d, re, MPFR_RNDN);
    switch (kind) {
    case NC_TAN:
    case NC_COT:
        mpfr_mul(re, sc, cc, MPFR_RNDN);
        nc_mul0(re, re, s2);                      /* sin c cos c s^2 */
        mpfr_set(im, t, MPFR_RNDN);
        if (kind == NC_COT) mpfr_neg(im, im, MPFR_RNDN);
        if (hyp) mpfr_swap(re, im);               /* tanh / coth */
        if (hyp && kind == NC_COT) { mpfr_neg(re, re, MPFR_RNDN); mpfr_neg(im, im, MPFR_RNDN); }
        break;
    case NC_SEC:
        nc_mul0(re, cc, s);                       /* cos c s */
        mpfr_mul(im, sc, t, MPFR_RNDN);
        nc_mul0(im, im, s);                       /* sin c t s */
        if (hyp) mpfr_neg(im, im, MPFR_RNDN);     /* sech */
        break;
    case NC_CSC:
        if (!hyp) {
            nc_mul0(re, sc, s);                   /* sin c s */
            mpfr_mul(im, cc, t, MPFR_RNDN);
            nc_mul0(im, im, s);
            mpfr_neg(im, im, MPFR_RNDN);          /* -cos c t s */
        } else {
            mpfr_mul(re, cc, t, MPFR_RNDN);
            nc_mul0(re, re, s);                   /* cos c t s */
            nc_mul0(im, sc, s);
            mpfr_neg(im, im, MPFR_RNDN);          /* -sin c s */
        }
        break;
    }
    mpfr_div(out_re, re, d, MPFR_RNDN);
    mpfr_div(out_im, im, d, MPFR_RNDN);
    mpfr_clears(sc, cc, s, t, s2, re, im, d, (mpfr_ptr)NULL);
}

void mpfr_complex_tan(mpfr_t out_re, mpfr_t out_im,
                      const mpfr_t a, const mpfr_t b) {
    complex_tan_family(out_re, out_im, a, b, false, NC_TAN);
}

void mpfr_complex_sinh(mpfr_t out_re, mpfr_t out_im,
                       const mpfr_t a, const mpfr_t b) {
    /* sinh(a+bi) = sinh(a) cos(b) + i cosh(a) sin(b) */
    mpfr_prec_t p = mpfr_get_prec(out_re);
    mpfr_t sa, ca, sb, cb;
    mpfr_inits2(p, sa, ca, sb, cb, (mpfr_ptr)NULL);
    mpfr_sinh_cosh(sa, ca, a, MPFR_RNDN);
    mpfr_sin_cos(sb, cb, b, MPFR_RNDN);
    nc_mul0(out_re, sa, cb);
    nc_mul0(out_im, ca, sb);
    mpfr_clears(sa, ca, sb, cb, (mpfr_ptr)NULL);
}

void mpfr_complex_cosh(mpfr_t out_re, mpfr_t out_im,
                       const mpfr_t a, const mpfr_t b) {
    /* cosh(a+bi) = cosh(a) cos(b) + i sinh(a) sin(b) */
    mpfr_prec_t p = mpfr_get_prec(out_re);
    mpfr_t sa, ca, sb, cb;
    mpfr_inits2(p, sa, ca, sb, cb, (mpfr_ptr)NULL);
    mpfr_sinh_cosh(sa, ca, a, MPFR_RNDN);
    mpfr_sin_cos(sb, cb, b, MPFR_RNDN);
    nc_mul0(out_re, ca, cb);
    nc_mul0(out_im, sa, sb);
    mpfr_clears(sa, ca, sb, cb, (mpfr_ptr)NULL);
}

void mpfr_complex_tanh(mpfr_t out_re, mpfr_t out_im,
                       const mpfr_t a, const mpfr_t b) {
    complex_tan_family(out_re, out_im, a, b, true, NC_TAN);
}

/* --------------------------------------------------------------------
 *  Inverse transcendentals
 *
 * The values are those of the principal-branch definitions
 *   asin z  = -i log(i z + sqrt(1 - z^2))       acos z = Pi/2 - asin z
 *   atan z  = (i/2)(log(1 - i z) - log(1 + i z))
 *   asinh z = log(z + sqrt(z^2 + 1))
 *   acosh z = log(z + sqrt(z - 1) sqrt(z + 1))
 *   atanh z = (log(1 + z) - log(1 - z))/2
 * (and the reciprocal functions f(1/z)) with principal Log and Sqrt and a
 * zero part read as an exact zero -- the convention the exact evaluator
 * uses, so N[f[p], n] and f[SetPrecision[p, n]] agree ON the cuts too.  They
 * are computed through cancellation-free forms (Kahan, "Branch cuts for
 * complex elementary functions", 1987):
 *
 *   asin z  = atan2(Re z, Re(s1 s2)) + i asinh(Im(conj(s1) s2))
 *   acos z  = 2 atan2(Re s1, Re s2)  + i asinh(Im(conj(s2) s1))
 *   acosh z = asinh(Re(conj(r1) r2))  + 2 i atan2(Im r1, Re r2)
 *   atanh z = log(((1+a)^2 + b^2) / ((1-a)^2 + b^2))/4
 *             + i atan2(2b, (1 - a)(1 + a) - b^2)/2
 *
 * with s1 = sqrt(1 - z), s2 = sqrt(1 + z), r1 = sqrt(z - 1), r2 = sqrt(z + 1),
 * and asinh / atan / acot / acoth reduced to these by z -> i z.  In every
 * real sum the two terms share a sign, so nothing cancels.  The sqrt is
 * mc_sqrt, which reads a zero imaginary part as exact (a negative real has
 * the root +i sqrt|x|): that, and nothing else, is what puts each on-cut
 * value on the side the definitions above give.  atanh's atan2 would read
 * the sign of a zero b, so the real cuts are set explicitly.
 *
 * The reciprocal functions take 1 -+ w and w -+ 1 for w = 1/z as
 * (z -+ 1)/z, each one complex division, not 1 -+ (1/z): ArcSec at
 * 1 + 10^-25 loses no digits to the rounding of 1/z.
 * ------------------------------------------------------------------ */

/* asin from x = Re w, om = 1 - w, op = 1 + w; outputs at their own
 * precision, internals at wp. */
static void nc_asin_core(mpfr_t out_re, mpfr_t out_im, const mpfr_t x,
                         const mpfr_t om_re, const mpfr_t om_im,
                         const mpfr_t op_re, const mpfr_t op_im, mpfr_prec_t wp) {
    mpfr_t s1r, s1i, s2r, s2i, u, v;
    mpfr_inits2(wp, s1r, s1i, s2r, s2i, u, v, (mpfr_ptr)NULL);
    mc_sqrt(s1r, s1i, om_re, om_im);
    mc_sqrt(s2r, s2i, op_re, op_im);
    /* Re(s1 s2) = s1r s2r - s1i s2i  (s1i s2i <= 0) */
    mpfr_mul(u, s1r, s2r, MPFR_RNDN);
    mpfr_mul(v, s1i, s2i, MPFR_RNDN);
    mpfr_sub(u, u, v, MPFR_RNDN);
    mpfr_atan2(out_re, x, u, MPFR_RNDN);
    /* Im(conj(s1) s2) = s1r s2i - s1i s2r  (both terms of b's sign) */
    mpfr_mul(u, s1r, s2i, MPFR_RNDN);
    mpfr_mul(v, s1i, s2r, MPFR_RNDN);
    mpfr_sub(u, u, v, MPFR_RNDN);
    mpfr_asinh(out_im, u, MPFR_RNDN);
    mpfr_clears(s1r, s1i, s2r, s2i, u, v, (mpfr_ptr)NULL);
}

static void nc_acos_core(mpfr_t out_re, mpfr_t out_im,
                         const mpfr_t om_re, const mpfr_t om_im,
                         const mpfr_t op_re, const mpfr_t op_im, mpfr_prec_t wp) {
    mpfr_t s1r, s1i, s2r, s2i, u, v;
    mpfr_inits2(wp, s1r, s1i, s2r, s2i, u, v, (mpfr_ptr)NULL);
    mc_sqrt(s1r, s1i, om_re, om_im);
    mc_sqrt(s2r, s2i, op_re, op_im);
    mpfr_atan2(u, s1r, s2r, MPFR_RNDN);
    mpfr_mul_2ui(out_re, u, 1, MPFR_RNDN);
    /* Im(conj(s2) s1) = s2r s1i - s2i s1r  (both terms of -b's sign) */
    mpfr_mul(u, s2r, s1i, MPFR_RNDN);
    mpfr_mul(v, s2i, s1r, MPFR_RNDN);
    mpfr_sub(u, u, v, MPFR_RNDN);
    mpfr_asinh(out_im, u, MPFR_RNDN);
    mpfr_clears(s1r, s1i, s2r, s2i, u, v, (mpfr_ptr)NULL);
}

/* acosh from zm = z - 1, zp = z + 1. */
static void nc_acosh_core(mpfr_t out_re, mpfr_t out_im,
                          const mpfr_t zm_re, const mpfr_t zm_im,
                          const mpfr_t zp_re, const mpfr_t zp_im, mpfr_prec_t wp) {
    mpfr_t r1r, r1i, r2r, r2i, u, v;
    mpfr_inits2(wp, r1r, r1i, r2r, r2i, u, v, (mpfr_ptr)NULL);
    mc_sqrt(r1r, r1i, zm_re, zm_im);
    mc_sqrt(r2r, r2i, zp_re, zp_im);
    /* Re(conj(r1) r2) = r1r r2r + r1i r2i  (r1i, r2i share b's sign) */
    mpfr_mul(u, r1r, r2r, MPFR_RNDN);
    mpfr_mul(v, r1i, r2i, MPFR_RNDN);
    mpfr_add(u, u, v, MPFR_RNDN);
    mpfr_asinh(out_re, u, MPFR_RNDN);
    mpfr_atan2(u, r1i, r2r, MPFR_RNDN);
    mpfr_mul_2ui(out_im, u, 1, MPFR_RNDN);
    mpfr_clears(r1r, r1i, r2r, r2i, u, v, (mpfr_ptr)NULL);
}

/* atanh (coth = false) or acoth (coth = true) of a + bi.  Both have the
 * real part log(((a+1)^2 + b^2)/((a-1)^2 + b^2))/4 -- log1p(4a/den)/4 when
 * the ratio num/den = 1 + 4a/den is in [1/2, 2], so a tiny z keeps its
 * relative accuracy, and the log of the ratio outside it, where 1 + 4a/den
 * would cancel (z next to -1: num -> 0) -- and the
 * imaginary parts
 *   atanh: atan2( 2b, (1 - a)(1 + a) - b^2)/2,  on the cut b = 0, |a| > 1
 *   acoth: atan2(-2b, (a - 1)(a + 1) + b^2)/2,  on the cut b = 0, |a| < 1
 * where the definition gives -sign(a) Pi/2. */
static void nc_atanh_core(mpfr_t out_re, mpfr_t out_im,
                          const mpfr_t a, const mpfr_t b, bool coth, mpfr_prec_t wp) {
    mpfr_t am1, ap1, b2, num, den, u, v;
    mpfr_inits2(wp, am1, ap1, b2, num, den, u, v, (mpfr_ptr)NULL);
    mpfr_sub_ui(am1, a, 1, MPFR_RNDN);
    mpfr_add_ui(ap1, a, 1, MPFR_RNDN);
    mpfr_sqr(b2, b, MPFR_RNDN);
    mpfr_sqr(num, ap1, MPFR_RNDN);
    mpfr_add(num, num, b2, MPFR_RNDN);
    mpfr_sqr(den, am1, MPFR_RNDN);
    mpfr_add(den, den, b2, MPFR_RNDN);
    /* u = 4a / den = num/den - 1 */
    mpfr_mul_2ui(u, a, 2, MPFR_RNDN);
    mpfr_div(u, u, den, MPFR_RNDN);
    if (mpfr_number_p(u) && mpfr_cmp_d(u, -0.5) >= 0 && mpfr_cmp_ui(u, 1) <= 0) {
        mpfr_log1p(u, u, MPFR_RNDN);
    } else {
        mpfr_div(u, num, den, MPFR_RNDN);
        mpfr_log(u, u, MPFR_RNDN);
    }
    mpfr_div_2ui(out_re, u, 2, MPFR_RNDN);

    bool on_cut = mpfr_zero_p(b) &&
                  (coth ? mpfr_cmpabs_ui(a, 1) < 0 : mpfr_cmpabs_ui(a, 1) > 0);
    if (on_cut && !mpfr_zero_p(a)) {
        mpfr_const_pi(u, MPFR_RNDN);
        mpfr_div_2ui(u, u, 1, MPFR_RNDN);
        if (mpfr_sgn(a) > 0) mpfr_neg(u, u, MPFR_RNDN);
        mpfr_set(out_im, u, MPFR_RNDN);
    } else if (on_cut) {
        /* acoth(0): the definition's limit along the imaginary axis from
         * the right, +i Pi/2 -- what the machine kernel (catanh(1/0)) and
         * the exact rule give. */
        mpfr_const_pi(u, MPFR_RNDN);
        mpfr_div_2ui(out_im, u, 1, MPFR_RNDN);
    } else {
        mpfr_mul(u, am1, ap1, MPFR_RNDN);          /* a^2 - 1 */
        mpfr_mul_2ui(v, b, 1, MPFR_RNDN);
        if (coth) {
            mpfr_add(u, u, b2, MPFR_RNDN);
            mpfr_neg(v, v, MPFR_RNDN);
        } else {
            mpfr_neg(u, u, MPFR_RNDN);
            mpfr_sub(u, u, b2, MPFR_RNDN);
        }
        mpfr_atan2(u, v, u, MPFR_RNDN);
        mpfr_div_2ui(out_im, u, 1, MPFR_RNDN);
    }
    mpfr_clears(am1, ap1, b2, num, den, u, v, (mpfr_ptr)NULL);
}

/* (out_re, out_im) = (q_re, q_im) rounded; clears nothing. */
static void nc_round_out(mpfr_t out_re, mpfr_t out_im, const mpfr_t r, const mpfr_t i) {
    mpfr_set(out_re, r, MPFR_RNDN);
    mpfr_set(out_im, i, MPFR_RNDN);
}

void mpfr_complex_asin(mpfr_t out_re, mpfr_t out_im,
                       const mpfr_t a, const mpfr_t b) {
    mpfr_prec_t wp = mpfr_get_prec(out_re) + NC_GUARD;
    mpfr_t omr, omi, opr, opi, rr, ri;
    mpfr_inits2(wp, omr, omi, opr, opi, rr, ri, (mpfr_ptr)NULL);
    mpfr_ui_sub(omr, 1, a, MPFR_RNDN);
    mpfr_neg(omi, b, MPFR_RNDN);
    mpfr_add_ui(opr, a, 1, MPFR_RNDN);
    mpfr_set(opi, b, MPFR_RNDN);
    nc_asin_core(rr, ri, a, omr, omi, opr, opi, wp);
    nc_round_out(out_re, out_im, rr, ri);
    mpfr_clears(omr, omi, opr, opi, rr, ri, (mpfr_ptr)NULL);
}

void mpfr_complex_acos(mpfr_t out_re, mpfr_t out_im,
                       const mpfr_t a, const mpfr_t b) {
    mpfr_prec_t wp = mpfr_get_prec(out_re) + NC_GUARD;
    mpfr_t omr, omi, opr, opi, rr, ri;
    mpfr_inits2(wp, omr, omi, opr, opi, rr, ri, (mpfr_ptr)NULL);
    mpfr_ui_sub(omr, 1, a, MPFR_RNDN);
    mpfr_neg(omi, b, MPFR_RNDN);
    mpfr_add_ui(opr, a, 1, MPFR_RNDN);
    mpfr_set(opi, b, MPFR_RNDN);
    nc_acos_core(rr, ri, omr, omi, opr, opi, wp);
    nc_round_out(out_re, out_im, rr, ri);
    mpfr_clears(omr, omi, opr, opi, rr, ri, (mpfr_ptr)NULL);
}

/* atan z = -i atanh(i z), i z = -b + a i. */
void mpfr_complex_atan(mpfr_t out_re, mpfr_t out_im,
                       const mpfr_t a, const mpfr_t b) {
    mpfr_prec_t wp = mpfr_get_prec(out_re) + NC_GUARD;
    mpfr_t wr, rr, ri;
    mpfr_inits2(wp, wr, rr, ri, (mpfr_ptr)NULL);
    mpfr_neg(wr, b, MPFR_RNDN);
    nc_atanh_core(rr, ri, wr, a, false, wp);
    mpfr_neg(rr, rr, MPFR_RNDN);
    nc_round_out(out_re, out_im, ri, rr);       /* (Im, -Re) */
    mpfr_clears(wr, rr, ri, (mpfr_ptr)NULL);
}

/* asinh z = -i asin(i z), i z = w = -b + a i: 1 - w = (1 + b) - a i,
 * 1 + w = (1 - b) + a i. */
void mpfr_complex_asinh(mpfr_t out_re, mpfr_t out_im,
                        const mpfr_t a, const mpfr_t b) {
    mpfr_prec_t wp = mpfr_get_prec(out_re) + NC_GUARD;
    mpfr_t x, omr, omi, opr, opi, rr, ri;
    mpfr_inits2(wp, x, omr, omi, opr, opi, rr, ri, (mpfr_ptr)NULL);
    mpfr_neg(x, b, MPFR_RNDN);
    mpfr_add_ui(omr, b, 1, MPFR_RNDN);
    mpfr_neg(omi, a, MPFR_RNDN);
    mpfr_ui_sub(opr, 1, b, MPFR_RNDN);
    mpfr_set(opi, a, MPFR_RNDN);
    nc_asin_core(rr, ri, x, omr, omi, opr, opi, wp);
    mpfr_neg(rr, rr, MPFR_RNDN);
    nc_round_out(out_re, out_im, ri, rr);       /* (Im, -Re) */
    mpfr_clears(x, omr, omi, opr, opi, rr, ri, (mpfr_ptr)NULL);
}

void mpfr_complex_acosh(mpfr_t out_re, mpfr_t out_im,
                        const mpfr_t a, const mpfr_t b) {
    mpfr_prec_t wp = mpfr_get_prec(out_re) + NC_GUARD;
    mpfr_t zmr, zmi, zpr, zpi, rr, ri;
    mpfr_inits2(wp, zmr, zmi, zpr, zpi, rr, ri, (mpfr_ptr)NULL);
    mpfr_sub_ui(zmr, a, 1, MPFR_RNDN);
    mpfr_set(zmi, b, MPFR_RNDN);
    mpfr_add_ui(zpr, a, 1, MPFR_RNDN);
    mpfr_set(zpi, b, MPFR_RNDN);
    nc_acosh_core(rr, ri, zmr, zmi, zpr, zpi, wp);
    nc_round_out(out_re, out_im, rr, ri);
    mpfr_clears(zmr, zmi, zpr, zpi, rr, ri, (mpfr_ptr)NULL);
}

void mpfr_complex_atanh(mpfr_t out_re, mpfr_t out_im,
                        const mpfr_t a, const mpfr_t b) {
    mpfr_prec_t wp = mpfr_get_prec(out_re) + NC_GUARD;
    mpfr_t rr, ri;
    mpfr_inits2(wp, rr, ri, (mpfr_ptr)NULL);
    nc_atanh_core(rr, ri, a, b, false, wp);
    nc_round_out(out_re, out_im, rr, ri);
    mpfr_clears(rr, ri, (mpfr_ptr)NULL);
}

/* --------------------------------------------------------------------
 *  Reciprocal trig / hyperbolic composites
 * ------------------------------------------------------------------ */

void mpfr_complex_cot (mpfr_t r, mpfr_t i, const mpfr_t a, const mpfr_t b) { complex_tan_family(r, i, a, b, false, NC_COT); }
void mpfr_complex_sec (mpfr_t r, mpfr_t i, const mpfr_t a, const mpfr_t b) { complex_tan_family(r, i, a, b, false, NC_SEC); }
void mpfr_complex_csc (mpfr_t r, mpfr_t i, const mpfr_t a, const mpfr_t b) { complex_tan_family(r, i, a, b, false, NC_CSC); }
void mpfr_complex_coth(mpfr_t r, mpfr_t i, const mpfr_t a, const mpfr_t b) { complex_tan_family(r, i, a, b, true,  NC_COT); }
void mpfr_complex_sech(mpfr_t r, mpfr_t i, const mpfr_t a, const mpfr_t b) { complex_tan_family(r, i, a, b, true,  NC_SEC); }
void mpfr_complex_csch(mpfr_t r, mpfr_t i, const mpfr_t a, const mpfr_t b) { complex_tan_family(r, i, a, b, true,  NC_CSC); }

/* (z + c)/z for c = c_re + c_im i (one of them zero), at the precision of
 * the outputs, as (z + c) conj(z) / |z|^2:
 *   Re = (a (a + c_re) + b (b + c_im)) / |z|^2,  Im = (c_im a - c_re b) / |z|^2.
 * Every part is then accurate relative to itself.  A general complex division
 * (Smith's) is accurate only relative to the larger part, and the small one
 * decides the side of a cut: for z = -10^100 + I, (z - 1)/z has Im = 10^-200,
 * which Smith's algorithm rounded to an exact 0, so ArcSech read the point as
 * on its cut and returned +Pi/2 I for -Pi/2 I.  Next to a branch point
 * (z = 1 + 10^-25) a - 1 is exact, so nothing is lost to rounding 1/z
 * either. */
static void nc_shift_over_z(mpfr_t out_re, mpfr_t out_im, long c_re, long c_im,
                            const mpfr_t a, const mpfr_t b) {
    mpfr_prec_t wp = mpfr_get_prec(out_re);
    mpfr_t t, u, n2;
    mpfr_inits2(wp, t, u, n2, (mpfr_ptr)NULL);
    mpfr_sqr(n2, a, MPFR_RNDN);
    mpfr_fma(n2, b, b, n2, MPFR_RNDN);            /* |z|^2 */
    mpfr_add_si(t, a, c_re, MPFR_RNDN);
    mpfr_mul(t, t, a, MPFR_RNDN);
    mpfr_add_si(u, b, c_im, MPFR_RNDN);
    mpfr_mul(u, u, b, MPFR_RNDN);
    mpfr_add(t, t, u, MPFR_RNDN);
    /* Im numerator: c_im a - c_re b, exact up to a sign (|c| = 1) */
    if (c_im != 0) mpfr_mul_si(u, a, c_im, MPFR_RNDN);
    else           mpfr_mul_si(u, b, -c_re, MPFR_RNDN);
    mpfr_div(out_re, t, n2, MPFR_RNDN);
    mpfr_div(out_im, u, n2, MPFR_RNDN);
    mpfr_clears(t, u, n2, (mpfr_ptr)NULL);
}

/* acot z = i acoth(i z), i z = -b + a i. */
void mpfr_complex_acot(mpfr_t out_re, mpfr_t out_im, const mpfr_t a, const mpfr_t b) {
    mpfr_prec_t wp = mpfr_get_prec(out_re) + NC_GUARD;
    mpfr_t wr, rr, ri;
    mpfr_inits2(wp, wr, rr, ri, (mpfr_ptr)NULL);
    mpfr_neg(wr, b, MPFR_RNDN);
    nc_atanh_core(rr, ri, wr, a, true, wp);
    mpfr_neg(ri, ri, MPFR_RNDN);
    nc_round_out(out_re, out_im, ri, rr);       /* (-Im, Re) */
    mpfr_clears(wr, rr, ri, (mpfr_ptr)NULL);
}

/* asec z = acos(1/z): 1 - 1/z = (z - 1)/z, 1 + 1/z = (z + 1)/z. */
void mpfr_complex_asec(mpfr_t out_re, mpfr_t out_im, const mpfr_t a, const mpfr_t b) {
    mpfr_prec_t wp = mpfr_get_prec(out_re) + NC_GUARD;
    mpfr_t omr, omi, opr, opi, rr, ri;
    mpfr_inits2(wp, omr, omi, opr, opi, rr, ri, (mpfr_ptr)NULL);
    nc_shift_over_z(omr, omi, -1, 0, a, b);
    nc_shift_over_z(opr, opi, 1, 0, a, b);
    nc_acos_core(rr, ri, omr, omi, opr, opi, wp);
    nc_round_out(out_re, out_im, rr, ri);
    mpfr_clears(omr, omi, opr, opi, rr, ri, (mpfr_ptr)NULL);
}

/* acsc z = asin(1/z), Re(1/z) = a / |z|^2. */
void mpfr_complex_acsc(mpfr_t out_re, mpfr_t out_im, const mpfr_t a, const mpfr_t b) {
    mpfr_prec_t wp = mpfr_get_prec(out_re) + NC_GUARD;
    mpfr_t x, y, one, zero, omr, omi, opr, opi, rr, ri;
    mpfr_inits2(wp, x, y, one, zero, omr, omi, opr, opi, rr, ri, (mpfr_ptr)NULL);
    mpfr_set_ui(one, 1, MPFR_RNDN);
    mpfr_set_zero(zero, 1);
    mpfr_complex_div(x, y, one, zero, a, b);
    nc_shift_over_z(omr, omi, -1, 0, a, b);
    nc_shift_over_z(opr, opi, 1, 0, a, b);
    nc_asin_core(rr, ri, x, omr, omi, opr, opi, wp);
    nc_round_out(out_re, out_im, rr, ri);
    mpfr_clears(x, y, one, zero, omr, omi, opr, opi, rr, ri, (mpfr_ptr)NULL);
}

/* acoth z, directly (not atanh(1/z), whose rounding of 1/z costs digits next
 * to the branch points +-1). */
void mpfr_complex_acoth(mpfr_t out_re, mpfr_t out_im, const mpfr_t a, const mpfr_t b) {
    mpfr_prec_t wp = mpfr_get_prec(out_re) + NC_GUARD;
    mpfr_t rr, ri;
    mpfr_inits2(wp, rr, ri, (mpfr_ptr)NULL);
    nc_atanh_core(rr, ri, a, b, true, wp);
    nc_round_out(out_re, out_im, rr, ri);
    mpfr_clears(rr, ri, (mpfr_ptr)NULL);
}

/* asech z = acosh(1/z): 1/z - 1 = (1 - z)/z = -(z - 1)/z, 1/z + 1 = (z + 1)/z. */
void mpfr_complex_asech(mpfr_t out_re, mpfr_t out_im, const mpfr_t a, const mpfr_t b) {
    mpfr_prec_t wp = mpfr_get_prec(out_re) + NC_GUARD;
    mpfr_t zmr, zmi, zpr, zpi, rr, ri;
    mpfr_inits2(wp, zmr, zmi, zpr, zpi, rr, ri, (mpfr_ptr)NULL);
    nc_shift_over_z(zmr, zmi, -1, 0, a, b);
    mpfr_neg(zmr, zmr, MPFR_RNDN);
    mpfr_neg(zmi, zmi, MPFR_RNDN);
    nc_shift_over_z(zpr, zpi, 1, 0, a, b);
    nc_acosh_core(rr, ri, zmr, zmi, zpr, zpi, wp);
    nc_round_out(out_re, out_im, rr, ri);
    mpfr_clears(zmr, zmi, zpr, zpi, rr, ri, (mpfr_ptr)NULL);
}

/* acsch z = asinh(1/z) = -i asin(i/z): w = i/z, Re w = b/|z|^2,
 * 1 - w = (z - i)/z, 1 + w = (z + i)/z. */
void mpfr_complex_acsch(mpfr_t out_re, mpfr_t out_im, const mpfr_t a, const mpfr_t b) {
    mpfr_prec_t wp = mpfr_get_prec(out_re) + NC_GUARD;
    mpfr_t x, y, one, zero, omr, omi, opr, opi, rr, ri;
    mpfr_inits2(wp, x, y, one, zero, omr, omi, opr, opi, rr, ri, (mpfr_ptr)NULL);
    mpfr_set_ui(one, 1, MPFR_RNDN);
    mpfr_set_zero(zero, 1);
    mpfr_complex_div(x, y, zero, one, a, b);     /* i / z */
    nc_shift_over_z(omr, omi, 0, -1, a, b);
    nc_shift_over_z(opr, opi, 0, 1, a, b);
    nc_asin_core(rr, ri, x, omr, omi, opr, opi, wp);
    mpfr_neg(rr, rr, MPFR_RNDN);
    nc_round_out(out_re, out_im, ri, rr);       /* (Im, -Re) */
    mpfr_clears(x, y, one, zero, omr, omi, opr, opi, rr, ri, (mpfr_ptr)NULL);
}

/* --------------------------------------------------------------------
 *  Public dispatcher
 * ------------------------------------------------------------------ */

Expr* numeric_mpfr_apply_complex_unary(const Expr* e, long default_bits,
                                       MpfrComplexUnaryOp op) {
    long bits = numeric_combined_bits(e, NULL, default_bits);
    mpfr_t in_re, in_im;
    mpfr_init2(in_re, bits);
    mpfr_init2(in_im, bits);
    bool ok = get_approx_mpfr(e, in_re, in_im, NULL);
    if (!ok) {
        mpfr_clear(in_re);
        mpfr_clear(in_im);
        return NULL;
    }
    Expr* re_expr = expr_new_mpfr_bits(bits);
    Expr* im_expr = expr_new_mpfr_bits(bits);
    op(re_expr->data.mpfr, im_expr->data.mpfr, in_re, in_im);
    bool im_zero = mpfr_zero_p(im_expr->data.mpfr);
    mpfr_clear(in_re);
    mpfr_clear(in_im);
    if (im_zero) {
        expr_free(im_expr);
        return re_expr;
    }
    return make_complex(re_expr, im_expr);
}

/* --------------------------------------------------------------------
 *  ncpx — raw complex-MPFR arithmetic (see numeric_complex.h).
 *
 *  Each op is alias-safe: every input component that must outlive a write
 *  to the output is read into scratch first. Binary ops that allocate take
 *  the working precision `wp` for that scratch.
 * ------------------------------------------------------------------ */

#define NRND MPFR_RNDN

void ncpx_init (ncpx* z, mpfr_prec_t p) { mpfr_init2(z->re, p); mpfr_init2(z->im, p); }
void ncpx_clear(ncpx* z)                { mpfr_clear(z->re);    mpfr_clear(z->im);    }

void ncpx_set(ncpx* d, const ncpx* s) {
    mpfr_set(d->re, s->re, NRND);
    mpfr_set(d->im, s->im, NRND);
}
void ncpx_set_d(ncpx* d, double re, double im) {
    mpfr_set_d(d->re, re, NRND);
    mpfr_set_d(d->im, im, NRND);
}
void ncpx_set_ui(ncpx* d, unsigned long re) {
    mpfr_set_ui(d->re, re, NRND);
    mpfr_set_ui(d->im, 0,  NRND);
}

void ncpx_add(ncpx* out, const ncpx* a, const ncpx* b) {
    mpfr_add(out->re, a->re, b->re, NRND);
    mpfr_add(out->im, a->im, b->im, NRND);
}
void ncpx_sub(ncpx* out, const ncpx* a, const ncpx* b) {
    mpfr_sub(out->re, a->re, b->re, NRND);
    mpfr_sub(out->im, a->im, b->im, NRND);
}
void ncpx_neg(ncpx* out, const ncpx* a) {
    mpfr_neg(out->re, a->re, NRND);
    mpfr_neg(out->im, a->im, NRND);
}

/* out = a * b. */
void ncpx_mul(ncpx* out, const ncpx* a, const ncpx* b, mpfr_prec_t wp) {
    mpfr_t ac, bd, ad, bc;
    mpfr_inits2(wp, ac, bd, ad, bc, (mpfr_ptr)0);
    mpfr_mul(ac, a->re, b->re, NRND);
    mpfr_mul(bd, a->im, b->im, NRND);
    mpfr_mul(ad, a->re, b->im, NRND);
    mpfr_mul(bc, a->im, b->re, NRND);
    mpfr_sub(out->re, ac, bd, NRND);
    mpfr_add(out->im, ad, bc, NRND);
    mpfr_clears(ac, bd, ad, bc, (mpfr_ptr)0);
}

/* out = a / b. */
void ncpx_div(ncpx* out, const ncpx* a, const ncpx* b, mpfr_prec_t wp) {
    mpfr_t ac, bd, ad, bc, den;
    mpfr_inits2(wp, ac, bd, ad, bc, den, (mpfr_ptr)0);
    mpfr_mul(ac, a->re, b->re, NRND);
    mpfr_mul(bd, a->im, b->im, NRND);
    mpfr_mul(ad, a->im, b->re, NRND);
    mpfr_mul(bc, a->re, b->im, NRND);
    mpfr_mul(den, b->re, b->re, NRND);
    mpfr_fma(den, b->im, b->im, den, NRND);   /* |b|^2 */
    mpfr_add(ac, ac, bd, NRND);               /* re num = ac + bd */
    mpfr_sub(ad, ad, bc, NRND);               /* im num = ad - bc */
    mpfr_div(out->re, ac, den, NRND);
    mpfr_div(out->im, ad, den, NRND);
    mpfr_clears(ac, bd, ad, bc, den, (mpfr_ptr)0);
}

void ncpx_scale(ncpx* out, const ncpx* z, const mpfr_t s) {
    mpfr_mul(out->re, z->re, s, NRND);
    mpfr_mul(out->im, z->im, s, NRND);
}

void ncpx_abs(mpfr_t mag, const ncpx* z) { mpfr_hypot(mag, z->re, z->im, NRND); }
void ncpx_arg(mpfr_t out, const ncpx* z) { mpfr_atan2(out, z->im, z->re, NRND); }

/* out = exp(z) = e^re (cos im + i sin im). */
void ncpx_exp(ncpx* out, const ncpx* z, mpfr_prec_t wp) {
    mpfr_t ea, c, s;
    mpfr_inits2(wp, ea, c, s, (mpfr_ptr)0);
    mpfr_exp(ea, z->re, NRND);
    mpfr_sin_cos(s, c, z->im, NRND);
    mpfr_mul(out->re, ea, c, NRND);
    mpfr_mul(out->im, ea, s, NRND);
    mpfr_clears(ea, c, s, (mpfr_ptr)0);
}

/* out = log(z) = log|z| + i arg(z). */
void ncpx_log(ncpx* out, const ncpx* z, mpfr_prec_t wp) {
    mpfr_t r, th;
    mpfr_inits2(wp, r, th, (mpfr_ptr)0);
    mpfr_hypot(r, z->re, z->im, NRND);
    mpfr_atan2(th, z->im, z->re, NRND);
    mpfr_log(out->re, r, NRND);
    mpfr_set(out->im, th, NRND);
    mpfr_clears(r, th, (mpfr_ptr)0);
}

/* out = sin(z) = sin(re) cosh(im) + i cos(re) sinh(im). */
void ncpx_sin(ncpx* out, const ncpx* z, mpfr_prec_t wp) {
    mpfr_t sa, ca, sh, ch;
    mpfr_inits2(wp, sa, ca, sh, ch, (mpfr_ptr)0);
    mpfr_sin_cos(sa, ca, z->re, NRND);
    mpfr_sinh_cosh(sh, ch, z->im, NRND);
    mpfr_mul(out->re, sa, ch, NRND);
    mpfr_mul(out->im, ca, sh, NRND);
    mpfr_clears(sa, ca, sh, ch, (mpfr_ptr)0);
}

/* out = cos(z) = cos(re) cosh(im) - i sin(re) sinh(im). */
void ncpx_cos(ncpx* out, const ncpx* z, mpfr_prec_t wp) {
    mpfr_t sa, ca, sh, ch;
    mpfr_inits2(wp, sa, ca, sh, ch, (mpfr_ptr)0);
    mpfr_sin_cos(sa, ca, z->re, NRND);
    mpfr_sinh_cosh(sh, ch, z->im, NRND);
    mpfr_mul(out->re, ca, ch, NRND);
    mpfr_mul(out->im, sa, sh, NRND);
    mpfr_neg(out->im, out->im, NRND);
    mpfr_clears(sa, ca, sh, ch, (mpfr_ptr)0);
}

/* out = tan, csc, sec or cot of z (which = 0, 1, 2, 3) as a quotient of the
 * componentwise-accurate sin and cos above (8 guard bits): tan = sin/cos,
 * csc = 1/sin, sec = 1/cos, cot = cos/sin.  A pole gives a zero divisor, hence
 * a non-finite result, which a compiled caller declines. */
void ncpx_trig_ratio(ncpx* out, const ncpx* z, int which, mpfr_prec_t wp) {
    ncpx s, c;
    mpfr_prec_t gp = wp + 8;
    ncpx_init(&s, gp); ncpx_init(&c, gp);
    /* one sin_cos of Re z and one sinh_cosh of Im z serve both:
     * sin z = sin a cosh b + i cos a sinh b, cos z = cos a cosh b - i sin a sinh b */
    mpfr_t sa, ca, sh, ch;
    mpfr_inits2(gp, sa, ca, sh, ch, (mpfr_ptr)0);
    mpfr_sin_cos(sa, ca, z->re, NRND);
    mpfr_sinh_cosh(sh, ch, z->im, NRND);
    mpfr_mul(s.re, sa, ch, NRND);
    mpfr_mul(s.im, ca, sh, NRND);
    mpfr_mul(c.re, ca, ch, NRND);
    mpfr_mul(c.im, sa, sh, NRND);
    mpfr_neg(c.im, c.im, NRND);
    mpfr_clears(sa, ca, sh, ch, (mpfr_ptr)0);
    switch (which) {
        case 0: ncpx_div(out, &s, &c, wp); break;
        case 2: { ncpx_set_ui(&s, 1); ncpx_div(out, &s, &c, wp); break; }
        case 3: ncpx_div(out, &c, &s, wp); break;
        default: { ncpx_set_ui(&c, 1); ncpx_div(out, &c, &s, wp); break; }
    }
    ncpx_clear(&s); ncpx_clear(&c);
}

/* out = sqrt(z), principal branch (Re >= 0), via polar form. */
void ncpx_sqrt(ncpx* out, const ncpx* z, mpfr_prec_t wp) {
    /* Algebraic form, no atan2 / sin_cos: with r = |z|, z = a + i b,
     *   a >= 0:  u = sqrt((r + a)/2),  v = b/(2u);
     *   a <  0:  w = sqrt((r - a)/2),  u = |b|/(2w),  v = sign(b) w,
     * each branch free of cancellation (r + a and r - a add like signs).  The
     * sign of b, a signed zero included, picks the side of the cut, as
     * atan2(b, a)/2 did: sqrt(-4 - 0 i) = -2 i. */
    mpfr_prec_t gp = wp + 4;
    mpfr_t r, u, w;
    mpfr_inits2(gp, r, u, w, (mpfr_ptr)0);
    mpfr_hypot(r, z->re, z->im, NRND);
    if (mpfr_zero_p(r)) {
        mpfr_set_zero(out->re, 1); mpfr_set_zero(out->im, 1);
    } else if (mpfr_sgn(z->re) >= 0) {
        mpfr_add(u, r, z->re, NRND);
        mpfr_div_2ui(u, u, 1, NRND);
        mpfr_sqrt(u, u, NRND);
        mpfr_div(w, z->im, u, NRND);
        mpfr_div_2ui(w, w, 1, NRND);
        mpfr_set(out->re, u, NRND);
        mpfr_set(out->im, w, NRND);
    } else {
        mpfr_sub(w, r, z->re, NRND);
        mpfr_div_2ui(w, w, 1, NRND);
        mpfr_sqrt(w, w, NRND);
        mpfr_abs(u, z->im, NRND);
        mpfr_div(u, u, w, NRND);
        mpfr_div_2ui(u, u, 1, NRND);
        mpfr_set(out->re, u, NRND);
        if (mpfr_signbit(z->im)) mpfr_neg(w, w, NRND);
        mpfr_set(out->im, w, NRND);
    }
    mpfr_clears(r, u, w, (mpfr_ptr)0);
}

/* out = z^e (principal branch) for a real exponent e:
 *   out = |z|^e (cos(e arg z) + i sin(e arg z)). */
void ncpx_pow_d(ncpx* out, const ncpx* z, double e, mpfr_prec_t wp) {
    mpfr_t r, th, c, s, ex;
    mpfr_inits2(wp, r, th, c, s, ex, (mpfr_ptr)0);
    mpfr_hypot(r, z->re, z->im, NRND);
    mpfr_atan2(th, z->im, z->re, NRND);
    mpfr_set_d(ex, e, NRND);
    mpfr_pow(r, r, ex, NRND);              /* |z|^e */
    mpfr_mul_d(th, th, e, NRND);           /* e arg z */
    mpfr_sin_cos(s, c, th, NRND);
    mpfr_mul(out->re, r, c, NRND);
    mpfr_mul(out->im, r, s, NRND);
    mpfr_clears(r, th, c, s, ex, (mpfr_ptr)0);
}

/* d = a * b on caller-owned scratch t[4] (no allocation per product); d may
 * alias a or b. */
static void ncpx_mul_scr(ncpx* d, const ncpx* a, const ncpx* b, mpfr_t* t) {
    mpfr_mul(t[0], a->re, b->re, NRND);
    mpfr_mul(t[1], a->im, b->im, NRND);
    mpfr_mul(t[2], a->re, b->im, NRND);
    mpfr_mul(t[3], a->im, b->re, NRND);
    mpfr_sub(d->re, t[0], t[1], NRND);
    mpfr_add(d->im, t[2], t[3], NRND);
}

/* out = z^n for an integer n, by binary powering: about 2 log2|n| complex
 * products instead of the polar form's hypot, atan2, pow and sin_cos -- the
 * dominant cost of a compiled polynomial factor (x + 2)^1500 sampled on a
 * contour at hundreds of digits.  It is also the more accurate route: the
 * polar angle n arg z carries n times the angle's rounding error, the products
 * only O(log n) roundings, which the guard bits below absorb.  n < 0 inverts
 * once at the end; 0^n for n < 0 is a division by zero (non-finite, so a
 * compiled caller falls back), and z^0 = 1. */
void ncpx_pow_si(ncpx* out, const ncpx* z, long n, mpfr_prec_t wp) {
    unsigned long m = (n < 0) ? 0UL - (unsigned long)n : (unsigned long)n;
    mpfr_prec_t gp = wp + 8;
    for (unsigned long k = m; k; k >>= 1) gp++;
    ncpx b, r;
    mpfr_t t[4];
    ncpx_init(&b, gp); ncpx_init(&r, gp);
    for (int i = 0; i < 4; i++) mpfr_init2(t[i], gp);
    ncpx_set(&b, z);
    mpfr_set_ui(r.re, 1, NRND); mpfr_set_ui(r.im, 0, NRND);
    bool first = true;
    for (; m; m >>= 1) {
        if (m & 1UL) {
            if (first) { ncpx_set(&r, &b); first = false; }
            else ncpx_mul_scr(&r, &r, &b, t);
        }
        if (m > 1UL) ncpx_mul_scr(&b, &b, &b, t);
    }
    if (n < 0) {
        ncpx one;
        ncpx_init(&one, gp);
        mpfr_set_ui(one.re, 1, NRND); mpfr_set_ui(one.im, 0, NRND);
        ncpx_div(&r, &one, &r, gp);
        ncpx_clear(&one);
    }
    mpfr_set(out->re, r.re, NRND);
    mpfr_set(out->im, r.im, NRND);
    for (int i = 0; i < 4; i++) mpfr_clear(t[i]);
    ncpx_clear(&b); ncpx_clear(&r);
}

/* out = z^w = exp(w log z), principal branch.  An integer exponent (a real
 * w with zero imaginary part and an integral value) is z^n exactly, taken by
 * binary powering: the polar route costs a complex log and exp (an atan, a
 * sqrt, a log, an exp, a sin_cos) per power -- most of the time of a compiled
 * Product[1 - z^k, {k, 1, 1000}] at 50 digits -- and rounds more. */
void ncpx_pow(ncpx* out, const ncpx* z, const ncpx* w, mpfr_prec_t wp) {
    if (mpfr_zero_p(w->im) && mpfr_integer_p(w->re)
        && mpfr_cmpabs_ui(w->re, 1UL << 20) <= 0
        && !(mpfr_zero_p(z->re) && mpfr_zero_p(z->im))) {
        ncpx_pow_si(out, z, mpfr_get_si(w->re, NRND), wp);
        return;
    }
    ncpx l, m;
    ncpx_init(&l, wp); ncpx_init(&m, wp);
    ncpx_log(&l, z, wp);
    ncpx_mul(&m, w, &l, wp);
    ncpx_exp(out, &m, wp);
    ncpx_clear(&l); ncpx_clear(&m);
}

#endif /* USE_MPFR */
