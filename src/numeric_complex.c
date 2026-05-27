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

/* Complex multiply (out = a * b). out_re may alias inputs; the intermediate
 * goes through scratch to keep aliasing safe. */
static void mc_mul(mpfr_t out_re, mpfr_t out_im,
                   const mpfr_t a_re, const mpfr_t a_im,
                   const mpfr_t b_re, const mpfr_t b_im) {
    mpfr_prec_t p = mpfr_get_prec(out_re);
    mpfr_t t1, t2, re_scratch;
    mpfr_inits2(p, t1, t2, re_scratch, (mpfr_ptr)NULL);
    mpfr_mul(t1, a_re, b_re, MPFR_RNDN);
    mpfr_mul(t2, a_im, b_im, MPFR_RNDN);
    mpfr_sub(re_scratch, t1, t2, MPFR_RNDN);
    mpfr_mul(t1, a_re, b_im, MPFR_RNDN);
    mpfr_mul(t2, a_im, b_re, MPFR_RNDN);
    mpfr_add(out_im, t1, t2, MPFR_RNDN);
    mpfr_set(out_re, re_scratch, MPFR_RNDN);
    mpfr_clears(t1, t2, re_scratch, (mpfr_ptr)NULL);
}

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
 * ------------------------------------------------------------------ */

void mpfr_complex_exp(mpfr_t out_re, mpfr_t out_im,
                      const mpfr_t a, const mpfr_t b) {
    mpfr_prec_t p = mpfr_get_prec(out_re);
    mpfr_t ea, cb, sb;
    mpfr_inits2(p, ea, cb, sb, (mpfr_ptr)NULL);
    mpfr_exp(ea, a, MPFR_RNDN);
    mpfr_sin_cos(sb, cb, b, MPFR_RNDN);
    mpfr_mul(out_re, ea, cb, MPFR_RNDN);
    mpfr_mul(out_im, ea, sb, MPFR_RNDN);
    mpfr_clears(ea, cb, sb, (mpfr_ptr)NULL);
}

void mpfr_complex_log(mpfr_t out_re, mpfr_t out_im,
                      const mpfr_t a, const mpfr_t b) {
    mpfr_prec_t p = mpfr_get_prec(out_re);
    mpfr_t r;
    mpfr_init2(r, p);
    mpfr_hypot(r, a, b, MPFR_RNDN);
    mpfr_log(out_re, r, MPFR_RNDN);
    mpfr_atan2(out_im, b, a, MPFR_RNDN);
    mpfr_clear(r);
}

void mpfr_complex_sin(mpfr_t out_re, mpfr_t out_im,
                      const mpfr_t a, const mpfr_t b) {
    mpfr_prec_t p = mpfr_get_prec(out_re);
    mpfr_t sa, ca, sb, cb;
    mpfr_inits2(p, sa, ca, sb, cb, (mpfr_ptr)NULL);
    mpfr_sin_cos(sa, ca, a, MPFR_RNDN);
    mpfr_sinh_cosh(sb, cb, b, MPFR_RNDN);
    mpfr_mul(out_re, sa, cb, MPFR_RNDN);
    mpfr_mul(out_im, ca, sb, MPFR_RNDN);
    mpfr_clears(sa, ca, sb, cb, (mpfr_ptr)NULL);
}

void mpfr_complex_cos(mpfr_t out_re, mpfr_t out_im,
                      const mpfr_t a, const mpfr_t b) {
    mpfr_prec_t p = mpfr_get_prec(out_re);
    mpfr_t sa, ca, sb, cb;
    mpfr_inits2(p, sa, ca, sb, cb, (mpfr_ptr)NULL);
    mpfr_sin_cos(sa, ca, a, MPFR_RNDN);
    mpfr_sinh_cosh(sb, cb, b, MPFR_RNDN);
    mpfr_mul(out_re, ca, cb, MPFR_RNDN);
    mpfr_mul(out_im, sa, sb, MPFR_RNDN);
    mpfr_neg(out_im, out_im, MPFR_RNDN);
    mpfr_clears(sa, ca, sb, cb, (mpfr_ptr)NULL);
}

/* tan(a+bi) = (sin(2a) + i sinh(2b)) / (cos(2a) + cosh(2b))
 *
 * The 2a / 2b form keeps the denominator strictly positive (cos+cosh
 * is bounded below by cosh(2b) - 1 > 0 for |b| > 0 and exactly 0 only at
 * the poles where the symbolic Sin/Cos paths have already handled the
 * input), and avoids the sin(a)cos(a)/cos(a)^2 cancellation that would
 * appear in the naive form for arguments near pi/2. */
void mpfr_complex_tan(mpfr_t out_re, mpfr_t out_im,
                      const mpfr_t a, const mpfr_t b) {
    mpfr_prec_t p = mpfr_get_prec(out_re);
    mpfr_t two_a, two_b, sin2a, cos2a, sinh2b, cosh2b, den;
    mpfr_inits2(p, two_a, two_b, sin2a, cos2a, sinh2b, cosh2b, den, (mpfr_ptr)NULL);
    mpfr_mul_2ui(two_a, a, 1, MPFR_RNDN);
    mpfr_mul_2ui(two_b, b, 1, MPFR_RNDN);
    mpfr_sin_cos(sin2a, cos2a, two_a, MPFR_RNDN);
    mpfr_sinh_cosh(sinh2b, cosh2b, two_b, MPFR_RNDN);
    mpfr_add(den, cos2a, cosh2b, MPFR_RNDN);
    mpfr_div(out_re, sin2a, den, MPFR_RNDN);
    mpfr_div(out_im, sinh2b, den, MPFR_RNDN);
    mpfr_clears(two_a, two_b, sin2a, cos2a, sinh2b, cosh2b, den, (mpfr_ptr)NULL);
}

void mpfr_complex_sinh(mpfr_t out_re, mpfr_t out_im,
                       const mpfr_t a, const mpfr_t b) {
    /* sinh(a+bi) = sinh(a) cos(b) + i cosh(a) sin(b) */
    mpfr_prec_t p = mpfr_get_prec(out_re);
    mpfr_t sa, ca, sb, cb;
    mpfr_inits2(p, sa, ca, sb, cb, (mpfr_ptr)NULL);
    mpfr_sinh_cosh(sa, ca, a, MPFR_RNDN);
    mpfr_sin_cos(sb, cb, b, MPFR_RNDN);
    mpfr_mul(out_re, sa, cb, MPFR_RNDN);
    mpfr_mul(out_im, ca, sb, MPFR_RNDN);
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
    mpfr_mul(out_re, ca, cb, MPFR_RNDN);
    mpfr_mul(out_im, sa, sb, MPFR_RNDN);
    mpfr_clears(sa, ca, sb, cb, (mpfr_ptr)NULL);
}

/* tanh(a+bi) = (sinh(2a) + i sin(2b)) / (cosh(2a) + cos(2b)) — the dual
 * of the tan identity, same stability rationale. */
void mpfr_complex_tanh(mpfr_t out_re, mpfr_t out_im,
                       const mpfr_t a, const mpfr_t b) {
    mpfr_prec_t p = mpfr_get_prec(out_re);
    mpfr_t two_a, two_b, sinh2a, cosh2a, sin2b, cos2b, den;
    mpfr_inits2(p, two_a, two_b, sinh2a, cosh2a, sin2b, cos2b, den, (mpfr_ptr)NULL);
    mpfr_mul_2ui(two_a, a, 1, MPFR_RNDN);
    mpfr_mul_2ui(two_b, b, 1, MPFR_RNDN);
    mpfr_sinh_cosh(sinh2a, cosh2a, two_a, MPFR_RNDN);
    mpfr_sin_cos(sin2b, cos2b, two_b, MPFR_RNDN);
    mpfr_add(den, cosh2a, cos2b, MPFR_RNDN);
    mpfr_div(out_re, sinh2a, den, MPFR_RNDN);
    mpfr_div(out_im, sin2b, den, MPFR_RNDN);
    mpfr_clears(two_a, two_b, sinh2a, cosh2a, sin2b, cos2b, den, (mpfr_ptr)NULL);
}

/* --------------------------------------------------------------------
 *  Inverse transcendentals (log-form definitions)
 * ------------------------------------------------------------------ */

void mpfr_complex_asin(mpfr_t out_re, mpfr_t out_im,
                       const mpfr_t a, const mpfr_t b) {
    /* asin(z) = -i log(i z + sqrt(1 - z^2)) */
    mpfr_prec_t p = mpfr_get_prec(out_re);
    mpfr_t z2_re, z2_im, sq_re, sq_im, inner_re, inner_im, lr, li, t;
    mpfr_inits2(p, z2_re, z2_im, sq_re, sq_im, inner_re, inner_im, lr, li, t, (mpfr_ptr)NULL);
    /* 1 - z^2 = (1 - a^2 + b^2) + (-2 a b) i */
    mpfr_sqr(z2_re, a, MPFR_RNDN);
    mpfr_sqr(t,    b, MPFR_RNDN);
    mpfr_sub(z2_re, z2_re, t, MPFR_RNDN);
    mpfr_ui_sub(z2_re, 1, z2_re, MPFR_RNDN);
    mpfr_mul(z2_im, a, b, MPFR_RNDN);
    mpfr_mul_2ui(z2_im, z2_im, 1, MPFR_RNDN);
    mpfr_neg(z2_im, z2_im, MPFR_RNDN);
    mc_sqrt(sq_re, sq_im, z2_re, z2_im);
    /* inner = i z + sqrt(1 - z^2) = (-b + sq_re) + (a + sq_im) i */
    mpfr_sub(inner_re, sq_re, b, MPFR_RNDN);
    mpfr_add(inner_im, sq_im, a, MPFR_RNDN);
    mpfr_complex_log(lr, li, inner_re, inner_im);
    /* -i (lr + li i) = li - lr i */
    mpfr_set(out_re, li, MPFR_RNDN);
    mpfr_neg(out_im, lr, MPFR_RNDN);
    mpfr_clears(z2_re, z2_im, sq_re, sq_im, inner_re, inner_im, lr, li, t, (mpfr_ptr)NULL);
}

void mpfr_complex_acos(mpfr_t out_re, mpfr_t out_im,
                       const mpfr_t a, const mpfr_t b) {
    /* acos(z) = pi/2 - asin(z) */
    mpfr_complex_asin(out_re, out_im, a, b);
    mpfr_prec_t p = mpfr_get_prec(out_re);
    mpfr_t pi_half;
    mpfr_init2(pi_half, p);
    mpfr_const_pi(pi_half, MPFR_RNDN);
    mpfr_div_2ui(pi_half, pi_half, 1, MPFR_RNDN);
    mpfr_sub(out_re, pi_half, out_re, MPFR_RNDN);
    mpfr_neg(out_im, out_im, MPFR_RNDN);
    mpfr_clear(pi_half);
}

void mpfr_complex_atan(mpfr_t out_re, mpfr_t out_im,
                       const mpfr_t a, const mpfr_t b) {
    /* atan(z) = (i/2)(log(1 - i z) - log(1 + i z))
     *   1 - i z = (1 + b) + (-a) i
     *   1 + i z = (1 - b) + ( a) i */
    mpfr_prec_t p = mpfr_get_prec(out_re);
    mpfr_t u_re, u_im, v_re, v_im, lr1, li1, lr2, li2;
    mpfr_inits2(p, u_re, u_im, v_re, v_im, lr1, li1, lr2, li2, (mpfr_ptr)NULL);
    mpfr_add_ui(u_re, b, 1, MPFR_RNDN);
    mpfr_neg(u_im, a, MPFR_RNDN);
    mpfr_ui_sub(v_re, 1, b, MPFR_RNDN);
    mpfr_set(v_im, a, MPFR_RNDN);
    mpfr_complex_log(lr1, li1, u_re, u_im);
    mpfr_complex_log(lr2, li2, v_re, v_im);
    /*  (i/2)((lr1 + li1 i) - (lr2 + li2 i))
     *  = (i/2)((lr1 - lr2) + (li1 - li2) i)
     *  =  -(li1 - li2)/2 + (lr1 - lr2)/2 i */
    mpfr_sub(out_re, li2, li1, MPFR_RNDN);
    mpfr_div_2ui(out_re, out_re, 1, MPFR_RNDN);
    mpfr_sub(out_im, lr1, lr2, MPFR_RNDN);
    mpfr_div_2ui(out_im, out_im, 1, MPFR_RNDN);
    mpfr_clears(u_re, u_im, v_re, v_im, lr1, li1, lr2, li2, (mpfr_ptr)NULL);
}

void mpfr_complex_asinh(mpfr_t out_re, mpfr_t out_im,
                        const mpfr_t a, const mpfr_t b) {
    /* asinh(z) = log(z + sqrt(z^2 + 1))
     *   z^2 + 1 = (a^2 - b^2 + 1) + (2 a b) i */
    mpfr_prec_t p = mpfr_get_prec(out_re);
    mpfr_t z2_re, z2_im, sq_re, sq_im, t;
    mpfr_inits2(p, z2_re, z2_im, sq_re, sq_im, t, (mpfr_ptr)NULL);
    mpfr_sqr(z2_re, a, MPFR_RNDN);
    mpfr_sqr(t,     b, MPFR_RNDN);
    mpfr_sub(z2_re, z2_re, t, MPFR_RNDN);
    mpfr_add_ui(z2_re, z2_re, 1, MPFR_RNDN);
    mpfr_mul(z2_im, a, b, MPFR_RNDN);
    mpfr_mul_2ui(z2_im, z2_im, 1, MPFR_RNDN);
    mc_sqrt(sq_re, sq_im, z2_re, z2_im);
    mpfr_add(sq_re, a, sq_re, MPFR_RNDN);
    mpfr_add(sq_im, b, sq_im, MPFR_RNDN);
    mpfr_complex_log(out_re, out_im, sq_re, sq_im);
    mpfr_clears(z2_re, z2_im, sq_re, sq_im, t, (mpfr_ptr)NULL);
}

void mpfr_complex_acosh(mpfr_t out_re, mpfr_t out_im,
                        const mpfr_t a, const mpfr_t b) {
    /* acosh(z) = log(z + sqrt(z - 1) * sqrt(z + 1))
     * The two-sqrt product form keeps the branch cut on (-inf, 1) — the
     * Mathematica convention. The product is computed via mc_mul on the
     * two independent sqrt evaluations. */
    mpfr_prec_t p = mpfr_get_prec(out_re);
    mpfr_t zm_re, zm_im, zp_re, zp_im;
    mpfr_t s1_re, s1_im, s2_re, s2_im;
    mpfr_t prod_re, prod_im, inner_re, inner_im;
    mpfr_inits2(p, zm_re, zm_im, zp_re, zp_im,
                   s1_re, s1_im, s2_re, s2_im,
                   prod_re, prod_im, inner_re, inner_im,
                   (mpfr_ptr)NULL);
    mpfr_sub_ui(zm_re, a, 1, MPFR_RNDN);
    mpfr_set   (zm_im, b,    MPFR_RNDN);
    mpfr_add_ui(zp_re, a, 1, MPFR_RNDN);
    mpfr_set   (zp_im, b,    MPFR_RNDN);
    mc_sqrt(s1_re, s1_im, zm_re, zm_im);
    mc_sqrt(s2_re, s2_im, zp_re, zp_im);
    mc_mul(prod_re, prod_im, s1_re, s1_im, s2_re, s2_im);
    mpfr_add(inner_re, a, prod_re, MPFR_RNDN);
    mpfr_add(inner_im, b, prod_im, MPFR_RNDN);
    mpfr_complex_log(out_re, out_im, inner_re, inner_im);
    mpfr_clears(zm_re, zm_im, zp_re, zp_im,
                s1_re, s1_im, s2_re, s2_im,
                prod_re, prod_im, inner_re, inner_im,
                (mpfr_ptr)NULL);
}

void mpfr_complex_atanh(mpfr_t out_re, mpfr_t out_im,
                        const mpfr_t a, const mpfr_t b) {
    /* atanh(z) = (1/2) log((1 + z) / (1 - z))
     *   1 + z = (1 + a) + b i
     *   1 - z = (1 - a) + (-b) i */
    mpfr_prec_t p = mpfr_get_prec(out_re);
    mpfr_t u_re, u_im, v_re, v_im, q_re, q_im, lr, li;
    mpfr_inits2(p, u_re, u_im, v_re, v_im, q_re, q_im, lr, li, (mpfr_ptr)NULL);
    mpfr_add_ui(u_re, a, 1, MPFR_RNDN);
    mpfr_set   (u_im, b,    MPFR_RNDN);
    mpfr_ui_sub(v_re, 1, a, MPFR_RNDN);
    mpfr_neg   (v_im, b,    MPFR_RNDN);
    mpfr_complex_div(q_re, q_im, u_re, u_im, v_re, v_im);
    mpfr_complex_log(lr, li, q_re, q_im);
    mpfr_div_2ui(out_re, lr, 1, MPFR_RNDN);
    mpfr_div_2ui(out_im, li, 1, MPFR_RNDN);
    mpfr_clears(u_re, u_im, v_re, v_im, q_re, q_im, lr, li, (mpfr_ptr)NULL);
}

/* --------------------------------------------------------------------
 *  Reciprocal trig / hyperbolic composites
 * ------------------------------------------------------------------ */

/* out = 1 / forward(in) */
static void reciprocal_of(mpfr_t out_re, mpfr_t out_im,
                          MpfrComplexUnaryOp forward,
                          const mpfr_t in_re, const mpfr_t in_im) {
    mpfr_prec_t p = mpfr_get_prec(out_re);
    mpfr_t f_re, f_im, one, zero;
    mpfr_inits2(p, f_re, f_im, one, zero, (mpfr_ptr)NULL);
    forward(f_re, f_im, in_re, in_im);
    mpfr_set_ui(one, 1, MPFR_RNDN);
    mpfr_set_zero(zero, +1);
    mpfr_complex_div(out_re, out_im, one, zero, f_re, f_im);
    mpfr_clears(f_re, f_im, one, zero, (mpfr_ptr)NULL);
}

/* out = inverse_of_forward(1 / in) — for ArcCot/ArcSec/... */
static void inverse_of_reciprocal(mpfr_t out_re, mpfr_t out_im,
                                  MpfrComplexUnaryOp inverse_of_forward,
                                  const mpfr_t in_re, const mpfr_t in_im) {
    mpfr_prec_t p = mpfr_get_prec(out_re);
    mpfr_t inv_re, inv_im, one, zero;
    mpfr_inits2(p, inv_re, inv_im, one, zero, (mpfr_ptr)NULL);
    mpfr_set_ui(one, 1, MPFR_RNDN);
    mpfr_set_zero(zero, +1);
    mpfr_complex_div(inv_re, inv_im, one, zero, in_re, in_im);
    inverse_of_forward(out_re, out_im, inv_re, inv_im);
    mpfr_clears(inv_re, inv_im, one, zero, (mpfr_ptr)NULL);
}

void mpfr_complex_cot (mpfr_t r, mpfr_t i, const mpfr_t a, const mpfr_t b) { reciprocal_of(r, i, mpfr_complex_tan,  a, b); }
void mpfr_complex_sec (mpfr_t r, mpfr_t i, const mpfr_t a, const mpfr_t b) { reciprocal_of(r, i, mpfr_complex_cos,  a, b); }
void mpfr_complex_csc (mpfr_t r, mpfr_t i, const mpfr_t a, const mpfr_t b) { reciprocal_of(r, i, mpfr_complex_sin,  a, b); }
void mpfr_complex_coth(mpfr_t r, mpfr_t i, const mpfr_t a, const mpfr_t b) { reciprocal_of(r, i, mpfr_complex_tanh, a, b); }
void mpfr_complex_sech(mpfr_t r, mpfr_t i, const mpfr_t a, const mpfr_t b) { reciprocal_of(r, i, mpfr_complex_cosh, a, b); }
void mpfr_complex_csch(mpfr_t r, mpfr_t i, const mpfr_t a, const mpfr_t b) { reciprocal_of(r, i, mpfr_complex_sinh, a, b); }

void mpfr_complex_acot (mpfr_t r, mpfr_t i, const mpfr_t a, const mpfr_t b) { inverse_of_reciprocal(r, i, mpfr_complex_atan,  a, b); }
void mpfr_complex_asec (mpfr_t r, mpfr_t i, const mpfr_t a, const mpfr_t b) { inverse_of_reciprocal(r, i, mpfr_complex_acos,  a, b); }
void mpfr_complex_acsc (mpfr_t r, mpfr_t i, const mpfr_t a, const mpfr_t b) { inverse_of_reciprocal(r, i, mpfr_complex_asin,  a, b); }
void mpfr_complex_acoth(mpfr_t r, mpfr_t i, const mpfr_t a, const mpfr_t b) { inverse_of_reciprocal(r, i, mpfr_complex_atanh, a, b); }
void mpfr_complex_asech(mpfr_t r, mpfr_t i, const mpfr_t a, const mpfr_t b) { inverse_of_reciprocal(r, i, mpfr_complex_acosh, a, b); }
void mpfr_complex_acsch(mpfr_t r, mpfr_t i, const mpfr_t a, const mpfr_t b) { inverse_of_reciprocal(r, i, mpfr_complex_asinh, a, b); }

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

#endif /* USE_MPFR */
