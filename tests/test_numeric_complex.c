/* Mathilda — unit tests for src/numeric_complex.{c,h}.
 *
 * These call the MPFR-complex ops directly (rather than through any
 * builtin) so the helpers are exercised independently of the per-
 * function wirings added in later phases. Each test:
 *
 *   - allocates input + output MPFRs at PREC bits,
 *   - runs the op,
 *   - checks the result against an MPFR oracle computed inline from the
 *     same identity, within 4 ulp at PREC bits (TOL_BITS below).
 *
 * The 4-ulp tolerance absorbs the rounding cost of going through 3-4
 * mpfr_* calls in the op vs the oracle's single closed-form arithmetic.
 */

#include "test_utils.h"

#ifdef USE_MPFR

#include "numeric_complex.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "common.h"

#include <mpfr.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

extern void symtab_init(void);
extern void core_init(void);
extern void files_init(void);

#define PREC 128
#define TOL_BITS 124   /* allow ~4 ulp = drop the bottom 4 bits */

/* Check two MPFRs agree to TOL_BITS bits of relative accuracy:
 *   |a - b| <= 2^-TOL_BITS * max(|a|, |b|, 1)
 * Falls back to absolute when both are near zero. */
static void assert_mpfr_close(const mpfr_t got, const mpfr_t expect,
                              const char* label) {
    mpfr_t diff, scale, tol;
    mpfr_inits2(PREC, diff, scale, tol, (mpfr_ptr)NULL);
    mpfr_sub(diff, got, expect, MPFR_RNDN);
    mpfr_abs(diff, diff, MPFR_RNDN);

    mpfr_t ag, ae;
    mpfr_inits2(PREC, ag, ae, (mpfr_ptr)NULL);
    mpfr_abs(ag, got,    MPFR_RNDN);
    mpfr_abs(ae, expect, MPFR_RNDN);
    mpfr_max(scale, ag, ae, MPFR_RNDN);
    if (mpfr_cmp_ui(scale, 1) < 0) mpfr_set_ui(scale, 1, MPFR_RNDN);

    mpfr_set_ui(tol, 1, MPFR_RNDN);
    mpfr_div_2ui(tol, tol, TOL_BITS, MPFR_RNDN);
    mpfr_mul(tol, tol, scale, MPFR_RNDN);

    if (mpfr_cmp(diff, tol) > 0) {
        mpfr_printf("FAIL: %s\n  got    = %.40Re\n  expect = %.40Re\n  diff   = %.10Re\n  tol    = %.10Re\n",
                    label, got, expect, diff, tol);
        mpfr_clears(diff, scale, tol, ag, ae, (mpfr_ptr)NULL);
        exit(1);
    }
    mpfr_clears(diff, scale, tol, ag, ae, (mpfr_ptr)NULL);
}

/* Set an mpfr_t to a fresh constant value (decimal string), at PREC. */
static void mset(mpfr_t x, const char* s) {
    mpfr_init2(x, PREC);
    mpfr_set_str(x, s, 10, MPFR_RNDN);
}

/* --------------------------------------------------------------------
 *  Construction
 * ------------------------------------------------------------------ */

static void test_make_complex_collapses_on_zero_im(void) {
    mpfr_t re, im;
    mset(re, "3.5");
    mset(im, "0");
    Expr* r = numeric_mpfr_make_complex(re, im);
    ASSERT(r != NULL);
    ASSERT(r->type == EXPR_MPFR);
    ASSERT(mpfr_cmp_d(r->data.mpfr, 3.5) == 0);
    expr_free(r);
    mpfr_clears(re, im, (mpfr_ptr)NULL);
}

static void test_make_complex_keeps_complex_when_im_nonzero(void) {
    mpfr_t re, im;
    mset(re, "3.5");
    mset(im, "1.25");
    Expr* r = numeric_mpfr_make_complex(re, im);
    ASSERT(r != NULL);
    ASSERT(r->type == EXPR_FUNCTION);
    ASSERT(r->data.function.arg_count == 2);
    ASSERT(r->data.function.args[0]->type == EXPR_MPFR);
    ASSERT(r->data.function.args[1]->type == EXPR_MPFR);
    ASSERT(mpfr_cmp_d(r->data.function.args[0]->data.mpfr, 3.5)  == 0);
    ASSERT(mpfr_cmp_d(r->data.function.args[1]->data.mpfr, 1.25) == 0);
    expr_free(r);
    mpfr_clears(re, im, (mpfr_ptr)NULL);
}

/* --------------------------------------------------------------------
 *  Forward transcendentals
 * ------------------------------------------------------------------ */

static void test_exp_euler_identity(void) {
    /* exp(i pi) == -1 + 0 i */
    mpfr_t a, b, out_re, out_im, expect_re, expect_im;
    mpfr_inits2(PREC, a, out_re, out_im, expect_re, expect_im, (mpfr_ptr)NULL);
    mpfr_init2(b, PREC);
    mpfr_set_zero(a, +1);
    mpfr_const_pi(b, MPFR_RNDN);
    mpfr_complex_exp(out_re, out_im, a, b);
    mpfr_set_si(expect_re, -1, MPFR_RNDN);
    mpfr_set_zero(expect_im, +1);
    assert_mpfr_close(out_re, expect_re, "exp(i pi).re");
    /* imag is ~0; allow the same absolute tolerance via the scale floor */
    assert_mpfr_close(out_im, expect_im, "exp(i pi).im");
    mpfr_clears(a, b, out_re, out_im, expect_re, expect_im, (mpfr_ptr)NULL);
}

static void test_log_inverse_of_exp(void) {
    /* log(exp(z)) == z (mod 2 pi in imag part); test a + b in (-pi, pi]. */
    mpfr_t a, b, e_re, e_im, l_re, l_im;
    mpfr_inits2(PREC, a, b, e_re, e_im, l_re, l_im, (mpfr_ptr)NULL);
    mpfr_set_d(a, 1.5,  MPFR_RNDN);
    mpfr_set_d(b, 0.75, MPFR_RNDN);
    mpfr_complex_exp(e_re, e_im, a, b);
    mpfr_complex_log(l_re, l_im, e_re, e_im);
    assert_mpfr_close(l_re, a, "log(exp(z)).re");
    assert_mpfr_close(l_im, b, "log(exp(z)).im");
    mpfr_clears(a, b, e_re, e_im, l_re, l_im, (mpfr_ptr)NULL);
}

static void test_sin_cos_identity(void) {
    /* sin^2(z) + cos^2(z) == 1 + 0 i */
    mpfr_t a, b, s_re, s_im, c_re, c_im;
    mpfr_t s2_re, s2_im, c2_re, c2_im, sum_re, sum_im;
    mpfr_inits2(PREC, a, b, s_re, s_im, c_re, c_im,
                      s2_re, s2_im, c2_re, c2_im, sum_re, sum_im,
                      (mpfr_ptr)NULL);
    mpfr_set_d(a, 0.7, MPFR_RNDN);
    mpfr_set_d(b, 1.3, MPFR_RNDN);
    mpfr_complex_sin(s_re, s_im, a, b);
    mpfr_complex_cos(c_re, c_im, a, b);
    /* s^2 = (s_re^2 - s_im^2) + (2 s_re s_im) i  (since (a+bi)^2 = a^2-b^2 + 2abi) */
    mpfr_sqr(s2_re, s_re, MPFR_RNDN);
    mpfr_sqr(sum_re, s_im, MPFR_RNDN);
    mpfr_sub(s2_re, s2_re, sum_re, MPFR_RNDN);
    mpfr_mul(s2_im, s_re, s_im, MPFR_RNDN);
    mpfr_mul_2ui(s2_im, s2_im, 1, MPFR_RNDN);
    mpfr_sqr(c2_re, c_re, MPFR_RNDN);
    mpfr_sqr(sum_re, c_im, MPFR_RNDN);
    mpfr_sub(c2_re, c2_re, sum_re, MPFR_RNDN);
    mpfr_mul(c2_im, c_re, c_im, MPFR_RNDN);
    mpfr_mul_2ui(c2_im, c2_im, 1, MPFR_RNDN);
    mpfr_add(sum_re, s2_re, c2_re, MPFR_RNDN);
    mpfr_add(sum_im, s2_im, c2_im, MPFR_RNDN);
    mpfr_t one, zero;
    mpfr_inits2(PREC, one, zero, (mpfr_ptr)NULL);
    mpfr_set_ui(one, 1, MPFR_RNDN);
    mpfr_set_zero(zero, +1);
    assert_mpfr_close(sum_re, one,  "sin^2 + cos^2 .re");
    assert_mpfr_close(sum_im, zero, "sin^2 + cos^2 .im");
    mpfr_clears(a, b, s_re, s_im, c_re, c_im,
                s2_re, s2_im, c2_re, c2_im, sum_re, sum_im, one, zero,
                (mpfr_ptr)NULL);
}

static void test_tan_matches_sin_over_cos(void) {
    /* tan(z) == sin(z) / cos(z) */
    mpfr_t a, b;
    mpfr_t s_re, s_im, c_re, c_im, t_re, t_im;
    mpfr_t q_re, q_im;
    mpfr_inits2(PREC, a, b, s_re, s_im, c_re, c_im, t_re, t_im,
                      q_re, q_im, (mpfr_ptr)NULL);
    mpfr_set_d(a, 0.8, MPFR_RNDN);
    mpfr_set_d(b, 0.4, MPFR_RNDN);
    mpfr_complex_sin(s_re, s_im, a, b);
    mpfr_complex_cos(c_re, c_im, a, b);
    mpfr_complex_div(q_re, q_im, s_re, s_im, c_re, c_im);
    mpfr_complex_tan(t_re, t_im, a, b);
    assert_mpfr_close(t_re, q_re, "tan vs sin/cos .re");
    assert_mpfr_close(t_im, q_im, "tan vs sin/cos .im");
    mpfr_clears(a, b, s_re, s_im, c_re, c_im, t_re, t_im, q_re, q_im,
                (mpfr_ptr)NULL);
}

static void test_sinh_cosh_identity(void) {
    /* cosh^2 - sinh^2 == 1 (complex form: same algebraic identity) */
    mpfr_t a, b;
    mpfr_t s_re, s_im, c_re, c_im;
    mpfr_t s2_re, s2_im, c2_re, c2_im, diff_re, diff_im, t;
    mpfr_inits2(PREC, a, b, s_re, s_im, c_re, c_im,
                      s2_re, s2_im, c2_re, c2_im, diff_re, diff_im, t,
                      (mpfr_ptr)NULL);
    mpfr_set_d(a, 0.5, MPFR_RNDN);
    mpfr_set_d(b, 0.9, MPFR_RNDN);
    mpfr_complex_sinh(s_re, s_im, a, b);
    mpfr_complex_cosh(c_re, c_im, a, b);
    /* cosh^2 = (c_re^2 - c_im^2) + 2 c_re c_im i */
    mpfr_sqr(c2_re, c_re, MPFR_RNDN);
    mpfr_sqr(t,     c_im, MPFR_RNDN);
    mpfr_sub(c2_re, c2_re, t, MPFR_RNDN);
    mpfr_mul(c2_im, c_re, c_im, MPFR_RNDN);
    mpfr_mul_2ui(c2_im, c2_im, 1, MPFR_RNDN);
    /* sinh^2 = (s_re^2 - s_im^2) + 2 s_re s_im i */
    mpfr_sqr(s2_re, s_re, MPFR_RNDN);
    mpfr_sqr(t,     s_im, MPFR_RNDN);
    mpfr_sub(s2_re, s2_re, t, MPFR_RNDN);
    mpfr_mul(s2_im, s_re, s_im, MPFR_RNDN);
    mpfr_mul_2ui(s2_im, s2_im, 1, MPFR_RNDN);
    mpfr_sub(diff_re, c2_re, s2_re, MPFR_RNDN);
    mpfr_sub(diff_im, c2_im, s2_im, MPFR_RNDN);
    mpfr_t one, zero;
    mpfr_inits2(PREC, one, zero, (mpfr_ptr)NULL);
    mpfr_set_ui(one, 1, MPFR_RNDN);
    mpfr_set_zero(zero, +1);
    assert_mpfr_close(diff_re, one,  "cosh^2 - sinh^2 .re");
    assert_mpfr_close(diff_im, zero, "cosh^2 - sinh^2 .im");
    mpfr_clears(a, b, s_re, s_im, c_re, c_im, s2_re, s2_im, c2_re, c2_im,
                diff_re, diff_im, t, one, zero, (mpfr_ptr)NULL);
}

static void test_tanh_matches_sinh_over_cosh(void) {
    mpfr_t a, b, s_re, s_im, c_re, c_im, t_re, t_im, q_re, q_im;
    mpfr_inits2(PREC, a, b, s_re, s_im, c_re, c_im, t_re, t_im, q_re, q_im, (mpfr_ptr)NULL);
    mpfr_set_d(a, 0.4,  MPFR_RNDN);
    mpfr_set_d(b, 0.65, MPFR_RNDN);
    mpfr_complex_sinh(s_re, s_im, a, b);
    mpfr_complex_cosh(c_re, c_im, a, b);
    mpfr_complex_div(q_re, q_im, s_re, s_im, c_re, c_im);
    mpfr_complex_tanh(t_re, t_im, a, b);
    assert_mpfr_close(t_re, q_re, "tanh vs sinh/cosh .re");
    assert_mpfr_close(t_im, q_im, "tanh vs sinh/cosh .im");
    mpfr_clears(a, b, s_re, s_im, c_re, c_im, t_re, t_im, q_re, q_im, (mpfr_ptr)NULL);
}

/* --------------------------------------------------------------------
 *  Inverse transcendentals
 * ------------------------------------------------------------------ */

static void test_asin_inverse_of_sin(void) {
    /* asin(sin(z)) == z for z with |re| <= pi/2 (principal branch). */
    mpfr_t a, b, s_re, s_im, r_re, r_im;
    mpfr_inits2(PREC, a, b, s_re, s_im, r_re, r_im, (mpfr_ptr)NULL);
    mpfr_set_d(a, 0.6, MPFR_RNDN);
    mpfr_set_d(b, 0.4, MPFR_RNDN);
    mpfr_complex_sin (s_re, s_im, a, b);
    mpfr_complex_asin(r_re, r_im, s_re, s_im);
    assert_mpfr_close(r_re, a, "asin(sin(z)).re");
    assert_mpfr_close(r_im, b, "asin(sin(z)).im");
    mpfr_clears(a, b, s_re, s_im, r_re, r_im, (mpfr_ptr)NULL);
}

static void test_atan_inverse_of_tan(void) {
    mpfr_t a, b, t_re, t_im, r_re, r_im;
    mpfr_inits2(PREC, a, b, t_re, t_im, r_re, r_im, (mpfr_ptr)NULL);
    mpfr_set_d(a, 0.3, MPFR_RNDN);
    mpfr_set_d(b, 0.5, MPFR_RNDN);
    mpfr_complex_tan (t_re, t_im, a, b);
    mpfr_complex_atan(r_re, r_im, t_re, t_im);
    assert_mpfr_close(r_re, a, "atan(tan(z)).re");
    assert_mpfr_close(r_im, b, "atan(tan(z)).im");
    mpfr_clears(a, b, t_re, t_im, r_re, r_im, (mpfr_ptr)NULL);
}

static void test_asinh_inverse_of_sinh(void) {
    mpfr_t a, b, s_re, s_im, r_re, r_im;
    mpfr_inits2(PREC, a, b, s_re, s_im, r_re, r_im, (mpfr_ptr)NULL);
    mpfr_set_d(a, 0.4, MPFR_RNDN);
    mpfr_set_d(b, 0.7, MPFR_RNDN);
    mpfr_complex_sinh (s_re, s_im, a, b);
    mpfr_complex_asinh(r_re, r_im, s_re, s_im);
    assert_mpfr_close(r_re, a, "asinh(sinh(z)).re");
    assert_mpfr_close(r_im, b, "asinh(sinh(z)).im");
    mpfr_clears(a, b, s_re, s_im, r_re, r_im, (mpfr_ptr)NULL);
}

static void test_atanh_inverse_of_tanh(void) {
    mpfr_t a, b, t_re, t_im, r_re, r_im;
    mpfr_inits2(PREC, a, b, t_re, t_im, r_re, r_im, (mpfr_ptr)NULL);
    mpfr_set_d(a, 0.25, MPFR_RNDN);
    mpfr_set_d(b, 0.5,  MPFR_RNDN);
    mpfr_complex_tanh (t_re, t_im, a, b);
    mpfr_complex_atanh(r_re, r_im, t_re, t_im);
    assert_mpfr_close(r_re, a, "atanh(tanh(z)).re");
    assert_mpfr_close(r_im, b, "atanh(tanh(z)).im");
    mpfr_clears(a, b, t_re, t_im, r_re, r_im, (mpfr_ptr)NULL);
}

static void test_acos_pi_half_minus_asin(void) {
    /* acos(z) + asin(z) == pi/2 + 0 i */
    mpfr_t a, b, ac_re, ac_im, as_re, as_im, sum_re, sum_im, pi_half, zero;
    mpfr_inits2(PREC, a, b, ac_re, ac_im, as_re, as_im, sum_re, sum_im,
                      pi_half, zero, (mpfr_ptr)NULL);
    mpfr_set_d(a, 0.4, MPFR_RNDN);
    mpfr_set_d(b, 0.7, MPFR_RNDN);
    mpfr_complex_acos(ac_re, ac_im, a, b);
    mpfr_complex_asin(as_re, as_im, a, b);
    mpfr_add(sum_re, ac_re, as_re, MPFR_RNDN);
    mpfr_add(sum_im, ac_im, as_im, MPFR_RNDN);
    mpfr_const_pi(pi_half, MPFR_RNDN);
    mpfr_div_2ui(pi_half, pi_half, 1, MPFR_RNDN);
    mpfr_set_zero(zero, +1);
    assert_mpfr_close(sum_re, pi_half, "acos + asin .re");
    assert_mpfr_close(sum_im, zero,    "acos + asin .im");
    mpfr_clears(a, b, ac_re, ac_im, as_re, as_im, sum_re, sum_im,
                pi_half, zero, (mpfr_ptr)NULL);
}

static void test_acosh_branch(void) {
    /* acosh(2 + 0 i) == log(2 + sqrt(3)) + 0 i (positive real). */
    mpfr_t a, b, r_re, r_im, expect_re, three, two, sq3, sum, zero;
    mpfr_inits2(PREC, a, b, r_re, r_im, expect_re, three, two, sq3, sum, zero,
                (mpfr_ptr)NULL);
    mpfr_set_ui(a, 2, MPFR_RNDN);
    mpfr_set_zero(b, +1);
    mpfr_complex_acosh(r_re, r_im, a, b);
    mpfr_set_ui(three, 3, MPFR_RNDN);
    mpfr_sqrt(sq3, three, MPFR_RNDN);
    mpfr_set_ui(two, 2, MPFR_RNDN);
    mpfr_add(sum, two, sq3, MPFR_RNDN);
    mpfr_log(expect_re, sum, MPFR_RNDN);
    mpfr_set_zero(zero, +1);
    assert_mpfr_close(r_re, expect_re, "acosh(2).re");
    assert_mpfr_close(r_im, zero,      "acosh(2).im");
    mpfr_clears(a, b, r_re, r_im, expect_re, three, two, sq3, sum, zero,
                (mpfr_ptr)NULL);
}

/* --------------------------------------------------------------------
 *  Dispatcher
 * ------------------------------------------------------------------ */

static void test_apply_complex_unary_through_expr(void) {
    /* Drive the dispatcher with a plain integer (which get_approx_mpfr
     * widens to (3, 0) at PREC bits) and an exp op. exp(3 + 0 i) == e^3,
     * which collapses to pure EXPR_MPFR. */
    Expr* in = expr_new_integer(3);
    Expr* out = numeric_mpfr_apply_complex_unary(in, PREC, mpfr_complex_exp);
    ASSERT(out != NULL);
    ASSERT(out->type == EXPR_MPFR);
    mpfr_t expect, three;
    mpfr_inits2(PREC, expect, three, (mpfr_ptr)NULL);
    mpfr_set_ui(three, 3, MPFR_RNDN);
    mpfr_exp(expect, three, MPFR_RNDN);
    assert_mpfr_close(out->data.mpfr, expect, "apply_complex_unary(exp, 3)");
    expr_free(in);
    expr_free(out);
    mpfr_clears(expect, three, (mpfr_ptr)NULL);
}

static void test_apply_complex_unary_keeps_complex_when_im_nonzero(void) {
    /* Sin of Complex[1, 1] at PREC: result has nonzero imag part, so the
     * dispatcher must return Complex[EXPR_MPFR, EXPR_MPFR]. */
    Expr* re = expr_new_mpfr_from_d(1.0, PREC);
    Expr* im = expr_new_mpfr_from_d(1.0, PREC);
    Expr* args[2] = { re, im };
    Expr* in = expr_new_function(expr_new_symbol("Complex"), args, 2);
    Expr* out = numeric_mpfr_apply_complex_unary(in, PREC, mpfr_complex_sin);
    ASSERT(out != NULL);
    ASSERT(out->type == EXPR_FUNCTION);
    ASSERT(out->data.function.arg_count == 2);
    ASSERT(out->data.function.args[0]->type == EXPR_MPFR);
    ASSERT(out->data.function.args[1]->type == EXPR_MPFR);
    expr_free(in);
    expr_free(out);
}

/* --------------------------------------------------------------------
 *  Reciprocals
 * ------------------------------------------------------------------ */

static void test_csc_is_one_over_sin(void) {
    mpfr_t a, b, s_re, s_im, csc_re, csc_im, one, zero, inv_re, inv_im;
    mpfr_inits2(PREC, a, b, s_re, s_im, csc_re, csc_im, one, zero, inv_re, inv_im,
                (mpfr_ptr)NULL);
    mpfr_set_d(a, 0.5, MPFR_RNDN);
    mpfr_set_d(b, 0.3, MPFR_RNDN);
    mpfr_complex_sin(s_re, s_im, a, b);
    mpfr_complex_csc(csc_re, csc_im, a, b);
    mpfr_set_ui(one, 1, MPFR_RNDN);
    mpfr_set_zero(zero, +1);
    mpfr_complex_div(inv_re, inv_im, one, zero, s_re, s_im);
    assert_mpfr_close(csc_re, inv_re, "csc vs 1/sin .re");
    assert_mpfr_close(csc_im, inv_im, "csc vs 1/sin .im");
    mpfr_clears(a, b, s_re, s_im, csc_re, csc_im, one, zero, inv_re, inv_im,
                (mpfr_ptr)NULL);
}

/* --------------------------------------------------------------------
 *  main
 * ------------------------------------------------------------------ */

int main(void) {
    symtab_init();
    core_init();
    files_init();

    TEST(test_make_complex_collapses_on_zero_im);
    TEST(test_make_complex_keeps_complex_when_im_nonzero);
    TEST(test_exp_euler_identity);
    TEST(test_log_inverse_of_exp);
    TEST(test_sin_cos_identity);
    TEST(test_tan_matches_sin_over_cos);
    TEST(test_sinh_cosh_identity);
    TEST(test_tanh_matches_sinh_over_cosh);
    TEST(test_asin_inverse_of_sin);
    TEST(test_atan_inverse_of_tan);
    TEST(test_asinh_inverse_of_sinh);
    TEST(test_atanh_inverse_of_tanh);
    TEST(test_acos_pi_half_minus_asin);
    TEST(test_acosh_branch);
    TEST(test_apply_complex_unary_through_expr);
    TEST(test_apply_complex_unary_keeps_complex_when_im_nonzero);
    TEST(test_csc_is_one_over_sin);
    printf("All numeric_complex tests passed.\n");
    return 0;
}

#else  /* USE_MPFR not defined */

int main(void) {
    printf("USE_MPFR not enabled; numeric_complex tests are no-ops.\n");
    return 0;
}

#endif
