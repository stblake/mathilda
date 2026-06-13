/* Tests for ND — numerical derivative.
 *
 * Cover: the default EulerSum method (forward-difference Richardson
 * extrapolation) for real, complex and non-analytic expressions; directional
 * derivatives via Scale (left/right/complex); the Scale and Terms options;
 * arbitrary precision via WorkingPrecision; the NIntegrate (Cauchy-integral
 * via NResidue) method including fractional and complex order; list-threading
 * over arg 0; the n = 0 case; option/argument-shape edge cases; the Protected
 * attribute; and memory hygiene.
 *
 * Expected values are taken from closed forms where available, otherwise from
 * Mathematica's NumericalCalculus`ND documentation. Numerical results are
 * compared *inside* the language (N[Abs[result - expected]] < tol) rather than
 * by parsing the printed form, which is rounded to ~6 significant figures.
 */

#include "core.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "test_utils.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char* eval_str(const char* input) {
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* e = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(e);
    expr_free(e);
    return s;
}

/* True if input is within `tol` of expected, evaluated at full internal
 * precision. Uses the L1 magnitude |Re d| + |Im d| of the difference: Abs of a
 * machine complex can carry a spurious tiny imaginary residual that leaves the
 * modulus complex, so a real comparand is built from Re/Im directly. */
static bool close_to(const char* input, const char* expected, double tol) {
    char buf[2048];
    snprintf(buf, sizeof buf,
             "N[Abs[Re[(%s) - (%s)]] + Abs[Im[(%s) - (%s)]]] < %.17g",
             input, expected, input, expected, tol);
    char* s = eval_str(buf);
    bool ok = (strcmp(s, "True") == 0);
    if (!ok) fprintf(stderr, "  close_to FAIL: |%s - (%s)| < %g  =>  %s\n",
                     input, expected, tol, s);
    free(s);
    return ok;
}

#define ASSERT_CLOSE(input, expected, tol)                                   \
    ASSERT_MSG(close_to((input), (expected), (tol)),                         \
               "%s ~= %s (tol %g)", (input), (expected), (tol))

/* ------------------------------------------------------------------------
 *  EulerSum (default) — elementary derivatives
 * ---------------------------------------------------------------------- */

static void test_basic_first_derivative(void) {
    ASSERT_CLOSE("ND[Exp[x],x,1]", "E", 1e-6);          /* d/dx e^x at 1 = e */
    ASSERT_CLOSE("ND[Sin[x],x,1]", "Cos[1]", 1e-6);
}

static void test_second_derivative(void) {
    ASSERT_CLOSE("ND[Cos[x]^3,{x,2},0]", "-3", 1e-5);   /* d2/dx2 cos^3 at 0 */
}

static void test_n_equals_zero(void) {
    ASSERT_CLOSE("ND[Exp[x],{x,0},2]", "Exp[2]", 1e-9); /* order 0 = f(x0) */
}

/* ------------------------------------------------------------------------
 *  Complex evaluation points
 * ---------------------------------------------------------------------- */

static void test_complex_point_real_result(void) {
    /* d/dx Sin[x] at Pi I = Cos[Pi I] = Cosh[Pi]. */
    ASSERT_CLOSE("ND[Sin[x],x,Pi I]", "Cosh[Pi]", 1e-4);
}

static void test_complex_point_complex_result(void) {
    /* d/dx Cos[I x] at 1+I = -I Sin[I(1+I)]. Doc: 0.634964 + 1.29846 I. */
    ASSERT_CLOSE("ND[Cos[I x],x,1+I]", "0.634964 + 1.29846 I", 1e-4);
}

/* ------------------------------------------------------------------------
 *  Non-analytic expressions — only EulerSum works
 * ---------------------------------------------------------------------- */

static void test_nonanalytic_eulersum(void) {
    /* Re[Cos[I y]] = Cosh[y]; derivative at 1 = Sinh[1]. EulerSum samples
     * along the real axis and is correct. */
    ASSERT_CLOSE("ND[Re[Cos[I y]],y,1]", "Sinh[1]", 1e-4);
}

/* ------------------------------------------------------------------------
 *  Directional derivatives via Scale
 * ---------------------------------------------------------------------- */

static void test_directional_abs(void) {
    ASSERT_CLOSE("ND[Abs[x],{x,1},0]", "1", 1e-9);            /* right deriv */
    ASSERT_CLOSE("ND[Abs[x],{x,1},0,Scale->-1]", "-1", 1e-9); /* left deriv */
}

static void test_directional_complex(void) {
    /* Complex direction: |(1+I)h|/((1+I)h) = (1 - I)/Sqrt[2]. */
    ASSERT_CLOSE("ND[Abs[x],{x,1},0,Scale->1+I]", "(1 - I)/Sqrt[2]", 1e-6);
}

/* ------------------------------------------------------------------------
 *  Scale captures the scale of variation; Terms controls accuracy
 * ---------------------------------------------------------------------- */

static void test_scale_variation(void) {
    /* d/dx Sin[100 x] at 0 = 100. The default Scale->1 mis-resolves the rapid
     * oscillation; Scale->1/100 (or more Terms) recovers it. */
    ASSERT_CLOSE("ND[Sin[100x],x,0,Scale->1/100]", "100", 1e-4);
    ASSERT_CLOSE("ND[Sin[100x],x,0,Terms->11]", "100", 1e-3);
}

/* ------------------------------------------------------------------------
 *  Arbitrary precision (MPFR) — fights subtractive cancellation
 * ---------------------------------------------------------------------- */

static void test_working_precision(void) {
    /* d^10/dx^10 Exp[x] at 0 = 1; machine EulerSum is swamped by cancellation,
     * higher WorkingPrecision + Terms recovers it. */
    ASSERT_CLOSE("ND[Exp[x],{x,10},0,WorkingPrecision->40,Terms->10]", "1", 1e-4);
    /* d^3/dx^3 Sin[x^2] at 1, to higher precision (doc value). */
    ASSERT_CLOSE("ND[Sin[x^2],{x,3},1,Terms->20,WorkingPrecision->40]",
                 "-14.420070264639875819038", 1e-9);
}

/* ------------------------------------------------------------------------
 *  NIntegrate method (Cauchy integral via NResidue)
 * ---------------------------------------------------------------------- */

static void test_nintegrate_analytic(void) {
    /* d^3/dx^3 Exp[x^2] at 1 = e^(x^2)(12x + 8x^3) | x=1 = 20 e. */
    ASSERT_CLOSE("ND[Exp[x^2],{x,3},1,Method->NIntegrate]", "20 E", 1e-4);
    /* d^4/dx^4 Exp[x^2] at 0 = 12. */
    ASSERT_CLOSE("ND[Exp[x^2],{x,4},0,Method->NIntegrate]", "12", 1e-6);
}

static void test_nintegrate_equals_factorial_nresidue(void) {
    /* ND[...,NIntegrate] == n! NResidue[expr/(x-x0)^(n+1), {x,x0}, Radius->1]. */
    ASSERT_CLOSE("ND[Exp[x^2],{x,4},0,Method->NIntegrate]",
                 "4! NResidue[Exp[x^2]/x^5,{x,0},Radius->1]", 1e-4);
}

static void test_nintegrate_fractional_order(void) {
    /* Half-integral of x at 1: Gamma generalization gives 4/(3 Sqrt[Pi]). */
    ASSERT_CLOSE("ND[x,{x,-1/2},1,Method->NIntegrate]", "4/(3 Sqrt[Pi])", 1e-4);
}

static void test_nintegrate_complex_order(void) {
    /* Complex order: doc value 0.0632028 + 1.11425 I. */
    ASSERT_CLOSE("ND[x^4,{x,I},1,Method->NIntegrate]", "0.0632028 + 1.11425 I", 1e-4);
}

/* ------------------------------------------------------------------------
 *  List-threading over arg 0
 * ---------------------------------------------------------------------- */

static void test_list_threading(void) {
    char* s = eval_str("ND[{Exp[x],Sin[x]},x,1]");
    ASSERT_MSG(s[0] == '{', "ND over a list should return a list: %s", s);
    free(s);
    ASSERT_CLOSE("ND[{Exp[x],Sin[x]},x,1][[1]]", "E", 1e-6);
    ASSERT_CLOSE("ND[{Exp[x],Sin[x]},x,1][[2]]", "Cos[1]", 1e-6);
}

/* ------------------------------------------------------------------------
 *  Edge cases / unevaluated forms
 * ---------------------------------------------------------------------- */

static void test_unevaluated_forms(void) {
    char* s1 = eval_str("ND[Exp[x],x,a]");      /* non-numeric point */
    ASSERT_MSG(strstr(s1, "ND[") != NULL, "symbolic x0 should stay unevaluated: %s", s1);
    free(s1);
    char* s2 = eval_str("ND[Exp[x],x]");        /* too few args */
    ASSERT_MSG(strstr(s2, "ND[") != NULL, "2-arg ND should stay unevaluated: %s", s2);
    free(s2);
    char* s3 = eval_str("ND[Exp[x],{1,2},0]");  /* non-symbol variable */
    ASSERT_MSG(strstr(s3, "ND[") != NULL, "bad spec should stay unevaluated: %s", s3);
    free(s3);
}

static void test_protected(void) {
    char* s = eval_str("MemberQ[Attributes[ND], Protected]");
    ASSERT_MSG(strcmp(s, "True") == 0, "ND should be Protected: %s", s);
    free(s);
}

static void test_memory_loop(void) {
    for (int i = 0; i < 25; i++) {
        const char* inputs[] = {
            "ND[Exp[x],x,1]",
            "ND[Sin[x^2],{x,2},1]",
            "ND[Exp[x^2],{x,3},1,Method->NIntegrate]",
            "ND[{Exp[x],Cos[x]},x,0]",
        };
        for (size_t k = 0; k < sizeof(inputs) / sizeof(inputs[0]); k++) {
            Expr* p = parse_expression(inputs[k]);
            ASSERT(p != NULL);
            Expr* v = evaluate(p);
            expr_free(p);
            expr_free(v);
        }
    }
}

/* ------------------------------------------------------------------------
 *  Main
 * ---------------------------------------------------------------------- */

int main(void) {
    symtab_init();
    core_init();

    TEST(test_basic_first_derivative);
    TEST(test_second_derivative);
    TEST(test_n_equals_zero);

    TEST(test_complex_point_real_result);
    TEST(test_complex_point_complex_result);

    TEST(test_nonanalytic_eulersum);

    TEST(test_directional_abs);
    TEST(test_directional_complex);

    TEST(test_scale_variation);

    TEST(test_working_precision);

    TEST(test_nintegrate_analytic);
    TEST(test_nintegrate_equals_factorial_nresidue);
    TEST(test_nintegrate_fractional_order);
    TEST(test_nintegrate_complex_order);

    TEST(test_list_threading);

    TEST(test_unevaluated_forms);
    TEST(test_protected);
    TEST(test_memory_loop);

    printf("All nderiv_tests passed.\n");
    return 0;
}
