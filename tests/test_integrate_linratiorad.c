/* test_integrate_linratiorad.c
 *
 * Tests for the linear-ratio-radicals (Möbius) substitution
 * (Integrate`LinearRatioRadicals, Method -> "LinearRatioRadicals"; added
 * 2026-06-06).  Correctness is asserted by differentiating the antiderivative
 * back to the integrand.  Where the residual carries cube-/sixth-root radicals
 * that Simplify cannot reduce symbolically, a numeric check at a sample point is
 * used instead (the substitution is an exact bijection, so no symbolic
 * differentiate-back gate is run in the implementation either).
 */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

/* Assert that the evaluated `input` is an unevaluated call with head `head`
 * (used for strict no-match cases under the forced method / explicit head). */
static void assert_head_unevaluated(const char* input, const char* head) {
    Expr* parsed = parse_expression(input);
    Expr* result = evaluate(parsed);
    expr_free(parsed);
    ASSERT(result != NULL);
    ASSERT_MSG(result->type == EXPR_FUNCTION
        && result->data.function.head
        && result->data.function.head->type == EXPR_SYMBOL
        && strcmp(result->data.function.head->data.symbol.name, head) == 0,
        "expected unevaluated %s[...] for: %s", head, input);
    expr_free(result);
}

/* The Automatic cascade picks up linear-ratio radicals after the quadratic
 * stage.  Each case differentiates back to the integrand. */
static void test_automatic_cascade(void) {
    /* Sqrt of a Möbius argument in the denominator. */
    assert_eval_eq(
        "Simplify[D[Integrate[1/Sqrt[(x + 1)/(x - 1)], x], x]"
        " - 1/Sqrt[(x + 1)/(x - 1)]]", "0", 0);
    /* Radical in the numerator. */
    assert_eval_eq(
        "Simplify[D[Integrate[Sqrt[(x - 1)/(x + 1)], x], x]"
        " - Sqrt[(x - 1)/(x + 1)]]", "0", 0);
    /* General Möbius base with a, b, c, d all distinct. */
    assert_eval_eq(
        "Simplify[D[Integrate[1/Sqrt[(2 x + 1)/(x + 3)], x], x]"
        " - 1/Sqrt[(2 x + 1)/(x + 3)]]", "0", 0);
}

/* Method plumbing: the option string and the explicit package head both route
 * to the routine and close the integrand. */
static void test_method_plumbing(void) {
    assert_eval_eq(
        "Simplify[D[Integrate[1/Sqrt[(x + 1)/(x - 1)], x,"
        " Method -> \"LinearRatioRadicals\"], x] - 1/Sqrt[(x + 1)/(x - 1)]]",
        "0", 0);
    assert_eval_eq(
        "Simplify[D[Integrate`LinearRatioRadicals[1/Sqrt[(x + 1)/(x - 1)], x], x]"
        " - 1/Sqrt[(x + 1)/(x - 1)]]", "0", 0);
}

/* Cube root and an LCM-of-denominators case: the residual carries radicals
 * Simplify cannot collapse, so verify numerically at a sample point. */
static void test_cuberoot_and_lcm_numeric(void) {
    /* Cube root of a Möbius argument times a rational factor (n = 3). */
    assert_eval_eq(
        "Abs[N[(D[Integrate[((x - 1)/(x + 1))^(1/3)/x, x], x]"
        " - ((x - 1)/(x + 1))^(1/3)/x) /. x -> 2.7]] < 10^-10",
        "True", 0);
    /* Two radicals (1/2 and 1/3) sharing one base: n = LCM[2,3] = 6. */
    assert_eval_eq(
        "Abs[N[(D[Integrate[Sqrt[(x + 1)/(x - 1)] + ((x + 1)/(x - 1))^(1/3), x], x]"
        " - (Sqrt[(x + 1)/(x - 1)] + ((x + 1)/(x - 1))^(1/3))) /. x -> 4.2]]"
        " < 10^-9", "True", 0);
}

/* Inexact (numeric) coefficients: the integrand is rationalised, integrated
 * exactly and numericalised back; the derivative matches at a sample point. */
static void test_inexact_coefficients(void) {
    assert_eval_eq(
        "Abs[N[(D[Integrate[1/Sqrt[(2.0 x + 1)/(x - 1)], x], x]"
        " - 1/Sqrt[(2.0 x + 1)/(x - 1)]) /. x -> 3.3]] < 10^-9", "True", 0);
}

/* Strict: the forced method / explicit head returns unevaluated when the
 * integrand is not a rational function of radicals of one Möbius argument. */
static void test_strict_no_match(void) {
    /* A bare linear radical belongs to LinearRadicals, not here. */
    assert_head_unevaluated(
        "Integrate`LinearRatioRadicals[Sqrt[x + 1], x]",
        "Integrate`LinearRatioRadicals");
    /* A quadratic radical belongs to QuadraticRadicals. */
    assert_head_unevaluated(
        "Integrate`LinearRatioRadicals[Sqrt[x^2 + 1], x]",
        "Integrate`LinearRatioRadicals");
    /* Bare Sqrt[x] (denominator is constant). */
    assert_head_unevaluated(
        "Integrate`LinearRatioRadicals[Sqrt[x], x]",
        "Integrate`LinearRatioRadicals");
    /* Two distinct Möbius bases. */
    assert_head_unevaluated(
        "Integrate`LinearRatioRadicals[Sqrt[(x + 1)/(x - 1)]"
        " + Sqrt[(x + 2)/(x - 3)], x]",
        "Integrate`LinearRatioRadicals");
    /* No radical at all. */
    assert_head_unevaluated(
        "Integrate`LinearRatioRadicals[Sin[x], x]",
        "Integrate`LinearRatioRadicals");
    /* Free of x: nothing to substitute. */
    assert_head_unevaluated(
        "Integrate`LinearRatioRadicals[5, x]",
        "Integrate`LinearRatioRadicals");
}

/* The cascade must not let this method steal the linear / quadratic cases. */
static void test_cascade_partition(void) {
    /* Linear radical still routes through LinearRadicals. */
    assert_eval_eq("Integrate[1/Sqrt[x + 1], x]", "2 Sqrt[1 + x]", 0);
    /* Quadratic radical still closes through QuadraticRadicals. */
    assert_eval_eq(
        "Simplify[D[Integrate[1/Sqrt[x^2 + 1], x], x] - 1/Sqrt[x^2 + 1]]",
        "0", 0);
}

void test_integrate_linratiorad(void) {
    symtab_init();
    core_init();

    TEST(test_automatic_cascade);
    TEST(test_method_plumbing);
    TEST(test_cuberoot_and_lcm_numeric);
    TEST(test_inexact_coefficients);
    TEST(test_strict_no_match);
    TEST(test_cascade_partition);

    printf("All Integrate LinearRatioRadicals tests passed!\n");
}

int main(void) {
    test_integrate_linratiorad();
    return 0;
}
