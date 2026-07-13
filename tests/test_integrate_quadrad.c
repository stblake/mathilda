/* test_integrate_quadrad.c
 *
 * Tests for the quadratic-radicals (Euler) substitution
 * (Integrate`QuadraticRadicals, Method -> "QuadraticRadicals"; added 2026-06-06).
 * Correctness is asserted by the universal predicate
 * Simplify[D[Integrate[f, x], x] - f] === 0 rather than by fixed output
 * strings, so the tests survive surface-form changes.
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
 * (used for strict no-match cases under the forced method). */
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

/* The Automatic cascade picks up quadratic radicals after the linear-radical
 * stage.  Each case differentiates back to the integrand. */
static void test_automatic_cascade(void) {
    /* a > 0 (Euler I): square-root in the denominator. */
    assert_eval_eq(
        "Simplify[D[Integrate[1/Sqrt[x^2 + 1], x], x] - 1/Sqrt[x^2 + 1]]", "0", 0);
    /* a < 0, c > 0, disc > 0 (Euler III): the arcsine family. */
    assert_eval_eq(
        "Simplify[D[Integrate[1/Sqrt[1 - x^2], x], x] - 1/Sqrt[1 - x^2]]", "0", 0);
    /* a > 0, real roots: the arccosh family. */
    assert_eval_eq(
        "Simplify[D[Integrate[1/Sqrt[x^2 - 1], x], x] - 1/Sqrt[x^2 - 1]]", "0", 0);
    /* Radical in the numerator. */
    assert_eval_eq(
        "Simplify[D[Integrate[Sqrt[x^2 + 1], x], x] - Sqrt[x^2 + 1]]", "0", 0);
    /* General irreducible quadratic with a non-unit, non-square leading term. */
    assert_eval_eq(
        "Simplify[D[Integrate[1/Sqrt[2 x^2 + 3 x + 5], x], x]"
        " - 1/Sqrt[2 x^2 + 3 x + 5]]", "0", 0);
}

/* Method plumbing: the option string and the explicit package head both route
 * to the routine and close the integrand. */
static void test_method_plumbing(void) {
    assert_eval_eq(
        "Simplify[D[Integrate[1/Sqrt[x^2 + x + 1], x,"
        " Method -> \"QuadraticRadicals\"], x] - 1/Sqrt[x^2 + x + 1]]", "0", 0);
    assert_eval_eq(
        "Simplify[D[Integrate`QuadraticRadicals[1/Sqrt[x^2 + x + 1], x], x]"
        " - 1/Sqrt[x^2 + x + 1]]", "0", 0);
    /* Rational factor multiplying the radical (a > 0). */
    assert_eval_eq(
        "Simplify[D[Integrate[1/(x Sqrt[x^2 + x + 1]), x,"
        " Method -> \"QuadraticRadicals\"], x] - 1/(x Sqrt[x^2 + x + 1])]", "0", 0);
}

/* Symbolic leading coefficient: Euler I with a Sqrt[a] constant. */
static void test_symbolic_coefficient(void) {
    assert_eval_eq(
        "Simplify[D[Integrate[1/Sqrt[a x^2 + 1], x], x] - 1/Sqrt[a x^2 + 1]]",
        "0", 0);
}

/* Strict: the forced method / explicit head returns unevaluated when the
 * integrand is not a rational function of square roots of one quadratic. */
static void test_strict_no_match(void) {
    /* Two distinct (linear) radicands — no degree-2 radical at all. */
    assert_head_unevaluated(
        "Integrate[Sqrt[x] + Sqrt[x + 1], x, Method -> \"QuadraticRadicals\"]",
        "Integrate");
    /* Cube root of a quadratic (denominator != 2). */
    assert_head_unevaluated(
        "Integrate`QuadraticRadicals[(x^2 + 1)^(1/3), x]",
        "Integrate`QuadraticRadicals");
    /* Square root of a cubic (radicand not degree 2). */
    assert_head_unevaluated(
        "Integrate`QuadraticRadicals[Sqrt[x^3 + 1], x]",
        "Integrate`QuadraticRadicals");
    /* A linear radical belongs to LinearRadicals, not here. */
    assert_head_unevaluated(
        "Integrate`QuadraticRadicals[Sqrt[x + 1], x]",
        "Integrate`QuadraticRadicals");
    /* No radical at all. */
    assert_head_unevaluated(
        "Integrate`QuadraticRadicals[Sin[x], x]",
        "Integrate`QuadraticRadicals");
    /* Free of x: nothing to substitute. */
    assert_head_unevaluated(
        "Integrate`QuadraticRadicals[5, x]",
        "Integrate`QuadraticRadicals");
}

void test_integrate_quadrad(void) {
    symtab_init();
    core_init();

    TEST(test_automatic_cascade);
    TEST(test_method_plumbing);
    TEST(test_symbolic_coefficient);
    TEST(test_strict_no_match);

    printf("All Integrate QuadraticRadicals tests passed!\n");
}

int main(void) {
    test_integrate_quadrad();
    return 0;
}
