/* test_integrate_linrad.c
 *
 * Tests for the linear-radicals rationalising substitution
 * (Integrate`LinearRadicals, Method -> "LinearRadicals"; added 2026-06-06).
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
        && strcmp(result->data.function.head->data.symbol, head) == 0,
        "expected unevaluated %s[...] for: %s", head, input);
    expr_free(result);
}

/* The Automatic cascade picks up linear radicals after the rational stage. */
static void test_automatic_cascade(void) {
    /* Single square-root radical of (x + 1). */
    assert_eval_eq(
        "Simplify[D[Integrate[1/Sqrt[x + 1], x], x] - 1/Sqrt[x + 1]]", "0", 0);
    /* Rational function of Sqrt[x]. */
    assert_eval_eq(
        "Simplify[D[Integrate[Sqrt[x]/(1 + Sqrt[x]), x], x]"
        " - Sqrt[x]/(1 + Sqrt[x])]", "0", 0);
    /* Bare Sqrt[x] (a = 1, b = 0). */
    assert_eval_eq(
        "Simplify[D[Integrate[Sqrt[x], x], x] - Sqrt[x]]", "0", 0);
}

/* Method plumbing: the option string and the explicit package head both route
 * to the routine and close the integrand. */
static void test_method_plumbing(void) {
    /* Cube-root radical: n = 3, mixed integer/log result. */
    assert_eval_eq(
        "Simplify[D[Integrate[1/(1 + x^(1/3)), x, Method -> \"LinearRadicals\"], x]"
        " - 1/(1 + x^(1/3))]", "0", 0);
    assert_eval_eq(
        "Simplify[D[Integrate`LinearRadicals[1/(1 + x^(1/3)), x], x]"
        " - 1/(1 + x^(1/3))]", "0", 0);
    /* Non-unit linear coefficients a = 2, b = 3. */
    assert_eval_eq(
        "Simplify[D[Integrate[Sqrt[2 x + 3]/x, x, Method -> \"LinearRadicals\"], x]"
        " - Sqrt[2 x + 3]/x]", "0", 0);
    /* Two different radical denominators share one base: n = LCM[2, 3] = 6. */
    assert_eval_eq(
        "Simplify[D[Integrate[1/(Sqrt[x] + x^(1/3)), x,"
        " Method -> \"LinearRadicals\"], x] - 1/(Sqrt[x] + x^(1/3))]", "0", 0);
}

/* Strict: the forced method returns unevaluated when the integrand is not a
 * rational function of radicals of a single linear argument. */
static void test_strict_no_match(void) {
    /* Two distinct linear bases (Union length > 1). */
    assert_head_unevaluated(
        "Integrate[Sqrt[x] + Sqrt[x + 1], x, Method -> \"LinearRadicals\"]",
        "Integrate");
    /* Radical of a non-linear argument. */
    assert_head_unevaluated(
        "Integrate`LinearRadicals[Sqrt[x^2 + 1], x]",
        "Integrate`LinearRadicals");
    /* No radical at all. */
    assert_head_unevaluated(
        "Integrate`LinearRadicals[Sin[x], x]",
        "Integrate`LinearRadicals");
    /* Free of x: nothing to substitute. */
    assert_head_unevaluated(
        "Integrate`LinearRadicals[5, x]",
        "Integrate`LinearRadicals");
}

void test_integrate_linrad(void) {
    symtab_init();
    core_init();

    TEST(test_automatic_cascade);
    TEST(test_method_plumbing);
    TEST(test_strict_no_match);

    printf("All Integrate LinearRadicals tests passed!\n");
}

int main(void) {
    test_integrate_linrad();
    return 0;
}
