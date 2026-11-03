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

/* Symbolic leading coefficient (Nasser, 2026-11-03).  Sqrt[(a x + c)^2 - 1] used
 * to come back WITHOUT its Log term (the Euler substitution's Sqrt[a^2] vs
 * (a^2)^(3/2) hid a repeated factor from the rational integrator), a
 * b-scaled radicand gave a bogus RootSum, and (k x + c)/Sqrt[1 - (k x + c)^2]
 * took 110 s to fail (radicand degree read as 0).  Checked numerically: a
 * dropped Log is exactly what a symbolic Simplify can be fooled about. */
static void test_symbolic_leading_coefficient(void) {
    const char* fs[] = {
        "Sqrt[a^2 x^2 + 2 a c x + c^2 - 1]",
        "Sqrt[a^2 b^4 x^2 + 2 a b^4 c x + b^4 c^2 - 1]",
        "Sqrt[4 a^2 b^6 x^2 + 8 a b^6 c x + 4 b^6 c^2 - 4 b^2]/(2 b^2)",
        "x Sqrt[-b^2 (4 - 4 a^2 b^4 x^2 - 8 a b^4 x c - 4 b^4 c^2)]",
        "(k x + c)/Sqrt[1 - (k x + c)^2]",
        "x^2 Sqrt[a x^2 + 1]",
        "Sqrt[a x^2 + b x + c]",
    };
    for (size_t i = 0; i < sizeof fs / sizeof fs[0]; i++) {
        char buf[1024];
        snprintf(buf, sizeof buf,
            "Module[{g = Integrate[%s, x]}, FreeQ[g, Integrate | RootSum] && "
            "Abs[N[(D[g, x] - (%s)) /. {a -> 13/10, b -> 7/10, c -> 2/5, k -> 9/10, "
            "x -> 21/100}]] < 10^-10]", fs[i], fs[i]);
        assert_eval_eq(buf, "True", 0);
    }
    /* the cancelled-Log closed form, directly */
    assert_eval_eq("Integrate[(k x + c)/Sqrt[1 - (k x + c)^2], x]",
                   "-Sqrt[1 - (c + k x)^2]/k", 0);
    /* textbook shape of the symbolic-a closed form */
    assert_eval_eq("Integrate[Sqrt[a x^2 + b x + c], x]",
        "((b + 2 a x) Sqrt[c + b x + a x^2])/(4 a) + ((-b^2 + 4 a c) "
        "Log[b + 2 a x + 2 Sqrt[a] Sqrt[c + b x + a x^2]])/(8 a^(3/2))", 0);
}

void test_integrate_quadrad(void) {
    symtab_init();
    core_init();

    TEST(test_automatic_cascade);
    TEST(test_method_plumbing);
    TEST(test_symbolic_coefficient);
    TEST(test_strict_no_match);
    TEST(test_symbolic_leading_coefficient);

    printf("All Integrate QuadraticRadicals tests passed!\n");
}

int main(void) {
    test_integrate_quadrad();
    return 0;
}
