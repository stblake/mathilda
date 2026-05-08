/* test_intrat.c — Phase 1 of the IntegrateRational port.
 *
 * Coverage:
 *   * Helpers: Content / Primitive / Monic / LeadingCoefficient.
 *   * IntegratePolynomial — term-by-term integration over Q[x].
 *   * HermiteReduce — Mack's linear version, on the four reference
 *     cases lifted from IntegrateRational.m:1331-1358 plus a couple
 *     of single-factor sanity checks.
 *   * Integrate dispatcher — polynomial / 1/(x-a) / 1/(x-a)^2 /
 *     c·D[d]/d^k cases that Phase 1 closes via derivative
 *     recognition.
 *
 * The universal correctness predicate is differentiation:
 *   ASSERT_INTEGRAL_OK(f, x)  ≡  Cancel[Together[D[Integrate[f,x],x]-f]] === 0
 * which catches algorithmic bugs even when the printed form differs
 * from the Mathematica reference.  HermiteReduce uses the analogous
 * `Cancel[Together[D[g, x] + h - f]] === 0` check.
 */

#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* Test helpers                                                        */
/* ------------------------------------------------------------------ */

/* Compare evaluated output to a string. */
static void run_eq(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* got = expr_to_string(res);
    if (strcmp(got, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n",
               input, expected, got);
        ASSERT_STR_EQ(got, expected);
    }
    free(got);
    expr_free(e);
    expr_free(res);
}

/* Check that ∫ integrand dx differentiates back to the integrand.
 * We Expand both sides before canceling so structurally-different
 * polynomial coefficients (e.g. Plus distributed differently) collapse
 * to the same canonical form. */
static void assert_integral_correct(const char* integrand) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
        "Cancel[Together[Expand[D[Integrate[%s, x], x] - (%s)]]]",
        integrand, integrand);
    Expr* e = parse_expression(buf);
    Expr* res = evaluate(e);
    char* got = expr_to_string(res);
    if (strcmp(got, "0") != 0) {
        printf("FAIL: D[Integrate[%s, x], x] - %s != 0\n  Got: %s\n",
               integrand, integrand, got);
        ASSERT_STR_EQ(got, "0");
    }
    free(got);
    expr_free(e);
    expr_free(res);
}

/* Check HermiteReduce[f, x] = {g, h} satisfies D[g, x] + h - f == 0. */
static void assert_hermite_correct(const char* integrand) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
        "With[{hr = Integrate`HermiteReduce[%s, x]}, "
        "Cancel[Together[D[hr[[1]], x] + hr[[2]] - (%s)]]]",
        integrand, integrand);
    Expr* e = parse_expression(buf);
    Expr* res = evaluate(e);
    char* got = expr_to_string(res);
    if (strcmp(got, "0") != 0) {
        printf("FAIL: HermiteReduce decomposition incorrect for %s\n  Got: %s\n",
               integrand, got);
        ASSERT_STR_EQ(got, "0");
    }
    free(got);
    expr_free(e);
    expr_free(res);
}

/* ------------------------------------------------------------------ */
/* Helpers                                                              */
/* ------------------------------------------------------------------ */

static void test_helpers_content(void) {
    run_eq("Integrate`Helpers`Content[3 + 6 x + 9 x^2, x]", "3");
    run_eq("Integrate`Helpers`Content[2 x^2 + 4 x, x]", "2");
    run_eq("Integrate`Helpers`Content[x^3 + x, x]", "1");
}

static void test_helpers_primitive(void) {
    run_eq("Integrate`Helpers`Primitive[3 + 6 x + 9 x^2, x]",
           "1 + 2 x + 3 x^2");
    run_eq("Integrate`Helpers`Primitive[2 x^2 + 4 x, x]",
           "2 x + x^2");
}

static void test_helpers_leading_coefficient(void) {
    run_eq("Integrate`Helpers`LeadingCoefficient[3 x^2 + 1, x]", "3");
    run_eq("Integrate`Helpers`LeadingCoefficient[5, x]", "5");
    run_eq("Integrate`Helpers`LeadingCoefficient[x^4 - x + 2, x]", "1");
}

static void test_helpers_monic(void) {
    run_eq("Integrate`Helpers`Monic[2 x^2 + 4 x + 6, x]",
           "3 + 2 x + x^2");
    /* Polynomial free of x ⇒ 1 (Mathematica convention). */
    run_eq("Integrate`Helpers`Monic[7, x]", "1");
}

/* ------------------------------------------------------------------ */
/* IntegratePolynomial                                                 */
/* ------------------------------------------------------------------ */

static void test_integrate_polynomial(void) {
    run_eq("Integrate`IntegratePolynomial[0, x]", "0");
    run_eq("Integrate`IntegratePolynomial[1, x]", "x");
    run_eq("Integrate`IntegratePolynomial[x, x]", "1/2 x^2");
    run_eq("Integrate`IntegratePolynomial[x^2, x]", "1/3 x^3");
    /* Linear combination. */
    run_eq("Integrate`IntegratePolynomial[3 + 5 x + 2 x^2, x]",
           "3 x + 5/2 x^2 + 2/3 x^3");
    /* Higher-order with parameter — Phase 1 leaves the coefficient
     * unchanged because IntegratePolynomial assumes a polynomial in x.
     * The y^2-3y+1 factor is treated as a constant. */
    assert_integral_correct("(y^2 - 3 y + 1)*(1 - x^2)^3");
}

/* ------------------------------------------------------------------ */
/* HermiteReduce                                                       */
/* ------------------------------------------------------------------ */

static void test_hermite_simple(void) {
    /* deg(D̄) starts at 0 ⇒ f stays as-is, g == 0. */
    run_eq("Integrate`HermiteReduce[1/x, x]", "{0, 1/x}");
    run_eq("Integrate`HermiteReduce[1/(x-a), x]", "{0, 1/(-a + x)}");

    /* 1/(x-a)^2: classical reduction down to g = -1/(x-a), h = 0. */
    assert_hermite_correct("1/(x-1)^2");
    assert_hermite_correct("1/(x-a)^2");

    /* Two repeated factors: the algorithm iterates twice. */
    assert_hermite_correct("1/(x^2 (x-1)^3)");
    assert_hermite_correct("(2 x + 3)/(x^2 (x+1)^4)");

    /* Squarefull denominator with an integer-content factor — caught
     * the scale-invariance bug in the gcd choice during Phase 1.
     * We use the monic gcd to keep the algorithm correct. */
    assert_hermite_correct("1/(4 x^2 (5 - 4 x)^2)");
}

static void test_hermite_reference(void) {
    /* The four cases lifted from IntegrateRational.m:1331-1358.
     * Verification is by differentiation. */
    assert_hermite_correct("(x^7 - 24 x^4 - 4 x^2 + 8 x - 8)/(x^8 + 6 x^6 + 12 x^4 + 8 x^2)");
    assert_hermite_correct("1/(x^4 + 1)^8");
    assert_hermite_correct("(-4 + x - 9 x^2)/(10 x - 8 x^2)^3");
    assert_hermite_correct("(28 + 110 x + 147 x^2 - 1131 x^3 - 945 x^4 - 189 x^5)/"
                            "(16 - 32 x + 100 x^2 - 56 x^3 + 296 x^4 + 216 x^5 + 36 x^6)");
}

/* ------------------------------------------------------------------ */
/* Integrate dispatcher                                                */
/* ------------------------------------------------------------------ */

static void test_integrate_polynomial_dispatch(void) {
    run_eq("Integrate[x^2, x]", "1/3 x^3");
    run_eq("Integrate[3 + 5 x + 2 x^2, x]", "3 x + 5/2 x^2 + 2/3 x^3");
    run_eq("Integrate[a x^3 + b, x]", "1/4 a x^4 + b x");
    /* Pure constant. */
    run_eq("Integrate[5, x]", "5 x");
}

static void test_integrate_log_recognition(void) {
    /* num is exactly D[den]/den ⇒ Log[den]. */
    run_eq("Integrate[1/x, x]", "Log[x]");
    run_eq("Integrate[2 x/(x^2 + 1), x]", "Log[1 + x^2]");
    run_eq("Integrate[1/(x - 1), x]", "Log[-1 + x]");
    /* General correctness checks — formatting may vary but the
     * differential check is the contract. */
    assert_integral_correct("(2 x + 3)/(x^2 + 3 x + 5)");
    assert_integral_correct("1/(x - a)");
}

static void test_integrate_rational_via_recognition(void) {
    /* num is exactly c·D[pol] with pol a squarefull factor of den. */
    assert_integral_correct("1/(x - a)^2");
    assert_integral_correct("(2 x + 3)/(x^2 + 3 x + 5)^2");
    assert_integral_correct("(x - 1)/(x + 1)^3");
    /* Improper rational: poly part + recognition. */
    assert_integral_correct("(x^3 + 2 x)/(x^2 + 1)");
}

static void test_integrate_unevaluated(void) {
    /* No closed form via Phase 1 (LRT is Phase 2). */
    run_eq("Integrate[1/(x^2 + 1), x]", "Integrate[1/(1 + x^2), x]");
    /* Non-rational integrand: stays as Integrate[...]. */
    run_eq("Integrate[Sin[x], x]", "Integrate[Sin[x], x]");
}

/* ------------------------------------------------------------------ */
/* PolynomialQuotientRemainder & SubresultantPolynomialRemainders      */
/* ------------------------------------------------------------------ */

static void test_polynomial_quotient_remainder(void) {
    run_eq("PolynomialQuotientRemainder[x^3 + x + 1, x^2 + 1, x]",
           "{x, 1}");
    run_eq("PolynomialQuotientRemainder[x^4 - 1, x^2 - 1, x]",
           "{1 + x^2, 0}");
    /* With the Extension option. */
    run_eq("PolynomialQuotientRemainder[x^2 - 2, x - Sqrt[2], x, Extension -> Sqrt[2]]",
           "{Sqrt[2] + x, 0}");
}

static void test_subresultant_chain(void) {
    /* Coprime inputs: chain ends at a constant. */
    run_eq("SubresultantPolynomialRemainders[x^4 + 1, 2 x^3, x]",
           "{1 + x^4, 2 x^3, 2}");
    /* Standard univariate example used to drive the LRT pipeline. */
    Expr* e = parse_expression(
        "SubresultantPolynomialRemainders[x^3 - 1, x^2 - 1, x]");
    Expr* res = evaluate(e);
    /* The chain has at least 3 elements — first two are the inputs
     * (with bigger-degree first), third is the pseudo-remainder. */
    ASSERT(res->type == EXPR_FUNCTION);
    ASSERT(res->data.function.arg_count >= 3);
    expr_free(e);
    expr_free(res);
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

int main(void) {
    symtab_init();
    core_init();

    TEST(test_helpers_content);
    TEST(test_helpers_primitive);
    TEST(test_helpers_leading_coefficient);
    TEST(test_helpers_monic);

    TEST(test_integrate_polynomial);

    TEST(test_hermite_simple);
    TEST(test_hermite_reference);

    TEST(test_integrate_polynomial_dispatch);
    TEST(test_integrate_log_recognition);
    TEST(test_integrate_rational_via_recognition);
    TEST(test_integrate_unevaluated);

    TEST(test_polynomial_quotient_remainder);
    TEST(test_subresultant_chain);

    printf("All Phase 1 (IntegrateRational scaffolding + Hermite) tests passed!\n");
    return 0;
}
