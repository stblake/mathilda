#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include <string.h>
#include <stdlib.h>

/* Shared driver: parse, evaluate, compare the printed form. */
static void check(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* res_str = expr_to_string(res);
    if (strcmp(res_str, expected) != 0) {
        printf("Cherry sigma test failed: %s\n  expected: %s\n  got:      %s\n",
               input, expected, res_str);
        ASSERT(0);
    }
    free(res_str);
    expr_free(e);
    expr_free(res);
}

/* ---- Integrate`SigmaDecomposition — Thm 4.4, all-equal degree-1 ---- */

/* Existence: Phi in K[P], P = Prod f_j. */
void test_sigma_exists() {
    /* Ex 5.1 reduced function x^2/2 = 1/2 + 1/2 (x^2-1). */
    check("Integrate`SigmaDecomposition[x^2/2, {x-1, x+1}, x]",
          "{{1/2, {0, 0}}, {1/2, {1, 1}}}");
    /* w^2 + w over Sigma = (x-1, x+1). */
    check("Integrate`SigmaDecomposition[(x^2-1)^2 + (x^2-1), {x-1, x+1}, x]",
          "{{1, {1, 1}}, {1, {2, 2}}}");
    /* A pure constant is a single 0-exponent term. */
    check("Integrate`SigmaDecomposition[1, {x-1, x+1}, x]", "{{1, {0, 0}}}");
    check("Integrate`SigmaDecomposition[7, {x-1, x+1}, x]", "{{7, {0, 0}}}");
    /* Three-factor product base P = x(x-1)(x+1) = x^3 - x. */
    check("Integrate`SigmaDecomposition[3 + 5 (x^3-x)^2, {x, x-1, x+1}, x]",
          "{{3, {0, 0, 0}}, {5, {2, 2, 2}}}");
}

/* Non-existence: Thm 4.4 termination proves no decomposition. */
void test_sigma_nonexistent() {
    /* Ex 5.2 reduced function x/2 (odd) has no even-power decomposition. */
    check("Integrate`SigmaDecomposition[x/2, {x-1, x+1}, x]", "$Failed");
    /* x^2/Log[x^3-x] reduced function x^2/(3x^2-1) — not a polynomial in P. */
    check("Integrate`SigmaDecomposition[x^2/(3x^2-1), {x, x-1, x+1}, x]", "$Failed");
    /* 1/(x^2+1): a proper rational, not in K[x^2-1]. */
    check("Integrate`SigmaDecomposition[1/(x^2+1), {x-1, x+1}, x]", "$Failed");
    /* A pole in a factor: negative multiplicity is out of T = Z_{>=0}. */
    check("Integrate`SigmaDecomposition[1/(x-1), {x-1, x+1}, x]", "$Failed");
}

/* ---- Integrate`LiElementaryQ — the li decision property ---- */

/* True: the Cherry engine exhibits an elementary + LogIntegral antiderivative. */
void test_lielementary_true() {
    check("Integrate`LiElementaryQ[x^3/Log[x^2-1], x]", "True");  /* Ex 5.1 */
    check("Integrate`LiElementaryQ[x/Log[x]^2, x]", "True");      /* d1 */
    check("Integrate`LiElementaryQ[x^2/Log[x+1], x]", "True");    /* d3 */
    check("Integrate`LiElementaryQ[1/(Log[x]+3), x]", "True");    /* d2 rescale */
    check("Integrate`LiElementaryQ[x/Log[x^2-1], x]", "True");    /* A/w' constant */
    check("Integrate`LiElementaryQ[1/Log[x], x]", "True");        /* = LogIntegral[x] */
}

/* False: the Sigma-decomposition non-existence certificate fires. */
void test_lielementary_false() {
    check("Integrate`LiElementaryQ[x^2/Log[x^2-1], x]", "False"); /* Ex 5.2 */
    check("Integrate`LiElementaryQ[x^2/Log[x^3-x], x]", "False"); /* multi-log */
    check("Integrate`LiElementaryQ[1/Log[x^2-1], x]", "False");   /* A/w'=1/(2x) not in K[w] */
}

/* Contrast with Risch`ElementaryIntegralQ: a li-integrable integrand is NOT
 * elementary, so the two predicates disagree on the li cases (this is the whole
 * point of a separate li decision). */
void test_lielementary_vs_elementary() {
    check("Integrate`LiElementaryQ[x^3/Log[x^2-1], x]", "True");
    check("Risch`ElementaryIntegralQ[x^3/Log[x^2-1], x]", "False");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_sigma_exists);
    TEST(test_sigma_nonexistent);
    TEST(test_lielementary_true);
    TEST(test_lielementary_false);
    TEST(test_lielementary_vs_elementary);

    printf("All Cherry Sigma-decomposition tests passed!\n");
    return 0;
}
