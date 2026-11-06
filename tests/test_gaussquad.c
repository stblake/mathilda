/* Tests for GaussianQuadratureWeights / GaussianQuadratureError.
 *
 * Expected values are independent of the implementation: closed-form
 * Gauss-Legendre nodes and weights for n = 2, 3, 5 (the zeros of P_2, P_3,
 * P_5 and w = 2/((1 - x^2) P_n'(x)^2) worked by hand), the defining exactness
 * of the rule on monomials of degree <= 2n - 1 against Integrate[x^k] =
 * 1/(k + 1) (and its failure at degree 2n), and the error-term constant
 * checked on f = x^(2n), where the error formula is exact.
 */

#include "core.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "test_utils.h"

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

/* Assert that `input` evaluates to True. */
static void assert_true(const char* input) {
    char* s = eval_str(input);
    ASSERT_MSG(strcmp(s, "True") == 0, "%s  ->  %s", input, s);
    free(s);
}

static void test_closed_form_n2_machine(void) {
    assert_true("g = GaussianQuadratureWeights[2, -1, 1]; "
                "Max[Abs[g - {{-1/Sqrt[3], 1}, {1/Sqrt[3], 1}}]] < 3*10^-16");
}

static void test_closed_form_n3_30digits(void) {
    assert_true("g = GaussianQuadratureWeights[3, -1, 1, 30]; "
                "Max[Abs[g - {{-Sqrt[3/5], 5/9}, {0, 8/9}, {Sqrt[3/5], 5/9}}]] < 10^-30");
    /* every entry carries the requested precision */
    assert_true("Min[Precision /@ Flatten[GaussianQuadratureWeights[3, -1, 1, 30]]] >= 30");
}

static void test_closed_form_n5_40digits(void) {
    /* P_5 zeros 0, +-(1/3) Sqrt[5 -+ 2 Sqrt[10/7]];
     * weights 128/225, (322 +- 13 Sqrt[70])/900. */
    assert_true("g = GaussianQuadratureWeights[5, -1, 1, 40]; "
                "e = {{-Sqrt[5 + 2 Sqrt[10/7]]/3, (322 - 13 Sqrt[70])/900}, "
                "     {-Sqrt[5 - 2 Sqrt[10/7]]/3, (322 + 13 Sqrt[70])/900}, "
                "     {0, 128/225}, "
                "     {Sqrt[5 - 2 Sqrt[10/7]]/3, (322 + 13 Sqrt[70])/900}, "
                "     {Sqrt[5 + 2 Sqrt[10/7]]/3, (322 - 13 Sqrt[70])/900}}; "
                "Max[Abs[g - e]] < 10^-40");
}

static void test_exactness_degree_2n_minus_1_50digits(void) {
    /* 12-point rule on [0, 1] at 50 digits: x^k integrated to 50 digits
     * for k <= 23, and visibly not for k = 24. */
    assert_true("r = GaussianQuadratureWeights[12, 0, 1, 50]; "
                "q[k_] := Total[#[[2]] #[[1]]^k & /@ r]; "
                "Max[Table[Abs[q[k] - 1/(k + 1)], {k, 0, 23}]] < 10^-49");
    assert_true("Abs[q[24] - 1/25] > 10^-20");
}

static void test_interval_mapping(void) {
    /* [a, b] = [-2, 3]: weights sum to 5, nodes increasing inside (a, b),
     * and Integrate[x^5 - x^2, {x, -2, 3}] = 665/6 - 35/3 is exact. */
    assert_true("r = GaussianQuadratureWeights[4, -2, 3, 30]; "
                "And[Abs[Total[r[[All, 2]]] - 5] < 10^-29, "
                "    OrderedQ[r[[All, 1]]], -2 < r[[1, 1]], r[[-1, 1]] < 3, "
                "    Abs[Total[#[[2]] (#[[1]]^5 - #[[1]]^2) & /@ r] - (665/6 - 35/3)] < 10^-26]");
    /* exact rational endpoints at machine precision */
    assert_true("r = GaussianQuadratureWeights[6, 1/3, 2]; "
                "Abs[Total[#[[2]] #[[1]]^11 & /@ r] - (2^12 - (1/3)^12)/12] < 10^-12");
}

static void test_odd_order_centre(void) {
    assert_true("GaussianQuadratureWeights[1, 0, 2] == {{1., 2.}}");
    /* n = 7: centre weight 2/P_7'(0)^2 = 2/(7 P_6(0))^2 = 2/(35/16)^2 = 512/1225 */
    assert_true("r = GaussianQuadratureWeights[7, -1, 1, 25]; r[[4, 1]] == 0 && "
                "Abs[r[[4, 2]] - 512/1225] < 10^-25");
}

static void test_large_n_machine(void) {
    /* n = 1000: symmetric, weights sum to 2, and every node and weight
     * agrees with the 30-digit rule (weights to 1e-12 relative). */
    assert_true("{t, g} = AbsoluteTiming[GaussianQuadratureWeights[1000, -1, 1]]; "
                "h = GaussianQuadratureWeights[1000, -1, 1, 30]; "
                "And[Abs[Total[g[[All, 2]]] - 2] < 10^-13, "
                "    g[[All, 1]] == -Reverse[g[[All, 1]]], "
                "    Max[Abs[g[[All, 1]] - h[[All, 1]]]] < 3*10^-16, "
                "    Max[Abs[(g[[All, 2]] - h[[All, 2]])/h[[All, 2]]]] < 10^-12]");
    /* milliseconds in practice; a generous bound for loaded CI machines */
    assert_true("t < 0.5");
}

static void test_error_term(void) {
    assert_true("GaussianQuadratureError[3, f, a, b] - (b - a)^7 Derivative[6][f]/2016000 // Simplify // # === 0 &");
    /* For f = x^6 the formula is exact: f^(6) = 720, so on [0, 1] the
     * difference Integrate - rule is 720/2016000 = 1/2800. */
    assert_true("r = GaussianQuadratureWeights[3, 0, 1, 40]; "
                "Abs[(1/7 - Total[#[[2]] #[[1]]^6 & /@ r]) - 1/2800] < 10^-38");
}

static void test_bad_arguments(void) {
    test_mute_messages();
    assert_true("Head[GaussianQuadratureWeights[0, 0, 1]] === GaussianQuadratureWeights");
    assert_true("Head[GaussianQuadratureWeights[3, a, 1]] === GaussianQuadratureWeights");
    assert_true("Head[GaussianQuadratureWeights[3, 0, 1, -2]] === GaussianQuadratureWeights");
    assert_true("MemberQ[Attributes[GaussianQuadratureWeights], Protected]");
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_closed_form_n2_machine);
    TEST(test_closed_form_n3_30digits);
    TEST(test_closed_form_n5_40digits);
    TEST(test_exactness_degree_2n_minus_1_50digits);
    TEST(test_interval_mapping);
    TEST(test_odd_order_centre);
    TEST(test_large_n_machine);
    TEST(test_error_term);
    TEST(test_bad_arguments);
    printf("All gaussquad tests passed.\n");
    return 0;
}
