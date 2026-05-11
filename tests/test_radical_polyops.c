/* Tests for the algebraic-generator (radical / exponential) pass shared
 * by Factor, Apart, Together, and Cancel.  When the input contains a
 * sub-expression u with a fractional rational exponent (radical case)
 * or with a rational-multiple-of-an-atom exponent (exponential case),
 * the four polynomial operations substitute Power[u, A/m] -> g
 * (a fresh symbol), run in g, then back-substitute.
 *
 * These tests pin the headline behaviours described in the user's
 * original request:
 *
 *   Factor[1 + r^(1/5) - 2 r^(2/5) - r^(3/5) - 2 r^(4/5)]
 *     = -(1 + 2 r^(1/5)) (-1 + r^(1/5) + r^(3/5))
 *
 *   Factor[<same with r replaced by Log[r]>]
 *     = -(1 + 2 Log[r]^(1/5)) (-1 + Log[r]^(1/5) + Log[r]^(3/5))
 *
 *   Apart[1 / (-1 + r^(3/7))]   gives partial fractions in r^(1/7)
 *
 *   Factor[Exp[2x] + 2 Exp[x] + 1] = (1 + E^x)^2
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

static void run_test(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* res_str = expr_to_string(res);
    if (strcmp(res_str, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, res_str);
        ASSERT_STR_EQ(res_str, expected);
    } else {
        printf("PASS: %s -> %s\n", input, res_str);
    }
    free(res_str);
    expr_free(e);
    expr_free(res);
}

static void test_factor_radical_symbol(void) {
    /* Polynomial in r^(1/5).  The substitution g = r^(1/5) maps
     * 1 + r^(1/5) - 2 r^(2/5) - r^(3/5) - 2 r^(4/5)
     *   -> 1 + g - 2 g^2 - g^3 - 2 g^4
     *   -> -(1+2g)(g^3+g-1) after Factor. */
    run_test("Factor[1 + r^(1/5) - 2 r^(2/5) - r^(3/5) - 2 r^(4/5)]",
             "-(-1 + r^(1/5) + r^(3/5)) (1 + 2 r^(1/5))");

    /* Same shape, simpler exponents. */
    run_test("Factor[r^(1/2) - 1]", "-1 + Sqrt[r]");
    run_test("Factor[r - 1]", "-1 + r");
}

static void test_factor_radical_compound_base(void) {
    /* Non-symbol base: the algebraic generator is a function call
     * (Log[r]) whose Power forms participate as Log[r]^(p/q). */
    run_test("Factor[1 + Log[r]^(1/5) - 2 Log[r]^(2/5) - Log[r]^(3/5) - 2 Log[r]^(4/5)]",
             "-(-1 + Log[r]^(1/5) + Log[r]^(3/5)) (1 + 2 Log[r]^(1/5))");
}

static void test_factor_exp(void) {
    /* Exponential case (atom != 1): treat E^x as the polynomial generator.
     * E^(2x) + 2 E^x + 1  ->  g^2 + 2g + 1 = (g+1)^2  ->  (1 + E^x)^2 */
    run_test("Factor[Exp[2x] + 2 Exp[x] + 1]", "(1 + E^x)^2");
    run_test("Factor[Exp[2x] - 1]", "(-1 + E^x) (1 + E^x)");
    run_test("Factor[Exp[2x] + 3 Exp[x] + 2]", "(1 + E^x) (2 + E^x)");
    run_test("Factor[E^(4x) - 1]", "(-1 + E^x) (1 + E^x) (1 + E^(2 x))");

    /* Mixed rational-coefficient exponents: atom = x, c values {1, 1/2}
     * give m = 2, g = E^(x/2), so E^x -> g^2 and E^(x/2) -> g. */
    run_test("Factor[Exp[x] - 2 Exp[x/2] + 1]", "(-1 + E^(1/2 x))^2");
}

static void test_apart_radical(void) {
    /* The headline regression: Apart should partial-fraction this in
     * r^(1/7) by substituting g = r^(1/7), then 1/(g^3-1) decomposes,
     * and back-substitution yields the radical form. */
    run_test("Apart[1/(-1 + r^(3/7))]",
             "1/3/(-1 + r^(1/7)) + (-2/3 - 1/3 r^(1/7))/(1 + r^(1/7) + r^(2/7))");
}

static void test_no_trigger_when_unhelpful(void) {
    /* Single Power[E, x] with no other matching site: the substitution
     * g = E^x would just rename the term, so the pass should NOT fire
     * and the input should come back unchanged. */
    run_test("Factor[Exp[x] + 1]", "1 + E^x");

    /* Pure-rational input with no fractional exponents: no radical pass. */
    run_test("Factor[r + r^2]", "r (1 + r)");
}

int main(void) {
    symtab_init();
    core_init();

    printf("Running radical polynomial-op tests...\n");
    TEST(test_factor_radical_symbol);
    TEST(test_factor_radical_compound_base);
    TEST(test_factor_exp);
    TEST(test_apart_radical);
    TEST(test_no_trigger_when_unhelpful);
    printf("All radical polynomial-op tests passed!\n");
    return 0;
}
