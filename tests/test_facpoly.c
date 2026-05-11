#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

void run_test(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    if (!e) {
        printf("Failed to parse: %s\n", input);
        ASSERT(0);
    }
    Expr* res = evaluate(e);
    char* res_str = expr_to_string_fullform(res);
    if (strcmp(res_str, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, res_str);
        ASSERT(0);
    } else {
        printf("PASS: %s -> %s\n", input, res_str);
    }
    free(res_str);
    expr_free(res);
    expr_free(e);
}

void test_factorsquarefree() {
    run_test("FactorSquareFree[x^5-x^3-x^2+1]", "Times[Power[Plus[-1, x], 2], Plus[1, Times[2, x], Times[2, Power[x, 2]], Power[x, 3]]]");
    run_test("FactorSquareFree[x^4-9x^3+29x^2-39x+18]", "Times[Power[Plus[-3, x], 2], Plus[2, Times[-3, x], Power[x, 2]]]");
    run_test("FactorSquareFree[x^5-x^3*y^2-x^2*y^3+y^5]", "Times[Power[Plus[Times[-1, x], y], 2], Plus[Power[x, 3], Times[2, Times[Power[x, 2], y]], Times[2, Times[x, Power[y, 2]]], Power[y, 3]]]");
    run_test("FactorSquareFree[{(x^2-1)(x-1),(x^4-1)(x^2-1)}]", "List[Times[Plus[1, x], Power[Plus[-1, x], 2]], Times[Plus[1, Power[x, 2]], Power[Plus[-1, Power[x, 2]], 2]]]");
    run_test("FactorSquareFree[Exp[2 x] - Log[x] - 1/x]", "Plus[Times[-1, Power[x, -1]], Times[-1, Log[x]], Power[E, Times[2, x]]]");
    run_test("FactorSquareFree[1 - x^3]", "Plus[1, Times[-1, Power[x, 3]]]");
    run_test("FactorSquareFree[(1 - x^3)^2]", "Power[Plus[1, Times[-1, Power[x, 3]]], 2]");
}

/* F4 Phase 1: cheap squarefree pre-check.  These cases exercise the
 * fast path (multivariate squarefree inputs) and ensure the slow path
 * still kicks in correctly when a real repeated factor is present.
 * Expected outputs match the pre-F4 baseline -- the optimisation is
 * meant to be observationally invisible. */
void test_factorsquarefree_f4_fastpath() {
    /* Multivariate squarefree -- the recursive content extraction
     * still pulls out the (xy-1) factor; the pre-check just skips
     * the expensive multivariate gcd inside Yun's algorithm. */
    run_test("FactorSquareFree[Expand[(x*y - 1)*(x + y + z)*(x*y*z + 1)]]",
             "Times[Plus[-1, Times[x, y]], Plus[x, y, z, Times[Power[x, 2], y, z], Times[x, Power[y, 2], z], Times[x, y, Power[z, 2]]]]");

    /* Multivariate non-squarefree -- slow path must catch (x*y-1)^2. */
    run_test("FactorSquareFree[Expand[(x*y - 1)^2*(x + y + z)]]",
             "Times[Power[Plus[-1, Times[x, y]], 2], Plus[x, y, z]]");

    /* Multivariate non-squarefree with two squared factors. */
    run_test("FactorSquareFree[Expand[(x + y)^2*(x - y)^2*(x*y - 1)]]",
             "Times[Plus[-1, Times[x, y]], Power[Plus[Times[-1, Power[x, 2]], Power[y, 2]], 2]]");

    /* Three-variable squarefree input -- this is the kind of case
     * the F4 fast path is targeting. */
    run_test("FactorSquareFree[Expand[(1 - x^3)*(1 + y^4)*(1 + z^5)]]",
             "Times[Plus[-1, Power[x, 3]], Plus[1, Power[y, 4]], Plus[-1, Times[-1, Power[z, 5]]]]");

    /* Repeated factor that is constant in the main variable (the
     * recursive content path catches it). */
    run_test("FactorSquareFree[Expand[(y + 1)^2*(x + y)]]",
             "Times[Power[Plus[1, y], 2], Plus[x, y]]");
}

/* Regression: x^p - c with p equal to the prime first tried by
 * factor_zassenhaus (p = 13) used to hang in the squarefree probe.
 * The derivative p*x^(p-1) is identically zero mod p, leaving an
 * untrimmed all-zero divisor that drove upoly_div_rem_mod into an
 * infinite loop.  These cases must complete; for x^13 - 1 the
 * factorisation is non-trivial; the others are irreducible. */
void test_factor_xp_minus_c_regression() {
    run_test("Factor[x^13 - 1]",
             "Times[Plus[-1, x], Plus[1, x, Power[x, 2], Power[x, 3], Power[x, 4], Power[x, 5], Power[x, 6], Power[x, 7], Power[x, 8], Power[x, 9], Power[x, 10], Power[x, 11], Power[x, 12]]]");
    run_test("Factor[x^13 - 2]",
             "Plus[-2, Power[x, 13]]");
    run_test("Factor[x^26 - 2]",
             "Plus[-2, Power[x, 26]]");
    /* End-to-end user case: the irreducibility short-circuit's
     * univariate images (y^13 - c, z^14 - c) used to trigger the
     * same hang.  Full factorisation completes in single-digit
     * hundreds of milliseconds. */
    run_test("Factor[Expand[x^2 (1 - x^12) (1 + x - y^13) (1 - y - z^14)]]",
             "Times[Power[x, 2], Plus[-1, x], Plus[1, x], Plus[1, Power[x, 2]], Plus[1, x, Power[x, 2]], Plus[1, Times[-1, x], Power[x, 2]], Plus[1, Times[-1, Power[x, 2]], Power[x, 4]], Plus[-1, Times[-1, x], Power[y, 13]], Plus[1, Times[-1, y], Times[-1, Power[z, 14]]]]");
}

void test_factor() {
    run_test("Factor[1 + 2 x + x^2]", "Power[Plus[1, x], 2]");
    run_test("Factor[x^10 - 1]", "Times[Plus[-1, x], Plus[1, x], Plus[1, x, Power[x, 2], Power[x, 3], Power[x, 4]], Plus[1, Times[-1, x], Power[x, 2], Times[-1, Power[x, 3]], Power[x, 4]]]");
    run_test("Factor[1 - x^3]", "Times[-1, Plus[-1, x], Plus[1, x, Power[x, 2]]]");
    run_test("Factor[x^10 - y^10]", "Times[Plus[x, y], Plus[x, Times[-1, y]], Plus[Power[x, 4], Times[-1, Times[Power[x, 3], y]], Times[Power[x, 2], Power[y, 2]], Times[-1, Times[x, Power[y, 3]]], Power[y, 4]], Plus[Power[x, 4], Times[Power[x, 3], y], Times[Power[x, 2], Power[y, 2]], Times[x, Power[y, 3]], Power[y, 4]]]");
    run_test("Factor[x^3 - 6x^2 + 11x - 6]", "Times[Plus[-3, x], Plus[-2, x], Plus[-1, x]]");
    run_test("Factor[2x^3 y - 2a^2 x y - 3a^2 x^2 + 3a^4]", "Times[Plus[a, x], Plus[Times[-1, a], x], Plus[Times[-3, Power[a, 2]], Times[2, Times[x, y]]]]");
    run_test("Factor[(x^3+2x^2)/(x^2-4y^2)-(x+2)/(x^2-4y^2)]", "Times[Plus[-1, x], Plus[1, x], Plus[2, x], Power[Plus[x, Times[-2, y]], -1], Power[Plus[x, Times[2, y]], -1]]");
    run_test("Factor[{x^2-1, x^4-1, x^8-1}]", "List[Times[Plus[-1, x], Plus[1, x]], Times[Plus[-1, x], Plus[1, x], Plus[1, Power[x, 2]]], Times[Plus[-1, x], Plus[1, x], Plus[1, Power[x, 2]], Plus[1, Power[x, 4]]]]");
    run_test("Factor[1 < 1 + 2 x + x^2 + 1/(1+x) < 2]", "Less[Less[1, Times[Plus[2, x], Power[Plus[1, x], -1], Plus[1, x, Power[x, 2]]]], 2]");
}

int main() {
    symtab_init();
    core_init();

    printf("Running facpoly tests...\n");
    TEST(test_factorsquarefree);
    TEST(test_factorsquarefree_f4_fastpath);
    TEST(test_factor_xp_minus_c_regression);
    TEST(test_factor);
    printf("All facpoly tests passed!\n");
    return 0;
}
