/*
 * test_solveradicals.c
 *
 * Unit tests for the radicals specialist `Solve`SolveRadicalsEquality`
 * (src/solverad.c).  Covers:
 *   - single-base substitution (x^(p/q) for a single base in x)
 *   - rational input that needs Together-then-numerator preprocessing
 *   - multi-base elimination via Resultant chain (Sqrt[x+1] + Sqrt[x-1])
 *   - extraneous-root filtering (numerical verification at the tail)
 *   - empty-solution detection (Sqrt[...] == -1 in Reals)
 *   - parametric inputs (Solve::nongen)
 *
 * Each test pins down the FullForm of the surviving solution list, so
 * any regression -- in the substitution, the resultant elimination, or
 * the verifier -- will fail it.
 */

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

static void run_test(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    if (!e) {
        printf("FAIL: failed to parse: %s\n", input);
        ASSERT(0);
        return;
    }
    Expr* res = evaluate(e);
    char* res_str = expr_to_string_fullform(res);
    if (strcmp(res_str, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n",
               input, expected, res_str);
        free(res_str);
        expr_free(res);
        expr_free(e);
        ASSERT(0);
        return;
    }
    printf("PASS: %s -> %s\n", input, res_str);
    free(res_str);
    expr_free(res);
    expr_free(e);
}

/* Sqrt[x] + 3 == 5  ->  x = 4.  Trivial: single Sqrt, no extraneous
 * roots.  Verifies the substitution u = Sqrt[x] -> Power[x, 1/2] is
 * detected and the resultant chain collapses to 4 - x. */
static void test_single_sqrt_constant_rhs(void) {
    run_test("Solve[Sqrt[x] + 3 == 5, x]",
             "List[List[Rule[x, 4]]]");
}

/* Sqrt[x] + 3 x == 5  ->  3 u^2 + u - 5 == 0 after u = Sqrt[x].
 * The negative-u branch back-substitutes to a spurious x; the verifier
 * must drop it.  Survivor: (31 - Sqrt[61]) / 18. */
static void test_sqrt_plus_linear_keeps_principal_branch(void) {
    run_test("Solve[Sqrt[x] + 3 x == 5, x]",
             "List[List[Rule[x, Times[Rational[1, 18], "
                  "Plus[31, Times[-1, Power[61, Rational[1, 2]]]]]]]]");
}

/* x - 8 Sqrt[x] + 15 == 0  ->  both u = 3 and u = 5 survive; x = 9
 * and x = 25 are both valid principal-branch roots. */
static void test_quadratic_in_sqrt(void) {
    run_test("Solve[x - 8 Sqrt[x] + 15 == 0, x]",
             "List[List[Rule[x, 9]], List[Rule[x, 25]]]");
}

/* Mixed exponents 1/2 and 1/4 of the same base x.  L = lcm(2, 4) = 4,
 * so u = x^(1/4).  Equation becomes u^2 + 3 u - 5 == 0.  One branch
 * survives verification: x = (311 - 57 Sqrt[29]) / 2. */
static void test_mixed_exponent_same_base(void) {
    run_test("Solve[Sqrt[x] + 3 x^(1/4) == 5, x]",
             "List[List[Rule[x, Times[Rational[1, 2], "
                  "Plus[311, Times[-57, Power[29, Rational[1, 2]]]]]]]]");
}

/* x^(2/3) - 3 x^(1/3) + 2 == 0.  Negative-q in the rational exponent;
 * u = x^(1/3) gives u^2 - 3 u + 2 == 0 -> u = 1 or 2 -> x = 1 or 8. */
static void test_cube_root_exponents(void) {
    run_test("Solve[x^(2/3) - 3 x^(1/3) + 2 == 0, x]",
             "List[List[Rule[x, 1]], List[Rule[x, 8]]]");
}

/* (x + Sqrt[x]) / Sqrt[x] + Sqrt[x] / (x + Sqrt[x]) == 4.  Rational
 * radical: Together-then-Numerator turns this into a polynomial in
 * x and u = Sqrt[x].  The resultant carries an extra x^2 factor;
 * verification drops x = 0 (denominator zero) and the conjugate
 * 4 - 2 Sqrt[3] root (residual ~2.31, not 0). */
static void test_rational_radical_extraneous_filter(void) {
    run_test("Solve[(x + Sqrt[x])/Sqrt[x] + Sqrt[x]/(x + Sqrt[x]) == 4, x]",
             "List[List[Rule[x, Times[Rational[1, 2], "
                  "Plus[8, Times[4, Power[3, Rational[1, 2]]]]]]]]");
}

/* x + Sqrt[x - 1] == 1.  Polynomial elimination yields x = 1 and x = 2;
 * only x = 1 survives back-substitution.  Confirms that the
 * verification step really runs against the *original* equation. */
static void test_extraneous_root_filter(void) {
    run_test("Solve[x + Sqrt[x - 1] == 1, x]",
             "List[List[Rule[x, 1]]]");
}

/* Two distinct radical bases.  The elimination chain processes one
 * Resultant per generator: first u_1 = Sqrt[x+1], then u_2 = Sqrt[x-1]. */
static void test_two_distinct_bases(void) {
    run_test("Solve[Sqrt[x + 1] + Sqrt[x - 1] == 3, x]",
             "List[List[Rule[x, Rational[85, 36]]]]");
}

/* Sqrt[x + 5] + Sqrt[x] == -1.  Both radicals are real-non-negative
 * on the real line, so the sum cannot equal -1.  Polynomial
 * elimination yields a candidate x = 4 (where Sqrt[9] + Sqrt[4] = 5);
 * verification drops it -- residual is +6, not 0.  Expected output:
 * the empty list. */
static void test_no_real_solution(void) {
    run_test("Solve[Sqrt[x + 5] + Sqrt[x] == -1, x]",
             "List[]");
}

/* Parametric input.  Solve cannot decide the sign of the residual
 * symbolically, so both branches survive and Solve::nongen is raised. */
static void test_parametric_nongen(void) {
    run_test("Solve[Sqrt[a x + c] + 3 x == 5, x]",
             "List["
             "List[Rule[x, Times[Rational[1, 18], "
                  "Plus[30, a, "
                       "Power[Plus[Power[Plus[-30, Times[-1, a]], 2], "
                                  "Times[-36, Plus[25, Times[-1, c]]]], "
                             "Rational[1, 2]]]]]], "
             "List[Rule[x, Times[Rational[1, 18], "
                  "Plus[30, a, "
                       "Times[-1, Power[Plus[Power[Plus[-30, Times[-1, a]], 2], "
                                              "Times[-36, Plus[25, Times[-1, c]]]], "
                                        "Rational[1, 2]]]]]]]]");
}

/* The qualified builtin is directly callable: a quick smoke test
 * verifies the dispatch wiring registered the right entry. */
static void test_qualified_builtin(void) {
    run_test("Solve`SolveRadicalsEquality[Sqrt[x] + 3 == 5, x]",
             "List[List[Rule[x, 4]]]");
}

int main(void) {
    symtab_init();
    core_init();

    printf("Running solveradicals tests...\n");
    TEST(test_single_sqrt_constant_rhs);
    TEST(test_sqrt_plus_linear_keeps_principal_branch);
    TEST(test_quadratic_in_sqrt);
    TEST(test_mixed_exponent_same_base);
    TEST(test_cube_root_exponents);
    TEST(test_rational_radical_extraneous_filter);
    TEST(test_extraneous_root_filter);
    TEST(test_two_distinct_bases);
    TEST(test_no_real_solution);
    TEST(test_parametric_nongen);
    TEST(test_qualified_builtin);

    printf("All solveradicals tests passed!\n");
    return 0;
}
