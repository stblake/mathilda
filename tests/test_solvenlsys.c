/*
 * test_solvenlsys.c
 *
 * Unit tests for the nonlinear polynomial-system specialist
 * `Solve`SolveNonlinearSystem` (src/solvenlsys.c), reached through the
 * public `Solve` router for non-affine systems of polynomial equations.
 *
 * Outputs are compared against FullForm strings, so the canonical
 * solution representation (List[List[Rule[var, val], ...], ...]) is
 * asserted exactly.  Solution order matches the lex-Gröbner triangular
 * back-substitution; when the order changes intentionally, update the
 * strings.
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

/* Symmetric 2x2: x y == 1, x + y == 3  ->  golden-ratio conjugate pair. */
static void test_symmetric_pair(void) {
    run_test("Solve[x y == 1 && x + y == 3, {x, y}]",
        "List["
          "List[Rule[x, Times[Rational[-1, 2], Plus[-3, Power[5, Rational[1, 2]]]]], "
               "Rule[y, Times[Rational[1, 2], Plus[3, Power[5, Rational[1, 2]]]]]], "
          "List[Rule[x, Times[Rational[-1, 2], Plus[-3, Times[-1, Power[5, Rational[1, 2]]]]]], "
               "Rule[y, Times[Rational[1, 2], Plus[3, Times[-1, Power[5, Rational[1, 2]]]]]]]]");
}

/* Circle meets the diagonal:  x^2 + y^2 == 1, x == y  ->  +-1/Sqrt[2]. */
static void test_circle_diagonal(void) {
    run_test("Solve[x^2 + y^2 == 1 && x == y, {x, y}]",
        "List["
          "List[Rule[x, Times[-1, Power[2, Rational[-1, 2]]]], "
               "Rule[y, Times[-1, Power[2, Rational[-1, 2]]]]], "
          "List[Rule[x, Power[2, Rational[-1, 2]]], "
               "Rule[y, Power[2, Rational[-1, 2]]]]]");
}

/* Inconsistent system  ->  {}. */
static void test_inconsistent(void) {
    run_test("Solve[x^2 + y^2 == 1 && x^2 + y^2 == 2, {x, y}]",
             "List[]");
}

/* Complex solutions over the default domain (Complexes). */
static void test_complex_default(void) {
    run_test("Solve[x^2 + y^2 == -1 && x == y, {x, y}]",
        "List["
          "List[Rule[x, Times[Complex[0, Rational[-1, 2]], Power[2, Rational[1, 2]]]], "
               "Rule[y, Times[Complex[0, Rational[-1, 2]], Power[2, Rational[1, 2]]]]], "
          "List[Rule[x, Times[Complex[0, Rational[1, 2]], Power[2, Rational[1, 2]]]], "
               "Rule[y, Times[Complex[0, Rational[1, 2]], Power[2, Rational[1, 2]]]]]]");
}

/* Reals domain prunes the complex branch to {}. */
static void test_reals_drop(void) {
    run_test("Solve[x^2 + y^2 == -1 && x == y, {x, y}, Reals]",
             "List[]");
}

/* Integers domain:  x y == 6, x + y == 5  ->  {3,2} and {2,3}. */
static void test_integers(void) {
    run_test("Solve[x y == 6 && x + y == 5, {x, y}, Integers]",
             "List[List[Rule[x, 3], Rule[y, 2]], "
                  "List[Rule[x, 2], Rule[y, 3]]]");
}

/* Mixed: quadratic in x, linear in y. */
static void test_quadratic_and_linear(void) {
    run_test("Solve[x^2 == 4 && y == x + 1, {x, y}]",
             "List[List[Rule[x, -2], Rule[y, -1]], "
                  "List[Rule[x, 2], Rule[y, 3]]]");
}

/* Circle and hyperbola: four real intersection points. */
static void test_circle_hyperbola(void) {
    run_test("Solve[x^2 + y^2 == 25 && x^2 - y^2 == 7, {x, y}]",
             "List[List[Rule[x, -4], Rule[y, -3]], "
                  "List[Rule[x, 4], Rule[y, -3]], "
                  "List[Rule[x, -4], Rule[y, 3]], "
                  "List[Rule[x, 4], Rule[y, 3]]]");
}

/* Three variables (elementary symmetric): permutations of {1,2,3}. */
static void test_three_var_symmetric(void) {
    run_test("Solve[x + y + z == 6 && x y z == 6 && x^2 + y^2 + z^2 == 14, {x, y, z}]",
             "List[List[Rule[x, 3], Rule[y, 2], Rule[z, 1]], "
                  "List[Rule[x, 2], Rule[y, 3], Rule[z, 1]], "
                  "List[Rule[x, 3], Rule[y, 1], Rule[z, 2]], "
                  "List[Rule[x, 1], Rule[y, 3], Rule[z, 2]], "
                  "List[Rule[x, 2], Rule[y, 1], Rule[z, 3]], "
                  "List[Rule[x, 1], Rule[y, 2], Rule[z, 3]]]");
}

/* List-form (not And) of equations is accepted too. */
static void test_list_form(void) {
    run_test("Solve[{x^2 + y^2 == 1, x - y == 0}, {x, y}]",
        "List["
          "List[Rule[x, Times[-1, Power[2, Rational[-1, 2]]]], "
               "Rule[y, Times[-1, Power[2, Rational[-1, 2]]]]], "
          "List[Rule[x, Power[2, Rational[-1, 2]]], "
               "Rule[y, Power[2, Rational[-1, 2]]]]]");
}

/* Positive-dimensional ideal (infinitely many solutions): left
 * unevaluated (head stays Solve), and Solve::nsdim is emitted. */
static void test_positive_dimensional(void) {
    run_test("Solve[x^2 - y^2 == 0, {x, y}]",
             "Solve[Equal[Plus[Power[x, 2], Times[-1, Power[y, 2]]], 0], "
                   "List[x, y]]");
}

/* Non-polynomial system (transcendental head): left unevaluated. */
static void test_non_polynomial(void) {
    run_test("Solve[Sin[x] + y == 0 && x + y == 1, {x, y}]",
             "Solve[And[Equal[Plus[Sin[x], y], 0], Equal[Plus[x, y], 1]], "
                   "List[x, y]]");
}

/* Regression: a purely linear system must keep taking the linear path
 * and produce the unique solution unchanged. */
static void test_linear_regression(void) {
    run_test("Solve[x + y == 3 && x - y == 1, {x, y}]",
             "List[List[Rule[x, 2], Rule[y, 1]]]");
    run_test("Solve[{2 x + y == 5, x - y == 1}, {x, y}]",
             "List[List[Rule[x, 2], Rule[y, 1]]]");
}

int main(void) {
    symtab_init();
    core_init();

    test_symmetric_pair();
    test_circle_diagonal();
    test_inconsistent();
    test_complex_default();
    test_reals_drop();
    test_integers();
    test_quadratic_and_linear();
    test_circle_hyperbola();
    test_three_var_symmetric();
    test_list_form();
    test_positive_dimensional();
    test_non_polynomial();
    test_linear_regression();

    printf("\nAll solvenlsys tests passed.\n");
    symtab_clear();
    return 0;
}
