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
          "List[Rule[x, Times[Rational[1, 2], Plus[3, Times[-1, Power[5, Rational[1, 2]]]]]], "
               "Rule[y, Times[Rational[1, 2], Plus[3, Power[5, Rational[1, 2]]]]]], "
          "List[Rule[x, Times[Rational[1, 2], Plus[3, Power[5, Rational[1, 2]]]]], "
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
    /* x = y = -+ I/Sqrt[2].  The eliminant returns x in the canonical
     * unrationalised form Complex[0, -+1] 2^(-1/2) and y in the equal
     * rationalised form Complex[0, -+1/2] 2^(1/2) (= -+ I Sqrt[2]/2); both are
     * the same value.  (A single global Times radical-canonicalisation to force
     * one convention breaks the unrationalised Normalize/Series forms, so the
     * two solution components are left in their equal natural forms.) */
    run_test("Solve[x^2 + y^2 == -1 && x == y, {x, y}]",
        "List["
          "List[Rule[x, Times[Complex[0, -1], Power[2, Rational[-1, 2]]]], "
               "Rule[y, Times[Complex[0, Rational[-1, 2]], Power[2, Rational[1, 2]]]]], "
          "List[Rule[x, Times[Complex[0, 1], Power[2, Rational[-1, 2]]]], "
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
             "List[List[Rule[x, 2], Rule[y, 3]], "
                  "List[Rule[x, 3], Rule[y, 2]]]");
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
                  "List[Rule[x, -4], Rule[y, 3]], "
                  "List[Rule[x, 4], Rule[y, -3]], "
                  "List[Rule[x, 4], Rule[y, 3]]]");
}

/* Three variables (elementary symmetric): permutations of {1,2,3}. */
static void test_three_var_symmetric(void) {
    run_test("Solve[x + y + z == 6 && x y z == 6 && x^2 + y^2 + z^2 == 14, {x, y, z}]",
             "List[List[Rule[x, 1], Rule[y, 2], Rule[z, 3]], List[Rule[x, 1], Rule[y, 3], Rule[z, 2]], List[Rule[x, 2], Rule[y, 1], Rule[z, 3]], List[Rule[x, 2], Rule[y, 3], Rule[z, 1]], List[Rule[x, 3], Rule[y, 1], Rule[z, 2]], List[Rule[x, 3], Rule[y, 2], Rule[z, 1]]]");
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

/* Positive-dimensional ideal from a multi-equation SYSTEM (infinitely
 * many solutions): left unevaluated (head stays Solve), Solve::nsdim
 * emitted.  A conjunction never takes the single-equation-multivariable
 * fallback (which handles a lone equation like x^2 - y^2 == 0, solving
 * it for one variable as {{x -> -Sqrt[y^2]}, {x -> Sqrt[y^2]}}). */
static void test_positive_dimensional(void) {
    run_test("Solve[x^2 - y^2 == 0 && x - y == 0, {x, y}]",
             "Solve[And[Equal[Plus[Power[x, 2], Times[-1, Power[y, 2]]], 0], "
                       "Equal[Plus[x, Times[-1, y]], 0]], "
                   "List[x, y]]");
}

/* Single equation in several variables: solved for the earliest-listed
 * variable it is polynomial in, treating the rest as parameters -- an
 * explicit rule (Solve's job), not a Reduce-style case split. */
static void test_single_eq_multivar(void) {
    /* x y == 1  ->  x -> 1/y. */
    run_test("Solve[x y == 1, {x, y}]",
             "List[List[Rule[x, Power[y, -1]]]]");
    /* Sin[x] + y == 0: x is not polynomial, so solve for the next
     * variable y  ->  y -> -Sin[x]. */
    run_test("Solve[Sin[x] + y == 0, {x, y}]",
             "List[List[Rule[y, Times[-1, Sin[x]]]]]");
}

/* Non-polynomial system (transcendental head): left unevaluated. */
static void test_non_polynomial(void) {
    run_test("Solve[Sin[x] + y == 0 && x + y == 1, {x, y}]",
             "Solve[And[Equal[Plus[y, Sin[x]], 0], Equal[Plus[x, y], 1]], List[x, y]]");
}

/* Regression: a purely linear system must keep taking the linear path
 * and produce the unique solution unchanged. */
/* Unexpanded / product forms: solvenlsys must Expand each equation before
 * GBPoly conversion (gb_from_expr's single-term parser cannot ingest
 * Power[Plus,k] / Times[Plus,...]).  Two circles meeting on their radical
 * axis: x = 1/2, y = +-Sqrt[3]/2. */
static void test_unexpanded_two_circles(void) {
    run_test("Solve[{x^2 + y^2 == 1, (x - 1)^2 + y^2 == 1}, {x, y}]",
        "List["
          "List[Rule[x, Rational[1, 2]], "
               "Rule[y, Times[Rational[-1, 2], Power[3, Rational[1, 2]]]]], "
          "List[Rule[x, Rational[1, 2]], "
               "Rule[y, Times[Rational[1, 2], Power[3, Rational[1, 2]]]]]]");
}

/* Product of sums == 0 combined with a linear constraint: the x+y=0 branch
 * contradicts x+y=2, leaving the single solution {1,1}. */
static void test_unexpanded_product(void) {
    run_test("Solve[{(x + y) (x - y) == 0, x + y == 2}, {x, y}]",
        "List[List[Rule[x, 1], Rule[y, 1]]]");
}

/* Parametric system: free symbols (here `a`) that are not solve variables are
 * coefficients in Q(params).  Routed through the RationalFunctions-coefficient
 * Gröbner engine and solved symbolically: x^2 == a, x + y == 0 -> x = +-Sqrt[a],
 * y = -/+Sqrt[a]. */
static void test_parametric_system(void) {
    run_test("Solve[{x^2 == a, x + y == 0}, {x, y}]",
        "List[List[Rule[x, Times[-1, Power[a, Rational[1, 2]]]], Rule[y, Power[a, Rational[1, 2]]]], List[Rule[x, Power[a, Rational[1, 2]]], Rule[y, Times[-1, Power[a, Rational[1, 2]]]]]]");
}

static void test_linear_regression(void) {
    run_test("Solve[x + y == 3 && x - y == 1, {x, y}]",
             "List[List[Rule[x, 2], Rule[y, 1]]]");
    run_test("Solve[{2 x + y == 5, x - y == 1}, {x, y}]",
             "List[List[Rule[x, 2], Rule[y, 1]]]");
}

/* A quadratic over Q(params) is factored before the formula.  On
 * a^2 b^3 (1 + 3c + 2c^2) the formula gave (-3a^2b^3 -+ Sqrt[a^4 b^6])/(4a^2b^3),
 * each -1 or -1/2 only by the branch of the square root, and the
 * triangular search then took n2 = 0 on both branches: the norm system
 * below (from ParallelMixedTower's NormSearchAll) came back {}. */
static void test_parametric_square_discriminant(void) {
    run_test("Solve[a^2 b^3 (1 + 3 c + 2 c^2) == 0, c]",
             "List[List[Rule[c, -1]], List[Rule[c, Rational[-1, 2]]]]");
    run_test("Solve[a + (2 a + b) c + 2 b c^2 == 0, c]",
             "List[List[Rule[c, Rational[-1, 2]]], List[Rule[c, Times[-1, a, Power[b, -1]]]]]");
    run_test("Solve[{n0^2 == b^2 (1 + c), 2 n0 n1 == 0, n1^2 + 2 n0 n2 == 2 a b c, "
             "2 n1 n2 == 0, n2^2 == a^2 (1 + c)}, {n0, n1, n2, c}]",
             "List[List[Rule[n0, 0], Rule[n1, Times[Rational[-1, 2], Power[Times[-8, a, b], Rational[1, 2]]]], Rule[n2, 0], Rule[c, -1]], "
             "List[Rule[n0, 0], Rule[n1, Times[Rational[1, 2], Power[Times[-8, a, b], Rational[1, 2]]]], Rule[n2, 0], Rule[c, -1]], "
             "List[Rule[n0, Times[Rational[-1, 4], Power[a, -1], Power[Times[8, Power[a, 2]], Rational[1, 2]], b]], Rule[n1, 0], "
             "Rule[n2, Times[Rational[1, 4], Power[Times[8, Power[a, 2]], Rational[1, 2]]]], Rule[c, Rational[-1, 2]]], "
             "List[Rule[n0, Times[Rational[1, 4], Power[a, -1], Power[Times[8, Power[a, 2]], Rational[1, 2]], b]], Rule[n1, 0], "
             "Rule[n2, Times[Rational[-1, 4], Power[Times[8, Power[a, 2]], Rational[1, 2]]]], Rule[c, Rational[-1, 2]]]]");
    /* irreducible over Q(a, b, c): the formula, unchanged */
    run_test("Solve[a x^2 + b x + c == 0, x]",
             "List[List[Rule[x, Times[Rational[1, 2], Power[a, -1], Plus[Times[-1, b], Times[-1, Power[Plus[Power[b, 2], Times[-4, a, c]], Rational[1, 2]]]]]]], "
             "List[Rule[x, Times[Rational[1, 2], Power[a, -1], Plus[Times[-1, b], Power[Plus[Power[b, 2], Times[-4, a, c]], Rational[1, 2]]]]]]]");
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
    test_single_eq_multivar();
    test_non_polynomial();
    test_unexpanded_two_circles();
    test_unexpanded_product();
    test_parametric_square_discriminant();
    test_parametric_system();
    test_linear_regression();

    printf("\nAll solvenlsys tests passed.\n");
    symtab_clear();
    return 0;
}
