/*
 * test_solvealways.c
 *
 * Unit tests for `SolveAlways` (src/solvealways.c).
 *
 * The result is the output of `Solve` on the collected coefficient
 * system; output ordering therefore follows whatever Solve produces
 * today.  When the order changes intentionally, update the strings.
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

/* Linear: a x + b == 0  ⇒  a = 0, b = 0 (works for every x).
 * Solve's rule-ordering puts the constant-coefficient equation first,
 * so the surfaced order is (b, a). */
static void test_linear(void) {
    run_test("SolveAlways[a x + b == 0, x]",
             "List[List[Rule[b, 0], Rule[a, 0]]]");
}

/* Vars as a one-element list -- same answer. */
static void test_vars_as_list(void) {
    run_test("SolveAlways[a x + b == 0, {x}]",
             "List[List[Rule[b, 0], Rule[a, 0]]]");
}

/* List of equations: every coefficient of x in every eq must vanish. */
static void test_list_of_equations(void) {
    run_test("SolveAlways[{a x + b == 0, c x + d == 0}, x]",
             "List[List[Rule[b, 0], Rule[a, 0], Rule[d, 0], Rule[c, 0]]]");
}

/* And of equations: identical to the list case. */
static void test_and_of_equations(void) {
    run_test("SolveAlways[a x + b == 0 && c x + d == 0, x]",
             "List[List[Rule[b, 0], Rule[a, 0], Rule[d, 0], Rule[c, 0]]]");
}

/* Bivariate: (a+b) x + (a-b) y == 0 ⇒ a = b = 0. */
static void test_bivariate(void) {
    run_test("SolveAlways[(a + b) x + (a - b) y == 0, {x, y}]",
             "List[List[Rule[a, 0], Rule[b, 0]]]");
}

/* No parameters present in eqns: result is {} per spec. */
static void test_no_parameters(void) {
    run_test("SolveAlways[2 x + 3 == 0, x]",
             "List[]");
}

/* Empty list of equations: no parameters present -> {}. */
static void test_empty_list_of_equations(void) {
    run_test("SolveAlways[{}, x]",
             "List[]");
}

/* Tautology x == x evaluates pre-builtin to True; no parameters -> {}.
 * (Equal[a, a] also reduces to True, so we cannot test a tautology
 * with a free parameter via that surface syntax -- the parameter is
 * stripped before SolveAlways sees it.) */
static void test_tautology(void) {
    run_test("SolveAlways[x == x, x]",
             "List[]");
}

/* One unknown relation: (a - b) x == 0 ⇒ a == b.  Solve picks `b` as
 * the dependent variable (Solve::svars warning is expected). */
static void test_one_unknown_relation(void) {
    run_test("SolveAlways[(a - b) x == 0, x]",
             "List[List[Rule[b, a]]]");
}

/* Rational (denominator-bearing) equations: the identity holds for all x iff the
 * NUMERATOR of the combined fraction vanishes, so SolveAlways must clear the
 * var-denominator (Together + Numerator) before coefficient extraction.  Feeding the
 * raw rational function to CoefficientList previously gave a spurious {} no-solution. */
static void test_rational_clear_denominator(void) {
    run_test("SolveAlways[a/(x - 1) + b/(x - 2) == (3 x - 5)/((x - 1)(x - 2)), x]",
             "List[List[Rule[a, 2], Rule[b, 1]]]");
    run_test("SolveAlways[A/(x - 1) + B/(x + 1) == 1/(x^2 - 1), x]",
             "List[List[Rule[A, Rational[1, 2]], Rule[B, Rational[-1, 2]]]]");
    /* A genuinely inconsistent rational identity is still {}. */
    run_test("SolveAlways[a/(x - 1) == 1/(x - 2), x]", "List[]");
}

int main(void) {
    symtab_init();
    core_init();

    printf("Running solvealways tests...\n");
    TEST(test_linear);
    TEST(test_vars_as_list);
    TEST(test_list_of_equations);
    TEST(test_and_of_equations);
    TEST(test_bivariate);
    TEST(test_no_parameters);
    TEST(test_empty_list_of_equations);
    TEST(test_tautology);
    TEST(test_one_unknown_relation);
    TEST(test_rational_clear_denominator);
    printf("All solvealways tests passed!\n");
    return 0;
}
