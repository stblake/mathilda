/* Unit tests for Minimize / Maximize (src/numerical_calculus/minimize.c).
 *
 * Minimize is the EXACT (symbolic) global optimizer. Exact input yields exact
 * output; inexact input delegates to NMinimize. The engine is sound by
 * construction: when a case cannot be decided exactly it leaves the expression
 * unevaluated rather than guessing, so several tests assert a *decline*
 * (Head[...] === Minimize) as the correct behaviour.
 *
 * Coverage:
 *   - Result shape / typing; {f_opt, {x -> x_opt, ...}}.
 *   - Univariate polynomial: quadratic, quartic, constant, high-degree, the
 *     degree-6 symmetric showcase, unbounded (odd degree / negative leading).
 *   - Maximize via objective negation + min/max duality.
 *   - Multivariate unconstrained polynomial with isolated minima.
 *   - Constrained polynomial over the Reals (disk, equality+inequality, LP).
 *   - Infeasible region -> {Infinity, {x -> Indeterminate, ...}}.
 *   - Inexact input -> numeric NMinimize fallback.
 *   - Sound declines: transcendental exact input, positive-dimensional minima,
 *     unbounded-constrained, bad arity.
 *   - Memory-hygiene smoke loop.
 *
 * Run binary directly: ./minimize_tests */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

/* Parse, evaluate, FullForm-compare. */
static void check_eq(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* got = expr_to_string_fullform(res);
    if (strcmp(got, expected) != 0) {
        fprintf(stdout, "FAIL: %s\n  expected: %s\n  got:      %s\n",
                input, expected, got);
        ASSERT_STR_EQ(got, expected);
    }
    free(got);
    expr_free(e);
    expr_free(res);
}

/* Parse+evaluate a boolean predicate and require it to be True. */
static void check_true(const char* input) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* got = expr_to_string_fullform(res);
    if (strcmp(got, "True") != 0) {
        fprintf(stdout, "FAIL: %s\n  expected: True\n  got:      %s\n",
                input, got);
        ASSERT_STR_EQ(got, "True");
    }
    free(got);
    expr_free(e);
    expr_free(res);
}

/* ------------------------------------------------------------------ */
/* 1. Result shape                                                     */
/* ------------------------------------------------------------------ */

static void test_result_shape(void) {
    check_eq("Head[Minimize[x^2, x]]", "List");
    check_eq("Length[Minimize[x^2, x]]", "2");
    check_eq("Head[Last[Minimize[x^2, x]]]", "List");
    check_eq("Head[First[Last[Minimize[x^2, x]]]]", "Rule");
}

/* ------------------------------------------------------------------ */
/* 2. Univariate unconstrained polynomial                              */
/* ------------------------------------------------------------------ */

static void test_quadratic(void) {
    check_eq("Minimize[2 x^2 - 3 x + 5, x]",
             "List[Rational[31, 8], List[Rule[x, Rational[3, 4]]]]");
}

static void test_quadratic_simple(void) {
    check_eq("Minimize[x^2, x]", "List[0, List[Rule[x, 0]]]");
}

static void test_shifted_parabola(void) {
    check_eq("Minimize[(x - 3)^2 + 1, x]", "List[1, List[Rule[x, 3]]]");
}

static void test_constant(void) {
    check_eq("Minimize[5, x]", "List[5, List[Rule[x, 0]]]");
}

static void test_quartic_value(void) {
    /* x^4 - 2 x^2 -> global min -1 at x = ±1 (one point returned). */
    check_true("First[Minimize[x^4 - 2 x^2, x]] == -1");
    check_true("(x /. Last[Minimize[x^4 - 2 x^2, x]])^2 == 1");
}

static void test_degree6_showcase(void) {
    /* The classic exact-global-min showcase; minimum -35721/64 at x = 11/2. */
    check_eq("Minimize[Expand[(x-1)(x-2)(x-4)(x-7)(x-9)(x-10)], x]",
             "List[Rational[-35721, 64], List[Rule[x, Rational[11, 2]]]]");
}

static void test_unbounded_odd(void) {
    check_true("Minimize[x^3, x] === {-Infinity, {x -> Indeterminate}}");
}

static void test_unbounded_negative_leading(void) {
    check_true("Minimize[-x^2, x] === {-Infinity, {x -> Indeterminate}}");
}

/* ------------------------------------------------------------------ */
/* 3. Maximize                                                         */
/* ------------------------------------------------------------------ */

static void test_maximize_parabola(void) {
    check_eq("Maximize[-x^2, x]", "List[0, List[Rule[x, 0]]]");
}

static void test_maximize_quadratic(void) {
    /* max of -(2x^2-3x+5) = -(min) at the same point. */
    check_eq("Maximize[-(2 x^2 - 3 x + 5), x]",
             "List[Rational[-31, 8], List[Rule[x, Rational[3, 4]]]]");
}

static void test_maximize_unbounded(void) {
    check_true("Maximize[x^2, x] === {Infinity, {x -> Indeterminate}}");
}

static void test_min_max_duality(void) {
    check_true("First[Maximize[-(x - 2)^2 - 1, x]] == -1");
    check_true("(x /. Last[Maximize[-(x - 2)^2 - 1, x]]) == 2");
}

/* ------------------------------------------------------------------ */
/* 4. Multivariate unconstrained polynomial                            */
/* ------------------------------------------------------------------ */

static void test_bowl_2d(void) {
    check_eq("Minimize[x^2 + y^2, {x, y}]",
             "List[0, List[Rule[x, 0], Rule[y, 0]]]");
}

static void test_shifted_bowl_2d(void) {
    check_eq("Minimize[(x - 1)^2 + (y + 2)^2, {x, y}]",
             "List[0, List[Rule[x, 1], Rule[y, -2]]]");
}

static void test_positive_dimensional_declines(void) {
    /* Minimizer set is the curve x^2 - 2y = 1/2 (positive-dimensional): sound
     * decline, left unevaluated. */
    check_eq("Head[Minimize[(x^2 - 2 y)^2 - x^2 + 2 y - 1, {x, y}]]", "Minimize");
}

/* ------------------------------------------------------------------ */
/* 5. Constrained polynomial over the Reals                            */
/* ------------------------------------------------------------------ */

static void test_disk_linear(void) {
    /* {x + y, x^2 + y^2 <= 1} -> -Sqrt[2] at (-1/Sqrt2, -1/Sqrt2). */
    check_true("Abs[N[First[Minimize[{x + y, x^2 + y^2 <= 1}, {x, y}]]] - N[-Sqrt[2]]] < 1*^-9");
    check_true("Abs[N[(x /. Last[Minimize[{x + y, x^2 + y^2 <= 1}, {x, y}]])] + N[1/Sqrt[2]]] < 1*^-9");
}

static void test_equality_inequality(void) {
    /* {x + 2y, x^2 + 2y^2 <= 3 && x + y == 2 && x >= 1} -> 7/3 at (5/3, 1/3). */
    check_eq("Minimize[{x + 2 y, x^2 + 2 y^2 <= 3 && x + y == 2 && x >= 1}, {x, y}]",
             "List[Rational[7, 3], List[Rule[x, Rational[5, 3]], Rule[y, Rational[1, 3]]]]");
}

static void test_linear_program(void) {
    /* {2x + 3y - z, 1<=x+y+z<=2 && 1<=x-y+z<=2 && x-y-z==3} -> 3. */
    check_true("First[Minimize[{2 x + 3 y - z, "
               "1 <= x + y + z <= 2 && 1 <= x - y + z <= 2 && x - y - z == 3}, "
               "{x, y, z}]] == 3");
}

static void test_infeasible(void) {
    check_true("Minimize[{x + y, x^2 < -1}, {x, y}] === "
               "{Infinity, {x -> Indeterminate, y -> Indeterminate}}");
}

static void test_univariate_constrained(void) {
    /* {3x^2 - x + 9, 2x^3 + 5x - 7 >= 0} -> 11 at x = 1. */
    check_eq("Minimize[{3 x^2 - x + 9, 2 x^3 + 5 x - 7 >= 0}, x]",
             "List[11, List[Rule[x, 1]]]");
}

/* ------------------------------------------------------------------ */
/* 6. Inexact input -> NMinimize fallback                              */
/* ------------------------------------------------------------------ */

static void test_inexact_fallback(void) {
    /* 2.5 x^2 - 3 x -> min -0.9 at x = 0.6 (numeric). */
    check_true("Abs[First[Minimize[2.5 x^2 - 3 x, x]] - (-0.9)] < 1*^-3");
    check_true("Head[First[Minimize[2.5 x^2 - 3 x, x]]] === Real");
}

/* ------------------------------------------------------------------ */
/* 7. Sound declines                                                   */
/* ------------------------------------------------------------------ */

static void test_transcendental_declines(void) {
    check_eq("Head[Minimize[Sin[x] + Cos[x], x]]", "Minimize");
}

static void test_unbounded_constrained_declines(void) {
    /* Unbounded below over the feasible region: general unboundedness needs QE
     * (deferred), so this is left unevaluated rather than guessed. */
    check_eq("Head[Minimize[{x + y, x <= y^2}, {x, y}]]", "Minimize");
}

static void test_bad_arity(void) {
    check_eq("Head[Minimize[x^2]]", "Minimize");
}

/* ------------------------------------------------------------------ */
/* 8. Memory smoke                                                     */
/* ------------------------------------------------------------------ */

static void test_memory_smoke(void) {
    for (int i = 0; i < 50; i++) {
        Expr* e = parse_expression("Minimize[2 x^2 - 3 x + 5, x]");
        Expr* r = evaluate(e);
        expr_free(e);
        expr_free(r);
    }
    for (int i = 0; i < 20; i++) {
        Expr* e = parse_expression("Minimize[{x + y, x^2 + y^2 <= 1}, {x, y}]");
        Expr* r = evaluate(e);
        expr_free(e);
        expr_free(r);
    }
}

int main(void) {
    symtab_init();
    core_init();

    /* Intended-error / decline tests write diagnostics to stderr; keep output
     * clean. Successful exact solves are silent. */
    freopen("/dev/null", "w", stderr);

    /* 1. Shape */
    TEST(test_result_shape);

    /* 2. Univariate */
    TEST(test_quadratic);
    TEST(test_quadratic_simple);
    TEST(test_shifted_parabola);
    TEST(test_constant);
    TEST(test_quartic_value);
    TEST(test_degree6_showcase);
    TEST(test_unbounded_odd);
    TEST(test_unbounded_negative_leading);

    /* 3. Maximize */
    TEST(test_maximize_parabola);
    TEST(test_maximize_quadratic);
    TEST(test_maximize_unbounded);
    TEST(test_min_max_duality);

    /* 4. Multivariate unconstrained */
    TEST(test_bowl_2d);
    TEST(test_shifted_bowl_2d);
    TEST(test_positive_dimensional_declines);

    /* 5. Constrained */
    TEST(test_disk_linear);
    TEST(test_equality_inequality);
    TEST(test_linear_program);
    TEST(test_infeasible);
    TEST(test_univariate_constrained);

    /* 6. Inexact fallback */
    TEST(test_inexact_fallback);

    /* 7. Declines */
    TEST(test_transcendental_declines);
    TEST(test_unbounded_constrained_declines);
    TEST(test_bad_arity);

    /* 8. Memory */
    TEST(test_memory_smoke);

    printf("All Minimize tests passed.\n");
    return 0;
}
