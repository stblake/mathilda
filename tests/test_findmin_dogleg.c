/* Unit tests for the dogleg method of FindMinimum / FindMaximum
 * (Method -> "Dogleg", alias "dogleg"; src/numerical_calculus/findmin.c).
 *
 * Dogleg is a TRUST-REGION method that builds the dense Hessian (symbolic, else
 * finite-differenced) and blends the Cauchy (steepest-descent) point with the
 * Newton point along the dogleg path, intersected with the trust radius. It
 * needs a positive-definite model; when the Hessian is indefinite it falls back
 * to the steepest-descent-to-boundary step, so it stays robust on convex and
 * near-convex problems (matching scipy.optimize.minimize(method="dogleg"), which
 * targets the same regime). UNCONSTRAINED: general constraints are rejected.
 * Every pinned target is an ANALYTIC optimum.
 *
 * Run binary directly: ./dogleg_tests */

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

static void check_true(const char* input) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* got = expr_to_string_fullform(res);
    if (strcmp(got, "True") != 0) {
        fprintf(stdout, "FAIL: %s\n  expected: True\n  got:      %s\n", input, got);
        ASSERT_STR_EQ(got, "True");
    }
    free(got);
    expr_free(e);
    expr_free(res);
}

static void check_eq(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* got = expr_to_string_fullform(res);
    if (strcmp(got, expected) != 0) {
        fprintf(stdout, "FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, got);
        ASSERT_STR_EQ(got, expected);
    }
    free(got);
    expr_free(e);
    expr_free(res);
}

/* ------------------------------------------------------------------ */
/* Quadratics — the Newton point lands inside the trust region         */
/* ------------------------------------------------------------------ */

static void test_quadratic_exact(void) {
    check_true("With[{r=FindMinimum[(x-1)^2+(y-2)^2, {{x,0},{y,0}}, Method->\"Dogleg\"]}, "
               "Abs[First[r]] < 1.*^-10 && Abs[(x/.Last[r])-1]+Abs[(y/.Last[r])-2] < 1.*^-6]");
}

static void test_matyas(void) {
    check_true("Abs[First[FindMinimum[0.26(x^2+y^2)-0.48 x y, {{x,3},{y,3}}, "
               "Method->\"Dogleg\", MaxIterations->2000]]] < 1.*^-8");
}

static void test_illcond_n10_c1e4(void) {
    /* condition 1e4 diagonal quadratic: dense Hessian is exact, Newton step. */
    check_true("Module[{v=Table[Symbol[\"z\"<>ToString[i]],{i,1,10}]}, "
               "Abs[First[FindMinimum[Evaluate[Sum[10^(4 (i-1)/9) v[[i]]^2,{i,1,10}]], "
               "Evaluate[Table[{v[[k]],1.0},{k,1,10}]], Method->\"Dogleg\", "
               "MaxIterations->3000]]] < 1.*^-6]");
}

/* ------------------------------------------------------------------ */
/* Curved valleys                                                      */
/* ------------------------------------------------------------------ */

static void test_rosenbrock_2d(void) {
    /* indefinite Hessian en route → the Cauchy fallback keeps progress. */
    check_true("Abs[First[FindMinimum[(1-x)^2+100(y-x^2)^2, {{x,-1.2},{y,1}}, "
               "Method->\"Dogleg\", MaxIterations->3000]]] < 1.*^-6");
    check_true("With[{r=Last[FindMinimum[(1-x)^2+100(y-x^2)^2, {{x,-1.2},{y,1}}, "
               "Method->\"Dogleg\", MaxIterations->3000]]}, "
               "Abs[(x/.r)-1]+Abs[(y/.r)-1] < 1.*^-3]");
}

static void test_booth(void) {
    check_true("With[{r=Last[FindMinimum[(x+2 y-7)^2+(2 x+y-5)^2, {{x,0},{y,0}}, "
               "Method->\"Dogleg\"]]}, Abs[(x/.r)-1]+Abs[(y/.r)-3] < 1.*^-5]");
}

static void test_trid_n6(void) {
    /* Trid n=6, min -50; constant positive-definite Hessian → clean dogleg. */
    check_true("Module[{v=Table[Symbol[\"z\"<>ToString[i]],{i,1,6}]}, "
               "Abs[First[FindMinimum["
               "Evaluate[Sum[(v[[i]]-1)^2,{i,1,6}] - Sum[v[[i]] v[[i-1]],{i,2,6}]], "
               "Evaluate[Table[{v[[k]],0.0},{k,1,6}]], Method->\"Dogleg\", "
               "MaxIterations->2000]] - (-50.0)] < 1.*^-4]");
}

/* ------------------------------------------------------------------ */
/* n==1, parity, maximum, alias, shape, constraint rejection           */
/* ------------------------------------------------------------------ */

static void test_n1(void) {
    /* Sin[x] from x=2: f'' < 0 near the start → Cauchy-to-boundary fallback
     * steers to the minimum 3 Pi/2 ~ 4.712, value -1. */
    check_true("With[{r=FindMinimum[Sin[x], {x,2}, Method->\"Dogleg\"]}, "
               "Abs[First[r]-(-1.0)] < 1.*^-6 && Abs[(x/.Last[r])-4.712388980] < 1.*^-4]");
}

static void test_parity_quasinewton(void) {
    check_true("With[{a=Last[FindMinimum[(x-1)^2+(y-2)^2+(x y-1)^2, {{x,0},{y,0}}, Method->\"Dogleg\"]], "
               "b=Last[FindMinimum[(x-1)^2+(y-2)^2+(x y-1)^2, {{x,0},{y,0}}, Method->\"QuasiNewton\"]]}, "
               "Abs[(x/.a)-(x/.b)]+Abs[(y/.a)-(y/.b)] < 1.*^-4]");
}

static void test_findmaximum(void) {
    check_true("With[{r=FindMaximum[-(x-3)^2-(y+1)^2, {{x,0},{y,0}}, Method->\"Dogleg\"]}, "
               "Abs[First[r]] < 1.*^-8 && Abs[(x/.Last[r])-3]+Abs[(y/.Last[r])+1] < 1.*^-4]");
}

static void test_alias(void) {
    check_true("First[FindMinimum[(x-1)^2+(y-2)^2, {{x,0},{y,0}}, Method->\"dogleg\"]] == "
               "First[FindMinimum[(x-1)^2+(y-2)^2, {{x,0},{y,0}}, Method->\"Dogleg\"]]");
}

static void test_shape(void) {
    check_eq("Length[FindMinimum[(x-1)^2+(y-2)^2, {{x,0},{y,0}}, Method->\"Dogleg\"]]", "2");
    check_eq("Length[Last[FindMinimum[(x-1)^2+(y-2)^2, {{x,0},{y,0}}, Method->\"Dogleg\"]]]", "2");
}

static void test_reject_general_constraint(void) {
    check_eq("Head[FindMinimum[{x^2+y^2, x+y>=1}, {{x,0},{y,0}}, Method->\"Dogleg\"]]",
             "FindMinimum");
}

/* ------------------------------------------------------------------ */
/* Memory hygiene smoke                                                */
/* ------------------------------------------------------------------ */

static void test_no_leak_many_calls(void) {
    for (int i = 0; i < 30; i++) {
        Expr* e = parse_expression(
            "FindMinimum[(x-2)^2+(y-3)^2+0.1(x y-1)^2, {{x,0},{y,0}}, Method->\"Dogleg\"]");
        Expr* r = evaluate(e);
        expr_free(e);
        expr_free(r);
    }
    check_true("True");
}

/* ------------------------------------------------------------------ */

int main(void) {
    symtab_init();
    core_init();
    freopen("/dev/null", "w", stderr);

    TEST(test_quadratic_exact);
    TEST(test_matyas);
    TEST(test_illcond_n10_c1e4);
    TEST(test_rosenbrock_2d);
    TEST(test_booth);
    TEST(test_trid_n6);
    TEST(test_n1);
    TEST(test_parity_quasinewton);
    TEST(test_findmaximum);
    TEST(test_alias);
    TEST(test_shape);
    TEST(test_reject_general_constraint);
    TEST(test_no_leak_many_calls);

    printf("All dogleg tests passed!\n");
    return 0;
}
