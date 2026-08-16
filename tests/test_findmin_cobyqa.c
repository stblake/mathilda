/* Unit tests for the COBYQA method of FindMinimum / FindMaximum
 * (Method -> "COBYQA"; src/numerical_calculus/findmin.c).
 *
 * COBYQA (Constrained Optimization BY Quadratic Approximations; Ragonneau &
 * Zhang 2023) is a DERIVATIVE-FREE trust-region SQP with QUADRATIC interpolation
 * models. Where COBYLA models the objective and constraints by linear
 * approximations, COBYQA models each by a full quadratic (built by finite
 * differences on a structured stencil), so it captures curvature: it converges
 * tighter on smooth problems and -- the discriminating case here -- navigates
 * the curved Rosenbrock valley to machine precision, where COBYLA's linear
 * models stall near f ~ 1e-3. It handles equality + inequality + bound
 * constraints natively (equalities are NOT split). The trust-region SQP step
 * reuses the SLSQP dual active-set QP; acceptance uses an L1 penalty merit.
 *
 * Every case fixes its start explicitly (COBYQA is deterministic) and pins an
 * analytic optimum. Run binary directly: ./cobyqa_tests */

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
/* Q1. Curvature -- quadratic models beat COBYLA's linear ones          */
/* ------------------------------------------------------------------ */

static void test_rosenbrock(void) {
    /* min 0 at (1,1). The curved valley COBYLA's LINEAR models cannot navigate
     * to machine precision -- COBYQA's quadratic models reach it. */
    check_true("Abs[First[FindMinimum[(1-x)^2+100(y-x^2)^2, {{x,-1.2},{y,1}}, "
               "Method->\"COBYQA\", MaxIterations->3000]]] < 1.*^-6");
    check_true("With[{r=Last[FindMinimum[(1-x)^2+100(y-x^2)^2, {{x,-1.2},{y,1}}, "
               "Method->\"COBYQA\", MaxIterations->3000]]}, "
               "Abs[(x/.r)-1]+Abs[(y/.r)-1] < 1.*^-3]");
}

static void test_uncon_quadratic(void) {
    check_true("Abs[First[FindMinimum[(x-1)^2+(y-2)^2, {{x,0},{y,0}}, "
               "Method->\"COBYQA\", MaxIterations->2000]]] < 1.*^-8");
}

/* ------------------------------------------------------------------ */
/* Q2. Inequality-constrained (analytic KKT)                           */
/* ------------------------------------------------------------------ */

static void test_ineq_corner(void) {
    /* min x+y s.t. 3x+2y>=7, x>=0, y>=0 -> 7/3 at (7/3,0). */
    check_true("Abs[First[FindMinimum[{x+y, 3 x+2 y>=7 && x>=0 && y>=0}, "
               "{{x,1},{y,1}}, Method->\"COBYQA\", MaxIterations->2000]] - 7/3] < 1.*^-4");
}

static void test_ineq_halfplane(void) {
    check_true("Abs[First[FindMinimum[{x^2+y^2, x>=1}, {{x,2},{y,1}}, "
               "Method->\"COBYQA\", MaxIterations->2000]] - 1.] < 1.*^-5");
}

static void test_qp_tutorial(void) {
    check_true("Abs[First[FindMinimum[{(x-1)^2+(y-2.5)^2, "
               "x-2 y+2>=0 && -x-2 y+6>=0 && -x+2 y+2>=0 && x>=0 && y>=0}, "
               "{{x,2},{y,0}}, Method->\"COBYQA\", MaxIterations->2000]] - 0.8] < 1.*^-5");
}

/* ------------------------------------------------------------------ */
/* Q3. Equality -- handled NATIVELY (not split)                        */
/* ------------------------------------------------------------------ */

static void test_eq_line(void) {
    check_true("Abs[First[FindMinimum[{x^2+y^2, x+y==1}, {{x,0},{y,0}}, "
               "Method->\"COBYQA\", MaxIterations->2000]] - 0.5] < 1.*^-5");
    check_true("With[{r=Last[FindMinimum[{x^2+y^2, x+y==1}, {{x,0},{y,0}}, "
               "Method->\"COBYQA\", MaxIterations->2000]]}, "
               "Abs[(x/.r)-0.5]+Abs[(y/.r)-0.5] < 1.*^-3]");
}

static void test_eq_circle(void) {
    /* nonlinear equality -- native handling with constraint curvature in H_L. */
    check_true("Abs[First[FindMinimum[{x+y, x^2+y^2==1}, {{x,-1},{y,-1}}, "
               "Method->\"COBYQA\", MaxIterations->2000]] + Sqrt[2.]] < 1.*^-5");
}

/* ------------------------------------------------------------------ */
/* Q4. Mixed eq + ineq + bounds -- HS71 (derivative-free)              */
/* ------------------------------------------------------------------ */

static void test_hs71(void) {
    /* HS71 = 17.0140 at (1, 4.743, 3.821, 1.379) -- reached WITHOUT a gradient. */
    check_true("Abs[First[FindMinimum[{x1 x4 (x1+x2+x3)+x3, "
               "x1 x2 x3 x4>=25 && x1^2+x2^2+x3^2+x4^2==40 && "
               "1<=x1<=5 && 1<=x2<=5 && 1<=x3<=5 && 1<=x4<=5}, "
               "{{x1,1},{x2,5},{x3,5},{x4,1}}, Method->\"COBYQA\", "
               "MaxIterations->3000]] - 17.0140173] < 1.*^-2");
}

/* ------------------------------------------------------------------ */
/* Q5. Box-only bounds                                                 */
/* ------------------------------------------------------------------ */

static void test_box_corner(void) {
    check_true("Abs[First[FindMinimum[{(x-2)^2+(y-3)^2, 0<=x<=1 && 0<=y<=1}, "
               "{{x,0.5},{y,0.5}}, Method->\"COBYQA\", MaxIterations->2000]] - 5.] < 1.*^-5");
}

/* ------------------------------------------------------------------ */
/* Q6. FindMaximum                                                     */
/* ------------------------------------------------------------------ */

static void test_findmaximum(void) {
    check_true("Abs[First[FindMaximum[{x y, x+y==10}, {{x,1},{y,1}}, "
               "Method->\"COBYQA\", MaxIterations->2000]] - 25.] < 1.*^-4");
}

/* ------------------------------------------------------------------ */
/* Q7. Edge cases                                                      */
/* ------------------------------------------------------------------ */

static void test_n1(void) {
    check_true("Abs[First[FindMinimum[(x-3)^2, {x,0}, Method->\"COBYQA\"]]] < 1.*^-8");
    check_true("Abs[First[FindMinimum[{(x-3)^2, x>=1}, {x,0}, "
               "Method->\"COBYQA\", MaxIterations->2000]]] < 1.*^-4");
}

static void test_inconsistent_no_crash(void) {
    check_eq("Head[FindMinimum[{x, x==0 && x==1}, {{x,0.5}}, Method->\"COBYQA\"]]",
             "List");
}

static void test_shape(void) {
    check_true("MatchQ[FindMinimum[{x^2+y^2, x+y==1}, {{x,0},{y,0}}, "
               "Method->\"COBYQA\"], {_?NumberQ, {(_->_)..}}]");
}

static void test_no_leak_many_calls(void) {
    check_true("Module[{s=0.}, Do[s += First[FindMinimum[{x^2+y^2, x+y==1}, "
               "{{x,0},{y,0}}, Method->\"COBYQA\", MaxIterations->1000]], {40}]; "
               "Abs[s - 40*0.5] < 1.*^-2]");
}

int main(void) {
    symtab_init();
    core_init();
    freopen("/dev/null", "w", stderr);

    /* Q1 curvature */
    TEST(test_rosenbrock);
    TEST(test_uncon_quadratic);

    /* Q2 inequality-constrained */
    TEST(test_ineq_corner);
    TEST(test_ineq_halfplane);
    TEST(test_qp_tutorial);

    /* Q3 equality (native) */
    TEST(test_eq_line);
    TEST(test_eq_circle);

    /* Q4 mixed HS71 */
    TEST(test_hs71);

    /* Q5 box-only */
    TEST(test_box_corner);

    /* Q6 FindMaximum */
    TEST(test_findmaximum);

    /* Q7 edge cases */
    TEST(test_n1);
    TEST(test_inconsistent_no_crash);
    TEST(test_shape);
    TEST(test_no_leak_many_calls);

    printf("All COBYQA tests passed!\n");
    return 0;
}
