/* Unit tests for the COBYLA method of FindMinimum / FindMaximum
 * (Method -> "COBYLA"; src/numerical_calculus/findmin.c).
 *
 * COBYLA (Powell's Constrained Optimization BY Linear Approximation) is a
 * DERIVATIVE-FREE trust-region method for CONSTRAINED optimization -- the first
 * derivative-free method in Mathilda to accept general (non-box) constraints
 * (Powell/NelderMead reject them; SLSQP/TNC need a gradient). It models the
 * objective and every constraint by a linear approximation (central differences
 * on a coordinate cross), solves a two-stage trust-region LP (feasibility, then
 * objective) via the same active-set solver as SLSQP, and accepts steps by an
 * L-infinity exact-penalty merit. Equalities are handled by splitting h==0 into
 * h<=0 and -h<=0 (so this is strictly more capable than scipy's inequality-only
 * COBYLA); box bounds enter as ordinary linear inequalities.
 *
 * COBYLA converges to ~rhoend = 10^-PrecisionGoal accuracy (a trust-region
 * radius, not a gradient test), so value pins use ~1e-4 tolerances. Its
 * distinguishing strength -- exercised here -- is NON-SMOOTH constrained
 * objectives (L1 / max-of-affine), where the gradient methods cannot help.
 * Every case fixes its start explicitly (COBYLA is deterministic) and pins an
 * analytic optimum.
 *
 * Run binary directly: ./cobyla_tests */

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
/* C1. Unconstrained reduction (smooth + non-smooth)                   */
/* ------------------------------------------------------------------ */

static void test_uncon_quadratic(void) {
    /* min 0 at (1,2). */
    check_true("Abs[First[FindMinimum[(x-1)^2+(y-2)^2, {{x,0},{y,0}}, "
               "Method->\"COBYLA\", MaxIterations->2000]]] < 1.*^-6");
    check_true("With[{r=Last[FindMinimum[(x-1)^2+(y-2)^2, {{x,0},{y,0}}, "
               "Method->\"COBYLA\", MaxIterations->2000]]}, "
               "Abs[(x/.r)-1]+Abs[(y/.r)-2] < 1.*^-3]");
}

static void test_uncon_nonsmooth(void) {
    /* min |x-3|+|y+2| -> 0 at (3,-2). Derivative methods can't touch this;
     * COBYLA's linear-approximation trust region does. */
    check_true("First[FindMinimum[Abs[x-3]+Abs[y+2], {{x,0},{y,0}}, "
               "Method->\"COBYLA\", MaxIterations->2000]] < 1.*^-3");
    check_true("With[{r=Last[FindMinimum[Abs[x-3]+Abs[y+2], {{x,0},{y,0}}, "
               "Method->\"COBYLA\", MaxIterations->2000]]}, "
               "Abs[(x/.r)-3]+Abs[(y/.r)+2] < 1.*^-2]");
}

/* ------------------------------------------------------------------ */
/* C2. Inequality-constrained (analytic KKT)                           */
/* ------------------------------------------------------------------ */

static void test_ineq_corner(void) {
    /* min x+y s.t. 3x+2y>=7, x>=0, y>=0 -> 7/3 at (7/3,0). */
    check_true("Abs[First[FindMinimum[{x+y, 3 x+2 y>=7 && x>=0 && y>=0}, "
               "{{x,1},{y,1}}, Method->\"COBYLA\", MaxIterations->2000]] - 7/3] < 1.*^-4");
}

static void test_ineq_halfplane(void) {
    /* min x^2+y^2 s.t. x>=1 -> 1 at (1,0). */
    check_true("Abs[First[FindMinimum[{x^2+y^2, x>=1}, {{x,2},{y,1}}, "
               "Method->\"COBYLA\", MaxIterations->2000]] - 1.] < 1.*^-4");
}

static void test_ineq_offset(void) {
    /* min (x-2)^2+(y-1)^2 s.t. x+y<=1 -> 2 at (1,0). */
    check_true("Abs[First[FindMinimum[{(x-2)^2+(y-1)^2, x+y<=1}, {{x,0},{y,0}}, "
               "Method->\"COBYLA\", MaxIterations->2000]] - 2.] < 1.*^-4");
}

static void test_qp_tutorial(void) {
    /* scipy's COBYLA/SLSQP tutorial QP: f* = 0.8 at (1.4,1.7). */
    check_true("Abs[First[FindMinimum[{(x-1)^2+(y-2.5)^2, "
               "x-2 y+2>=0 && -x-2 y+6>=0 && -x+2 y+2>=0 && x>=0 && y>=0}, "
               "{{x,2},{y,0}}, Method->\"COBYLA\", MaxIterations->2000]] - 0.8] < 1.*^-4");
}

/* ------------------------------------------------------------------ */
/* C3. Equality (handled by the h<=0 && -h<=0 split)                   */
/* ------------------------------------------------------------------ */

static void test_eq_line(void) {
    /* min x^2+y^2 s.t. x+y==1 -> 1/2 at (1/2,1/2). */
    check_true("Abs[First[FindMinimum[{x^2+y^2, x+y==1}, {{x,0},{y,0}}, "
               "Method->\"COBYLA\", MaxIterations->2000]] - 0.5] < 1.*^-4");
    check_true("With[{r=Last[FindMinimum[{x^2+y^2, x+y==1}, {{x,0},{y,0}}, "
               "Method->\"COBYLA\", MaxIterations->2000]]}, "
               "Abs[(x/.r)-0.5]+Abs[(y/.r)-0.5] < 1.*^-2]");
}

static void test_eq_circle(void) {
    /* min x+y s.t. x^2+y^2==1 -> -Sqrt[2]. */
    check_true("Abs[First[FindMinimum[{x+y, x^2+y^2==1}, {{x,-1},{y,-1}}, "
               "Method->\"COBYLA\", MaxIterations->2000]] + Sqrt[2.]] < 1.*^-4");
}

/* ------------------------------------------------------------------ */
/* C4. Non-smooth CONSTRAINED -- COBYLA's niche                        */
/* ------------------------------------------------------------------ */

static void test_nonsmooth_constrained(void) {
    /* min |x-2|+|y+1| s.t. x+y>=0 -> unconstrained min (2,-1) is feasible
     * (2-1=1>=0), so f*=0 at (2,-1). The kink at the optimum is why central
     * differences (not forward) are needed. */
    check_true("First[FindMinimum[{Abs[x-2]+Abs[y+1], x+y>=0}, {{x,0},{y,0}}, "
               "Method->\"COBYLA\", MaxIterations->2000]] < 1.*^-3");
}

static void test_nonsmooth_minimax(void) {
    /* min max(x, y, -x-y) s.t. (unconstrained) -> the Chebyshev center 0 at
     * (0,0)? min over R^2 of max(x,y,-x-y): the min of the max-of-affine is 0
     * at the origin (each affine is 0 there and they can't all be negative).
     * Expressed with Max. */
    check_true("First[FindMinimum[Max[x, y, -x-y], {{x,1},{y,1}}, "
               "Method->\"COBYLA\", MaxIterations->2000]] < 1.*^-3");
}

/* ------------------------------------------------------------------ */
/* C5. Box-only bounds                                                 */
/* ------------------------------------------------------------------ */

static void test_box_corner(void) {
    /* unconstrained min (2,3) outside [0,1]^2 -> corner (1,1), value 5. */
    check_true("Abs[First[FindMinimum[{(x-2)^2+(y-3)^2, 0<=x<=1 && 0<=y<=1}, "
               "{{x,0.5},{y,0.5}}, Method->\"COBYLA\", MaxIterations->2000]] - 5.] < 1.*^-4");
}

/* ------------------------------------------------------------------ */
/* C6. FindMaximum                                                     */
/* ------------------------------------------------------------------ */

static void test_findmaximum(void) {
    /* max x y s.t. x+y==10 -> 25 at (5,5). */
    check_true("Abs[First[FindMaximum[{x y, x+y==10}, {{x,1},{y,1}}, "
               "Method->\"COBYLA\", MaxIterations->2000]] - 25.] < 1.*^-3");
}

/* ------------------------------------------------------------------ */
/* C7. Edge cases                                                      */
/* ------------------------------------------------------------------ */

static void test_n1(void) {
    /* n==1 unconstrained delegates to Brent: min (x-3)^2 -> 0 at 3. */
    check_true("Abs[First[FindMinimum[(x-3)^2, {x,0}, Method->\"COBYLA\"]]] < 1.*^-8");
    /* n==1 WITH a constraint uses the full COBYLA machinery. */
    check_true("Abs[First[FindMinimum[{(x-3)^2, x>=1}, {x,0}, "
               "Method->\"COBYLA\", MaxIterations->2000]]] < 1.*^-4");
}

static void test_inconsistent_no_crash(void) {
    check_eq("Head[FindMinimum[{x, x==0 && x==1}, {{x,0.5}}, Method->\"COBYLA\"]]",
             "List");
}

static void test_shape(void) {
    check_true("MatchQ[FindMinimum[{x^2+y^2, x+y==1}, {{x,0},{y,0}}, "
               "Method->\"COBYLA\"], {_?NumberQ, {(_->_)..}}]");
}

static void test_no_leak_many_calls(void) {
    check_true("Module[{s=0.}, Do[s += First[FindMinimum[{x+y, 3 x+2 y>=7 && x>=0 && y>=0}, "
               "{{x,1},{y,1}}, Method->\"COBYLA\", MaxIterations->1000]], {40}]; "
               "Abs[s - 40*7/3] < 1.*^-2]");
}

int main(void) {
    symtab_init();
    core_init();
    freopen("/dev/null", "w", stderr);

    /* C1 unconstrained */
    TEST(test_uncon_quadratic);
    TEST(test_uncon_nonsmooth);

    /* C2 inequality-constrained */
    TEST(test_ineq_corner);
    TEST(test_ineq_halfplane);
    TEST(test_ineq_offset);
    TEST(test_qp_tutorial);

    /* C3 equality via split */
    TEST(test_eq_line);
    TEST(test_eq_circle);

    /* C4 non-smooth constrained */
    TEST(test_nonsmooth_constrained);
    TEST(test_nonsmooth_minimax);

    /* C5 box-only */
    TEST(test_box_corner);

    /* C6 FindMaximum */
    TEST(test_findmaximum);

    /* C7 edge cases */
    TEST(test_n1);
    TEST(test_inconsistent_no_crash);
    TEST(test_shape);
    TEST(test_no_leak_many_calls);

    printf("All COBYLA tests passed!\n");
    return 0;
}
