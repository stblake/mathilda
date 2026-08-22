/* Unit tests for the Newton-CG method of FindMinimum / FindMaximum
 * (Method -> "NewtonCG", alias "Newton-CG"; src/numerical_calculus/findmin.c).
 *
 * Newton-CG is a LINE-SEARCH truncated-Newton method: an inner conjugate-gradient
 * loop approximately solves the Newton system H*p = -g using Hessian-vector
 * products (finite differences of the exact compiled gradient) with the
 * Eisenstat-Walker forcing sequence and negative-curvature truncation, followed
 * by a unit-step-first Wolfe line search. It is UNCONSTRAINED, mirroring
 * scipy.optimize.minimize(method="Newton-CG"): general (non-box) constraints are
 * rejected, like Powell/NelderMead. Every pinned target is an ANALYTIC optimum.
 *
 * Large objectives use genuine scalar symbols generated with
 * Symbol["z"<>ToString[i]], wrapped in Evaluate[...] past FindMinimum's HoldAll.
 *
 * Run binary directly: ./newtoncg_tests */

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
/* Curved valleys — true Newton curvature via Hv                       */
/* ------------------------------------------------------------------ */

static void test_quadratic_exact(void) {
    /* min 0 at (1,2); the inner CG solves the 2x2 Newton system exactly. */
    check_true("With[{r=FindMinimum[(x-1)^2+(y-2)^2, {{x,0},{y,0}}, Method->\"NewtonCG\"]}, "
               "Abs[First[r]] < 1.*^-10 && Abs[(x/.Last[r])-1]+Abs[(y/.Last[r])-2] < 1.*^-6]");
}

static void test_rosenbrock_2d(void) {
    check_true("Abs[First[FindMinimum[(1-x)^2+100(y-x^2)^2, {{x,-1.2},{y,1}}, "
               "Method->\"NewtonCG\", MaxIterations->2000]]] < 1.*^-8");
    check_true("With[{r=Last[FindMinimum[(1-x)^2+100(y-x^2)^2, {{x,-1.2},{y,1}}, "
               "Method->\"NewtonCG\", MaxIterations->2000]]}, "
               "Abs[(x/.r)-1]+Abs[(y/.r)-1] < 1.*^-4]");
}

static void test_rosenbrock_extended_n10(void) {
    check_true("Module[{v=Table[Symbol[\"z\"<>ToString[i]],{i,1,10}]}, "
               "Abs[First[FindMinimum["
               "Evaluate[Sum[100(v[[j+1]]-v[[j]]^2)^2+(1-v[[j]])^2,{j,1,9}]], "
               "Evaluate[Table[{v[[k]],-1.2},{k,1,10}]], Method->\"NewtonCG\", "
               "MaxIterations->4000]]] < 1.*^-4]");
}

/* ------------------------------------------------------------------ */
/* Ill-conditioned quadratics                                          */
/* ------------------------------------------------------------------ */

static void test_illcond_n10_c1e4(void) {
    check_true("Module[{v=Table[Symbol[\"z\"<>ToString[i]],{i,1,10}]}, "
               "Abs[First[FindMinimum[Evaluate[Sum[10^(4 (i-1)/9) v[[i]]^2,{i,1,10}]], "
               "Evaluate[Table[{v[[k]],1.0},{k,1,10}]], Method->\"NewtonCG\", "
               "MaxIterations->3000]]] < 1.*^-6]");
}

/* ------------------------------------------------------------------ */
/* Standard smooth functions                                           */
/* ------------------------------------------------------------------ */

static void test_booth(void) {
    check_true("With[{r=Last[FindMinimum[(x+2 y-7)^2+(2 x+y-5)^2, {{x,0},{y,0}}, "
               "Method->\"NewtonCG\"]]}, Abs[(x/.r)-1]+Abs[(y/.r)-3] < 1.*^-5]");
}

static void test_matyas(void) {
    check_true("Abs[First[FindMinimum[0.26(x^2+y^2)-0.48 x y, {{x,3},{y,3}}, "
               "Method->\"NewtonCG\", MaxIterations->2000]]] < 1.*^-8");
}

static void test_trid_n6(void) {
    /* Trid n=6, analytic global min -50 (Hessian constant, positive definite). */
    check_true("Module[{v=Table[Symbol[\"z\"<>ToString[i]],{i,1,6}]}, "
               "Abs[First[FindMinimum["
               "Evaluate[Sum[(v[[i]]-1)^2,{i,1,6}] - Sum[v[[i]] v[[i-1]],{i,2,6}]], "
               "Evaluate[Table[{v[[k]],0.0},{k,1,6}]], Method->\"NewtonCG\", "
               "MaxIterations->2000]] - (-50.0)] < 1.*^-4]");
}

static void test_highdim_quadratic_n20(void) {
    /* min 0 at z_i = i; Hessian-free method scales without forming H. */
    check_true("Module[{v=Table[Symbol[\"z\"<>ToString[i]],{i,1,20}]}, "
               "Abs[First[FindMinimum[Evaluate[Sum[(v[[i]]-i)^2,{i,1,20}]], "
               "Evaluate[Table[{v[[k]],0.0},{k,1,20}]], Method->\"NewtonCG\"]]] < 1.*^-8]");
}

/* ------------------------------------------------------------------ */
/* n==1, parity, maximum, alias, shape, constraint rejection           */
/* ------------------------------------------------------------------ */

static void test_n1(void) {
    /* Sin[x] from x=2 → 3 Pi/2 ~ 4.712, value -1 (runs the general machinery,
     * no Brent delegation when the method is requested explicitly). */
    check_true("With[{r=FindMinimum[Sin[x], {x,2}, Method->\"NewtonCG\"]}, "
               "Abs[First[r]-(-1.0)] < 1.*^-6 && Abs[(x/.Last[r])-4.712388980] < 1.*^-4]");
}

static void test_parity_quasinewton(void) {
    check_true("With[{a=Last[FindMinimum[(x-1)^2+(y-2)^2+(x y-1)^2, {{x,0},{y,0}}, Method->\"NewtonCG\"]], "
               "b=Last[FindMinimum[(x-1)^2+(y-2)^2+(x y-1)^2, {{x,0},{y,0}}, Method->\"QuasiNewton\"]]}, "
               "Abs[(x/.a)-(x/.b)]+Abs[(y/.a)-(y/.b)] < 1.*^-4]");
}

static void test_findmaximum(void) {
    check_true("With[{r=FindMaximum[-(x-3)^2-(y+1)^2, {{x,0},{y,0}}, Method->\"NewtonCG\"]}, "
               "Abs[First[r]] < 1.*^-8 && Abs[(x/.Last[r])-3]+Abs[(y/.Last[r])+1] < 1.*^-4]");
}

static void test_alias(void) {
    /* "Newton-CG" (scipy spelling) routes to the identical solver. */
    check_true("First[FindMinimum[(1-x)^2+100(y-x^2)^2, {{x,-1.2},{y,1}}, "
               "Method->\"Newton-CG\", MaxIterations->2000]] == "
               "First[FindMinimum[(1-x)^2+100(y-x^2)^2, {{x,-1.2},{y,1}}, "
               "Method->\"NewtonCG\", MaxIterations->2000]]");
}

static void test_shape(void) {
    check_eq("Length[FindMinimum[(x-1)^2+(y-2)^2, {{x,0},{y,0}}, Method->\"NewtonCG\"]]", "2");
    check_eq("Length[Last[FindMinimum[(x-1)^2+(y-2)^2, {{x,0},{y,0}}, Method->\"NewtonCG\"]]]", "2");
}

static void test_reject_general_constraint(void) {
    /* Unconstrained method (mirrors scipy): a general constraint is rejected and
     * the expression is returned unevaluated (head stays FindMinimum). */
    check_eq("Head[FindMinimum[{x^2+y^2, x+y>=1}, {{x,0},{y,0}}, Method->\"NewtonCG\"]]",
             "FindMinimum");
}

/* ------------------------------------------------------------------ */
/* Memory hygiene smoke                                                */
/* ------------------------------------------------------------------ */

static void test_no_leak_many_calls(void) {
    for (int i = 0; i < 30; i++) {
        Expr* e = parse_expression(
            "FindMinimum[(x-2)^2+(y-3)^2+(x y-1)^2, {{x,0},{y,0}}, Method->\"NewtonCG\"]");
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
    TEST(test_rosenbrock_2d);
    TEST(test_rosenbrock_extended_n10);
    TEST(test_illcond_n10_c1e4);
    TEST(test_booth);
    TEST(test_matyas);
    TEST(test_trid_n6);
    TEST(test_highdim_quadratic_n20);
    TEST(test_n1);
    TEST(test_parity_quasinewton);
    TEST(test_findmaximum);
    TEST(test_alias);
    TEST(test_shape);
    TEST(test_reject_general_constraint);
    TEST(test_no_leak_many_calls);

    printf("All Newton-CG tests passed!\n");
    return 0;
}
