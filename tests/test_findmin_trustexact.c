/* Unit tests for the trust-exact method of FindMinimum / FindMaximum
 * (Method -> "TrustExact", alias "trust-exact"; src/numerical_calculus/findmin.c).
 *
 * trust-exact solves the trust-region subproblem NEARLY EXACTLY by the
 * Moré-Sorensen algorithm: it builds the dense Hessian and iterates on the
 * Levenberg shift λ so that (B+λI) is positive semidefinite and ||p(λ)|| ≈ Δ,
 * with a dedicated hard-case branch (g ⟂ the least eigenvector) that adds the
 * least-eigenvector direction out to the boundary. It therefore handles
 * indefinite Hessians that dogleg cannot, matching
 * scipy.optimize.minimize(method="trust-exact"). UNCONSTRAINED: general
 * constraints are rejected. Every pinned target is an ANALYTIC optimum.
 *
 * Run binary directly: ./trustexact_tests */

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
/* Quadratics — the exact solver reaches machine precision             */
/* ------------------------------------------------------------------ */

static void test_quadratic_exact(void) {
    check_true("With[{r=FindMinimum[(x-1)^2+(y-2)^2, {{x,0},{y,0}}, Method->\"TrustExact\"]}, "
               "Abs[First[r]] < 1.*^-12 && Abs[(x/.Last[r])-1]+Abs[(y/.Last[r])-2] < 1.*^-7]");
}

static void test_illcond_n10_c1e6(void) {
    /* condition 1e6: the exact subproblem solve reaches ~0 where a low-rank
     * quasi-Newton model would crawl. */
    check_true("Module[{v=Table[Symbol[\"z\"<>ToString[i]],{i,1,10}]}, "
               "Abs[First[FindMinimum[Evaluate[Sum[10^(6 (i-1)/9) v[[i]]^2,{i,1,10}]], "
               "Evaluate[Table[{v[[k]],1.0},{k,1,10}]], Method->\"TrustExact\", "
               "MaxIterations->3000]]] < 1.*^-5]");
}

/* ------------------------------------------------------------------ */
/* Indefinite Hessians — where dogleg's SPD assumption fails            */
/* ------------------------------------------------------------------ */

static void test_double_well_negcurv(void) {
    /* (x^2-1)^2 + y^2 from x≈0: indefinite Hessian at the start. The exact
     * subproblem shifts λ past -λmin and steps to the boundary. min 0 at x=±1. */
    check_true("Abs[First[FindMinimum[(x^2-1)^2+y^2, {{x,0.1},{y,0.5}}, "
               "Method->\"TrustExact\", MaxIterations->2000]]] < 1.*^-10");
}

static void test_beale(void) {
    /* Beale: multimodal with indefinite regions, min 0 at (3, 0.5). */
    check_true("With[{r=FindMinimum["
               "(1.5-x+x y)^2+(2.25-x+x y^2)^2+(2.625-x+x y^3)^2, {{x,1},{y,1}}, "
               "Method->\"TrustExact\", MaxIterations->2000]}, "
               "Abs[First[r]] < 1.*^-8 && Abs[(x/.Last[r])-3]+Abs[(y/.Last[r])-0.5] < 1.*^-4]");
}

/* ------------------------------------------------------------------ */
/* Curved valleys and standard functions                               */
/* ------------------------------------------------------------------ */

static void test_rosenbrock_2d(void) {
    check_true("Abs[First[FindMinimum[(1-x)^2+100(y-x^2)^2, {{x,-1.2},{y,1}}, "
               "Method->\"TrustExact\", MaxIterations->2000]]] < 1.*^-8");
    check_true("With[{r=Last[FindMinimum[(1-x)^2+100(y-x^2)^2, {{x,-1.2},{y,1}}, "
               "Method->\"TrustExact\", MaxIterations->2000]]}, "
               "Abs[(x/.r)-1]+Abs[(y/.r)-1] < 1.*^-4]");
}

static void test_booth(void) {
    check_true("With[{r=Last[FindMinimum[(x+2 y-7)^2+(2 x+y-5)^2, {{x,0},{y,0}}, "
               "Method->\"TrustExact\"]]}, Abs[(x/.r)-1]+Abs[(y/.r)-3] < 1.*^-5]");
}

static void test_trid_n6(void) {
    check_true("Module[{v=Table[Symbol[\"z\"<>ToString[i]],{i,1,6}]}, "
               "Abs[First[FindMinimum["
               "Evaluate[Sum[(v[[i]]-1)^2,{i,1,6}] - Sum[v[[i]] v[[i-1]],{i,2,6}]], "
               "Evaluate[Table[{v[[k]],0.0},{k,1,6}]], Method->\"TrustExact\", "
               "MaxIterations->2000]] - (-50.0)] < 1.*^-4]");
}

/* ------------------------------------------------------------------ */
/* n==1, parity, maximum, alias, shape, constraint rejection           */
/* ------------------------------------------------------------------ */

static void test_n1(void) {
    check_true("With[{r=FindMinimum[Sin[x], {x,2}, Method->\"TrustExact\"]}, "
               "Abs[First[r]-(-1.0)] < 1.*^-6 && Abs[(x/.Last[r])-4.712388980] < 1.*^-4]");
}

static void test_parity_quasinewton(void) {
    check_true("With[{a=Last[FindMinimum[(x-1)^2+(y-2)^2+(x y-1)^2, {{x,0},{y,0}}, Method->\"TrustExact\"]], "
               "b=Last[FindMinimum[(x-1)^2+(y-2)^2+(x y-1)^2, {{x,0},{y,0}}, Method->\"QuasiNewton\"]]}, "
               "Abs[(x/.a)-(x/.b)]+Abs[(y/.a)-(y/.b)] < 1.*^-4]");
}

static void test_findmaximum(void) {
    check_true("With[{r=FindMaximum[-(x-3)^2-(y+1)^2, {{x,0},{y,0}}, Method->\"TrustExact\"]}, "
               "Abs[First[r]] < 1.*^-8 && Abs[(x/.Last[r])-3]+Abs[(y/.Last[r])+1] < 1.*^-4]");
}

static void test_alias(void) {
    check_true("First[FindMinimum[(x-1)^2+(y-2)^2, {{x,0},{y,0}}, Method->\"trust-exact\"]] == "
               "First[FindMinimum[(x-1)^2+(y-2)^2, {{x,0},{y,0}}, Method->\"TrustExact\"]]");
}

static void test_shape(void) {
    check_eq("Length[FindMinimum[(x-1)^2+(y-2)^2, {{x,0},{y,0}}, Method->\"TrustExact\"]]", "2");
    check_eq("Length[Last[FindMinimum[(x-1)^2+(y-2)^2, {{x,0},{y,0}}, Method->\"TrustExact\"]]]", "2");
}

static void test_reject_general_constraint(void) {
    check_eq("Head[FindMinimum[{x^2+y^2, x+y>=1}, {{x,0},{y,0}}, Method->\"TrustExact\"]]",
             "FindMinimum");
}

/* ------------------------------------------------------------------ */
/* Memory hygiene smoke                                                */
/* ------------------------------------------------------------------ */

static void test_no_leak_many_calls(void) {
    for (int i = 0; i < 30; i++) {
        Expr* e = parse_expression(
            "FindMinimum[(x^2-1)^2+y^2, {{x,0.1},{y,0.5}}, Method->\"TrustExact\"]");
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
    TEST(test_illcond_n10_c1e6);
    TEST(test_double_well_negcurv);
    TEST(test_beale);
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

    printf("All trust-exact tests passed!\n");
    return 0;
}
