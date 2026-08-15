/* Cross-method extensive test battery for FindMinimum / FindMaximum
 * (src/numerical_calculus/findmin.c).
 *
 * The per-method files (test_findmin.c, test_lbfgsb.c, test_findmin_powell.c,
 * test_findmin_neldermead.c) each exercise one method in depth. This file does
 * the complementary sweep: it runs a SHARED battery of standard test functions
 * with known analytic optima across EVERY applicable method, so a regression in
 * any one method on any one function surfaces here. Each n-D case is solved by
 * all six n-D methods
 *     QuasiNewton, ConjugateGradient, Newton, LBFGSB, Powell, NelderMead
 * and every 1-D case by Brent plus the two derivative-free methods that delegate
 * n==1 to Brent. Failures are labelled with the offending method.
 *
 * All targets are analytic optima (0, 5, ...) or the exact optimiser location,
 * never a snapshot of the current solver output. MaxIterations is set generously
 * so the slower derivative-free methods also converge; the gradient methods stop
 * earlier of their own accord.
 *
 * Run binary directly: ./findmin_methods_tests */

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

/* Every method that operates on n>=2 variables. */
static const char* ND_METHODS[] = {
    "QuasiNewton", "ConjugateGradient", "Newton", "LBFGSB", "Powell", "NelderMead", "TNC"
};
#define N_ND_METHODS ((int)(sizeof(ND_METHODS) / sizeof(ND_METHODS[0])))

/* Evaluate `input`, require the fullform to equal "True", label failures. */
static void check_true_labeled(const char* input, const char* label) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* got = expr_to_string_fullform(res);
    if (strcmp(got, "True") != 0) {
        fprintf(stdout, "FAIL [%s]: %s\n  expected: True\n  got:      %s\n",
                label, input, got);
        ASSERT_STR_EQ(got, "True");
    }
    free(got);
    expr_free(e);
    expr_free(res);
}

/* Solve `obj` over `spec` with EVERY n-D method and assert `assertion` (which
 * refers to the result as `r`). */
static void check_all_nd(const char* obj, const char* spec, const char* assertion) {
    for (int i = 0; i < N_ND_METHODS; i++) {
        char buf[1500];
        snprintf(buf, sizeof buf,
                 "With[{r=FindMinimum[%s, %s, Method->\"%s\", MaxIterations->5000]}, %s]",
                 obj, spec, ND_METHODS[i], assertion);
        check_true_labeled(buf, ND_METHODS[i]);
    }
}

/* ------------------------------------------------------------------ */
/* M1. Smooth unconstrained battery — every n-D method                 */
/* ------------------------------------------------------------------ */

static void test_quadratic_bowl(void) {
    /* min 0 at (1,2). */
    check_all_nd("(x-1)^2+(y-2)^2", "{{x,0},{y,0}}",
                 "Abs[First[r]] < 1.*^-6 && Abs[(x/.Last[r])-1]+Abs[(y/.Last[r])-2] < 1.*^-3");
}

static void test_rosenbrock(void) {
    /* curved valley; min 0 at (1,1). Point tol is loose (5e-3) because the
     * derivative-free methods stop with f~1e-8 while still ~1e-3 from (1,1). */
    check_all_nd("(1-x)^2+100(y-x^2)^2", "{{x,-1.2},{y,1}}",
                 "Abs[First[r]] < 1.*^-5 && Abs[(x/.Last[r])-1]+Abs[(y/.Last[r])-1] < 5.*^-3");
}

static void test_booth(void) {
    /* min 0 at (1,3). */
    check_all_nd("(x+2 y-7)^2+(2 x+y-5)^2", "{{x,0},{y,0}}",
                 "Abs[First[r]] < 1.*^-6 && Abs[(x/.Last[r])-1]+Abs[(y/.Last[r])-3] < 1.*^-3");
}

static void test_himmelblau(void) {
    /* four equal minima at f=0; from (1,1) every method descends to (3,2). */
    check_all_nd("(x^2+y-11)^2+(x+y^2-7)^2", "{{x,1},{y,1}}",
                 "Abs[First[r]] < 1.*^-5 && Abs[(x/.Last[r])-3]+Abs[(y/.Last[r])-2] < 1.*^-3");
}

static void test_sphere_n5(void) {
    /* separable quadratic n=5, min 0 at origin, over every n-D method. */
    for (int i = 0; i < N_ND_METHODS; i++) {
        char buf[1500];
        snprintf(buf, sizeof buf,
                 "Module[{v=Table[Symbol[\"z\"<>ToString[i]],{i,1,5}]}, "
                 "Abs[First[FindMinimum[Evaluate[Total[Table[v[[i]]^2,{i,1,5}]]], "
                 "Evaluate[Table[{v[[k]],1.0},{k,1,5}]], Method->\"%s\", "
                 "MaxIterations->5000]]] < 1.*^-7]", ND_METHODS[i]);
        check_true_labeled(buf, ND_METHODS[i]);
    }
}

/* ------------------------------------------------------------------ */
/* M2. Box bounds — every n-D method                                   */
/* ------------------------------------------------------------------ */

static void test_bound_corner(void) {
    /* unconstrained min (2,3) is outside [0,1]^2; optimum on the corner (1,1),
     * value 5. Every method honours the box (projection / active set / clamp). */
    check_all_nd("{(x-2)^2+(y-3)^2, 0<=x<=1 && 0<=y<=1}", "{{x,0.5},{y,0.5}}",
                 "Abs[First[r]-5.0] < 1.*^-4 && Abs[(x/.Last[r])-1]+Abs[(y/.Last[r])-1] < 1.*^-3");
}

/* ------------------------------------------------------------------ */
/* M3. FindMaximum — every n-D method                                  */
/* ------------------------------------------------------------------ */

static void test_findmaximum_all(void) {
    /* concave paraboloid, max 0 at (1,2). */
    for (int i = 0; i < N_ND_METHODS; i++) {
        char buf[1500];
        snprintf(buf, sizeof buf,
                 "With[{r=FindMaximum[-(x-1)^2-(y-2)^2, {{x,0},{y,0}}, "
                 "Method->\"%s\", MaxIterations->5000]}, "
                 "Abs[First[r]] < 1.*^-6 && Abs[(x/.Last[r])-1]+Abs[(y/.Last[r])-2] < 1.*^-3]",
                 ND_METHODS[i]);
        check_true_labeled(buf, ND_METHODS[i]);
    }
}

/* ------------------------------------------------------------------ */
/* M4. Result shape — every n-D method                                 */
/* ------------------------------------------------------------------ */

static void test_shape_all(void) {
    for (int i = 0; i < N_ND_METHODS; i++) {
        char buf[900];
        snprintf(buf, sizeof buf,
                 "Length[FindMinimum[(x-1)^2+(y-2)^2, {{x,0},{y,0}}, Method->\"%s\"]] == 2 && "
                 "Length[Last[FindMinimum[(x-1)^2+(y-2)^2, {{x,0},{y,0}}, Method->\"%s\"]]] == 2",
                 ND_METHODS[i], ND_METHODS[i]);
        check_true_labeled(buf, ND_METHODS[i]);
    }
}

/* ------------------------------------------------------------------ */
/* M5. 1-D methods                                                     */
/* ------------------------------------------------------------------ */

static void test_1d_methods(void) {
    /* (x-3)^2 min at 3. Brent is the native 1-D method; Powell and NelderMead
     * delegate n==1 to Brent, so all three agree. */
    const char* m1d[] = {"Brent", "Powell", "NelderMead"};
    for (int i = 0; i < 3; i++) {
        char buf[600];
        snprintf(buf, sizeof buf,
                 "Abs[(x /. Last[FindMinimum[(x-3)^2, {x, 0}, Method->\"%s\"]]) - 3.0] < 1.*^-5",
                 m1d[i]);
        check_true_labeled(buf, m1d[i]);
    }
    /* A transcendental 1-D min: x Cos[x] near x=3 → ~3.4256. */
    const char* m1d2[] = {"Brent", "Powell", "NelderMead"};
    for (int i = 0; i < 3; i++) {
        char buf[600];
        snprintf(buf, sizeof buf,
                 "Abs[First[FindMinimum[x Cos[x], {x, 3}, Method->\"%s\"]] - (-3.28837)] < 1.*^-3",
                 m1d2[i]);
        check_true_labeled(buf, m1d2[i]);
    }
}

/* ------------------------------------------------------------------ */
/* M6. Automatic method selection                                      */
/* ------------------------------------------------------------------ */

static void test_automatic(void) {
    /* Automatic: Brent for n==1, QuasiNewton for n>=2. Both reach the optimum. */
    check_true_labeled("Abs[(x /. Last[FindMinimum[(x-3)^2, {x, 0}]]) - 3.0] < 1.*^-5",
                       "Automatic-1D");
    check_true_labeled("With[{r=Last[FindMinimum[(x-1)^2+(y-2)^2, {{x,0},{y,0}}]]}, "
                       "Abs[(x/.r)-1]+Abs[(y/.r)-2] < 1.*^-4]", "Automatic-nD");
}

/* ------------------------------------------------------------------ */

int main(void) {
    symtab_init();
    core_init();
    freopen("/dev/null", "w", stderr);

    TEST(test_quadratic_bowl);
    TEST(test_rosenbrock);
    TEST(test_booth);
    TEST(test_himmelblau);
    TEST(test_sphere_n5);
    TEST(test_bound_corner);
    TEST(test_findmaximum_all);
    TEST(test_shape_all);
    TEST(test_1d_methods);
    TEST(test_automatic);

    printf("All FindMinimum cross-method tests passed!\n");
    return 0;
}
