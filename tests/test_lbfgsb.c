/* Unit tests for the L-BFGS-B method of FindMinimum / FindMaximum
 * (Method -> "LBFGSB", src/numerical_calculus/findmin.c).
 *
 * L-BFGS-B is a Mathilda extension to FindMinimum's method set: limited-memory
 * BFGS (O(m*n) per iteration, vs the full-memory QuasiNewton's O(n^2)) with box
 * constraints handled by an active-set projection, plus a strong-Wolfe line
 * search. It is deterministic given a start point (no RNG), so every case fixes
 * its start explicitly. Every pinned target is an ANALYTIC optimum (0, 0.25, 5,
 * n, 3, ...) or a relation between two solves — never a snapshot of the
 * solver's current output — so a future accuracy change cannot silently ratify
 * a wrong answer.
 *
 * Large-n objectives use genuine scalar symbols generated with
 * Symbol["z"<>ToString[i]] (FindMinimum does not auto-expand indexed x[i]
 * specs the way NMinimize does), wrapped in Evaluate[...] to expand past
 * FindMinimum's HoldAll. Hard problems (extended Rosenbrock at large n, high
 * conditioning) are given an explicit MaxIterations budget, as L-BFGS needs
 * many cheap iterations to resolve disparate scales (scipy's L-BFGS-B likewise
 * defaults to 15000).
 *
 * Run binary directly: ./lbfgsb_tests */

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
        fprintf(stdout, "FAIL: %s\n  expected: True\n  got:      %s\n",
                input, got);
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
        fprintf(stdout, "FAIL: %s\n  expected: %s\n  got:      %s\n",
                input, expected, got);
        ASSERT_STR_EQ(got, expected);
    }
    free(got);
    expr_free(e);
    expr_free(res);
}

/* ------------------------------------------------------------------ */
/* G1. Unconstrained smooth, scaling in dimension                      */
/* ------------------------------------------------------------------ */

static void test_rosenbrock_2d(void) {
    /* global min 0 at (1,1). */
    check_true("Abs[First[FindMinimum[(1-x)^2+100(y-x^2)^2, {{x,-1.2},{y,1}}, "
               "Method->\"LBFGSB\"]]] < 1.*^-6");
    check_true("With[{r=Last[FindMinimum[(1-x)^2+100(y-x^2)^2, {{x,-1.2},{y,1}}, "
               "Method->\"LBFGSB\"]]}, Abs[(x/.r)-1]+Abs[(y/.r)-1] < 1.*^-4]");
}

static void test_rosenbrock_n10(void) {
    check_true("Module[{v=Table[Symbol[\"z\"<>ToString[i]],{i,1,10}]}, "
               "Abs[First[FindMinimum["
               "Evaluate[Sum[100(v[[j+1]]-v[[j]]^2)^2+(1-v[[j]])^2,{j,1,9}]], "
               "Evaluate[Table[{v[[k]],-1.2},{k,1,10}]], Method->\"LBFGSB\", "
               "MaxIterations->3000]]] < 1.*^-5]");
}

static void test_rosenbrock_n50(void) {
    check_true("Module[{v=Table[Symbol[\"z\"<>ToString[i]],{i,1,50}]}, "
               "First[FindMinimum["
               "Evaluate[Sum[100(v[[j+1]]-v[[j]]^2)^2+(1-v[[j]])^2,{j,1,49}]], "
               "Evaluate[Table[{v[[k]],-1.2},{k,1,50}]], Method->\"LBFGSB\", "
               "MaxIterations->5000]] < 1.*^-4]");
}

static void test_illcond_quadratic_n200_scaling(void) {
    /* Unimodal ill-conditioned (cond 1e4) quadratic at n=200 → 0. Completing at
     * this size within the budget and the alarm is portable evidence of
     * scaling: a full dense n*n Hessian would be 320 KB and O(n^2) per step.
     * (Extended Rosenbrock is deliberately NOT used at large n — its sum form
     * is multimodal for n>=4, so the reachable optimum is start-dependent and a
     * poor cross-method reference.) */
    check_true("Module[{v=Table[Symbol[\"z\"<>ToString[i]],{i,1,200}]}, "
               "First[FindMinimum[Evaluate[Sum[10^(4 (i-1)/199) v[[i]]^2,{i,1,200}]], "
               "Evaluate[Table[{v[[k]],1.0},{k,1,200}]], Method->\"LBFGSB\", "
               "MaxIterations->3000]] < 1.*^-6]");
}

static void test_illcond_quadratic_n100(void) {
    /* diagonal quadratic, condition 1e4, min 0 at origin. */
    check_true("Module[{v=Table[Symbol[\"z\"<>ToString[i]],{i,1,100}]}, "
               "First[FindMinimum[Evaluate[Sum[10^(4 (i-1)/99) v[[i]]^2,{i,1,100}]], "
               "Evaluate[Table[{v[[k]],1.0},{k,1,100}]], Method->\"LBFGSB\", "
               "MaxIterations->2000]] < 1.*^-6]");
}

static void test_wellcond_quadratic_n1000_scaling(void) {
    /* n=1000 separable quadratic → 0. Large-n completion is the scaling proof;
     * a full-memory method would allocate 8 MB for the Hessian at this size. */
    check_true("Module[{v=Table[Symbol[\"z\"<>ToString[i]],{i,1,1000}]}, "
               "First[FindMinimum[Evaluate[Sum[v[[i]]^2,{i,1,1000}]], "
               "Evaluate[Table[{v[[k]],1.0},{k,1,1000}]], Method->\"LBFGSB\"]] < 1.*^-8]");
}

static void test_beale(void) {
    /* Beale, min 0 at (3, 0.5). */
    check_true("With[{r=Last[FindMinimum["
               "(1.5-x+x y)^2+(2.25-x+x y^2)^2+(2.625-x+x y^3)^2, {{x,1},{y,1}}, "
               "Method->\"LBFGSB\"]]}, Abs[(x/.r)-3.0]+Abs[(y/.r)-0.5] < 1.*^-3]");
}

static void test_wood(void) {
    /* Wood/Colville 4D, min 0 at (1,1,1,1). */
    check_true("Abs[First[FindMinimum["
               "100(y-x^2)^2+(1-x)^2+90(w-z^2)^2+(1-z)^2+10.1((y-1)^2+(w-1)^2)"
               "+19.8(y-1)(w-1), {{x,-3},{y,-1},{z,-3},{w,-1}}, Method->\"LBFGSB\", "
               "MaxIterations->2000]]] < 1.*^-4");
}

static void test_powell_singular(void) {
    /* Powell singular 4D: min 0 at origin, Hessian singular there. */
    check_true("First[FindMinimum["
               "(x+10 y)^2+5(z-w)^2+(y-2 z)^4+10(x-w)^4, "
               "{{x,3},{y,-1},{z,0},{w,1}}, Method->\"LBFGSB\", MaxIterations->2000]] < 1.*^-6");
}

static void test_booth(void) {
    /* Booth, min 0 at (1,3). */
    check_true("With[{r=Last[FindMinimum[(x+2 y-7)^2+(2 x+y-5)^2, {{x,0},{y,0}}, "
               "Method->\"LBFGSB\"]]}, Abs[(x/.r)-1]+Abs[(y/.r)-3] < 1.*^-5]");
}

/* ------------------------------------------------------------------ */
/* G2. Bound constraints active at the solution (the "-B" path)        */
/* ------------------------------------------------------------------ */

static void test_bound_rosenbrock_on_face(void) {
    /* x <= 0.5 cuts off the unconstrained min (1,1); the constrained optimum
     * sits ON the bound at (0.5, 0.25), f = 0.25. This is the case a plain
     * projected-gradient L-BFGS gets wrong. */
    check_true("Abs[First[FindMinimum[{(1-x)^2+100(y-x^2)^2, x<=0.5}, "
               "{{x,0},{y,0}}, Method->\"LBFGSB\"]] - 0.25] < 1.*^-3");
    check_true("(x <= 0.5000001) /. Last[FindMinimum[{(1-x)^2+100(y-x^2)^2, x<=0.5}, "
               "{{x,0},{y,0}}, Method->\"LBFGSB\"]]");
}

static void test_bound_corner(void) {
    /* Unconstrained min (2,3) lies outside [0,1]^2 → constrained optimum is the
     * corner (1,1), f = 5, both coordinates on the upper bound. */
    check_true("Abs[First[FindMinimum[(x-2)^2+(y-3)^2, {{x,0,0,1},{y,0,0,1}}, "
               "Method->\"LBFGSB\"]] - 5.0] < 1.*^-4");
    check_true("With[{r=Last[FindMinimum[(x-2)^2+(y-3)^2, {{x,0,0,1},{y,0,0,1}}, "
               "Method->\"LBFGSB\"]]}, Abs[(x/.r)-1]+Abs[(y/.r)-1] < 1.*^-4]");
}

static void test_bound_many_active(void) {
    /* Sum (x_i-2)^2 on [0,1]^10 → every x_i = 1, f = 10. Assert the value AND
     * that every coordinate respects its box. */
    check_true("Module[{v=Table[Symbol[\"z\"<>ToString[i]],{i,1,10}]}, "
               "Abs[First[FindMinimum[Evaluate[Sum[(v[[i]]-2)^2,{i,1,10}]], "
               "Evaluate[Table[{v[[k]],0,0,1},{k,1,10}]], Method->\"LBFGSB\"]] - 10.0] < 1.*^-4]");
    check_true("Module[{v=Table[Symbol[\"z\"<>ToString[i]],{i,1,10}]}, "
               "And @@ ((0 <= # <= 1.000001) & /@ (v /. Last[FindMinimum["
               "Evaluate[Sum[(v[[i]]-2)^2,{i,1,10}]], "
               "Evaluate[Table[{v[[k]],0,0,1},{k,1,10}]], Method->\"LBFGSB\"]]))]");
}

/* ------------------------------------------------------------------ */
/* G3. Correctness / robustness                                        */
/* ------------------------------------------------------------------ */

static void test_start_at_optimum(void) {
    /* Starting at the optimum converges immediately (no descent step). */
    check_true("Abs[First[FindMinimum[x^2+y^2+z^2, {{x,0},{y,0},{z,0}}, "
               "Method->\"LBFGSB\"]]] < 1.*^-12");
    check_true("(lbCnt=0; FindMinimum[x^2+y^2+z^2, {{x,0},{y,0},{z,0}}, "
               "Method->\"LBFGSB\", StepMonitor:>(lbCnt=lbCnt+1)]; lbCnt <= 1)");
}

static void test_start_on_bound(void) {
    /* min (x-0.5)^2 s.t. x >= 1 → optimum x = 1 (on the bound), f = 0.25. */
    check_true("Abs[First[FindMinimum[{(x-0.5)^2, x>=1}, {x,1}, Method->\"LBFGSB\"]] - 0.25] < 1.*^-6");
}

static void test_many_curvature_updates(void) {
    /* cond 1e4, n=60: the memory ring must cycle correctly to still reach 0. */
    check_true("Module[{v=Table[Symbol[\"z\"<>ToString[i]],{i,1,60}]}, "
               "First[FindMinimum[Evaluate[Sum[10^(4 (i-1)/59) v[[i]]^2,{i,1,60}]], "
               "Evaluate[Table[{v[[k]],1.0},{k,1,60}]], Method->\"LBFGSB\", "
               "MaxIterations->2000]] < 1.*^-6]");
}

static void test_nonconvex_separable(void) {
    /* Sum (x_i^4 - 3 x_i^2), min -2.25 per coord at x_i = Sqrt[3/2]; from
     * x_i = 1 it descends into the positive well. n=5 → 5*(-2.25). */
    check_true("Module[{v=Table[Symbol[\"z\"<>ToString[i]],{i,1,5}]}, "
               "Abs[First[FindMinimum[Evaluate[Sum[v[[i]]^4-3 v[[i]]^2,{i,1,5}]], "
               "Evaluate[Table[{v[[k]],1.0},{k,1,5}]], Method->\"LBFGSB\"]] - (5 (-2.25))] < 1.*^-3]");
}

static void test_gradient_supplied_matches_fd(void) {
    /* User Gradient and the finite-difference fallback reach the same optimum. */
    check_true("Abs[First[FindMinimum["
               "(1.5-x+x y)^2+(2.25-x+x y^2)^2+(2.625-x+x y^3)^2, {{x,1},{y,1}}, "
               "Method->\"LBFGSB\", Gradient->{"
               "2(1.5-x+x y)(y-1)+2(2.25-x+x y^2)(y^2-1)+2(2.625-x+x y^3)(y^3-1), "
               "2(1.5-x+x y)x+2(2.25-x+x y^2)(2 x y)+2(2.625-x+x y^3)(3 x y^2)}]]] < 1.*^-6");
    check_true("Abs[First[FindMinimum["
               "(1.5-x+x y)^2+(2.25-x+x y^2)^2+(2.625-x+x y^3)^2, {{x,1},{y,1}}, "
               "Method->\"LBFGSB\"]]] < 1.*^-5");
}

/* ------------------------------------------------------------------ */
/* G4. Parity with the existing full-memory method                     */
/* ------------------------------------------------------------------ */

static void test_parity_rosenbrock(void) {
    /* LBFGSB and QuasiNewton both reach the analytic optimum 0. */
    check_true("Module[{v=Table[Symbol[\"z\"<>ToString[i]],{i,1,10}]}, "
               "First[FindMinimum[Evaluate[Sum[100(v[[j+1]]-v[[j]]^2)^2+(1-v[[j]])^2,{j,1,9}]], "
               "Evaluate[Table[{v[[k]],-1.2},{k,1,10}]], Method->\"LBFGSB\", MaxIterations->3000]] < 1.*^-5]");
    check_true("Module[{v=Table[Symbol[\"z\"<>ToString[i]],{i,1,10}]}, "
               "First[FindMinimum[Evaluate[Sum[100(v[[j+1]]-v[[j]]^2)^2+(1-v[[j]])^2,{j,1,9}]], "
               "Evaluate[Table[{v[[k]],-1.2},{k,1,10}]], Method->\"QuasiNewton\", MaxIterations->3000]] < 1.*^-5]");
}

static void test_parity_illcond(void) {
    check_true("Module[{v=Table[Symbol[\"z\"<>ToString[i]],{i,1,30}]}, "
               "First[FindMinimum[Evaluate[Sum[10^(4 (i-1)/29) v[[i]]^2,{i,1,30}]], "
               "Evaluate[Table[{v[[k]],1.0},{k,1,30}]], Method->\"LBFGSB\", MaxIterations->1000]] < 1.*^-6]");
    check_true("Module[{v=Table[Symbol[\"z\"<>ToString[i]],{i,1,30}]}, "
               "First[FindMinimum[Evaluate[Sum[10^(4 (i-1)/29) v[[i]]^2,{i,1,30}]], "
               "Evaluate[Table[{v[[k]],1.0},{k,1,30}]], Method->\"QuasiNewton\", MaxIterations->1000]] < 1.*^-6]");
}

/* ------------------------------------------------------------------ */
/* G5. General-constraint interplay (L-BFGS-B as the inner AL solver)  */
/* ------------------------------------------------------------------ */

static void test_general_equality_box(void) {
    /* min x^2+y^2+z^2 s.t. x+y+z==3, 0<=x,y,z<=2 → KKT optimum (1,1,1), f=3. */
    check_true("Abs[First[FindMinimum[{x^2+y^2+z^2, x+y+z==3 && 0<=x<=2 && 0<=y<=2 && 0<=z<=2}, "
               "{{x,0},{y,0},{z,0}}, Method->\"LBFGSB\"]] - 3.0] < 1.*^-3");
    check_true("Abs[(x+y+z /. Last[FindMinimum[{x^2+y^2+z^2, x+y+z==3 && 0<=x<=2 && 0<=y<=2 && 0<=z<=2}, "
               "{{x,0},{y,0},{z,0}}, Method->\"LBFGSB\"]]) - 3.0] < 1.*^-3");
}

static void test_general_equality_binds_bound(void) {
    /* min (x-2)^2+(y-2)^2 on the unit circle with x,y>=0 → (1/Sqrt2, 1/Sqrt2). */
    check_true("Abs[First[FindMinimum[{(x-2)^2+(y-2)^2, x^2+y^2==1 && x>=0 && y>=0}, "
               "{{x,0.6},{y,0.6}}, Method->\"LBFGSB\"]] - 2 (1/Sqrt[2]-2)^2] < 1.*^-2");
}

/* ------------------------------------------------------------------ */
/* G6. Shape / diagnostics                                             */
/* ------------------------------------------------------------------ */

static void test_shape(void) {
    check_eq("Head[FindMinimum[x^2+y^2, {{x,1},{y,1}}, Method->\"LBFGSB\"]]", "List");
    check_eq("Length[FindMinimum[x^2+y^2, {{x,1},{y,1}}, Method->\"LBFGSB\"]]", "2");
    /* A too-small MaxIterations still returns a well-formed result. */
    check_eq("Head[FindMinimum[(1-x)^2+100(y-x^2)^2, {{x,0},{y,0}}, "
             "Method->\"LBFGSB\", MaxIterations->3]]", "List");
}

static void test_aliases(void) {
    /* "LBFGS" and "LimitedMemoryBFGS" select the same method. */
    check_true("Abs[First[FindMinimum[(1-x)^2+100(y-x^2)^2, {{x,-1.2},{y,1}}, "
               "Method->\"LBFGS\"]]] < 1.*^-6");
    check_true("Abs[First[FindMinimum[(1-x)^2+100(y-x^2)^2, {{x,-1.2},{y,1}}, "
               "Method->\"LimitedMemoryBFGS\"]]] < 1.*^-6");
}

static void test_findmaximum(void) {
    /* FindMaximum of an inverted paraboloid with a box → (1,1), value 0. */
    check_true("Abs[First[FindMaximum[-(x-1)^2-(y-1)^2, {{x,0,0,2},{y,0,0,2}}, "
               "Method->\"LBFGSB\"]]] < 1.*^-8");
}

/* ------------------------------------------------------------------ */
/* G7. Memory hygiene smoke                                            */
/* ------------------------------------------------------------------ */

static void test_no_leak_many_calls(void) {
    for (int i = 0; i < 30; i++) {
        Expr* e = parse_expression(
            "FindMinimum[{(x-2)^2+(y-3)^2, x+y<=1}, {{x,0,-5,5},{y,0,-5,5}}, Method->\"LBFGSB\"]");
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
    /* Numeric routines emit expected diagnostics (lstol at a hard bound, etc.);
     * keep them off stdout so only FAIL lines show. */
    freopen("/dev/null", "w", stderr);

    /* G1 unconstrained + scaling */
    TEST(test_rosenbrock_2d);
    TEST(test_rosenbrock_n10);
    TEST(test_rosenbrock_n50);
    TEST(test_illcond_quadratic_n200_scaling);
    TEST(test_illcond_quadratic_n100);
    TEST(test_wellcond_quadratic_n1000_scaling);
    TEST(test_beale);
    TEST(test_wood);
    TEST(test_powell_singular);
    TEST(test_booth);

    /* G2 bounds active */
    TEST(test_bound_rosenbrock_on_face);
    TEST(test_bound_corner);
    TEST(test_bound_many_active);

    /* G3 robustness */
    TEST(test_start_at_optimum);
    TEST(test_start_on_bound);
    TEST(test_many_curvature_updates);
    TEST(test_nonconvex_separable);
    TEST(test_gradient_supplied_matches_fd);

    /* G4 parity */
    TEST(test_parity_rosenbrock);
    TEST(test_parity_illcond);

    /* G5 general constraints */
    TEST(test_general_equality_box);
    TEST(test_general_equality_binds_bound);

    /* G6 shape / diagnostics */
    TEST(test_shape);
    TEST(test_aliases);
    TEST(test_findmaximum);

    /* G7 memory */
    TEST(test_no_leak_many_calls);

    printf("All L-BFGS-B tests passed!\n");
    return 0;
}
