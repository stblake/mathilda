/* Unit tests for the Nelder-Mead method of FindMinimum / FindMaximum
 * (Method -> "NelderMead"; src/numerical_calculus/findmin.c).
 *
 * Nelder-Mead is FindMinimum's downhill-simplex derivative-free method (the
 * LOCAL simplex, matching scipy's _minimize_neldermead — distinct from the
 * restarted, box-sampling simplex the global NMinimize driver uses). It is
 * deterministic given a start point (no RNG), so every case fixes its start
 * explicitly. Every pinned target is an ANALYTIC optimum (0, 18, -50, ...) or a
 * relation between two solves — never a snapshot of the solver's current output.
 *
 * Nelder-Mead is strong on SMOOTH derivative-free problems; it is famously weak
 * on NON-SMOOTH ones (McKinnon 1998 — the simplex can converge to a
 * non-stationary point). The non-smooth L1 case that Powell solves is therefore
 * NOT asserted here; use Method -> "Powell" for non-smooth objectives.
 *
 * Large objectives use genuine scalar symbols generated with
 * Symbol["z"<>ToString[i]], wrapped in Evaluate[...] past FindMinimum's HoldAll.
 *
 * Run binary directly: ./neldermead_tests */

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
/* N1. Unconstrained smooth, standard test functions                   */
/* ------------------------------------------------------------------ */

static void test_rosenbrock_2d(void) {
    /* curved valley; global min 0 at (1,1). */
    check_true("Abs[First[FindMinimum[(1-x)^2+100(y-x^2)^2, {{x,-1.2},{y,1}}, "
               "Method->\"NelderMead\", MaxIterations->2000]]] < 1.*^-6");
    check_true("With[{r=Last[FindMinimum[(1-x)^2+100(y-x^2)^2, {{x,-1.2},{y,1}}, "
               "Method->\"NelderMead\", MaxIterations->2000]]}, "
               "Abs[(x/.r)-1]+Abs[(y/.r)-1] < 1.*^-3]");
}

static void test_beale(void) {
    /* Beale, min 0 at (3, 0.5). */
    check_true("With[{r=Last[FindMinimum["
               "(1.5-x+x y)^2+(2.25-x+x y^2)^2+(2.625-x+x y^3)^2, {{x,1},{y,1}}, "
               "Method->\"NelderMead\", MaxIterations->2000]]}, "
               "Abs[(x/.r)-3.0]+Abs[(y/.r)-0.5] < 1.*^-3]");
}

static void test_booth(void) {
    /* Booth, min 0 at (1,3). */
    check_true("With[{r=Last[FindMinimum[(x+2 y-7)^2+(2 x+y-5)^2, {{x,0},{y,0}}, "
               "Method->\"NelderMead\", MaxIterations->2000]]}, "
               "Abs[(x/.r)-1]+Abs[(y/.r)-3] < 1.*^-4]");
}

static void test_matyas(void) {
    /* Matyas, min 0 at (0,0); nearly singular Hessian. */
    check_true("Abs[First[FindMinimum[0.26(x^2+y^2)-0.48 x y, {{x,3},{y,3}}, "
               "Method->\"NelderMead\", MaxIterations->2000]]] < 1.*^-8");
}

static void test_powell_singular(void) {
    /* Powell's singular quartic 4D, min 0 at the origin. */
    check_true("Abs[First[FindMinimum["
               "(x+10 y)^2+5(z-w)^2+(y-2 z)^4+10(x-w)^4, "
               "{{x,3},{y,-1},{z,0},{w,1}}, Method->\"NelderMead\", "
               "MaxIterations->4000]]] < 1.*^-6");
}

static void test_wood(void) {
    /* Wood/Colville 4D, min 0 at (1,1,1,1). */
    check_true("Abs[First[FindMinimum["
               "100(y-x^2)^2+(1-x)^2+90(w-z^2)^2+(1-z)^2+10.1((y-1)^2+(w-1)^2)"
               "+19.8(y-1)(w-1), {{x,-3},{y,-1},{z,-3},{w,-1}}, "
               "Method->\"NelderMead\", MaxIterations->4000]]] < 1.*^-6");
}

static void test_himmelblau(void) {
    /* Himmelblau: four equal global minima at f=0; from (1,1) the simplex
     * descends into the (3,2) basin. */
    check_true("With[{r=FindMinimum[(x^2+y-11)^2+(x+y^2-7)^2, {{x,1},{y,1}}, "
               "Method->\"NelderMead\", MaxIterations->3000]}, "
               "Abs[First[r]] < 1.*^-6 && Abs[(x/.Last[r])-3]+Abs[(y/.Last[r])-2] < 1.*^-3]");
}

static void test_sphere_n5(void) {
    /* separable quadratic n=5, min 0 at origin. */
    check_true("Module[{v=Table[Symbol[\"z\"<>ToString[i]],{i,1,5}]}, "
               "Abs[First[FindMinimum[Evaluate[Total[Table[v[[i]]^2,{i,1,5}]]], "
               "Evaluate[Table[{v[[k]],1.0},{k,1,5}]], Method->\"NelderMead\", "
               "MaxIterations->3000]]] < 1.*^-8]");
}

/* ------------------------------------------------------------------ */
/* N2. Box bounds                                                       */
/* ------------------------------------------------------------------ */

static void test_bound_corner(void) {
    /* Unconstrained min (5,5) is outside x<=2, y<=2; the constrained optimum
     * sits in the box corner (2,2), value 18. */
    check_true("With[{r=FindMinimum[{(x-5)^2+(y-5)^2, x<=2 && y<=2}, "
               "{{x,0},{y,0}}, Method->\"NelderMead\", MaxIterations->2000]}, "
               "Abs[First[r]-18.0] < 1.*^-5 && Abs[(x/.Last[r])-2]+Abs[(y/.Last[r])-2] < 1.*^-4]");
}

/* ------------------------------------------------------------------ */
/* N3. n == 1 delegates to the exact Brent path                        */
/* ------------------------------------------------------------------ */

static void test_n1_matches_brent(void) {
    /* A 1-variable Nelder-Mead solve delegates to the same Brent 1-D minimiser
     * as Method -> "Brent" (a 2-vertex simplex is a crude line search), so the
     * located minimiser must agree. Sin[x] from x=2 → 3 Pi/2 ~ 4.712. */
    check_true("With[{a=x/.Last[FindMinimum[Sin[x], {x,2}, Method->\"NelderMead\"]], "
               "b=x/.Last[FindMinimum[Sin[x], {x,2}, Method->\"Brent\"]]}, "
               "Abs[a-b] < 1.*^-6]");
}

/* ------------------------------------------------------------------ */
/* N4. FindMaximum, shape, diagnostics                                 */
/* ------------------------------------------------------------------ */

static void test_findmaximum(void) {
    /* FindMaximum minimises -f; concave paraboloid, max 0 at (3,-1). */
    check_true("With[{r=FindMaximum[-(x-3)^2-(y+1)^2, {{x,0},{y,0}}, "
               "Method->\"NelderMead\", MaxIterations->2000]}, "
               "Abs[First[r]] < 1.*^-6 && Abs[(x/.Last[r])-3]+Abs[(y/.Last[r])+1] < 1.*^-3]");
}

static void test_shape(void) {
    check_eq("Length[FindMinimum[(x-1)^2+(y-2)^2, {{x,0},{y,0}}, Method->\"NelderMead\"]]", "2");
    check_eq("Length[Last[FindMinimum[(x-1)^2+(y-2)^2, {{x,0},{y,0}}, Method->\"NelderMead\"]]]", "2");
}

static void test_general_constraint_rejected(void) {
    /* Like scipy's Nelder-Mead, general (non-box) constraints are unsupported —
     * FindMinimum::nimpl, and the call returns unevaluated. */
    check_eq("Head[FindMinimum[{x^2+y^2, x^2+y^2>=1}, {{x,2},{y,2}}, Method->\"NelderMead\"]]",
             "FindMinimum");
}

/* ------------------------------------------------------------------ */
/* N5. Memory hygiene smoke                                            */
/* ------------------------------------------------------------------ */

static void test_no_leak_many_calls(void) {
    for (int i = 0; i < 30; i++) {
        Expr* e = parse_expression(
            "FindMinimum[(1-x)^2+100(y-x^2)^2, {{x,-1.2},{y,1}}, Method->\"NelderMead\"]");
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
    /* Numeric routines emit expected diagnostics (nimpl for the rejected
     * general constraint, etc.); keep them off stdout so only FAIL lines show. */
    freopen("/dev/null", "w", stderr);

    /* N1 unconstrained smooth */
    TEST(test_rosenbrock_2d);
    TEST(test_beale);
    TEST(test_booth);
    TEST(test_matyas);
    TEST(test_powell_singular);
    TEST(test_wood);
    TEST(test_himmelblau);
    TEST(test_sphere_n5);

    /* N2 box bounds */
    TEST(test_bound_corner);

    /* N3 n == 1 parity */
    TEST(test_n1_matches_brent);

    /* N4 maximum / shape / diagnostics */
    TEST(test_findmaximum);
    TEST(test_shape);
    TEST(test_general_constraint_rejected);

    /* N5 memory */
    TEST(test_no_leak_many_calls);

    printf("All Nelder-Mead tests passed!\n");
    return 0;
}
