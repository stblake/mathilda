/* Unit tests for the Powell method of FindMinimum / FindMaximum
 * (Method -> "Powell", alias "PrincipalAxis"; src/numerical_calculus/findmin.c).
 *
 * Powell is FindMinimum's derivative-free local method: a conjugate-direction
 * search that minimises using function values only (no gradient, analytic or
 * finite-difference), with the same Brent parabolic-interpolation line search
 * as Method -> "Brent" restricted to phi(t) = f(p + t d). It is deterministic
 * given a start point (no RNG), so every case fixes its start explicitly. Every
 * pinned target is an ANALYTIC optimum (0, 18, -50, (3,0.5), ...) or a relation
 * between two solves — never a snapshot of the solver's current output — so a
 * future accuracy change cannot silently ratify a wrong answer.
 *
 * Large objectives use genuine scalar symbols generated with
 * Symbol["z"<>ToString[i]] (FindMinimum does not auto-expand indexed x[i] specs
 * the way NMinimize does), wrapped in Evaluate[...] to expand past FindMinimum's
 * HoldAll.
 *
 * Run binary directly: ./powell_tests */

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
/* P1. Unconstrained smooth, standard test functions                   */
/* ------------------------------------------------------------------ */

static void test_rosenbrock_2d(void) {
    /* curved valley; global min 0 at (1,1). */
    check_true("Abs[First[FindMinimum[(1-x)^2+100(y-x^2)^2, {{x,-1.2},{y,1}}, "
               "Method->\"Powell\"]]] < 1.*^-6");
    /* Point tolerance is 1e-3, not 1e-4: on Rosenbrock's flat curved valley a
     * derivative-free method drives f below 1e-9 while the point is still ~1e-4
     * from (1,1) (scipy's Powell lands the same distance out). The objective
     * value above is the tight, meaningful check. */
    check_true("With[{r=Last[FindMinimum[(1-x)^2+100(y-x^2)^2, {{x,-1.2},{y,1}}, "
               "Method->\"Powell\"]]}, Abs[(x/.r)-1]+Abs[(y/.r)-1] < 1.*^-3]");
}

static void test_beale(void) {
    /* Beale, min 0 at (3, 0.5). Regression for the fm_bracket_line strict-`<`
     * growth test: Beale is FLAT in x at y==1 (the start row), and a `<=` test
     * chased x to ~1e37 where floating-point cancellation faked a spurious
     * minimum. Strict `<` stops at the flat region and the search reaches
     * (3, 0.5). */
    check_true("Abs[First[FindMinimum["
               "(1.5-x+x y)^2+(2.25-x+x y^2)^2+(2.625-x+x y^3)^2, {{x,1},{y,1}}, "
               "Method->\"Powell\"]]] < 1.*^-6");
    check_true("With[{r=Last[FindMinimum["
               "(1.5-x+x y)^2+(2.25-x+x y^2)^2+(2.625-x+x y^3)^2, {{x,1},{y,1}}, "
               "Method->\"Powell\"]]}, Abs[(x/.r)-3.0]+Abs[(y/.r)-0.5] < 1.*^-3]");
}

static void test_booth(void) {
    /* Booth, min 0 at (1,3). */
    check_true("With[{r=Last[FindMinimum[(x+2 y-7)^2+(2 x+y-5)^2, {{x,0},{y,0}}, "
               "Method->\"Powell\"]]}, Abs[(x/.r)-1]+Abs[(y/.r)-3] < 1.*^-5]");
}

static void test_matyas(void) {
    /* Matyas, min 0 at (0,0); nearly singular Hessian. */
    check_true("Abs[First[FindMinimum[0.26(x^2+y^2)-0.48 x y, {{x,-3},{y,2}}, "
               "Method->\"Powell\"]]] < 1.*^-8");
}

static void test_powell_singular(void) {
    /* Powell's own singular quartic 4D, min 0 at the origin (Hessian singular
     * there — the namesake hard case for the method). */
    check_true("Abs[First[FindMinimum["
               "(x+10 y)^2+5(z-w)^2+(y-2 z)^4+10(x-w)^4, "
               "{{x,3},{y,-1},{z,0},{w,1}}, Method->\"Powell\", "
               "MaxIterations->2000]]] < 1.*^-6");
}

static void test_wood(void) {
    /* Wood/Colville 4D, min 0 at (1,1,1,1). */
    check_true("Abs[First[FindMinimum["
               "100(y-x^2)^2+(1-x)^2+90(w-z^2)^2+(1-z)^2+10.1((y-1)^2+(w-1)^2)"
               "+19.8(y-1)(w-1), {{x,-3},{y,-1},{z,-3},{w,-1}}, "
               "Method->\"Powell\", MaxIterations->3000]]] < 1.*^-6");
}

static void test_trid_n6(void) {
    /* Trid function n=6: Sum[(x_i-1)^2] - Sum[x_i x_{i-1}]; analytic global
     * min -50. Exercises the conjugate-direction build-up over n=6 coupled
     * variables (dense quadratic, off-diagonal coupling). */
    check_true("Module[{v=Table[Symbol[\"z\"<>ToString[i]],{i,1,6}]}, "
               "Abs[First[FindMinimum["
               "Evaluate[Sum[(v[[i]]-1)^2,{i,1,6}] - Sum[v[[i]] v[[i-1]],{i,2,6}]], "
               "Evaluate[Table[{v[[k]],1.0},{k,1,6}]], Method->\"Powell\", "
               "MaxIterations->2000]] - (-50.0)] < 1.*^-4]");
}

/* ------------------------------------------------------------------ */
/* P2. Derivative-free advantage: non-smooth objective                 */
/* ------------------------------------------------------------------ */

static void test_nonsmooth_l1(void) {
    /* Sum of absolute values, min 0 at (1,2,3). The objective is
     * non-differentiable at the optimum, so the gradient methods stall; Powell
     * (function values only, Brent line search over the piecewise-linear
     * profile) reaches it. This is the case the whole method exists for. */
    check_true("Abs[First[FindMinimum[Abs[x-1]+Abs[y-2]+Abs[w-3], "
               "{{x,0},{y,0},{w,0}}, Method->\"Powell\"]]] < 1.*^-4");
    check_true("With[{r=Last[FindMinimum[Abs[x-1]+Abs[y-2]+Abs[w-3], "
               "{{x,0},{y,0},{w,0}}, Method->\"Powell\"]]}, "
               "Abs[(x/.r)-1]+Abs[(y/.r)-2]+Abs[(w/.r)-3] < 1.*^-4]");
}

/* ------------------------------------------------------------------ */
/* P3. Box bounds                                                       */
/* ------------------------------------------------------------------ */

static void test_bound_corner(void) {
    /* Unconstrained min of (x-5)^2+(y-5)^2 is (5,5), outside x<=2, y<=2; the
     * constrained optimum sits in the box corner (2,2) with value 18. Tests the
     * feasible-t-interval clamp of the line search on two active bounds. */
    check_true("With[{r=FindMinimum[{(x-5)^2+(y-5)^2, x<=2 && y<=2}, "
               "{{x,0},{y,0}}, Method->\"Powell\"]}, "
               "Abs[First[r]-18.0] < 1.*^-6 && Abs[(x/.Last[r])-2]+Abs[(y/.Last[r])-2] < 1.*^-6]");
}

static void test_bound_single(void) {
    /* One active bound: (x-5)^2 over x<=2 → x=2, value 9. n==1 delegates to the
     * Brent path, which stops within its line-search tolerance of the boundary
     * minimum (residual ~2e-7), so the check is 1e-6. */
    check_true("Abs[First[FindMinimum[{(x-5)^2, x<=2}, {x,0}, Method->\"Powell\"]] - 9.0] < 1.*^-6");
}

/* ------------------------------------------------------------------ */
/* P4. n == 1 delegates to the exact Brent path                        */
/* ------------------------------------------------------------------ */

static void test_n1_matches_brent(void) {
    /* A 1-variable Powell solve delegates to the same Brent 1-D minimiser as
     * Method -> "Brent", so the located minimiser must agree. Sin[x] from x=2
     * → 3 Pi/2 ~ 4.712. */
    check_true("With[{a=x/.Last[FindMinimum[Sin[x], {x,2}, Method->\"Powell\"]], "
               "b=x/.Last[FindMinimum[Sin[x], {x,2}, Method->\"Brent\"]]}, "
               "Abs[a-b] < 1.*^-6]");
}

/* ------------------------------------------------------------------ */
/* P5. Alias, FindMaximum, shape                                       */
/* ------------------------------------------------------------------ */

static void test_principalaxis_alias(void) {
    /* "PrincipalAxis" (Mathematica's name) routes to the identical solver, so
     * the deterministic result is bit-for-bit the same as "Powell" — and it no
     * longer errors FindMinimum::nimpl the way it did before this method landed. */
    check_true("First[FindMinimum[(1-x)^2+100(y-x^2)^2, {{x,-1.2},{y,1}}, "
               "Method->\"PrincipalAxis\"]] == "
               "First[FindMinimum[(1-x)^2+100(y-x^2)^2, {{x,-1.2},{y,1}}, "
               "Method->\"Powell\"]]");
}

static void test_findmaximum(void) {
    /* FindMaximum minimises -f; concave paraboloid, max 0 at (3,-1). */
    check_true("With[{r=FindMaximum[-(x-3)^2-(y+1)^2, {{x,0},{y,0}}, "
               "Method->\"Powell\"]}, "
               "Abs[First[r]] < 1.*^-6 && Abs[(x/.Last[r])-3]+Abs[(y/.Last[r])+1] < 1.*^-4]");
}

static void test_shape(void) {
    /* Result is {fmin, {rules...}}: a 2-element list. */
    check_eq("Length[FindMinimum[(x-1)^2+(y-2)^2, {{x,0},{y,0}}, Method->\"Powell\"]]", "2");
    check_eq("Length[Last[FindMinimum[(x-1)^2+(y-2)^2, {{x,0},{y,0}}, Method->\"Powell\"]]]", "2");
}

static void test_general_constraint_rejected(void) {
    /* scipy's Powell has no constrained form; Mathilda mirrors that — a general
     * (non-box) constraint emits FindMinimum::nimpl and the call returns
     * unevaluated (head stays FindMinimum). */
    check_eq("Head[FindMinimum[{x^2+y^2, x^2+y^2>=1}, {{x,2},{y,2}}, Method->\"Powell\"]]",
             "FindMinimum");
}

/* ------------------------------------------------------------------ */
/* P6. Memory hygiene smoke                                            */
/* ------------------------------------------------------------------ */

static void test_no_leak_many_calls(void) {
    for (int i = 0; i < 30; i++) {
        Expr* e = parse_expression(
            "FindMinimum[(1-x)^2+100(y-x^2)^2, {{x,-1.2},{y,1}}, Method->\"Powell\"]");
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

    /* P1 unconstrained smooth */
    TEST(test_rosenbrock_2d);
    TEST(test_beale);
    TEST(test_booth);
    TEST(test_matyas);
    TEST(test_powell_singular);
    TEST(test_wood);
    TEST(test_trid_n6);

    /* P2 non-smooth */
    TEST(test_nonsmooth_l1);

    /* P3 box bounds */
    TEST(test_bound_corner);
    TEST(test_bound_single);

    /* P4 n == 1 parity */
    TEST(test_n1_matches_brent);

    /* P5 alias / maximum / shape / diagnostics */
    TEST(test_principalaxis_alias);
    TEST(test_findmaximum);
    TEST(test_shape);
    TEST(test_general_constraint_rejected);

    /* P6 memory */
    TEST(test_no_leak_many_calls);

    printf("All Powell tests passed!\n");
    return 0;
}
