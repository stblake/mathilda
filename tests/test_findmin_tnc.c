/* Unit tests for the TNC method of FindMinimum / FindMaximum
 * (Method -> "TNC", alias "TruncatedNewton"; src/numerical_calculus/findmin.c).
 *
 * TNC is a Hessian-FREE truncated Newton method: an inner conjugate-gradient
 * loop approximately solves the Newton system H*p = -g using Hessian-vector
 * products (finite differences of the exact compiled gradient), with active-set
 * bound handling. It is gradient-based and deterministic given a start point, so
 * every case fixes its start explicitly. Every pinned target is an ANALYTIC
 * optimum (0, 18, 80, -50, 0.5, ...) or a relation between two solves.
 *
 * TNC's distinguishing strengths, exercised here: true Newton curvature on a
 * curved valley (Rosenbrock -> machine precision) and on ill-conditioned
 * quadratics (cond up to 1e6, where L-BFGS's low-rank model struggles). Unlike
 * the derivative-free Powell/NelderMead, general (non-box) constraints ARE
 * supported (via the augmented-Lagrangian penalty wrapper, like L-BFGS-B).
 *
 * Large objectives use genuine scalar symbols generated with
 * Symbol["z"<>ToString[i]], wrapped in Evaluate[...] past FindMinimum's HoldAll.
 *
 * Run binary directly: ./tnc_tests */

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
/* T1. Curved valleys — true Newton curvature via Hv                   */
/* ------------------------------------------------------------------ */

static void test_rosenbrock_2d(void) {
    /* min 0 at (1,1). TNC's inner CG solves the 2x2 Newton system exactly, so
     * unlike a steepest/quasi-Newton step it reaches (1,1) to near machine
     * precision — the point tolerance can be tight. */
    check_true("Abs[First[FindMinimum[(1-x)^2+100(y-x^2)^2, {{x,-1.2},{y,1}}, "
               "Method->\"TNC\", MaxIterations->2000]]] < 1.*^-8");
    check_true("With[{r=Last[FindMinimum[(1-x)^2+100(y-x^2)^2, {{x,-1.2},{y,1}}, "
               "Method->\"TNC\", MaxIterations->2000]]}, "
               "Abs[(x/.r)-1]+Abs[(y/.r)-1] < 1.*^-4]");
}

static void test_rosenbrock_extended_n10(void) {
    /* extended Rosenbrock n=10 from -1.2 → all-ones, min 0. */
    check_true("Module[{v=Table[Symbol[\"z\"<>ToString[i]],{i,1,10}]}, "
               "Abs[First[FindMinimum["
               "Evaluate[Sum[100(v[[j+1]]-v[[j]]^2)^2+(1-v[[j]])^2,{j,1,9}]], "
               "Evaluate[Table[{v[[k]],-1.2},{k,1,10}]], Method->\"TNC\", "
               "MaxIterations->3000]]] < 1.*^-5]");
}

/* ------------------------------------------------------------------ */
/* T2. Ill-conditioned quadratics — where TNC beats L-BFGS             */
/* ------------------------------------------------------------------ */

static void test_illcond_n10_c1e4(void) {
    check_true("Module[{v=Table[Symbol[\"z\"<>ToString[i]],{i,1,10}]}, "
               "Abs[First[FindMinimum[Evaluate[Sum[10^(4 (i-1)/9) v[[i]]^2,{i,1,10}]], "
               "Evaluate[Table[{v[[k]],1.0},{k,1,10}]], Method->\"TNC\", "
               "MaxIterations->2000]]] < 1.*^-6]");
}

static void test_illcond_n10_c1e6(void) {
    /* condition 1e6 — the case a low-rank quasi-Newton model resolves slowly;
     * TNC's true curvature drives it to ~0. */
    check_true("Module[{v=Table[Symbol[\"z\"<>ToString[i]],{i,1,10}]}, "
               "Abs[First[FindMinimum[Evaluate[Sum[10^(6 (i-1)/9) v[[i]]^2,{i,1,10}]], "
               "Evaluate[Table[{v[[k]],1.0},{k,1,10}]], Method->\"TNC\", "
               "MaxIterations->3000]]] < 1.*^-5]");
}

static void test_illcond_n50_c1e4(void) {
    check_true("Module[{v=Table[Symbol[\"z\"<>ToString[i]],{i,1,50}]}, "
               "Abs[First[FindMinimum[Evaluate[Sum[10^(4 (i-1)/49) v[[i]]^2,{i,1,50}]], "
               "Evaluate[Table[{v[[k]],1.0},{k,1,50}]], Method->\"TNC\", "
               "MaxIterations->2000]]] < 1.*^-6]");
}

/* ------------------------------------------------------------------ */
/* T3. Other standard smooth functions                                 */
/* ------------------------------------------------------------------ */

static void test_booth(void) {
    check_true("With[{r=Last[FindMinimum[(x+2 y-7)^2+(2 x+y-5)^2, {{x,0},{y,0}}, "
               "Method->\"TNC\"]]}, Abs[(x/.r)-1]+Abs[(y/.r)-3] < 1.*^-5]");
}

static void test_matyas(void) {
    check_true("Abs[First[FindMinimum[0.26(x^2+y^2)-0.48 x y, {{x,3},{y,3}}, "
               "Method->\"TNC\", MaxIterations->2000]]] < 1.*^-8");
}

static void test_wood(void) {
    check_true("Abs[First[FindMinimum["
               "100(y-x^2)^2+(1-x)^2+90(w-z^2)^2+(1-z)^2+10.1((y-1)^2+(w-1)^2)"
               "+19.8(y-1)(w-1), {{x,-3},{y,-1},{z,-3},{w,-1}}, "
               "Method->\"TNC\", MaxIterations->2000]]] < 1.*^-6");
}

static void test_trid_n6(void) {
    /* Trid n=6, analytic global min -50. */
    check_true("Module[{v=Table[Symbol[\"z\"<>ToString[i]],{i,1,6}]}, "
               "Abs[First[FindMinimum["
               "Evaluate[Sum[(v[[i]]-1)^2,{i,1,6}] - Sum[v[[i]] v[[i-1]],{i,2,6}]], "
               "Evaluate[Table[{v[[k]],0.0},{k,1,6}]], Method->\"TNC\", "
               "MaxIterations->2000]] - (-50.0)] < 1.*^-4]");
}

/* ------------------------------------------------------------------ */
/* T4. Box bounds (active-set)                                         */
/* ------------------------------------------------------------------ */

static void test_bound_corner(void) {
    /* unconstrained min (5,5) outside x<=2,y<=2; optimum at the corner (2,2)=18. */
    check_true("With[{r=FindMinimum[{(x-5)^2+(y-5)^2, x<=2 && y<=2}, "
               "{{x,0},{y,0}}, Method->\"TNC\"]}, "
               "Abs[First[r]-18.0] < 1.*^-6 && Abs[(x/.Last[r])-2]+Abs[(y/.Last[r])-2] < 1.*^-6]");
}

static void test_bound_many_active(void) {
    /* 5 variables each capped at 1, unconstrained min at 5 → all bounds active,
     * value 5*(5-1)^2 = 80 at all-ones. Stresses the active-set projection. */
    check_true("With[{r=FindMinimum[{(z1-5)^2+(z2-5)^2+(z3-5)^2+(z4-5)^2+(z5-5)^2, "
               "z1<=1 && z2<=1 && z3<=1 && z4<=1 && z5<=1}, "
               "{{z1,0},{z2,0},{z3,0},{z4,0},{z5,0}}, Method->\"TNC\"]}, "
               "Abs[First[r]-80.0] < 1.*^-6]");
}

/* ------------------------------------------------------------------ */
/* T5. General (non-box) constraints via the penalty wrapper           */
/* ------------------------------------------------------------------ */

static void test_general_constraint(void) {
    /* Unlike Powell/NelderMead (which reject general constraints), TNC — a
     * gradient method — routes them through the augmented-Lagrangian wrapper,
     * like L-BFGS-B. min x^2+y^2 s.t. x+y>=1 → (0.5,0.5), value 0.5. */
    check_true("With[{r=FindMinimum[{x^2+y^2, x+y>=1}, {{x,2},{y,2}}, Method->\"TNC\"]}, "
               "Abs[First[r]-0.5] < 1.*^-4 && Abs[(x/.Last[r])-0.5]+Abs[(y/.Last[r])-0.5] < 1.*^-3]");
}

/* ------------------------------------------------------------------ */
/* T6. n==1 (no delegation), parity, maximum, shape                    */
/* ------------------------------------------------------------------ */

static void test_n1(void) {
    /* TNC runs its general machinery at n==1 (no Brent delegation): one CG step
     * solves the 1x1 Newton system, and the negative-curvature guard steers it
     * off the maximum near the start. Sin[x] from x=2 → 3 Pi/2 ~ 4.712. */
    check_true("With[{r=FindMinimum[Sin[x], {x,2}, Method->\"TNC\"]}, "
               "Abs[First[r]-(-1.0)] < 1.*^-6 && Abs[(x/.Last[r])-4.712388980] < 1.*^-4]");
}

static void test_parity_quasinewton(void) {
    /* TNC and QuasiNewton reach the same minimiser on a smooth problem. */
    check_true("With[{a=Last[FindMinimum[(x-1)^2+(y-2)^2+(x y-1)^2, {{x,0},{y,0}}, Method->\"TNC\"]], "
               "b=Last[FindMinimum[(x-1)^2+(y-2)^2+(x y-1)^2, {{x,0},{y,0}}, Method->\"QuasiNewton\"]]}, "
               "Abs[(x/.a)-(x/.b)]+Abs[(y/.a)-(y/.b)] < 1.*^-4]");
}

static void test_findmaximum(void) {
    /* concave paraboloid, max 0 at (3,-1). */
    check_true("With[{r=FindMaximum[-(x-3)^2-(y+1)^2, {{x,0},{y,0}}, Method->\"TNC\"]}, "
               "Abs[First[r]] < 1.*^-8 && Abs[(x/.Last[r])-3]+Abs[(y/.Last[r])+1] < 1.*^-4]");
}

static void test_alias(void) {
    /* "TruncatedNewton" routes to the identical solver → bit-identical result. */
    check_true("First[FindMinimum[(1-x)^2+100(y-x^2)^2, {{x,-1.2},{y,1}}, "
               "Method->\"TruncatedNewton\", MaxIterations->2000]] == "
               "First[FindMinimum[(1-x)^2+100(y-x^2)^2, {{x,-1.2},{y,1}}, "
               "Method->\"TNC\", MaxIterations->2000]]");
}

static void test_shape(void) {
    check_eq("Length[FindMinimum[(x-1)^2+(y-2)^2, {{x,0},{y,0}}, Method->\"TNC\"]]", "2");
    check_eq("Length[Last[FindMinimum[(x-1)^2+(y-2)^2, {{x,0},{y,0}}, Method->\"TNC\"]]]", "2");
}

/* ------------------------------------------------------------------ */
/* T7. Memory hygiene smoke                                            */
/* ------------------------------------------------------------------ */

static void test_no_leak_many_calls(void) {
    for (int i = 0; i < 30; i++) {
        Expr* e = parse_expression(
            "FindMinimum[{(x-2)^2+(y-3)^2, x+y<=1}, {{x,0},{y,0}}, Method->\"TNC\"]");
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
    /* Numeric routines emit expected diagnostics (penalty chatter, lstol, ...);
     * keep them off stdout so only FAIL lines show. */
    freopen("/dev/null", "w", stderr);

    /* T1 curved valleys */
    TEST(test_rosenbrock_2d);
    TEST(test_rosenbrock_extended_n10);

    /* T2 ill-conditioned */
    TEST(test_illcond_n10_c1e4);
    TEST(test_illcond_n10_c1e6);
    TEST(test_illcond_n50_c1e4);

    /* T3 standard smooth */
    TEST(test_booth);
    TEST(test_matyas);
    TEST(test_wood);
    TEST(test_trid_n6);

    /* T4 box bounds */
    TEST(test_bound_corner);
    TEST(test_bound_many_active);

    /* T5 general constraints via penalty */
    TEST(test_general_constraint);

    /* T6 n==1 / parity / maximum / alias / shape */
    TEST(test_n1);
    TEST(test_parity_quasinewton);
    TEST(test_findmaximum);
    TEST(test_alias);
    TEST(test_shape);

    /* T7 memory */
    TEST(test_no_leak_many_calls);

    printf("All TNC tests passed!\n");
    return 0;
}
