/* Unit tests for the SLSQP method of FindMinimum / FindMaximum
 * (Method -> "SLSQP", alias "SequentialQuadraticProgramming";
 * src/numerical_calculus/findmin.c).
 *
 * SLSQP is a Han-Powell sequential quadratic programming method for smooth
 * CONSTRAINED optimization (equality + inequality + bounds). Each outer step
 * solves a QP built from the objective's quadratic model (a BFGS approximation
 * to the Hessian of the Lagrangian, kept SPD by Powell damping) and the
 * linearized constraints, via a Goldfarb-Idnani dual active-set QP; the step
 * is accepted by an L1 exact-penalty merit line search. Unlike the other
 * gradient methods, SLSQP consumes general (non-box) constraints DIRECTLY
 * through the QP rather than the augmented-Lagrangian penalty wrapper, so it
 * reaches the true constrained optimum (accurate constraint satisfaction) with
 * super-linear local convergence. With no constraints it degrades to a
 * damped-BFGS line-search solve.
 *
 * SLSQP is deterministic given a start point, so every case fixes its start
 * explicitly and pins an ANALYTIC optimum (KKT solution or an exact objective
 * value). The discriminating constrained cases -- the corner minimum at 7/3
 * (which the penalty method converges to only approximately), HS71 = 17.014,
 * the box-constrained Rosenbrock at (0.5, 0.25) -- exercise the parts that
 * distinguish SLSQP from the penalty-based path.
 *
 * Run binary directly: ./slsqp_tests */

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
/* S1. Unconstrained reduction (SLSQP with no constraints == BFGS)     */
/* ------------------------------------------------------------------ */

static void test_rosenbrock_2d(void) {
    /* min 0 at (1,1); with no constraints the QP is d = -B^-1 g and the merit
     * is f, so this is a damped-BFGS line-search solve. */
    check_true("Abs[First[FindMinimum[100 (y-x^2)^2+(1-x)^2, {{x,-1.2},{y,1}}, "
               "Method->\"SLSQP\", MaxIterations->2000]]] < 1.*^-8");
    check_true("With[{r=Last[FindMinimum[100 (y-x^2)^2+(1-x)^2, {{x,-1.2},{y,1}}, "
               "Method->\"SLSQP\", MaxIterations->2000]]}, "
               "Abs[(x/.r)-1]+Abs[(y/.r)-1] < 1.*^-4]");
}

static void test_illcond_quadratic(void) {
    /* condition 1e4 diagonal quadratic, min 0 at the origin. */
    check_true("Module[{v=Table[Symbol[\"z\"<>ToString[i]],{i,1,8}]}, "
               "Abs[First[FindMinimum[Evaluate[Sum[10^(4 (i-1)/7) v[[i]]^2,{i,1,8}]], "
               "Evaluate[Table[{v[[k]],1.0},{k,1,8}]], Method->\"SLSQP\", "
               "MaxIterations->2000]]] < 1.*^-6]");
}

static void test_sphere(void) {
    check_true("Abs[First[FindMinimum[x^2+y^2+z^2, {{x,3},{y,-2},{z,1}}, "
               "Method->\"SLSQP\"]]] < 1.*^-10");
}

/* ------------------------------------------------------------------ */
/* S2. Equality-constrained (analytic KKT)                             */
/* ------------------------------------------------------------------ */

static void test_eq_linear_line(void) {
    /* min x^2+y^2 s.t. x+y==1 -> 1/2 at (1/2,1/2). */
    check_true("Abs[First[FindMinimum[{x^2+y^2, x+y==1}, {{x,0},{y,0}}, "
               "Method->\"SLSQP\"]] - 0.5] < 1.*^-8");
    check_true("With[{r=Last[FindMinimum[{x^2+y^2, x+y==1}, {{x,0},{y,0}}, "
               "Method->\"SLSQP\"]]}, Abs[(x/.r)-0.5]+Abs[(y/.r)-0.5] < 1.*^-6]");
}

static void test_eq_circle(void) {
    /* min x+y s.t. x^2+y^2==1 -> -Sqrt[2] at (-1/Sqrt2,-1/Sqrt2). */
    check_true("Abs[First[FindMinimum[{x+y, x^2+y^2==1}, {{x,-1},{y,-1}}, "
               "Method->\"SLSQP\"]] + Sqrt[2.]] < 1.*^-7");
}

static void test_eq_offset_line(void) {
    /* min (x-2)^2+(y-1)^2 s.t. x+y==1 -> 2 at (1,0). */
    check_true("Abs[First[FindMinimum[{(x-2)^2+(y-1)^2, x+y==1}, {{x,0},{y,0}}, "
               "Method->\"SLSQP\"]] - 2.] < 1.*^-8");
}

/* ------------------------------------------------------------------ */
/* S3. Inequality-constrained (analytic active set)                    */
/* ------------------------------------------------------------------ */

static void test_ineq_corner(void) {
    /* min x+y s.t. 3x+2y>=7, x>=0, y>=0 -> 7/3 at (7/3,0). This is the exact
     * case the penalty method reaches only approximately (it lands near
     * x~2.45 once feasible); SLSQP drives to the corner via the active set. */
    check_true("Abs[First[FindMinimum[{x+y, 3 x+2 y>=7 && x>=0 && y>=0}, "
               "{{x,1},{y,1}}, Method->\"SLSQP\"]] - 7/3] < 1.*^-6");
    check_true("With[{r=Last[FindMinimum[{x+y, 3 x+2 y>=7 && x>=0 && y>=0}, "
               "{{x,1},{y,1}}, Method->\"SLSQP\"]]}, "
               "Abs[(x/.r)-7/3] + Abs[(y/.r)] < 1.*^-5]");
}

static void test_ineq_halfplane(void) {
    /* min x^2+y^2 s.t. x>=1 -> 1 at (1,0). */
    check_true("Abs[First[FindMinimum[{x^2+y^2, x>=1}, {{x,2},{y,1}}, "
               "Method->\"SLSQP\"]] - 1.] < 1.*^-8");
}

static void test_ineq_offset(void) {
    /* min (x-2)^2+(y-1)^2 s.t. x+y<=1 -> 2 at (1,0) (unconstrained min (2,1)
     * is infeasible, so the constraint is active on the line x+y==1). */
    check_true("Abs[First[FindMinimum[{(x-2)^2+(y-1)^2, x+y<=1}, {{x,0},{y,0}}, "
               "Method->\"SLSQP\"]] - 2.] < 1.*^-8");
}

/* ------------------------------------------------------------------ */
/* S4. Mixed equality + inequality + bounds (Hock-Schittkowski)        */
/* ------------------------------------------------------------------ */

static void test_hs71(void) {
    /* HS71 -- the canonical SLSQP test problem:
     *   min x1 x4 (x1+x2+x3) + x3
     *   s.t. x1 x2 x3 x4 >= 25, x1^2+x2^2+x3^2+x4^2 == 40, 1 <= xi <= 5.
     * Optimum f = 17.0140 at (1, 4.743, 3.821, 1.379). */
    check_true("Abs[First[FindMinimum[{x1 x4 (x1+x2+x3)+x3, "
               "x1 x2 x3 x4>=25 && x1^2+x2^2+x3^2+x4^2==40 && "
               "1<=x1<=5 && 1<=x2<=5 && 1<=x3<=5 && 1<=x4<=5}, "
               "{{x1,1},{x2,5},{x3,5},{x4,1}}, Method->\"SLSQP\", "
               "MaxIterations->500]] - 17.0140173] < 1.*^-3");
    /* First coordinate pins to the lower bound x1==1. */
    check_true("With[{r=Last[FindMinimum[{x1 x4 (x1+x2+x3)+x3, "
               "x1 x2 x3 x4>=25 && x1^2+x2^2+x3^2+x4^2==40 && "
               "1<=x1<=5 && 1<=x2<=5 && 1<=x3<=5 && 1<=x4<=5}, "
               "{{x1,1},{x2,5},{x3,5},{x4,1}}, Method->\"SLSQP\"]]}, "
               "Abs[(x1/.r)-1] < 1.*^-4]");
}

static void test_mixed_eq_ineq(void) {
    /* min (x-1)^2+(y-2)^2 s.t. y-x==1 && x+y<=2 -> 0.5 at (0.5,1.5)
     * (equality active always, inequality active at the optimum). */
    check_true("Abs[First[FindMinimum[{(x-1)^2+(y-2)^2, y-x==1 && x+y<=2}, "
               "{{x,0},{y,1}}, Method->\"SLSQP\"]] - 0.5] < 1.*^-7");
}

static void test_portfolio_simplex(void) {
    /* min w1^2+w2^2+w3^2 s.t. w1+w2+w3==1, wi>=0 -> 1/3 at (1/3,1/3,1/3). */
    check_true("Abs[First[FindMinimum[{w1^2+w2^2+w3^2, "
               "w1+w2+w3==1 && w1>=0 && w2>=0 && w3>=0}, "
               "{{w1,0.5},{w2,0.3},{w3,0.2}}, Method->\"SLSQP\"]] - 1/3] < 1.*^-7");
}

/* ------------------------------------------------------------------ */
/* S5. Box-only bounds                                                 */
/* ------------------------------------------------------------------ */

static void test_box_rosenbrock(void) {
    /* Rosenbrock s.t. x<=0.5 -> 0.25 at (0.5,0.25); the LBFGSB discriminating
     * case, cross-checking that SLSQP's bound rows reach the same optimum. */
    check_true("Abs[First[FindMinimum[{100 (y-x^2)^2+(1-x)^2, x<=0.5}, "
               "{{x,0},{y,0}}, Method->\"SLSQP\", MaxIterations->2000]] - 0.25] < 1.*^-6");
    check_true("With[{r=Last[FindMinimum[{100 (y-x^2)^2+(1-x)^2, x<=0.5}, "
               "{{x,0},{y,0}}, Method->\"SLSQP\", MaxIterations->2000]]}, "
               "Abs[(x/.r)-0.5]+Abs[(y/.r)-0.25] < 1.*^-4]");
}

static void test_box_bracket_spec(void) {
    /* {x, x0, xmin, xmax} bracket spec: min x^2 on [1,3] -> 1 at x==1. */
    check_true("Abs[First[FindMinimum[x^2, {x, 2, 1, 3}, Method->\"SLSQP\"]] - 1.] "
               "< 1.*^-8");
}

/* ------------------------------------------------------------------ */
/* S6. FindMaximum (objective negation)                                */
/* ------------------------------------------------------------------ */

static void test_findmaximum_eq(void) {
    /* max x y s.t. x+y==10 -> 25 at (5,5). */
    check_true("Abs[First[FindMaximum[{x y, x+y==10}, {{x,1},{y,1}}, "
               "Method->\"SLSQP\"]] - 25.] < 1.*^-6");
}

static void test_findmaximum_ineq(void) {
    /* max -(x-3)^2-(y+1)^2 s.t. x+y<=0 -> constrained max. Unconstrained max at
     * (3,-1) has x+y=2>0 infeasible; active on x+y==0. Maximize
     * -(x-3)^2-(y+1)^2 s.t. x+y=0 -> x-3=-(y+1)=-(-x+1)=x-1?? solve:
     * grad(-(x-3)^2-(y+1)^2)=(-2(x-3),-2(y+1))=mu(1,1) with y=-x
     * -> -(x-3) = -(y+1) -> x-3 = y+1 = -x+1 -> 2x=4 -> x=2, y=-2, f=-(−1)^2-(−1)^2=-2. */
    check_true("Abs[First[FindMaximum[{-(x-3)^2-(y+1)^2, x+y<=0}, {{x,0},{y,0}}, "
               "Method->\"SLSQP\"]] + 2.] < 1.*^-7");
}

/* ------------------------------------------------------------------ */
/* S7. Edge cases / options                                            */
/* ------------------------------------------------------------------ */

static void test_n1(void) {
    /* n==1 explicit SLSQP: min (x-3)^2 -> 0 at 3. */
    check_true("Abs[First[FindMinimum[(x-3)^2, {x,0}, Method->\"SLSQP\"]]] < 1.*^-10");
    check_true("With[{r=Last[FindMinimum[(x-3)^2, {x,0}, Method->\"SLSQP\"]]}, "
               "Abs[(x/.r)-3] < 1.*^-5]");
}

static void test_gradient_override(void) {
    /* User-supplied gradient reaches the same optimum. */
    check_true("Abs[First[FindMinimum[x^2+y^2, {{x,1},{y,1}}, Method->\"SLSQP\", "
               "Gradient->{2 x, 2 y}]]] < 1.*^-12");
}

static void test_alias(void) {
    /* "SequentialQuadraticProgramming" is the same method as "SLSQP". */
    check_true("Abs[First[FindMinimum[{x^2+y^2, x+y==1}, {{x,0},{y,0}}, "
               "Method->\"SequentialQuadraticProgramming\"]] - 0.5] < 1.*^-8");
}

static void test_parity_quasinewton(void) {
    /* Unconstrained: SLSQP and QuasiNewton agree to rounding at the optimum. */
    check_true("With[{a=First[FindMinimum[100 (y-x^2)^2+(1-x)^2, {{x,-1.2},{y,1}}, "
               "Method->\"SLSQP\", MaxIterations->2000]], "
               "b=First[FindMinimum[100 (y-x^2)^2+(1-x)^2, {{x,-1.2},{y,1}}, "
               "Method->\"QuasiNewton\", MaxIterations->2000]]}, Abs[a-b] < 1.*^-8]");
}

static void test_inconsistent_no_crash(void) {
    /* Contradictory equalities: SLSQP must return a well-formed List result
     * (best-effort / infeasible sentinel), not crash or hang. */
    check_eq("Head[FindMinimum[{x, x==0 && x==1}, {{x,0.5}}, Method->\"SLSQP\"]]",
             "List");
}

static void test_shape(void) {
    /* Result shape: {fmin, {rules...}} with a numeric first element. */
    check_true("MatchQ[FindMinimum[{x^2+y^2, x+y==1}, {{x,0},{y,0}}, "
               "Method->\"SLSQP\"], {_?NumberQ, {(_->_)..}}]");
}

static void test_no_leak_many_calls(void) {
    /* Exercise many solves (valgrind checks leaks separately); asserts the
     * repeated constrained solves stay correct. */
    check_true("Module[{s=0.}, Do[s += First[FindMinimum[{x^2+y^2, x+y==1}, "
               "{{x,0},{y,0}}, Method->\"SLSQP\"]], {50}]; Abs[s - 25.] < 1.*^-4]");
}

int main(void) {
    symtab_init();
    core_init();
    /* Numeric routines emit expected diagnostics (infeas, lstol, ...); keep
     * them off stdout so only FAIL lines show. */
    freopen("/dev/null", "w", stderr);

    /* S1 unconstrained reduction */
    TEST(test_rosenbrock_2d);
    TEST(test_illcond_quadratic);
    TEST(test_sphere);

    /* S2 equality-constrained */
    TEST(test_eq_linear_line);
    TEST(test_eq_circle);
    TEST(test_eq_offset_line);

    /* S3 inequality-constrained */
    TEST(test_ineq_corner);
    TEST(test_ineq_halfplane);
    TEST(test_ineq_offset);

    /* S4 mixed eq + ineq + bounds */
    TEST(test_hs71);
    TEST(test_mixed_eq_ineq);
    TEST(test_portfolio_simplex);

    /* S5 box-only */
    TEST(test_box_rosenbrock);
    TEST(test_box_bracket_spec);

    /* S6 FindMaximum */
    TEST(test_findmaximum_eq);
    TEST(test_findmaximum_ineq);

    /* S7 edge cases / options */
    TEST(test_n1);
    TEST(test_gradient_override);
    TEST(test_alias);
    TEST(test_parity_quasinewton);
    TEST(test_inconsistent_no_crash);
    TEST(test_shape);
    TEST(test_no_leak_many_calls);

    printf("All SLSQP tests passed!\n");
    return 0;
}
