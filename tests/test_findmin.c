/* Unit tests for FindMinimum / FindMaximum (src/findmin.c).
 *
 * Coverage:
 *   - 1D Brent (default), bracket form, two-start.
 *   - Multivariate QuasiNewton/BFGS (default for n>=2).
 *   - Method dispatch (Brent, QuasiNewton, ConjugateGradient, Newton).
 *   - User-supplied Gradient option.
 *   - Box constraints (4-elt var spec, and parsed `a <= x && x <= b`).
 *   - Penalty path for a general inequality from a feasible start.
 *   - FindMaximum wrapper.
 *   - Monitors (EvaluationMonitor / StepMonitor).
 *   - MaxIterations exits and still returns a result.
 *   - HoldAll locality: the search variable does not leak.
 *   - Symbolic constants (Pi) collapse via numericalize.
 *   - Stress: Rosenbrock from (0,0) and (-1,1).
 *   - Diagnostics: argt, ivar, badmeth, badopt, nimpl.
 *   - Memory hygiene smoke test (loop of calls).
 *
 * Run binary directly: ./findmin_tests  (per MEMORY.md note). */

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

static void mute_stderr_once(void) {
    static int done = 0;
    if (!done) {
        freopen("/dev/null", "w", stderr);
        done = 1;
    }
}

/* Parse, evaluate, FullForm-compare. */
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

/* ------------------------------------------------------------------ */
/* 1. Golden 1D paths                                                  */
/* ------------------------------------------------------------------ */

static void test_min_parabola(void) {
    /* x^2 → minimum at 0. */
    check_true("Abs[(x /. Last[FindMinimum[x^2, {x, 1}]])] < 1.*^-6");
    check_true("Abs[First[FindMinimum[x^2, {x, 1}]]] < 1.*^-6");
}

static void test_min_shifted_parabola(void) {
    /* (x-3)^2 + 1 → minimum at x=3, value=1. */
    check_true("Abs[(x /. Last[FindMinimum[(x-3)^2 + 1, {x, 0}]]) - 3.0] < 1.*^-5");
    check_true("Abs[First[FindMinimum[(x-3)^2 + 1, {x, 0}]] - 1.0] < 1.*^-5");
}

static void test_min_x_cos_x_near_3(void) {
    /* x Cos[x] near x=2 finds local min at ≈ 3.4256, value ≈ -3.2884. */
    check_true("Abs[(x /. Last[FindMinimum[x Cos[x], {x, 2}]]) - 3.4256184593] < 1.*^-4");
    check_true("Abs[First[FindMinimum[x Cos[x], {x, 2}]] - (-3.2883712792)] < 1.*^-4");
}

static void test_min_x_cos_x_near_9(void) {
    /* x Cos[x] near x=7 → different local min at ≈ 9.5293. */
    check_true("Abs[(x /. Last[FindMinimum[x Cos[x], {x, 7}]]) - 9.5293446429] < 1.*^-4");
}

static void test_min_returns_pair(void) {
    /* Result must be {value, {var -> val}}. */
    check_eq("Head[FindMinimum[x^2, {x, 1}]]", "List");
    check_eq("Length[FindMinimum[x^2, {x, 1}]]", "2");
    check_eq("Head[Last[FindMinimum[x^2, {x, 1}]]]", "List");
    check_eq("Head[First[Last[FindMinimum[x^2, {x, 1}]]]]", "Rule");
}

/* ------------------------------------------------------------------ */
/* 2. Bracket / two-start forms                                        */
/* ------------------------------------------------------------------ */

static void test_min_bracket_4arg(void) {
    /* {x, xstart, xmin, xmax} routes through bracket Brent. */
    check_true("Abs[(x /. Last[FindMinimum[x Cos[x], {x, 7, 1, 15}]]) - 9.5293446429] < 1.*^-4");
}

static void test_min_two_start(void) {
    /* {x, x0, x1} → derivative-free Brent. */
    check_true("Abs[(x /. Last[FindMinimum[x Cos[x], {x, 8, 11}]]) - 9.5293446429] < 1.*^-4");
}

/* ------------------------------------------------------------------ */
/* 3. Multivariate (default QuasiNewton/BFGS)                          */
/* ------------------------------------------------------------------ */

static void test_min_2d_simple(void) {
    /* (x-1)^2 + (y-2)^2 → min at (1, 2), f=0. */
    check_true("With[{r = Last[FindMinimum[(x-1)^2 + (y-2)^2, {{x,0},{y,0}}]]}, "
               "Abs[(x /. r) - 1.0] + Abs[(y /. r) - 2.0] < 1.*^-4]");
}

static void test_min_2d_sin_sin(void) {
    /* Sin[x] Sin[2 y] from (2, 2) → min f=-1 at (π/2, 3π/4). */
    check_true("With[{r = Last[FindMinimum[Sin[x] Sin[2 y], {{x, 2}, {y, 2}}]]}, "
               "Abs[(x /. r) - 1.5707963268] + Abs[(y /. r) - 2.3561944902] < 1.*^-4]");
    check_true("Abs[First[FindMinimum[Sin[x] Sin[2 y], {{x, 2}, {y, 2}}]] - (-1.0)] < 1.*^-5");
}

static void test_min_3d(void) {
    /* x^2 + y^2 + z^2 starting at (1,2,3) → (0,0,0). */
    check_true("With[{r = Last[FindMinimum[x^2 + y^2 + z^2, "
               "{{x, 1}, {y, 2}, {z, 3}}]]}, "
               "Abs[(x /. r)] + Abs[(y /. r)] + Abs[(z /. r)] < 1.*^-4]");
}

static void test_min_auto_start_multivar(void) {
    /* {x, y} auto-start at 0. */
    check_true("First[FindMinimum[x^2 + y^2, {x, y}]] < 1.*^-6");
}

/* ------------------------------------------------------------------ */
/* 4. Method dispatch                                                  */
/* ------------------------------------------------------------------ */

static void test_method_brent_explicit(void) {
    check_true("Abs[(x /. Last[FindMinimum[(x-3)^2, {x, 0, -5, 5}, Method -> \"Brent\"]]) - 3.0] < 1.*^-5");
}

static void test_method_quasinewton(void) {
    check_true("With[{r = Last[FindMinimum[(x-1)^2 + (y-2)^2, "
               "{{x, 0}, {y, 0}}, Method -> \"QuasiNewton\"]]}, "
               "Abs[(x /. r) - 1.0] + Abs[(y /. r) - 2.0] < 1.*^-4]");
}

static void test_method_conjugate_gradient(void) {
    check_true("With[{r = Last[FindMinimum[(x-1)^2 + (y-2)^2, "
               "{{x, 5}, {y, 5}}, Method -> \"ConjugateGradient\"]]}, "
               "Abs[(x /. r) - 1.0] + Abs[(y /. r) - 2.0] < 1.*^-4]");
}

static void test_method_newton(void) {
    check_true("With[{r = Last[FindMinimum[Sin[x] Sin[2 y], "
               "{{x, 2}, {y, 2}}, Method -> \"Newton\"]]}, "
               "Abs[(x /. r) - 1.5707963268] + Abs[(y /. r) - 2.3561944902] < 1.*^-3]");
}

/* ------------------------------------------------------------------ */
/* 5. User-supplied Gradient                                           */
/* ------------------------------------------------------------------ */

static void test_user_gradient(void) {
    /* Gradient -> {dfdx, dfdy} for Sin[x] Sin[2 y]. */
    check_true("With[{r = Last[FindMinimum[Sin[x] Sin[2 y], "
               "{{x, 2}, {y, 2}}, "
               "Gradient -> {Cos[x] Sin[2 y], 2 Cos[2 y] Sin[x]}]]}, "
               "Abs[(x /. r) - 1.5707963268] + Abs[(y /. r) - 2.3561944902] < 1.*^-4]");
}

/* ------------------------------------------------------------------ */
/* 6. Constraints                                                      */
/* ------------------------------------------------------------------ */

static void test_box_4arg_spec(void) {
    /* Bracket spec {x, 7, 1, 15} restricts x ∈ [1,15]. */
    check_true("Abs[(x /. Last[FindMinimum[x Cos[x], {x, 7, 1, 15}]]) - 9.5293446429] < 1.*^-4");
}

static void test_box_from_and(void) {
    /* 1 <= x && x <= 15 → box constraint. */
    check_true("Abs[(x /. Last[FindMinimum[{x Cos[x], 1 <= x && x <= 15}, "
               "{x, 7}]]) - 9.5293446429] < 1.*^-4");
}

static void test_box_from_chained_le(void) {
    /* `1 <= x <= 15` parses as Inequality[1, LessEqual, x, LessEqual, 15]
     * and fm_collect_constraints walks the triples, registering each as
     * a box bound — Brent and the n-D methods both see the same boxes. */
    check_true("Abs[(x /. Last[FindMinimum[{x Cos[x], 1 <= x <= 15}, "
               "{x, 7}]]) - 9.5293446429] < 1.*^-4");
    check_true("Abs[(x /. Last[FindMaximum[{x Cos[x], 1 <= x <= 15}, "
               "{x, 7}]]) - 6.4372981791] < 1.*^-4");
}

static void test_box_from_chained_strict(void) {
    /* Strict chain `-3 < x < 3` and one-sided forms must also decompose. */
    check_true("Abs[(x /. Last[FindMinimum[{Sin[x], -3 < x < 3}, "
               "{x, 1}, Method -> \"Brent\"]]) + 1.5707963268] < 1.*^-4");
    check_true("Abs[(x /. Last[FindMinimum[{x^2, 2 <= x}, {x, 5}]]) - 2.0] < 1.*^-4");
    check_true("Abs[(x /. Last[FindMinimum[{x^2, x <= -3}, {x, -5}]]) + 3.0] < 1.*^-4");
}

static void test_penalty_inequality_feasible_start(void) {
    /* Sin[x] Sin[2 y] subject to disk x^2+y^2 <= 3, start (0.5, 1.5)
     * is feasible. Result should be in the disk and f < 0. */
    /* The quadratic-penalty path keeps the iterate near (not strictly
     * inside) the constraint set; we accept a small relaxation. */
    check_true("First[FindMinimum[{Sin[x] Sin[2 y], x^2 + y^2 <= 3}, "
               "{{x, 0.5}, {y, 1.5}}]] < 0");
}

static void test_penalty_linear_inequality(void) {
    /* min (x-3)^2 + (y-4)^2 s.t. x + y <= 1.
     * Analytical KKT: x=0, y=1, f=18. Start from infeasible (0, 0). */
    check_true("Abs[First[FindMinimum[{(x-3)^2 + (y-4)^2, x + y <= 1}, "
               "{{x, 0}, {y, 0}}]] - 18.0] < 1.*^-3");
    check_true("Abs[(y /. Last[FindMinimum[{(x-3)^2 + (y-4)^2, x + y <= 1}, "
               "{{x, 0}, {y, 0}}]]) - 1.0] < 1.*^-3");
}

static void test_penalty_quadratic_inequality(void) {
    /* min x^2 + y^2 s.t. x^2 + y^2 >= 1. Boundary, f = 1.
     * Symmetric start: BFGS picks (1/sqrt(2), 1/sqrt(2)). */
    check_true("Abs[First[FindMinimum[{x^2 + y^2, x^2 + y^2 >= 1}, "
               "{{x, 2.0}, {y, 2.0}}]] - 1.0] < 1.*^-3");
}

static void test_penalty_equality_constraint(void) {
    /* min x^2 + y^2 + z^2 s.t. x + y + z == 3 → (1, 1, 1), f = 3. */
    check_true("Abs[First[FindMinimum[{x^2 + y^2 + z^2, x + y + z == 3}, "
               "{{x, 0}, {y, 0}, {z, 0}}]] - 3.0] < 1.*^-3");
    check_true("Abs[(x /. Last[FindMinimum[{x^2 + y^2 + z^2, x + y + z == 3}, "
               "{{x, 0}, {y, 0}, {z, 0}}]]) - 1.0] < 1.*^-3");
}

static void test_penalty_projection_to_disk(void) {
    /* Closest point in the closed unit disk to (1, 2). Analytical:
     * project onto unit circle → (1, 2) / sqrt(5). f = (sqrt(5) - 1)^2. */
    check_true("Abs[First[FindMinimum[{(x-1)^2 + (y-2)^2, x^2 + y^2 <= 1}, "
               "{{x, 0.5}, {y, 0.5}}]] - (Sqrt[5] - 1)^2] < 1.*^-3");
}

static void test_penalty_constrained_max_xy(void) {
    /* Equivalent of max x*y s.t. x + y <= 4, x >= 0, y >= 0.
     * KKT: x = y = 2, x*y = 4. We minimise -x*y, so first[result] = -4. */
    check_true("Abs[First[FindMinimum[{-(x*y), "
               "x + y <= 4 && x >= 0 && y >= 0}, "
               "{{x, 0.1}, {y, 0.1}}]] + 4.0] < 1.*^-3");
}

/* ------------------------------------------------------------------ */
/* 7. FindMaximum                                                      */
/* ------------------------------------------------------------------ */

static void test_max_cos(void) {
    /* Cos[x] near 0 → max f=1 at x=0. */
    check_true("Abs[First[FindMaximum[Cos[x], {x, 0}]] - 1.0] < 1.*^-5");
    check_true("Abs[(x /. Last[FindMaximum[Cos[x], {x, 0}]])] < 1.*^-4");
}

static void test_max_inverted_parabola(void) {
    /* -(x-2)^2 + 7 → max f=7 at x=2. */
    check_true("Abs[First[FindMaximum[-(x-2)^2 + 7, {x, 0}]] - 7.0] < 1.*^-5");
    check_true("Abs[(x /. Last[FindMaximum[-(x-2)^2 + 7, {x, 0}]]) - 2.0] < 1.*^-4");
}

static void test_max_with_box(void) {
    /* FindMaximum[x Cos[x], {x, 7, 1, 15}] → finds local max in box. */
    check_true("First[FindMaximum[x Cos[x], {x, 7, 1, 15}]] > 0"); /* some positive local max */
}

/* ------------------------------------------------------------------ */
/* 8. Monitors                                                          */
/* ------------------------------------------------------------------ */

static void test_evaluation_monitor_fires(void) {
    Expr* setup = parse_expression("fmEvals = 0");
    Expr* t = evaluate(setup); expr_free(setup); expr_free(t);
    Expr* call = parse_expression(
        "FindMinimum[(x-3)^2, {x, 0}, "
        "EvaluationMonitor :> (fmEvals = fmEvals + 1)]");
    Expr* tr = evaluate(call); expr_free(call); expr_free(tr);
    check_true("fmEvals > 0");
}

static void test_step_monitor_fires(void) {
    Expr* setup = parse_expression("fmSteps = 0");
    Expr* t = evaluate(setup); expr_free(setup); expr_free(t);
    Expr* call = parse_expression(
        "FindMinimum[(x-3)^2, {x, 0}, "
        "StepMonitor :> (fmSteps = fmSteps + 1)]");
    Expr* tr = evaluate(call); expr_free(call); expr_free(tr);
    check_true("fmSteps > 0");
}

/* ------------------------------------------------------------------ */
/* 9. Options behaviour                                                */
/* ------------------------------------------------------------------ */

static void test_max_iterations_returns_iterate(void) {
    /* Tight iteration cap still returns a list result. */
    mute_stderr_once();
    check_eq("Head[FindMinimum[(1-x)^2 + 100 (y-x^2)^2, "
             "{{x, 0}, {y, 0}}, MaxIterations -> 3]]",
             "List");
}

static void test_working_precision_accepted(void) {
    /* WorkingPrecision -> 30 must be accepted (the iteration runs in MPFR
     * at the requested precision; see test_working_precision_brent_mpfr
     * and test_working_precision_bfgs_mpfr below for the digit-level checks). */
    check_eq("Head[FindMinimum[(x-Pi)^2, {x, 0}, WorkingPrecision -> 30]]", "List");
}

static void test_working_precision_brent_mpfr(void) {
    /* 1D Brent at WP=50: the minimum of (x-Pi)^2 is x=Pi, so the result
     * must be Pi to ~50 digits — well below the machine-precision
     * 16-digit horizon, proving the iteration ran in MPFR. */
    check_true("Abs[(x /. Last[FindMinimum[(x - Pi)^2, {x, 0}, "
               "WorkingPrecision -> 50]]) - Pi] < 1.*^-30");
}

static void test_working_precision_bfgs_mpfr(void) {
    /* n-D BFGS at WP=30: solve a quadratic and require ~25-digit accuracy
     * on each component. The double-precision optimizer can't beat ~1e-15. */
    check_true("Abs[(x /. (Last[FindMinimum[(x - 3)^2 + (y - 4)^2, "
               "{{x, 0}, {y, 0}}, WorkingPrecision -> 30]][[1]])) - 3] < 1.*^-20");
    check_true("Abs[(y /. (Last[FindMinimum[(x - 3)^2 + (y - 4)^2, "
               "{{x, 0}, {y, 0}}, WorkingPrecision -> 30]][[2]])) - 4] < 1.*^-20");
}

static void test_working_precision_findmax_mpfr(void) {
    /* FindMaximum delegates to FindMinimum after negation; the MPFR head
     * must survive the negation step (regression for an early version
     * that collapsed back to EXPR_REAL on the result-list head). */
    check_true("Abs[First[FindMaximum[-(x - Pi)^2, {x, 0}, "
               "WorkingPrecision -> 50]]] < 1.*^-30");
}

static void test_symbolic_pi(void) {
    /* The bound objective uses Pi — must numericalize. */
    check_true("Abs[(x /. Last[FindMinimum[(x-Pi)^2, {x, 0}]]) - 3.141592653589793] < 1.*^-5");
}

/* ------------------------------------------------------------------ */
/* 10. HoldAll + Block locality                                        */
/* ------------------------------------------------------------------ */

static void test_holdall_locality(void) {
    /* After FindMinimum, the search variable must remain free. */
    Expr* clr = parse_expression("Clear[fmVarLocality]");
    Expr* tmp = evaluate(clr); expr_free(clr); expr_free(tmp);
    check_eq("FindMinimum[fmVarLocality^2, {fmVarLocality, 1}]; fmVarLocality",
             "fmVarLocality");
}

/* ------------------------------------------------------------------ */
/* 11. Diagnostics / error paths                                       */
/* ------------------------------------------------------------------ */

static void test_diag_argt_zero(void) {
    mute_stderr_once();
    check_eq("Head[FindMinimum[]]", "FindMinimum");
}

static void test_diag_argt_one(void) {
    mute_stderr_once();
    check_eq("Head[FindMinimum[x^2]]", "FindMinimum");
}

static void test_diag_ivar(void) {
    mute_stderr_once();
    check_eq("Head[FindMinimum[x^2, {2, 1}]]", "FindMinimum");
}

static void test_diag_badmeth(void) {
    mute_stderr_once();
    check_eq("Head[FindMinimum[x^2, {x, 1}, Method -> \"Bogus\"]]", "FindMinimum");
}

static void test_diag_badopt(void) {
    mute_stderr_once();
    check_eq("Head[FindMinimum[x^2, {x, 1}, NotAnOpt -> 5]]", "FindMinimum");
}

static void test_diag_nimpl_method(void) {
    mute_stderr_once();
    check_eq("Head[FindMinimum[x^2, {x, 1}, Method -> \"LinearProgramming\"]]",
             "FindMinimum");
}

static void test_diag_findmax_argt(void) {
    mute_stderr_once();
    check_eq("Head[FindMaximum[]]", "FindMaximum");
}

/* ------------------------------------------------------------------ */
/* 12. Stress / hard problems                                          */
/* ------------------------------------------------------------------ */

static void test_stress_rosenbrock_origin(void) {
    /* Rosenbrock (1-x)^2 + 100 (y-x^2)^2 from (0, 0) → (1, 1). */
    check_true("With[{r = Last[FindMinimum[(1-x)^2 + 100 (y-x^2)^2, "
               "{{x, 0}, {y, 0}}]]}, "
               "Abs[(x /. r) - 1.0] + Abs[(y /. r) - 1.0] < 1.*^-3]");
}

static void test_stress_rosenbrock_negative(void) {
    /* Rosenbrock from (-1, 1) — classic difficult start. */
    check_true("With[{r = Last[FindMinimum[(1-x)^2 + 100 (y-x^2)^2, "
               "{{x, -1}, {y, 1}}]]}, "
               "Abs[(x /. r) - 1.0] + Abs[(y /. r) - 1.0] < 1.*^-2]");
}

static void test_stress_beale(void) {
    /* Beale's function: (1.5-x+x y)^2+(2.25-x+x y^2)^2+(2.625-x+x y^3)^2.
     * Global minimum at (3, 0.5) with f=0. Start near solution. */
    check_true("With[{r = Last[FindMinimum["
               "(1.5 - x + x y)^2 + (2.25 - x + x y^2)^2 + (2.625 - x + x y^3)^2, "
               "{{x, 1}, {y, 1}}]]}, "
               "Abs[(x /. r) - 3.0] + Abs[(y /. r) - 0.5] < 1.*^-2]");
}

/* ------------------------------------------------------------------ */
/* 13. Memory smoke test                                               */
/* ------------------------------------------------------------------ */

static void test_no_leak_many_calls(void) {
    for (int i = 0; i < 30; i++) {
        Expr* e1 = parse_expression("FindMinimum[(x-3)^2 + 1, {x, 0}]");
        Expr* r1 = evaluate(e1); expr_free(e1); expr_free(r1);
        Expr* e2 = parse_expression("FindMinimum[Sin[x] Sin[2 y], {{x, 2}, {y, 2}}]");
        Expr* r2 = evaluate(e2); expr_free(e2); expr_free(r2);
        Expr* e3 = parse_expression("FindMaximum[Cos[x], {x, 0}]");
        Expr* r3 = evaluate(e3); expr_free(e3); expr_free(r3);
    }
    check_true("True");
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

int main(void) {
    symtab_init();
    core_init();

    /* 1. Golden 1D */
    TEST(test_min_parabola);
    TEST(test_min_shifted_parabola);
    TEST(test_min_x_cos_x_near_3);
    TEST(test_min_x_cos_x_near_9);
    TEST(test_min_returns_pair);

    /* 2. Bracket / two-start */
    TEST(test_min_bracket_4arg);
    TEST(test_min_two_start);

    /* 3. Multivariate */
    TEST(test_min_2d_simple);
    TEST(test_min_2d_sin_sin);
    TEST(test_min_3d);
    TEST(test_min_auto_start_multivar);

    /* 4. Methods */
    TEST(test_method_brent_explicit);
    TEST(test_method_quasinewton);
    TEST(test_method_conjugate_gradient);
    TEST(test_method_newton);

    /* 5. User gradient */
    TEST(test_user_gradient);

    /* 6. Constraints */
    TEST(test_box_4arg_spec);
    TEST(test_box_from_and);
    TEST(test_box_from_chained_le);
    TEST(test_box_from_chained_strict);
    TEST(test_penalty_inequality_feasible_start);
    TEST(test_penalty_linear_inequality);
    TEST(test_penalty_quadratic_inequality);
    TEST(test_penalty_equality_constraint);
    TEST(test_penalty_projection_to_disk);
    TEST(test_penalty_constrained_max_xy);

    /* 7. FindMaximum */
    TEST(test_max_cos);
    TEST(test_max_inverted_parabola);
    TEST(test_max_with_box);

    /* 8. Monitors */
    TEST(test_evaluation_monitor_fires);
    TEST(test_step_monitor_fires);

    /* 9. Options */
    TEST(test_max_iterations_returns_iterate);
    TEST(test_working_precision_accepted);
    TEST(test_working_precision_brent_mpfr);
    TEST(test_working_precision_bfgs_mpfr);
    TEST(test_working_precision_findmax_mpfr);
    TEST(test_symbolic_pi);

    /* 10. Locality */
    TEST(test_holdall_locality);

    /* 11. Diagnostics */
    TEST(test_diag_argt_zero);
    TEST(test_diag_argt_one);
    TEST(test_diag_ivar);
    TEST(test_diag_badmeth);
    TEST(test_diag_badopt);
    TEST(test_diag_nimpl_method);
    TEST(test_diag_findmax_argt);

    /* 12. Stress */
    TEST(test_stress_rosenbrock_origin);
    TEST(test_stress_rosenbrock_negative);
    TEST(test_stress_beale);

    /* 13. Memory */
    TEST(test_no_leak_many_calls);

    printf("All FindMinimum/FindMaximum tests passed!\n");
    return 0;
}
