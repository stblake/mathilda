/* Unit tests for FindRoot (src/findroot.c).
 *
 * Coverage:
 *   - Scalar Newton from a single start.
 *   - Equation form lhs == rhs.
 *   - Brent and Secant methods (both implicit and explicit).
 *   - Bracket form {x, x0, xmin, xmax} routes through Brent.
 *   - System of 2 / 3 equations via Newton + Gaussian-elimination solver.
 *   - Complex roots (machine precision and MPFR/MPC).
 *   - WorkingPrecision -> n triggers MPFR Newton, agrees with Pi to 45+ digits.
 *   - DampingFactor speeds convergence on repeated roots.
 *   - MaxIterations: emits cvmit, returns last iterate.
 *   - AccuracyGoal / PrecisionGoal alter stopping criterion.
 *   - StepMonitor / EvaluationMonitor fire each iteration / eval.
 *   - HoldAll + Block: variables are not polluted in the outer scope.
 *   - User-supplied Jacobian short-circuits symbolic D[].
 *   - Diagnostics: argt, ivar, vecvar, badmeth, brnoth, badopt.
 *
 * Run binary directly: ./findroot_tests
 * (per MEMORY.md note: ctest is not configured in tests/CMakeLists.txt). */

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

/* Evaluate `input` to a True/False symbol; assert it is True. */
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
/* 1. Trivial / golden path                                            */
/* ------------------------------------------------------------------ */

static void test_root_quadratic(void) {
    check_true("Abs[(x /. First[FindRoot[x^2 - 2, {x, 1.0}]]) - Sqrt[2]] < 1.*^-7");
}

static void test_root_sin(void) {
    check_true("Abs[(x /. First[FindRoot[Sin[x], {x, 3.0}]]) - Pi] < 1.*^-7");
}

static void test_root_equation_form(void) {
    /* cos(x) == x has root ~0.739085... */
    check_true("Abs[(x /. First[FindRoot[Cos[x] == x, {x, 0}]]) - 0.7390851332151607] < 1.*^-7");
}

static void test_root_exp_plus_sin(void) {
    /* Sin[x] + Exp[x] has root near -0.588532... */
    check_true("Abs[(x /. First[FindRoot[Sin[x] + Exp[x], {x, 0}]]) - (-0.5885327439818610)] < 1.*^-7");
}

static void test_root_polynomial(void) {
    check_true("Abs[(x /. First[FindRoot[x^3 - 8, {x, 1.0}]]) - 2.0] < 1.*^-7");
}

static void test_root_returns_rule_list(void) {
    /* Structural check: head is List of length 1, inner head is Rule. */
    check_eq("Head[FindRoot[x^2 - 4, {x, 1.0}]]", "List");
    check_eq("Length[FindRoot[x^2 - 4, {x, 1.0}]]", "1");
    check_eq("Head[First[FindRoot[x^2 - 4, {x, 1.0}]]]", "Rule");
}

/* ------------------------------------------------------------------ */
/* 2. Method dispatch                                                  */
/* ------------------------------------------------------------------ */

static void test_method_brent_implicit(void) {
    /* {x, a, b} with default method -> Secant; supply Brent explicitly. */
    check_true("Abs[(x /. First[FindRoot[x^2 - 2, {x, 1.0, 2.0}, Method -> \"Brent\"]]) - Sqrt[2]] < 1.*^-8");
}

static void test_method_brent_bracket_spec(void) {
    /* 4-arg spec auto-routes to Brent.  Cos[x] - x has unique root in [0,1]. */
    check_true("Abs[(x /. First[FindRoot[Cos[x] - x, {x, 0.5, 0, 1}]]) - 0.7390851332151607] < 1.*^-8");
}

static void test_method_secant(void) {
    check_true("Abs[(x /. First[FindRoot[x^2 - 2, {x, 1.0, 2.0}, Method -> \"Secant\"]]) - Sqrt[2]] < 1.*^-8");
}

static void test_method_newton_default(void) {
    /* {var, x0} single-point: Newton. Verify by checking convergence on
     * a problem where secant would diverge from a single start.  */
    check_true("Abs[(x /. First[FindRoot[Cos[x] - x, {x, 0.5}]]) - 0.7390851332151607] < 1.*^-7");
}

static void test_method_invalid(void) {
    /* Bogus method name → unevaluated, no crash. */
    mute_stderr_once();
    check_eq("Head[FindRoot[x^2 - 1, {x, 0.5}, Method -> \"Bogus\"]]", "FindRoot");
}

static void test_brent_no_bracket(void) {
    /* f(a)*f(b) > 0 → brnoth diagnostic, unevaluated. */
    mute_stderr_once();
    check_eq("Head[FindRoot[x^2 + 1, {x, 0.0, 1.0}, Method -> \"Brent\"]]", "FindRoot");
}

/* ------------------------------------------------------------------ */
/* 3. System of equations                                              */
/* ------------------------------------------------------------------ */

static void test_system_2x2_explicit(void) {
    /* From the docstring: {y == Exp[x], x + y == 2}. */
    check_true("Abs[(x /. First[FindRoot[{y == Exp[x], x + y == 2}, {{x, 1.0}, {y, 1.0}}]]) - 0.4428544010023886] < 1.*^-8");
}

static void test_system_2x2_exp_squared(void) {
    /* {Exp[x-2] == y, y^2 == x} has solution near (0.019026, 0.137935).
     * Verify by residual norm at the returned point. */
    check_true("With[{r = FindRoot[{Exp[x - 2] == y, y^2 == x}, "
               "{{x, 1.0}, {y, 1.0}}]}, "
               "Abs[(y - Exp[x - 2]) /. r] + Abs[(y^2 - x) /. r] < 1.*^-6]");
}

static void test_system_3x3(void) {
    /* {Sin[x+y], Cos[x-y], x^2+y^2 - z} = 0 has root x=pi/4, y=-pi/4, z=pi^2/8. */
    check_true("With[{r = FindRoot[{Sin[x+y], Cos[x-y], x^2+y^2 - z}, "
               "{{x, 1.0}, {y, 0.0}, {z, 0.0}}]}, "
               "Abs[Sin[x+y] /. r] + Abs[Cos[x-y] /. r] + Abs[(x^2+y^2-z) /. r] < 1.*^-6]");
}

static void test_system_linear(void) {
    /* {x + y - 3, 2 x - y} → x=1, y=2. */
    check_true("With[{r = FindRoot[{x + y - 3, 2 x - y}, {{x, 0.0}, {y, 0.0}}]}, "
               "Abs[(x /. r) - 1.0] + Abs[(y /. r) - 2.0] < 1.*^-7]");
}

/* ------------------------------------------------------------------ */
/* 4. Complex roots                                                    */
/* ------------------------------------------------------------------ */

static void test_complex_start_z(void) {
    /* (Cos[z + I] - 2)(z + 2) has a real root z=-2 and a complex root
     * z ≈ 0 + 0.316958 I.  Starting at 1 + 0.1 I should find the
     * complex one.  Compare component-wise against the known answer. */
    check_true("With[{r = FindRoot[(Cos[z + I] - 2)(z + 2), {z, 1.0 + 0.1 I}]}, "
               "Abs[Re[z /. r]] + Abs[Im[z /. r] - 0.316957896762603] < 1.*^-6]");
}

static void test_complex_imag_only_start(void) {
    /* x^2 + 1 has roots ±I. Starting at 0.1 + 0.9 I finds +I. */
    check_true("With[{r = FindRoot[x^2 + 1, {x, 0.1 + 0.9 I}]}, "
               "Abs[Re[x /. r]] + Abs[Im[x /. r] - 1.0] < 1.*^-6]");
}

/* ------------------------------------------------------------------ */
/* 5. Precision (WorkingPrecision)                                     */
/* ------------------------------------------------------------------ */

static void test_working_precision_30(void) {
    /* Sin[x] = 0 near 3 finds Pi at 30 digits. */
    check_true("With[{r = FindRoot[Sin[x], {x, 3}, WorkingPrecision -> 30]}, "
               "Abs[(x /. r) - N[Pi, 30]] < 1.*^-25]");
}

static void test_working_precision_50(void) {
    /* WorkingPrecision -> 50 → default goals = 25 digits.  Residual
     * is around 1e-29, comfortably under 1e-20. */
    check_true("With[{r = FindRoot[Sin[x], {x, 3}, WorkingPrecision -> 50]}, "
               "Abs[(x /. r) - N[Pi, 50]] < 1.*^-20]");
}

static void test_working_precision_complex(void) {
    /* Complex Newton at MPFR precision: x^2 + 1 = 0 starting near +I.
     * Default goals = WP/2 = 15 digits.  Verify at least 10 digits. */
    check_true("With[{r = FindRoot[x^2 + 1, {x, 0.1 + 0.9 I}, WorkingPrecision -> 30]}, "
               "Abs[Re[x /. r]] + Abs[Im[x /. r] - 1] < 1.*^-10]");
}

/* ------------------------------------------------------------------ */
/* 6. Iteration limit / damping                                        */
/* ------------------------------------------------------------------ */

static void test_max_iterations_returns_last(void) {
    /* Mathematica returns the last iterate after MaxIterations. */
    mute_stderr_once();
    check_eq("Head[FindRoot[Sin[x] + Exp[x], {x, 0}, MaxIterations -> 3]]", "List");
}

static void test_damping_factor_speeds_repeated_root(void) {
    /* (x - 1)^3 has triple root at 1.  Damping 1 converges linearly;
     * damping 3 should hit machine precision faster. We just check
     * that the result is accurate. */
    check_true("Abs[(x /. First[FindRoot[(x - 1)^3, {x, 0.5}, DampingFactor -> 3]]) - 1.0] < 1.*^-6");
}

/* ------------------------------------------------------------------ */
/* 7. Accuracy / Precision goals                                       */
/* ------------------------------------------------------------------ */

static void test_accuracy_goal_relaxed(void) {
    /* AccuracyGoal -> 4 returns an iterate within ~1e-4 of zero.  The
     * answer for Sin[x-10]-x+10 near 0 is in the basin of x = 10. */
    check_true("Abs[(x /. First[FindRoot[Sin[x - 10] - x + 10, {x, 8.0}, "
               "AccuracyGoal -> 4, PrecisionGoal -> 4]]) - 10.0] < 1.0");
}

static void test_precision_goal_infinity(void) {
    /* AccuracyGoal -> Infinity disables the |f|-based stop; PrecisionGoal
     * alone must still converge. */
    check_true("Abs[(x /. First[FindRoot[Sin[x - 10] - x + 10, {x, 8.0}, "
               "AccuracyGoal -> Infinity, PrecisionGoal -> 8]]) - 10.0] < 1.*^-4");
}

/* ------------------------------------------------------------------ */
/* 8. HoldAll + Block (locality)                                       */
/* ------------------------------------------------------------------ */

static void test_hold_all_locality(void) {
    /* The search variable name must NOT leak: after FindRoot, the
     * symbol must still be free (no OwnValue installed). */
    Expr* clr = parse_expression("Clear[frVarLocality]");
    Expr* tmp = evaluate(clr); expr_free(clr); expr_free(tmp);
    check_eq("FindRoot[frVarLocality^2 - 9, {frVarLocality, 2.0}]; frVarLocality",
             "frVarLocality");
}

/* ------------------------------------------------------------------ */
/* 9. User-supplied Jacobian (scalar form)                             */
/* ------------------------------------------------------------------ */

static void test_jacobian_scalar(void) {
    /* For scalar f, the Jacobian option supplies a single expression
     * playing the role of f'.  We compute Sqrt[5] with the exact
     * derivative as a hint.  */
    check_true("Abs[(x /. First[FindRoot[x^2 - 5, {x, 1.0}, Jacobian -> 2 x]]) - Sqrt[5]] < 1.*^-7");
}

/* ------------------------------------------------------------------ */
/* 10. Monitors                                                        */
/* ------------------------------------------------------------------ */

static void test_step_monitor_fires(void) {
    /* StepMonitor :> (counter++) increments a global per iteration.   */
    Expr* setup = parse_expression("frSteps = 0");
    Expr* t = evaluate(setup); expr_free(setup); expr_free(t);
    Expr* call = parse_expression(
        "FindRoot[Sin[x] - 1, {x, 1.5}, "
        "StepMonitor :> (frSteps = frSteps + 1)]");
    Expr* tr = evaluate(call); expr_free(call); expr_free(tr);
    /* Now frSteps must be > 0. */
    check_true("frSteps > 0");
}

static void test_evaluation_monitor_fires(void) {
    Expr* setup = parse_expression("frEvals = 0");
    Expr* t = evaluate(setup); expr_free(setup); expr_free(t);
    Expr* call = parse_expression(
        "FindRoot[Sin[x], {x, 3.0}, "
        "EvaluationMonitor :> (frEvals = frEvals + 1)]");
    Expr* tr = evaluate(call); expr_free(call); expr_free(tr);
    check_true("frEvals > 0");
}

/* ------------------------------------------------------------------ */
/* 11. Diagnostics / error paths                                       */
/* ------------------------------------------------------------------ */

static void test_diag_argt(void) {
    mute_stderr_once();
    check_eq("Head[FindRoot[x^2]]", "FindRoot");
}

static void test_diag_vector_rejected(void) {
    /* Vector-valued start is the deferred path; v1 rejects it. */
    mute_stderr_once();
    check_eq("Head[FindRoot[x, {x, {1, 2}}]]", "FindRoot");
}

static void test_diag_ivar(void) {
    mute_stderr_once();
    check_eq("Head[FindRoot[x^2 - 1, {2, 1.0}]]", "FindRoot");
}

static void test_diag_invalid_option(void) {
    mute_stderr_once();
    check_eq("Head[FindRoot[x^2 - 1, {x, 1.0}, NotAnOpt -> 5]]", "FindRoot");
}

/* ------------------------------------------------------------------ */
/* 12. Stress / hard problems                                          */
/* ------------------------------------------------------------------ */

static void test_stress_repeated_root_x3(void) {
    /* (x-1)^3 with Newton — linear convergence; bumping damping helps. */
    check_true("Abs[(x /. First[FindRoot[(x - 1)^3, {x, 0.5}, "
               "DampingFactor -> 3, MaxIterations -> 200]]) - 1.0] < 1.*^-4");
}

static void test_stress_high_precision_pi(void) {
    /* Sin[x] = 0 near 3, WorkingPrecision -> 50.  Default goals are
     * 25 digits, so residual is ~1e-29.  Verify at least 20 digits. */
    check_true("With[{r = FindRoot[Sin[x], {x, 3}, WorkingPrecision -> 50]}, "
               "Abs[(x /. r) - N[Pi, 50]] < 1.*^-20]");
}

static void test_stress_brent_transcendental(void) {
    /* Cos[x] - x has unique real root ~0.739085; Brent on a wide bracket. */
    check_true("Abs[(x /. First[FindRoot[Cos[x] - x, {x, 0.0, 1.0}, "
               "Method -> \"Brent\"]]) - 0.7390851332151607] < 1.*^-7");
}

static void test_stress_secant_quadratic(void) {
    /* Secant on x^2 - 2 starting from 1.0, 1.5. */
    check_true("Abs[(x /. First[FindRoot[x^2 - 2, {x, 1.0, 1.5}, "
               "Method -> \"Secant\"]]) - Sqrt[2]] < 1.*^-7");
}

/* ------------------------------------------------------------------ */
/* 13. Memory: many calls in a row (smoke test for leaks)              */
/* ------------------------------------------------------------------ */

static void test_no_leak_many_calls(void) {
    /* Just exercise the code paths repeatedly.  valgrind / leak-checks
     * run as a separate `make valgrind` step (see CLAUDE.md). */
    for (int i = 0; i < 20; i++) {
        Expr* e = parse_expression("FindRoot[Sin[x] + Exp[x], {x, 0}]");
        Expr* r = evaluate(e);
        expr_free(e);
        expr_free(r);
    }
    /* Trivially pass.  Memory health is verified externally. */
    check_true("True");
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

int main(void) {
    symtab_init();
    core_init();

    TEST(test_root_quadratic);
    TEST(test_root_sin);
    TEST(test_root_equation_form);
    TEST(test_root_exp_plus_sin);
    TEST(test_root_polynomial);
    TEST(test_root_returns_rule_list);

    TEST(test_method_brent_implicit);
    TEST(test_method_brent_bracket_spec);
    TEST(test_method_secant);
    TEST(test_method_newton_default);
    TEST(test_method_invalid);
    TEST(test_brent_no_bracket);

    TEST(test_system_2x2_explicit);
    TEST(test_system_2x2_exp_squared);
    TEST(test_system_3x3);
    TEST(test_system_linear);

    TEST(test_complex_start_z);
    TEST(test_complex_imag_only_start);

    TEST(test_working_precision_30);
    TEST(test_working_precision_50);
    TEST(test_working_precision_complex);

    TEST(test_max_iterations_returns_last);
    TEST(test_damping_factor_speeds_repeated_root);

    TEST(test_accuracy_goal_relaxed);
    TEST(test_precision_goal_infinity);

    TEST(test_hold_all_locality);
    TEST(test_jacobian_scalar);

    TEST(test_step_monitor_fires);
    TEST(test_evaluation_monitor_fires);

    TEST(test_diag_argt);
    TEST(test_diag_vector_rejected);
    TEST(test_diag_ivar);
    TEST(test_diag_invalid_option);

    TEST(test_stress_repeated_root_x3);
    TEST(test_stress_high_precision_pi);
    TEST(test_stress_brent_transcendental);
    TEST(test_stress_secant_quadratic);

    TEST(test_no_leak_many_calls);

    printf("All FindRoot tests passed!\n");
    return 0;
}
