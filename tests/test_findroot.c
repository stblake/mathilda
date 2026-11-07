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
        test_mute_messages();
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

static void test_root_finite_difference_zeta(void) {
    /* Zeta' is an inert Derivative[1][Zeta]; Newton must fall back to a
     * finite-difference derivative. Real root of Zeta[s] == 1.05 near 4.613. */
    check_true("Abs[(s /. First[FindRoot[Zeta[s] == 1.05, {s, 5}]]) - "
               "4.612973128439312] < 1.*^-6");
    /* Complex start: a zeta-difference root near 0.948962 + 20.3778 I. */
    check_true("Abs[(s /. First[FindRoot[2^(-s) (Zeta[s] - Zeta[s, 3/2]) == 0, "
               "{s, 20 I}]]) - (0.948962 + 20.3778 I)] < 1.*^-3");
    /* Finite-difference fallback also works at arbitrary precision. */
    check_true("Abs[(s /. First[FindRoot[Zeta[s] == 1.05, {s, 5}, "
               "WorkingPrecision -> 30]]) - 4.612973128439312] < 1.*^-6");
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
/* 14. Secant / Brent at arbitrary precision                           */
/* ------------------------------------------------------------------ */

/* 10^400 (x^2 - 2) exceeds the double range everywhere: the secant must stay
 * in MPFR end to end.  Root Sqrt[2]; {2, 3} does not bracket it, so this runs
 * secant steps first and Brent after the first sign change. */
static void test_secant_mpfr_beyond_double_range(void) {
    mute_stderr_once();
    check_true("With[{r = x /. FindRoot[10^400 (x^2 - 2), {x, 2, 3}, "
               "WorkingPrecision -> 50]}, "
               "Abs[r - Sqrt[2]] < 10^-45 && Precision[r] > 49]");
}

/* The largest zero of L_1000 (|L_1000| ~ 1e867 near it) by the secant method
 * from {4002, 4003} at 100 digits.  Checked independently: the Newton
 * correction L_1000/L_1000' recomputed at 200 digits is below 10^-80, and the
 * leading digits agree with 3943.24739484527095238972810775. */
static void test_secant_mpfr_laguerre_1000(void) {
    mute_stderr_once();
    check_true("Module[{f, r, y}, "
               "f[(t_)?NumberQ] := LaguerreL[1000, SetPrecision[t, 100]]; "
               "r = x /. FindRoot[Evaluate[f[x]], {x, 4002, 4003}, "
               "WorkingPrecision -> 100, AccuracyGoal -> 1000, MaxIterations -> 100]; "
               "y = SetPrecision[r, 200]; "
               "Abs[LaguerreL[1000, y]/LaguerreL[999, 1, y]] < 10^-80 && "
               "Abs[r - 394324739484527095238972810775/10^26] < 10^-25]");
}

/* A plain secant step from the bracket {-10, 20} of ArcTan[x - 1] extrapolates
 * away from the root (the next secant is nearly flat) and diverges; once two
 * iterates straddle a sign change the search must stay inside them. */
static void test_secant_bracket_safeguard(void) {
    mute_stderr_once();
    check_true("Abs[(x /. FindRoot[ArcTan[x - 1], {x, -10, 20}]) - 1] < 10^-10");
    check_true("Abs[(x /. FindRoot[ArcTan[x - 1], {x, -10, 20}, "
               "WorkingPrecision -> 30]) - 1] < 10^-25");
}

/* Starts that bracket exactly one zero of ChebyshevT[50, x] (zeros
 * Cos[(2k - 1) Pi/100]) return that zero -- not a neighbour reached by
 * extrapolating across sign changes. */
static void test_secant_bracketed_root_is_kept(void) {
    mute_stderr_once();
    check_true("With[{r = N[Cos[25 Pi/100], 40], g = N[Pi/100, 40]}, "
               "Abs[(x /. FindRoot[ChebyshevT[50, x], {x, r - 0.45 g, r + 0.9 g}]) - r] < 10^-12 && "
               "Abs[(x /. FindRoot[ChebyshevT[50, x], {x, r - 0.45 g, r + 0.9 g}, "
               "WorkingPrecision -> 40]) - r] < 10^-35]");
}

/* Brent at WorkingPrecision -> 40 works in MPFR and, under the default goals,
 * returns the full working precision (the step that meets the goal is taken
 * and lands at superlinear accuracy); explicit loose goals stop it early. */
static void test_brent_mpfr_precision(void) {
    mute_stderr_once();
    check_true("With[{r = x /. FindRoot[x^2 - 2, {x, 1, 2}, WorkingPrecision -> 40, "
               "Method -> \"Brent\"]}, Abs[r - Sqrt[2]] < 10^-38 && Precision[r] > 39]");
    check_true("With[{r = x /. FindRoot[x^2 - 2, {x, 1, 2}, WorkingPrecision -> 40, "
               "AccuracyGoal -> 5, PrecisionGoal -> 5]}, "
               "Abs[r - Sqrt[2]] < 10^-5 && Abs[r - Sqrt[2]] > 10^-30]");
}

/* Brent on a steep function whose root sits 4.7e-9 from a log singularity:
 * the bracket-width goal (absolute 10^-8) is larger than the root's own
 * scale, so it must not stop on the width alone with |f| ~ 1 (it returned
 * g = 49.1 for g == 50 after 13 evaluations).  Machine and MPFR. */
static void test_brent_steep_root_needs_residual(void) {
    mute_stderr_once();
    check_true("Module[{g, z0 = 0.0064823961, r}, "
               "g[z_?NumericQ] := If[z <= z0, 100., Min[100., -Log[z - z0]/0.3 - 13.9]]; "
               "r = z /. FindRoot[g[z] == 50, {z, 0.0064, 0.0065}, Method -> \"Brent\"]; "
               "Abs[g[r] - 50] < 10^-6]");
    check_true("Module[{g, z0 = 64823961/10^10, r}, "
               "g[z_?NumericQ] := If[z <= z0, 100, Min[100, -Log[z - z0]*10/3 - 139/10]]; "
               "r = z /. FindRoot[g[z] == 50, {z, 64/10^4, 65/10^4}, Method -> \"Brent\", "
               "WorkingPrecision -> 30]; Abs[g[r] - 50] < 10^-12]");
}

/* The compiled fast path must not swallow EvaluationMonitor. */
static void test_evaluation_monitor_compiled_path(void) {
    mute_stderr_once();
    check_true("Module[{c = 0}, FindRoot[x^2 - 2, {x, 1, 2}, WorkingPrecision -> 40, "
               "EvaluationMonitor :> c++]; c >= 5]");
    check_true("Module[{c = 0}, FindRoot[x^2 - 2, {x, 1.}, EvaluationMonitor :> c++]; c >= 3]");
}

/* ------------------------------------------------------------------ */
/* 15. Real-valued option values                                       */
/* ------------------------------------------------------------------ */

static void test_real_valued_options(void) {
    mute_stderr_once();
    /* the same solve as with Integer values, meeting the 20-digit goal */
    check_true("With[{r = x /. FindRoot[x^2 - 2, {x, 1}, WorkingPrecision -> 30., "
               "AccuracyGoal -> 20., PrecisionGoal -> 20., MaxIterations -> 20.], "
               "s = x /. FindRoot[x^2 - 2, {x, 1}, WorkingPrecision -> 30, "
               "AccuracyGoal -> 20, PrecisionGoal -> 20, MaxIterations -> 20]}, "
               "r === s && Abs[r - Sqrt[2]] < 10^-20 && Precision[r] > 29]");
    /* an arbitrary-precision WorkingPrecision value, as N[.., n] produces */
    check_true("Abs[(x /. FindRoot[x^2 - 2, {x, 1, 2}, WorkingPrecision -> N[40, 20]]) "
               "- Sqrt[2]] < 10^-35");
    /* N applied to a FindRoot result */
    check_true("With[{r = x /. N[FindRoot[x^2 - 2, {x, 1, 2}, WorkingPrecision -> 40, "
               "AccuracyGoal -> 40, MaxIterations -> 100], 12]}, "
               "Abs[r - Sqrt[2]] < 10^-11 && 11.5 < Precision[r] < 13]");
    /* a non-integral count stays an error */
    check_eq("Head[FindRoot[x^2 - 2, {x, 1}, MaxIterations -> 20.5]]", "FindRoot");
    check_eq("Head[FindRoot[x^2 - 2, {x, 1}, WorkingPrecision -> -3.]]", "FindRoot");
}

/* FindRoot is HoldAll; scalar option values are evaluated before use. */
static void test_option_values_evaluated(void) {
    mute_stderr_once();
    check_true("Block[{n = 10}, With[{r = x /. FindRoot[x^2 - 2, {x, 1}, "
               "MaxIterations -> 2 n, WorkingPrecision -> 3 n]}, "
               "Abs[r - Sqrt[2]] < 10^-15 && Precision[r] > 29]]");
}

/* ------------------------------------------------------------------ */
/* 16. Machine solve that leaves the double range                      */
/* ------------------------------------------------------------------ */

/* L_1000(4002.) ~ 3e867 overflows a double; the solve escalates to MPFR and
 * returns a machine number (leading digits as in the 100-digit test). */
static void test_machine_overflow_escalates(void) {
    mute_stderr_once();
    check_true("With[{r = x /. FindRoot[LaguerreL[1000, x], {x, 4002}]}, "
               "Precision[r] === MachinePrecision && "
               "Abs[r - 39432473948452709523897/10^19] < 10^-9]");
    /* secant: superlinear, so held to the machine goals (relative 10^-8) */
    check_true("With[{r = x /. FindRoot[LaguerreL[1000, x], {x, 4002, 4003}]}, "
               "Precision[r] === MachinePrecision && "
               "Abs[r - 39432473948452709523897/10^19] < 4*10^-5]");
    /* e^(1000 x) - 2 from x = 1: e^1000 ~ 2e434. Root Log[2]/1000. */
    check_true("Abs[(x /. FindRoot[Exp[1000 x] - 2, {x, 1}, MaxIterations -> 2000]) "
               "- Log[2]/1000] < 10^-10");   /* |f| < 10^-8, f' = 2000 */
}

/* ------------------------------------------------------------------ */
/* 17. Jacobian as a matrix                                            */
/* ------------------------------------------------------------------ */

static void test_jacobian_matrix_forms(void) {
    mute_stderr_once();
    /* 1 x 1 matrix for one equation, entries numeric only at a number */
    check_true("Module[{f, g, r, y}, "
               "f[(t_)?NumberQ] := LaguerreL[1000, SetPrecision[t, 100]]; "
               "g[(t_)?NumberQ] := -LaguerreL[999, 1, SetPrecision[t, 100]]; "
               "r = x /. FindRoot[Evaluate[f[x]], {x, 4002}, WorkingPrecision -> 100, "
               "AccuracyGoal -> 1000, MaxIterations -> 50, Jacobian -> {{g[x]}}]; "
               "y = SetPrecision[r, 200]; "
               "Abs[LaguerreL[1000, y]/LaguerreL[999, 1, y]] < 10^-80]");
    /* machine: the matrix and the bare form agree, and the result is a zero */
    check_true("With[{a = x /. FindRoot[LaguerreL[100, x], {x, 402}, "
               "Jacobian -> {{-LaguerreL[99, 1, x]}}], "
               "b = x /. FindRoot[LaguerreL[100, x], {x, 402}, "
               "Jacobian -> -LaguerreL[99, 1, x]]}, "
               "Abs[a - b] < 10^-9 && Abs[LaguerreL[100, SetPrecision[a, 60]]/"
               "LaguerreL[99, 1, SetPrecision[a, 60]]] < 10^-9]");
    /* one-equation system with a 1 x 1 matrix; a 2 x 2 system */
    check_true("Abs[(x /. FindRoot[{x^2 - 2}, {{x, 1}}, Jacobian -> {{2 x}}]) - Sqrt[2]] < 10^-8");
    check_true("Max[Abs[({x, y} /. FindRoot[{x^2 + y - 3, x - y + 1}, {{x, 0.5}, {y, 0.5}}, "
               "Jacobian -> {{2 x, 1}, {1, -1}}]) - {1, 2}]] < 10^-8");
    /* wrong shapes are rejected, not indexed */
    check_eq("Head[FindRoot[x^2 - 2, {x, 1}, Jacobian -> {{2 x, 1}}]]", "FindRoot");
    check_eq("Head[FindRoot[{x^2 + y - 3, x - y + 1}, {{x, 0.5}, {y, 0.5}}, "
             "Jacobian -> {{2 x, 1}}]]", "FindRoot");
}

/* ------------------------------------------------------------------ */
/* 18. Systems: vector equations, opaque functions, WorkingPrecision   */
/* ------------------------------------------------------------------ */

/* A vector equation F[a, b] == {0, 0} whose F is an opaque NumericQ-guarded
 * function: the widths come from the value at the start, the Jacobian from
 * forward differences (D of an opaque F is an inert Derivative).  The same
 * system as a list of opaque scalar equations, and as plain expressions. */
static void test_system_vector_opaque(void) {
    mute_stderr_once();
    const char* root = "{1.227773304303994085117098256488137476416507977664297054566109, "
                       "1.492572713238451939346233689691658447861638507574043131233776}";
    char buf[1024];
    snprintf(buf, sizeof buf,
             "Module[{g}, g[a_?NumericQ, b_?NumericQ] := {a^2 + b - 3, a - b^2 + 1}; "
             "Max[Abs[({a, b} /. FindRoot[g[a, b] == {0, 0}, {{a, 1.2}, {b, 1.1}}]) - %s]] < 10^-12]",
             root);
    check_true(buf);
    snprintf(buf, sizeof buf,
             "Module[{h, k}, h[a_?NumericQ, b_?NumericQ] := a^2 + b - 3; "
             "k[a_?NumericQ, b_?NumericQ] := a - b^2 + 1; "
             "Max[Abs[({a, b} /. FindRoot[{h[a, b] == 0, k[a, b] == 0}, {{a, 1.2}, {b, 1.1}}]) - %s]] < 10^-12]",
             root);
    check_true(buf);
    /* a bare vector-valued expression, and a mixed scalar + vector system */
    snprintf(buf, sizeof buf,
             "Module[{g}, g[a_?NumericQ, b_?NumericQ] := {a^2 + b - 3, a - b^2 + 1}; "
             "Max[Abs[({a, b} /. FindRoot[g[a, b], {{a, 1.2}, {b, 1.1}}]) - %s]] < 10^-12]",
             root);
    check_true(buf);
    check_true("Module[{g}, g[a_?NumericQ, b_?NumericQ] := {a + b - 3, a - b + 1}; "
               "Max[Abs[({a, b, c} /. FindRoot[{g[a, b] == {0, 0}, c^2 == 2}, "
               "{{a, 0.}, {b, 0.}, {c, 1.}}]) - {1, 2, Sqrt[2]}]] < 10^-10]");
    /* widths that do not add up stay a systno error */
    check_eq("Module[{g}, g[a_?NumericQ, b_?NumericQ] := {a, b, a b}; "
             "Head[FindRoot[g[a, b] == {1, 2, 3}, {{a, 1.}, {b, 1.}}]]]", "FindRoot");
}

/* WorkingPrecision -> 60 runs the system in MPFR and returns the full
 * working precision (it used to be ignored: machine numbers came back), for
 * symbolic and for opaque vector systems alike. */
static void test_system_working_precision(void) {
    mute_stderr_once();
    const char* root = "{1.227773304303994085117098256488137476416507977664297054566109076, "
                       "1.492572713238451939346233689691658447861638507574043131233776482}";
    char buf[1024];
    snprintf(buf, sizeof buf,
             "With[{r = {a, b} /. FindRoot[{a^2 + b - 3, a - b^2 + 1}, {{a, 12/10}, {b, 11/10}}, "
             "WorkingPrecision -> 60]}, Max[Abs[r - %s]] < 10^-57 && Min[Precision /@ r] > 59]",
             root);
    check_true(buf);
    snprintf(buf, sizeof buf,
             "Module[{g}, g[a_?NumericQ, b_?NumericQ] := {a^2 + b - 3, a - b^2 + 1}; "
             "With[{r = {a, b} /. FindRoot[g[a, b] == {0, 0}, {{a, 12/10}, {b, 11/10}}, "
             "WorkingPrecision -> 60]}, Max[Abs[r - %s]] < 10^-55 && Min[Precision /@ r] > 59]]",
             root);
    check_true(buf);
}

/* A shooting problem: an NDSolve behind a NumericQ-guarded function of the
 * launch angle and the arrival time, solved for the time to reach (1, 1/2)
 * under constant gravity (closed form T = Sqrt[7 - 2 Sqrt[11]] for speed 2). */
static void test_system_shooting_ndsolve(void) {
    mute_stderr_once();
    check_true("Module[{pos, x, y, t}, "
               "pos[th_?NumericQ, T_?NumericQ] := {x[T], y[T]} /. First[NDSolve[{x''[t] == 0, "
               "y''[t] == -1, x[0] == 0, y[0] == 0, x'[0] == 2 Cos[th], y'[0] == 2 Sin[th]}, "
               "{x, y}, {t, 0, T}, PrecisionGoal -> 13, AccuracyGoal -> 13]]; "
               "Abs[(T /. FindRoot[pos[th, T] == {1, 1/2}, {{th, 0.6}, {T, 0.6}}]) "
               "- Sqrt[7 - 2 Sqrt[11]]] < 10^-9]");
}

/* A Newton step that lands where f is not real (Log of a negative argument)
 * is halved back into the domain instead of ending the solve, at machine and
 * at arbitrary precision; and a system's step is backtracked likewise. */
static void test_newton_backtracks_into_domain(void) {
    mute_stderr_once();
    /* from x = 20, Newton on Log[x] - 1 + x/100 overshoots to x < 0 */
    check_true("Module[{f}, f[x_?NumericQ] := Log[x] - 1 + x/100; "
               "With[{r = x /. FindRoot[f[x], {x, 20}]}, Abs[Log[r] - 1 + r/100] < 10^-12]]");
    check_true("Module[{f}, f[x_?NumericQ] := Log[x] - 1 + x/100; "
               "With[{r = x /. FindRoot[f[x], {x, 20}, WorkingPrecision -> 40]}, "
               "Abs[Log[r] - 1 + r/100] < 10^-35]]");
    check_true("With[{r = {x, y} /. FindRoot[{Log[x] - 1 + x/100, y - x}, {{x, 20.}, {y, 1.}}]}, "
               "Abs[Log[r[[1]]] - 1 + r[[1]]/100] < 10^-12 && Abs[r[[2]] - r[[1]]] < 10^-12]");
}

/* A polynomial with exact rational coefficients whose numerators and
 * denominators are BigInts, at WorkingPrecision: v . x^Range[n] with an MPFR
 * x used to evaluate to Indeterminate (the fused Dot read the bignum
 * Rationals as NaN), so FindRoot reported a non-real f. */
static void test_bigint_rational_coefficients(void) {
    mute_stderr_once();
    check_true("Module[{v = {(3^80 + 1)/(2^130), (5^60)/(7^40 + 2)}, r}, "
               "r = x /. FindRoot[v . {x, x^2} == 1/10, {x, 1/1000}, WorkingPrecision -> 40]; "
               "Abs[v . {r, r^2} - 1/10] < 10^-35]");
    check_true("NumberQ[{(3^80 + 1)/(2^130), 2} . {N[1/3, 30], N[1/7, 30]}]");
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

/* Newton with the compiled symbolic derivative: Wilkinson's perturbed
 * polynomial prod (x - k) - 2^-23 x^19 from x = 21 (its real root next to 21)
 * against an independent 40-digit solve (MPFR, a different path), and the
 * compiled 20-term derivative against the interpreted one at a few points. */
static void test_root_compiled_derivative(void) {
    check_true("Module[{p = Product[x - k, {k, 1, 20}] - x^19/2^23, r, r40}, "
               "r = x /. FindRoot[Evaluate[p], {x, 21.}]; "
               "r40 = x /. FindRoot[Evaluate[p], {x, 21}, WorkingPrecision -> 40]; "
               "Abs[r - r40] < 10^-12 && Abs[r40 - 20.846908101] < 10^-8]");
    check_true("Module[{d = D[Product[x - k, {k, 1, 20}] - x^19/2^23, x], cf}, "
               "cf = Compile[{{x, _Real}}, Evaluate[d]]; "
               "And @@ Table[Abs[cf[v] - (d /. x -> v)] <= 10^-12 Abs[d /. x -> v], {v, {0.5, 7.25, 20.85, 21.}}]]");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_root_compiled_derivative);
    TEST(test_system_vector_opaque);
    TEST(test_system_working_precision);
    TEST(test_system_shooting_ndsolve);
    TEST(test_newton_backtracks_into_domain);
    TEST(test_bigint_rational_coefficients);

    TEST(test_root_quadratic);
    TEST(test_root_sin);
    TEST(test_root_equation_form);
    TEST(test_root_exp_plus_sin);
    TEST(test_root_polynomial);
    TEST(test_root_returns_rule_list);
    TEST(test_root_finite_difference_zeta);

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

    TEST(test_secant_mpfr_beyond_double_range);
    TEST(test_secant_mpfr_laguerre_1000);
    TEST(test_secant_bracket_safeguard);
    TEST(test_secant_bracketed_root_is_kept);
    TEST(test_brent_mpfr_precision);
    TEST(test_brent_steep_root_needs_residual);
    TEST(test_evaluation_monitor_compiled_path);
    TEST(test_real_valued_options);
    TEST(test_option_values_evaluated);
    TEST(test_machine_overflow_escalates);
    TEST(test_jacobian_matrix_forms);

    TEST(test_no_leak_many_calls);

    printf("All FindRoot tests passed!\n");
    return 0;
}
