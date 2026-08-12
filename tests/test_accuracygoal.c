/* Tests for the AccuracyGoal / PrecisionGoal contract shared by the numerical
 * calculus operations (nc_accuracy.{c,h}).
 *
 * Covers:
 *   - the four goal-value forms parse (a number, Automatic, Infinity,
 *     MachinePrecision) and an invalid value leaves the call unevaluated;
 *   - the default AccuracyGoal -> MachinePrecision appears in Options[];
 *   - MachinePrecision tracks WorkingPrecision (a WP-30 sum reaches ~28 digits,
 *     not a fixed 16), and the historical machine-precision default accuracy is
 *     preserved (regression guard for the geometric series);
 *   - the adaptive routines still return the correct value (and a finite number,
 *     never $Failed / unevaluated) on a cancellation-heavy limit whose accuracy
 *     goal cannot be met -- the warn-and-return-best contract.
 *
 * Numerical results are compared inside the language to avoid depending on the
 * printed form; the warning text itself goes to stderr and is not asserted.
 */

#include "core.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "test_utils.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char* eval_str(const char* input) {
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* e = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(e);
    expr_free(e);
    return s;
}

/* True iff `input` evaluates to the symbol True. */
static bool is_true(const char* input) {
    char* s = eval_str(input);
    bool ok = (strcmp(s, "True") == 0);
    if (!ok) fprintf(stderr, "  is_true FAIL: %s  =>  %s\n", input, s);
    free(s);
    return ok;
}

/* |input - expected| < tol, evaluated at full internal precision. */
static bool close_to(const char* input, const char* expected, double tol) {
    char buf[4096];
    snprintf(buf, sizeof buf,
             "N[Abs[Re[(%s) - (%s)]] + Abs[Im[(%s) - (%s)]]] < %.17g",
             input, expected, input, expected, tol);
    char* s = eval_str(buf);
    bool ok = (strcmp(s, "True") == 0);
    if (!ok) fprintf(stderr, "  close_to FAIL: |%s - (%s)| < %g  =>  %s\n",
                     input, expected, tol, s);
    free(s);
    return ok;
}

#define ASSERT_TRUE(input)  ASSERT_MSG(is_true((input)), "%s => True", (input))
#define ASSERT_CLOSE(input, expected, tol)                                   \
    ASSERT_MSG(close_to((input), (expected), (tol)),                         \
               "%s ~= %s (tol %g)", (input), (expected), (tol))

/* ---------------------------------------------------------------------- */

/* All four goal-value forms are accepted (the result is a number). */
static void test_goal_forms_parse(void) {
    ASSERT_TRUE("NumberQ[NLimit[Sin[x]/x, x -> 0, AccuracyGoal -> 20]]");
    ASSERT_TRUE("NumberQ[NLimit[Sin[x]/x, x -> 0, AccuracyGoal -> Automatic]]");
    ASSERT_TRUE("NumberQ[NLimit[Sin[x]/x, x -> 0, AccuracyGoal -> Infinity, "
                "PrecisionGoal -> Infinity]]");
    ASSERT_TRUE("NumberQ[NLimit[Sin[x]/x, x -> 0, AccuracyGoal -> MachinePrecision]]");
    ASSERT_TRUE("NumberQ[NSum[1/n^2, {n, 1, Infinity}, PrecisionGoal -> 10]]");
    ASSERT_TRUE("NumberQ[ND[Sin[x], x, 0, AccuracyGoal -> 12, PrecisionGoal -> 12]]");
}

/* An invalid goal value leaves the call unevaluated (the head survives). */
static void test_bad_goal_unevaluated(void) {
    char* s = eval_str("NLimit[Sin[x]/x, x -> 0, AccuracyGoal -> \"nonsense\"]");
    ASSERT_MSG(strstr(s, "NLimit") != NULL,
               "bad AccuracyGoal stays unevaluated, got %s", s);
    free(s);
}

/* The default AccuracyGoal is MachinePrecision on every wired-up head. */
static void test_default_is_machineprecision(void) {
    ASSERT_TRUE("MemberQ[Options[NLimit],   AccuracyGoal -> MachinePrecision]");
    ASSERT_TRUE("MemberQ[Options[NSum],     AccuracyGoal -> MachinePrecision]");
    ASSERT_TRUE("MemberQ[Options[ND],       AccuracyGoal -> MachinePrecision]");
    ASSERT_TRUE("MemberQ[Options[NSeries],  AccuracyGoal -> MachinePrecision]");
    ASSERT_TRUE("MemberQ[Options[NResidue], AccuracyGoal -> MachinePrecision]");
    ASSERT_TRUE("MemberQ[Options[NIntegrate], AccuracyGoal -> MachinePrecision]");
    /* The stale NLimit Terms default was corrected to 13. */
    ASSERT_TRUE("MemberQ[Options[NLimit], Terms -> 13]");
}

/* MachinePrecision tracks WorkingPrecision: a WP-30 sum reaches far past 16
 * digits, so the default goal is NOT a fixed machine-precision absolute floor. */
static void test_machineprecision_tracks_wp(void) {
    /* Sum_{n>=1} Log[1 + 1/n^2] = Log[Sinh[Pi]/Pi]; needs > 20 digits. */
    ASSERT_CLOSE("NSum[Log[1 + 1/n^2], {n, 1, Infinity}, WorkingPrecision -> 30]",
                 "Log[Sinh[Pi]/Pi]", 1e-25);
    /* Zeta(2) at high precision. */
    ASSERT_CLOSE("NSum[1/n^2, {n, 1, Infinity}, WorkingPrecision -> 30]",
                 "Pi^2/6", 1e-25);
}

/* The machine-precision default still delivers its historical accuracy
 * (regression guard: a loose PrecisionGoal must not cap it). */
static void test_machine_default_accuracy(void) {
    ASSERT_CLOSE("NSum[1/2^i, {i, 0, Infinity}]", "2", 1e-12);
    ASSERT_CLOSE("NSum[1/n^2, {n, 1, Infinity}]", "Pi^2/6", 1e-12);
    ASSERT_CLOSE("NLimit[(1 + 1/n)^n, n -> Infinity]", "E", 1e-9);
}

/* Warn-and-return-best: a cancellation-heavy limit cannot reach the machine
 * accuracy goal, but must still return the correct value as a finite number
 * (never $Failed / unevaluated). The NLimit::accgl warning goes to stderr. */
static void test_warn_returns_best(void) {
    ASSERT_TRUE("NumberQ[NLimit[(Tan[x] - Sin[x])/x^3, x -> 0]]");
    ASSERT_CLOSE("NLimit[(Tan[x] - Sin[x])/x^3, x -> 0]", "1/2", 1e-4);
    /* Cancellation-heavy ND: adaptive refinement must not follow a spuriously
     * small tableau residual into the round-off regime (ND[Exp,{x,2},0] -> 1,
     * not ~0). Value-consistency guard regression. */
    ASSERT_CLOSE("ND[Exp[x], {x, 2}, 0]", "1", 1e-6);
    ASSERT_CLOSE("ND[Cos[x], {x, 2}, 0]", "-1", 1e-6);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_goal_forms_parse);
    TEST(test_bad_goal_unevaluated);
    TEST(test_default_is_machineprecision);
    TEST(test_machineprecision_tracks_wp);
    TEST(test_machine_default_accuracy);
    TEST(test_warn_returns_best);

    printf("All accuracygoal_tests passed.\n");
    return 0;
}
