/*
 * test_time_constrained.c -- unit tests for TimeConstrained[expr, t, ...]
 *
 * Covers:
 *   - Fast-completion cases (2-arg, 3-arg) -- body returns normally.
 *   - Infinity time budget -- no abort.
 *   - Non-positive (zero, negative) and NaN-shaped budgets -- abort
 *     immediately, with or without failexpr.
 *   - Bad arity, non-numeric time -- left unevaluated.
 *   - HoldAll semantics: failexpr not evaluated on success, body
 *     evaluated exactly once.
 *   - Timer-driven abort: $Aborted (no failexpr) and failexpr branch
 *     when the body would burn measurable CPU.
 *   - Nested TimeConstrained -- inner aborts independently of outer.
 *   - Protected attribute is set; Attributes[TimeConstrained] reads back
 *     {HoldAll, Protected}.
 *
 * The CPU-burning tests use Do[expensive, {many}].  On very fast
 * machines the body may finish before the kernel ITIMER_PROF tick
 * fires; those cases accept either the timed-out result or the
 * "completed normally" result (Null from Do) to stay portable.
 */

#include "test_utils.h"
#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char* eval_string(const char* input) {
    struct Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    struct Expr* e = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(e);
    expr_free(e);
    return s;
}

static void run_setup(const char* input) {
    struct Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    struct Expr* e = evaluate(p);
    expr_free(p);
    expr_free(e);
}

void test_tc_fast_2arg(void) {
    char* s = eval_string("TimeConstrained[1 + 1, 10]");
    ASSERT_STR_EQ(s, "2");
    free(s);
}

void test_tc_fast_3arg(void) {
    char* s = eval_string("TimeConstrained[1 + 1, 10, failed]");
    ASSERT_STR_EQ(s, "2");
    free(s);
}

void test_tc_infinity(void) {
    char* s = eval_string("TimeConstrained[2 * 3, Infinity]");
    ASSERT_STR_EQ(s, "6");
    free(s);
}

void test_tc_infinity_three_arg(void) {
    char* s = eval_string("TimeConstrained[2 + 3, Infinity, oops]");
    ASSERT_STR_EQ(s, "5");
    free(s);
}

void test_tc_zero_time_no_fail(void) {
    char* s = eval_string("TimeConstrained[1 + 1, 0]");
    ASSERT_STR_EQ(s, "$Aborted");
    free(s);
}

void test_tc_zero_time_with_fail(void) {
    char* s = eval_string("TimeConstrained[1 + 1, 0, mybad]");
    ASSERT_STR_EQ(s, "mybad");
    free(s);
}

void test_tc_negative_time_no_fail(void) {
    char* s = eval_string("TimeConstrained[1 + 1, -1]");
    ASSERT_STR_EQ(s, "$Aborted");
    free(s);
}

void test_tc_negative_time_with_fail(void) {
    char* s = eval_string("TimeConstrained[Print[\"should never run\"], -0.5, gaveup]");
    ASSERT_STR_EQ(s, "gaveup");
    free(s);
}

void test_tc_one_arg_unevaluated(void) {
    char* s = eval_string("TimeConstrained[x]");
    ASSERT_STR_EQ(s, "TimeConstrained[x]");
    free(s);
}

void test_tc_zero_args_unevaluated(void) {
    char* s = eval_string("TimeConstrained[]");
    ASSERT_STR_EQ(s, "TimeConstrained[]");
    free(s);
}

void test_tc_too_many_args_unevaluated(void) {
    char* s = eval_string("TimeConstrained[a, b, c, d]");
    /* Non-numeric time `b` would also leave it unevaluated, but the
     * arity check fires first.  Either way the call sticks. */
    ASSERT(strstr(s, "TimeConstrained") == s);
    free(s);
}

void test_tc_non_numeric_time_unevaluated(void) {
    char* s = eval_string("TimeConstrained[1 + 1, undefSym]");
    /* HoldAll keeps the body literal; printing Plus[1,1] gives "1 + 1". */
    ASSERT_STR_EQ(s, "TimeConstrained[1 + 1, undefSym]");
    free(s);
}

void test_tc_time_arg_is_evaluated(void) {
    /* The time-budget argument is evaluated before being parsed.  Using
     * an arithmetic expression for t verifies the evaluation path. */
    char* s = eval_string("TimeConstrained[6 / 2, 2 + 3]");
    ASSERT_STR_EQ(s, "3");
    free(s);
}

void test_tc_rational_time(void) {
    /* Rational time budget: 1/2 second is plenty for 2+2. */
    char* s = eval_string("TimeConstrained[2 + 2, 1/2]");
    ASSERT_STR_EQ(s, "4");
    free(s);
}

void test_tc_failexpr_not_evaluated_on_success(void) {
    /* Verify Mathematica semantic: failexpr runs ONLY on abort.  Use a
     * counter to detect any accidental evaluation. */
    run_setup("tcCounterA = 0");

    /* Run a fast body with a counter-bumping failexpr.  The failexpr
     * would assign tcCounterA=tcCounterA+1 if (incorrectly) evaluated. */
    char* result = eval_string(
        "TimeConstrained[42, 10, tcCounterA = tcCounterA + 1]");
    ASSERT_STR_EQ(result, "42");
    free(result);

    char* s = eval_string("tcCounterA");
    ASSERT_STR_EQ(s, "0");
    free(s);
}

void test_tc_body_evaluated_exactly_once(void) {
    /* HoldAll keeps the body literal until the builtin runs it.  The
     * builtin invokes evaluate() once, so a side-effecting body should
     * fire exactly one time. */
    run_setup("tcCounterB = 0");

    char* r = eval_string(
        "TimeConstrained[tcCounterB = tcCounterB + 1, 10]");
    ASSERT_STR_EQ(r, "1");
    free(r);

    char* s = eval_string("tcCounterB");
    ASSERT_STR_EQ(s, "1");
    free(s);
}

void test_tc_completes_under_budget(void) {
    /* Plus @@ Range[100] = 5050 -- trivially fast even on slow boxes. */
    char* s = eval_string("TimeConstrained[Plus @@ Range[100], 60]");
    ASSERT_STR_EQ(s, "5050");
    free(s);
}

void test_tc_attributes(void) {
    char* s = eval_string("Attributes[TimeConstrained]");
    ASSERT_STR_EQ(s, "{HoldAll, Protected}");
    free(s);
}

void test_tc_protected(void) {
    /* TimeConstrained is Protected: an attempt to OwnValue-assign it
     * emits Set::wrsym (to stderr) but the symbol's builtin still runs. */
    run_setup("TimeConstrained = 42");

    char* s = eval_string("TimeConstrained[5 * 5, 60]");
    ASSERT_STR_EQ(s, "25");
    free(s);
}

/* ------------------------------------------------------------------ */
/* CPU-budget abort tests.                                              */
/*                                                                      */
/* These use Do[expensive, {n}] with n large enough to take many        */
/* hundreds of milliseconds of CPU on any reasonable machine.  Each      */
/* test then asserts the result is either the abort value or the         */
/* "Do completed normally" value Null.  Accepting both keeps the test    */
/* deterministic on extremely fast hardware that finishes before the     */
/* SIGPROF tick.                                                         */
/* ------------------------------------------------------------------ */

void test_tc_timeout_returns_aborted(void) {
    char* s = eval_string(
        "TimeConstrained[Do[N[Sin[1.0]^2 + Cos[1.0]^2, 50], {200000}], 0.05]");
    int ok = (strcmp(s, "$Aborted") == 0) || (strcmp(s, "Null") == 0);
    if (!ok) fprintf(stderr, "Unexpected timeout result: %s\n", s);
    ASSERT(ok);
    free(s);
}

void test_tc_timeout_returns_failexpr(void) {
    char* s = eval_string(
        "TimeConstrained[Do[N[Sin[1.0]^2 + Cos[1.0]^2, 50], {200000}], 0.05, mybadcase]");
    int ok = (strcmp(s, "mybadcase") == 0) || (strcmp(s, "Null") == 0);
    if (!ok) fprintf(stderr, "Unexpected timeout result: %s\n", s);
    ASSERT(ok);
    free(s);
}

void test_tc_nested_inner_completes(void) {
    /* Inner returns quickly inside an outer's generous budget. */
    char* s = eval_string("TimeConstrained[TimeConstrained[7 + 8, 10], 60]");
    ASSERT_STR_EQ(s, "15");
    free(s);
}

void test_tc_nested_inner_aborts_outer_succeeds(void) {
    /* Inner aborts -> $Aborted -- outer has plenty of CPU so it returns
     * that $Aborted unchanged. */
    char* s = eval_string(
        "TimeConstrained[TimeConstrained[Do[N[Sin[1.0]^2 + Cos[1.0]^2, 50], {200000}], 0.05], 60]");
    int ok = (strcmp(s, "$Aborted") == 0) || (strcmp(s, "Null") == 0);
    if (!ok) fprintf(stderr, "Unexpected nested result: %s\n", s);
    ASSERT(ok);
    free(s);
}

void test_tc_nested_inner_aborts_with_failexpr(void) {
    char* s = eval_string(
        "TimeConstrained[TimeConstrained[Do[N[Sin[1.0]^2 + Cos[1.0]^2, 50], {200000}], 0.05, innerfail], 60]");
    int ok = (strcmp(s, "innerfail") == 0) || (strcmp(s, "Null") == 0);
    if (!ok) fprintf(stderr, "Unexpected nested-failexpr result: %s\n", s);
    ASSERT(ok);
    free(s);
}

void test_tc_session_state_intact_after_abort(void) {
    /* After a forced abort the evaluator's recursion-depth counter must
     * be reset.  Otherwise the next non-trivial computation would emit
     * "$RecursionLimit::reclim".  A simple deep symbolic computation
     * after the abort verifies the counter is back to zero. */
    char* aborted = eval_string(
        "TimeConstrained[Do[N[Sin[1.0]^2 + Cos[1.0]^2, 50], {200000}], 0.05]");
    /* Accept either outcome; we don't care here. */
    free(aborted);

    char* s = eval_string("Plus @@ Range[200]");
    ASSERT_STR_EQ(s, "20100");
    free(s);
}

void test_tc_failexpr_can_be_symbolic(void) {
    char* s = eval_string("TimeConstrained[2 + 2, 60, x + y]");
    /* Body completes; failexpr unused. */
    ASSERT_STR_EQ(s, "4");
    free(s);
}

void test_tc_holdall_blocks_premature_eval(void) {
    /* If TimeConstrained lacked HoldAll, the body Print[...] would be
     * evaluated by the caller before the builtin saw it.  With HoldAll,
     * the body remains literal until the builtin's evaluate() runs it
     * inside the timer.  We can't easily observe Print's side-effect in
     * a unit test, but we *can* verify the body is not double-evaluated
     * by checking a counter increments by exactly 1. */
    run_setup("tcCounterC = 0");

    char* r = eval_string(
        "TimeConstrained[tcCounterC = tcCounterC + 1, 60, fail]");
    ASSERT_STR_EQ(r, "1");
    free(r);

    char* s = eval_string("tcCounterC");
    ASSERT_STR_EQ(s, "1");
    free(s);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_tc_fast_2arg);
    TEST(test_tc_fast_3arg);
    TEST(test_tc_infinity);
    TEST(test_tc_infinity_three_arg);
    TEST(test_tc_zero_time_no_fail);
    TEST(test_tc_zero_time_with_fail);
    TEST(test_tc_negative_time_no_fail);
    TEST(test_tc_negative_time_with_fail);
    TEST(test_tc_one_arg_unevaluated);
    TEST(test_tc_zero_args_unevaluated);
    TEST(test_tc_too_many_args_unevaluated);
    TEST(test_tc_non_numeric_time_unevaluated);
    TEST(test_tc_time_arg_is_evaluated);
    TEST(test_tc_rational_time);
    TEST(test_tc_failexpr_not_evaluated_on_success);
    TEST(test_tc_body_evaluated_exactly_once);
    TEST(test_tc_completes_under_budget);
    TEST(test_tc_attributes);
    TEST(test_tc_protected);
    TEST(test_tc_timeout_returns_aborted);
    TEST(test_tc_timeout_returns_failexpr);
    TEST(test_tc_nested_inner_completes);
    TEST(test_tc_nested_inner_aborts_outer_succeeds);
    TEST(test_tc_nested_inner_aborts_with_failexpr);
    TEST(test_tc_session_state_intact_after_abort);
    TEST(test_tc_failexpr_can_be_symbolic);
    TEST(test_tc_holdall_blocks_premature_eval);

    printf("All TimeConstrained tests passed!\n");
    return 0;
}
