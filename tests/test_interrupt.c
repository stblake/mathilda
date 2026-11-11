/*
 * test_interrupt.c -- the user-interrupt flag (src/term/term_interrupt.h) as
 * the evaluator sees it: a pending interrupt becomes the ordinary Abort[]
 * sentinel at the next fixed-point step, CheckAbort catches it, AbortProtect
 * defers it, TimeConstrained does not mistake it for a timeout, Pause and the
 * heavy C kernels (through tc_deadline_passed) stop early, messages are muted
 * while it unwinds, and the flag never outlives the evaluation it aborted.
 *
 * TestInterrupt[] is a test-only builtin that raises the interrupt from inside
 * an evaluation, the way the REPL's SIGINT handler does from outside one.
 */
#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include "eval.h"
#include "sym_names.h"
#include "term/term_interrupt.h"

#include <stdio.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

bool tc_deadline_passed(void);

static Expr* builtin_test_interrupt(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 0) return NULL;
    mth_interrupt_request();
    return expr_new_symbol(SYM_Null);
}

/* TestDeadline[] -- tc_deadline_passed() as seen from inside an evaluation. */
static Expr* builtin_test_deadline(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 0) return NULL;
    return expr_new_symbol(tc_deadline_passed() ? SYM_True : SYM_False);
}

/* TestChildInterrupted[] -- what a subprocess builtin does after system()
 * reports a child killed by SIGINT: it requests an interrupt and returns. */
static Expr* builtin_test_child_interrupted(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 0) return NULL;
    volatile int child_killed_by_sigint = 1;  /* stands in for WIFSIGNALED && WTERMSIG == SIGINT */
    if (child_killed_by_sigint) mth_interrupt_request();
    return expr_new_integer(0);
}

static double now_s(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

/* Evaluate `input` as a top-level input; check the printed result and that
 * the interrupt did not leak past it. Returns the elapsed seconds. */
static double run_top(const char* input, const char* expected) {
    (void)eval_take_top_abort();
    double t0 = now_s();
    assert_eval_eq(input, expected, 0);
    double dt = now_s() - t0;
    ASSERT_MSG(!mth_interrupt_pending(), "interrupt still pending after %s", input);
    assert_eval_eq("1 + 1", "2", 0);   /* case 12: the flag never leaks */
    return dt;
}

static void test_loop(void) {
    /* 1. Raised inside, a 10^9-step loop ends at once, and the top level saw
     * the abort. */
    (void)eval_take_top_abort();
    double t0 = now_s();
    assert_eval_eq("TestInterrupt[]; Do[i, {i, 10^9}]", "$Aborted", 0);
    double dt = now_s() - t0;
    ASSERT_MSG(dt < 1.0, "Do took %.2fs", dt);
    ASSERT(eval_take_top_abort());
    ASSERT(!mth_interrupt_pending());
    assert_eval_eq("1 + 1", "2", 0);

    /* 2. Preset before evaluation, with an inert body that only ever takes the
     * stamped fixed-point fast path. */
    mth_interrupt_request();
    dt = run_top("Do[f[i], {i, 10^9}]", "$Aborted");
    ASSERT_MSG(dt < 1.0, "inert Do took %.2fs", dt);
}

static void test_checkabort_abortprotect(void) {
    /* 3. CheckAbort catches it and clears the flag. */
    run_top("CheckAbort[TestInterrupt[]; While[True], caught]", "caught");

    /* 4. AbortProtect defers it to the end of its body. */
    run_top("x = 0; AbortProtect[TestInterrupt[]; Do[x++, {1000}]]", "$Aborted");
    assert_eval_eq("x", "1000", 0);

    /* Nested: it takes effect at the outermost exit only. */
    run_top("y = 0; AbortProtect[AbortProtect[TestInterrupt[]; y = 1]; y = 2]", "$Aborted");
    assert_eval_eq("y", "2", 0);

    /* A CheckAbort around the AbortProtect catches it on the way out. */
    run_top("CheckAbort[AbortProtect[TestInterrupt[]; 1], caught]", "caught");
}

static void test_timeconstrained(void) {
    /* 5. An interrupt is not a timeout. */
    run_top("TimeConstrained[TestInterrupt[]; While[True], 100, fail]", "$Aborted");
}

static void test_pause(void) {
    /* 6. Pause ends early. */
    run_top("t = First[AbsoluteTiming[CheckAbort[TestInterrupt[]; Pause[5], c]]]; t < 0.5", "True");

    /* 7. ... but not inside AbortProtect. */
    double dt = run_top("AbortProtect[TestInterrupt[]; Pause[0.2]; 7]", "$Aborted");
    ASSERT_MSG(dt >= 0.2, "AbortProtect'd Pause cut short: %.3fs", dt);
}

static void test_deadline_poll(void) {
    /* 8. tc_deadline_passed() is the kernels' poll. */
    ASSERT(!tc_deadline_passed());
    mth_interrupt_request();
    ASSERT(tc_deadline_passed());
    mth_interrupt_clear();
    ASSERT(!tc_deadline_passed());
    /* Inside AbortProtect it stays false (checked from inside a builtin). */
    run_top("AbortProtect[TestInterrupt[]; r = TestDeadline[]]; r", "$Aborted");
    assert_eval_eq("r", "False", 0);
    run_top("CheckAbort[TestInterrupt[]; TestDeadline[], c]", "c");
}

static void test_callbacks_and_messages(void) {
    /* 9. Sticky through a callback whose value a builtin consumes. */
    run_top("Select[{1, 2, 3}, (TestInterrupt[]; True) &]", "$Aborted");

    /* 10. Messages are muted while it unwinds. */
    run_top("Quiet[TestInterrupt[]; 1/0; Do[i, {10^9}]]", "$Aborted");
    assert_eval_eq("$MessageList", "{}", 0);
}

static void test_heavy_kernel(void) {
    /* 11. A C kernel that polls tc_deadline_passed() declines promptly. */
    double dt = run_top("TestInterrupt[]; FactorInteger[NextPrime[10^30] NextPrime[10^31]]",
                        "$Aborted");
    ASSERT_MSG(dt < 2.0, "FactorInteger took %.2fs", dt);
}

static void test_public_contract(void) {
    /* 13. The public API from plain C, with no evaluation running ... */
    ASSERT(!mth_interrupt_pending());
    mth_interrupt_request();
    ASSERT(mth_interrupt_pending());
    double dt = run_top("Do[i, {i, 10^9}]", "$Aborted");
    ASSERT_MSG(dt < 1.0, "Do took %.2fs", dt);
    ASSERT(!mth_interrupt_pending());

    /* ... and from a builtin that saw its child die of SIGINT. */
    run_top("r = TestChildInterrupted[]; Do[i, {i, 10^9}]; r", "$Aborted");
    run_top("CheckAbort[TestChildInterrupted[]; Pause[5], interrupted]", "interrupted");
}

int main(void) {
    symtab_init();
    core_init();
    symtab_add_builtin("TestInterrupt", builtin_test_interrupt);
    symtab_add_builtin("TestDeadline", builtin_test_deadline);
    symtab_add_builtin("TestChildInterrupted", builtin_test_child_interrupted);
    TEST(test_loop);
    TEST(test_checkabort_abortprotect);
    TEST(test_timeconstrained);
    TEST(test_pause);
    TEST(test_deadline_poll);
    TEST(test_callbacks_and_messages);
    TEST(test_heavy_kernel);
    TEST(test_public_contract);
    printf("All interrupt tests passed.\n");
    return 0;
}
