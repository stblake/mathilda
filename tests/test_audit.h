/*
 * test_audit.h -- shared harness for the builtin-audit regression suites.
 *
 * The 2026-11-14 builtin audit (review/report.html, review/issues/<Head>.md)
 * confirmed 2,475 defects.  Every finding that is fixed becomes one check in
 * tests/test_audit_<category>.c, keyed by its audit id "<Head>#<n>" so a
 * failure names the write-up that explains it.  tools/audit_backlog.py tracks
 * which findings are still open.
 *
 * Expected values come from the mathematics (identities, SameQ / PossibleZeroQ
 * against a closed form, numeric agreement at a sample point) or from the
 * reference's documented message contract -- never from another system's
 * output.  Messages are asserted in the language itself:
 *     Check[expr, "fired", Head::tag]      -> "fired"    (routed, Check sees it)
 *     CheckAbort[expr, "aborted"]          -> "aborted"  (an Abort[] sentinel)
 * so each check also proves the diagnostic goes through mth_message (the
 * Quiet / Check surface), not a raw stderr write.
 *
 * Unlike assert_eval_eq, a failed check does not abort: every failure in the
 * suite is printed, and audit_finish() makes the exit status nonzero.  A crash
 * still kills the suite, which is the point for the crash-class findings.
 */
#ifndef MATHILDA_TEST_AUDIT_H
#define MATHILDA_TEST_AUDIT_H

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include "message.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_audit_fail = 0;
static int g_audit_pass = 0;

/* Evaluate `input`; its InputForm-style printed form must equal `expected`.
 * `id` is the audit finding ("Head#n"). */
static void audit_chk(const char* id, const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    if (!e) {
        fprintf(stderr, "FAIL [%s] does not parse: %s\n", id, input);
        g_audit_fail++;
        return;
    }
    Expr* r = evaluate(e);
    expr_free(e);
    char* got = expr_to_string(r);
    if (!got || strcmp(got, expected) != 0) {
        fprintf(stderr, "FAIL [%s] %s\n  expected: %s\n  got:      %s\n",
                id, input, expected, got ? got : "(null)");
        g_audit_fail++;
    } else {
        g_audit_pass++;
    }
    free(got);
    expr_free(r);
}

/* `input` must evaluate to True (for checks phrased as a predicate:
 * SameQ, PossibleZeroQ[a - b], Abs[a - b] < tol, ...). */
static void audit_true(const char* id, const char* input) {
    audit_chk(id, input, "True");
}

/* Evaluate `input` with every message muted (for side-effect statements
 * that set up a later check). */
static void audit_run(const char* input) {
    Expr* e = parse_expression(input);
    if (!e) return;
    mth_msg_suppress_push();
    Expr* r = evaluate(e);
    mth_msg_suppress_pop();
    expr_free(e);
    expr_free(r);
}

static int audit_finish(const char* suite) {
    printf("%s: %d passed, %d failed\n", suite, g_audit_pass, g_audit_fail);
    return g_audit_fail ? 1 : 0;
}

/* Common prologue for each suite's main(). */
static void audit_init(void) {
    symtab_init();
    core_init();
    test_load_init_m();
}

#endif /* MATHILDA_TEST_AUDIT_H */
