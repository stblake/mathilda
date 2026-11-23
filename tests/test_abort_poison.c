/*
 * test_abort_poison.c -- an aborted computation must not leave the Abort[]
 * sentinel behind as DATA (EVAL-abort-1).
 *
 * The bug: an abort escaped every TimeConstrained and aborted a LATER
 * computation that had no time limit.  A timeout siglongjmp'd out of a Simplify
 * inside ParallelMixedSpecial and skipped its factor_memo_pop, so that
 * Simplify's Factor/Together/Cancel memo stayed active for the rest of the
 * session.  It held entries made while the timeout was in flight -- when the
 * evaluator answers every call with Abort[], so a Together whose inner
 * evaluation aborted had stored a value with Abort[] coefficients.  A later,
 * unconstrained Integrate read one back through PolynomialExtendedGCD,
 * evaluated the Abort[] and aborted with no constraint left to own it.
 *
 * Checked here:
 *   - the memo stack depth is restored by a timeout jump, so a jumped-out
 *     memo is never consulted again (test_factor_memo_unwind);
 *   - no memo entry is stored while an abort is due (test_memo_store_gate);
 *   - the two evaluator paths that turned a sentinel into a datum: a head
 *     that evaluates to Abort[] / Throw (Abort[]["E"]), and a Module / With /
 *     Block initializer that aborts (the body ran with the local bound to
 *     Abort[]).
 */
#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include "eval.h"
#include "memacct.h"
#include "mtime.h"
#include "expr.h"
#include "sym_names.h"
#include "poly/facpoly.h"

#include <stdio.h>
#include <string.h>

static Expr* mk_abort(void) {
    return expr_new_function(expr_new_symbol(SYM_Abort), NULL, 0);
}

/* TestMemoSpin[] -- what a Simplify unwound by a timeout looks like to the
 * memo stack: push a memo holding a poisoned Together entry (the shape a
 * Together interrupted by the timeout used to store), then run a pure C loop
 * that never returns to the evaluator, so the timeout is taken by the
 * SIGPROF backstop's siglongjmp -- past the pop that would follow. */
static Expr* builtin_test_memo_spin(Expr* res) {
    if (res->data.function.arg_count != 0) return NULL;
    FactorMemo* m = factor_memo_new();
    factor_memo_push(m);
    Expr* key = parse_expression("Together[1/(xq + 1) + 1/(xq - 1)]");
    Expr* val = expr_new_function(expr_new_symbol(SYM_Plus),
                                  (Expr*[]){ mk_abort(), expr_new_symbol("xq") }, 2);
    factor_memo_store(m, key, val);
    expr_free(key);
    expr_free(val);
    double t0 = mth_mono_seconds();
    volatile unsigned long spin = 0;
    for (;;) {                     /* bounded only as a safety net: 20 s */
        if ((++spin & 0xFFFFFFul) == 0 && mth_mono_seconds() - t0 > 20.0) break;
    }
    factor_memo_pop();
    factor_memo_free(m);
    return expr_new_symbol("spunOut");
}

static void test_factor_memo_unwind(void) {
    symtab_add_builtin("TestMemoSpin", builtin_test_memo_spin);
    ASSERT(factor_memo_active() == NULL);
    assert_eval_eq("TimeConstrained[TestMemoSpin[], 0.2, \"to\"]", "\"to\"", 0);
    /* The jump put the stack back: nothing is active outside a Simplify. */
    ASSERT(factor_memo_active() == NULL);
    /* The poisoned entry is unreachable: the unconstrained computation gets
     * its answer, no abort escapes. */
    assert_eval_eq("CheckAbort[Together[1/(xq + 1) + 1/(xq - 1)], \"ESCAPED\"]",
                   "(2 xq)/((-1 + xq) (1 + xq))", 0);
    /* Twice in a row, and nested inside an outer budget. */
    assert_eval_eq("TimeConstrained[TestMemoSpin[], 0.2, \"to\"]", "\"to\"", 0);
    assert_eval_eq("TimeConstrained[{TimeConstrained[TestMemoSpin[], 0.2, \"in\"], 1}, 10, \"out\"]",
                   "{\"in\", 1}", 0);
    ASSERT(factor_memo_active() == NULL);
    /* A Simplify afterwards pushes into the abandoned slot (freeing that memo)
     * and pops cleanly. */
    assert_eval_eq("Simplify[(xq^2 - 1)/(xq - 1)]", "1 + xq", 0);
    ASSERT(factor_memo_active() == NULL);
    assert_eval_eq("CheckAbort[Together[1/(xq + 1) + 1/(xq - 1)], \"ESCAPED\"]",
                   "(2 xq)/((-1 + xq) (1 + xq))", 0);
}

/* While a TIME frame is tripped (an abort is due), a memo store is skipped:
 * whatever the caller computed may carry the sentinel. */
static void test_memo_store_gate(void) {
    Expr* key = parse_expression("Factor[yq^2 - 1]");
    Expr* val = parse_expression("(yq - 1) (yq + 1)");
    FactorMemo* m = factor_memo_new();

    ASSERT(!eval_abort_raised());
    double now = mth_mono_seconds();
    int me = cf_time_push(now - 1.0);        /* a deadline already past */
    ASSERT(me >= 0);
    cf_time_note(now);                       /* ...and recorded: tripped */
    ASSERT(eval_abort_raised());
    factor_memo_store(m, key, val);
    ASSERT(factor_memo_lookup(m, key) == NULL);
    cf_unwind_to(me);                        /* the frame is gone */
    ASSERT(!eval_abort_raised());
    factor_memo_store(m, key, val);
    ASSERT(factor_memo_lookup(m, key) != NULL);

    factor_memo_free(m);
    expr_free(key);
    expr_free(val);
}

/* An evaluated head that is an in-flight sentinel unwinds the call, as an
 * evaluated argument does; it is never kept as the head of a datum. */
static void test_sentinel_in_head(void) {
    assert_eval_eq("CheckAbort[(Abort[])[\"F\"], \"ab\"]", "\"ab\"", 0);
    assert_eval_eq("CheckAbort[{(Abort[])[\"E\"], 2}, \"ab\"]", "\"ab\"", 0);
    assert_eval_eq("CheckAbort[(Abort[])[1][2], \"ab\"]", "\"ab\"", 0);
    assert_eval_eq("ph1[x_] := Abort[]; CheckAbort[ph1[1][2], \"ab\"]", "\"ab\"", 0);
    assert_eval_eq("Catch[(Throw[7])[x]]", "7", 0);
    assert_eval_eq("Catch[{1, (Throw[8, t])[y]}, t]", "8", 0);
    /* An association read through a symbol bound to an aborted value. */
    assert_eval_eq("CheckAbort[tq = <|\"E\" -> 1|>; tq = Abort[]; tq[\"E\"], \"ab\"]", "\"ab\"", 0);
    /* Held, Abort[][x] is an ordinary expression. */
    assert_eval_eq("Hold[(Abort[])[x]]", "Hold[Abort[][x]]", 0);
    assert_eval_eq("f9[g9][h9]", "f9[g9][h9]", 0);
}

/* An initializer that aborts leaves the scoping construct at once: the body
 * does not run with the local bound to Abort[]. */
static void test_initializer_abort(void) {
    assert_eval_eq("ran = {}; {CheckAbort[Module[{a = Abort[]}, AppendTo[ran, 1]; a], \"ab\"], ran}",
                   "{\"ab\", {}}", 0);
    assert_eval_eq("ran = {}; {CheckAbort[With[{a = Abort[]}, AppendTo[ran, 1]; a[\"F\"]], \"ab\"], ran}",
                   "{\"ab\", {}}", 0);
    assert_eval_eq("ran = {}; {CheckAbort[Block[{a = Abort[]}, AppendTo[ran, 1]; a], \"ab\"], ran}",
                   "{\"ab\", {}}", 0);
    /* Block put back what it had localised before the abort. */
    assert_eval_eq("bq = 5; {CheckAbort[Block[{bq = 1, cq = Abort[]}, bq], \"ab\"], bq}",
                   "{\"ab\", 5}", 0);
    /* Normal initializers are untouched. */
    assert_eval_eq("Module[{a = 2}, a^2]", "4", 0);
    assert_eval_eq("With[{a = 3}, a + 1]", "4", 0);
    assert_eval_eq("Catch[Module[{a = Throw[5]}, 1]]", "5", 0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_memo_store_gate);
    TEST(test_sentinel_in_head);
    TEST(test_initializer_abort);
    TEST(test_factor_memo_unwind);
    printf("All abort-poison tests passed.\n");
    return 0;
}
