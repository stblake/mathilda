#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Default OwnValue is seeded by eval_init() called from core_init(). */
static void test_default_visible(void) {
    assert_eval_eq("$RecursionLimit", "1024", 0);
}

static void test_set_updates_c_state(void) {
    assert_eval_eq("$RecursionLimit = 2048", "2048", 0);
    assert_eval_eq("$RecursionLimit", "2048", 0);
    ASSERT(eval_get_recursion_limit() == 2048);
}

static void test_setdelayed_evaluates_at_assign(void) {
    assert_eval_eq("$RecursionLimit := 200 + 800", "Null", 0);
    /* OwnValue retrieves and evaluates the body each time. */
    assert_eval_eq("$RecursionLimit", "1000", 0);
    ASSERT(eval_get_recursion_limit() == 1000);
}

static void test_below_minimum_is_rejected(void) {
    /* Get into a known state first. */
    eval_set_recursion_limit(512);
    Expr* sym  = expr_new_symbol("$RecursionLimit");
    Expr* curr = expr_new_integer(512);
    symtab_add_own_value("$RecursionLimit", sym, curr);
    expr_free(sym);
    expr_free(curr);

    /* Below MIN_RECURSION_LIMIT (20). The OwnValue must be rolled back
     * to the active C-side limit (512), not left at the rejected value. */
    assert_eval_eq("$RecursionLimit = 5", "5", 0); /* Set returns the RHS */
    assert_eval_eq("$RecursionLimit", "512", 0);    /* but storage rolled back */
    ASSERT(eval_get_recursion_limit() == 512);
}

static void test_guard_fires_at_user_limit(void) {
    /* A self-referential rule blows the recursion. With a low limit we get
     * a Hold[]-wrapped sub-expression back instead of running indefinitely. */
    assert_eval_eq("$RecursionLimit = 50", "50", 0);
    assert_eval_eq("Clear[g]", "Null", 0);
    assert_eval_eq("g[x_] := g[x] + 1", "Null", 0);
    /* It terminates with the recursive CALL held: the guard stops the
     * g[1] that would re-enter the evaluator, so the result is
     * k + Hold[g[1]].  A valueless symbol never re-enters, so the head g on
     * its own is never what gets held (Hold[g][1] was an artifact of the
     * head being evaluated through the guarded path).  The count of
     * accumulated +1 terms depends on the precise unwind path and is not
     * pinned. */
    Expr* parsed = parse_expression("g[1]");
    Expr* evaluated = evaluate(parsed);
    expr_free(parsed);
    char* str = expr_to_string(evaluated);
    ASSERT_MSG(strstr(str, "Hold[g[1]]") != NULL,
               "expected Hold[g[1]] in overflow output, got: %s", str);
    ASSERT_MSG(strstr(str, "Hold[g]") == NULL,
               "a bare head must not be held, got: %s", str);
    free(str);
    expr_free(evaluated);
}

/* Under Trace the valueless-symbol fast path declines (every evaluation gets
 * a trace frame), so the head g reaches evaluate()'s recursion guard itself.
 * The guard must give the same answer the fast path does -- the symbol, not
 * Hold[g] -- or the overflow's shape would depend on whether a trace is
 * active: the traced run must hold the call g[1], never the bare head. */
static void test_guard_agrees_under_trace(void) {
    assert_eval_eq("$RecursionLimit = 50", "50", 0);
    assert_eval_eq("Clear[g]", "Null", 0);
    assert_eval_eq("g[x_] := g[x] + 1", "Null", 0);
    assert_eval_eq("t = ToString[InputForm[Trace[g[1]]]]; "
                   "{StringContainsQ[t, \"Hold[g[1]]\"], StringContainsQ[t, \"Hold[g][\"]}",
                   "{True, False}", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_default_visible);
    TEST(test_set_updates_c_state);
    TEST(test_setdelayed_evaluates_at_assign);
    TEST(test_below_minimum_is_rejected);
    TEST(test_guard_fires_at_user_limit);
    TEST(test_guard_agrees_under_trace);

    printf("All $RecursionLimit tests passed!\n");
    symtab_clear();
    return 0;
}
