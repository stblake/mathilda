/* Unit tests for the Trace collector (beads-planning-b3k, Phase 2).
 *
 * These exercise the evaluator-level entry point eval_collect_trace() directly
 * -- the user-facing Trace[] builtin lands in Phase 3 (src/trace.c). The
 * collector evaluates a BORROWED held expression to a fixed point and returns a
 * freshly-owned flat List of the top-level intermediate expressions:
 *   - >=1 top-level rewrite -> {e0, e1, ..., eN}
 *   - no top-level rewrite (inert atom / normal form) -> {} (empty List)
 * The returned expression always has head List.
 *
 * Run under valgrind to confirm the ownership contract: each stored step is
 * inc-ref'd (expr_copy) and the discarded final value is freed exactly once.
 */
#include "expr.h"
#include "eval.h"
#include "print.h"
#include "symtab.h"
#include "core.h"
#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Trace a freshly-parsed (unevaluated) expression. The parsed input is
 * borrowed by eval_collect_trace, so we own and free it here. Returns the
 * collector's List result (caller frees). */
static Expr* trace_of(const char* input) {
    Expr* parsed = parse_expression(input);
    ASSERT(parsed != NULL);
    Expr* trace = eval_collect_trace(parsed);
    expr_free(parsed);
    ASSERT(trace != NULL);
    return trace;
}

/* Assert the trace result is a List head with exactly `n` elements. */
static void assert_list_len(Expr* trace, size_t n) {
    ASSERT(trace->type == EXPR_FUNCTION);
    ASSERT(trace->data.function.head->type == EXPR_SYMBOL);
    ASSERT(strcmp(trace->data.function.head->data.symbol.name, "List") == 0);
    ASSERT_MSG(trace->data.function.arg_count == n,
               "expected %zu steps, got %zu", n, trace->data.function.arg_count);
}

/* 1 + 1 takes one top-level rewrite: {1 + 1, 2}. */
static void test_arithmetic_step(void) {
    Expr* t = trace_of("1 + 1");
    assert_list_len(t, 2);
    char* s = expr_to_string(t);
    ASSERT_STR_EQ(s, "{1 + 1, 2}");
    free(s);
    /* Last element is the final result 2. */
    Expr* last = t->data.function.args[1];
    ASSERT(last->type == EXPR_INTEGER && last->data.integer == 2);
    expr_free(t);
}

/* An inert atom triggers no top-level rewrite -> empty List. */
static void test_inert_atom(void) {
    Expr* t = trace_of("5");
    assert_list_len(t, 0);
    char* s = expr_to_string(t);
    ASSERT_STR_EQ(s, "{}");   /* Wolfram Trace[5] -> {} */
    free(s);
    /* FullForm of the empty result confirms the List head is present. */
    char* ff = expr_to_string_fullform(t);
    ASSERT_STR_EQ(ff, "List[]");
    free(ff);
    expr_free(t);
}

/* A bare undefined symbol is already in normal form -> empty List. */
static void test_inert_symbol(void) {
    Expr* t = trace_of("zzz");
    assert_list_len(t, 0);
    expr_free(t);
}

/* A multi-step reduction records the held input first and the final value
 * last, with no duplicate trailing entry. The *exact* intermediate set for
 * x + 1 is validated in Phase 3 (builtin_trace); here we assert only the
 * invariants the collector guarantees regardless of rewrite granularity. */
static void test_multistep(void) {
    assert_eval_eq("x = 3", "3", 0);
    Expr* t = trace_of("x + 1");
    /* Nested trace of Plus[x, 1] with x = 3: the x -> 3 argument
     * sub-evaluation nests as {x, 3}, then the reassembled 3 + 1, then 4:
     * {{x, 3}, 3 + 1, 4}. */
    ASSERT(t->type == EXPR_FUNCTION);
    ASSERT(t->data.function.arg_count >= 2);
    /* First element is the nested x -> 3 sub-evaluation, not the held input. */
    char* first = expr_to_string(t->data.function.args[0]);
    ASSERT_STR_EQ(first, "{x, 3}");
    free(first);
    /* Last element is the final result 4, and it appears exactly once at the
     * end (the fixed point is not double-recorded). */
    Expr* last = t->data.function.args[t->data.function.arg_count - 1];
    ASSERT(last->type == EXPR_INTEGER && last->data.integer == 4);
    Expr* penult = t->data.function.args[t->data.function.arg_count - 2];
    ASSERT(!expr_eq(penult, last));
    expr_free(t);
    assert_eval_eq("Clear[x]", "Null", 0);
}

/* Reentrancy: eval_collect_trace called back-to-back leaves no residual global
 * collector state, and each call produces an independent, correct trace. This
 * also guards the Phase-1 clock mitigation -- a second trace of the same input
 * must not be emptied by the evaluation-cache early-exit. */
static void test_reentrancy_and_repeat(void) {
    Expr* a = trace_of("2 + 3");
    assert_list_len(a, 2);
    char* sa = expr_to_string(a);
    ASSERT_STR_EQ(sa, "{2 + 3, 5}");
    free(sa);
    expr_free(a);

    /* Repeat the identical trace: the clock bump inside eval_collect_trace must
     * defeat the timestamp early-exit so we still see both steps. */
    Expr* b = trace_of("2 + 3");
    assert_list_len(b, 2);
    char* sb = expr_to_string(b);
    ASSERT_STR_EQ(sb, "{2 + 3, 5}");
    free(sb);
    expr_free(b);

    /* An inert trace between two live ones must still be empty (no leakage of
     * the previous collector's seeded state). */
    Expr* c = trace_of("7");
    assert_list_len(c, 0);
    expr_free(c);
}

/* The user-facing Trace[] builtin (Phase 3, src/trace.c). Unlike the direct
 * eval_collect_trace tests above, these go through the evaluator, so they also
 * exercise the HoldForm wrapping that keeps the returned list inert on the
 * evaluator's re-evaluation pass (without it, {1 + 1, 2} would collapse to
 * {2, 2}). */
static void test_builtin_basic(void) {
    assert_eval_eq("Trace[1 + 1]", "{1 + 1, 2}", 0);
    assert_eval_eq("Trace[5]", "{}", 0);              /* inert atom */
    assert_eval_eq("Trace[zzz]", "{}", 0);            /* normal-form symbol */
    /* HoldAll: the argument reaches the builtin unevaluated. */
    assert_eval_eq("Attributes[Trace]", "{HoldAll, Protected}", 0);
}

static void test_builtin_holdform_idempotent(void) {
    /* Re-tracing the same input in one session must still yield the full
     * trace (guards the clock mitigation); and the held first element must not
     * collapse under re-evaluation (guards the HoldForm wrapping). */
    assert_eval_eq("Trace[1 + 1]", "{1 + 1, 2}", 0);
    assert_eval_eq("Trace[1 + 1]", "{1 + 1, 2}", 0);
}

/* Nested traces: argument sub-evaluations appear as sublists and the
 * reassembled intermediate form appears as a step (matching Mathematica). */
static void test_builtin_nested(void) {
    /* Each Power argument reduces first (nested), then the reassembled sum,
     * then Plus fires. */
    assert_eval_eq("Trace[2^3 + 4^2 + 1]",
                   "{{2^3, 8}, {4^2, 16}, 8 + 16 + 1, 25}", 0);
    /* Range[10] reduces (nested); reassembled x^{...}; Listable result. The
     * builtin-internal Range folding and the threaded element evaluations are
     * suppressed, so each shows as a single rewrite. */
    assert_eval_eq("Trace[x^Range[10]]",
                   "{{Range[10], {1, 2, 3, 4, 5, 6, 7, 8, 9, 10}}, "
                   "x^{1, 2, 3, 4, 5, 6, 7, 8, 9, 10}, "
                   "{x, x^2, x^3, x^4, x^5, x^6, x^7, x^8, x^9, x^10}}", 0);
    /* Times over a reduced sum. */
    assert_eval_eq("Trace[2 (3 + 4)]", "{{3 + 4, 7}, 2 7, 14}", 0);
    /* Argument changes but the head has no rule: reassembled f[2] recorded,
     * result f[2] deduped away. */
    assert_eval_eq("Trace[f[1 + 1]]", "{{1 + 1, 2}, f[2]}", 0);
    /* A DownValue firing is a traced step, and its RHS re-evaluates at the
     * outer level. */
    assert_eval_eq("gg[y_] := y^2", "Null", 0);
    assert_eval_eq("Trace[gg[1 + 1]]", "{{1 + 1, 2}, gg[2], 2^2, 4}", 0);
    assert_eval_eq("Clear[gg]", "Null", 0);
    /* A no-sub-evaluation expression is unchanged from the flat form. */
    assert_eval_eq("Trace[1 + 1]", "{1 + 1, 2}", 0);
}

/* The two-argument form Trace[expr, form] (h4u): the trace filtered to the
 * step leaves whose expression matches the pattern form, flattening away
 * nesting. HoldAll delivers form to the builtin unevaluated, so pattern
 * literals like _Integer work directly. */
static void test_builtin_two_arg_form(void) {
    /* _Integer keeps only the integer results, dropping the unevaluated Plus. */
    assert_eval_eq("Trace[1 + 2 + 3, _Integer]", "{6}", 0);
    /* 2*3 is now a (nested) sub-evaluation, so its result 6 is a matched
     * integer step too; matches are flattened out of the nesting. */
    assert_eval_eq("Trace[2*3 + 1, _Integer]", "{6, 7}", 0);
    /* Blank matches every step -> identical to the one-arg form. */
    assert_eval_eq("Trace[1 + 2 + 3, _]", "{1 + 2 + 3, 6}", 0);
    /* A form nothing matches -> {} (there are steps, but none is a Real). */
    assert_eval_eq("Trace[1 + 2 + 3, _Real]", "{}", 0);
    /* A bare symbol form matches only that literal symbol, not the steps. */
    assert_eval_eq("Trace[1 + 1, x]", "{}", 0);
    /* Head pattern: _f keeps the f[...] result of a nested application. */
    assert_eval_eq("Trace[Nest[f, x, 3], _f]", "{f[f[f[x]]]}", 0);
    /* No steps at all (normal form) stays {} whatever the form is. */
    assert_eval_eq("Trace[zzz, _]", "{}", 0);
    /* Arities other than 1/2 stay unevaluated (returns NULL, HoldAll holds). */
    assert_eval_eq("Trace[1 + 1, x, y]", "Trace[1 + 1, x, y]", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_arithmetic_step);
    TEST(test_inert_atom);
    TEST(test_inert_symbol);
    TEST(test_multistep);
    TEST(test_reentrancy_and_repeat);
    TEST(test_builtin_basic);
    TEST(test_builtin_holdform_idempotent);
    TEST(test_builtin_nested);
    TEST(test_builtin_two_arg_form);

    printf("All trace tests passed!\n");
    symtab_clear();
    return 0;
}
