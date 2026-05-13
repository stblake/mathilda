/*
 * test_repl_hooks.c
 *
 * Coverage for the REPL session hooks $PreRead, $Pre, $Post, $PrePrint
 * and $Epilog. Each test exercises both the unset (pass-through) and
 * set (transforming) path through the helpers in repl_hooks.{c,h},
 * plus a handful of compose tests that mirror the order in which the
 * REPL itself runs the hooks during process_input.
 */

#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "repl_hooks.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/* helpers                                                            */
/* ------------------------------------------------------------------ */

/* Evaluate `src` for side effects (e.g. assigning a hook). Discards
 * the result. */
static void run(const char* src) {
    Expr* parsed = parse_expression(src);
    ASSERT(parsed != NULL);
    Expr* result = evaluate(parsed);
    expr_free(parsed);
    expr_free(result);
}

/* Wipe every REPL hook so subsequent tests start from a known state.
 * Uses Clear[] to mimic the user-facing way to disarm a hook. */
static void clear_all_hooks(void) {
    run("Clear[$PreRead]");
    run("Clear[$Pre]");
    run("Clear[$Post]");
    run("Clear[$PrePrint]");
    run("Clear[$Epilog]");
}

/* Convenience: parse `src` and return the evaluated Expr. Caller must
 * expr_free the result. */
static Expr* parse_eval(const char* src) {
    Expr* parsed = parse_expression(src);
    ASSERT(parsed != NULL);
    Expr* result = evaluate(parsed);
    expr_free(parsed);
    return result;
}

/* Convenience: compare an Expr's printed form against a literal. */
static void assert_print_eq(Expr* e, const char* expected) {
    char* s = expr_to_string(e);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: expected '%s', got '%s'\n", expected, s);
    }
    ASSERT(strcmp(s, expected) == 0);
    free(s);
}

/* ------------------------------------------------------------------ */
/* init                                                               */
/* ------------------------------------------------------------------ */

static void test_hooks_init_creates_symbols(void) {
    clear_all_hooks();
    /* Each hook symbol should resolve via symtab without crashing,
     * even though no OwnValue is installed. */
    ASSERT(symtab_lookup("$PreRead") != NULL);
    ASSERT(symtab_lookup("$Pre") != NULL);
    ASSERT(symtab_lookup("$Post") != NULL);
    ASSERT(symtab_lookup("$PrePrint") != NULL);
    ASSERT(symtab_lookup("$Epilog") != NULL);

    /* No default OwnValue means hooks default to no-op. */
    ASSERT(symtab_get_own_values("$PreRead") == NULL);
    ASSERT(symtab_get_own_values("$Pre") == NULL);
    ASSERT(symtab_get_own_values("$Post") == NULL);
    ASSERT(symtab_get_own_values("$PrePrint") == NULL);
    ASSERT(symtab_get_own_values("$Epilog") == NULL);
}

static void test_hooks_have_docstrings(void) {
    clear_all_hooks();
    const char* doc;
    doc = symtab_get_docstring("$PreRead");
    ASSERT(doc != NULL && strstr(doc, "input") != NULL);
    doc = symtab_get_docstring("$Pre");
    ASSERT(doc != NULL && strstr(doc, "input") != NULL);
    doc = symtab_get_docstring("$Post");
    ASSERT(doc != NULL && strstr(doc, "output") != NULL);
    doc = symtab_get_docstring("$PrePrint");
    ASSERT(doc != NULL && strstr(doc, "printed") != NULL);
    doc = symtab_get_docstring("$Epilog");
    ASSERT(doc != NULL && strstr(doc, "terminat") != NULL);
}

/* ------------------------------------------------------------------ */
/* $PreRead                                                            */
/* ------------------------------------------------------------------ */

static void test_pre_read_unset_passes_through(void) {
    clear_all_hooks();
    char* out = repl_apply_pre_read("Plus[1, 2]");
    ASSERT(out != NULL);
    ASSERT_STR_EQ(out, "Plus[1, 2]");
    free(out);
}

static void test_pre_read_unset_handles_empty_input(void) {
    clear_all_hooks();
    char* out = repl_apply_pre_read("");
    ASSERT(out != NULL);
    ASSERT_STR_EQ(out, "");
    free(out);
}

static void test_pre_read_set_to_identity(void) {
    clear_all_hooks();
    /* Identity function leaves the string unchanged. */
    run("$PreRead = Identity");
    char* out = repl_apply_pre_read("foo bar");
    ASSERT(out != NULL);
    ASSERT_STR_EQ(out, "foo bar");
    free(out);
    clear_all_hooks();
}

static void test_pre_read_set_to_string_join(void) {
    clear_all_hooks();
    /* Wrap input string with a prefix/suffix using StringJoin. */
    run("$PreRead = Function[s, StringJoin[\"<<\", s, \">>\"]]");
    char* out = repl_apply_pre_read("hi");
    ASSERT(out != NULL);
    ASSERT_STR_EQ(out, "<<hi>>");
    free(out);
    clear_all_hooks();
}

static void test_pre_read_can_rewrite_to_valid_expression(void) {
    clear_all_hooks();
    /* The hook should be able to rewrite "x" -> "x + 1" so the parser
     * sees a different expression. We verify this by composing the
     * hook with parse_expression + evaluate, mirroring the REPL. */
    run("$PreRead = Function[s, StringJoin[s, \" + 1\"]]");
    char* cooked = repl_apply_pre_read("3");
    ASSERT(cooked != NULL);
    ASSERT_STR_EQ(cooked, "3 + 1");
    Expr* parsed = parse_expression(cooked);
    free(cooked);
    Expr* result = evaluate(parsed);
    expr_free(parsed);
    assert_print_eq(result, "4");
    expr_free(result);
    clear_all_hooks();
}

static void test_pre_read_non_string_return_falls_back(void) {
    clear_all_hooks();
    /* Pre-read that returns an integer instead of a string: the
     * helper should warn (to stderr) and fall back to the original
     * input. We can only test the fallback value here; the warning
     * path is exercised but not asserted on. */
    run("$PreRead = Function[s, 42]");
    /* Re-open stderr at /dev/null while the diagnostic fires so it
     * does not clutter the test transcript. freopen is portable C99. */
    fflush(stderr);
    FILE* saved_err = freopen("/dev/null", "w", stderr);
    char* out = repl_apply_pre_read("untouched");
    if (saved_err) {
        /* Restore. On macOS dup2-style restore needs a fresh fd; the
         * simplest portable thing is to reopen stderr at /dev/tty. If
         * /dev/tty is unavailable (CI), fall back to leaving it
         * pointed at /dev/null -- the rest of the test does not write
         * to stderr. */
        FILE* restored = freopen("/dev/tty", "w", stderr);
        (void)restored;
    }
    ASSERT(out != NULL);
    ASSERT_STR_EQ(out, "untouched");
    free(out);
    clear_all_hooks();
}

static void test_pre_read_clear_disables_hook(void) {
    clear_all_hooks();
    run("$PreRead = Function[s, StringJoin[\"X\", s]]");
    char* out = repl_apply_pre_read("a");
    ASSERT_STR_EQ(out, "Xa");
    free(out);

    run("Clear[$PreRead]");
    out = repl_apply_pre_read("a");
    ASSERT_STR_EQ(out, "a");
    free(out);
}

/* ------------------------------------------------------------------ */
/* $Pre                                                                */
/* ------------------------------------------------------------------ */

static void test_pre_unset_passes_through(void) {
    clear_all_hooks();
    Expr* in = parse_expression("Plus[1, 2]");
    Expr* out = repl_apply_pre(in);
    /* Unset hook returns the input pointer unchanged. */
    ASSERT(out == in);
    expr_free(out);
}

static void test_pre_set_to_identity(void) {
    clear_all_hooks();
    run("$Pre = Identity");
    Expr* in = parse_expression("Plus[1, 2]");
    Expr* out = repl_apply_pre(in);
    /* Identity[Plus[1,2]] evaluates to 3. */
    assert_print_eq(out, "3");
    expr_free(out);
    clear_all_hooks();
}

static void test_pre_doubles_input(void) {
    clear_all_hooks();
    run("$Pre = Function[x, 2 x]");
    Expr* in = parse_expression("5");
    Expr* out = repl_apply_pre(in);
    /* The wrapped value is the integer 5; $Pre[5] -> 10. */
    assert_print_eq(out, "10");
    expr_free(out);
    clear_all_hooks();
}

static void test_pre_evaluates_input_first_when_not_held(void) {
    clear_all_hooks();
    /* Without HoldAll, the wrapped expression Plus[1,2] is evaluated
     * to 3 BEFORE $Pre receives it; the hook then doubles 3 -> 6. */
    run("$Pre = Function[x, 2 x]");
    Expr* in = parse_expression("Plus[1, 2]");
    Expr* out = repl_apply_pre(in);
    assert_print_eq(out, "6");
    expr_free(out);
    clear_all_hooks();
}

static void test_pre_with_holdall_sees_unevaluated(void) {
    clear_all_hooks();
    /* A HoldAll-attributed pure function captures the input
     * unevaluated. Wrap with HoldForm so we can inspect it. */
    run("$Pre = Function[x, Length[Unevaluated[x]], HoldAll]");
    Expr* in = parse_expression("Plus[a, b, c]");
    Expr* out = repl_apply_pre(in);
    /* Length of the held Plus expression is 3. */
    assert_print_eq(out, "3");
    expr_free(out);
    clear_all_hooks();
}

static void test_pre_clear_disables_hook(void) {
    clear_all_hooks();
    run("$Pre = Function[x, x + 100]");
    Expr* in = parse_expression("1");
    Expr* out = repl_apply_pre(in);
    assert_print_eq(out, "101");
    expr_free(out);

    run("Clear[$Pre]");
    in = parse_expression("1");
    out = repl_apply_pre(in);
    /* Pass-through; the integer 1 survives unchanged. */
    assert_print_eq(out, "1");
    expr_free(out);
}

/* ------------------------------------------------------------------ */
/* $Post                                                               */
/* ------------------------------------------------------------------ */

static void test_post_unset_passes_through(void) {
    clear_all_hooks();
    Expr* in = parse_expression("42");
    Expr* out = repl_apply_post(in);
    ASSERT(out == in);
    expr_free(out);
}

static void test_post_wraps_result(void) {
    clear_all_hooks();
    run("$Post = Function[x, wrapper[x]]");
    Expr* in = parse_eval("1 + 2");
    Expr* out = repl_apply_post(in);
    assert_print_eq(out, "wrapper[3]");
    expr_free(out);
    clear_all_hooks();
}

static void test_post_chains_on_evaluator_output(void) {
    clear_all_hooks();
    run("$Post = Function[x, x + 10]");
    Expr* in = parse_eval("Plus[1, 2]"); /* yields 3 */
    Expr* out = repl_apply_post(in);
    assert_print_eq(out, "13");
    expr_free(out);
    clear_all_hooks();
}

/* ------------------------------------------------------------------ */
/* $PrePrint                                                           */
/* ------------------------------------------------------------------ */

static void test_pre_print_unset_passes_through(void) {
    clear_all_hooks();
    Expr* in = parse_expression("foo");
    Expr* out = repl_apply_pre_print(in);
    ASSERT(out == in);
    expr_free(out);
}

static void test_pre_print_can_rewrite_for_display(void) {
    clear_all_hooks();
    run("$PrePrint = Function[x, displayWrap[x]]");
    Expr* in = parse_eval("1 + 2");
    Expr* out = repl_apply_pre_print(in);
    assert_print_eq(out, "displayWrap[3]");
    expr_free(out);
    clear_all_hooks();
}

/* ------------------------------------------------------------------ */
/* $Epilog                                                             */
/* ------------------------------------------------------------------ */

static void test_epilog_unset_is_noop(void) {
    clear_all_hooks();
    /* Should not crash, no observable side effect. */
    repl_apply_epilog();
    /* Probe a sentinel symbol that the noop must NOT have touched. */
    Expr* probe = parse_eval("epilogTouched");
    assert_print_eq(probe, "epilogTouched");
    expr_free(probe);
}

static void test_epilog_evaluates_assigned_value(void) {
    clear_all_hooks();
    /* Use SetDelayed so the body fires each time $Epilog evaluates. */
    run("epilogRan = False");
    run("$Epilog := (epilogRan = True; doneSentinel)");

    /* Sanity check: the side effect has not yet occurred. */
    Expr* probe = parse_eval("epilogRan");
    assert_print_eq(probe, "False");
    expr_free(probe);

    repl_apply_epilog();

    probe = parse_eval("epilogRan");
    assert_print_eq(probe, "True");
    expr_free(probe);

    /* Also verify the result of the epilog is allowed to be any
     * expression -- we don't print it but we do free it. */
    clear_all_hooks();
    run("Clear[epilogRan]");
}

static void test_epilog_can_run_multiple_times(void) {
    clear_all_hooks();
    run("epilogCount = 0");
    run("$Epilog := (epilogCount = epilogCount + 1)");

    repl_apply_epilog();
    repl_apply_epilog();
    repl_apply_epilog();

    Expr* probe = parse_eval("epilogCount");
    assert_print_eq(probe, "3");
    expr_free(probe);

    clear_all_hooks();
    run("Clear[epilogCount]");
}

static void test_epilog_with_set_uses_static_value(void) {
    clear_all_hooks();
    /* Plain Set assigns the *evaluated* RHS once. The OwnValue then
     * just resolves to that value on each call; the side effect
     * (incrementing the counter) only fires at assignment time. */
    run("staticCount = 0");
    run("$Epilog = (staticCount = staticCount + 1)");

    repl_apply_epilog();
    repl_apply_epilog();
    repl_apply_epilog();

    Expr* probe = parse_eval("staticCount");
    /* Expected: counter incremented exactly once at the Set. */
    assert_print_eq(probe, "1");
    expr_free(probe);

    clear_all_hooks();
    run("Clear[staticCount]");
}

/* ------------------------------------------------------------------ */
/* compose: simulate the REPL's read -> $PreRead -> parse -> $Pre   */
/* -> evaluate -> $Post -> $PrePrint -> print pipeline                */
/* ------------------------------------------------------------------ */

static void test_pipeline_no_hooks_set(void) {
    clear_all_hooks();
    char* cooked = repl_apply_pre_read("1 + 2 + 3");
    Expr* parsed = parse_expression(cooked);
    free(cooked);
    Expr* pre = repl_apply_pre(expr_copy(parsed));
    Expr* evaluated = evaluate(pre);
    expr_free(pre);
    Expr* post = repl_apply_post(evaluated);
    Expr* to_print = repl_apply_pre_print(expr_copy(post));
    assert_print_eq(to_print, "6");
    expr_free(to_print);
    expr_free(parsed);
    expr_free(post);
}

static void test_pipeline_pre_post_compose(void) {
    clear_all_hooks();
    run("$Pre = Function[x, x + 1]");
    run("$Post = Function[x, x * 10]");

    char* cooked = repl_apply_pre_read("5");
    Expr* parsed = parse_expression(cooked);
    free(cooked);
    Expr* pre = repl_apply_pre(expr_copy(parsed));
    Expr* evaluated = evaluate(pre);
    expr_free(pre);
    Expr* post = repl_apply_post(evaluated);

    /* 5 -> $Pre -> 6 -> evaluate -> 6 -> $Post -> 60. */
    assert_print_eq(post, "60");

    expr_free(parsed);
    expr_free(post);
    clear_all_hooks();
}

static void test_pipeline_full_stack(void) {
    clear_all_hooks();
    /* PreRead doubles the digit string by appending another copy.
     * Pre adds 1, Post multiplies by 10, PrePrint wraps in printed[]. */
    run("$PreRead = Function[s, StringJoin[s, \" + \", s]]");
    run("$Pre = Function[x, x + 1]");
    run("$Post = Function[x, x * 10]");
    run("$PrePrint = Function[x, printed[x]]");

    char* cooked = repl_apply_pre_read("4");
    /* "4" -> "4 + 4". */
    ASSERT_STR_EQ(cooked, "4 + 4");
    Expr* parsed = parse_expression(cooked);
    free(cooked);
    /* parsed = Plus[4, 4]. */

    Expr* pre = repl_apply_pre(expr_copy(parsed));
    /* $Pre[Plus[4,4]] -> 8 + 1 = 9 (Plus evaluated first as 8). */
    Expr* evaluated = evaluate(pre);
    expr_free(pre);
    /* evaluated = 9 (already a fixed point). */

    Expr* post = repl_apply_post(evaluated);
    /* $Post[9] -> 90. */
    assert_print_eq(post, "90");

    Expr* to_print = repl_apply_pre_print(expr_copy(post));
    /* $PrePrint[90] -> printed[90]. */
    assert_print_eq(to_print, "printed[90]");

    expr_free(to_print);
    expr_free(parsed);
    expr_free(post);
    clear_all_hooks();
}

static void test_pipeline_pre_print_does_not_affect_post(void) {
    clear_all_hooks();
    run("$Post = Function[x, x + 1000]");
    run("$PrePrint = Function[x, x + 7]");

    Expr* parsed = parse_expression("5");
    Expr* pre = repl_apply_pre(expr_copy(parsed));
    Expr* evaluated = evaluate(pre);
    expr_free(pre);

    Expr* post = repl_apply_post(evaluated);
    /* $Post yields 1005. */
    assert_print_eq(post, "1005");

    /* $PrePrint[1005] yields 1012, but post itself stays at 1005 so
     * Out[n] would still be assigned 1005 by the REPL. */
    Expr* to_print = repl_apply_pre_print(expr_copy(post));
    assert_print_eq(to_print, "1012");
    assert_print_eq(post, "1005");

    expr_free(to_print);
    expr_free(parsed);
    expr_free(post);
    clear_all_hooks();
}

/* Set $Pre to a head with HoldAll so it inspects the un-evaluated
 * parsed input. Mirrors the documented use case of intercepting the
 * raw input expression. */
static void test_pre_holdall_intercepts_raw_input(void) {
    clear_all_hooks();
    run("SetAttributes[capture, HoldAll]");
    run("capture[x_] := holdResult[Length[Unevaluated[x]]]");
    run("$Pre = capture");

    Expr* parsed = parse_expression("Plus[a, b, c, d]");
    Expr* out = repl_apply_pre(parsed);
    /* capture sees Plus[a,b,c,d] unevaluated; Length is 4. */
    assert_print_eq(out, "holdResult[4]");
    expr_free(out);

    clear_all_hooks();
    run("Clear[capture]");
}

/* Reassigning a hook while it's set should swap the behaviour without
 * leaking or stacking the previous transformation. */
static void test_hook_reassignment_replaces_behaviour(void) {
    clear_all_hooks();
    run("$Post = Function[x, first[x]]");
    Expr* a = repl_apply_post(parse_eval("1"));
    assert_print_eq(a, "first[1]");
    expr_free(a);

    run("$Post = Function[x, second[x]]");
    Expr* b = repl_apply_post(parse_eval("1"));
    assert_print_eq(b, "second[1]");
    expr_free(b);

    clear_all_hooks();
}

/* The $PreRead text hook should not see any expression evaluation:
 * its input is the literal user-typed string, not a parsed Expr. */
static void test_pre_read_passes_raw_string(void) {
    clear_all_hooks();
    /* Capture the input string into an OwnValue we can inspect. */
    run("$PreRead = Function[s, (lastInput = s; s)]");
    char* out = repl_apply_pre_read("Plus[1, 2 + 3]");
    ASSERT_STR_EQ(out, "Plus[1, 2 + 3]");
    free(out);

    Expr* probe = parse_eval("lastInput");
    assert_print_eq(probe, "\"Plus[1, 2 + 3]\"");
    expr_free(probe);

    clear_all_hooks();
    run("Clear[lastInput]");
}

/* Reading a hook back via OwnValues should expose the user's
 * assignment, confirming the hook lives in the regular symbol table
 * rather than a parallel mechanism. */
static void test_hook_visible_via_ownvalues(void) {
    clear_all_hooks();
    run("$Pre = Function[x, x]");
    Rule* r = symtab_get_own_values("$Pre");
    ASSERT(r != NULL);
    /* OwnValues for $Pre returns at least one rule whose pattern is
     * the bare symbol $Pre. */
    ASSERT(r->pattern != NULL);
    ASSERT(r->pattern->type == EXPR_SYMBOL);
    ASSERT_STR_EQ(r->pattern->data.symbol, "$Pre");
    clear_all_hooks();
}

/* ------------------------------------------------------------------ */
/* main                                                                */
/* ------------------------------------------------------------------ */

int main(void) {
    symtab_init();
    core_init();

    TEST(test_hooks_init_creates_symbols);
    TEST(test_hooks_have_docstrings);

    TEST(test_pre_read_unset_passes_through);
    TEST(test_pre_read_unset_handles_empty_input);
    TEST(test_pre_read_set_to_identity);
    TEST(test_pre_read_set_to_string_join);
    TEST(test_pre_read_can_rewrite_to_valid_expression);
    TEST(test_pre_read_non_string_return_falls_back);
    TEST(test_pre_read_clear_disables_hook);

    TEST(test_pre_unset_passes_through);
    TEST(test_pre_set_to_identity);
    TEST(test_pre_doubles_input);
    TEST(test_pre_evaluates_input_first_when_not_held);
    TEST(test_pre_with_holdall_sees_unevaluated);
    TEST(test_pre_clear_disables_hook);

    TEST(test_post_unset_passes_through);
    TEST(test_post_wraps_result);
    TEST(test_post_chains_on_evaluator_output);

    TEST(test_pre_print_unset_passes_through);
    TEST(test_pre_print_can_rewrite_for_display);

    TEST(test_epilog_unset_is_noop);
    TEST(test_epilog_evaluates_assigned_value);
    TEST(test_epilog_can_run_multiple_times);
    TEST(test_epilog_with_set_uses_static_value);

    TEST(test_pipeline_no_hooks_set);
    TEST(test_pipeline_pre_post_compose);
    TEST(test_pipeline_full_stack);
    TEST(test_pipeline_pre_print_does_not_affect_post);

    TEST(test_pre_holdall_intercepts_raw_input);
    TEST(test_hook_reassignment_replaces_behaviour);
    TEST(test_pre_read_passes_raw_string);
    TEST(test_hook_visible_via_ownvalues);

    printf("All REPL hook tests passed.\n");
    return 0;
}
