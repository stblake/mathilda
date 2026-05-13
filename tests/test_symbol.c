/*
 * Tests for the Symbol[] builtin in core.c.
 *
 * Symbol["name"] refers to a symbol with the specified name. It validates
 * the name, resolves it through the context system, and returns an
 * EXPR_SYMBOL.
 */

#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "core.h"
#include "symtab.h"
#include "context.h"
#include <stdio.h>
#include <string.h>

/* ---- Basic construction ---- */

static void test_symbol_bare_returns_symbol(void) {
    /* Symbol["foo"] in Global` context: legacy bare-name behavior. */
    assert_eval_eq("Symbol[\"foo\"]", "foo", 0);
}

static void test_symbol_head_is_symbol(void) {
    /* The constructed value has Head Symbol. */
    assert_eval_eq("Head[Symbol[\"xyzzy\"]]", "Symbol", 0);
}

static void test_head_of_plain_symbol(void) {
    /* Every symbol -- not just ones built via Symbol[] -- has head Symbol. */
    assert_eval_eq("Head[zorp]", "Symbol", 0);
}

/* ---- Pattern matching: x_Symbol ---- */

static void test_symbol_pattern_match(void) {
    /* The canonical example from the Mathematica docs. */
    assert_eval_eq(
        "{f[xx], f[\"xx\"], f[2]} /. f[s_Symbol] :> g[s]",
        "{g[xx], f[\"xx\"], f[2]}",
        0);
}

static void test_match_q_on_symbol(void) {
    assert_eval_eq("MatchQ[abc, _Symbol]", "True", 0);
    assert_eval_eq("MatchQ[\"abc\", _Symbol]", "False", 0);
    assert_eval_eq("MatchQ[42, _Symbol]", "False", 0);
}

/* ---- Context handling ---- */

static void test_symbol_qualified_absolute(void) {
    /* "a`x" is an absolutely-qualified name; should round-trip verbatim. */
    assert_eval_eq("Symbol[\"a`x\"]", "a`x", 0);
}

static void test_symbol_relative_in_global(void) {
    /* In Global`, a leading backtick yields a name relative to Global`.
     * Because Global` is the current context the short form just prints "x". */
    assert_eval_eq("Symbol[\"`relx\"]", "relx", 0);
}

static void test_symbol_relative_inside_package(void) {
    /* Symbol["`name"] inside Begin["Pkg$Test`"] should produce Pkg$Test`name.
     * Inside that context, Pkg$Test` matches $Context so the short form
     * "rel" is what gets printed. After End[], Pkg$Test` is no longer the
     * current context (nor on $ContextPath), so the same symbol now prints
     * with its fully-qualified name. */
    assert_eval_eq("Begin[\"Pkg$Test`\"]", "\"Pkg$Test`\"", 0);
    assert_eval_eq("Symbol[\"`rel\"]", "rel", 0);
    assert_eval_eq("End[]", "\"Pkg$Test`\"", 0);
    /* The earlier Symbol[] call has interned Pkg$Test`rel; build it again
     * from outside and confirm it now prints fully qualified. */
    assert_eval_eq("Symbol[\"Pkg$Test`rel\"]", "Pkg$Test`rel", 0);
}

/* ---- Symbol assignment & retrieval ---- */

static void test_symbol_value_retrieval(void) {
    /* Setting a value on a symbol and then fetching it via Symbol[]. */
    assert_eval_eq("varQ = 99", "99", 0);
    /* Symbol["varQ"] returns the symbol varQ, which evaluates to 99. */
    assert_eval_eq("Symbol[\"varQ\"]", "99", 0);
    /* Clear it so this test does not pollute subsequent tests. */
    assert_eval_eq("Clear[varQ]", "Null", 0);
}

/* ---- Invalid names: must be returned unevaluated, error printed ---- */

static void test_symbol_invalid_starts_with_digit(void) {
    /* "1x" is invalid: a name cannot start with a digit. The call must
     * return the original expression unevaluated. */
    Expr* parsed = parse_expression("Symbol[\"1x\"]");
    ASSERT(parsed != NULL);
    /* Silence the symname error message on stderr for clean test output. */
    FILE* saved_stderr = stderr;
    (void)saved_stderr;
    freopen("/dev/null", "w", stderr);
    Expr* result = evaluate(parsed);
    /* Reopen stderr to the terminal. */
    freopen("/dev/tty", "w", stderr);

    expr_free(parsed);
    char* s = expr_to_string(result);
    ASSERT_STR_EQ(s, "Symbol[\"1x\"]");
    free(s);
    expr_free(result);
}

static void test_symbol_invalid_empty_string(void) {
    Expr* parsed = parse_expression("Symbol[\"\"]");
    ASSERT(parsed != NULL);
    freopen("/dev/null", "w", stderr);
    Expr* result = evaluate(parsed);
    freopen("/dev/tty", "w", stderr);

    expr_free(parsed);
    char* s = expr_to_string(result);
    ASSERT_STR_EQ(s, "Symbol[\"\"]");
    free(s);
    expr_free(result);
}

static void test_symbol_invalid_punctuation(void) {
    Expr* parsed = parse_expression("Symbol[\"a-b\"]");
    ASSERT(parsed != NULL);
    freopen("/dev/null", "w", stderr);
    Expr* result = evaluate(parsed);
    freopen("/dev/tty", "w", stderr);

    expr_free(parsed);
    char* s = expr_to_string(result);
    ASSERT_STR_EQ(s, "Symbol[\"a-b\"]");
    free(s);
    expr_free(result);
}

/* ---- Non-string argument: leave unevaluated ---- */

static void test_symbol_non_string_arg(void) {
    /* Symbol[123] does not satisfy the input type and must return the
     * call unevaluated. */
    assert_eval_eq("Symbol[123]", "Symbol[123]", 0);
}

static void test_symbol_wrong_arity(void) {
    /* 0-arg and 2-arg forms are not defined; remain unevaluated. */
    assert_eval_eq("Symbol[]", "Symbol[]", 0);
    assert_eval_eq("Symbol[\"a\", \"b\"]", "Symbol[\"a\", \"b\"]", 0);
}

/* ---- Names with digits and $ ---- */

static void test_symbol_with_digit_after_letter(void) {
    assert_eval_eq("Symbol[\"x1\"]", "x1", 0);
    assert_eval_eq("Symbol[\"x123y\"]", "x123y", 0);
}

static void test_symbol_with_dollar(void) {
    /* '$' is a valid identifier character (start or middle). */
    assert_eval_eq("Symbol[\"$foo\"]", "$foo", 0);
    assert_eval_eq("Symbol[\"a$b\"]", "a$b", 0);
}

/* ---- Attributes ---- */

static void test_symbol_is_protected(void) {
    /* Verify Symbol carries the Protected attribute. */
    assert_eval_eq("MemberQ[Attributes[Symbol], Protected]", "True", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_symbol_bare_returns_symbol);
    TEST(test_symbol_head_is_symbol);
    TEST(test_head_of_plain_symbol);

    TEST(test_symbol_pattern_match);
    TEST(test_match_q_on_symbol);

    TEST(test_symbol_qualified_absolute);
    TEST(test_symbol_relative_in_global);
    TEST(test_symbol_relative_inside_package);

    TEST(test_symbol_value_retrieval);

    TEST(test_symbol_invalid_starts_with_digit);
    TEST(test_symbol_invalid_empty_string);
    TEST(test_symbol_invalid_punctuation);

    TEST(test_symbol_non_string_arg);
    TEST(test_symbol_wrong_arity);

    TEST(test_symbol_with_digit_after_letter);
    TEST(test_symbol_with_dollar);

    TEST(test_symbol_is_protected);

    printf("All Symbol tests passed.\n");
    return 0;
}
