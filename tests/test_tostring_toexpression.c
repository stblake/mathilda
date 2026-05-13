/*
 * Tests for ToString and ToExpression builtins.
 *
 * ToString[expr]            -> string in InputForm
 * ToString[expr, form]      -> string in InputForm / FullForm / TeXForm
 * ToExpression[input]       -> parse and evaluate
 * ToExpression[in, form]    -> parse with form, evaluate
 * ToExpression[in, form, h] -> wrap parsed expression with head h before
 *                              evaluation (e.g. h = Hold)
 *
 * Each test runs the call through evaluate() so we exercise the same
 * code path the REPL uses, then compares the printed Out[...] form.
 */

#include "core.h"
#include "expr.h"
#include "symtab.h"
#include "eval.h"
#include "test_utils.h"
#include "print.h"
#include "parse.h"
#include <stdio.h>
#include <string.h>

void test_tostring_default(void) {
    assert_eval_eq("ToString[x^2 + y^3]", "\"x^2 + y^3\"", 0);
    assert_eval_eq("ToString[1 + 1]", "\"2\"", 0);
    assert_eval_eq("ToString[123]", "\"123\"", 0);
    /* InputForm of a String includes the surrounding quotes, so the
     * stored content is the 7-byte string `"hello"`. When the test
     * framework prints that EXPR_STRING for comparison it adds another
     * pair of display quotes, giving "" + hello + "" = 9 visible chars. */
    assert_eval_eq("ToString[\"hello\"]", "\"\"hello\"\"", 0);
}

void test_tostring_inputform(void) {
    assert_eval_eq("ToString[x^2 + y^3, InputForm]", "\"x^2 + y^3\"", 0);
    assert_eval_eq("ToString[Sin[Pi/4], InputForm]", "\"1/Sqrt[2]\"", 0);
    assert_eval_eq("ToString[1/2, InputForm]", "\"1/2\"", 0);
}

void test_tostring_fullform(void) {
    assert_eval_eq("ToString[x^2 + y^3, FullForm]",
                   "\"Plus[Power[x, 2], Power[y, 3]]\"", 0);
    assert_eval_eq("ToString[1 + 1, FullForm]", "\"2\"", 0);
    assert_eval_eq("ToString[{a, b, c}, FullForm]",
                   "\"List[a, b, c]\"", 0);
}

void test_tostring_texform(void) {
    /* TeXForm of x^2+y^3 should produce a TeX-shaped string. We test
     * a known render to lock the contract; if the TeX printer changes
     * its precise output the test will need to follow. */
    assert_eval_eq("ToString[x^2 + y^3, TeXForm]", "\"x^{2}+y^{3}\"", 0);
    assert_eval_eq("ToString[1/2, TeXForm]", "\"\\frac{1}{2}\"", 0);
}

void test_tostring_invalid_form(void) {
    /* Unsupported forms should leave the call unevaluated rather than
     * silently fall through to InputForm. */
    assert_eval_eq("ToString[x, MadeUpForm]", "ToString[x, MadeUpForm]", 0);
}

void test_toexpression_default(void) {
    assert_eval_eq("ToExpression[\"1+1\"]", "2", 0);
    assert_eval_eq("ToExpression[\"x + y\"]", "x + y", 0);
    assert_eval_eq("ToExpression[\"{1, 2, 3}\"]", "{1, 2, 3}", 0);
    assert_eval_eq("ToExpression[\"Sin[Pi/2]\"]", "1", 0);
}

void test_toexpression_form(void) {
    /* InputForm and FullForm both currently route through parse_expression. */
    assert_eval_eq("ToExpression[\"2+3\", InputForm]", "5", 0);
    assert_eval_eq("ToExpression[\"Plus[2, 3]\", FullForm]", "5", 0);
}

void test_toexpression_with_hold(void) {
    /* The 3-arg form wraps the parsed expression with head h before the
     * outer evaluator re-enters it. Hold prevents the inner evaluation. */
    assert_eval_eq("ToExpression[\"1+1\", InputForm, Hold]", "Hold[1 + 1]", 0);
    assert_eval_eq("ToExpression[\"x*y\", InputForm, Hold]", "Hold[x y]", 0);
}

void test_toexpression_parse_failure(void) {
    /* parse_expression returns NULL on syntax error; we return $Failed. */
    assert_eval_eq("ToExpression[\"1+\"]", "$Failed", 0);
    assert_eval_eq("ToExpression[\"]\"]", "$Failed", 0);
}

void test_toexpression_non_string(void) {
    /* Non-string input leaves the call unevaluated. */
    assert_eval_eq("ToExpression[123]", "ToExpression[123]", 0);
}

void test_toexpression_listable(void) {
    /* ToExpression is Listable: threads over List arguments. */
    assert_eval_eq("ToExpression[{\"1+1\", \"2+2\", \"3+3\"}]",
                   "{2, 4, 6}", 0);
}

void test_roundtrip(void) {
    /* Compositional contract: ToExpression[ToString[expr]] reproduces
     * expr for the small set of expressions where the printed form is
     * canonical. */
    assert_eval_eq("ToExpression[ToString[x^2 + y^3]]", "x^2 + y^3", 0);
    assert_eval_eq("ToExpression[ToString[{1, 2, 3}]]", "{1, 2, 3}", 0);
    assert_eval_eq("ToExpression[ToString[Sin[Pi/4]]]", "1/Sqrt[2]", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_tostring_default);
    TEST(test_tostring_inputform);
    TEST(test_tostring_fullform);
    TEST(test_tostring_texform);
    TEST(test_tostring_invalid_form);
    TEST(test_toexpression_default);
    TEST(test_toexpression_form);
    TEST(test_toexpression_with_hold);
    TEST(test_toexpression_parse_failure);
    TEST(test_toexpression_non_string);
    TEST(test_toexpression_listable);
    TEST(test_roundtrip);

    printf("All ToString/ToExpression tests passed!\n");
    return 0;
}
