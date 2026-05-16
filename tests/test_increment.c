#include "core.h"
#include "expr.h"
#include "symtab.h"
#include "eval.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
 * Tests for Increment, Decrement, PreIncrement, PreDecrement, AddTo,
 * SubtractFrom and their parser aliases ++, --, +=, -=.
 *
 * Each REPL transcript from the feature request is translated into a test.
 * Because Mathilda tests run in a shared symbol table, each test clears the
 * symbols it touches at entry to make the assertions independent.
 */

static void reset(const char* name) {
    char buf[128];
    snprintf(buf, sizeof(buf), "Clear[%s]", name);
    Expr* p = parse_expression(buf);
    if (p) { Expr* r = evaluate(p); expr_free(p); expr_free(r); }
}

/* ---------- Parser smoke tests: all four tokens now parse ---------- */

static void test_parser_accepts_new_tokens(void) {
    reset("x");
    assert_eval_eq("x = 10", "10", 0);
    assert_eval_eq("x++", "10", 0);
    assert_eval_eq("x", "11", 0);
    assert_eval_eq("x--", "11", 0);
    assert_eval_eq("x", "10", 0);
    assert_eval_eq("x += 1", "11", 0);
    assert_eval_eq("x", "11", 0);
    assert_eval_eq("x -= 1", "10", 0);
    assert_eval_eq("x", "10", 0);
    reset("x");
}

/* ---------- Increment ---------- */

static void test_increment_returns_old_value(void) {
    reset("k");
    assert_eval_eq("k = 1", "1", 0);
    assert_eval_eq("k++", "1", 0);
    assert_eval_eq("k", "2", 0);
    reset("k");
}

static void test_increment_numerical_value(void) {
    reset("x");
    assert_eval_eq("x = 1.5", "1.5", 0);
    assert_eval_eq("x++", "1.5", 0);
    assert_eval_eq("x", "2.5", 0);
    reset("x");
}

static void test_increment_symbolic_value(void) {
    reset("v");
    reset("a");
    assert_eval_eq("v = a", "a", 0);
    assert_eval_eq("v++", "a", 0);
    assert_eval_eq("v", "1 + a", 0);
    reset("v");
}

static void test_increment_list(void) {
    reset("x");
    assert_eval_eq("x = {1, 2, 3}", "{1, 2, 3}", 0);
    assert_eval_eq("x++", "{1, 2, 3}", 0);
    assert_eval_eq("x", "{2, 3, 4}", 0);
    reset("x");
}

static void test_increment_part(void) {
    reset("list");
    assert_eval_eq("list = {1, 2, 3}", "{1, 2, 3}", 0);
    assert_eval_eq("list[[2]]++", "2", 0);
    assert_eval_eq("list", "{1, 3, 3}", 0);
    reset("list");
}

static void test_increment_rvalue(void) {
    reset("undefSym");
    /* Argument has no value: emits rvalue message, stays unevaluated. */
    assert_eval_eq("undefSym++", "Increment[undefSym]", 0);
    reset("undefSym");
}

/* ---------- PreIncrement ---------- */

static void test_preincrement_returns_new_value(void) {
    reset("k");
    assert_eval_eq("k = 1", "1", 0);
    assert_eval_eq("++k", "2", 0);
    assert_eval_eq("k", "2", 0);
    reset("k");
}

static void test_preincrement_numerical_value(void) {
    reset("x");
    assert_eval_eq("x = 1.5", "1.5", 0);
    assert_eval_eq("++x", "2.5", 0);
    assert_eval_eq("x", "2.5", 0);
    reset("x");
}

static void test_preincrement_symbolic_value(void) {
    reset("v");
    reset("a");
    assert_eval_eq("v = a", "a", 0);
    assert_eval_eq("++v", "1 + a", 0);
    assert_eval_eq("v", "1 + a", 0);
    reset("v");
}

static void test_preincrement_list(void) {
    reset("x");
    assert_eval_eq("x = {1, 2, 3}", "{1, 2, 3}", 0);
    assert_eval_eq("++x", "{2, 3, 4}", 0);
    assert_eval_eq("x", "{2, 3, 4}", 0);
    reset("x");
}

static void test_preincrement_part(void) {
    reset("list");
    assert_eval_eq("list = {1, 2, 3}", "{1, 2, 3}", 0);
    assert_eval_eq("++list[[2]]", "3", 0);
    assert_eval_eq("list", "{1, 3, 3}", 0);
    reset("list");
}

static void test_preincrement_rvalue(void) {
    reset("undefSym");
    assert_eval_eq("++undefSym", "PreIncrement[undefSym]", 0);
    reset("undefSym");
}

/* ---------- Decrement ---------- */

static void test_decrement_returns_old_value(void) {
    reset("k");
    assert_eval_eq("k = 1", "1", 0);
    assert_eval_eq("k--", "1", 0);
    assert_eval_eq("k", "0", 0);
    reset("k");
}

static void test_decrement_numerical_value(void) {
    reset("x");
    assert_eval_eq("x = 1.5", "1.5", 0);
    assert_eval_eq("x--", "1.5", 0);
    assert_eval_eq("x", "0.5", 0);
    reset("x");
}

static void test_decrement_symbolic_value(void) {
    reset("v");
    reset("a");
    assert_eval_eq("v = a", "a", 0);
    assert_eval_eq("v--", "a", 0);
    assert_eval_eq("v", "-1 + a", 0);
    reset("v");
}

static void test_decrement_list(void) {
    reset("x");
    assert_eval_eq("x = {1, 2, 3}", "{1, 2, 3}", 0);
    assert_eval_eq("x--", "{1, 2, 3}", 0);
    assert_eval_eq("x", "{0, 1, 2}", 0);
    reset("x");
}

static void test_decrement_part(void) {
    reset("list");
    assert_eval_eq("list = {1, 2, 3}", "{1, 2, 3}", 0);
    assert_eval_eq("list[[2]]--", "2", 0);
    assert_eval_eq("list", "{1, 1, 3}", 0);
    reset("list");
}

static void test_decrement_rvalue(void) {
    reset("undefSym");
    assert_eval_eq("undefSym--", "Decrement[undefSym]", 0);
    reset("undefSym");
}

/* ---------- PreDecrement ---------- */

static void test_predecrement_returns_new_value(void) {
    reset("k");
    assert_eval_eq("k = 1", "1", 0);
    assert_eval_eq("--k", "0", 0);
    assert_eval_eq("k", "0", 0);
    reset("k");
}

static void test_predecrement_numerical_value(void) {
    reset("x");
    assert_eval_eq("x = 1.5", "1.5", 0);
    assert_eval_eq("--x", "0.5", 0);
    assert_eval_eq("x", "0.5", 0);
    reset("x");
}

static void test_predecrement_symbolic_value(void) {
    reset("v");
    reset("a");
    assert_eval_eq("v = a", "a", 0);
    assert_eval_eq("--v", "-1 + a", 0);
    assert_eval_eq("v", "-1 + a", 0);
    reset("v");
}

static void test_predecrement_list(void) {
    reset("x");
    assert_eval_eq("x = {1, 2, 3}", "{1, 2, 3}", 0);
    assert_eval_eq("--x", "{0, 1, 2}", 0);
    assert_eval_eq("x", "{0, 1, 2}", 0);
    reset("x");
}

static void test_predecrement_part(void) {
    reset("list");
    assert_eval_eq("list = {1, 2, 3}", "{1, 2, 3}", 0);
    assert_eval_eq("--list[[2]]", "1", 0);
    assert_eval_eq("list", "{1, 1, 3}", 0);
    reset("list");
}

static void test_predecrement_rvalue(void) {
    reset("undefSym");
    assert_eval_eq("--undefSym", "PreDecrement[undefSym]", 0);
    reset("undefSym");
}

/* ---------- AddTo ---------- */

static void test_addto_returns_new_value(void) {
    reset("k");
    assert_eval_eq("k = 1", "1", 0);
    assert_eval_eq("k += 5", "6", 0);
    assert_eval_eq("k", "6", 0);
    reset("k");
}

static void test_addto_numerical_value(void) {
    reset("x");
    assert_eval_eq("x = 1.5", "1.5", 0);
    assert_eval_eq("x += 3.75", "5.25", 0);
    assert_eval_eq("x", "5.25", 0);
    reset("x");
}

static void test_addto_symbolic_value(void) {
    reset("v");
    reset("a");
    reset("b");
    assert_eval_eq("v = a", "a", 0);
    assert_eval_eq("v += b", "a + b", 0);
    assert_eval_eq("v", "a + b", 0);
    reset("v");
}

static void test_addto_list_scalar(void) {
    reset("x");
    assert_eval_eq("x = {1, 2, 3}", "{1, 2, 3}", 0);
    assert_eval_eq("x += 17", "{18, 19, 20}", 0);
    assert_eval_eq("x", "{18, 19, 20}", 0);
    reset("x");
}

static void test_addto_list_list(void) {
    reset("x");
    assert_eval_eq("x = {18, 19, 20}", "{18, 19, 20}", 0);
    assert_eval_eq("x += {20, 21, 22}", "{38, 40, 42}", 0);
    assert_eval_eq("x", "{38, 40, 42}", 0);
    reset("x");
}

static void test_addto_part(void) {
    reset("list");
    reset("c");
    assert_eval_eq("list = {1, 2, 3}", "{1, 2, 3}", 0);
    assert_eval_eq("list[[2]] += c", "2 + c", 0);
    assert_eval_eq("list", "{1, 2 + c, 3}", 0);
    reset("list");
}

static void test_addto_rvalue(void) {
    reset("undefSym");
    assert_eval_eq("undefSym += 2", "AddTo[undefSym, 2]", 0);
    reset("undefSym");
}

/* ---------- SubtractFrom ---------- */

static void test_subtractfrom_returns_new_value(void) {
    reset("k");
    assert_eval_eq("k = 1", "1", 0);
    assert_eval_eq("k -= 5", "-4", 0);
    assert_eval_eq("k", "-4", 0);
    reset("k");
}

static void test_subtractfrom_numerical_value(void) {
    reset("x");
    assert_eval_eq("x = 1.5", "1.5", 0);
    assert_eval_eq("x -= 0.75", "0.75", 0);
    assert_eval_eq("x", "0.75", 0);
    reset("x");
}

static void test_subtractfrom_symbolic_value(void) {
    reset("v");
    reset("a");
    reset("b");
    assert_eval_eq("v = a", "a", 0);
    assert_eval_eq("v -= b", "a - b", 0);
    assert_eval_eq("v", "a - b", 0);
    reset("v");
}

static void test_subtractfrom_list_scalar(void) {
    reset("x");
    assert_eval_eq("x = {100, 200, 300}", "{100, 200, 300}", 0);
    assert_eval_eq("x -= 17", "{83, 183, 283}", 0);
    assert_eval_eq("x", "{83, 183, 283}", 0);
    reset("x");
}

static void test_subtractfrom_list_list(void) {
    reset("x");
    assert_eval_eq("x = {83, 183, 283}", "{83, 183, 283}", 0);
    assert_eval_eq("x -= {20, 21, 22}", "{63, 162, 261}", 0);
    assert_eval_eq("x", "{63, 162, 261}", 0);
    reset("x");
}

static void test_subtractfrom_part(void) {
    reset("list");
    reset("c");
    assert_eval_eq("list = {1, 2, 3}", "{1, 2, 3}", 0);
    assert_eval_eq("list[[2]] -= c", "2 - c", 0);
    assert_eval_eq("list", "{1, 2 - c, 3}", 0);
    reset("list");
}

static void test_subtractfrom_rvalue(void) {
    reset("undefSym");
    assert_eval_eq("undefSym -= 2", "SubtractFrom[undefSym, 2]", 0);
    reset("undefSym");
}

/* ---------- Parser edge cases ---------- */

static void test_postfix_and_binary_coexist(void) {
    /* a++ b is Increment[a] * b (implicit multiplication after postfix ++) */
    reset("a");
    reset("b");
    assert_eval_eq("a = 10; b = 3", "3", 0);
    assert_eval_eq("a++ b", "30", 0);
    assert_eval_eq("a", "11", 0);
    reset("a");
    reset("b");
}

static void test_prefix_binds_tighter_than_plus(void) {
    /* ++a + 1 is (++a) + 1 = new_a + 1 */
    reset("a");
    assert_eval_eq("a = 10", "10", 0);
    assert_eval_eq("++a + 1", "12", 0);
    assert_eval_eq("a", "11", 0);
    reset("a");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_parser_accepts_new_tokens);

    TEST(test_increment_returns_old_value);
    TEST(test_increment_numerical_value);
    TEST(test_increment_symbolic_value);
    TEST(test_increment_list);
    TEST(test_increment_part);
    TEST(test_increment_rvalue);

    TEST(test_preincrement_returns_new_value);
    TEST(test_preincrement_numerical_value);
    TEST(test_preincrement_symbolic_value);
    TEST(test_preincrement_list);
    TEST(test_preincrement_part);
    TEST(test_preincrement_rvalue);

    TEST(test_decrement_returns_old_value);
    TEST(test_decrement_numerical_value);
    TEST(test_decrement_symbolic_value);
    TEST(test_decrement_list);
    TEST(test_decrement_part);
    TEST(test_decrement_rvalue);

    TEST(test_predecrement_returns_new_value);
    TEST(test_predecrement_numerical_value);
    TEST(test_predecrement_symbolic_value);
    TEST(test_predecrement_list);
    TEST(test_predecrement_part);
    TEST(test_predecrement_rvalue);

    TEST(test_addto_returns_new_value);
    TEST(test_addto_numerical_value);
    TEST(test_addto_symbolic_value);
    TEST(test_addto_list_scalar);
    TEST(test_addto_list_list);
    TEST(test_addto_part);
    TEST(test_addto_rvalue);

    TEST(test_subtractfrom_returns_new_value);
    TEST(test_subtractfrom_numerical_value);
    TEST(test_subtractfrom_symbolic_value);
    TEST(test_subtractfrom_list_scalar);
    TEST(test_subtractfrom_list_list);
    TEST(test_subtractfrom_part);
    TEST(test_subtractfrom_rvalue);

    TEST(test_postfix_and_binary_coexist);
    TEST(test_prefix_binds_tighter_than_plus);

    printf("All increment/decrement/addto/subtractfrom tests passed!\n");
    return 0;
}
