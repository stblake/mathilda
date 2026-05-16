#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void run_full(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    char* s = expr_to_string_fullform(r);
    ASSERT_MSG(strcmp(s, expected) == 0,
               "FoldList %s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

static void run_infix(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    char* s = expr_to_string(r);
    ASSERT_MSG(strcmp(s, expected) == 0,
               "FoldList %s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

/* ---------- Spec examples ---------- */

static void test_foldlist_basic_symbolic(void) {
    /* FoldList[f, x, {a, b, c, d}]
       -> {x, f[x, a], f[f[x, a], b], f[f[f[x, a], b], c], f[f[f[f[x, a], b], c], d]} */
    run_full("FoldList[f, x, {a, b, c, d}]",
             "List[x, f[x, a], f[f[x, a], b], f[f[f[x, a], b], c], f[f[f[f[x, a], b], c], d]]");
}

static void test_foldlist_no_seed(void) {
    /* FoldList[f, {a, b, c, d}]
       -> {a, f[a, b], f[f[a, b], c], f[f[f[a, b], c], d]} */
    run_full("FoldList[f, {a, b, c, d}]",
             "List[a, f[a, b], f[f[a, b], c], f[f[f[a, b], c], d]]");
}

static void test_foldlist_head_preserved(void) {
    /* FoldList[f, x, p[a, b, c, d]]
       -> p[x, f[x, a], f[f[x, a], b], f[f[f[x, a], b], c], f[f[f[f[x, a], b], c], d]] */
    run_full("FoldList[f, x, p[a, b, c, d]]",
             "p[x, f[x, a], f[f[x, a], b], f[f[f[x, a], b], c], f[f[f[f[x, a], b], c], d]]");
}

static void test_foldlist_fold_right(void) {
    /* FoldList[g[#2, #1] &, x, {a, b, c, d}]
       -> {x, g[a, x], g[b, g[a, x]], g[c, g[b, g[a, x]]], g[d, g[c, g[b, g[a, x]]]]} */
    run_full("FoldList[g[#2, #1] &, x, {a, b, c, d}]",
             "List[x, g[a, x], g[b, g[a, x]], g[c, g[b, g[a, x]]], g[d, g[c, g[b, g[a, x]]]]]");
}

static void test_foldlist_factorials(void) {
    /* FoldList[Times, 1, Range[10]]
       -> {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800} */
    run_full("FoldList[Times, 1, Range[10]]",
             "List[1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800]");
}

static void test_foldlist_digits(void) {
    /* FoldList[10 #1 + #2 &, 0, {4, 5, 1, 6, 7, 8}]
       -> {0, 4, 45, 451, 4516, 45167, 451678} */
    run_full("FoldList[10 #1 + #2 &, 0, {4, 5, 1, 6, 7, 8}]",
             "List[0, 4, 45, 451, 4516, 45167, 451678]");
}

static void test_foldlist_horner_numeric(void) {
    /* FoldList[2 #1 + #2 &, 0, {1, 0, 1, 1}]
       -> {0, 1, 2, 5, 11} */
    run_full("FoldList[2 #1 + #2 &, 0, {1, 0, 1, 1}]",
             "List[0, 1, 2, 5, 11]");
}

static void test_foldlist_cumulative_sum_numeric(void) {
    /* Cumulative sums 1..5 -> {0, 1, 3, 6, 10, 15} */
    run_full("FoldList[Plus, 0, Range[5]]",
             "List[0, 1, 3, 6, 10, 15]");
}

static void test_foldlist_alternating_sums_numeric(void) {
    /* Numeric variant of the alternating-sum example: Mathilda does not
       distribute unary minus through Plus, so symbolic simplification
       would not collapse to the Mathematica "a - b + c - ..." form. With
       explicit numbers the arithmetic fully evaluates: at each step we
       subtract the previous partial result from the next input. */
    run_full("FoldList[#2 - #1 &, 0, {1, 2, 3, 4, 5}]",
             "List[0, 1, 1, 2, 2, 3]");
}

static void test_foldlist_continued_fraction_infix(void) {
    /* FoldList[1/(#2 + #1) &, x, Reverse[{a, b, c}]]
       -> {x, 1/(c + x), 1/(b + 1/(c + x)), 1/(a + 1/(b + 1/(c + x)))} */
    run_infix("FoldList[1/(#2 + #1) &, x, Reverse[{a, b, c}]]",
              "{x, 1/(c + x), 1/(b + 1/(c + x)), 1/(a + 1/(b + 1/(c + x)))}");
}

static void test_foldlist_empty_list(void) {
    /* FoldList[f, {}] -> {} (empty list, no seed) */
    run_full("FoldList[f, {}]", "List[]");
}

static void test_foldlist_empty_list_with_seed(void) {
    /* FoldList[f, x, {}] -> {x} */
    run_full("FoldList[f, x, {}]", "List[x]");
}

static void test_foldlist_length_property(void) {
    /* FoldList[f, x, list] has length n+1 where n = Length[list] */
    run_full("Length[FoldList[f, x, Range[7]]]", "8");
    run_full("Length[FoldList[f, Range[7]]]", "7");
}

static void test_foldlist_first_element(void) {
    /* First[FoldList[f, x, list]] is always x (with seed) */
    run_full("First[FoldList[f, x, {a, b, c, d}]]", "x");
    /* First[FoldList[f, list]] is always list[[1]] (no seed) */
    run_full("First[FoldList[f, {a, b, c, d}]]", "a");
}

static void test_foldlist_last_equals_fold(void) {
    /* Last[FoldList[f, x, list]] == Fold[f, x, list] */
    run_full("Last[FoldList[f, x, {a, b, c}]]",
             "f[f[f[x, a], b], c]");
    run_full("Fold[f, x, {a, b, c}]",
             "f[f[f[x, a], b], c]");
}

static void test_foldlist_empty_head_preserved(void) {
    /* FoldList[f, p[]] -> p[] (empty, preserving head) */
    run_full("FoldList[f, p[]]", "p[]");
}

static void test_foldlist_nonlist_head_empty_with_seed(void) {
    /* FoldList[f, x, p[]] -> p[x] */
    run_full("FoldList[f, x, p[]]", "p[x]");
}

/* ---------- Unevaluated / error cases ---------- */

static void test_foldlist_wrong_arg_count(void) {
    run_full("FoldList[f]", "FoldList[f]");
    run_full("FoldList[]", "FoldList[]");
    run_full("FoldList[f, x, {a}, extra]", "FoldList[f, x, List[a], extra]");
}

static void test_foldlist_list_is_symbol(void) {
    run_full("FoldList[f, x, y]", "FoldList[f, x, y]");
}

static void test_foldlist_attributes_protected(void) {
    run_full("MemberQ[Attributes[FoldList], Protected]", "True");
}

int main(void) {
    symtab_init();
    core_init();

    /* Spec examples */
    TEST(test_foldlist_basic_symbolic);
    TEST(test_foldlist_no_seed);
    TEST(test_foldlist_head_preserved);
    TEST(test_foldlist_fold_right);
    TEST(test_foldlist_factorials);
    TEST(test_foldlist_digits);
    TEST(test_foldlist_horner_numeric);
    TEST(test_foldlist_cumulative_sum_numeric);
    TEST(test_foldlist_alternating_sums_numeric);
    TEST(test_foldlist_continued_fraction_infix);
    TEST(test_foldlist_empty_list);
    TEST(test_foldlist_empty_list_with_seed);
    TEST(test_foldlist_length_property);
    TEST(test_foldlist_first_element);
    TEST(test_foldlist_last_equals_fold);
    TEST(test_foldlist_empty_head_preserved);
    TEST(test_foldlist_nonlist_head_empty_with_seed);

    /* Unevaluated cases */
    TEST(test_foldlist_wrong_arg_count);
    TEST(test_foldlist_list_is_symbol);
    TEST(test_foldlist_attributes_protected);

    printf("All FoldList tests passed!\n");
    return 0;
}
