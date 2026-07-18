/*
 * test_complement.c -- unit tests for the Complement builtin.
 *
 * Complement[eall, e1, e2, ...] gives the sorted list of distinct elements of
 * eall that appear in none of the ei, using the head of the first argument.
 * With SameTest -> f a custom equivalence relation is used; Wolfram keeps the
 * canonically-greatest member of each equivalence class as the representative
 * (the SameTest cases below pin that behaviour).
 */

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

/* Evaluate input and return its printed form; caller must free. */
static char* eval_to_string(const char* src) {
    Expr* parsed = parse_expression(src);
    assert(parsed != NULL);
    Expr* evaluated = evaluate(parsed);
    expr_free(parsed);
    char* s = expr_to_string(evaluated);
    expr_free(evaluated);
    return s;
}

/* --- Basic set difference ------------------------------------------ */

static void test_basic_multilist(void) {
    assert_eval_eq("Complement[{a, b, c, d, e}, {a, c}, {d}]", "{b, e}", 0);
}

static void test_two_list(void) {
    assert_eval_eq("Complement[{1, 2, 3, 4, 5}, {2, 4}]", "{1, 3, 5}", 0);
}

/* Result is sorted and deduplicated regardless of input order/multiplicity. */
static void test_sorted_and_deduped(void) {
    assert_eval_eq("Complement[{b, e, d, a, b, c, d}, {b, c}]", "{a, d, e}", 0);
}

/* Single list: sorted distinct elements (like Union[list]). */
static void test_single_list(void) {
    assert_eval_eq("Complement[{c, a, b, a}]", "{a, b, c}", 0);
}

/* Everything removed yields {}. */
static void test_empty_result(void) {
    assert_eval_eq("Complement[{a, b}, {a, b, c}]", "{}", 0);
}

/* Subtracting an empty list is identity (sorted + deduped). */
static void test_with_empty_list(void) {
    assert_eval_eq("Complement[{3, 1, 2, 1}, {}]", "{1, 2, 3}", 0);
}

/* Subtracting disjoint lists leaves the first (sorted + deduped). */
static void test_disjoint(void) {
    assert_eval_eq("Complement[{a, b, c}, {d, e, f}]", "{a, b, c}", 0);
}

/* --- Arbitrary heads (not just List) ------------------------------- */

static void test_arbitrary_head(void) {
    assert_eval_eq("Complement[f[a, b, c, d], f[c, a], f[b, b, a]]", "f[d]", 0);
}

/* Mismatched heads leave the expression unevaluated. */
static void test_mismatched_heads_unevaluated(void) {
    char* s = eval_to_string("Complement[f[a, b], g[a]]");
    ASSERT_STR_EQ(s, "Complement[f[a, b], g[a]]");
    free(s);
}

/* --- Integration with other builtins ------------------------------- */

static void test_range_odds(void) {
    assert_eval_eq("Complement[Range[10], Range[2, 10, 2]]", "{1, 3, 5, 7, 9}", 0);
}

static void test_divisors(void) {
    assert_eval_eq("Complement[Divisors[12], Divisors[6]]", "{4, 12}", 0);
}

/* --- SameTest equivalence classes ---------------------------------- */

static void test_sametest_abs(void) {
    assert_eval_eq(
        "Complement[{2, -2, 1, 3}, {2, 1, -2, -1}, SameTest -> (Abs[#1] == Abs[#2] &)]",
        "{3}", 0);
}

static void test_sametest_floor(void) {
    assert_eval_eq(
        "Complement[{1.1, 3.4, .5, 7.6, 7.1, 1.9}, {1.2, 3.3, 1.3}, "
        "SameTest -> (Floor[#1] == Floor[#2] &)]",
        "{0.5, 7.1}", 0);
}

static void test_sametest_total(void) {
    assert_eval_eq(
        "Complement[{{1, 2}, {3}, {4, 5, 6}, {9, 5}}, {{2, 1}, {8, 4, 3}}, "
        "SameTest -> (Total[#1] == Total[#2] &)]",
        "{{9, 5}}", 0);
}

/* Explicit SameTest -> Automatic behaves as the default path. */
static void test_sametest_automatic(void) {
    assert_eval_eq("Complement[{3, 1, 2, 1}, {2}, SameTest -> Automatic]", "{1, 3}", 0);
}

/* --- Options ------------------------------------------------------- */

static void test_options(void) {
    assert_eval_eq("Options[Complement]", "{SameTest -> Automatic}", 0);
}

/* --- Attributes ---------------------------------------------------- */

static void test_attributes(void) {
    assert_eval_eq("MemberQ[Attributes[Complement], Protected]", "True", 0);
    /* Complement is order-sensitive, so unlike Intersection it is NOT Flat. */
    assert_eval_eq("MemberQ[Attributes[Complement], Flat]", "False", 0);
    assert_eval_eq("MemberQ[Attributes[Complement], OneIdentity]", "False", 0);
}

/* Zero-argument call reports an error and returns unevaluated (the error
 * message goes to stderr; we only assert the observable return form). */
static void test_no_args_unevaluated(void) {
    char* s = eval_to_string("Complement[]");
    ASSERT_STR_EQ(s, "Complement[]");
    free(s);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_basic_multilist);
    TEST(test_two_list);
    TEST(test_sorted_and_deduped);
    TEST(test_single_list);
    TEST(test_empty_result);
    TEST(test_with_empty_list);
    TEST(test_disjoint);
    TEST(test_arbitrary_head);
    TEST(test_mismatched_heads_unevaluated);
    TEST(test_range_odds);
    TEST(test_divisors);
    TEST(test_sametest_abs);
    TEST(test_sametest_floor);
    TEST(test_sametest_total);
    TEST(test_sametest_automatic);
    TEST(test_options);
    TEST(test_attributes);
    TEST(test_no_args_unevaluated);

    printf("All tests passed!\n");
    return 0;
}
