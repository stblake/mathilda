/*
 * test_intersection.c -- unit tests for the Intersection builtin.
 *
 * Intersection[l1, l2, ...] gives the sorted list of elements common to every
 * li, deduplicated, using the head of the first argument. With SameTest -> f
 * a custom equivalence relation is used; Wolfram keeps the canonically-greatest
 * member of each equivalence class as the representative (the SameTest cases
 * below pin that behaviour).
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

/* --- Basic set intersection ---------------------------------------- */

static void test_basic_multilist(void) {
    assert_eval_eq("Intersection[{1, 1, 2, 3}, {3, 1, 4}, {4, 1, 3, 3}]", "{1, 3}", 0);
}

static void test_two_list_overlap(void) {
    assert_eval_eq("Intersection[{1, 2, 3, 4}, {2, 4, 6}]", "{2, 4}", 0);
}

static void test_sorted_and_deduped(void) {
    /* Result is sorted and distinct regardless of input order/multiplicity. */
    assert_eval_eq("Intersection[{5, 3, 3, 1, 5}, {1, 3, 5, 5}]", "{1, 3, 5}", 0);
}

/* Single list: sorted distinct elements (like Union[list]). */
static void test_single_list(void) {
    assert_eval_eq("Intersection[{3, 1, 2, 1}]", "{1, 2, 3}", 0);
}

/* Empty intersection yields {}. */
static void test_empty_result(void) {
    assert_eval_eq("Intersection[{a, b, c}, {d, e, f}]", "{}", 0);
}

/* Intersecting with an empty list gives {}. */
static void test_with_empty_list(void) {
    assert_eval_eq("Intersection[{1, 2, 3}, {}]", "{}", 0);
}

/* Symbolic elements sort canonically. */
static void test_symbolic_elements(void) {
    assert_eval_eq("Intersection[{a, b, c, d}, {b, d, e}]", "{b, d}", 0);
}

/* --- Arbitrary heads (not just List) ------------------------------- */

static void test_arbitrary_head(void) {
    assert_eval_eq("Intersection[f[a, b], f[c, a], f[b, b, a]]", "f[a]", 0);
}

static void test_arbitrary_head_multi(void) {
    assert_eval_eq("Intersection[g[3, 1, 2], g[2, 3]]", "g[2, 3]", 0);
}

/* Mismatched heads leave the expression unevaluated. */
static void test_mismatched_heads_unevaluated(void) {
    char* s = eval_to_string("Intersection[f[a, b], g[a, b]]");
    ASSERT_STR_EQ(s, "Intersection[f[a, b], g[a, b]]");
    free(s);
}

/* --- Integration with other builtins ------------------------------- */

static void test_divisors(void) {
    assert_eval_eq("Intersection[Divisors[45], Divisors[78]]", "{1, 3}", 0);
}

/* --- SameTest equivalence classes ---------------------------------- */

static void test_sametest_abs(void) {
    assert_eval_eq(
        "Intersection[{2, -2, 1, 3, 1}, {2, 1, -2, -1}, SameTest -> (Abs[#1] == Abs[#2] &)]",
        "{1, 2}", 0);
}

static void test_sametest_floor(void) {
    assert_eval_eq(
        "Intersection[{1.1, 3.4, .5, 7.6, 7.1, 1.9}, {1.2, 3.3, 7.7, 1.3}, "
        "SameTest -> (Floor[#1] == Floor[#2] &)]",
        "{1.9, 3.4, 7.6}", 0);
}

static void test_sametest_total(void) {
    assert_eval_eq(
        "Intersection[{{1, 2}, {3}, {4, 5, 6}, {9, 6}}, {{2, 1}, {8, 4, 3}}, "
        "SameTest -> (Total[#1] == Total[#2] &)]",
        "{{1, 2}, {4, 5, 6}}", 0);
}

/* SameTest with a single list: greatest representative per class. */
static void test_sametest_single_list(void) {
    assert_eval_eq(
        "Intersection[{2, -2, 1, 3}, SameTest -> (Abs[#1] == Abs[#2] &)]",
        "{1, 2, 3}", 0);
}

/* --- Attributes ---------------------------------------------------- */

static void test_attributes(void) {
    assert_eval_eq("MemberQ[Attributes[Intersection], Flat]", "True", 0);
    assert_eval_eq("MemberQ[Attributes[Intersection], OneIdentity]", "True", 0);
    assert_eval_eq("MemberQ[Attributes[Intersection], Protected]", "True", 0);
}

/* Zero-argument call reports an error and returns unevaluated (the error
 * message goes to stderr; we only assert the observable return form). */
static void test_no_args_unevaluated(void) {
    char* s = eval_to_string("Intersection[]");
    ASSERT_STR_EQ(s, "Intersection[]");
    free(s);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_basic_multilist);
    TEST(test_two_list_overlap);
    TEST(test_sorted_and_deduped);
    TEST(test_single_list);
    TEST(test_empty_result);
    TEST(test_with_empty_list);
    TEST(test_symbolic_elements);
    TEST(test_arbitrary_head);
    TEST(test_arbitrary_head_multi);
    TEST(test_mismatched_heads_unevaluated);
    TEST(test_divisors);
    TEST(test_sametest_abs);
    TEST(test_sametest_floor);
    TEST(test_sametest_total);
    TEST(test_sametest_single_list);
    TEST(test_attributes);
    TEST(test_no_args_unevaluated);

    printf("All tests passed!\n");
    return 0;
}
