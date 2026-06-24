/* Unit tests for IntegerPartitions (src/partitions.c).
 *
 * Note: assert_eval_eq uses libc assert(), which CMake's Release build
 * compiles out under -DNDEBUG. It always prints "FAIL: ..." to stderr on a
 * mismatch, so grep the output for FAIL even when the exit code is 0. */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "parse.h"
#include "test_utils.h"
#include <stdio.h>

/* ----- documented examples (exact WL output) --------------------------- */

static void test_all_of_5(void) {
    assert_eval_eq("IntegerPartitions[5]",
        "{{5}, {4, 1}, {3, 2}, {3, 1, 1}, {2, 2, 1}, {2, 1, 1, 1}, {1, 1, 1, 1, 1}}", 0);
}

static void test_at_most_k(void) {
    assert_eval_eq("IntegerPartitions[8, 3]",
        "{{8}, {7, 1}, {6, 2}, {6, 1, 1}, {5, 3}, {5, 2, 1}, "
        "{4, 4}, {4, 3, 1}, {4, 2, 2}, {3, 3, 2}}", 0);
}

static void test_exactly_k(void) {
    assert_eval_eq("IntegerPartitions[8, {3}]",
        "{{6, 1, 1}, {5, 2, 1}, {4, 3, 1}, {4, 2, 2}, {3, 3, 2}}", 0);
}

static void test_stepped_length(void) {
    /* even length only */
    assert_eval_eq("IntegerPartitions[6, {2, Infinity, 2}]",
        "{{5, 1}, {4, 2}, {3, 3}, {3, 1, 1, 1}, {2, 2, 1, 1}, {1, 1, 1, 1, 1, 1}}", 0);
}

static void test_restricted_parts(void) {
    assert_eval_eq("IntegerPartitions[8, All, {1, 2, 5}]",
        "{{5, 2, 1}, {5, 1, 1, 1}, {2, 2, 2, 2}, {2, 2, 2, 1, 1}, "
        "{2, 2, 1, 1, 1, 1}, {2, 1, 1, 1, 1, 1, 1}, {1, 1, 1, 1, 1, 1, 1, 1}}", 0);
}

static void test_rational_parts(void) {
    assert_eval_eq("IntegerPartitions[3, 10, {1, 1/3, 3/4}]",
        "{{3/4, 3/4, 3/4, 3/4}, "
        "{1/3, 1/3, 1/3, 1/3, 1/3, 1/3, 1/3, 1/3, 1/3}, "
        "{1/3, 1/3, 1/3, 1/3, 1/3, 1/3, 1}, "
        "{1/3, 1/3, 1/3, 1, 1}, {1, 1, 1}}", 0);
}

static void test_negative_parts(void) {
    assert_eval_eq("IntegerPartitions[5, 10, {1, -1}]",
        "{{-1, -1, 1, 1, 1, 1, 1, 1, 1}, {-1, 1, 1, 1, 1, 1, 1}, {1, 1, 1, 1, 1}}", 0);
}

static void test_first_m(void) {
    assert_eval_eq("IntegerPartitions[15, All, All, 10]",
        "{{15}, {14, 1}, {13, 2}, {13, 1, 1}, {12, 3}, {12, 2, 1}, "
        "{12, 1, 1, 1}, {11, 4}, {11, 3, 1}, {11, 2, 2}}", 0);
}

static void test_last_m(void) {
    assert_eval_eq("IntegerPartitions[15, All, All, -3]",
        "{{2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, "
        "{2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}, "
        "{1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}}", 0);
}

static void test_coin_change(void) {
    assert_eval_eq("IntegerPartitions[156, 10, {1, 5, 10, 25}]",
        "{{25, 25, 25, 25, 25, 25, 5, 1}, "
        "{25, 25, 25, 25, 25, 10, 10, 10, 1}, "
        "{25, 25, 25, 25, 25, 10, 10, 5, 5, 1}}", 0);
}

static void test_mcnugget(void) {
    assert_eval_eq("IntegerPartitions[50, All, {6, 9, 20}]",
        "{{20, 9, 9, 6, 6}, {20, 6, 6, 6, 6, 6}}", 0);
}

static void test_mcnugget_counts(void) {
    /* number of McNugget partitions for i = 1..20 */
    assert_eval_eq("Table[Length[IntegerPartitions[i, All, {6, 9, 20}]], {i, 20}]",
        "{0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 0, 2, 0, 1}", 0);
}

/* ----- edge cases ------------------------------------------------------- */

static void test_rational_n_empty(void) {
    assert_eval_eq("IntegerPartitions[1/2]", "{}", 0);
}

static void test_rational_n_parts(void) {
    assert_eval_eq("IntegerPartitions[1/2, All, {1/6, 1/3}]",
        "{{1/3, 1/6}, {1/6, 1/6, 1/6}}", 0);
}

static void test_zero(void) {
    assert_eval_eq("IntegerPartitions[0]", "{{}}", 0);
}

static void test_one(void) {
    assert_eval_eq("IntegerPartitions[1]", "{{1}}", 0);
}

static void test_kmin_kmax_range(void) {
    /* between 2 and 3 parts */
    assert_eval_eq("IntegerPartitions[5, {2, 3}]",
        "{{4, 1}, {3, 2}, {3, 1, 1}, {2, 2, 1}}", 0);
}

static void test_length_is_partitionsp(void) {
    assert_eval_eq("Length[IntegerPartitions[10]]", "42", 0);
    assert_eval_eq("Length[IntegerPartitions[20]]", "627", 0);
}

/* ----- diagnostics: result is left unevaluated ------------------------- */

static void test_undef_unevaluated(void) {
    /* infinite part-set with no length / count bound stays unevaluated */
    assert_eval_eq("IntegerPartitions[5, All, {1, -1}]",
        "IntegerPartitions[5, All, {1, -1}]", 0);
}

static void test_argb_unevaluated(void) {
    assert_eval_eq("IntegerPartitions[]", "IntegerPartitions[]", 0);
}

static void test_symbolic_unevaluated(void) {
    assert_eval_eq("IntegerPartitions[x]", "IntegerPartitions[x]", 0);
}

static void test_take_too_many_returns_all(void) {
    /* requesting more than available emits ::take and returns the lot */
    assert_eval_eq("IntegerPartitions[3, All, All, 7]",
        "{{3}, {2, 1}, {1, 1, 1}}", 0);
}

static void test_m_zero_empty(void) {
    assert_eval_eq("IntegerPartitions[5, All, All, 0]", "{}", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_all_of_5);
    TEST(test_at_most_k);
    TEST(test_exactly_k);
    TEST(test_stepped_length);
    TEST(test_restricted_parts);
    TEST(test_rational_parts);
    TEST(test_negative_parts);
    TEST(test_first_m);
    TEST(test_last_m);
    TEST(test_coin_change);
    TEST(test_mcnugget);
    TEST(test_mcnugget_counts);
    TEST(test_rational_n_empty);
    TEST(test_rational_n_parts);
    TEST(test_zero);
    TEST(test_one);
    TEST(test_kmin_kmax_range);
    TEST(test_length_is_partitionsp);
    TEST(test_undef_unevaluated);
    TEST(test_argb_unevaluated);
    TEST(test_symbolic_unevaluated);
    TEST(test_take_too_many_returns_all);
    TEST(test_m_zero_empty);

    printf("All IntegerPartitions tests passed!\n");
    return 0;
}
