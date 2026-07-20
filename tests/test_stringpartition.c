/*
 * test_stringpartition.c - unit tests for the StringPartition builtin.
 *
 * StringPartition has no regex dependency, so these tests are unconditional
 * (unlike test_stringfns.c, which is USE_REGEX-gated).
 */

#include "test_utils.h"

extern void symtab_init(void);
extern void core_init(void);

/* StringPartition["s", n] - non-overlapping length-n blocks (offset d = n). */
static void test_basic(void) {
    assert_eval_eq("StringPartition[\"123456789123456789\", 9]",
                   "{\"123456789\", \"123456789\"}", 0);
    /* Trailing chars that cannot fill a full block are dropped. */
    assert_eval_eq("StringPartition[\"ababababab\", 3]",
                   "{\"aba\", \"bab\", \"aba\"}", 0);
    /* Exact multiple. */
    assert_eval_eq("StringPartition[\"abcdef\", 2]",
                   "{\"ab\", \"cd\", \"ef\"}", 0);
    /* n == length. */
    assert_eval_eq("StringPartition[\"abc\", 3]", "{\"abc\"}", 0);
    /* n larger than the string, plain form -> no block fits. */
    assert_eval_eq("StringPartition[\"abc\", 5]", "{}", 0);
}

/* StringPartition["s", n, d] - length-n blocks starting every d chars. */
static void test_offset(void) {
    /* d < n: overlapping windows. */
    assert_eval_eq("StringPartition[\"123456789\", 2, 1]",
                   "{\"12\", \"23\", \"34\", \"45\", \"56\", \"67\", \"78\", \"89\"}", 0);
    assert_eval_eq("StringPartition[\"12345\", 3, 1]",
                   "{\"123\", \"234\", \"345\"}", 0);
    /* d == n reproduces the two-argument form. */
    assert_eval_eq("StringPartition[\"abcdef\", 2, 2]",
                   "{\"ab\", \"cd\", \"ef\"}", 0);
    /* d > n: characters in the middle are skipped. */
    assert_eval_eq("StringPartition[\"123456789\", 2, 3]",
                   "{\"12\", \"45\", \"78\"}", 0);
}

/* StringPartition["s", UpTo[n]] - allow a shorter final block. */
static void test_upto(void) {
    assert_eval_eq("StringPartition[\"123456789\", UpTo[2]]",
                   "{\"12\", \"34\", \"56\", \"78\", \"9\"}", 0);
    /* Exact fit: no short tail. */
    assert_eval_eq("StringPartition[\"abcdef\", UpTo[2]]",
                   "{\"ab\", \"cd\", \"ef\"}", 0);
    /* n larger than the string -> the whole string as one short block. */
    assert_eval_eq("StringPartition[\"abc\", UpTo[5]]", "{\"abc\"}", 0);
    /* UpTo combined with an explicit offset. */
    assert_eval_eq("StringPartition[\"1234567\", UpTo[3], 3]",
                   "{\"123\", \"456\", \"7\"}", 0);
}

/* Threading over a list of strings. */
static void test_list_threading(void) {
    assert_eval_eq("StringPartition[{\"abcd\", \"ef\"}, 2]",
                   "{{\"ab\", \"cd\"}, {\"ef\"}}", 0);
    assert_eval_eq("StringPartition[{\"abcd\"}, 2, 1]",
                   "{{\"ab\", \"bc\", \"cd\"}}", 0);
    /* Empty list threads to an empty list. */
    assert_eval_eq("StringPartition[{}, 2]", "{}", 0);
}

/* Boundary conditions. */
static void test_edges(void) {
    /* Empty string yields no blocks. */
    assert_eval_eq("StringPartition[\"\", 2]", "{}", 0);
    assert_eval_eq("StringPartition[\"\", UpTo[2]]", "{}", 0);
    /* Single character. */
    assert_eval_eq("StringPartition[\"a\", 1]", "{\"a\"}", 0);
}

/* Invalid arguments leave the call unevaluated. */
static void test_unevaluated(void) {
    /* Non-string first argument. */
    assert_eval_eq("StringPartition[xyz, 2]", "StringPartition[xyz, 2]", 0);
    /* Non-positive length or offset (strings print with their quotes). */
    assert_eval_eq("StringPartition[\"abc\", 0]", "StringPartition[\"abc\", 0]", 0);
    assert_eval_eq("StringPartition[\"abc\", -2]", "StringPartition[\"abc\", -2]", 0);
    assert_eval_eq("StringPartition[\"abc\", 2, 0]", "StringPartition[\"abc\", 2, 0]", 0);
    /* Non-integer spec. */
    assert_eval_eq("StringPartition[\"abc\", x]", "StringPartition[\"abc\", x]", 0);
    /* Bad arity -> ::argt message, call unchanged. */
    assert_eval_eq("StringPartition[]", "StringPartition[]", 0);
    assert_eval_eq("StringPartition[\"abc\"]", "StringPartition[\"abc\"]", 0);
    assert_eval_eq("StringPartition[\"abc\", 2, 1, 5]",
                   "StringPartition[\"abc\", 2, 1, 5]", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_basic);
    TEST(test_offset);
    TEST(test_upto);
    TEST(test_list_threading);
    TEST(test_edges);
    TEST(test_unevaluated);

    printf("All StringPartition tests passed!\n");
    return 0;
}
