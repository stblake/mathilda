/*
 * test_stringrepeat.c - unit tests for the StringRepeat builtin.
 *
 * StringRepeat has no regex dependency, so these tests are unconditional.
 */

#include "test_utils.h"

extern void symtab_init(void);
extern void core_init(void);

/* StringRepeat["s", n] - "s" concatenated n times. */
static void test_basic(void) {
    assert_eval_eq("StringRepeat[\"a\", 50]",
                   "\"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\"", 0);
    assert_eval_eq("StringRepeat[\"abc\", 10]",
                   "\"abcabcabcabcabcabcabcabcabcabc\"", 0);
    /* n == 1 reproduces the input. */
    assert_eval_eq("StringRepeat[\"ab\", 1]", "\"ab\"", 0);
    /* Multi-character block, small count. */
    assert_eval_eq("StringRepeat[\"xy\", 3]", "\"xyxyxy\"", 0);
}

/* StringRepeat["s", n, max] - up to n copies, truncated to length <= max. */
static void test_truncation(void) {
    /* 10 copies of "ab" is 20 chars; capped to 19 -> partial final copy. */
    assert_eval_eq("StringRepeat[\"ab\", 10, 19]",
                   "\"abababababababababa\"", 0);
    /* max exactly equals the full length -> untruncated. */
    assert_eval_eq("StringRepeat[\"ab\", 10, 20]",
                   "\"abababababababababab\"", 0);
    /* max larger than the full length -> untruncated (only n copies). */
    assert_eval_eq("StringRepeat[\"ab\", 3, 100]", "\"ababab\"", 0);
    /* max an exact multiple below n copies. */
    assert_eval_eq("StringRepeat[\"abc\", 10, 6]", "\"abcabc\"", 0);
    /* Mid-copy cutoff lands inside a block. */
    assert_eval_eq("StringRepeat[\"abc\", 10, 7]", "\"abcabca\"", 0);
    /* Single-character block: max simply bounds the count. */
    assert_eval_eq("StringRepeat[\"a\", 50, 5]", "\"aaaaa\"", 0);
}

/* Boundary conditions. */
static void test_edges(void) {
    /* Zero copies -> empty string. */
    assert_eval_eq("StringRepeat[\"abc\", 0]", "\"\"", 0);
    /* max == 0 -> empty string. */
    assert_eval_eq("StringRepeat[\"abc\", 10, 0]", "\"\"", 0);
    /* Empty base string -> empty string regardless of n or max. */
    assert_eval_eq("StringRepeat[\"\", 5]", "\"\"", 0);
    assert_eval_eq("StringRepeat[\"\", 5, 10]", "\"\"", 0);
    /* Single character. */
    assert_eval_eq("StringRepeat[\"a\", 1]", "\"a\"", 0);
    /* n == 0 with a max still gives the empty string. */
    assert_eval_eq("StringRepeat[\"ab\", 0, 10]", "\"\"", 0);
}

/* Invalid arguments leave the call unevaluated. */
static void test_unevaluated(void) {
    /* Non-string first argument. */
    assert_eval_eq("StringRepeat[xyz, 3]", "StringRepeat[xyz, 3]", 0);
    /* Non-integer count. */
    assert_eval_eq("StringRepeat[\"a\", x]", "StringRepeat[\"a\", x]", 0);
    /* Non-integer max. */
    assert_eval_eq("StringRepeat[\"a\", 3, y]", "StringRepeat[\"a\", 3, y]", 0);
    /* Negative count. */
    assert_eval_eq("StringRepeat[\"a\", -1]", "StringRepeat[\"a\", -1]", 0);
    /* Negative max. */
    assert_eval_eq("StringRepeat[\"a\", 3, -2]", "StringRepeat[\"a\", 3, -2]", 0);
    /* Bad arity -> ::argt message, call unchanged. */
    assert_eval_eq("StringRepeat[]", "StringRepeat[]", 0);
    assert_eval_eq("StringRepeat[\"a\"]", "StringRepeat[\"a\"]", 0);
    assert_eval_eq("StringRepeat[\"a\", 2, 5, 9]",
                   "StringRepeat[\"a\", 2, 5, 9]", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_basic);
    TEST(test_truncation);
    TEST(test_edges);
    TEST(test_unevaluated);

    printf("All StringRepeat tests passed!\n");
    return 0;
}
