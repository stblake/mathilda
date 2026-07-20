/*
 * test_stringpad.c - unit tests for StringPadLeft and StringPadRight.
 *
 * These builtins have no regex dependency, so the tests are unconditional.
 */

#include "test_utils.h"

extern void symtab_init(void);
extern void core_init(void);

/* Basic padding with spaces on a single string. */
static void test_basic(void) {
    assert_eval_eq("StringPadLeft[\"abcde\", 10]",  "\"     abcde\"", 0);
    assert_eval_eq("StringPadRight[\"abcde\", 10]", "\"abcde     \"", 0);
    /* n == length -> unchanged. */
    assert_eval_eq("StringPadLeft[\"abcde\", 5]",  "\"abcde\"", 0);
    assert_eval_eq("StringPadRight[\"abcde\", 5]", "\"abcde\"", 0);
    /* Single-argument single string -> unchanged. */
    assert_eval_eq("StringPadLeft[\"abcde\"]",  "\"abcde\"", 0);
    assert_eval_eq("StringPadRight[\"abcde\"]", "\"abcde\"", 0);
}

/* Padding with an explicit (single-character) pad string. */
static void test_pad_char(void) {
    assert_eval_eq("StringPadLeft[\"abcde\", 10, \".\"]",  "\".....abcde\"", 0);
    assert_eval_eq("StringPadRight[\"abcde\", 10, \".\"]", "\"abcde.....\"", 0);
    assert_eval_eq("StringPadLeft[\"100\", 6, \"0\"]",  "\"000100\"", 0);
}

/* Multi-character pad string: cyclic copies read left-to-right (p[i % plen]). */
static void test_pad_multichar(void) {
    /* pad width 5, "ab" -> "ababa". */
    assert_eval_eq("StringPadLeft[\"x\", 6, \"ab\"]",  "\"ababax\"", 0);
    assert_eval_eq("StringPadRight[\"x\", 6, \"ab\"]", "\"xababa\"", 0);
    /* pad width 7, "xyz" -> "xyzxyzx". */
    assert_eval_eq("StringPadLeft[\"cat\", 10, \"xyz\"]",  "\"xyzxyzxcat\"", 0);
    assert_eval_eq("StringPadRight[\"cat\", 10, \"xyz\"]", "\"catxyzxyzx\"", 0);
}

/* Truncation when n < length: left keeps the last n, right keeps the first n. */
static void test_truncation(void) {
    assert_eval_eq("StringPadLeft[\"abcde\", 3]",  "\"cde\"", 0);
    assert_eval_eq("StringPadRight[\"abcde\", 3]", "\"abc\"", 0);
    /* Truncation ignores the pad string. */
    assert_eval_eq("StringPadLeft[\"abcde\", 2, \".\"]",  "\"de\"", 0);
    assert_eval_eq("StringPadRight[\"abcde\", 2, \".\"]", "\"ab\"", 0);
}

/* List, single argument: pad each to the length of the longest. */
static void test_list_maxlen(void) {
    assert_eval_eq(
        "StringPadLeft[{\"a\", \"ab\", \"abc\", \"abcd\", \"abcde\"}]",
        "{\"    a\", \"   ab\", \"  abc\", \" abcd\", \"abcde\"}", 0);
    assert_eval_eq(
        "StringPadRight[{\"a\", \"ab\", \"abc\", \"abcd\", \"abcde\"}]",
        "{\"a    \", \"ab   \", \"abc  \", \"abcd \", \"abcde\"}", 0);
}

/* List with an explicit length: pad or truncate each element. */
static void test_list_n(void) {
    assert_eval_eq(
        "StringPadLeft[{\"a\", \"ab\", \"abc\", \"abcd\", \"abcde\"}, 3]",
        "{\"  a\", \" ab\", \"abc\", \"bcd\", \"cde\"}", 0);
    assert_eval_eq(
        "StringPadRight[{\"a\", \"ab\", \"abc\", \"abcd\", \"abcde\"}, 3]",
        "{\"a  \", \"ab \", \"abc\", \"abc\", \"abc\"}", 0);
}

/* List with a length and a pad string. */
static void test_list_n_pad(void) {
    assert_eval_eq(
        "StringPadLeft[{\"a\", \"bb\", \"ccc\"}, 4, \"-\"]",
        "{\"---a\", \"--bb\", \"-ccc\"}", 0);
    assert_eval_eq(
        "StringPadRight[{\"a\", \"bb\", \"ccc\"}, 4, \"-\"]",
        "{\"a---\", \"bb--\", \"ccc-\"}", 0);
    /* Empty list stays empty. */
    assert_eval_eq("StringPadLeft[{}, 3]", "{}", 0);
    assert_eval_eq("StringPadRight[{}]", "{}", 0);
}

/* Boundary conditions. */
static void test_edges(void) {
    /* n == 0 -> empty string. */
    assert_eval_eq("StringPadLeft[\"abc\", 0]",  "\"\"", 0);
    assert_eval_eq("StringPadRight[\"abc\", 0]", "\"\"", 0);
    /* Empty input string padded up. */
    assert_eval_eq("StringPadLeft[\"\", 3]",  "\"   \"", 0);
    assert_eval_eq("StringPadRight[\"\", 3, \"x\"]", "\"xxx\"", 0);
    /* Empty input, empty pad, no width needed -> empty (w == 0 is fine). */
    assert_eval_eq("StringPadLeft[\"\", 0, \"\"]", "\"\"", 0);
    /* Empty input, 1-arg -> empty. */
    assert_eval_eq("StringPadRight[\"\"]", "\"\"", 0);
}

/* Invalid arguments leave the call unevaluated. */
static void test_unevaluated(void) {
    /* Non-string, non-list first argument. */
    assert_eval_eq("StringPadLeft[xyz, 3]", "StringPadLeft[xyz, 3]", 0);
    /* Non-integer length. */
    assert_eval_eq("StringPadLeft[\"a\", x]", "StringPadLeft[\"a\", x]", 0);
    /* Negative length. */
    assert_eval_eq("StringPadLeft[\"a\", -1]", "StringPadLeft[\"a\", -1]", 0);
    /* Non-string pad. */
    assert_eval_eq("StringPadLeft[\"a\", 5, 3]", "StringPadLeft[\"a\", 5, 3]", 0);
    /* List pad (Wolfram's list form is not supported). */
    assert_eval_eq("StringPadLeft[\"a\", 5, {\"x\", \"y\"}]",
                   "StringPadLeft[\"a\", 5, {\"x\", \"y\"}]", 0);
    /* Empty pad while padding is required. */
    assert_eval_eq("StringPadLeft[\"a\", 5, \"\"]",
                   "StringPadLeft[\"a\", 5, \"\"]", 0);
    /* Non-string element in the list. */
    assert_eval_eq("StringPadLeft[{\"a\", 7}, 3]", "StringPadLeft[{\"a\", 7}, 3]", 0);
    /* Bad arity -> ::argb message, call unchanged. */
    assert_eval_eq("StringPadLeft[]", "StringPadLeft[]", 0);
    assert_eval_eq("StringPadRight[]", "StringPadRight[]", 0);
    assert_eval_eq("StringPadLeft[\"a\", 2, \"x\", 9]",
                   "StringPadLeft[\"a\", 2, \"x\", 9]", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_basic);
    TEST(test_pad_char);
    TEST(test_pad_multichar);
    TEST(test_truncation);
    TEST(test_list_maxlen);
    TEST(test_list_n);
    TEST(test_list_n_pad);
    TEST(test_edges);
    TEST(test_unevaluated);

    printf("All StringPad tests passed!\n");
    return 0;
}
