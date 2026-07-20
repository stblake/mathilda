/*
 * test_stringriffle.c - unit tests for StringRiffle.
 *
 * StringRiffle has no regex dependency, so the tests are unconditional.
 *
 * Note on newlines: the Mathilda lexer does not decode backslash escapes (an
 * input "\n" parses to the literal string "n"), so the default newline scheme
 * is exercised via real newline bytes in the *expected* C literals (our C code
 * emits real newlines), while explicit-separator tests use printable
 * separators.
 */

#include "test_utils.h"

extern void symtab_init(void);
extern void core_init(void);

/* Flat list: default separator is a single space. */
static void test_default_flat(void) {
    assert_eval_eq("StringRiffle[{\"a\", \"b\", \"c\", \"d\", \"e\"}]",
                   "\"a b c d e\"", 0);
    /* Single element -> no separator. */
    assert_eval_eq("StringRiffle[{\"a\"}]", "\"a\"", 0);
    /* Non-string leaves are converted with ToString. */
    assert_eval_eq("StringRiffle[{1, 2, 3}]", "\"1 2 3\"", 0);
}

/* An explicit string separator replaces the top-level default. */
static void test_explicit_sep(void) {
    assert_eval_eq("StringRiffle[{\"a\", \"b\", \"c\", \"d\", \"e\"}, \", \"]",
                   "\"a, b, c, d, e\"", 0);
    assert_eval_eq("StringRiffle[{1, 2, 3}, \"-\"]", "\"1-2-3\"", 0);
}

/* A 3-string list separator is a {left, sep, right} delimiter triple. */
static void test_delimiter_triple(void) {
    assert_eval_eq(
        "StringRiffle[{\"a\", \"b\", \"c\", \"d\", \"e\"}, {\"(\", \" \", \")\"}]",
        "\"(a b c d e)\"", 0);
}

/* Default nested scheme: spaces at the bottom, newlines going up. The expected
 * literals contain real newline bytes. */
static void test_default_nested(void) {
    /* 2-D: rows separated by one newline, cells by a space. */
    assert_eval_eq(
        "StringRiffle[{{\"a\", \"b\", \"c\"}, {\"d\", \"e\", \"f\"}}]",
        "\"a b c\nd e f\"", 0);
    /* 3-D: blocks by two newlines, rows by one, cells by a space. */
    assert_eval_eq(
        "StringRiffle[{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}}]",
        "\"1 2\n3 4\n\n5 6\n7 8\"", 0);
}

/* Per-level separators (level 1 = outermost). Uses printable separators. */
static void test_multilevel_sep(void) {
    assert_eval_eq(
        "StringRiffle[{{\"a\", \"b\", \"c\"}, {\"d\", \"e\", \"f\"}}, \"|\", \"-\"]",
        "\"a-b-c|d-e-f\"", 0);
    /* JSON-like: a delimiter triple at level 1, a plain sep at level 2, with
     * integer leaves converted via ToString. */
    assert_eval_eq(
        "StringRiffle[{{\"a\", 27}, {\"b\", 28}, {\"c\", 29}}, "
        "{\"{\", \", \", \"}\"}, \": \"]",
        "\"{a: 27, b: 28, c: 29}\"", 0);
}

/* Deeper levels than supplied separators fall back to the default scheme. */
static void test_partial_sep_defaults(void) {
    /* Level-1 sep "|" between rows; the cells use the default space. */
    assert_eval_eq(
        "StringRiffle[{{\"a\", \"b\"}, {\"c\", \"d\"}}, \"|\"]",
        "\"a b|c d\"", 0);
}

/* Boundary conditions. */
static void test_edges(void) {
    /* Empty list -> "". */
    assert_eval_eq("StringRiffle[{}]", "\"\"", 0);
    /* Empty list with a delimiter triple -> just the delimiters. */
    assert_eval_eq("StringRiffle[{}, {\"[\", \",\", \"]\"}]", "\"[]\"", 0);
    /* A bare string argument is returned unchanged. */
    assert_eval_eq("StringRiffle[\"abc\"]", "\"abc\"", 0);
    /* Nested empty rows: two empty rows joined by a newline. */
    assert_eval_eq("StringRiffle[{{}, {}}]", "\"\n\"", 0);
}

/* Invalid arguments leave the call unevaluated. */
static void test_unevaluated(void) {
    /* No arguments -> StringRiffle::argm, call unchanged. */
    assert_eval_eq("StringRiffle[]", "StringRiffle[]", 0);
    /* First argument neither a list nor a string. */
    assert_eval_eq("StringRiffle[x]", "StringRiffle[x]", 0);
    /* Malformed separator: not a string, not a 3-string list. */
    assert_eval_eq("StringRiffle[{\"a\", \"b\"}, 5]",
                   "StringRiffle[{\"a\", \"b\"}, 5]", 0);
    /* A 2-element list is not a valid delimiter triple. */
    assert_eval_eq("StringRiffle[{\"a\", \"b\"}, {\"(\", \")\"}]",
                   "StringRiffle[{\"a\", \"b\"}, {\"(\", \")\"}]", 0);
    /* A 3-list with a non-string element is not a valid triple. */
    assert_eval_eq("StringRiffle[{\"a\", \"b\"}, {\"(\", 5, \")\"}]",
                   "StringRiffle[{\"a\", \"b\"}, {\"(\", 5, \")\"}]", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_default_flat);
    TEST(test_explicit_sep);
    TEST(test_delimiter_triple);
    TEST(test_default_nested);
    TEST(test_multilevel_sep);
    TEST(test_partial_sep_defaults);
    TEST(test_edges);
    TEST(test_unevaluated);

    printf("All StringRiffle tests passed!\n");
    return 0;
}
