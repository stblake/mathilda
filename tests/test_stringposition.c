/*
 * test_stringposition.c - unit tests for the StringPosition builtin.
 *
 * StringPosition is backed by the PCRE2 regex engine, so the whole suite is a
 * no-op when built without PCRE2 (USE_REGEX undefined), mirroring
 * test_stringfns.c.  Cases are taken directly from the StringPosition
 * specification.  Mathilda's lexer drops one backslash per escape, so a regex
 * backslash is written "\\" in Mathilda source (four in this C source); none of
 * these cases need one.
 */

#include "test_utils.h"

extern void symtab_init(void);
extern void core_init(void);

#ifdef USE_REGEX

/* Basic literal substring positions, in StringTake form. */
static void test_literal(void) {
    assert_eval_eq("StringPosition[\"abXYZaaabXYZaaaaXYZXYZ\", \"XYZ\"]",
                   "{{3, 5}, {10, 12}, {17, 19}, {20, 22}}", 0);
    assert_eval_eq("StringPosition[\"XYZabc\", \"XYZ\"]", "{{1, 3}}", 0);
    /* No match -> empty list. */
    assert_eval_eq("StringPosition[\"abc\", \"z\"]", "{}", 0);
}

/* Third positional argument caps the number of occurrences. */
static void test_count(void) {
    assert_eval_eq("StringPosition[\"abXYZaaabXYZaaaaXYZXYZ\", \"XYZ\", 1]",
                   "{{3, 5}}", 0);
    assert_eval_eq("StringPosition[\"abXYZaaabXYZaaaaXYZXYZ\", \"XYZ\", 2]",
                   "{{3, 5}, {10, 12}}", 0);
    /* n larger than the number of matches returns them all. */
    assert_eval_eq("StringPosition[\"XYZabc\", \"XYZ\", 5]", "{{1, 3}}", 0);
}

/* x_ ~~ x_ : two equal characters (named-pattern backreference). */
static void test_repeated_pattern(void) {
    assert_eval_eq("StringPosition[\"AABBBAABABBCCCBAAA\", x_ ~~ x_]",
                   "{{1, 2}, {3, 4}, {4, 5}, {6, 7}, {10, 11}, {12, 13}, "
                   "{13, 14}, {16, 17}, {17, 18}}", 0);
    assert_eval_eq("StringPosition[\"AABBBAABABBCCCBAAA\", x_ ~~ x_, "
                   "Overlaps -> False]",
                   "{{1, 2}, {3, 4}, {6, 7}, {10, 11}, {12, 13}, {16, 17}}", 0);
}

/* Overlaps -> True (default) versus False for a plain literal pattern. */
static void test_overlaps_literal(void) {
    assert_eval_eq("StringPosition[\"AAAAA\", \"AA\"]",
                   "{{1, 2}, {2, 3}, {3, 4}, {4, 5}}", 0);
    assert_eval_eq("StringPosition[\"AAAAA\", \"AA\", Overlaps -> False]",
                   "{{1, 2}, {3, 4}}", 0);
    /* Overlaps -> True is the explicit default. */
    assert_eval_eq("StringPosition[\"AAAAA\", \"AA\", Overlaps -> True]",
                   "{{1, 2}, {2, 3}, {3, 4}, {4, 5}}", 0);
}

/* A list of patterns: matches of each are merged, sorted by start position. */
static void test_pattern_list(void) {
    assert_eval_eq("StringPosition[\"ABAABBAABABB\", {\"ABA\", \"AA\"}]",
                   "{{1, 3}, {3, 4}, {7, 8}, {8, 10}}", 0);
    assert_eval_eq("StringPosition[\"ABAABBAABABB\", {\"ABA\", \"AA\"}, "
                   "Overlaps -> False]",
                   "{{1, 3}, {7, 8}}", 0);
}

/* IgnoreCase folds case for both literal and named patterns. */
static void test_ignorecase(void) {
    assert_eval_eq("StringPosition[\"abAB\", \"a\", IgnoreCase -> True]",
                   "{{1, 1}, {3, 3}}", 0);
    assert_eval_eq("StringPosition[\"abAB\", \"a\", IgnoreCase -> False]",
                   "{{1, 1}}", 0);
    assert_eval_eq("StringPosition[\"aAbaBabBBaABAaBa\", x_ ~~ x_, "
                   "IgnoreCase -> True]",
                   "{{1, 2}, {7, 8}, {8, 9}, {10, 11}, {13, 14}}", 0);
    assert_eval_eq("StringPosition[\"aAbaBabBBaABAaBa\", x_ ~~ x_, "
                   "IgnoreCase -> False]",
                   "{{8, 9}}", 0);
}

/* BlankSequence x__ : default keeps one (longest) substring per start. */
static void test_blanksequence(void) {
    assert_eval_eq("StringPosition[\"AAAA\", x__]",
                   "{{1, 4}, {2, 4}, {3, 4}, {4, 4}}", 0);
    /* Overlaps -> All includes every matching substring at every start. */
    assert_eval_eq("StringPosition[\"AAAA\", x__, Overlaps -> All]",
                   "{{1, 4}, {1, 3}, {1, 2}, {1, 1}, {2, 4}, {2, 3}, {2, 2}, "
                   "{3, 4}, {3, 3}, {4, 4}}", 0);
}

/* A list of patterns differs from Alternatives / a single StringExpression. */
static void test_list_vs_alternatives(void) {
    assert_eval_eq("StringPosition[\"AAAA\", {x__, \"A\"}]",
                   "{{1, 4}, {1, 1}, {2, 4}, {2, 2}, {3, 4}, {3, 3}, "
                   "{4, 4}, {4, 4}}", 0);
    assert_eval_eq("StringPosition[\"AAAA\", StringExpression[{x__, \"A\"}]]",
                   "{{1, 4}, {2, 4}, {3, 4}, {4, 4}}", 0);
    assert_eval_eq("StringPosition[\"AAAA\", x__ | \"A\"]",
                   "{{1, 4}, {2, 4}, {3, 4}, {4, 4}}", 0);
}

/* RegularExpression and character-class patterns. */
static void test_regex_and_classes(void) {
    assert_eval_eq("StringPosition[\"a12b345\", DigitCharacter]",
                   "{{2, 2}, {3, 3}, {5, 5}, {6, 6}, {7, 7}}", 0);
    /* Greedy [0-9]+ under the default Overlaps -> True reports a (longest) match
     * starting at each digit position, exactly like the x__ example. */
    assert_eval_eq("StringPosition[\"a12b345\", RegularExpression[\"[0-9]+\"]]",
                   "{{2, 3}, {3, 3}, {5, 7}, {6, 7}, {7, 7}}", 0);
    /* Overlaps -> False collapses each run to one non-overlapping match. */
    assert_eval_eq("StringPosition[\"a12b345\", RegularExpression[\"[0-9]+\"], "
                   "Overlaps -> False]",
                   "{{2, 3}, {5, 7}}", 0);
}

/* A list of subject strings threads to a list of results. */
static void test_subject_threading(void) {
    assert_eval_eq("StringPosition[{\"abc\", \"XYZabc\"}, \"a\"]",
                   "{{{1, 1}}, {{4, 4}}}", 0);
    /* Empty subject list threads to an empty list. */
    assert_eval_eq("StringPosition[{}, \"a\"]", "{}", 0);
    /* Non-string elements are passed through unchanged. */
    assert_eval_eq("StringPosition[{\"aa\", 7}, \"a\"]",
                   "{{{1, 1}, {2, 2}}, 7}", 0);
}

/* Compatibility with StringTake: the pair really selects the matched text. */
static void test_stringtake_roundtrip(void) {
    assert_eval_eq("StringTake[\"XYZabc\", StringPosition[\"XYZabc\", \"XYZ\"][[1]]]",
                   "\"XYZ\"", 0);
}

/* Empty string / empty match edges. */
static void test_edges(void) {
    assert_eval_eq("StringPosition[\"\", \"a\"]", "{}", 0);
    assert_eval_eq("StringPosition[\"aaa\", \"a\", 0]",
                   "{{1, 1}, {2, 2}, {3, 3}}", 0);   /* n <= 0 imposes no cap */
}

/* Options[StringPosition] exposes defaults; SetOptions changes the behaviour. */
static void test_options(void) {
    assert_eval_eq("Options[StringPosition]",
                   "{IgnoreCase -> False, Overlaps -> True}", 0);
    /* SetOptions redefines the default used when no explicit option is given. */
    assert_eval_eq("SetOptions[StringPosition, Overlaps -> False]",
                   "{IgnoreCase -> False, Overlaps -> False}", 0);
    assert_eval_eq("StringPosition[\"AAAAA\", \"AA\"]", "{{1, 2}, {3, 4}}", 0);
    /* An explicit option still overrides the changed default. */
    assert_eval_eq("StringPosition[\"AAAAA\", \"AA\", Overlaps -> True]",
                   "{{1, 2}, {2, 3}, {3, 4}, {4, 5}}", 0);
    /* Restore the default so test order does not matter. */
    assert_eval_eq("SetOptions[StringPosition, Overlaps -> True]",
                   "{IgnoreCase -> False, Overlaps -> True}", 0);
    assert_eval_eq("StringPosition[\"AAAAA\", \"AA\"]",
                   "{{1, 2}, {2, 3}, {3, 4}, {4, 5}}", 0);
}

/* Invalid arguments leave the call unevaluated (message emitted). */
static void test_unevaluated(void) {
    assert_eval_eq("StringPosition[xyz, \"a\"]", "StringPosition[xyz, \"a\"]", 0);
    assert_eval_eq("StringPosition[]", "StringPosition[]", 0);
    /* A single argument is below the required arity. */
    assert_eval_eq("StringPosition[\"abc\"]", "StringPosition[\"abc\"]", 0);
    /* Non-integer third argument. */
    assert_eval_eq("StringPosition[\"abc\", \"a\", x]",
                   "StringPosition[\"abc\", \"a\", x]", 0);
}

#endif /* USE_REGEX */

int main(void) {
    symtab_init();
    core_init();

#ifdef USE_REGEX
    TEST(test_literal);
    TEST(test_count);
    TEST(test_repeated_pattern);
    TEST(test_overlaps_literal);
    TEST(test_pattern_list);
    TEST(test_ignorecase);
    TEST(test_blanksequence);
    TEST(test_list_vs_alternatives);
    TEST(test_regex_and_classes);
    TEST(test_subject_threading);
    TEST(test_stringtake_roundtrip);
    TEST(test_edges);
    TEST(test_options);
    TEST(test_unevaluated);
#else
    printf("USE_REGEX not defined; skipping StringPosition tests\n");
#endif

    printf("All StringPosition tests passed!\n");
    return 0;
}
