/*
 * test_stringcount.c - unit tests for the StringCount builtin.
 *
 * StringCount is backed by the PCRE2 regex engine, so the whole suite is a
 * no-op when built without PCRE2 (USE_REGEX undefined), mirroring
 * test_stringposition.c.  Mathilda's lexer drops one backslash per escape, so a
 * regex backslash is written "\\" in Mathilda source and therefore "\\\\" in
 * this C source.
 *
 * StringCount shares its match enumeration (regex_scan, regex_common.c) with
 * StringCases and StringPosition; test_invariants below is the regression net
 * that keeps the three in agreement.
 */

#include "test_utils.h"

extern void symtab_init(void);
extern void core_init(void);

#ifdef USE_REGEX

/* ======================= literal substring counts ======================= */

static void test_literal(void) {
    assert_eval_eq("StringCount[\"the cat sat on the mat\", \"at\"]", "3", 0);
    assert_eval_eq("StringCount[\"abcabc\", \"a\"]", "2", 0);
    assert_eval_eq("StringCount[\"XYZabc\", \"XYZ\"]", "1", 0);
    /* Whole subject is one match. */
    assert_eval_eq("StringCount[\"XYZ\", \"XYZ\"]", "1", 0);
    /* No match. */
    assert_eval_eq("StringCount[\"abc\", \"z\"]", "0", 0);
    /* Empty subject. */
    assert_eval_eq("StringCount[\"\", \"a\"]", "0", 0);
    /* Single-character subject. */
    assert_eval_eq("StringCount[\"a\", \"a\"]", "1", 0);
    /* A bare string is a LITERAL pattern: "." matches only a dot. */
    assert_eval_eq("StringCount[\"a.b.c\", \".\"]", "2", 0);
}

/* ====================== general string expressions ====================== */

static void test_string_patterns(void) {
    assert_eval_eq("StringCount[\"a12b345\", DigitCharacter]", "5", 0);
    assert_eval_eq("StringCount[\"a12b345\", LetterCharacter]", "2", 0);
    assert_eval_eq("StringCount[\"a1 b2\", WordCharacter]", "4", 0);
    /* Whitespace matches a whole RUN of whitespace, so "a b  c" has two. */
    assert_eval_eq("StringCount[\"a b  c\", Whitespace]", "2", 0);
    assert_eval_eq("StringCount[\"x=1.5, y=-2\", NumberString]", "2", 0);
    assert_eval_eq("StringCount[\"abcde\", Except[\"a\"]]", "4", 0);
    /* StringExpression concatenation (the ~~ operator). */
    assert_eval_eq("StringCount[\"aaab aab ab\", \"a\" ~~ \"b\"]", "3", 0);
    /* Named pattern with a backreference: two equal adjacent characters. */
    assert_eval_eq("StringCount[\"AABBBAABABB\", x_ ~~ x_]", "4", 0);
    /* Greedy BlankSequence swallows the whole subject in one match; the
     * null-sequence form additionally matches empty at the end. */
    assert_eval_eq("StringCount[\"AAAA\", x__]", "1", 0);
    assert_eval_eq("StringCount[\"AAAA\", x___]", "2", 0);
}

static void test_regex_patterns(void) {
    assert_eval_eq("StringCount[\"a23b4222c63333d80\", RegularExpression[\"\\\\d+\"]]",
                   "4", 0);
    assert_eval_eq("StringCount[\"adefgh12c34\", RegularExpression[\"[a-e]+\"]]",
                   "2", 0);
    assert_eval_eq("StringCount[\"aabc1aaaagh2ade\", RegularExpression[\"a{2,3}\"]]",
                   "2", 0);
    /* Non-greedy quantifier. */
    assert_eval_eq("StringCount[\"abc1agh2cde\", RegularExpression[\"a.+?\\\\d\"]]",
                   "2", 0);
    /* Inline case-insensitivity modifier. */
    assert_eval_eq("StringCount[\"AaBbCc\", RegularExpression[\"(?i)[a-c]\"]]",
                   "6", 0);
    /* RegularExpression["."] matches any character, unlike the literal ".". */
    assert_eval_eq("StringCount[\"a1b2\", RegularExpression[\".\"]]", "4", 0);
}

/* ============================ pattern lists ============================ */

static void test_pattern_list(void) {
    assert_eval_eq("StringCount[\"a1b2\", {\"a\", DigitCharacter}]", "3", 0);
    assert_eval_eq("StringCount[\"ABAABBAABABB\", {\"ABA\", \"AA\"}]", "2", 0);
    assert_eval_eq("StringCount[\"ABAABBAABABB\", {\"ABA\", \"AA\"}, Overlaps -> True]",
                   "4", 0);
    /* At a given start the earlier list element wins, so {x__, "A"} counts the
     * one greedy x__ match, not the single "A" as well. */
    assert_eval_eq("StringCount[\"AAAA\", {x__, \"A\"}]", "1", 0);
    assert_eval_eq("StringCount[\"AAAA\", x__ | \"A\"]", "1", 0);
}

/* A Rule/RuleDelayed pattern counts matches of the LHS; the RHS is irrelevant
 * to a count and must not change the answer. */
static void test_rule_pattern(void) {
    assert_eval_eq("StringCount[\"a1b2\", DigitCharacter -> \"x\"]", "2", 0);
    assert_eval_eq("StringCount[\"a1b2\", DigitCharacter :> \"x\"]", "2", 0);
    assert_eval_eq("StringCount[\"a1b2\", DigitCharacter -> \"x\"] == "
                   "StringCount[\"a1b2\", DigitCharacter]", "True", 0);
}

/* ========================= subject threading ========================= */

static void test_subject_threading(void) {
    assert_eval_eq("StringCount[{\"a1\", \"b22\", \"ccc\"}, DigitCharacter]",
                   "{1, 2, 0}", 0);
    assert_eval_eq("StringCount[{\"abc\", \"XYZabc\"}, \"a\"]", "{1, 1}", 0);
    /* Empty subject list threads to an empty list. */
    assert_eval_eq("StringCount[{}, \"a\"]", "{}", 0);
    /* Non-string elements pass through unchanged. */
    assert_eval_eq("StringCount[{\"aa\", 7}, \"a\"]", "{2, 7}", 0);
    assert_eval_eq("StringCount[{\"a1\", 5, \"b22\"}, DigitCharacter]",
                   "{1, 5, 2}", 0);
}

/* ============================== Overlaps ============================== */

static void test_overlaps(void) {
    /* Default is Overlaps -> False: overlapping substrings are not separate. */
    assert_eval_eq("StringCount[\"AAAA\", \"AA\"]", "2", 0);
    assert_eval_eq("StringCount[\"AAAA\", \"AA\", Overlaps -> False]", "2", 0);
    /* True: overlaps count, one substring per start position. */
    assert_eval_eq("StringCount[\"AAAA\", \"AA\", Overlaps -> True]", "3", 0);
    assert_eval_eq("StringCount[\"AAAAA\", \"AA\", Overlaps -> True]", "4", 0);
    /* All: every matching substring at every start.  For a fixed-length
     * pattern that agrees with True; for a variable-length one it does not. */
    assert_eval_eq("StringCount[\"AAAA\", \"AA\", Overlaps -> All]", "3", 0);
    assert_eval_eq("StringCount[\"AAAA\", x__, Overlaps -> All]", "10", 0);
    assert_eval_eq("StringCount[\"abc\", x__, Overlaps -> All]", "6", 0);
    assert_eval_eq("StringCount[\"aaa\", \"a\", Overlaps -> All]", "3", 0);
    /* Named-pattern backreference under each policy. */
    assert_eval_eq("StringCount[\"AABBBAABABB\", x_ ~~ x_, Overlaps -> True]",
                   "5", 0);
}

/* ============================= IgnoreCase ============================= */

static void test_ignorecase(void) {
    assert_eval_eq("StringCount[\"aAbB\", \"a\", IgnoreCase -> True]", "2", 0);
    assert_eval_eq("StringCount[\"aAbB\", \"a\", IgnoreCase -> False]", "1", 0);
    /* Default is IgnoreCase -> False. */
    assert_eval_eq("StringCount[\"aAbB\", \"a\"]", "1", 0);
    assert_eval_eq("StringCount[\"AaBbCc\", LetterCharacter, IgnoreCase -> True]",
                   "6", 0);
    /* Both options together, in either order. */
    assert_eval_eq("StringCount[\"AaAa\", \"a\", IgnoreCase -> True, Overlaps -> True]",
                   "4", 0);
    assert_eval_eq("StringCount[\"AaAa\", \"a\", Overlaps -> True, IgnoreCase -> True]",
                   "4", 0);
    /* RuleDelayed (:>) is accepted as an option form too. */
    assert_eval_eq("StringCount[\"aAbB\", \"a\", IgnoreCase :> True]", "2", 0);
    assert_eval_eq("StringCount[\"AAAA\", \"AA\", Overlaps :> True]", "3", 0);
}

/* ==================== Options / SetOptions defaults ==================== */

static void test_options(void) {
    assert_eval_eq("Options[StringCount]",
                   "{IgnoreCase -> False, Overlaps -> False}", 0);
    /* SetOptions redefines the default used when no explicit option is given. */
    assert_eval_eq("SetOptions[StringCount, Overlaps -> True]",
                   "{IgnoreCase -> False, Overlaps -> True}", 0);
    assert_eval_eq("StringCount[\"AAAA\", \"AA\"]", "3", 0);
    /* An explicit option still overrides the changed default. */
    assert_eval_eq("StringCount[\"AAAA\", \"AA\", Overlaps -> False]", "2", 0);
    /* Restore the default so test order does not matter. */
    assert_eval_eq("SetOptions[StringCount, Overlaps -> False]",
                   "{IgnoreCase -> False, Overlaps -> False}", 0);
    assert_eval_eq("StringCount[\"AAAA\", \"AA\"]", "2", 0);

    /* The same for IgnoreCase. */
    assert_eval_eq("SetOptions[StringCount, IgnoreCase -> True]",
                   "{IgnoreCase -> True, Overlaps -> False}", 0);
    assert_eval_eq("StringCount[\"aAbB\", \"a\"]", "2", 0);
    assert_eval_eq("SetOptions[StringCount, IgnoreCase -> False]",
                   "{IgnoreCase -> False, Overlaps -> False}", 0);
    assert_eval_eq("StringCount[\"aAbB\", \"a\"]", "1", 0);
}

/* ===================== cross-builtin consistency ===================== */

/*
 * StringCount, StringCases and StringPosition run the same scan, so the counts
 * must agree for every pattern shape and every Overlaps policy.  If the shared
 * scanner is ever changed for one of them, this is what catches it.
 */
static void test_invariants(void) {
    /* StringCount == Length[StringCases] at the default options. */
    assert_eval_eq("StringCount[\"a12b345\", DigitCharacter] == "
                   "Length[StringCases[\"a12b345\", DigitCharacter]]", "True", 0);
    assert_eval_eq("StringCount[\"the cat sat on the mat\", \"at\"] == "
                   "Length[StringCases[\"the cat sat on the mat\", \"at\"]]",
                   "True", 0);
    assert_eval_eq("StringCount[\"ABAABBAABABB\", {\"ABA\", \"AA\"}] == "
                   "Length[StringCases[\"ABAABBAABABB\", {\"ABA\", \"AA\"}]]",
                   "True", 0);

    /* ... and under each explicit Overlaps setting. */
    assert_eval_eq("StringCount[\"AAAA\", \"AA\", Overlaps -> False] == "
                   "Length[StringCases[\"AAAA\", \"AA\", Overlaps -> False]]",
                   "True", 0);
    assert_eval_eq("StringCount[\"AAAA\", \"AA\", Overlaps -> True] == "
                   "Length[StringCases[\"AAAA\", \"AA\", Overlaps -> True]]",
                   "True", 0);
    assert_eval_eq("StringCount[\"AAAA\", x__, Overlaps -> All] == "
                   "Length[StringCases[\"AAAA\", x__, Overlaps -> All]]",
                   "True", 0);

    /* StringCount == Length[StringPosition] for each Overlaps policy. */
    assert_eval_eq("StringCount[\"AAAAA\", \"AA\", Overlaps -> False] == "
                   "Length[StringPosition[\"AAAAA\", \"AA\", Overlaps -> False]]",
                   "True", 0);
    assert_eval_eq("StringCount[\"AAAAA\", \"AA\", Overlaps -> True] == "
                   "Length[StringPosition[\"AAAAA\", \"AA\", Overlaps -> True]]",
                   "True", 0);
    assert_eval_eq("StringCount[\"AAAA\", x__, Overlaps -> All] == "
                   "Length[StringPosition[\"AAAA\", x__, Overlaps -> All]]",
                   "True", 0);
    assert_eval_eq("StringCount[\"ABAABBAABABB\", {\"ABA\", \"AA\"}, Overlaps -> All] == "
                   "Length[StringPosition[\"ABAABBAABABB\", {\"ABA\", \"AA\"}, "
                   "Overlaps -> All]]", "True", 0);
    assert_eval_eq("StringCount[\"AABBBAABABB\", x_ ~~ x_, Overlaps -> True] == "
                   "Length[StringPosition[\"AABBBAABABB\", x_ ~~ x_, "
                   "Overlaps -> True]]", "True", 0);
    /* IgnoreCase agrees across all three as well. */
    assert_eval_eq("StringCount[\"aAbaBabB\", \"a\", IgnoreCase -> True] == "
                   "Length[StringPosition[\"aAbaBabB\", \"a\", IgnoreCase -> True, "
                   "Overlaps -> False]]", "True", 0);
}

/* ============================ zero-width ============================ */

/*
 * A pattern that can match the empty string yields one zero-width match at each
 * position the scan reaches, including one past the last character.  These pin
 * the measured behaviour so a future scanner change cannot silently alter it.
 */
static void test_zero_width(void) {
    /* "a*" matches "" at 0, "aa" at 1, then "" at the end. */
    assert_eval_eq("StringCount[\"baa\", RegularExpression[\"a*\"]]", "3", 0);
    assert_eval_eq("StringCases[\"baa\", RegularExpression[\"a*\"]]",
                   "{\"\", \"aa\", \"\"}", 0);
    /* The empty pattern matches at every position plus the end. */
    assert_eval_eq("StringCount[\"abc\", \"\"]", "4", 0);
    assert_eval_eq("StringCount[\"\", \"\"]", "1", 0);
    /* Consistent with StringCases, and the scan always terminates. */
    assert_eval_eq("StringCount[\"abc\", \"\"] == Length[StringCases[\"abc\", \"\"]]",
                   "True", 0);
}

/* ====================== unevaluated / error paths ====================== */

static void test_unevaluated(void) {
    /* Non-string, non-list subject leaves the call unevaluated. */
    assert_eval_eq("StringCount[123, \"a\"]", "StringCount[123, \"a\"]", 0);
    assert_eval_eq("StringCount[xyz, \"a\"]", "StringCount[xyz, \"a\"]", 0);
    /* Unsupported pattern leaves the call unevaluated. */
    assert_eval_eq("StringCount[\"abc\", 5]", "StringCount[\"abc\", 5]", 0);
    /* Wrong arity emits StringCount::argrx and leaves the call unevaluated. */
    assert_eval_eq("StringCount[]", "StringCount[]", 0);
    assert_eval_eq("StringCount[\"abc\"]", "StringCount[\"abc\"]", 0);
}

/* ================= StringCases gained the same options ================= */

static void test_stringcases_options(void) {
    assert_eval_eq("Options[StringCases]",
                   "{IgnoreCase -> False, Overlaps -> False}", 0);
    assert_eval_eq("StringCases[\"AAAA\", \"AA\", Overlaps -> True]",
                   "{\"AA\", \"AA\", \"AA\"}", 0);
    assert_eval_eq("StringCases[\"AAA\", x__, Overlaps -> All]",
                   "{\"AAA\", \"AA\", \"A\", \"AA\", \"A\", \"A\"}", 0);
    assert_eval_eq("StringCases[\"aAbB\", \"a\", IgnoreCase -> True]",
                   "{\"a\", \"A\"}", 0);
    /* $n expansion still works when the scan carries capture offsets, both in
     * the overlapping and the all-substrings policies. */
    assert_eval_eq("StringCases[\"AAAA\", \"AA\" -> \"$0!\", Overlaps -> True]",
                   "{\"AA!\", \"AA!\", \"AA!\"}", 0);
    assert_eval_eq("StringCases[\"a1b22\", "
                   "RegularExpression[\"([a-z])(\\\\d+)\"] -> \"$2$1\"]",
                   "{\"1a\", \"22b\"}", 0);
    /* A rule whose LHS folds case still replaces. */
    assert_eval_eq("StringCases[\"ABC\", \"b\" -> \"X\", IgnoreCase -> True]",
                   "{\"X\"}", 0);
}

#endif /* USE_REGEX */

int main(void) {
    symtab_init();
    core_init();

#ifdef USE_REGEX
    TEST(test_literal);
    TEST(test_string_patterns);
    TEST(test_regex_patterns);
    TEST(test_pattern_list);
    TEST(test_rule_pattern);
    TEST(test_subject_threading);
    TEST(test_overlaps);
    TEST(test_ignorecase);
    TEST(test_options);
    TEST(test_invariants);
    TEST(test_zero_width);
    TEST(test_unevaluated);
    TEST(test_stringcases_options);
#else
    printf("USE_REGEX not defined; skipping StringCount tests\n");
#endif

    printf("All StringCount tests passed!\n");
    return 0;
}
