/*
 * test_stringcontainsq.c - unit tests for the substring-predicate family:
 * StringContainsQ, StringFreeQ, StringStartsQ and StringEndsQ.
 *
 * All four are backed by the PCRE2 regex engine, so the whole suite is a no-op
 * when built without PCRE2 (USE_REGEX undefined), mirroring test_stringcount.c.
 * Mathilda's lexer drops one backslash per escape, so a regex backslash is
 * written "\\" in Mathilda source and therefore "\\\\" in this C source.
 *
 * The four share one core (src/strings/regex/stringcontainsq.c) and are related
 * by the Wolfram Language identities
 *
 *   StringFreeQ  [s, p] == !StringContainsQ[s, p]
 *   StringStartsQ[s, p] ==  StringContainsQ[s, StartOfString ~~ p]
 *   StringEndsQ  [s, p] ==  StringContainsQ[s, p ~~ EndOfString]
 *   StringContainsQ[s, p] == StringMatchQ[s, ___ ~~ p ~~ ___]
 *
 * test_duality below asserts those identities directly over a table of
 * subjects and patterns; it is the regression net that keeps the four in
 * agreement when the shared core changes.
 *
 * StringStartsQ / StringEndsQ anchor by wrapping the pattern in a synthesised
 * StringExpression[StartOfString, patt] rather than by touching the shared rule
 * builder.  test_anchor_grouping guards the subtle failure that design could
 * have: an alternation-bearing pattern must be parenthesised, or `\Aa|b` would
 * make StringStartsQ["zb", RegularExpression["a|b"]] wrongly True.
 */

#include "test_utils.h"

extern void symtab_init(void);
extern void core_init(void);

#ifdef USE_REGEX

/* ========================= literal substrings ========================= */

static void test_contains_literal(void) {
    assert_eval_eq("StringContainsQ[\"abcd\", \"b\"]", "True", 0);
    assert_eval_eq("StringContainsQ[\"abcd\", \"a\"]", "True", 0);   /* at the start */
    assert_eval_eq("StringContainsQ[\"abcd\", \"d\"]", "True", 0);   /* at the end */
    assert_eval_eq("StringContainsQ[\"abcd\", \"bc\"]", "True", 0);
    assert_eval_eq("StringContainsQ[\"abcd\", \"abcd\"]", "True", 0); /* whole subject */
    assert_eval_eq("StringContainsQ[\"abcd\", \"z\"]", "False", 0);
    assert_eval_eq("StringContainsQ[\"abcd\", \"abcde\"]", "False", 0);
    /* Empty subject / empty pattern. */
    assert_eval_eq("StringContainsQ[\"\", \"a\"]", "False", 0);
    assert_eval_eq("StringContainsQ[\"abc\", \"\"]", "True", 0);
    assert_eval_eq("StringContainsQ[\"\", \"\"]", "True", 0);
    /* A bare string is a LITERAL pattern: "." matches only a dot. */
    assert_eval_eq("StringContainsQ[\"a.b\", \".\"]", "True", 0);
    assert_eval_eq("StringContainsQ[\"ab\", \".\"]", "False", 0);
}

static void test_free_literal(void) {
    assert_eval_eq("StringFreeQ[\"abcd\", \"a\"]", "False", 0);
    assert_eval_eq("StringFreeQ[\"abcd\", \"z\"]", "True", 0);
    assert_eval_eq("StringFreeQ[\"\", \"a\"]", "True", 0);
    assert_eval_eq("StringFreeQ[\"\", \"\"]", "False", 0);
    assert_eval_eq("StringFreeQ[\"ab\", \".\"]", "True", 0);
    assert_eval_eq("StringFreeQ[\"commit\", \"co\"]", "False", 0);
    assert_eval_eq("StringContainsQ[\"commit\", \"co\"]", "True", 0);
}

/* ====================== general string expressions ====================== */

static void test_string_patterns(void) {
    assert_eval_eq("StringFreeQ[\"bcde\", \"b\" ~~ __ ~~ \"e\"]", "False", 0);
    assert_eval_eq("StringFreeQ[\"bcde\", \"c\" ~~ __ ~~ \"t\"]", "True", 0);
    assert_eval_eq("StringContainsQ[\"bcde\", \"b\" ~~ __ ~~ \"e\"]", "True", 0);
    assert_eval_eq("StringContainsQ[\"bcde\", \"c\" ~~ __ ~~ \"t\"]", "False", 0);

    /* Character classes. */
    assert_eval_eq("StringFreeQ[\"a1 and a2\", DigitCharacter ..]", "False", 0);
    assert_eval_eq("StringFreeQ[\"abc\", DigitCharacter ..]", "True", 0);
    assert_eval_eq("StringContainsQ[\"a1 and a2\", DigitCharacter ..]", "True", 0);
    assert_eval_eq("StringContainsQ[\"hello world\", Whitespace]", "True", 0);
    assert_eval_eq("StringContainsQ[\"helloworld\", Whitespace]", "False", 0);
    assert_eval_eq("StringContainsQ[\"ab3\", LetterCharacter ~~ DigitCharacter]", "True", 0);
    assert_eval_eq("StringContainsQ[\"x 12.5\", NumberString]", "True", 0);

    /* Except: a single character that does not begin an "a"-match. */
    assert_eval_eq("StringContainsQ[\"abc\", Except[\"a\"]]", "True", 0);
    assert_eval_eq("StringContainsQ[\"aaa\", Except[\"a\"]]", "False", 0);

    /* A repeated named pattern is a backreference: x_ ~~ x_ is two equal
     * adjacent characters. */
    assert_eval_eq("StringFreeQ[\"abcade\", x_ ~~ x_]", "True", 0);
    assert_eval_eq("StringContainsQ[\"abcade\", x_ ~~ x_]", "False", 0);
    assert_eval_eq("StringContainsQ[\"abbc\", x_ ~~ x_]", "True", 0);
}

static void test_regex_patterns(void) {
    assert_eval_eq("StringFreeQ[\"abcde\", RegularExpression[\"b.*d\"]]", "False", 0);
    assert_eval_eq("StringContainsQ[\"abcde\", RegularExpression[\"b.*d\"]]", "True", 0);
    assert_eval_eq("StringContainsQ[\"abc\", RegularExpression[\"\\\\d+\"]]", "False", 0);
    assert_eval_eq("StringContainsQ[\"a1\", RegularExpression[\"\\\\d+\"]]", "True", 0);

    /* Regular expressions mixed with string patterns. */
    assert_eval_eq("StringFreeQ[\"bac 123\", "
                   "RegularExpression[\"a.*\"] ~~ DigitCharacter ..]", "False", 0);
    assert_eval_eq("StringContainsQ[\"bac 123\", "
                   "RegularExpression[\"a.*\"] ~~ DigitCharacter ..]", "True", 0);

    /* An alternation-bearing RegularExpression must be grouped when it is
     * concatenated, or `a|bz` would match a bare "a". */
    assert_eval_eq("StringContainsQ[\"bz\", RegularExpression[\"a|b\"] ~~ \"z\"]", "True", 0);
    assert_eval_eq("StringContainsQ[\"zzz\", RegularExpression[\"a|b\"] ~~ \"z\"]", "False", 0);
}

/* ========================== lists of patterns ========================== */

static void test_pattern_list(void) {
    assert_eval_eq("StringFreeQ[\"abcdabcdcd\", {\"abc\", \"cd\"}]", "False", 0);
    assert_eval_eq("StringFreeQ[\"abcdabcdcd\", \"abc\" | \"cd\"]", "False", 0);
    assert_eval_eq("StringContainsQ[\"abcdabcdcd\", {\"abc\", \"cd\"}]", "True", 0);
    assert_eval_eq("StringContainsQ[\"abcdabcdcd\", \"abc\" | \"cd\"]", "True", 0);

    /* A list of patterns is Alternatives: any one of them suffices. */
    assert_eval_eq("StringFreeQ[\"cad\", {\"a\", \"b\"}]", "False", 0);
    assert_eval_eq("StringFreeQ[\"cad\", \"a\" | \"b\"]", "False", 0);
    assert_eval_eq("StringContainsQ[\"cad\", {\"a\", \"b\"}]", "True", 0);
    assert_eval_eq("StringContainsQ[\"cad\", {\"q\", \"d\"}]", "True", 0);  /* second one */
    assert_eval_eq("StringFreeQ[\"xyz\", {\"a\", \"b\"}]", "True", 0);

    /* The list and the Alternatives spelling must agree everywhere. */
    assert_eval_eq("StringContainsQ[\"cad\", {\"a\", \"b\"}] == "
                   "StringContainsQ[\"cad\", \"a\" | \"b\"]", "True", 0);

    /* An empty list of patterns matches nothing. */
    assert_eval_eq("StringContainsQ[\"abc\", {}]", "False", 0);
    assert_eval_eq("StringFreeQ[\"abc\", {}]", "True", 0);
    assert_eval_eq("StringStartsQ[\"abc\", {}]", "False", 0);
    assert_eval_eq("StringEndsQ[\"abc\", {}]", "False", 0);
    assert_eval_eq("StringFreeQ[{\"a\", \"b\"}, {}]", "{True, True}", 0);
}

/* ====================== threading over subject lists ==================== */

static void test_subject_threading(void) {
    assert_eval_eq("StringFreeQ[{\"a\", \"b\", \"ab\", \"abcd\", \"bcde\"}, \"a\"]",
                   "{False, True, False, False, True}", 0);
    assert_eval_eq("StringContainsQ[{\"a\", \"b\", \"ab\", \"abcd\", \"bcde\"}, \"a\"]",
                   "{True, False, True, True, False}", 0);
    assert_eval_eq("StringFreeQ[{\"aba\", \"bcd\", \"cea\"}, \"c\"]",
                   "{True, False, False}", 0);
    assert_eval_eq("StringContainsQ[{\"aba\", \"bcd\", \"cea\"}, \"a\"]",
                   "{True, False, True}", 0);
    assert_eval_eq("StringFreeQ[{\"ability\", \"listable\", \"argument\"}, "
                   "\"a\" ~~ __ ~~ \"t\" ~~ ___]", "{False, True, False}", 0);
    assert_eval_eq("StringContainsQ[{\"ability\", \"listable\", \"argument\"}, "
                   "\"a\" ~~ __ ~~ \"t\" ~~ ___]", "{True, False, True}", 0);

    /* An empty subject list threads to an empty list. */
    assert_eval_eq("StringContainsQ[{}, \"a\"]", "{}", 0);
    assert_eval_eq("StringFreeQ[{}, \"a\"]", "{}", 0);
    assert_eval_eq("StringStartsQ[{}, \"a\"]", "{}", 0);

    /* Threading is over the SUBJECT only: a list pattern stays alternatives. */
    assert_eval_eq("StringContainsQ[{\"ax\", \"by\", \"cz\"}, {\"a\", \"b\"}]",
                   "{True, True, False}", 0);
}

/* ============================== IgnoreCase ============================== */

static void test_ignorecase(void) {
    assert_eval_eq("StringFreeQ[\"BACCD\", \"ac\", IgnoreCase -> False]", "True", 0);
    assert_eval_eq("StringFreeQ[\"BACCD\", \"ac\", IgnoreCase -> True]", "False", 0);
    assert_eval_eq("StringContainsQ[\"abcd\", \"BC\", IgnoreCase -> False]", "False", 0);
    assert_eval_eq("StringContainsQ[\"abcd\", \"BC\", IgnoreCase -> True]", "True", 0);

    /* False is the default. */
    assert_eval_eq("StringContainsQ[\"abcd\", \"BC\"]", "False", 0);

    /* IgnoreCase also applies to the anchored faces and to list patterns. */
    assert_eval_eq("StringStartsQ[\"Commit\", \"co\"]", "False", 0);
    assert_eval_eq("StringStartsQ[\"Commit\", \"co\", IgnoreCase -> True]", "True", 0);
    assert_eval_eq("StringEndsQ[\"commIT\", \"it\"]", "False", 0);
    assert_eval_eq("StringEndsQ[\"commIT\", \"it\", IgnoreCase -> True]", "True", 0);
    assert_eval_eq("StringContainsQ[\"ABC\", {\"x\", \"b\"}, IgnoreCase -> True]", "True", 0);
    assert_eval_eq("StringContainsQ[{\"AB\", \"cd\"}, \"ab\", IgnoreCase -> True]",
                   "{True, False}", 0);
}

static void test_options(void) {
    assert_eval_eq("Options[StringContainsQ]", "{IgnoreCase -> False}", 0);
    assert_eval_eq("Options[StringFreeQ]", "{IgnoreCase -> False}", 0);
    assert_eval_eq("Options[StringStartsQ]", "{IgnoreCase -> False}", 0);
    assert_eval_eq("Options[StringEndsQ]", "{IgnoreCase -> False}", 0);

    /* SetOptions changes the default; an explicit option still overrides it. */
    assert_eval_eq("SetOptions[StringContainsQ, IgnoreCase -> True]; "
                   "StringContainsQ[\"abcd\", \"BC\"]", "True", 0);
    assert_eval_eq("StringContainsQ[\"abcd\", \"BC\", IgnoreCase -> False]", "False", 0);
    assert_eval_eq("SetOptions[StringContainsQ, IgnoreCase -> False]; "
                   "StringContainsQ[\"abcd\", \"BC\"]", "False", 0);
}

/* ============================ operator form ============================ */

static void test_operator_form(void) {
    assert_eval_eq("StringContainsQ[\"a\"][\"abc\"]", "True", 0);
    assert_eval_eq("StringContainsQ[\"a\"][\"xyz\"]", "False", 0);
    assert_eval_eq("StringFreeQ[\"a\"][\"bcd\"]", "True", 0);
    assert_eval_eq("StringStartsQ[\"co\"][\"commit\"]", "True", 0);
    assert_eval_eq("StringEndsQ[\"it\"][\"commit\"]", "True", 0);

    /* A curried option is carried through to the eventual application. */
    assert_eval_eq("StringFreeQ[\"ac\", IgnoreCase -> True][\"BACCD\"]", "False", 0);
    assert_eval_eq("StringFreeQ[\"ac\", IgnoreCase -> False][\"BACCD\"]", "True", 0);

    /* The operator threads over a subject list just like the two-argument form. */
    assert_eval_eq("StringContainsQ[\"a\"][{\"ab\", \"cd\"}]", "{True, False}", 0);

    /* Usable wherever a predicate is expected. */
    assert_eval_eq("StringContainsQ[\"a\"] /@ {\"abc\", \"xyz\"}", "{True, False}", 0);
    assert_eval_eq("Select[{\"abc\", \"xyz\", \"bat\"}, StringContainsQ[\"a\"]]",
                   "{\"abc\", \"bat\"}", 0);
    assert_eval_eq("Select[{\"abc\", \"xyz\", \"bat\"}, StringFreeQ[\"a\"]]",
                   "{\"xyz\"}", 0);
}

/* ======================= StringStartsQ / StringEndsQ ==================== */

static void test_startsq_endsq(void) {
    assert_eval_eq("StringStartsQ[\"commit\", \"co\"]", "True", 0);
    assert_eval_eq("StringStartsQ[\"commit\", \"om\"]", "False", 0);
    assert_eval_eq("StringStartsQ[\"commit\", \"it\"]", "False", 0);
    assert_eval_eq("StringEndsQ[\"commit\", \"it\"]", "True", 0);
    assert_eval_eq("StringEndsQ[\"commit\", \"mi\"]", "False", 0);
    assert_eval_eq("StringEndsQ[\"commit\", \"co\"]", "False", 0);

    /* The whole string is both a prefix and a suffix; "" is too. */
    assert_eval_eq("StringStartsQ[\"abc\", \"abc\"]", "True", 0);
    assert_eval_eq("StringEndsQ[\"abc\", \"abc\"]", "True", 0);
    assert_eval_eq("StringStartsQ[\"abc\", \"\"]", "True", 0);
    assert_eval_eq("StringEndsQ[\"abc\", \"\"]", "True", 0);
    assert_eval_eq("StringStartsQ[\"abc\", \"abcd\"]", "False", 0);
    assert_eval_eq("StringEndsQ[\"\", \"a\"]", "False", 0);

    /* Lists of patterns and lists of subjects. */
    assert_eval_eq("StringStartsQ[\"abc\", {\"x\", \"a\"}]", "True", 0);
    assert_eval_eq("StringEndsQ[\"abc\", {\"x\", \"c\"}]", "True", 0);
    assert_eval_eq("StringStartsQ[\"abc\", {\"x\", \"y\"}]", "False", 0);
    assert_eval_eq("StringStartsQ[{\"apple\", \"banana\"}, \"a\"]", "{True, False}", 0);
    assert_eval_eq("StringEndsQ[{\"apple\", \"banana\"}, \"a\"]", "{False, True}", 0);

    /* General string patterns, not just literals. */
    assert_eval_eq("StringStartsQ[\"a123\", LetterCharacter ~~ DigitCharacter ..]",
                   "True", 0);
    assert_eval_eq("StringEndsQ[\"a123\", DigitCharacter ..]", "True", 0);
    assert_eval_eq("StringEndsQ[\"123a\", DigitCharacter ..]", "False", 0);
}

/*
 * The anchor must bind to the whole pattern, not to its first alternative.
 * These four cases fail loudly if the synthesised StringExpression wrapper ever
 * stops being parenthesised (`\Aa|b` instead of `\A(?:a|b)`).
 */
static void test_anchor_grouping(void) {
    assert_eval_eq("StringStartsQ[\"bz\", RegularExpression[\"a|b\"]]", "True", 0);
    assert_eval_eq("StringStartsQ[\"zb\", RegularExpression[\"a|b\"]]", "False", 0);
    assert_eval_eq("StringEndsQ[\"za\", RegularExpression[\"a|b\"]]", "True", 0);
    assert_eval_eq("StringEndsQ[\"az\", RegularExpression[\"a|b\"]]", "False", 0);

    /* Likewise the anchor must not be satisfied by a later start offset. */
    assert_eval_eq("StringStartsQ[\"xa123\", LetterCharacter ~~ DigitCharacter ..]",
                   "False", 0);
    assert_eval_eq("StringStartsQ[\"abcabc\", \"bc\"]", "False", 0);
    assert_eval_eq("StringEndsQ[\"abcabc\", \"ab\"]", "False", 0);
}

/* ================ StartOfString / EndOfString in the translator ========= */

static void test_anchors_in_patterns(void) {
    assert_eval_eq("StringContainsQ[\"commit\", StartOfString ~~ \"co\"]", "True", 0);
    assert_eval_eq("StringContainsQ[\"commit\", StartOfString ~~ \"om\"]", "False", 0);
    assert_eval_eq("StringContainsQ[\"commit\", \"it\" ~~ EndOfString]", "True", 0);
    assert_eval_eq("StringContainsQ[\"commit\", \"mi\" ~~ EndOfString]", "False", 0);
    assert_eval_eq("StringFreeQ[\"commit\", StartOfString ~~ \"om\"]", "True", 0);

    /* The two heads are shared with the rest of the regex string family. */
    assert_eval_eq("StringMatchQ[\"abc\", StartOfString ~~ \"abc\" ~~ EndOfString]",
                   "True", 0);
    assert_eval_eq("StringCases[\"abcabc\", StartOfString ~~ \"abc\"]", "{\"abc\"}", 0);
    assert_eval_eq("StringCount[\"abcabc\", StartOfString ~~ \"abc\"]", "1", 0);
    assert_eval_eq("StringCount[\"abcabc\", \"abc\"]", "2", 0);
    assert_eval_eq("StringPosition[\"abcabc\", \"abc\" ~~ EndOfString]", "{{4, 6}}", 0);
    assert_eval_eq("StringReplace[\"abcabc\", StartOfString ~~ \"abc\" -> \"X\"]",
                   "\"Xabc\"", 0);
}

/* ========================== the four identities ========================= */

static void test_duality(void) {
    /* StringFreeQ is exactly the negation of StringContainsQ, over a table of
     * subjects crossed with patterns (literal, absent, empty, and general). */
    assert_eval_eq(
        "And @@ Flatten[Table[StringFreeQ[s, p] == !StringContainsQ[s, p], "
        "{s, {\"abc\", \"\", \"xyz\", \"aabbcc\"}}, "
        "{p, {\"a\", \"z\", \"\", \"a\" ~~ __ ~~ \"c\", DigitCharacter}}]]",
        "True", 0);

    /* StringStartsQ[s, p] == StringContainsQ[s, StartOfString ~~ p]. */
    assert_eval_eq(
        "And @@ Flatten[Table[StringStartsQ[s, p] == "
        "StringContainsQ[s, StartOfString ~~ p], "
        "{s, {\"commit\", \"abc\", \"\"}}, {p, {\"co\", \"om\", \"a\", \"\"}}]]",
        "True", 0);

    /* StringEndsQ[s, p] == StringContainsQ[s, p ~~ EndOfString]. */
    assert_eval_eq(
        "And @@ Flatten[Table[StringEndsQ[s, p] == "
        "StringContainsQ[s, p ~~ EndOfString], "
        "{s, {\"commit\", \"abc\", \"\"}}, {p, {\"it\", \"mi\", \"c\", \"\"}}]]",
        "True", 0);

    /* StringContainsQ[s, p] == StringMatchQ[s, ___ ~~ p ~~ ___]. */
    assert_eval_eq(
        "And @@ Flatten[Table[StringContainsQ[s, p] == "
        "StringMatchQ[s, ___ ~~ p ~~ ___], "
        "{s, {\"commit\", \"abc\", \"xyz\"}}, {p, {\"om\", \"a\", \"q\"}}]]",
        "True", 0);

    /* And it agrees with the counting builtin it does not share code with. */
    assert_eval_eq(
        "And @@ Flatten[Table[StringContainsQ[s, p] == (StringCount[s, p] > 0), "
        "{s, {\"abcabc\", \"xyz\", \"\"}}, {p, {\"abc\", \"c\", \"q\"}}]]",
        "True", 0);
}

/* ============================= unevaluated ============================= */

static void test_unevaluated(void) {
    /* Zero arguments: an argt message, expression returned unchanged. */
    assert_eval_eq("StringFreeQ[]", "StringFreeQ[]", 0);
    assert_eval_eq("StringContainsQ[]", "StringContainsQ[]", 0);
    assert_eval_eq("StringStartsQ[]", "StringStartsQ[]", 0);
    assert_eval_eq("StringEndsQ[]", "StringEndsQ[]", 0);

    /* Too many positional arguments. */
    assert_eval_eq("StringContainsQ[\"a\", \"b\", \"c\"]",
                   "StringContainsQ[\"a\", \"b\", \"c\"]", 0);

    /* A non-string subject leaves the call alone rather than answering. */
    assert_eval_eq("StringFreeQ[5, \"a\"]", "StringFreeQ[5, \"a\"]", 0);
    assert_eval_eq("StringContainsQ[x, \"a\"]", "StringContainsQ[x, \"a\"]", 0);

    /* So does a subject list holding a non-string: a Boolean predicate must not
     * pass the offending element through as if it were a result. */
    assert_eval_eq("StringContainsQ[{\"a\", 7}, \"a\"]",
                   "StringContainsQ[{\"a\", 7}, \"a\"]", 0);

    /* An unsupported pattern is left for whatever might define it later. */
    assert_eval_eq("StringContainsQ[\"abc\", 5]", "StringContainsQ[\"abc\", 5]", 0);
    assert_eval_eq("StringEndsQ[\"abc\", f[x]]", "StringEndsQ[\"abc\", f[x]]", 0);

    /* A single argument is the operator form, not an error. */
    assert_eval_eq("StringContainsQ[\"a\"]", "StringContainsQ[#1, \"a\"] &", 0);
}

/* ============================== attributes ============================= */

static void test_attributes(void) {
    assert_eval_eq("Attributes[StringContainsQ]", "{Protected}", 0);
    assert_eval_eq("Attributes[StringFreeQ]", "{Protected}", 0);
    assert_eval_eq("Attributes[StringStartsQ]", "{Protected}", 0);
    assert_eval_eq("Attributes[StringEndsQ]", "{Protected}", 0);
    assert_eval_eq("Attributes[StartOfString]", "{Protected}", 0);
    assert_eval_eq("Attributes[EndOfString]", "{Protected}", 0);
}

#endif /* USE_REGEX */

int main(void) {
    symtab_init();
    core_init();

#ifdef USE_REGEX
    TEST(test_contains_literal);
    TEST(test_free_literal);
    TEST(test_string_patterns);
    TEST(test_regex_patterns);
    TEST(test_pattern_list);
    TEST(test_subject_threading);
    TEST(test_ignorecase);
    TEST(test_options);
    TEST(test_operator_form);
    TEST(test_startsq_endsq);
    TEST(test_anchor_grouping);
    TEST(test_anchors_in_patterns);
    TEST(test_duality);
    TEST(test_unevaluated);
    TEST(test_attributes);
#else
    printf("USE_REGEX not defined; skipping StringContainsQ tests\n");
#endif

    printf("All StringContainsQ / StringFreeQ / StringStartsQ / StringEndsQ tests passed!\n");
    return 0;
}
