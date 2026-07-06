/*
 * test_stringfns.c - regression tests for the regex-aware string functions
 * StringMatchQ / StringCases / StringReplace / StringSplit, driven directly
 * from the worked examples in the RegularExpression specification.
 *
 * NOTE on escaping: Mathilda's string lexer consumes one backslash from every
 * "\x" escape, so a regex backslash must be written as "\\" in Mathilda source
 * and therefore as four backslashes ("\\\\") in this C source.  Real newlines
 * are embedded as literal '\n' bytes (Mathilda has no \n string escape).
 *
 * The whole suite is a no-op when built without PCRE2 (USE_REGEX undefined).
 */

#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"

#ifdef USE_REGEX

/* ============================= StringCases ============================= */

static void test_cases_charclass(void) {
    assert_eval_eq("StringCases[\"adefgh12c34\", RegularExpression[\"[a-e]+\"]]",
                   "{\"ade\", \"c\"}", 0);
    assert_eval_eq("StringCases[\"a13b12c1da32efg\", RegularExpression[\"[ab]\"]]",
                   "{\"a\", \"b\", \"a\"}", 0);
    assert_eval_eq("StringCases[\"adefgh12c34\", RegularExpression[\"[a-e]\"]]",
                   "{\"a\", \"d\", \"e\", \"c\"}", 0);
    assert_eval_eq("StringCases[\"a13b12c17a32\", RegularExpression[\"[^a1]\"]]",
                   "{\"3\", \"b\", \"2\", \"c\", \"7\", \"3\", \"2\"}", 0);
    assert_eval_eq("StringCases[\"aaaa bbbb 1234\", RegularExpression[\"[a-z]+\"]]",
                   "{\"aaaa\", \"bbbb\"}", 0);
}

static void test_cases_quantifiers(void) {
    assert_eval_eq("StringCases[\"a23b4222c63333d80\", RegularExpression[\"\\\\d+\"]]",
                   "{\"23\", \"4222\", \"63333\", \"80\"}", 0);
    assert_eval_eq("StringCases[\"aabc1aaaagh2ade\", RegularExpression[\"a{2,3}\"]]",
                   "{\"aa\", \"aaa\"}", 0);
    /* non-greedy */
    assert_eval_eq("StringCases[\"abc1agh2cde\", RegularExpression[\"a.+?\\\\d\"]]",
                   "{\"abc1\", \"agh2\"}", 0);
}

static void test_cases_shorthand_classes(void) {
    assert_eval_eq("StringCases[\"a2322c63333d80\", RegularExpression[\"\\\\d\"]]",
                   "{\"2\", \"3\", \"2\", \"2\", \"6\", \"3\", \"3\", \"3\", \"3\", \"8\", \"0\"}", 0);
    assert_eval_eq("StringCases[\"a2322c63333d80\", RegularExpression[\"\\\\D\"]]",
                   "{\"a\", \"c\", \"d\"}", 0);
    assert_eval_eq("StringCases[\"a23b42c63,d80\", RegularExpression[\"\\\\w\"]]",
                   "{\"a\", \"2\", \"3\", \"b\", \"4\", \"2\", \"c\", \"6\", \"3\", \"d\", \"8\", \"0\"}", 0);
    assert_eval_eq("StringCases[\"a23b:42c63;d80\", RegularExpression[\"\\\\W\"]]",
                   "{\":\", \";\"}", 0);
}

static void test_cases_named_class(void) {
    assert_eval_eq("StringCases[\"AaBBccDDeefG\", RegularExpression[\"[[:upper:]]+\"]]",
                   "{\"A\", \"BB\", \"DD\", \"G\"}", 0);
}

static void test_cases_literal_and_dot(void) {
    /* A bare string is a LITERAL pattern (dot matches only a dot). */
    assert_eval_eq("StringCases[\"a.b.c\", \".\"]", "{\".\", \".\"}", 0);
    /* RegularExpression["."] matches any character. */
    assert_eval_eq("StringCases[\"a1b2\", RegularExpression[\".\"]]",
                   "{\"a\", \"1\", \"b\", \"2\"}", 0);
}

static void test_cases_inline_option(void) {
    assert_eval_eq("StringCases[\"AaBbCc\", RegularExpression[\"(?i)[a-c]\"]]",
                   "{\"A\", \"a\", \"B\", \"b\", \"C\", \"c\"}", 0);
}

static void test_cases_empty(void) {
    assert_eval_eq("StringCases[\"\", RegularExpression[\"\\\\d\"]]", "{}", 0);
}

static void test_cases_thread(void) {
    assert_eval_eq("StringCases[{\"a1\", \"b22\"}, RegularExpression[\"\\\\d+\"]]",
                   "{{\"1\"}, {\"22\"}}", 0);
}

/* ============================ StringMatchQ ============================= */

static void test_matchq_basic(void) {
    assert_eval_eq("StringMatchQ[\"12345\", RegularExpression[\"\\\\d+\"]]", "True", 0);
    assert_eval_eq("StringMatchQ[\"12a45\", RegularExpression[\"\\\\d+\"]]", "False", 0);
    /* whole-string anchoring: a partial match is not enough */
    assert_eval_eq("StringMatchQ[\"abc\", \"ab\"]", "False", 0);
    assert_eval_eq("StringMatchQ[\"abc\", \"abc\"]", "True", 0);
}

static void test_matchq_newline(void) {
    /* Real newline bytes; \\w|\\s covers letters, digits and whitespace. */
    assert_eval_eq("StringMatchQ[\"abcd\nefgh\n1234\", RegularExpression[\"(\\\\w|\\\\s)*\"]]",
                   "True", 0);
    assert_eval_eq("StringMatchQ[\"abcd\nefgh\n1234\", RegularExpression[\"(.*|\\\\s*)*\"]]",
                   "True", 0);
}

static void test_matchq_thread(void) {
    assert_eval_eq("StringMatchQ[{\"12\", \"x\"}, RegularExpression[\"\\\\d+\"]]",
                   "{True, False}", 0);
}

static void test_matchq_unevaluated(void) {
    /* Non-string subject leaves the call unevaluated. */
    assert_eval_eq("StringMatchQ[xyz, RegularExpression[\"a\"]]",
                   "StringMatchQ[xyz, RegularExpression[\"a\"]]", 0);
}

/* ============================ StringReplace =========================== */

static void test_replace_backref(void) {
    assert_eval_eq("StringReplace[\"a13b12c1da32efg\", RegularExpression[\"(\\\\d+)\"] -> \"[$1]\"]",
                   "\"a[13]b[12]c[1]da[32]efg\"", 0);
}

static void test_replace_zero_width(void) {
    /* \\b is a zero-width word boundary; :> keeps the RHS a constant. */
    assert_eval_eq("StringReplace[\"123 45 6 789\", RegularExpression[\"\\\\b\"] :> \"X\"]",
                   "\"X123X X45X X6X X789X\"", 0);
}

static void test_replace_literal(void) {
    assert_eval_eq("StringReplace[\"hello\", \"l\" -> \"L\"]", "\"heLLo\"", 0);
    assert_eval_eq("StringReplace[\"a+b\", \"+\" -> \"plus\"]", "\"aplusb\"", 0);
}

static void test_replace_multi_rule(void) {
    assert_eval_eq("StringReplace[\"abcde\", {\"a\" -> \"1\", \"c\" -> \"3\"}]",
                   "\"1b3de\"", 0);
}

/* ============================= StringSplit ============================ */

static void test_split_delimiter(void) {
    assert_eval_eq("StringSplit[\"1.23, 4.56  7.89\", RegularExpression[\"(\\\\s|,)+\"]]",
                   "{\"1.23\", \"4.56\", \"7.89\"}", 0);
    assert_eval_eq("StringSplit[\"a,b,c\", \",\"]", "{\"a\", \"b\", \"c\"}", 0);
    /* Literal delimiter: the dot is not a metacharacter here. */
    assert_eval_eq("StringSplit[\"a.b.c\", \".\"]", "{\"a\", \"b\", \"c\"}", 0);
}

static void test_split_zero_width(void) {
    /* Zero-width delimiters at line starts / ends; result printed via InputForm.
     * Mathilda prints newlines raw, so the expected strings embed real '\n'. */
    assert_eval_eq("InputForm[StringSplit[\"line1\nline2\nline3\", RegularExpression[\"(?m)^\"]]]",
                   "{\"line1\n\", \"line2\n\", \"line3\"}", 0);
    assert_eval_eq("InputForm[StringSplit[\"line1\nline2\nline3\", RegularExpression[\"(?m)$\"]]]",
                   "{\"line1\", \"\nline2\", \"\nline3\"}", 0);
}

static void test_split_no_match(void) {
    /* Delimiter absent -> single piece; empty subject -> empty list. */
    assert_eval_eq("StringSplit[\"abc\", RegularExpression[\"x\"]]", "{\"abc\"}", 0);
    assert_eval_eq("StringSplit[\"\", RegularExpression[\"x\"]]", "{}", 0);
}

#endif /* USE_REGEX */

int main(void) {
    symtab_init();
    core_init();

#ifdef USE_REGEX
    TEST(test_cases_charclass);
    TEST(test_cases_quantifiers);
    TEST(test_cases_shorthand_classes);
    TEST(test_cases_named_class);
    TEST(test_cases_literal_and_dot);
    TEST(test_cases_inline_option);
    TEST(test_cases_empty);
    TEST(test_cases_thread);

    TEST(test_matchq_basic);
    TEST(test_matchq_newline);
    TEST(test_matchq_thread);
    TEST(test_matchq_unevaluated);

    TEST(test_replace_backref);
    TEST(test_replace_zero_width);
    TEST(test_replace_literal);
    TEST(test_replace_multi_rule);

    TEST(test_split_delimiter);
    TEST(test_split_zero_width);
    TEST(test_split_no_match);
#else
    printf("USE_REGEX not defined; skipping regex string-function tests\n");
#endif

    printf("All string-function tests passed!\n");
    return 0;
}
