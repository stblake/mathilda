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

static void test_split_whitespace_default(void) {
    /* One-argument form splits at runs of whitespace. */
    assert_eval_eq("StringSplit[\"a bbb  cccc aa   d\"]",
                   "{\"a\", \"bbb\", \"cccc\", \"aa\", \"d\"}", 0);
    assert_eval_eq("StringSplit[\"This is a sentence--which goes on.\"]",
                   "{\"This\", \"is\", \"a\", \"sentence--which\", \"goes\", \"on.\"}", 0);
}

static void test_split_literal_alternatives(void) {
    /* Literal delimiter (dot is not a metacharacter). */
    assert_eval_eq("StringSplit[\"192.168.0.1\", \".\"]",
                   "{\"192\", \"168\", \"0\", \"1\"}", 0);
    /* A list of delimiters and the equivalent Alternatives (|) form. */
    assert_eval_eq("StringSplit[\"a-b:c-d:e-f-g\", {\":\", \"-\"}]",
                   "{\"a\", \"b\", \"c\", \"d\", \"e\", \"f\", \"g\"}", 0);
    assert_eval_eq("StringSplit[\"a-b:c-d:e-f-g\", \":\" | \"-\"]",
                   "{\"a\", \"b\", \"c\", \"d\", \"e\", \"f\", \"g\"}", 0);
}

static void test_split_charclass(void) {
    assert_eval_eq("StringSplit[\"123  2.3  4  6\", WhitespaceCharacter ..]",
                   "{\"123\", \"2.3\", \"4\", \"6\"}", 0);
    assert_eval_eq("StringSplit[\"11a22b3\", _?LetterQ]",
                   "{\"11\", \"22\", \"3\"}", 0);
    assert_eval_eq("StringSplit[\"This is a sentence, which goes on.\", Except[WordCharacter] ..]",
                   "{\"This\", \"is\", \"a\", \"sentence\", \"which\", \"goes\", \"on\"}", 0);
}

static void test_split_stringexpression(void) {
    /* ~~ concatenates string patterns: whitespace, a digit, whitespace. */
    assert_eval_eq("StringSplit[\"primes: 2 two 3 three 5 five ...\", "
                   "Whitespace ~~ RegularExpression[\"\\\\d\"] ~~ Whitespace]",
                   "{\"primes:\", \"two\", \"three\", \"five ...\"}", 0);
}

static void test_split_rules(void) {
    /* Insert a value at each delimiter (patt -> val). */
    assert_eval_eq("StringSplit[\"a b::c d::e f g\", \"::\" -> \"--\"]",
                   "{\"a b\", \"--\", \"c d\", \"--\", \"e f g\"}", 0);
    /* :> with a named pattern keeps the delimiter itself. */
    assert_eval_eq("StringSplit[\"a--b c--d e\", x : \"--\" :> x]",
                   "{\"a\", \"--\", \"b c\", \"--\", \"d e\"}", 0);
    /* Equivalent to StringReplace when the pieces are rejoined. */
    assert_eval_eq("StringJoin[StringSplit[\"ab::c::d::ef\", \"::\" -> \"X\"]]",
                   "\"abXcXdXef\"", 0);
}

static void test_split_n_and_all(void) {
    /* At most n substrings; the remainder is the last piece. */
    assert_eval_eq("StringSplit[\"a b c d\", \" \", 2]", "{\"a\", \"b c d\"}", 0);
    /* All keeps the leading/trailing empty substrings. */
    assert_eval_eq("InputForm[StringSplit[\":a:b:c:\", \":\", All]]",
                   "{\"\", \"a\", \"b\", \"c\", \"\"}", 0);
    /* Interior empties are kept by default; boundary empties dropped. */
    assert_eval_eq("InputForm[StringSplit[\"a,,b\", \",\"]]",
                   "{\"a\", \"\", \"b\"}", 0);
    assert_eval_eq("InputForm[StringSplit[\":a:b:c:\", \":\"]]",
                   "{\"a\", \"b\", \"c\"}", 0);
}

static void test_split_null_and_ignorecase(void) {
    /* A null delimiter splits at every character. */
    assert_eval_eq("StringSplit[\"abcdefg\", \"\"]",
                   "{\"a\", \"b\", \"c\", \"d\", \"e\", \"f\", \"g\"}", 0);
    /* IgnoreCase treats c and C as equivalent. */
    assert_eval_eq("StringSplit[\"cat Cat hat CAT\", \"c\", IgnoreCase -> True]",
                   "{\"at \", \"at hat \", \"AT\"}", 0);
}

static void test_split_thread_and_regex(void) {
    assert_eval_eq("StringSplit[{\"a:b:c:d\", \"listable:element\"}, \":\"]",
                   "{{\"a\", \"b\", \"c\", \"d\"}, {\"listable\", \"element\"}}", 0);
    assert_eval_eq("StringSplit[\"A tree, an apple, four pears. And more: two sacks\", "
                   "RegularExpression[\"\\\\W+\"]]",
                   "{\"A\", \"tree\", \"an\", \"apple\", \"four\", \"pears\", "
                   "\"And\", \"more\", \"two\", \"sacks\"}", 0);
}

static void test_split_unevaluated(void) {
    /* A non-string subject leaves the call unevaluated. */
    assert_eval_eq("StringSplit[xyz, \":\"]", "StringSplit[xyz, \":\"]", 0);
}

/* The WL string-pattern vocabulary is shared: StringCases gains it too. */
static void test_cases_wl_pattern(void) {
    assert_eval_eq("StringCases[\"a bbb  cccc aa   d\", Except[WhitespaceCharacter] ..]",
                   "{\"a\", \"bbb\", \"cccc\", \"aa\", \"d\"}", 0);
    assert_eval_eq("StringReplace[\"a1b2c3\", DigitCharacter -> \"#\"]",
                   "\"a#b#c#\"", 0);
}

/* StringTrim - trims matching substrings from both ends of a string. */
static void test_trim_whitespace_default(void) {
    /* Default pattern trims leading/trailing whitespace runs. */
    assert_eval_eq("StringTrim[\"   aaa bbb ccc   \"]",
                   "\"aaa bbb ccc\"", 0);
    /* No trimmable ends: string is returned unchanged. */
    assert_eval_eq("StringTrim[\"hello\"]", "\"hello\"", 0);
    /* Tabs and real newlines count as whitespace (embedded as literal bytes). */
    assert_eval_eq("StringTrim[\"\t\n  hi there  \n\t\"]", "\"hi there\"", 0);
    /* An all-whitespace string trims to empty; empty input stays empty. */
    assert_eval_eq("StringTrim[\"   \"]", "\"\"", 0);
    assert_eval_eq("StringTrim[\"\"]", "\"\"", 0);
}

static void test_trim_pattern(void) {
    /* Multi-character alternative repeated at both ends: ("+" | "-") ... */
    assert_eval_eq("StringTrim[\"++++aaa bbb ccc----\", (\"+\" | \"-\") ...]",
                   "\"aaa bbb ccc\"", 0);
    /* A literal single-character pattern is stripped to a fixed point. */
    assert_eval_eq("StringTrim[\"xxabcxx\", \"x\"]", "\"abc\"", 0);
    /* Character-class Repeated strips digit runs from both ends. */
    assert_eval_eq("StringTrim[\"007bond007\", DigitCharacter ..]", "\"bond\"", 0);
    /* Interior matches are left untouched; only the ends are trimmed. */
    assert_eval_eq("StringTrim[\"xxaxbxx\", \"x\"]", "\"axb\"", 0);
}

static void test_trim_regex_front_only(void) {
    /* RegularExpression["^ *"] anchors to the start, so only the front trims. */
    assert_eval_eq("StringTrim[\"   aaa bbb ccc   \", RegularExpression[\"^ *\"]]",
                   "\"aaa bbb ccc   \"", 0);
}

static void test_trim_thread(void) {
    /* Threads over a list of strings; non-string elements pass through. */
    assert_eval_eq("StringTrim[{\"  a  \", \"  b  \"}]", "{\"a\", \"b\"}", 0);
    assert_eval_eq("StringTrim[{\"  a  \", 5, x}]", "{\"a\", 5, x}", 0);
}

static void test_trim_unevaluated(void) {
    /* A non-string subject leaves the call unevaluated. */
    assert_eval_eq("StringTrim[x]", "StringTrim[x]", 0);
    /* Zero arguments: unevaluated (a StringTrim::argt message is written). */
    assert_eval_eq("StringTrim[]", "StringTrim[]", 0);
}

/* ============================ StringExtract =========================== */

/* Single-level whitespace position specs: n, -n, {..}, span, All. */
static void test_extract_positions(void) {
    assert_eval_eq("StringExtract[\"a bbb  cccc aa   d\", 2]", "\"bbb\"", 0);
    assert_eval_eq("StringExtract[\"a bbb  cccc aa   d\", -1]", "\"d\"", 0);
    assert_eval_eq("StringExtract[\"a bbb  cccc aa   d\", 2 ;; 4]",
                   "{\"bbb\", \"cccc\", \"aa\"}", 0);
    assert_eval_eq("StringExtract[\"a bbb  cccc aa   d\", {1, 3}]",
                   "{\"a\", \"cccc\"}", 0);
    assert_eval_eq("StringExtract[\"a bbb  cccc aa   d\", {1, 4, 5}]",
                   "{\"a\", \"aa\", \"d\"}", 0);
}

/* {pos..} is equivalent to Part[StringSplit[..], {pos..}]. */
static void test_extract_part_equivalence(void) {
    assert_eval_eq("StringExtract[\"a bbb  cccc aa   d\", {1, 4, 5}]",
                   "Part[StringSplit[\"a bbb  cccc aa   d\"], {1, 4, 5}]", 0);
}

/* All extracts every whitespace-delimited token (== StringSplit). */
static void test_extract_all(void) {
    assert_eval_eq("StringExtract[\"A tree, an apple, four pears. And more: two sacks\", All]",
                   "{\"A\", \"tree,\", \"an\", \"apple,\", \"four\", \"pears.\", "
                   "\"And\", \"more:\", \"two\", \"sacks\"}", 0);
    assert_eval_eq("StringExtract[\"A tree, an apple, four pears. And more: two sacks\", All]",
                   "StringSplit[\"A tree, an apple, four pears. And more: two sacks\"]", 0);
}

/* sep -> pos with a literal separator and a pattern separator. */
static void test_extract_separator(void) {
    assert_eval_eq("StringExtract[\"a--bbb--ccc--dddd\", \"--\" -> 3]", "\"ccc\"", 0);
    /* sep -> All is exactly StringSplit[string, sep]. */
    assert_eval_eq("StringExtract[\"a--bbb---ccc--dddd\", \"--\" -> All]",
                   "StringSplit[\"a--bbb---ccc--dddd\", \"--\"]", 0);
    /* A string pattern as separator: runs of non-word characters. */
    assert_eval_eq("StringExtract[\"A tree, an apple, four pears. And more: two sacks\", "
                   "Except[WordCharacter] .. -> All]",
                   "{\"A\", \"tree\", \"an\", \"apple\", \"four\", \"pears\", "
                   "\"And\", \"more\", \"two\", \"sacks\"}", 0);
}

/* Multi-level extraction with whitespace/newline depth defaults. */
static void test_extract_multilevel(void) {
    /* table = "a 1\nb 2\nc 3 x" with real newline bytes. */
    assert_eval_eq("StringExtract[\"a 1\nb 2\nc 3 x\", All, 1]",
                   "{\"a\", \"b\", \"c\"}", 0);
    assert_eval_eq("StringExtract[\"a 1\nb 2\nc 3 x\", 3, 1]", "\"c\"", 0);
    assert_eval_eq("StringExtract[\"a 1\nb 2\nc 3 x\", 3, All]",
                   "{\"c\", \"3\", \"x\"}", 0);
    /* Three-level: "\n\n" -> 2, "\n" -> 2, Whitespace -> 3. */
    assert_eval_eq("StringExtract[\"a b c\nd e f\ng h i\n\nj k l\nm n o\np q r s"
                   "\n\nt u v\nw x y z\", 2, 2, 3]", "\"o\"", 0);
}

/* Out-of-range blocks yield Missing["PartAbsent", pos]. */
static void test_extract_missing(void) {
    assert_eval_eq("StringExtract[\"a b c\", 5]", "Missing[\"PartAbsent\", 5]", 0);
    assert_eval_eq("StringExtract[\"a b c\", -9]", "Missing[\"PartAbsent\", -9]", 0);
    /* Mapped across lines: absent third word becomes Missing per line. */
    assert_eval_eq("StringExtract[\"a 1\nb 2\nc 3 x\", All, 3]",
                   "{Missing[\"PartAbsent\", 3], Missing[\"PartAbsent\", 3], \"x\"}", 0);
}

/* A list of strings threads. */
static void test_extract_thread(void) {
    assert_eval_eq("StringExtract[{\"a b c\", \"x y z\"}, 2]", "{\"b\", \"y\"}", 0);
    assert_eval_eq("StringExtract[{\"a b c\", \"x y z\"}, {1, 3}]",
                   "{{\"a\", \"c\"}, {\"x\", \"z\"}}", 0);
}

static void test_extract_unevaluated(void) {
    /* A non-string, non-list subject leaves the call unevaluated. */
    assert_eval_eq("StringExtract[xyz, 2]", "StringExtract[xyz, 2]", 0);
    /* Zero / one argument: unevaluated (a StringExtract::argm message is written). */
    assert_eval_eq("StringExtract[]", "StringExtract[]", 0);
    assert_eval_eq("StringExtract[\"abc\"]", "StringExtract[\"abc\"]", 0);
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
    TEST(test_split_whitespace_default);
    TEST(test_split_literal_alternatives);
    TEST(test_split_charclass);
    TEST(test_split_stringexpression);
    TEST(test_split_rules);
    TEST(test_split_n_and_all);
    TEST(test_split_null_and_ignorecase);
    TEST(test_split_thread_and_regex);
    TEST(test_split_unevaluated);
    TEST(test_cases_wl_pattern);

    TEST(test_trim_whitespace_default);
    TEST(test_trim_pattern);
    TEST(test_trim_regex_front_only);
    TEST(test_trim_thread);
    TEST(test_trim_unevaluated);

    TEST(test_extract_positions);
    TEST(test_extract_part_equivalence);
    TEST(test_extract_all);
    TEST(test_extract_separator);
    TEST(test_extract_multilevel);
    TEST(test_extract_missing);
    TEST(test_extract_thread);
    TEST(test_extract_unevaluated);
#else
    printf("USE_REGEX not defined; skipping regex string-function tests\n");
#endif

    printf("All string-function tests passed!\n");
    return 0;
}
