/*
 * test_strings.c - Unit tests for StringLength and Characters builtins
 */

#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"

/* ===== StringLength tests ===== */

void test_stringlength_basic() {
    assert_eval_eq("StringLength[\"tiger\"]", "5", 0);
}

void test_stringlength_empty() {
    assert_eval_eq("StringLength[\"\"]", "0", 0);
}

void test_stringlength_single_char() {
    assert_eval_eq("StringLength[\"a\"]", "1", 0);
}

void test_stringlength_spaces() {
    assert_eval_eq("StringLength[\"hello world\"]", "11", 0);
}

void test_stringlength_special_chars() {
    assert_eval_eq("StringLength[\"a\\nb\"]", "3", 0);
}

void test_stringlength_digits() {
    assert_eval_eq("StringLength[\"12345\"]", "5", 0);
}

void test_stringlength_mixed() {
    assert_eval_eq("StringLength[\"abc 123!@#\"]", "10", 0);
}

void test_stringlength_symbolic() {
    /* Non-string argument: should remain unevaluated */
    assert_eval_eq("StringLength[x]", "StringLength[x]", 0);
}

void test_stringlength_integer_arg() {
    /* Integer argument: should remain unevaluated */
    assert_eval_eq("StringLength[42]", "StringLength[42]", 0);
}

void test_stringlength_no_args() {
    /* No arguments: should remain unevaluated */
    assert_eval_eq("StringLength[]", "StringLength[]", 0);
}

void test_stringlength_two_args() {
    /* Two arguments: should remain unevaluated */
    assert_eval_eq("StringLength[\"a\", \"b\"]", "StringLength[\"a\", \"b\"]", 0);
}

void test_stringlength_listable() {
    /* Listable: threads over lists */
    assert_eval_eq("StringLength[{\"ABC\", \"DE\", \"F\"}]", "{3, 2, 1}", 0);
}

void test_stringlength_listable_with_empty() {
    assert_eval_eq("StringLength[{\"hello\", \"\", \"world\"}]", "{5, 0, 5}", 0);
}

void test_stringlength_long_string() {
    assert_eval_eq("StringLength[\"abcdefghijklmnopqrstuvwxyz\"]", "26", 0);
}

void test_stringlength_tab() {
    assert_eval_eq("StringLength[\"a\\tb\"]", "3", 0);
}

/* ===== Characters tests ===== */

void test_characters_basic() {
    assert_eval_eq("Characters[\"ABC\"]", "{\"A\", \"B\", \"C\"}", 0);
}

void test_characters_empty() {
    assert_eval_eq("Characters[\"\"]", "{}", 0);
}

void test_characters_single_char() {
    assert_eval_eq("Characters[\"x\"]", "{\"x\"}", 0);
}

void test_characters_with_spaces() {
    assert_eval_eq("Characters[\"A string.\"]", "{\"A\", \" \", \"s\", \"t\", \"r\", \"i\", \"n\", \"g\", \".\"}", 0);
}

void test_characters_digits() {
    assert_eval_eq("Characters[\"123\"]", "{\"1\", \"2\", \"3\"}", 0);
}

void test_characters_symbolic() {
    /* Non-string argument: should remain unevaluated */
    assert_eval_eq("Characters[x]", "Characters[x]", 0);
}

void test_characters_integer_arg() {
    /* Integer argument: should remain unevaluated */
    assert_eval_eq("Characters[42]", "Characters[42]", 0);
}

void test_characters_no_args() {
    assert_eval_eq("Characters[]", "Characters[]", 0);
}

void test_characters_two_args() {
    assert_eval_eq("Characters[\"a\", \"b\"]", "Characters[\"a\", \"b\"]", 0);
}

void test_characters_listable() {
    /* Listable: threads over lists */
    assert_eval_eq("Characters[{\"ABC\", \"DEF\", \"XYZ\"}]",
                   "{{\"A\", \"B\", \"C\"}, {\"D\", \"E\", \"F\"}, {\"X\", \"Y\", \"Z\"}}", 0);
}

void test_characters_listable_mixed_lengths() {
    assert_eval_eq("Characters[{\"ab\", \"c\"}]",
                   "{{\"a\", \"b\"}, {\"c\"}}", 0);
}

void test_characters_special_chars() {
    assert_eval_eq("Characters[\"!@#\"]", "{\"!\", \"@\", \"#\"}", 0);
}

/* ===== StringJoin tests ===== */

void test_stringjoin_basic() {
    assert_eval_eq("StringJoin[\"abcd\", \"ABCD\", \"xyz\"]", "\"abcdABCDxyz\"", 0);
}

void test_stringjoin_two_args() {
    assert_eval_eq("StringJoin[\"hello\", \" world\"]", "\"hello world\"", 0);
}

void test_stringjoin_single_arg() {
    assert_eval_eq("StringJoin[\"abc\"]", "\"abc\"", 0);
}

void test_stringjoin_empty_args() {
    assert_eval_eq("StringJoin[]", "\"\"", 0);
}

void test_stringjoin_empty_strings() {
    assert_eval_eq("StringJoin[\"\", \"\", \"\"]", "\"\"", 0);
}

void test_stringjoin_mixed_empty() {
    assert_eval_eq("StringJoin[\"a\", \"\", \"b\"]", "\"ab\"", 0);
}

void test_stringjoin_flatten_list() {
    /* All lists are recursively flattened */
    assert_eval_eq("StringJoin[{{\"AB\", \"CD\"}, \"XY\"}]", "\"ABCDXY\"", 0);
}

void test_stringjoin_nested_lists() {
    assert_eval_eq("StringJoin[{{{\"a\"}, \"b\"}, \"c\"}]", "\"abc\"", 0);
}

void test_stringjoin_list_with_multiple_args() {
    assert_eval_eq("StringJoin[{\"A\", \"B\"}, \"C\"]", "\"ABC\"", 0);
}

void test_stringjoin_symbolic_unevaluated() {
    /* Non-string argument: remain unevaluated */
    assert_eval_eq("StringJoin[\"a\", x]", "StringJoin[\"a\", x]", 0);
}

void test_stringjoin_integer_unevaluated() {
    assert_eval_eq("StringJoin[\"a\", 42]", "StringJoin[\"a\", 42]", 0);
}

void test_stringjoin_newline() {
    /* Parser treats \n as literal 'n' (no escape interpretation) */
    assert_eval_eq("StringJoin[\"one\", \"\\n\", \"two\"]", "\"onentwo\"", 0);
}

void test_stringjoin_many_args() {
    assert_eval_eq("StringJoin[\"a\", \"b\", \"c\", \"d\", \"e\"]", "\"abcde\"", 0);
}

/* ===== <> infix parser tests ===== */

void test_infix_stringjoin_basic() {
    assert_eval_eq("\"abcd\" <> \"ABCD\" <> \"xyz\"", "\"abcdABCDxyz\"", 0);
}

void test_infix_stringjoin_two() {
    assert_eval_eq("\"hello\" <> \" world\"", "\"hello world\"", 0);
}

void test_infix_stringjoin_empty() {
    assert_eval_eq("\"\" <> \"\"", "\"\"", 0);
}

void test_infix_stringjoin_chain_four() {
    assert_eval_eq("\"a\" <> \"b\" <> \"c\" <> \"d\"", "\"abcd\"", 0);
}

void test_infix_stringjoin_symbolic() {
    assert_eval_eq("\"a\" <> x", "StringJoin[\"a\", x]", 0);
}

void test_infix_stringjoin_newline() {
    /* Parser treats \n as literal 'n' (no escape interpretation) */
    assert_eval_eq("\"one\" <> \"\\n\" <> \"two\"", "\"onentwo\"", 0);
}

void test_infix_stringjoin_fullform() {
    /* FullForm of string result is just the string */
    assert_eval_eq("\"a\" <> \"b\"", "\"ab\"", 1);
}

void test_infix_stringjoin_symbolic_fullform() {
    /* Verify parser produces StringJoin head for unevaluated cases */
    assert_eval_eq("\"a\" <> x", "StringJoin[\"a\", x]", 1);
}

void test_infix_stringjoin_three_fullform() {
    /* Flat attribute should flatten chained <> into one StringJoin */
    assert_eval_eq("\"a\" <> x <> \"c\"", "StringJoin[\"a\", x, \"c\"]", 1);
}

void test_infix_stringjoin_with_parens() {
    assert_eval_eq("(\"a\" <> \"b\") <> \"c\"", "\"abc\"", 0);
}

void test_infix_stringjoin_spaces() {
    assert_eval_eq("\"hello\"  <>  \"world\"", "\"helloworld\"", 0);
}

/* ===== StringPart tests ===== */

/* Single index */
void test_stringpart_single_index() {
    assert_eval_eq("StringPart[\"abcdefghijklm\", 6]", "\"f\"", 0);
}

void test_stringpart_first_char() {
    assert_eval_eq("StringPart[\"abcdef\", 1]", "\"a\"", 0);
}

void test_stringpart_last_char() {
    assert_eval_eq("StringPart[\"abcdef\", 6]", "\"f\"", 0);
}

/* Negative index */
void test_stringpart_negative_index() {
    assert_eval_eq("StringPart[\"abcdefghijklm\", -4]", "\"j\"", 0);
}

void test_stringpart_negative_one() {
    assert_eval_eq("StringPart[\"abcdef\", -1]", "\"f\"", 0);
}

void test_stringpart_negative_full() {
    assert_eval_eq("StringPart[\"abcdef\", -6]", "\"a\"", 0);
}

/* List of indices */
void test_stringpart_list_indices() {
    assert_eval_eq("StringPart[\"abcdefghijklm\", {1, 3, 5}]", "{\"a\", \"c\", \"e\"}", 0);
}

void test_stringpart_list_single() {
    assert_eval_eq("StringPart[\"abcdef\", {3}]", "{\"c\"}", 0);
}

void test_stringpart_list_with_negative() {
    assert_eval_eq("StringPart[\"abcdef\", {1, -1}]", "{\"a\", \"f\"}", 0);
}

/* Span: m;;n */
void test_stringpart_span_basic() {
    assert_eval_eq("StringPart[\"abcdefghijklm\", 1;;6]",
                   "{\"a\", \"b\", \"c\", \"d\", \"e\", \"f\"}", 0);
}

void test_stringpart_span_middle() {
    assert_eval_eq("StringPart[\"abcdef\", 2;;4]",
                   "{\"b\", \"c\", \"d\"}", 0);
}

/* Span: m;;n;;s */
void test_stringpart_span_step() {
    assert_eval_eq("StringPart[\"abcdefghijklm\", 1;;-1;;2]",
                   "{\"a\", \"c\", \"e\", \"g\", \"i\", \"k\", \"m\"}", 0);
}

void test_stringpart_span_negative_step() {
    assert_eval_eq("StringPart[\"abcdefghijklm\", -1;;1;;-2]",
                   "{\"m\", \"k\", \"i\", \"g\", \"e\", \"c\", \"a\"}", 0);
}

void test_stringpart_span_step_3() {
    assert_eval_eq("StringPart[\"abcdefghijklm\", 1;;-1;;3]",
                   "{\"a\", \"d\", \"g\", \"j\", \"m\"}", 0);
}

/* List of strings */
void test_stringpart_list_strings_single_index() {
    assert_eval_eq("StringPart[{\"abcd\", \"efgh\", \"ijklm\"}, 1]",
                   "{\"a\", \"e\", \"i\"}", 0);
}

void test_stringpart_list_strings_list_indices() {
    assert_eval_eq("StringPart[{\"abcd\", \"efgh\", \"ijklm\"}, {1, -1}]",
                   "{{\"a\", \"d\"}, {\"e\", \"h\"}, {\"i\", \"m\"}}", 0);
}

/* Edge cases */
void test_stringpart_unevaluated_no_args() {
    assert_eval_eq("StringPart[]", "StringPart[]", 0);
}

void test_stringpart_unevaluated_one_arg() {
    assert_eval_eq("StringPart[\"abc\"]", "StringPart[\"abc\"]", 0);
}

void test_stringpart_unevaluated_three_args() {
    assert_eval_eq("StringPart[\"abc\", 1, 2]", "StringPart[\"abc\", 1, 2]", 0);
}

void test_stringpart_unevaluated_symbolic() {
    assert_eval_eq("StringPart[x, 1]", "StringPart[x, 1]", 0);
}

void test_stringpart_unevaluated_symbolic_index() {
    assert_eval_eq("StringPart[\"abc\", x]", "StringPart[\"abc\", x]", 0);
}

void test_stringpart_out_of_bounds() {
    /* Out of bounds should remain unevaluated */
    assert_eval_eq("StringPart[\"abc\", 4]", "StringPart[\"abc\", 4]", 0);
}

void test_stringpart_out_of_bounds_negative() {
    assert_eval_eq("StringPart[\"abc\", -4]", "StringPart[\"abc\", -4]", 0);
}

void test_stringpart_zero_index() {
    /* 0 is not a valid index */
    assert_eval_eq("StringPart[\"abc\", 0]", "StringPart[\"abc\", 0]", 0);
}

void test_stringpart_span_all() {
    assert_eval_eq("StringPart[\"abc\", 1;;3]",
                   "{\"a\", \"b\", \"c\"}", 0);
}

void test_stringpart_span_reverse() {
    assert_eval_eq("StringPart[\"abcde\", 5;;1;;-1]",
                   "{\"e\", \"d\", \"c\", \"b\", \"a\"}", 0);
}

void test_stringpart_span_single_element() {
    assert_eval_eq("StringPart[\"abcdef\", 3;;3]", "{\"c\"}", 0);
}

void test_stringpart_span_negative_endpoints() {
    assert_eval_eq("StringPart[\"abcdef\", -3;;-1]",
                   "{\"d\", \"e\", \"f\"}", 0);
}

/* ===== StringTake tests ===== */

/* StringTake["string", n] - first n characters */
void test_stringtake_first_n() {
    assert_eval_eq("StringTake[\"abcdefghijklm\", 6]", "\"abcdef\"", 0);
}

void test_stringtake_first_1() {
    assert_eval_eq("StringTake[\"abcdef\", 1]", "\"a\"", 0);
}

void test_stringtake_first_all() {
    assert_eval_eq("StringTake[\"abc\", 3]", "\"abc\"", 0);
}

void test_stringtake_zero() {
    assert_eval_eq("StringTake[\"abc\", 0]", "\"\"", 0);
}

/* StringTake["string", -n] - last n characters */
void test_stringtake_last_n() {
    assert_eval_eq("StringTake[\"abcdefghijklm\", -4]", "\"jklm\"", 0);
}

void test_stringtake_last_1() {
    assert_eval_eq("StringTake[\"abcdef\", -1]", "\"f\"", 0);
}

void test_stringtake_last_all() {
    assert_eval_eq("StringTake[\"abc\", -3]", "\"abc\"", 0);
}

/* StringTake["string", {n}] - nth character */
void test_stringtake_single_char() {
    assert_eval_eq("StringTake[\"abcdefghijklm\", {6}]", "\"f\"", 0);
}

void test_stringtake_single_first() {
    assert_eval_eq("StringTake[\"abcdef\", {1}]", "\"a\"", 0);
}

void test_stringtake_single_last() {
    assert_eval_eq("StringTake[\"abcdef\", {-1}]", "\"f\"", 0);
}

void test_stringtake_single_negative() {
    assert_eval_eq("StringTake[\"abcdef\", {-3}]", "\"d\"", 0);
}

/* StringTake["string", {m, n}] - characters m through n */
void test_stringtake_range() {
    assert_eval_eq("StringTake[\"abcdefghijklm\", {5, 10}]", "\"efghij\"", 0);
}

void test_stringtake_range_from_start() {
    assert_eval_eq("StringTake[\"abcdef\", {1, 3}]", "\"abc\"", 0);
}

void test_stringtake_range_to_end() {
    assert_eval_eq("StringTake[\"abcdef\", {4, 6}]", "\"def\"", 0);
}

void test_stringtake_range_single() {
    assert_eval_eq("StringTake[\"abcdef\", {3, 3}]", "\"c\"", 0);
}

void test_stringtake_range_negative() {
    assert_eval_eq("StringTake[\"abcdef\", {-3, -1}]", "\"def\"", 0);
}

/* StringTake["string", {m, n, s}] - range with step */
void test_stringtake_step() {
    assert_eval_eq("StringTake[\"abcdefghijklm\", {1, -1, 2}]", "\"acegikm\"", 0);
}

void test_stringtake_step_3() {
    assert_eval_eq("StringTake[\"abcdefghijklm\", {1, -1, 3}]", "\"adgjm\"", 0);
}

void test_stringtake_step_negative() {
    assert_eval_eq("StringTake[\"abcde\", {5, 1, -1}]", "\"edcba\"", 0);
}

void test_stringtake_step_neg2() {
    assert_eval_eq("StringTake[\"abcdefghijklm\", {-1, 1, -2}]", "\"mkigeca\"", 0);
}

/* StringTake["string", UpTo[n]] */
void test_stringtake_upto_enough() {
    assert_eval_eq("StringTake[\"abcdef\", UpTo[4]]", "\"abcd\"", 0);
}

void test_stringtake_upto_more_than_available() {
    assert_eval_eq("StringTake[\"abc\", UpTo[4]]", "\"abc\"", 0);
}

void test_stringtake_upto_exact() {
    assert_eval_eq("StringTake[\"abc\", UpTo[3]]", "\"abc\"", 0);
}

void test_stringtake_upto_one() {
    assert_eval_eq("StringTake[\"abc\", UpTo[1]]", "\"a\"", 0);
}

void test_stringtake_upto_large() {
    assert_eval_eq("StringTake[\"abc\", UpTo[100]]", "\"abc\"", 0);
}

/* StringTake[{s1, s2, ...}, spec] - list of strings */
void test_stringtake_list_strings() {
    assert_eval_eq("StringTake[{\"abcdef\", \"stuv\", \"xyzw\"}, -2]",
                   "{\"ef\", \"uv\", \"zw\"}", 0);
}

void test_stringtake_list_strings_first() {
    assert_eval_eq("StringTake[{\"abcdef\", \"xyz\"}, 2]",
                   "{\"ab\", \"xy\"}", 0);
}

void test_stringtake_list_strings_range() {
    assert_eval_eq("StringTake[{\"abcdef\", \"ghijkl\"}, {2, 4}]",
                   "{\"bcd\", \"hij\"}", 0);
}

/* Edge cases - unevaluated */
void test_stringtake_unevaluated_no_args() {
    assert_eval_eq("StringTake[]", "StringTake[]", 0);
}

void test_stringtake_unevaluated_one_arg() {
    assert_eval_eq("StringTake[\"abc\"]", "StringTake[\"abc\"]", 0);
}

void test_stringtake_unevaluated_three_args() {
    assert_eval_eq("StringTake[\"abc\", 1, 2]", "StringTake[\"abc\", 1, 2]", 0);
}

void test_stringtake_unevaluated_symbolic() {
    assert_eval_eq("StringTake[x, 1]", "StringTake[x, 1]", 0);
}

void test_stringtake_unevaluated_symbolic_spec() {
    assert_eval_eq("StringTake[\"abc\", x]", "StringTake[\"abc\", x]", 0);
}

void test_stringtake_out_of_bounds() {
    assert_eval_eq("StringTake[\"abc\", 4]", "StringTake[\"abc\", 4]", 0);
}

void test_stringtake_out_of_bounds_negative() {
    assert_eval_eq("StringTake[\"abc\", -4]", "StringTake[\"abc\", -4]", 0);
}

void test_stringtake_empty_string_zero() {
    assert_eval_eq("StringTake[\"\", 0]", "\"\"", 0);
}

void test_stringtake_empty_upto() {
    assert_eval_eq("StringTake[\"\", UpTo[5]]", "\"\"", 0);
}

/* ===== Combined / integration tests ===== */

void test_stringjoin_with_characters() {
    /* Join the characters back together */
    assert_eval_eq("StringJoin[Characters[\"hello\"]]", "\"hello\"", 0);
}

void test_stringlength_of_stringjoin() {
    assert_eval_eq("StringLength[\"abc\" <> \"def\"]", "6", 0);
}

void test_stringlength_of_characters() {
    /* StringLength applied to the result of Characters should give list of 1s */
    assert_eval_eq("StringLength /@ Characters[\"hello\"]", "{1, 1, 1, 1, 1}", 0);
}

void test_length_of_characters() {
    /* Length of Characters result should equal StringLength */
    assert_eval_eq("Length[Characters[\"tiger\"]]", "5", 0);
}

/* ===== StringReverse tests ===== */

void test_stringreverse_basic() {
    assert_eval_eq("StringReverse[\"abcdef\"]", "\"fedcba\"", 0);
}

void test_stringreverse_single_char() {
    assert_eval_eq("StringReverse[\"a\"]", "\"a\"", 0);
}

void test_stringreverse_empty() {
    assert_eval_eq("StringReverse[\"\"]", "\"\"", 0);
}

void test_stringreverse_palindrome() {
    assert_eval_eq("StringReverse[\"racecar\"]", "\"racecar\"", 0);
}

void test_stringreverse_with_spaces() {
    assert_eval_eq("StringReverse[\"a b c\"]", "\"c b a\"", 0);
}

void test_stringreverse_digits() {
    assert_eval_eq("StringReverse[\"12345\"]", "\"54321\"", 0);
}

void test_stringreverse_special_chars() {
    /* Parser treats \n as the literal char 'n' (no escape interpretation),
     * so "a\nb" is the 3-char string "anb", which reverses to "bna". */
    assert_eval_eq("StringReverse[\"a\\nb\"]", "\"bna\"", 0);
}

void test_stringreverse_long() {
    assert_eval_eq("StringReverse[\"abcdefghijklmnopqrstuvwxyz\"]",
                   "\"zyxwvutsrqponmlkjihgfedcba\"", 0);
}

/* Involution: reversing twice gives the original. */
void test_stringreverse_involution() {
    assert_eval_eq("StringReverse[StringReverse[\"coelenterate\"]]",
                   "\"coelenterate\"", 0);
}

/* Listable: threads automatically over a list of strings. */
void test_stringreverse_listable() {
    assert_eval_eq("StringReverse[{\"cat\", \"dog\", \"fish\", \"coelenterate\"}]",
                   "{\"tac\", \"god\", \"hsif\", \"etaretneleoc\"}", 0);
}

void test_stringreverse_listable_empty_elements() {
    assert_eval_eq("StringReverse[{\"\", \"ab\", \"\"}]",
                   "{\"\", \"ba\", \"\"}", 0);
}

/* Composes with StringJoin. */
void test_stringreverse_of_stringjoin() {
    assert_eval_eq("StringReverse[\"abc\" <> \"def\"]", "\"fedcba\"", 0);
}

/* Non-string argument: left unevaluated (symbolic flow-through). */
void test_stringreverse_symbolic() {
    assert_eval_eq("StringReverse[x]", "StringReverse[x]", 0);
}

void test_stringreverse_integer_arg() {
    assert_eval_eq("StringReverse[12345]", "StringReverse[12345]", 0);
}

/* Wrong argument counts: StringReverse::argx, call left unevaluated. */
void test_stringreverse_no_args() {
    assert_eval_eq("StringReverse[]", "StringReverse[]", 0);
}

void test_stringreverse_two_args() {
    assert_eval_eq("StringReverse[\"a\", \"b\"]", "StringReverse[\"a\", \"b\"]", 0);
}

/* ===== StringDrop tests ===== */

/* StringDrop["string", n] - drop the first n characters */
void test_stringdrop_first_n() {
    assert_eval_eq("StringDrop[\"abcdefghijklm\", 4]", "\"efghijklm\"", 0);
}

void test_stringdrop_first_1() {
    assert_eval_eq("StringDrop[\"abcdef\", 1]", "\"bcdef\"", 0);
}

void test_stringdrop_first_all() {
    /* Dropping every character yields the empty string. */
    assert_eval_eq("StringDrop[\"abc\", 3]", "\"\"", 0);
}

void test_stringdrop_zero() {
    /* Dropping zero characters returns the whole string. */
    assert_eval_eq("StringDrop[\"abc\", 0]", "\"abc\"", 0);
}

/* StringDrop["string", -n] - drop the last n characters */
void test_stringdrop_last_n() {
    assert_eval_eq("StringDrop[\"abcdefghijklm\", -4]", "\"abcdefghi\"", 0);
}

void test_stringdrop_last_1() {
    assert_eval_eq("StringDrop[\"abcdef\", -1]", "\"abcde\"", 0);
}

void test_stringdrop_last_all() {
    assert_eval_eq("StringDrop[\"abc\", -3]", "\"\"", 0);
}

/* StringDrop["string", {n}] - drop the nth character */
void test_stringdrop_single_char() {
    assert_eval_eq("StringDrop[\"abcdefghijklm\", {3}]", "\"abdefghijklm\"", 0);
}

void test_stringdrop_single_first() {
    assert_eval_eq("StringDrop[\"abcdef\", {1}]", "\"bcdef\"", 0);
}

void test_stringdrop_single_last() {
    assert_eval_eq("StringDrop[\"abcdef\", {-1}]", "\"abcde\"", 0);
}

void test_stringdrop_single_negative() {
    assert_eval_eq("StringDrop[\"abcdef\", {-3}]", "\"abcef\"", 0);
}

/* StringDrop["string", {m, n}] - drop characters m through n */
void test_stringdrop_range() {
    assert_eval_eq("StringDrop[\"abcdefghijklm\", {5, 10}]", "\"abcdklm\"", 0);
}

void test_stringdrop_range_from_start() {
    assert_eval_eq("StringDrop[\"abcdef\", {1, 3}]", "\"def\"", 0);
}

void test_stringdrop_range_to_end() {
    assert_eval_eq("StringDrop[\"abcdef\", {4, 6}]", "\"abc\"", 0);
}

void test_stringdrop_range_single() {
    assert_eval_eq("StringDrop[\"abcdef\", {3, 3}]", "\"abdef\"", 0);
}

void test_stringdrop_range_negative() {
    assert_eval_eq("StringDrop[\"abcdef\", {-3, -1}]", "\"abc\"", 0);
}

void test_stringdrop_range_mixed_signs() {
    /* {5, -4} on a 13-char string drops positions 5..10. */
    assert_eval_eq("StringDrop[\"abcdefghijklm\", {5, -4}]", "\"abcdklm\"", 0);
}

void test_stringdrop_range_all() {
    assert_eval_eq("StringDrop[\"abcdef\", {1, 6}]", "\"\"", 0);
}

void test_stringdrop_range_empty() {
    /* A decreasing range is empty and drops nothing. */
    assert_eval_eq("StringDrop[\"abc\", {5, 4}]", "\"abc\"", 0);
}

void test_stringdrop_range_empty_in_bounds() {
    assert_eval_eq("StringDrop[\"abcdef\", {4, 2}]", "\"abcdef\"", 0);
}

/* StringDrop["string", {m, n, s}] - drop range with step */
void test_stringdrop_step() {
    assert_eval_eq("StringDrop[\"abcdefghijklm\", {1, -1, 2}]", "\"bdfhjl\"", 0);
}

void test_stringdrop_step_3() {
    assert_eval_eq("StringDrop[\"abcdefghijklm\", {1, -1, 3}]", "\"bcefhikl\"", 0);
}

void test_stringdrop_step_negative() {
    /* Dropping positions 5,4,3,2,1 leaves the tail. */
    assert_eval_eq("StringDrop[\"abcde\", {5, 1, -1}]", "\"\"", 0);
}

void test_stringdrop_step_range() {
    assert_eval_eq("StringDrop[\"abcdefghij\", {2, 8, 2}]", "\"acegij\"", 0);
}

/* StringDrop["string", UpTo[n]] */
void test_stringdrop_upto_enough() {
    assert_eval_eq("StringDrop[\"abcdef\", UpTo[4]]", "\"ef\"", 0);
}

void test_stringdrop_upto_more_than_available() {
    assert_eval_eq("StringDrop[\"abc\", UpTo[4]]", "\"\"", 0);
}

void test_stringdrop_upto_exact() {
    assert_eval_eq("StringDrop[\"abc\", UpTo[3]]", "\"\"", 0);
}

void test_stringdrop_upto_one() {
    assert_eval_eq("StringDrop[\"abc\", UpTo[1]]", "\"bc\"", 0);
}

void test_stringdrop_upto_large() {
    assert_eval_eq("StringDrop[\"abc\", UpTo[100]]", "\"\"", 0);
}

/* StringDrop[{s1, s2, ...}, spec] - list of strings */
void test_stringdrop_list_strings() {
    assert_eval_eq("StringDrop[{\"abcdef\", \"xyzw\", \"stuv\"}, -2]",
                   "{\"abcd\", \"xy\", \"st\"}", 0);
}

void test_stringdrop_list_strings_first() {
    assert_eval_eq("StringDrop[{\"abcdef\", \"xyz\"}, 2]",
                   "{\"cdef\", \"z\"}", 0);
}

void test_stringdrop_list_strings_range() {
    assert_eval_eq("StringDrop[{\"abcdef\", \"ghijkl\"}, {2, 4}]",
                   "{\"aef\", \"gkl\"}", 0);
}

/* Newline (parsed as literal 'n') counts as a single character. */
void test_stringdrop_newline() {
    assert_eval_eq("StringDrop[\"abc\\ndef\", 4]", "\"def\"", 0);
}

/* Consistency with Drop over Characters (both sides are evaluated). */
void test_stringdrop_matches_drop_characters() {
    assert_eval_eq("StringDrop[\"abcdefghijklm\", {5, -4}] === "
                   "StringJoin[Drop[Characters[\"abcdefghijklm\"], {5, -4}]]",
                   "True", 0);
}

/* Edge cases - out of bounds / unevaluated */
void test_stringdrop_out_of_bounds() {
    assert_eval_eq("StringDrop[\"abc\", 5]", "StringDrop[\"abc\", 5]", 0);
}

void test_stringdrop_out_of_bounds_negative() {
    assert_eval_eq("StringDrop[\"abc\", -5]", "StringDrop[\"abc\", -5]", 0);
}

void test_stringdrop_single_out_of_bounds() {
    assert_eval_eq("StringDrop[\"abc\", {5}]", "StringDrop[\"abc\", {5}]", 0);
}

void test_stringdrop_range_out_of_bounds() {
    assert_eval_eq("StringDrop[\"abc\", {2, 9}]", "StringDrop[\"abc\", {2, 9}]", 0);
}

void test_stringdrop_step_zero() {
    assert_eval_eq("StringDrop[\"abc\", {1, 3, 0}]", "StringDrop[\"abc\", {1, 3, 0}]", 0);
}

void test_stringdrop_unevaluated_symbolic() {
    assert_eval_eq("StringDrop[x, 1]", "StringDrop[x, 1]", 0);
}

void test_stringdrop_unevaluated_symbolic_spec() {
    assert_eval_eq("StringDrop[\"abc\", x]", "StringDrop[\"abc\", x]", 0);
}

void test_stringdrop_no_args() {
    /* StringDrop::argrx, left unevaluated. */
    assert_eval_eq("StringDrop[]", "StringDrop[]", 0);
}

void test_stringdrop_one_arg() {
    assert_eval_eq("StringDrop[\"abc\"]", "StringDrop[\"abc\"]", 0);
}

void test_stringdrop_three_args() {
    assert_eval_eq("StringDrop[\"abc\", 1, 2]", "StringDrop[\"abc\", 1, 2]", 0);
}

void test_stringdrop_empty_string_zero() {
    assert_eval_eq("StringDrop[\"\", 0]", "\"\"", 0);
}

void test_stringdrop_empty_string_upto() {
    assert_eval_eq("StringDrop[\"\", UpTo[5]]", "\"\"", 0);
}

int main() {
    symtab_init();
    core_init();

    /* StringLength tests */
    TEST(test_stringlength_basic);
    TEST(test_stringlength_empty);
    TEST(test_stringlength_single_char);
    TEST(test_stringlength_spaces);
    TEST(test_stringlength_special_chars);
    TEST(test_stringlength_digits);
    TEST(test_stringlength_mixed);
    TEST(test_stringlength_symbolic);
    TEST(test_stringlength_integer_arg);
    TEST(test_stringlength_no_args);
    TEST(test_stringlength_two_args);
    TEST(test_stringlength_listable);
    TEST(test_stringlength_listable_with_empty);
    TEST(test_stringlength_long_string);
    TEST(test_stringlength_tab);

    /* Characters tests */
    TEST(test_characters_basic);
    TEST(test_characters_empty);
    TEST(test_characters_single_char);
    TEST(test_characters_with_spaces);
    TEST(test_characters_digits);
    TEST(test_characters_symbolic);
    TEST(test_characters_integer_arg);
    TEST(test_characters_no_args);
    TEST(test_characters_two_args);
    TEST(test_characters_listable);
    TEST(test_characters_listable_mixed_lengths);
    TEST(test_characters_special_chars);

    /* StringJoin tests */
    TEST(test_stringjoin_basic);
    TEST(test_stringjoin_two_args);
    TEST(test_stringjoin_single_arg);
    TEST(test_stringjoin_empty_args);
    TEST(test_stringjoin_empty_strings);
    TEST(test_stringjoin_mixed_empty);
    TEST(test_stringjoin_flatten_list);
    TEST(test_stringjoin_nested_lists);
    TEST(test_stringjoin_list_with_multiple_args);
    TEST(test_stringjoin_symbolic_unevaluated);
    TEST(test_stringjoin_integer_unevaluated);
    TEST(test_stringjoin_newline);
    TEST(test_stringjoin_many_args);

    /* <> infix parser tests */
    TEST(test_infix_stringjoin_basic);
    TEST(test_infix_stringjoin_two);
    TEST(test_infix_stringjoin_empty);
    TEST(test_infix_stringjoin_chain_four);
    TEST(test_infix_stringjoin_symbolic);
    TEST(test_infix_stringjoin_newline);
    TEST(test_infix_stringjoin_fullform);
    TEST(test_infix_stringjoin_symbolic_fullform);
    TEST(test_infix_stringjoin_three_fullform);
    TEST(test_infix_stringjoin_with_parens);
    TEST(test_infix_stringjoin_spaces);

    /* StringPart tests */
    TEST(test_stringpart_single_index);
    TEST(test_stringpart_first_char);
    TEST(test_stringpart_last_char);
    TEST(test_stringpart_negative_index);
    TEST(test_stringpart_negative_one);
    TEST(test_stringpart_negative_full);
    TEST(test_stringpart_list_indices);
    TEST(test_stringpart_list_single);
    TEST(test_stringpart_list_with_negative);
    TEST(test_stringpart_span_basic);
    TEST(test_stringpart_span_middle);
    TEST(test_stringpart_span_step);
    TEST(test_stringpart_span_negative_step);
    TEST(test_stringpart_span_step_3);
    TEST(test_stringpart_list_strings_single_index);
    TEST(test_stringpart_list_strings_list_indices);
    TEST(test_stringpart_unevaluated_no_args);
    TEST(test_stringpart_unevaluated_one_arg);
    TEST(test_stringpart_unevaluated_three_args);
    TEST(test_stringpart_unevaluated_symbolic);
    TEST(test_stringpart_unevaluated_symbolic_index);
    TEST(test_stringpart_out_of_bounds);
    TEST(test_stringpart_out_of_bounds_negative);
    TEST(test_stringpart_zero_index);
    TEST(test_stringpart_span_all);
    TEST(test_stringpart_span_reverse);
    TEST(test_stringpart_span_single_element);
    TEST(test_stringpart_span_negative_endpoints);

    /* StringTake tests */
    TEST(test_stringtake_first_n);
    TEST(test_stringtake_first_1);
    TEST(test_stringtake_first_all);
    TEST(test_stringtake_zero);
    TEST(test_stringtake_last_n);
    TEST(test_stringtake_last_1);
    TEST(test_stringtake_last_all);
    TEST(test_stringtake_single_char);
    TEST(test_stringtake_single_first);
    TEST(test_stringtake_single_last);
    TEST(test_stringtake_single_negative);
    TEST(test_stringtake_range);
    TEST(test_stringtake_range_from_start);
    TEST(test_stringtake_range_to_end);
    TEST(test_stringtake_range_single);
    TEST(test_stringtake_range_negative);
    TEST(test_stringtake_step);
    TEST(test_stringtake_step_3);
    TEST(test_stringtake_step_negative);
    TEST(test_stringtake_step_neg2);
    TEST(test_stringtake_upto_enough);
    TEST(test_stringtake_upto_more_than_available);
    TEST(test_stringtake_upto_exact);
    TEST(test_stringtake_upto_one);
    TEST(test_stringtake_upto_large);
    TEST(test_stringtake_list_strings);
    TEST(test_stringtake_list_strings_first);
    TEST(test_stringtake_list_strings_range);
    TEST(test_stringtake_unevaluated_no_args);
    TEST(test_stringtake_unevaluated_one_arg);
    TEST(test_stringtake_unevaluated_three_args);
    TEST(test_stringtake_unevaluated_symbolic);
    TEST(test_stringtake_unevaluated_symbolic_spec);
    TEST(test_stringtake_out_of_bounds);
    TEST(test_stringtake_out_of_bounds_negative);
    TEST(test_stringtake_empty_string_zero);
    TEST(test_stringtake_empty_upto);

    /* Combined tests */
    TEST(test_stringjoin_with_characters);
    TEST(test_stringlength_of_stringjoin);
    TEST(test_stringlength_of_characters);
    TEST(test_length_of_characters);

    /* StringDrop tests */
    TEST(test_stringdrop_first_n);
    TEST(test_stringdrop_first_1);
    TEST(test_stringdrop_first_all);
    TEST(test_stringdrop_zero);
    TEST(test_stringdrop_last_n);
    TEST(test_stringdrop_last_1);
    TEST(test_stringdrop_last_all);
    TEST(test_stringdrop_single_char);
    TEST(test_stringdrop_single_first);
    TEST(test_stringdrop_single_last);
    TEST(test_stringdrop_single_negative);
    TEST(test_stringdrop_range);
    TEST(test_stringdrop_range_from_start);
    TEST(test_stringdrop_range_to_end);
    TEST(test_stringdrop_range_single);
    TEST(test_stringdrop_range_negative);
    TEST(test_stringdrop_range_mixed_signs);
    TEST(test_stringdrop_range_all);
    TEST(test_stringdrop_range_empty);
    TEST(test_stringdrop_range_empty_in_bounds);
    TEST(test_stringdrop_step);
    TEST(test_stringdrop_step_3);
    TEST(test_stringdrop_step_negative);
    TEST(test_stringdrop_step_range);
    TEST(test_stringdrop_upto_enough);
    TEST(test_stringdrop_upto_more_than_available);
    TEST(test_stringdrop_upto_exact);
    TEST(test_stringdrop_upto_one);
    TEST(test_stringdrop_upto_large);
    TEST(test_stringdrop_list_strings);
    TEST(test_stringdrop_list_strings_first);
    TEST(test_stringdrop_list_strings_range);
    TEST(test_stringdrop_newline);
    TEST(test_stringdrop_matches_drop_characters);
    TEST(test_stringdrop_out_of_bounds);
    TEST(test_stringdrop_out_of_bounds_negative);
    TEST(test_stringdrop_single_out_of_bounds);
    TEST(test_stringdrop_range_out_of_bounds);
    TEST(test_stringdrop_step_zero);
    TEST(test_stringdrop_unevaluated_symbolic);
    TEST(test_stringdrop_unevaluated_symbolic_spec);
    TEST(test_stringdrop_no_args);
    TEST(test_stringdrop_one_arg);
    TEST(test_stringdrop_three_args);
    TEST(test_stringdrop_empty_string_zero);
    TEST(test_stringdrop_empty_string_upto);

    /* StringReverse tests */
    TEST(test_stringreverse_basic);
    TEST(test_stringreverse_single_char);
    TEST(test_stringreverse_empty);
    TEST(test_stringreverse_palindrome);
    TEST(test_stringreverse_with_spaces);
    TEST(test_stringreverse_digits);
    TEST(test_stringreverse_special_chars);
    TEST(test_stringreverse_long);
    TEST(test_stringreverse_involution);
    TEST(test_stringreverse_listable);
    TEST(test_stringreverse_listable_empty_elements);
    TEST(test_stringreverse_of_stringjoin);
    TEST(test_stringreverse_symbolic);
    TEST(test_stringreverse_integer_arg);
    TEST(test_stringreverse_no_args);
    TEST(test_stringreverse_two_args);

    printf("All string tests passed!\n");
    return 0;
}
