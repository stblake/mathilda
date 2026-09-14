/*
 * test_readlist.c — Unit tests for ReadList (src/io/readlist.c).
 *
 * Covers every read type (Byte, Character, Expression, Number, Real, Record,
 * String, Word), type-sequence grouping with EndOfFile padding, the n-object
 * limit, the separator options (RecordSeparators, WordSeparators, TokenWords,
 * NullRecords, NullWords), SetOptions integration, and the error/edge paths
 * (missing file, malformed number, empty file).
 *
 * Assertions compare with === (SameQ) against a parsed literal so the checks
 * are independent of how reals/strings print. assert_eval_eq frees the parse
 * tree, the evaluated result, and the printed string, so the suite is
 * leak-clean by construction.
 */

#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "parse.h"
#include "print.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/* Unique-ish per-test scratch path under /tmp (pid + tag), pre-unlinked. */
static void scratch_path(char* out, size_t cap, const char* tag) {
    snprintf(out, cap, "/tmp/Mathilda_readlist_%d_%s.txt", (int)getpid(), tag);
    unlink(out);
}

static void write_file(const char* path, const char* content) {
    FILE* fp = fopen(path, "wb");
    ASSERT(fp != NULL);
    if (content && *content) fputs(content, fp);
    fclose(fp);
}

/* Write `content` to a fresh scratch file, evaluate `expr_fmt` (which contains
 * one %s for the path), and assert it equals `expected`. */
static void check(const char* tag, const char* content,
                  const char* expr_fmt, const char* expected) {
    char path[256];
    scratch_path(path, sizeof(path), tag);
    write_file(path, content);
    char input[1024];
    snprintf(input, sizeof(input), expr_fmt, path);
    assert_eval_eq(input, expected, 0);
    unlink(path);
}

/* ===== Expression (default type) ===== */

void test_default_reads_expressions(void) {
    check("expr", "1+1\n2*3\nx+x\n",
          "ReadList[\"%s\"] === {2, 6, 2 x}", "True");
}

void test_explicit_expression_type(void) {
    check("expr2", "a; b;\nc\n",
          "ReadList[\"%s\", Expression] === {a, b, c}", "True");
}

void test_empty_file_default(void) {
    check("empty", "", "ReadList[\"%s\"] === {}", "True");
}

void test_empty_file_number(void) {
    check("emptyn", "", "ReadList[\"%s\", Number] === {}", "True");
}

/* ===== Number / Real ===== */

void test_number_int_and_real(void) {
    check("nums", "1 2 3\n4.5 6\n",
          "ReadList[\"%s\", Number] === {1, 2, 3, 4.5, 6}", "True");
}

void test_number_scientific_and_fortran(void) {
    /* 2e5 / 2.e5 have E-notation -> real; 2.*10^5 evaluates to a real. */
    check("sci", "2e5 2.e5 1.5e-3\n",
          "ReadList[\"%s\", Number] === {200000., 200000., 0.0015}", "True");
    check("fort", "2.*10^5\n",
          "ReadList[\"%s\", Number] === {200000.}", "True");
}

void test_real_always_approximate(void) {
    check("reals", "1 2 3\n",
          "ReadList[\"%s\", Real] === {1., 2., 3.}", "True");
}

void test_number_head_is_integer(void) {
    check("nhead", "42\n",
          "Head[First[ReadList[\"%s\", Number]]] === Integer", "True");
}

void test_real_head_is_real(void) {
    check("rhead", "42\n",
          "Head[First[ReadList[\"%s\", Real]]] === Real", "True");
}

void test_malformed_number_yields_failed(void) {
    /* Prints ReadList::readn and inserts $Failed, then continues. */
    check("bad", "1 abc 3\n",
          "ReadList[\"%s\", Number] === {1, $Failed, 3}", "True");
}

/* ===== Word / Character / Byte ===== */

void test_word_splits_on_whitespace(void) {
    check("words", "the quick brown\nfox jumps\n",
          "ReadList[\"%s\", Word] === {\"the\", \"quick\", \"brown\", \"fox\", \"jumps\"}",
          "True");
}

void test_word_head_is_string(void) {
    check("whead", "alpha beta\n",
          "Head[First[ReadList[\"%s\", Word]]] === String", "True");
}

void test_character_reads_one_char_strings(void) {
    check("chars", "abc",
          "ReadList[\"%s\", Character] === {\"a\", \"b\", \"c\"}", "True");
}

void test_byte_reads_integer_codes(void) {
    check("bytes", "ABC",
          "ReadList[\"%s\", Byte] === {65, 66, 67}", "True");
}

/* ===== Record / String ===== */

void test_string_reads_lines(void) {
    check("lines", "1+1\n2*3\nx+x\n",
          "ReadList[\"%s\", String] === {\"1+1\", \"2*3\", \"x+x\"}", "True");
}

void test_record_default_skips_empty(void) {
    check("recs", "x\n\ny\n",
          "ReadList[\"%s\", Record] === {\"x\", \"y\"}", "True");
}

void test_record_nullrecords_keeps_empty(void) {
    check("recsn", "x\n\ny\n",
          "ReadList[\"%s\", Record, NullRecords -> True] === {\"x\", \"\", \"y\"}",
          "True");
}

/* ===== Type-sequence grouping ===== */

void test_group_pairs(void) {
    check("pairs", "a 1\nb 2\nc 3\n",
          "ReadList[\"%s\", {Word, Number}] === {{\"a\", 1}, {\"b\", 2}, {\"c\", 3}}",
          "True");
}

void test_group_endoffile_padding(void) {
    check("three", "1 2 3\n",
          "ReadList[\"%s\", {Number, Number}] === {{1, 2}, {3, EndOfFile}}", "True");
}

/* ===== n limit ===== */

void test_n_limits_flat(void) {
    check("nlim", "1 2 3 4 5\n",
          "ReadList[\"%s\", Number, 2] === {1, 2}", "True");
}

void test_n_limits_groups(void) {
    check("nglim", "a 1\nb 2\nc 3\n",
          "ReadList[\"%s\", {Word, Number}, 2] === {{\"a\", 1}, {\"b\", 2}}", "True");
}

/* ===== Options ===== */

void test_wordseparators_option(void) {
    check("csv", "a,b,c\n",
          "ReadList[\"%s\", Word, WordSeparators -> {\",\"}] === {\"a\", \"b\", \"c\"}",
          "True");
}

void test_nullwords_option(void) {
    check("dbl", "a  b\n",
          "ReadList[\"%s\", Word, NullWords -> True] === {\"a\", \"\", \"b\"}", "True");
}

void test_tokenwords_option(void) {
    check("plus", "a+b\n",
          "ReadList[\"%s\", Word, TokenWords -> {\"+\"}] === {\"a\", \"+\", \"b\"}",
          "True");
}

void test_recordseparators_option(void) {
    check("semi", "x;y;z",
          "ReadList[\"%s\", Record, RecordSeparators -> {\";\"}] === {\"x\", \"y\", \"z\"}",
          "True");
}

/* ===== SetOptions integration ===== */

void test_setoptions_default(void) {
    char path[256];
    scratch_path(path, sizeof(path), "setopt");
    write_file(path, "a,b,c\n");
    char input[1024];
    /* Change the default WordSeparators, read without an inline option, reset. */
    assert_eval_eq("SetOptions[ReadList, WordSeparators -> {\",\"}]; 42", "42", 0);
    snprintf(input, sizeof(input),
             "ReadList[\"%s\", Word] === {\"a\", \"b\", \"c\"}", path);
    assert_eval_eq(input, "True", 0);
    assert_eval_eq("SetOptions[ReadList, WordSeparators -> {\" \", \"\\t\"}]; 42",
                   "42", 0);
    unlink(path);
}

/* ===== Error paths ===== */

void test_missing_file_returns_failed(void) {
    /* Prints ReadList::noopen and returns $Failed. */
    assert_eval_eq("ReadList[\"/no/such/dir/Mathilda_readlist_missing.txt\"] === $Failed",
                   "True", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_default_reads_expressions);
    TEST(test_explicit_expression_type);
    TEST(test_empty_file_default);
    TEST(test_empty_file_number);

    TEST(test_number_int_and_real);
    TEST(test_number_scientific_and_fortran);
    TEST(test_real_always_approximate);
    TEST(test_number_head_is_integer);
    TEST(test_real_head_is_real);
    TEST(test_malformed_number_yields_failed);

    TEST(test_word_splits_on_whitespace);
    TEST(test_word_head_is_string);
    TEST(test_character_reads_one_char_strings);
    TEST(test_byte_reads_integer_codes);

    TEST(test_string_reads_lines);
    TEST(test_record_default_skips_empty);
    TEST(test_record_nullrecords_keeps_empty);

    TEST(test_group_pairs);
    TEST(test_group_endoffile_padding);

    TEST(test_n_limits_flat);
    TEST(test_n_limits_groups);

    TEST(test_wordseparators_option);
    TEST(test_nullwords_option);
    TEST(test_tokenwords_option);
    TEST(test_recordseparators_option);

    TEST(test_setoptions_default);

    TEST(test_missing_file_returns_failed);

    printf("All ReadList tests passed.\n");
    return 0;
}
