/*
 * test_streams.c — Unit tests for the src/io stream layer: Read, OpenRead,
 * OpenWrite, OpenAppend, Write, WriteString, Close, Streams, StreamPosition,
 * SetStreamPosition, plus the InputStream/OutputStream/File objects.
 *
 * Each test writes a scratch file, then evaluates one program that binds the
 * path to the symbol `p` and drives the stream through `p`.  The stream
 * registry and symbol table are process-global, so successive Read/Write calls
 * inside one program share a current point.  Assertions compare with === (SameQ)
 * so they do not depend on how reals/strings print; assert_eval_eq frees the
 * parse tree, the result, and the printed string, so the suite is leak-clean by
 * construction (open streams are freed by Close, or by the atexit teardown).
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
    snprintf(out, cap, "/tmp/Mathilda_streams_%d_%s.txt", (int)getpid(), tag);
    unlink(out);
}

static void write_file(const char* path, const char* content) {
    FILE* fp = fopen(path, "wb");
    ASSERT(fp != NULL);
    if (content && *content) fputs(content, fp);
    fclose(fp);
}

/* Write `content` to a fresh scratch file, evaluate `prog` (which contains one
 * %s for the path — used to bind the Mathilda symbol p), assert == `expected`. */
static void check(const char* tag, const char* content,
                  const char* prog, const char* expected) {
    char path[256];
    scratch_path(path, sizeof(path), tag);
    write_file(path, content);
    char input[2048];
    snprintf(input, sizeof(input), prog, path);
    assert_eval_eq(input, expected, 0);
    unlink(path);
}

/* ===== OpenRead / stream objects ===== */

void test_openread_returns_inputstream(void) {
    check("or", "1 2 3\n",
          "p=\"%s\"; s=OpenRead[p]; h=Head[s]; Close[s]; h === InputStream", "True");
}

void test_openread_keeps_filename(void) {
    check("orname", "x\n",
          "p=\"%s\"; s=OpenRead[p]; q=StringQ[First[s]]; Close[s]; q", "True");
}

/* ===== Successive Read advances the current point ===== */

void test_read_advances(void) {
    check("adv", "1 2 3\n",
          "p=\"%s\"; s=OpenRead[p]; "
          "r={Read[s,Number],Read[s,Number],Read[s,Number],Read[s,Number]}; "
          "Close[s]; r === {1, 2, 3, EndOfFile}", "True");
}

void test_read_default_is_expression(void) {
    check("def", "1+1\n2*3\n",
          "p=\"%s\"; s=OpenRead[p]; r={Read[s],Read[s]}; Close[s]; r === {2, 6}",
          "True");
}

/* ===== Type-sequence grouping and nested structures ===== */

void test_read_group_pair(void) {
    check("pair", "1 2 3\n",
          "p=\"%s\"; s=OpenRead[p]; "
          "r={Read[s,{Number,Number}], Read[s,{Number,Number}]}; Close[s]; "
          "r === {{1, 2}, {3, EndOfFile}}", "True");
}

void test_read_nested_matrix(void) {
    check("mat", "1 2 3 4\n",
          "p=\"%s\"; s=OpenRead[p]; r=Read[s,{{Number,Number},{Number,Number}}]; "
          "Close[s]; r === {{1, 2}, {3, 4}}", "True");
}

void test_read_hold_expression_stays_held(void) {
    check("hold", "1+1\n",
          "p=\"%s\"; s=OpenRead[p]; r=Read[s,Hold[Expression]]; Close[s]; "
          "r === Hold[1 + 1]", "True");
}

void test_read_expression_evaluates(void) {
    check("evalx", "1+1\n",
          "p=\"%s\"; s=OpenRead[p]; r=Read[s,Expression]; Close[s]; r === 2", "True");
}

void test_read_arbitrary_head(void) {
    check("head", "3 4\n",
          "p=\"%s\"; s=OpenRead[p]; r=Read[s,ftest[Number,Number]]; Close[s]; "
          "r === ftest[3, 4]", "True");
}

/* ===== Per-type single reads ===== */

void test_read_byte(void) {
    check("byte", "ABC",
          "p=\"%s\"; s=OpenRead[p]; r={Read[s,Byte],Read[s,Byte],Read[s,Byte]}; "
          "Close[s]; r === {65, 66, 67}", "True");
}

void test_read_character(void) {
    check("char", "ab",
          "p=\"%s\"; s=OpenRead[p]; r={Read[s,Character],Read[s,Character]}; "
          "Close[s]; r === {\"a\", \"b\"}", "True");
}

void test_read_word(void) {
    check("word", "foo bar\n",
          "p=\"%s\"; s=OpenRead[p]; r={Read[s,Word],Read[s,Word]}; Close[s]; "
          "r === {\"foo\", \"bar\"}", "True");
}

void test_read_string_line(void) {
    check("line", "hello\nworld\n",
          "p=\"%s\"; s=OpenRead[p]; r={Read[s,String],Read[s,String]}; Close[s]; "
          "r === {\"hello\", \"world\"}", "True");
}

void test_read_record_with_option(void) {
    check("rec", "x;y;z",
          "p=\"%s\"; s=OpenRead[p]; "
          "r={Read[s,Record,RecordSeparators->{\";\"}], "
          "   Read[s,Record,RecordSeparators->{\";\"}]}; Close[s]; "
          "r === {\"x\", \"y\"}", "True");
}

void test_read_number_scientific(void) {
    check("sci", "2e5 2.e5\n",
          "p=\"%s\"; s=OpenRead[p]; r={Read[s,Number],Read[s,Number]}; Close[s]; "
          "r === {200000., 200000.}", "True");
}

void test_read_real_always_approximate(void) {
    check("real", "3\n",
          "p=\"%s\"; s=OpenRead[p]; r=Read[s,Real]; Close[s]; r === 3.", "True");
}

/* ===== EOF and error paths ===== */

void test_read_eof_returns_endoffile(void) {
    check("eof", "",
          "p=\"%s\"; s=OpenRead[p]; r=Read[s,Number]; Close[s]; r === EndOfFile",
          "True");
}

void test_read_group_eof_bare(void) {
    check("eofg", "",
          "p=\"%s\"; s=OpenRead[p]; r=Read[s,{Number,Number}]; Close[s]; "
          "r === EndOfFile", "True");
}

void test_read_malformed_number(void) {
    /* Prints Read::readn and returns $Failed. */
    check("bad", "abc\n",
          "p=\"%s\"; s=OpenRead[p]; r=Read[s,Number]; Close[s]; r === $Failed",
          "True");
}

void test_read_after_close_fails(void) {
    /* Prints Read::openx and returns $Failed. */
    check("closed", "1\n",
          "p=\"%s\"; s=OpenRead[p]; Close[s]; Read[s,Number] === $Failed", "True");
}

void test_read_missing_file(void) {
    assert_eval_eq(
        "Read[\"/no/such/dir/Mathilda_streams_missing.txt\", Number] === $Failed",
        "True", 0);
}

/* ===== Options ===== */

void test_read_wordseparators_option(void) {
    check("wsep", "a,b\n",
          "p=\"%s\"; s=OpenRead[p]; "
          "r={Read[s,Word,WordSeparators->{\",\"}], "
          "   Read[s,Word,WordSeparators->{\",\"}]}; Close[s]; "
          "r === {\"a\", \"b\"}", "True");
}

void test_read_setoptions(void) {
    char path[256];
    scratch_path(path, sizeof(path), "setopt");
    write_file(path, "a,b\n");
    char input[1024];
    assert_eval_eq("SetOptions[Read, WordSeparators -> {\",\"}]; 42", "42", 0);
    snprintf(input, sizeof(input),
             "p=\"%s\"; s=OpenRead[p]; r={Read[s,Word],Read[s,Word]}; Close[s]; "
             "r === {\"a\", \"b\"}", path);
    assert_eval_eq(input, "True", 0);
    assert_eval_eq("SetOptions[Read, WordSeparators -> {\" \", \"\\t\"}]; 42",
                   "42", 0);
    unlink(path);
}

/* ===== File["..."] form and filename persistence ===== */

void test_read_file_wrapper(void) {
    check("filew", "1 2\n",
          "p=\"%s\"; s=OpenRead[File[p]]; r=Read[s,Number]; Close[s]; r === 1",
          "True");
}

void test_read_filename_persists(void) {
    /* No explicit OpenRead: the named file is auto-opened and stays open, so a
     * second Read advances past the first. */
    check("persist", "1 2 3\n",
          "p=\"%s\"; a=Read[p,Number]; b=Read[p,Number]; Close[p]; {a,b} === {1, 2}",
          "True");
}

/* ===== Output: Write / WriteString / OpenAppend ===== */

void test_write_roundtrip(void) {
    check("wrt", "",
          "p=\"%s\"; s=OpenWrite[p]; Write[s,1+1]; Write[s,a+b]; Close[s]; "
          "ReadList[p] === {2, a + b}", "True");
}

void test_writestring_raw(void) {
    check("wstr", "",
          "p=\"%s\"; s=OpenWrite[p]; WriteString[s,\"hello\",\" \",\"world\"]; "
          "Close[s]; ReadList[p,String] === {\"hello world\"}", "True");
}

void test_openappend(void) {
    check("app", "",
          "p=\"%s\"; s=OpenWrite[p]; WriteString[s,\"a\\n\"]; Close[s]; "
          "t=OpenAppend[p]; WriteString[t,\"b\\n\"]; Close[t]; "
          "ReadList[p,String] === {\"a\", \"b\"}", "True");
}

/* ===== Streams / positioning ===== */

void test_streams_lists_open(void) {
    check("slist", "x\n",
          "p=\"%s\"; s=OpenRead[p]; n=Length[Streams[p]]; Close[s]; n === 1",
          "True");
}

void test_streams_empty_after_close(void) {
    check("sclose", "x\n",
          "p=\"%s\"; s=OpenRead[p]; Close[s]; Streams[p] === {}", "True");
}

void test_streamposition_input(void) {
    check("spos", "abcdef",
          "p=\"%s\"; s=OpenRead[p]; Read[s,Character]; k=StreamPosition[s]; "
          "Close[s]; k === 1", "True");
}

void test_setstreamposition_input(void) {
    check("sset", "abcdef",
          "p=\"%s\"; s=OpenRead[p]; SetStreamPosition[s,4]; c=Read[s,Character]; "
          "Close[s]; c === \"e\"", "True");
}

void test_setstreamposition_infinity(void) {
    check("sinf", "abcdef",
          "p=\"%s\"; s=OpenRead[p]; SetStreamPosition[s,Infinity]; "
          "r=Read[s,Character]; Close[s]; r === EndOfFile", "True");
}

void test_streamposition_output(void) {
    check("opos", "",
          "p=\"%s\"; s=OpenWrite[p]; WriteString[s,\"abc\"]; k=StreamPosition[s]; "
          "Close[s]; k === 3", "True");
}

/* ===== ReadList / Read consistency ===== */

void test_readlist_still_flat(void) {
    check("rl", "1 2 3\n", "ReadList[\"%s\", Number] === {1, 2, 3}", "True");
}

void test_readlist_from_open_stream(void) {
    /* One Read, then ReadList continues from the same current point. */
    check("rlstream", "1 2 3\n",
          "p=\"%s\"; s=OpenRead[p]; a=Read[s,Number]; rest=ReadList[s,Number]; "
          "Close[s]; {a, rest} === {1, {2, 3}}", "True");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_openread_returns_inputstream);
    TEST(test_openread_keeps_filename);

    TEST(test_read_advances);
    TEST(test_read_default_is_expression);

    TEST(test_read_group_pair);
    TEST(test_read_nested_matrix);
    TEST(test_read_hold_expression_stays_held);
    TEST(test_read_expression_evaluates);
    TEST(test_read_arbitrary_head);

    TEST(test_read_byte);
    TEST(test_read_character);
    TEST(test_read_word);
    TEST(test_read_string_line);
    TEST(test_read_record_with_option);
    TEST(test_read_number_scientific);
    TEST(test_read_real_always_approximate);

    TEST(test_read_eof_returns_endoffile);
    TEST(test_read_group_eof_bare);
    TEST(test_read_malformed_number);
    TEST(test_read_after_close_fails);
    TEST(test_read_missing_file);

    TEST(test_read_wordseparators_option);
    TEST(test_read_setoptions);

    TEST(test_read_file_wrapper);
    TEST(test_read_filename_persists);

    TEST(test_write_roundtrip);
    TEST(test_writestring_raw);
    TEST(test_openappend);

    TEST(test_streams_lists_open);
    TEST(test_streams_empty_after_close);
    TEST(test_streamposition_input);
    TEST(test_setstreamposition_input);
    TEST(test_setstreamposition_infinity);
    TEST(test_streamposition_output);

    TEST(test_readlist_still_flat);
    TEST(test_readlist_from_open_stream);

    printf("All stream tests passed.\n");
    return 0;
}
