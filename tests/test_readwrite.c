/*
 * test_readwrite.c — Unit tests for Get / Put / PutAppend.
 *
 * Covers:
 *   - parsing of expr >> filename, expr >> "filename", expr >>> filename
 *   - that bare and quoted filenames lower to the same Put/PutAppend AST
 *   - generation of files via Put (truncating semantics)
 *   - empty-file creation via Put["filename"]
 *   - append semantics of PutAppend
 *   - Get round-tripping a value written by Put
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

/* Build a unique-ish per-test scratch path under /tmp.  The CTest
 * working directory varies, and tests may run concurrently, so we
 * compose pid + a per-test tag to avoid collisions. */
static void scratch_path(char* out, size_t cap, const char* tag) {
    snprintf(out, cap, "/tmp/Mathilda_rw_%d_%s.txt", (int)getpid(), tag);
    unlink(out);  /* ignore errors — file may not exist */
}

static char* read_whole_file(const char* path) {
    FILE* fp = fopen(path, "rb");
    if (!fp) return NULL;
    fseek(fp, 0, SEEK_END);
    long sz = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    if (sz < 0) { fclose(fp); return NULL; }
    char* buf = malloc((size_t)sz + 1);
    if (!buf) { fclose(fp); return NULL; }
    size_t n = fread(buf, 1, (size_t)sz, fp);
    buf[n] = '\0';
    fclose(fp);
    return buf;
}

/* ===== Parser tests ===== */

/* `expr >> "file.txt"` parses to Put[expr, "file.txt"]. FullForm gives
 * a stable text representation we can compare. */
void test_parse_put_quoted_filename() {
    assert_eval_eq("FullForm[Hold[2 + 3 >> \"out.m\"]]",
                   "Hold[Put[Plus[2, 3], \"out.m\"]]", 0);
}

/* Bare filename desugars to a quoted-string second argument. */
void test_parse_put_bare_filename() {
    assert_eval_eq("FullForm[Hold[x >> myfile]]",
                   "Hold[Put[x, \"myfile\"]]", 0);
}

/* Filenames may contain dots, slashes, dashes, underscores. */
void test_parse_put_path_chars() {
    assert_eval_eq("FullForm[Hold[x >> ./sub_dir/file-1.m]]",
                   "Hold[Put[x, \"./sub_dir/file-1.m\"]]", 0);
}

void test_parse_putappend_quoted() {
    assert_eval_eq("FullForm[Hold[a >>> \"log.m\"]]",
                   "Hold[PutAppend[a, \"log.m\"]]", 0);
}

void test_parse_putappend_bare() {
    assert_eval_eq("FullForm[Hold[a >>> log]]",
                   "Hold[PutAppend[a, \"log\"]]", 0);
}

/* `;` has lower precedence (10) than `>>` (30), so semicolons split
 * statements normally even when one of them uses `>>`. */
void test_parse_put_in_compound() {
    assert_eval_eq("FullForm[Hold[a = 5; a >> \"f\"]]",
                   "Hold[CompoundExpression[Set[a, 5], Put[a, \"f\"]]]", 0);
}

/* `+` (310) binds tighter than `>>` (30): `a + b >> "f"` is
 * `Put[a + b, "f"]`, NOT `a + Put[b, "f"]`. */
void test_parse_put_lower_than_plus() {
    assert_eval_eq("FullForm[Hold[a + b >> \"f\"]]",
                   "Hold[Put[Plus[a, b], \"f\"]]", 0);
}

/* `>>>` must match before `>>`, which must match before `>`. Verify by
 * giving each a chance to misparse the leading characters. */
void test_parse_putappend_vs_put_lex() {
    /* `>>>` first */
    assert_eval_eq("FullForm[Hold[1 >>> file]]",
                   "Hold[PutAppend[1, \"file\"]]", 0);
    /* `>>` second */
    assert_eval_eq("FullForm[Hold[1 >> file]]",
                   "Hold[Put[1, \"file\"]]", 0);
    /* `>` still works as Greater for arithmetic comparison. */
    assert_eval_eq("FullForm[Hold[1 > 2]]",
                   "Hold[Greater[1, 2]]", 0);
}

/* ===== File generation tests ===== */

void test_put_writes_integer() {
    char path[256];
    scratch_path(path, sizeof(path), "writeint");

    char input[512];
    snprintf(input, sizeof(input), "Put[42, \"%s\"]", path);
    assert_eval_eq(input, "Null", 0);

    char* got = read_whole_file(path);
    ASSERT(got != NULL);
    ASSERT_STR_EQ(got, "42\n");
    free(got);
    unlink(path);
}

void test_put_writes_big_integer() {
    char path[256];
    scratch_path(path, sizeof(path), "bigint");

    char input[512];
    snprintf(input, sizeof(input), "Put[47!, \"%s\"]", path);
    assert_eval_eq(input, "Null", 0);

    char* got = read_whole_file(path);
    ASSERT(got != NULL);
    ASSERT_STR_EQ(got,
        "258623241511168180642964355153611979969197632389120000000000\n");
    free(got);
    unlink(path);
}

void test_put_overwrites_existing() {
    char path[256];
    scratch_path(path, sizeof(path), "overwrite");

    /* Pre-populate with junk. */
    FILE* fp = fopen(path, "w");
    ASSERT(fp != NULL);
    fputs("garbage that must be deleted\nmore garbage\n", fp);
    fclose(fp);

    char input[512];
    snprintf(input, sizeof(input), "Put[7, \"%s\"]", path);
    assert_eval_eq(input, "Null", 0);

    char* got = read_whole_file(path);
    ASSERT(got != NULL);
    ASSERT_STR_EQ(got, "7\n");
    free(got);
    unlink(path);
}

void test_put_multiple_expressions() {
    char path[256];
    scratch_path(path, sizeof(path), "multi");

    char input[512];
    snprintf(input, sizeof(input), "Put[1, 2, 3, \"%s\"]", path);
    assert_eval_eq(input, "Null", 0);

    char* got = read_whole_file(path);
    ASSERT(got != NULL);
    ASSERT_STR_EQ(got, "1\n2\n3\n");
    free(got);
    unlink(path);
}

void test_put_empty_file() {
    char path[256];
    scratch_path(path, sizeof(path), "empty");

    char input[512];
    snprintf(input, sizeof(input), "Put[\"%s\"]", path);
    assert_eval_eq(input, "Null", 0);

    char* got = read_whole_file(path);
    ASSERT(got != NULL);
    ASSERT_STR_EQ(got, "");
    free(got);
    unlink(path);
}

/* `expr >> filename` form must produce the same file as
 * `Put[expr, "filename"]`. */
void test_put_shorthand_writes_file() {
    char path[256];
    scratch_path(path, sizeof(path), "shorthand");

    char input[512];
    snprintf(input, sizeof(input), "100 >> \"%s\"", path);
    assert_eval_eq(input, "Null", 0);

    char* got = read_whole_file(path);
    ASSERT(got != NULL);
    ASSERT_STR_EQ(got, "100\n");
    free(got);
    unlink(path);
}

/* ===== PutAppend tests ===== */

void test_putappend_creates_file() {
    char path[256];
    scratch_path(path, sizeof(path), "appendnew");

    char input[512];
    snprintf(input, sizeof(input), "PutAppend[\"hi\", \"%s\"]", path);
    assert_eval_eq(input, "Null", 0);

    char* got = read_whole_file(path);
    ASSERT(got != NULL);
    ASSERT_STR_EQ(got, "\"hi\"\n");
    free(got);
    unlink(path);
}

/* PutAppend preserves prior contents, unlike Put. */
void test_putappend_preserves_existing() {
    char path[256];
    scratch_path(path, sizeof(path), "appendkeep");

    char put_input[512];
    snprintf(put_input, sizeof(put_input), "Put[1, \"%s\"]", path);
    assert_eval_eq(put_input, "Null", 0);

    char app_input[512];
    snprintf(app_input, sizeof(app_input), "PutAppend[2, \"%s\"]", path);
    assert_eval_eq(app_input, "Null", 0);

    char* got = read_whole_file(path);
    ASSERT(got != NULL);
    ASSERT_STR_EQ(got, "1\n2\n");
    free(got);
    unlink(path);
}

void test_putappend_shorthand() {
    char path[256];
    scratch_path(path, sizeof(path), "appendshort");

    /* Sequence mirrors the Mathematica reference example:
     *   FactorInteger[40320]  >>  factorizations
     *   FactorInteger[479001600] >>> factorizations
     * Use small literals here so we don't rely on factorisation output
     * formatting beyond the basic List printing. */
    char a[512], b[512];
    snprintf(a, sizeof(a), "{{2,7},{3,2}} >> \"%s\"",  path);
    snprintf(b, sizeof(b), "{{2,10},{3,5}} >>> \"%s\"", path);
    assert_eval_eq(a, "Null", 0);
    assert_eval_eq(b, "Null", 0);

    char* got = read_whole_file(path);
    ASSERT(got != NULL);
    ASSERT_STR_EQ(got, "{{2, 7}, {3, 2}}\n{{2, 10}, {3, 5}}\n");
    free(got);
    unlink(path);
}

/* PutAppend with multiple expressions writes them all, in order, to the
 * end of the file. */
void test_putappend_multiple_expressions() {
    char path[256];
    scratch_path(path, sizeof(path), "appendmulti");

    char p[512];
    snprintf(p, sizeof(p), "Put[0, \"%s\"]", path);
    assert_eval_eq(p, "Null", 0);

    char a[512];
    snprintf(a, sizeof(a), "PutAppend[10, 20, 30, \"%s\"]", path);
    assert_eval_eq(a, "Null", 0);

    char* got = read_whole_file(path);
    ASSERT(got != NULL);
    ASSERT_STR_EQ(got, "0\n10\n20\n30\n");
    free(got);
    unlink(path);
}

/* ===== Get round-trip test ===== */

/* What Put writes Get must be able to read back. */
void test_put_then_get_round_trip() {
    char path[256];
    scratch_path(path, sizeof(path), "roundtrip");

    char input[512];
    snprintf(input, sizeof(input), "Put[1 + 2 + 3, \"%s\"]", path);
    assert_eval_eq(input, "Null", 0);

    char get_input[512];
    snprintf(get_input, sizeof(get_input), "Get[\"%s\"]", path);
    /* 1 + 2 + 3 evaluates to 6 before Put writes it. */
    assert_eval_eq(get_input, "6", 0);

    unlink(path);
}

/* ===== Attribute checks ===== */

void test_put_putappend_protected() {
    /* AttributeQ / Attributes returns the list; we just check it
     * contains Protected. */
    assert_eval_eq("MemberQ[Attributes[Put], Protected]", "True", 0);
    assert_eval_eq("MemberQ[Attributes[PutAppend], Protected]", "True", 0);
    assert_eval_eq("MemberQ[Attributes[Get], Protected]", "True", 0);
}

int main(void) {
    symtab_init();
    core_init();

    /* Parsing */
    TEST(test_parse_put_quoted_filename);
    TEST(test_parse_put_bare_filename);
    TEST(test_parse_put_path_chars);
    TEST(test_parse_putappend_quoted);
    TEST(test_parse_putappend_bare);
    TEST(test_parse_put_in_compound);
    TEST(test_parse_put_lower_than_plus);
    TEST(test_parse_putappend_vs_put_lex);

    /* File generation — Put */
    TEST(test_put_writes_integer);
    TEST(test_put_writes_big_integer);
    TEST(test_put_overwrites_existing);
    TEST(test_put_multiple_expressions);
    TEST(test_put_empty_file);
    TEST(test_put_shorthand_writes_file);

    /* File generation — PutAppend */
    TEST(test_putappend_creates_file);
    TEST(test_putappend_preserves_existing);
    TEST(test_putappend_shorthand);
    TEST(test_putappend_multiple_expressions);

    /* Round trip with Get */
    TEST(test_put_then_get_round_trip);

    /* Attributes */
    TEST(test_put_putappend_protected);

    printf("All tests passed!\n");
    return 0;
}
