/*
 * test_files.c — Unit tests for the filesystem builtins in
 * src/files.c.
 *
 * Covers:
 *   FileExistsQ  — true/false on a freshly created scratch file,
 *                  on directories, on dangling symlinks, plus the
 *                  unevaluated-on-bad-args contract.
 *   FileExtension — every Mathematica edge case from the docstring
 *                   (no extension, ends-with-dot, leading-dot,
 *                   nested extensions, directory specifications,
 *                   pure-directory inputs).
 *   FileBaseName  — symmetric coverage with FileExtension, plus the
 *                   "split off last extension only" behaviour.
 *
 * All scratch paths live under /tmp keyed by pid + tag so concurrent
 * test binaries don't collide; every test cleans up the files it
 * creates, and bad-arity / wrong-type inputs are checked to verify
 * the unevaluated path is taken (and therefore that res is not
 * freed twice).
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
#include <sys/stat.h>

/* Unique per-test scratch path under /tmp.  We pre-unlink so a
 * previous abnormal run doesn't leave state behind. */
static void scratch_path(char* out, size_t cap, const char* tag) {
    snprintf(out, cap, "/tmp/Mathilda_files_%d_%s", (int)getpid(), tag);
    unlink(out);
    rmdir(out);
}

/* ===== FileExistsQ ===== */

void test_fileexistsq_true_for_existing_file(void) {
    char path[256];
    scratch_path(path, sizeof(path), "exists");

    FILE* fp = fopen(path, "w");
    ASSERT(fp != NULL);
    fputs("x", fp);
    fclose(fp);

    char input[512];
    snprintf(input, sizeof(input), "FileExistsQ[\"%s\"]", path);
    assert_eval_eq(input, "True", 0);

    unlink(path);
}

void test_fileexistsq_false_for_missing_file(void) {
    char path[256];
    scratch_path(path, sizeof(path), "missing");
    /* scratch_path already unlinks it; we just trust it's gone. */

    char input[512];
    snprintf(input, sizeof(input), "FileExistsQ[\"%s\"]", path);
    assert_eval_eq(input, "False", 0);
}

/* "FileExistsQ tests for files, directories, or any other filesystem
 * objects." — verify directories count. */
void test_fileexistsq_true_for_directory(void) {
    char path[256];
    scratch_path(path, sizeof(path), "dir");
    ASSERT(mkdir(path, 0700) == 0);

    char input[512];
    snprintf(input, sizeof(input), "FileExistsQ[\"%s\"]", path);
    assert_eval_eq(input, "True", 0);

    rmdir(path);
}

/* A dangling symlink IS a filesystem object — lstat() succeeds, so
 * FileExistsQ must return True even though the target is missing. */
void test_fileexistsq_true_for_dangling_symlink(void) {
    char link_path[256];
    char target[256];
    scratch_path(link_path, sizeof(link_path), "danglink");
    scratch_path(target,    sizeof(target),    "danglinktgt");

    /* Make sure the target does NOT exist, then create the link. */
    if (symlink(target, link_path) != 0) {
        /* symlink() can fail on filesystems that don't support it
         * (rare on /tmp).  Skip rather than fail — the surrounding
         * coverage is plenty. */
        fprintf(stderr, "skip: symlink not supported on /tmp\n");
        return;
    }

    char input[512];
    snprintf(input, sizeof(input), "FileExistsQ[\"%s\"]", link_path);
    assert_eval_eq(input, "True", 0);

    unlink(link_path);
}

/* "" is not a valid path; lstat returns ENOENT.  False, no crash. */
void test_fileexistsq_empty_string(void) {
    assert_eval_eq("FileExistsQ[\"\"]", "False", 0);
}

/* Bad-args contract: the call must remain unevaluated (NULL return
 * from the builtin) when given the wrong arity or wrong types.  This
 * also catches a double-free bug — if the builtin freed res before
 * returning NULL the next reference to it (the resulting symbolic
 * form) would crash. */
void test_fileexistsq_unevaluated_on_bad_args(void) {
    assert_eval_eq("FileExistsQ[]",        "FileExistsQ[]",        0);
    assert_eval_eq("FileExistsQ[42]",      "FileExistsQ[42]",      0);
    assert_eval_eq("FileExistsQ[a, b]",    "FileExistsQ[a, b]",    0);
    /* Symbolic argument flows through — important for callers that
     * build FileExistsQ[x] expressions before x is bound. */
    assert_eval_eq("FileExistsQ[x]",       "FileExistsQ[x]",       0);
}

void test_fileexistsq_protected(void) {
    assert_eval_eq("MemberQ[Attributes[FileExistsQ], Protected]", "True", 0);
}

/* ===== FileExtension ===== */

void test_fileextension_simple(void) {
    assert_eval_eq("FileExtension[\"file.txt\"]", "\"txt\"", 0);
}

void test_fileextension_only_last(void) {
    /* "If there are multiple endings to a file name, separated by .,
     * FileExtension gives only the last one." */
    assert_eval_eq("FileExtension[\"file.tar.gz\"]", "\"gz\"", 0);
}

void test_fileextension_none(void) {
    assert_eval_eq("FileExtension[\"file\"]", "\"\"", 0);
}

/* "...gives \"\" if ... ends with a . character." */
void test_fileextension_trailing_dot(void) {
    assert_eval_eq("FileExtension[\"file.\"]", "\"\"", 0);
}

/* Hidden-style names with a leading dot have no extension; the
 * leading dot is part of the base name. */
void test_fileextension_leading_dot_only(void) {
    assert_eval_eq("FileExtension[\".bashrc\"]", "\"\"", 0);
}

void test_fileextension_ignores_directory(void) {
    /* "FileExtension ignores any directory specification." */
    assert_eval_eq("FileExtension[\"/foo/bar/file.txt\"]", "\"txt\"", 0);
    /* Dots in the directory part don't count. */
    assert_eval_eq("FileExtension[\"/a.b/c.d/file.txt\"]", "\"txt\"", 0);
    /* No extension on the leaf even though the directory has dots. */
    assert_eval_eq("FileExtension[\"/a.b/file\"]", "\"\"", 0);
}

void test_fileextension_directory_name(void) {
    /* "...if the file name has the form of a directory name" — the
     * leaf component is empty, so there is no extension. */
    assert_eval_eq("FileExtension[\"/foo/bar/\"]", "\"\"", 0);
    assert_eval_eq("FileExtension[\"\"]",        "\"\"", 0);
}

void test_fileextension_unevaluated_on_bad_args(void) {
    assert_eval_eq("FileExtension[]",     "FileExtension[]",     0);
    assert_eval_eq("FileExtension[42]",   "FileExtension[42]",   0);
    assert_eval_eq("FileExtension[a, b]", "FileExtension[a, b]", 0);
    assert_eval_eq("FileExtension[x]",    "FileExtension[x]",    0);
}

void test_fileextension_protected(void) {
    assert_eval_eq("MemberQ[Attributes[FileExtension], Protected]", "True", 0);
}

/* ===== FileBaseName ===== */

void test_filebasename_simple(void) {
    assert_eval_eq("FileBaseName[\"file.txt\"]", "\"file\"", 0);
}

/* "FileBaseName[\"file.tar.gz\"] gives \"file.tar\"" — only the
 * final extension is split off. */
void test_filebasename_drops_only_last_extension(void) {
    assert_eval_eq("FileBaseName[\"file.tar.gz\"]", "\"file.tar\"", 0);
}

void test_filebasename_no_extension(void) {
    assert_eval_eq("FileBaseName[\"file\"]", "\"file\"", 0);
}

/* "file." has no extension, so the trailing dot remains. */
void test_filebasename_trailing_dot(void) {
    assert_eval_eq("FileBaseName[\"file.\"]", "\"file.\"", 0);
}

/* ".bashrc" has no extension, so it's its own base name. */
void test_filebasename_leading_dot_only(void) {
    assert_eval_eq("FileBaseName[\".bashrc\"]", "\".bashrc\"", 0);
}

void test_filebasename_drops_directory(void) {
    /* "FileBaseName drops all directory specifications." */
    assert_eval_eq("FileBaseName[\"/foo/bar/file.txt\"]", "\"file\"", 0);
    assert_eval_eq("FileBaseName[\"/foo/bar/file\"]",     "\"file\"", 0);
    /* Directory dots don't bleed into the leaf base name. */
    assert_eval_eq("FileBaseName[\"/a.b/c.d/file.txt\"]", "\"file\"", 0);
}

/* Pure-directory inputs have an empty leaf, so the base name is "". */
void test_filebasename_directory_input(void) {
    assert_eval_eq("FileBaseName[\"/foo/bar/\"]", "\"\"", 0);
    assert_eval_eq("FileBaseName[\"\"]",         "\"\"", 0);
}

void test_filebasename_unevaluated_on_bad_args(void) {
    assert_eval_eq("FileBaseName[]",     "FileBaseName[]",     0);
    assert_eval_eq("FileBaseName[42]",   "FileBaseName[42]",   0);
    assert_eval_eq("FileBaseName[a, b]", "FileBaseName[a, b]", 0);
    assert_eval_eq("FileBaseName[x]",    "FileBaseName[x]",    0);
}

void test_filebasename_protected(void) {
    assert_eval_eq("MemberQ[Attributes[FileBaseName], Protected]", "True", 0);
}

/* Cross-check: FileBaseName + "." + FileExtension reconstructs the
 * filename component when an extension is present. */
void test_filebasename_extension_round_trip(void) {
    assert_eval_eq("StringJoin[FileBaseName[\"a/b/c.d\"], \".\", FileExtension[\"a/b/c.d\"]]",
                   "\"c.d\"", 0);
}

int main(void) {
    symtab_init();
    core_init();

    /* FileExistsQ */
    TEST(test_fileexistsq_true_for_existing_file);
    TEST(test_fileexistsq_false_for_missing_file);
    TEST(test_fileexistsq_true_for_directory);
    TEST(test_fileexistsq_true_for_dangling_symlink);
    TEST(test_fileexistsq_empty_string);
    TEST(test_fileexistsq_unevaluated_on_bad_args);
    TEST(test_fileexistsq_protected);

    /* FileExtension */
    TEST(test_fileextension_simple);
    TEST(test_fileextension_only_last);
    TEST(test_fileextension_none);
    TEST(test_fileextension_trailing_dot);
    TEST(test_fileextension_leading_dot_only);
    TEST(test_fileextension_ignores_directory);
    TEST(test_fileextension_directory_name);
    TEST(test_fileextension_unevaluated_on_bad_args);
    TEST(test_fileextension_protected);

    /* FileBaseName */
    TEST(test_filebasename_simple);
    TEST(test_filebasename_drops_only_last_extension);
    TEST(test_filebasename_no_extension);
    TEST(test_filebasename_trailing_dot);
    TEST(test_filebasename_leading_dot_only);
    TEST(test_filebasename_drops_directory);
    TEST(test_filebasename_directory_input);
    TEST(test_filebasename_unevaluated_on_bad_args);
    TEST(test_filebasename_protected);
    TEST(test_filebasename_extension_round_trip);

    printf("All tests passed!\n");
    return 0;
}
