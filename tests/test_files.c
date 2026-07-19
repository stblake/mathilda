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

#include <fcntl.h>
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

/* ===== FilePrint =====
 *
 * FilePrint writes to stdout rather than returning the content, so each
 * test redirects fd 1 to a scratch capture file, runs the Mathilda
 * expression, restores stdout, and compares the captured bytes against
 * an expected string.  The redirection has to live in a helper so the
 * file descriptors get cleaned up even when an assertion fails (we
 * keep that risk small by asserting only after restoration).
 */

/* Write `content` to a fresh path keyed by `tag`.  Returns a pointer
 * into a static buffer holding the resulting path — fine for tests, not
 * thread-safe (we don't run tests concurrently). */
static const char* fileprint_make_input(const char* tag, const char* content) {
    static char path[256];
    scratch_path(path, sizeof(path), tag);
    FILE* fp = fopen(path, "wb");
    ASSERT(fp != NULL);
    if (*content) {
        size_t n = strlen(content);
        ASSERT(fwrite(content, 1, n, fp) == n);
    }
    fclose(fp);
    return path;
}

/* Run a Mathilda expression with stdout redirected to a scratch file,
 * then return the captured bytes (caller frees).  Any non-NULL evaluator
 * result is also freed so the leak checker stays happy. */
static char* fileprint_eval_capture(const char* mathilda_expr) {
    char capture[256];
    scratch_path(capture, sizeof(capture), "stdout");

    fflush(stdout);
    int saved = dup(STDOUT_FILENO);
    ASSERT(saved >= 0);
    int fd = open(capture, O_WRONLY | O_CREAT | O_TRUNC, 0600);
    ASSERT(fd >= 0);
    ASSERT(dup2(fd, STDOUT_FILENO) >= 0);
    close(fd);

    Expr* parsed = parse_expression(mathilda_expr);
    Expr* result = (parsed != NULL) ? evaluate(parsed) : NULL;

    fflush(stdout);
    ASSERT(dup2(saved, STDOUT_FILENO) >= 0);
    close(saved);

    if (parsed) expr_free(parsed);
    if (result) expr_free(result);

    FILE* fp = fopen(capture, "rb");
    ASSERT(fp != NULL);
    fseek(fp, 0, SEEK_END);
    long sz = ftell(fp);
    fseek(fp, 0, SEEK_SET);
    char* buf = (char*)malloc((size_t)sz + 1);
    ASSERT(buf != NULL);
    size_t got = fread(buf, 1, (size_t)sz, fp);
    buf[got] = '\0';
    fclose(fp);
    unlink(capture);
    return buf;
}

/* Five-line fixture used by most of the selector tests. */
static const char* FP_FIVE = "alpha\nbeta\ngamma\ndelta\nepsilon\n";

void test_fileprint_full_file(void) {
    const char* path = fileprint_make_input("fp_full", FP_FIVE);
    char expr[512];
    snprintf(expr, sizeof(expr), "FilePrint[\"%s\"]", path);
    char* out = fileprint_eval_capture(expr);
    ASSERT_STR_EQ(out, FP_FIVE);
    free(out);
    unlink(path);
}

void test_fileprint_first_n(void) {
    const char* path = fileprint_make_input("fp_first", FP_FIVE);
    char expr[512];
    snprintf(expr, sizeof(expr), "FilePrint[\"%s\", 2]", path);
    char* out = fileprint_eval_capture(expr);
    ASSERT_STR_EQ(out, "alpha\nbeta\n");
    free(out);
    unlink(path);
}

/* n larger than the file size clamps to "print everything". */
void test_fileprint_first_n_clamps(void) {
    const char* path = fileprint_make_input("fp_first_big", FP_FIVE);
    char expr[512];
    snprintf(expr, sizeof(expr), "FilePrint[\"%s\", 100]", path);
    char* out = fileprint_eval_capture(expr);
    ASSERT_STR_EQ(out, FP_FIVE);
    free(out);
    unlink(path);
}

/* n == 0 is a valid call that simply emits nothing. */
void test_fileprint_zero_is_empty(void) {
    const char* path = fileprint_make_input("fp_zero", FP_FIVE);
    char expr[512];
    snprintf(expr, sizeof(expr), "FilePrint[\"%s\", 0]", path);
    char* out = fileprint_eval_capture(expr);
    ASSERT_STR_EQ(out, "");
    free(out);
    unlink(path);
}

void test_fileprint_last_n(void) {
    const char* path = fileprint_make_input("fp_last", FP_FIVE);
    char expr[512];
    snprintf(expr, sizeof(expr), "FilePrint[\"%s\", -2]", path);
    char* out = fileprint_eval_capture(expr);
    ASSERT_STR_EQ(out, "delta\nepsilon\n");
    free(out);
    unlink(path);
}

void test_fileprint_last_n_clamps(void) {
    const char* path = fileprint_make_input("fp_last_big", FP_FIVE);
    char expr[512];
    snprintf(expr, sizeof(expr), "FilePrint[\"%s\", -100]", path);
    char* out = fileprint_eval_capture(expr);
    ASSERT_STR_EQ(out, FP_FIVE);
    free(out);
    unlink(path);
}

void test_fileprint_span_basic(void) {
    const char* path = fileprint_make_input("fp_span", FP_FIVE);
    char expr[512];
    snprintf(expr, sizeof(expr), "FilePrint[\"%s\", 2;;4]", path);
    char* out = fileprint_eval_capture(expr);
    ASSERT_STR_EQ(out, "beta\ngamma\ndelta\n");
    free(out);
    unlink(path);
}

/* Negative endpoints inside the Span count from the end. */
void test_fileprint_span_negative_end(void) {
    const char* path = fileprint_make_input("fp_span_neg", FP_FIVE);
    char expr[512];
    snprintf(expr, sizeof(expr), "FilePrint[\"%s\", 2;;-2]", path);
    char* out = fileprint_eval_capture(expr);
    ASSERT_STR_EQ(out, "beta\ngamma\ndelta\n");
    free(out);
    unlink(path);
}

/* Span[m, n] with m > n and an implicit +1 step yields no output. */
void test_fileprint_span_empty(void) {
    const char* path = fileprint_make_input("fp_span_empty", FP_FIVE);
    char expr[512];
    snprintf(expr, sizeof(expr), "FilePrint[\"%s\", 4;;2]", path);
    char* out = fileprint_eval_capture(expr);
    ASSERT_STR_EQ(out, "");
    free(out);
    unlink(path);
}

void test_fileprint_span_step(void) {
    const char* path = fileprint_make_input("fp_span_step", FP_FIVE);
    char expr[512];
    snprintf(expr, sizeof(expr), "FilePrint[\"%s\", 1;;5;;2]", path);
    char* out = fileprint_eval_capture(expr);
    ASSERT_STR_EQ(out, "alpha\ngamma\nepsilon\n");
    free(out);
    unlink(path);
}

/* Negative step walks backwards.  Verifies output order matches the
 * iteration order, not the source order. */
void test_fileprint_span_negative_step(void) {
    const char* path = fileprint_make_input("fp_span_revstep", FP_FIVE);
    char expr[512];
    snprintf(expr, sizeof(expr), "FilePrint[\"%s\", 5;;1;;-1]", path);
    char* out = fileprint_eval_capture(expr);
    ASSERT_STR_EQ(out, "epsilon\ndelta\ngamma\nbeta\nalpha\n");
    free(out);
    unlink(path);
}

void test_fileprint_span_negative_step_partial(void) {
    const char* path = fileprint_make_input("fp_span_revpart", FP_FIVE);
    char expr[512];
    snprintf(expr, sizeof(expr), "FilePrint[\"%s\", 5;;1;;-2]", path);
    char* out = fileprint_eval_capture(expr);
    ASSERT_STR_EQ(out, "epsilon\ngamma\nalpha\n");
    free(out);
    unlink(path);
}

/* Empty input file: every selector form is a no-op. */
void test_fileprint_empty_file(void) {
    const char* path = fileprint_make_input("fp_empty", "");
    char expr[512];
    snprintf(expr, sizeof(expr), "FilePrint[\"%s\"]", path);
    char* out = fileprint_eval_capture(expr);
    ASSERT_STR_EQ(out, "");
    free(out);
    snprintf(expr, sizeof(expr), "FilePrint[\"%s\", 3]", path);
    out = fileprint_eval_capture(expr);
    ASSERT_STR_EQ(out, "");
    free(out);
    snprintf(expr, sizeof(expr), "FilePrint[\"%s\", -3]", path);
    out = fileprint_eval_capture(expr);
    ASSERT_STR_EQ(out, "");
    free(out);
    unlink(path);
}

/* Files without a trailing newline get one synthesised so the next
 * REPL prompt isn't appended to the last line.  This applies only when
 * we actually emit the unterminated tail. */
void test_fileprint_no_trailing_newline(void) {
    const char* path = fileprint_make_input("fp_unterm", "one\ntwo\nthree");
    char expr[512];
    snprintf(expr, sizeof(expr), "FilePrint[\"%s\"]", path);
    char* out = fileprint_eval_capture(expr);
    ASSERT_STR_EQ(out, "one\ntwo\nthree\n");
    free(out);
    /* Print only the first two lines — both end with '\n' in the file,
     * so we should NOT inject an extra newline. */
    snprintf(expr, sizeof(expr), "FilePrint[\"%s\", 2]", path);
    out = fileprint_eval_capture(expr);
    ASSERT_STR_EQ(out, "one\ntwo\n");
    free(out);
    unlink(path);
}

/* A non-existent file should print the diagnostic and yield $Failed
 * without leaving stdout corrupted (i.e. fflush/close work cleanly). */
void test_fileprint_missing_file_returns_failed(void) {
    char missing[256];
    scratch_path(missing, sizeof(missing), "fp_missing");
    /* scratch_path already unlinks it. */

    /* The diagnostic line ends up in our capture; we don't pin its
     * exact text but do verify that the evaluator returned the
     * sentinel symbol. */
    char expr[512];
    snprintf(expr, sizeof(expr), "FilePrint[\"%s\"]", missing);
    char* out = fileprint_eval_capture(expr);
    free(out);
    assert_eval_eq(expr, "$Failed", 0);
}

void test_fileprint_unevaluated_on_bad_args(void) {
    /* Wrong arity. */
    assert_eval_eq("FilePrint[]",                "FilePrint[]",                0);
    assert_eval_eq("FilePrint[a, b, c]",         "FilePrint[a, b, c]",         0);
    /* Non-string filename. */
    assert_eval_eq("FilePrint[42]",              "FilePrint[42]",              0);
    /* Bad selector type — neither integer nor Span. */
    assert_eval_eq("FilePrint[\"/tmp/x\", x]",   "FilePrint[\"/tmp/x\", x]",   0);
    /* Symbolic filename flows through (the caller might bind it later). */
    assert_eval_eq("FilePrint[x]",               "FilePrint[x]",               0);
}

void test_fileprint_zero_step_is_unevaluated(void) {
    /* Span with step 0 is mathematically ill-formed; leaving the call
     * unevaluated is safer than infinite-looping.  The printer falls
     * back to FullForm-ish Span[...] notation for the zero-step case,
     * so we pin that exact rendering. */
    const char* path = fileprint_make_input("fp_zero_step", FP_FIVE);
    char expr[512];
    char expected[512];
    snprintf(expr,     sizeof(expr),     "FilePrint[\"%s\", 1;;5;;0]", path);
    snprintf(expected, sizeof(expected), "FilePrint[\"%s\", Span[1, 5, 0]]", path);
    assert_eval_eq(expr, expected, 0);
    unlink(path);
}

void test_fileprint_protected(void) {
    assert_eval_eq("MemberQ[Attributes[FilePrint], Protected]", "True", 0);
}

/* ===== FileNameJoin ===== */

void test_filenamejoin_basic(void) {
    assert_eval_eq("FileNameJoin[{\"dir1\", \"dir2\", \"file\"}]",
                   "\"dir1/dir2/file\"", 0);
}

void test_filenamejoin_prejoined_component(void) {
    /* A component may itself contain separators; it is split and rejoined. */
    assert_eval_eq("FileNameJoin[{\"dir1/dir2\", \"file\"}]",
                   "\"dir1/dir2/file\"", 0);
}

void test_filenamejoin_absolute_from_empty_leading(void) {
    /* An empty leading component produces a leading separator. */
    assert_eval_eq("FileNameJoin[{\"\", \"usr\", \"bin\"}]", "\"/usr/bin\"", 0);
}

void test_filenamejoin_absolute_from_leading_separator(void) {
    assert_eval_eq("FileNameJoin[{\"/usr\", \"bin\"}]", "\"/usr/bin\"", 0);
}

void test_filenamejoin_single_name_canonicalizes(void) {
    /* FileNameJoin["name"] canonicalizes a lone name. */
    assert_eval_eq("FileNameJoin[\"a/b/c\"]", "\"a/b/c\"", 0);
}

void test_filenamejoin_collapses_duplicate_separators(void) {
    assert_eval_eq("FileNameJoin[{\"dir1/\", \"file\"}]", "\"dir1/file\"", 0);
    assert_eval_eq("FileNameJoin[{\"a//b\", \"c\"}]",     "\"a/b/c\"",     0);
}

void test_filenamejoin_empty_list(void) {
    assert_eval_eq("FileNameJoin[{}]", "\"\"", 0);
}

void test_filenamejoin_windows_separator(void) {
    /* Backslash separators when OperatingSystem->"Windows". */
    assert_eval_eq(
        "FileNameJoin[{\"dir1\", \"dir2\", \"file\"}, OperatingSystem->\"Windows\"]",
        "\"dir1\\dir2\\file\"", 0);
}

void test_filenamejoin_windows_unc_share(void) {
    /* A leading \\server\share is preserved verbatim as a UNC prefix. */
    assert_eval_eq(
        "FileNameJoin[{\"\\\\\\\\server\\\\share\", \"path\", \"file\"}, OperatingSystem->\"Windows\"]",
        "\"\\\\server\\share\\path\\file\"", 0);
}

void test_filenamejoin_unix_option(void) {
    /* OperatingSystem->"Unix" forces '/' regardless of host. */
    assert_eval_eq(
        "FileNameJoin[{\"dir1\", \"dir2\"}, OperatingSystem->\"Unix\"]",
        "\"dir1/dir2\"", 0);
}

void test_filenamejoin_macosx_option(void) {
    /* "MacOSX" behaves like "Unix". */
    assert_eval_eq(
        "FileNameJoin[{\"dir1\", \"dir2\"}, OperatingSystem->\"MacOSX\"]",
        "\"dir1/dir2\"", 0);
}

void test_filenamejoin_unevaluated_on_bad_args(void) {
    assert_eval_eq("FileNameJoin[42]",       "FileNameJoin[42]",       0);
    assert_eval_eq("FileNameJoin[{1, 2}]",   "FileNameJoin[{1, 2}]",   0);
    assert_eval_eq("FileNameJoin[x]",        "FileNameJoin[x]",        0);
}

void test_filenamejoin_unknown_os_is_unevaluated(void) {
    assert_eval_eq(
        "FileNameJoin[{\"a\"}, OperatingSystem->\"VMS\"]",
        "FileNameJoin[{\"a\"}, OperatingSystem -> \"VMS\"]", 0);
}

void test_filenamejoin_zero_args_is_unevaluated(void) {
    /* Prints FileNameJoin::argx to stderr and leaves the call unevaluated. */
    assert_eval_eq("FileNameJoin[]", "FileNameJoin[]", 0);
}

void test_filenamejoin_protected(void) {
    assert_eval_eq("MemberQ[Attributes[FileNameJoin], Protected]", "True", 0);
}

/* ===== FileNameSplit ===== */

void test_filenamesplit_basic(void) {
    assert_eval_eq("FileNameSplit[\"a/b/c\"]", "{\"a\", \"b\", \"c\"}", 0);
}

/* The docstring's own example: leading separator yields a leading "",
 * and the trailing separator produces no empty final element. */
void test_filenamesplit_absolute_and_trailing(void) {
    assert_eval_eq("FileNameSplit[\"/home/sb/mathilda/examples/\"]",
                   "{\"\", \"home\", \"sb\", \"mathilda\", \"examples\"}", 0);
}

/* Duplicate and trailing separators collapse away. */
void test_filenamesplit_collapses_separators(void) {
    assert_eval_eq("FileNameSplit[\"a//b/\"]", "{\"a\", \"b\"}", 0);
}

/* A lone name splits to a single part. */
void test_filenamesplit_single_name(void) {
    assert_eval_eq("FileNameSplit[\"file\"]", "{\"file\"}", 0);
}

/* Empty input yields the empty list; "/" yields just the leading "". */
void test_filenamesplit_edge_empty_and_root(void) {
    assert_eval_eq("FileNameSplit[\"\"]", "{}", 0);
    assert_eval_eq("FileNameSplit[\"/\"]", "{\"\"}", 0);
}

/* FileNameJoin is the inverse: joining a split reconstructs the name. */
void test_filenamesplit_join_round_trip(void) {
    assert_eval_eq("FileNameJoin[FileNameSplit[\"/a/b/c\"]]", "\"/a/b/c\"", 0);
    assert_eval_eq("FileNameJoin[FileNameSplit[\"a/b/c\"]]",  "\"a/b/c\"",  0);
}

/* Under Windows a drive letter falls out as an ordinary first part. */
void test_filenamesplit_windows_drive(void) {
    assert_eval_eq(
        "FileNameSplit[\"C:\\\\path\\\\file\", OperatingSystem->\"Windows\"]",
        "{\"C:\", \"path\", \"file\"}", 0);
}

/* Under Windows a leading \\server\share UNC prefix is kept as one part. */
void test_filenamesplit_windows_unc_share(void) {
    assert_eval_eq(
        "FileNameSplit[\"\\\\\\\\server\\\\share\\\\path\\\\file\", OperatingSystem->\"Windows\"]",
        "{\"\\\\server\\share\", \"path\", \"file\"}", 0);
}

/* OperatingSystem->"Unix"/"MacOSX" force '/' regardless of host; a bare
 * backslash is then an ordinary character, not a separator. */
void test_filenamesplit_unix_and_macosx_option(void) {
    assert_eval_eq(
        "FileNameSplit[\"dir1/dir2\", OperatingSystem->\"Unix\"]",
        "{\"dir1\", \"dir2\"}", 0);
    assert_eval_eq(
        "FileNameSplit[\"dir1/dir2\", OperatingSystem->\"MacOSX\"]",
        "{\"dir1\", \"dir2\"}", 0);
}

void test_filenamesplit_unknown_os_is_unevaluated(void) {
    assert_eval_eq(
        "FileNameSplit[\"a\", OperatingSystem->\"VMS\"]",
        "FileNameSplit[\"a\", OperatingSystem -> \"VMS\"]", 0);
}

void test_filenamesplit_unevaluated_on_bad_args(void) {
    assert_eval_eq("FileNameSplit[42]",     "FileNameSplit[42]",     0);
    assert_eval_eq("FileNameSplit[{1, 2}]", "FileNameSplit[{1, 2}]", 0);
    /* Symbolic argument flows through unevaluated (the caller may bind it). */
    assert_eval_eq("FileNameSplit[x]",      "FileNameSplit[x]",      0);
}

void test_filenamesplit_zero_args_is_unevaluated(void) {
    /* Prints FileNameSplit::argx to stderr and leaves the call unevaluated. */
    assert_eval_eq("FileNameSplit[]", "FileNameSplit[]", 0);
}

void test_filenamesplit_protected(void) {
    assert_eval_eq("MemberQ[Attributes[FileNameSplit], Protected]", "True", 0);
}

void test_filenamesplit_options(void) {
    assert_eval_eq("Options[FileNameSplit]",
                   "{OperatingSystem -> \"Unix\"}", 0);
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

    /* FilePrint */
    TEST(test_fileprint_full_file);
    TEST(test_fileprint_first_n);
    TEST(test_fileprint_first_n_clamps);
    TEST(test_fileprint_zero_is_empty);
    TEST(test_fileprint_last_n);
    TEST(test_fileprint_last_n_clamps);
    TEST(test_fileprint_span_basic);
    TEST(test_fileprint_span_negative_end);
    TEST(test_fileprint_span_empty);
    TEST(test_fileprint_span_step);
    TEST(test_fileprint_span_negative_step);
    TEST(test_fileprint_span_negative_step_partial);
    TEST(test_fileprint_empty_file);
    TEST(test_fileprint_no_trailing_newline);
    TEST(test_fileprint_missing_file_returns_failed);
    TEST(test_fileprint_unevaluated_on_bad_args);
    TEST(test_fileprint_zero_step_is_unevaluated);
    TEST(test_fileprint_protected);

    /* FileNameJoin */
    TEST(test_filenamejoin_basic);
    TEST(test_filenamejoin_prejoined_component);
    TEST(test_filenamejoin_absolute_from_empty_leading);
    TEST(test_filenamejoin_absolute_from_leading_separator);
    TEST(test_filenamejoin_single_name_canonicalizes);
    TEST(test_filenamejoin_collapses_duplicate_separators);
    TEST(test_filenamejoin_empty_list);
    TEST(test_filenamejoin_windows_separator);
    TEST(test_filenamejoin_windows_unc_share);
    TEST(test_filenamejoin_unix_option);
    TEST(test_filenamejoin_macosx_option);
    TEST(test_filenamejoin_unevaluated_on_bad_args);
    TEST(test_filenamejoin_unknown_os_is_unevaluated);
    TEST(test_filenamejoin_zero_args_is_unevaluated);
    TEST(test_filenamejoin_protected);

    /* FileNameSplit */
    TEST(test_filenamesplit_basic);
    TEST(test_filenamesplit_absolute_and_trailing);
    TEST(test_filenamesplit_collapses_separators);
    TEST(test_filenamesplit_single_name);
    TEST(test_filenamesplit_edge_empty_and_root);
    TEST(test_filenamesplit_join_round_trip);
    TEST(test_filenamesplit_windows_drive);
    TEST(test_filenamesplit_windows_unc_share);
    TEST(test_filenamesplit_unix_and_macosx_option);
    TEST(test_filenamesplit_unknown_os_is_unevaluated);
    TEST(test_filenamesplit_unevaluated_on_bad_args);
    TEST(test_filenamesplit_zero_args_is_unevaluated);
    TEST(test_filenamesplit_protected);
    TEST(test_filenamesplit_options);

    printf("All tests passed!\n");
    return 0;
}
