/*
 * test_files.c — Unit tests for the filesystem builtins in
 * src/files.c.
 *
 * Covers:
 *   FileExistsQ  — true/false on a freshly created scratch file,
 *                  on directories, on dangling symlinks, plus the
 *                  unevaluated-on-bad-args contract.
 *   FileSize     — exact byte count of a scratch file, zero for an
 *                  empty file, an Integer (not a Quantity) head, the
 *                  $Failed / Quiet[] / Check[] behaviour on a missing
 *                  file, plus the unevaluated-on-bad-args contract.
 *   FileExtension — every Mathematica edge case from the docstring
 *                   (no extension, ends-with-dot, leading-dot,
 *                   nested extensions, directory specifications,
 *                   pure-directory inputs).
 *   FileBaseName  — symmetric coverage with FileExtension, plus the
 *                   "split off last extension only" behaviour.
 *   Environment   — set/unset variables and the bad-argument messages.
 *   DirectoryName, FileNames, $InputFileName — reference cases checked
 *                   against wolframscript, on a scratch directory tree.
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

/* ===== FileSize ===== */

/* Exact byte count of a scratch file whose contents we control. */
void test_filesize_matches_written_bytes(void) {
    char path[256];
    scratch_path(path, sizeof(path), "size");

    FILE* fp = fopen(path, "wb");
    ASSERT(fp != NULL);
    /* Write exactly 11 bytes — no trailing newline from fputs. */
    ASSERT(fwrite("hello world", 1, 11, fp) == 11);
    fclose(fp);

    char input[512];
    snprintf(input, sizeof(input), "FileSize[\"%s\"]", path);
    assert_eval_eq(input, "11", 0);

    unlink(path);
}

/* An empty file has size zero. */
void test_filesize_zero_for_empty_file(void) {
    char path[256];
    scratch_path(path, sizeof(path), "size_empty");

    FILE* fp = fopen(path, "wb");
    ASSERT(fp != NULL);
    fclose(fp);

    char input[512];
    snprintf(input, sizeof(input), "FileSize[\"%s\"]", path);
    assert_eval_eq(input, "0", 0);

    unlink(path);
}

/* The result is a plain Integer, not a Quantity. */
void test_filesize_head_is_integer(void) {
    char path[256];
    scratch_path(path, sizeof(path), "size_head");

    FILE* fp = fopen(path, "wb");
    ASSERT(fp != NULL);
    ASSERT(fwrite("abc", 1, 3, fp) == 3);
    fclose(fp);

    char input[512];
    snprintf(input, sizeof(input), "Head[FileSize[\"%s\"]]", path);
    assert_eval_eq(input, "Integer", 0);

    unlink(path);
}

/* A missing file yields $Failed.  Wrapped in Quiet[] so the
 * FileSize::nffil diagnostic is suppressed — which also verifies the
 * Quiet-aware message helper actually silences it. */
void test_filesize_missing_returns_failed(void) {
    char missing[256];
    scratch_path(missing, sizeof(missing), "size_missing");
    /* scratch_path already unlinks it. */

    char input[512];
    snprintf(input, sizeof(input), "Quiet[FileSize[\"%s\"]]", missing);
    assert_eval_eq(input, "$Failed", 0);
}

/* Check[] must see that FileSize fired a message on the missing file
 * (mth_msg_note_fired wiring) and return the failure expression.  The
 * outer Quiet[] keeps the message off stderr while Check[] still detects
 * it — the idiomatic Quiet[Check[...]] pairing. */
void test_filesize_check_sees_message(void) {
    char missing[256];
    scratch_path(missing, sizeof(missing), "size_check");

    char input[512];
    snprintf(input, sizeof(input),
             "Quiet[Check[FileSize[\"%s\"], \"caught\"]]", missing);
    assert_eval_eq(input, "\"caught\"", 0);
}

/* Bad-args contract: wrong arity or non-string argument leaves the call
 * unevaluated (NULL return, res not freed), and a symbolic argument flows
 * through so callers can build FileSize[x] before x is bound. */
void test_filesize_unevaluated_on_bad_args(void) {
    assert_eval_eq("FileSize[]",     "FileSize[]",     0);
    assert_eval_eq("FileSize[42]",   "FileSize[42]",   0);
    assert_eval_eq("FileSize[a, b]", "FileSize[a, b]", 0);
    assert_eval_eq("FileSize[x]",    "FileSize[x]",    0);
}

void test_filesize_protected(void) {
    assert_eval_eq("MemberQ[Attributes[FileSize], Protected]", "True", 0);
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
        "\"dir1\\\\dir2\\\\file\"", 0);
}

void test_filenamejoin_windows_unc_share(void) {
    /* A leading \\server\share is preserved verbatim as a UNC prefix. */
    assert_eval_eq(
        "FileNameJoin[{\"\\\\\\\\server\\\\share\", \"path\", \"file\"}, OperatingSystem->\"Windows\"]",
        "\"\\\\\\\\server\\\\share\\\\path\\\\file\"", 0);
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
        "{\"\\\\\\\\server\\\\share\", \"path\", \"file\"}", 0);
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

/* ===== DirectoryName ===== */

void test_directoryname_cases(void) {
    /* Reference values from wolframscript. */
    assert_eval_eq("DirectoryName[\"/a/b/c.txt\"]", "\"/a/b/\"", 0);
    assert_eval_eq("DirectoryName[\"/a/b/\"]", "\"/a/\"", 0);
    assert_eval_eq("DirectoryName[\"c.txt\"]", "\"\"", 0);
    assert_eval_eq("DirectoryName[\"/\"]", "\"\"", 0);
    assert_eval_eq("DirectoryName[\"/a\"]", "\"/\"", 0);
    assert_eval_eq("DirectoryName[\"a/b\"]", "\"a/\"", 0);
    assert_eval_eq("DirectoryName[\"a//b\"]", "\"a/\"", 0);
    assert_eval_eq("DirectoryName[\"///a//b//\"]", "\"/a/\"", 0);
    assert_eval_eq("DirectoryName[\"/a/b/c\", 2]", "\"/a/\"", 0);
    assert_eval_eq("DirectoryName[\"/a/b/c\", 5]", "\"\"", 0);
}

void test_directoryname_unevaluated_on_bad_args(void) {
    assert_eval_eq("Quiet[DirectoryName[3]]", "DirectoryName[3]", 0);
    assert_eval_eq("Quiet[DirectoryName[\"/a/b\", 0]]", "DirectoryName[\"/a/b\", 0]", 0);
    assert_eval_eq("Check[DirectoryName[3], \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("Attributes[DirectoryName]", "{Protected}", 0);
}

/* ===== FileNames ===== */

/* Build  root/{a.m, b.wl, c.txt, .hid, sub/{d.m, deep/e.m}}. */
static void fn_tree(char* root, size_t cap) {
    char p[512];
    scratch_path(root, cap, "fntree");
    ASSERT(mkdir(root, 0700) == 0);
    const char* files[] = { "a.m", "b.wl", "c.txt", ".hid" };
    for (size_t i = 0; i < 4; i++) {
        snprintf(p, sizeof p, "%s/%s", root, files[i]);
        FILE* fp = fopen(p, "w"); ASSERT(fp != NULL); fclose(fp);
    }
    snprintf(p, sizeof p, "%s/sub", root);       ASSERT(mkdir(p, 0700) == 0);
    snprintf(p, sizeof p, "%s/sub/deep", root);  ASSERT(mkdir(p, 0700) == 0);
    snprintf(p, sizeof p, "%s/sub/d.m", root);
    FILE* fp = fopen(p, "w"); ASSERT(fp != NULL); fclose(fp);
    snprintf(p, sizeof p, "%s/sub/deep/e.m", root);
    fp = fopen(p, "w"); ASSERT(fp != NULL); fclose(fp);
}

static void fn_tree_rm(const char* root) {
    char p[512];
    const char* rel[] = { "sub/deep/e.m", "sub/d.m", "a.m", "b.wl", "c.txt", ".hid" };
    for (size_t i = 0; i < 6; i++) {
        snprintf(p, sizeof p, "%s/%s", root, rel[i]);
        unlink(p);
    }
    snprintf(p, sizeof p, "%s/sub/deep", root); rmdir(p);
    snprintf(p, sizeof p, "%s/sub", root);      rmdir(p);
    rmdir(root);
}

/* Evaluate `fmt` with every %s replaced by the tree root, expect True. */
static void fn_check(const char* root, const char* fmt) {
    char src[2048], cur[2048];
    snprintf(cur, sizeof cur, "%s", fmt);
    for (;;) {
        char* at = strstr(cur, "%s");
        if (!at) break;
        *at = '\0';
        snprintf(src, sizeof src, "%s%s%s", cur, root, at + 2);
        snprintf(cur, sizeof cur, "%s", src);
    }
    assert_eval_eq(cur, "True", 0);
}

void test_filenames_tree(void) {
    char root[256];
    fn_tree(root, sizeof root);
    /* Depth 1, a dir prefix on every name, hidden files and dirs included. */
    fn_check(root, "FileNames[\"*\", \"%s\"] === Sort[{\"%s/a.m\", \"%s/b.wl\", "
                   "\"%s/c.txt\", \"%s/.hid\", \"%s/sub\"}]");
    fn_check(root, "FileNames[\"*.m\", \"%s\"] === {\"%s/a.m\"}");
    /* A trailing separator on the dir is not doubled. */
    fn_check(root, "FileNames[\"*.m\", \"%s/\"] === {\"%s/a.m\"}");
    /* Depth 2 and Infinity. */
    fn_check(root, "FileNames[\"*.m\", \"%s\", 2] === Sort[{\"%s/a.m\", \"%s/sub/d.m\"}]");
    fn_check(root, "FileNames[\"*.m\", \"%s\", Infinity] === "
                   "Sort[{\"%s/a.m\", \"%s/sub/d.m\", \"%s/sub/deep/e.m\"}]");
    /* Several dirs; the directory part of a form moves onto the dir. */
    fn_check(root, "FileNames[\"*.m\", {\"%s\", \"%s/sub\"}] === "
                   "Sort[{\"%s/a.m\", \"%s/sub/d.m\"}]");
    fn_check(root, "FileNames[\"sub/*.m\", \"%s\"] === {\"%s/sub/d.m\"}");
    /* List forms (deduplicated), StringExpression and RegularExpression. */
    fn_check(root, "FileNames[{\"*.m\", \"a*\"}, \"%s\"] === {\"%s/a.m\"}");
    fn_check(root, "FileNames[__ ~~ \".wl\", \"%s\"] === {\"%s/b.wl\"}");
    fn_check(root, "FileNames[RegularExpression[\"[ab]\\\\..*\"], \"%s\"] === "
                   "Sort[{\"%s/a.m\", \"%s/b.wl\"}]");
    /* Case: explicit IgnoreCase overrides the platform default. */
    fn_check(root, "FileNames[\"*.M\", \"%s\", IgnoreCase -> False] === {}");
    fn_check(root, "FileNames[\"*.M\", \"%s\", IgnoreCase -> True] === {\"%s/a.m\"}");
    /* A missing directory or a plain file contributes nothing. */
    fn_check(root, "FileNames[\"*\", \"%s/nonexistent\"] === {}");
    fn_check(root, "FileNames[\"*\", \"%s/a.m\"] === {}");
    fn_tree_rm(root);
}

void test_filenames_bad_args(void) {
    assert_eval_eq("Quiet[FileNames[\"*\", \".\", -1]]", "FileNames[\"*\", \".\", -1]", 0);
    assert_eval_eq("Check[FileNames[\"*\", \".\", -1], \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("FileNames[\"*\", 3]", "FileNames[\"*\", 3]", 0);
    assert_eval_eq("Attributes[FileNames]", "{Protected}", 0);
    assert_eval_eq("Options[FileNames]", "{IgnoreCase -> Automatic}", 0);
}

/* ===== Directory / SetDirectory / ResetDirectory / DirectoryStack ===== */

void test_directory_stack(void) {
    char root[256], src[1024];
    fn_tree(root, sizeof root);
    assert_eval_eq("d0 = Directory[]; StringQ[d0] && DirectoryStack[] === {}", "True", 0);

    /* SetDirectory returns the canonical new directory and pushes the old. */
    snprintf(src, sizeof src, "d1 = SetDirectory[\"%s\"]; d1 === Directory[]", root);
    assert_eval_eq(src, "True", 0);
    assert_eval_eq("DirectoryStack[] === {d0}", "True", 0);
    /* Relative paths and FileNames now see the new directory. */
    assert_eval_eq("d2 = SetDirectory[\"./sub/../sub/\"]; "
                   "{d2 === d1 <> \"/sub\", DirectoryStack[] === {d1, d0}}",
                   "{True, True}", 0);
    assert_eval_eq("FileNames[\"*.m\"]", "{\"d.m\"}", 0);
    /* Pop back, then the empty-stack message. */
    assert_eval_eq("{ResetDirectory[] === d1, ResetDirectory[] === d0, Directory[] === d0}",
                   "{True, True, True}", 0);
    assert_eval_eq("Quiet[ResetDirectory[]] === d0", "True", 0);
    assert_eval_eq("Check[ResetDirectory[], \"msg\"]", "\"msg\"", 0);

    /* Failures leave the directory and the stack unchanged. */
    snprintf(src, sizeof src, "Quiet[SetDirectory[\"%s/nonexistent\"]]", root);
    assert_eval_eq(src, "$Failed", 0);
    snprintf(src, sizeof src, "Quiet[SetDirectory[\"%s/a.m\"]]", root);
    assert_eval_eq(src, "$Failed", 0);
    assert_eval_eq("Quiet[SetDirectory[3]]", "SetDirectory[3]", 0);
    assert_eval_eq("{Directory[] === d0, DirectoryStack[]}", "{True, {}}", 0);
    assert_eval_eq("Quiet[Directory[1]]", "Directory[1]", 0);

    /* SetDirectory[] goes home. */
    assert_eval_eq("StringQ[$HomeDirectory] && SetDirectory[] === Directory[] && "
                   "ResetDirectory[] === d0", "True", 0);
    assert_eval_eq("{Attributes[SetDirectory], Attributes[ResetDirectory], "
                   "Attributes[Directory], Attributes[DirectoryStack]}",
                   "{{Protected}, {Protected}, {Protected}, {Protected}}", 0);
    /* $HomeDirectory / $InstallationDirectory: Protected strings; the
     * installation root is absolute and actually holds the bootstrap tree
     * (from tests/build it is reached through the relative CWD ladder). */
    assert_eval_eq("{Attributes[$HomeDirectory], Attributes[$InstallationDirectory]}",
                   "{{Protected}, {Protected}}", 0);
    assert_eval_eq("StringQ[$InstallationDirectory] && "
                   "(FileExistsQ[FileNameJoin[{$InstallationDirectory, \"src\", \"internal\", \"init.m\"}]] || "
                   "FileExistsQ[FileNameJoin[{$InstallationDirectory, \"share\", \"mathilda\", \"internal\", \"init.m\"}]])",
                   "True", 0);
#ifndef _WIN32
    assert_eval_eq("StringTake[$InstallationDirectory, 1]", "\"/\"", 0);
#endif
    fn_tree_rm(root);
}

/* ===== Environment ===== */

void test_environment(void) {
    /* A set variable comes back as its string value (PATH is set in every
     * environment the suite runs in); an unset one is $Failed. */
    assert_eval_eq("StringQ[Environment[\"PATH\"]] && StringLength[Environment[\"PATH\"]] > 0",
                   "True", 0);
    assert_eval_eq("Environment[\"MATHILDA_SURELY_UNSET_VARIABLE_42\"]", "$Failed", 0);
    /* Non-string argument: Environment::string, unevaluated; wrong arity: argx. */
    assert_eval_eq("Quiet[Environment[3]]", "Environment[3]", 0);
    assert_eval_eq("Check[Environment[3], \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("Check[Environment[], \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("Quiet[Environment[\"PATH\", 1]]", "Environment[\"PATH\", 1]", 0);
    assert_eval_eq("Attributes[Environment]", "{Protected}", 0);
}

/* ===== ParentDirectory / FileNameTake / AbsoluteFileName / FileType /
 *       ToFileName / DeleteFile ===== */

void test_parentdirectory(void) {
    char root[256], src[1200];
    fn_tree(root, sizeof root);
    /* Canonical absolute form of the scratch root (on macOS /tmp is a symlink). */
    snprintf(src, sizeof src, "pdRoot = AbsoluteFileName[\"%s\"]; StringQ[pdRoot]", root);
    assert_eval_eq(src, "True", 0);
    snprintf(src, sizeof src, "ParentDirectory[\"%s/sub\"] === pdRoot", root);
    assert_eval_eq(src, "True", 0);
    snprintf(src, sizeof src, "ParentDirectory[\"%s/sub/\"] === pdRoot", root);
    assert_eval_eq(src, "True", 0);
    snprintf(src, sizeof src, "ParentDirectory[\"%s/sub\", 2] === ParentDirectory[pdRoot]", root);
    assert_eval_eq(src, "True", 0);
    assert_eval_eq("{ParentDirectory[\"/\"], ParentDirectory[\"/tmp\", 50]}", "{\"/\", \"/\"}", 0);
    assert_eval_eq("ParentDirectory[] === ParentDirectory[Directory[]]", "True", 0);
    /* A relative directory resolves against the working directory. */
    snprintf(src, sizeof src, "SetDirectory[\"%s\"]; pd = ParentDirectory[\"sub\"]; ResetDirectory[]; pd === pdRoot", root);
    assert_eval_eq(src, "True", 0);
    /* Not an existing directory: dirnex, unevaluated. */
    snprintf(src, sizeof src, "Quiet[ParentDirectory[\"%s/a.m\"]] === ParentDirectory[\"%s/a.m\"]", root, root);
    assert_eval_eq(src, "True", 0);
    assert_eval_eq("Check[ParentDirectory[\"/no/such/dir\"], \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("Quiet[ParentDirectory[3]]", "ParentDirectory[3]", 0);
    assert_eval_eq("Quiet[ParentDirectory[\"/\", 0]]", "ParentDirectory[\"/\", 0]", 0);
    assert_eval_eq("Attributes[ParentDirectory]", "{Protected}", 0);
    fn_tree_rm(root);
}

void test_filenametake(void) {
    assert_eval_eq("FileNameTake[\"/a/b/c.txt\"]", "\"c.txt\"", 0);
    assert_eval_eq("{FileNameTake[\"/a/b/c.txt\", 1], FileNameTake[\"/a/b/c.txt\", 2]}", "{\"/\", \"/a\"}", 0);
    assert_eval_eq("FileNameTake[\"/a/b/c.txt\", -2]", "\"b/c.txt\"", 0);
    assert_eval_eq("{FileNameTake[\"/a/b/c.txt\", {2}], FileNameTake[\"/a/b/c.txt\", {-1}], "
                   "FileNameTake[\"/a/b/c.txt\", {1}]}", "{\"a\", \"c.txt\", \"/\"}", 0);
    assert_eval_eq("FileNameTake[\"/a/b/c.txt\", {2, 3}]", "\"a/b\"", 0);
    assert_eval_eq("{FileNameTake[\"/a/b/c.txt\", 10], FileNameTake[\"/a/b/c.txt\", -10]}",
                   "{\"/a/b/c.txt\", \"/a/b/c.txt\"}", 0);
    assert_eval_eq("{FileNameTake[\"/a/b/c.txt\", 0], FileNameTake[\"/a/b/c.txt\", {5}], "
                   "FileNameTake[\"\", 1]}", "{\"\", \"\", \"\"}", 0);
    assert_eval_eq("{FileNameTake[\"a/b/\"], FileNameTake[\"a/b/c.txt\", 1]}", "{\"b\", \"a\"}", 0);
    assert_eval_eq("Quiet[FileNameTake[\"/a\", All]]", "FileNameTake[\"/a\", All]", 0);
    assert_eval_eq("Check[FileNameTake[\"/a\", x], \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("Attributes[FileNameTake]", "{Protected}", 0);
}

void test_absolutefilename_filetype(void) {
    char root[256], src[1200];
    fn_tree(root, sizeof root);
    snprintf(src, sizeof src, "afn = AbsoluteFileName[\"%s/sub/../a.m\"]; "
             "{StringMatchQ[afn, \"/*/a.m\"], StringFreeQ[afn, \"..\"]}", root);
    assert_eval_eq(src, "{True, True}", 0);
    snprintf(src, sizeof src, "Quiet[AbsoluteFileName[\"%s/nope\"]]", root);
    assert_eval_eq(src, "$Failed", 0);
    snprintf(src, sizeof src, "Check[AbsoluteFileName[\"%s/nope\"], \"msg\"]", root);
    assert_eval_eq(src, "\"msg\"", 0);
    snprintf(src, sizeof src, "{FileType[\"%s/a.m\"], FileType[\"%s/sub\"], FileType[\"%s/nope\"]}",
             root, root, root);
    assert_eval_eq(src, "{File, Directory, None}", 0);
    assert_eval_eq("Quiet[FileType[3]]", "FileType[3]", 0);
    assert_eval_eq("Check[FileType[3], \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("{Attributes[AbsoluteFileName], Attributes[FileType]}", "{{Protected}, {Protected}}", 0);
    fn_tree_rm(root);
}

void test_tofilename(void) {
    assert_eval_eq("{ToFileName[{\"a\", \"b\"}, \"c.m\"], ToFileName[\"a\", \"c.m\"], ToFileName[{\"a\", \"b\"}]}",
                   "{\"a/b/c.m\", \"a/c.m\", \"a/b/\"}", 0);
    assert_eval_eq("{ToFileName[\"\", \"c.m\"], ToFileName[{}, \"c.m\"]}", "{\"c.m\", \"c.m\"}", 0);
    assert_eval_eq("{ToFileName[{\"/a\", \"b\"}, \"c\"], ToFileName[{\"a/\", \"b\"}, \"c\"], ToFileName[{\"/\", \"a\"}, \"c\"]}",
                   "{\"/a/b/c\", \"a/b/c\", \"/a/c\"}", 0);
    assert_eval_eq("ToFileName[3, \"c\"]", "ToFileName[3, \"c\"]", 0);
    assert_eval_eq("Attributes[ToFileName]", "{Protected}", 0);
}

void test_deletefile(void) {
    char root[256], src[1200];
    fn_tree(root, sizeof root);
    snprintf(src, sizeof src, "Put[1, \"%s/del.txt\"]; {DeleteFile[\"%s/del.txt\"], FileExistsQ[\"%s/del.txt\"]}",
             root, root, root);
    assert_eval_eq(src, "{Null, False}", 0);
    /* Missing file and directory: a message and $Failed. */
    snprintf(src, sizeof src, "Quiet[{DeleteFile[\"%s/del.txt\"], DeleteFile[\"%s/sub\"]}]", root, root);
    assert_eval_eq(src, "{$Failed, $Failed}", 0);
    snprintf(src, sizeof src, "Check[DeleteFile[\"%s/del.txt\"], \"msg\"]", root);
    assert_eval_eq(src, "\"msg\"", 0);
    snprintf(src, sizeof src, "FileType[\"%s/sub\"]", root);
    assert_eval_eq(src, "Directory", 0);
    /* List form. */
    snprintf(src, sizeof src, "Put[1, \"%s/d1\"]; Put[2, \"%s/d2\"]; "
             "{DeleteFile[{\"%s/d1\", \"%s/d2\"}], FileExistsQ[\"%s/d1\"], FileExistsQ[\"%s/d2\"], DeleteFile[{}]}",
             root, root, root, root, root, root);
    assert_eval_eq(src, "{Null, False, False, Null}", 0);
    assert_eval_eq("Quiet[DeleteFile[3]]", "DeleteFile[3]", 0);
    assert_eval_eq("Attributes[DeleteFile]", "{Protected}", 0);
    fn_tree_rm(root);
}

void test_syntaxq(void) {
    assert_eval_eq("{SyntaxQ[\"1+2\"], SyntaxQ[\"a=1\\nb=2\"], SyntaxQ[\"a;;b\"], SyntaxQ[\"(*c*) x\"], "
                   "SyntaxQ[\"(*c*)\"], SyntaxQ[\" \"], SyntaxQ[\"a;\"]}",
                   "{True, True, True, True, True, True, True}", 0);
    assert_eval_eq("{SyntaxQ[\"1+\"], SyntaxQ[\"f[x\"], SyntaxQ[\"x)\"], SyntaxQ[\"\"], SyntaxQ[\"\\\"abc\"]}",
                   "{False, False, False, False, False}", 0);
    /* Nothing is evaluated. */
    assert_eval_eq("synQ = 0; {SyntaxQ[\"synQ = 5\"], synQ}", "{True, 0}", 0);
    assert_eval_eq("Quiet[SyntaxQ[3]]", "SyntaxQ[3]", 0);
    assert_eval_eq("Check[SyntaxQ[3], \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("Attributes[SyntaxQ]", "{Protected}", 0);
}

/* ===== Import text formats / ReadString ===== */

/* Write `bytes` (exact, binary) to a scratch file named by tag. */
static void write_bytes(char* path, size_t cap, const char* tag, const char* bytes) {
    scratch_path(path, cap, tag);
    FILE* fp = fopen(path, "wb");
    ASSERT(fp != NULL);
    fputs(bytes, fp);
    fclose(fp);
}

void test_import_text(void) {
    char t[256], u[256], e[256], v[256], x[256], txt[300], src[2048];
    write_bytes(t, sizeof t, "it_t", "line1\nline2\n\nline4");
    write_bytes(u, sizeof u, "it_u", "abc\n");
    write_bytes(e, sizeof e, "it_e", "");
    write_bytes(v, sizeof v, "it_v", "abc\n\n\n");
    write_bytes(x, sizeof x, "it_x", "a\r\nb\r\n");
    /* "Text" drops ONE trailing line break; "String" is the bytes exactly;
     * "Lines" drops trailing empty lines.  Text and Lines read CRLF and a lone
     * CR as one line break; String and ReadString keep the CR bytes. */
    snprintf(src, sizeof src, "{Import[\"%s\", \"Text\"], Import[\"%s\", \"Text\"], Import[\"%s\", \"Text\"], "
             "Import[\"%s\", \"Text\"], Import[\"%s\", \"Text\"]}", t, u, e, v, x);
    assert_eval_eq(src, "{\"line1\\nline2\\n\\nline4\", \"abc\", \"\", \"abc\\n\\n\", \"a\\nb\"}", 0);
    snprintf(src, sizeof src, "{Import[\"%s\", \"String\"], Import[\"%s\", \"String\"], Import[\"%s\", {\"Text\"}]}",
             u, x, u);
    assert_eval_eq(src, "{\"abc\\n\", \"a\\r\\nb\\r\\n\", \"abc\"}", 0);
    snprintf(src, sizeof src, "{Import[\"%s\", \"Lines\"], Import[\"%s\", \"Lines\"], Import[\"%s\", \"Lines\"], "
             "Import[\"%s\", \"Lines\"]}", t, e, v, x);
    assert_eval_eq(src, "{{\"line1\", \"line2\", \"\", \"line4\"}, {}, {\"abc\"}, {\"a\", \"b\"}}", 0);
    /* A .txt file imports as Text by default. */
    snprintf(txt, sizeof txt, "%s.txt", u);
    rename(u, txt);
    snprintf(src, sizeof src, "Import[\"%s\"]", txt);
    assert_eval_eq(src, "\"abc\"", 0);
    /* Missing file: Import::nffil, $Failed. */
    snprintf(src, sizeof src, "Quiet[Import[\"%s.nope\", \"Text\"]]", t);
    assert_eval_eq(src, "$Failed", 0);
    snprintf(src, sizeof src, "Check[Import[\"%s.nope\", \"Text\"], \"msg\"]", t);
    assert_eval_eq(src, "\"msg\"", 0);

    /* ReadString: whole file, terminator, EndOfFile, streams. */
    snprintf(src, sizeof src, "{ReadString[\"%s\"], ReadString[\"%s\"], ReadString[\"%s\"], "
             "ReadString[\"%s\", \"\\n\"], ReadString[\"%s\", {\"2\", \"1\"}]}", txt, x, e, t, t);
    assert_eval_eq(src, "{\"abc\\n\", \"a\\r\\nb\\r\\n\", EndOfFile, \"line1\", \"line\"}", 0);
    /* A lone CR is a line break for Text, Lines and String records. */
    {
        char y[256];
        write_bytes(y, sizeof y, "it_y", "x\ry\r");
        snprintf(src, sizeof src, "{Import[\"%s\", \"Text\"], Import[\"%s\", \"Lines\"], "
                 "ToCharacterCode[Import[\"%s\", \"String\"]], ReadList[\"%s\", String]}", y, y, y, y);
        assert_eval_eq(src, "{\"x\\ny\", {\"x\", \"y\"}, {120, 13, 121, 13}, {\"x\", \"y\"}}", 0);
        unlink(y);
        /* Leading empty lines are dropped as trailing ones are; interior stay. */
        write_bytes(y, sizeof y, "it_y", "\n\nab\n\ncd\n\n");
        snprintf(src, sizeof src, "Import[\"%s\", \"Lines\"]", y);
        assert_eval_eq(src, "{\"ab\", \"\", \"cd\"}", 0);
        unlink(y);
    }
    /* On a stream the terminator is left unread, and the stream advances. */
    snprintf(src, sizeof src, "rs = OpenRead[\"%s\"]; r = {ReadString[rs, \"\\n\"], ReadString[rs], ReadString[rs]}; "
             "Close[rs]; r", t);
    assert_eval_eq(src, "{\"line1\", \"\\nline2\\n\\nline4\", EndOfFile}", 0);
    snprintf(src, sizeof src, "Quiet[ReadString[\"%s.nope\"]]", t);
    assert_eval_eq(src, "$Failed", 0);
    snprintf(src, sizeof src, "Check[ReadString[\"%s.nope\"], \"msg\"]", t);
    assert_eval_eq(src, "\"msg\"", 0);
    assert_eval_eq("ReadString[3]", "ReadString[3]", 0);
    assert_eval_eq("Attributes[ReadString]", "{Protected}", 0);
    unlink(t); unlink(txt); unlink(e); unlink(v); unlink(x);
}

/* ===== $InputFileName ===== */

void test_inputfilename_get(void) {
    char path[256], src[512];
    scratch_path(path, sizeof(path), "ifn.m");
    FILE* fp = fopen(path, "w");
    ASSERT(fp != NULL);
    fputs("ifnSeen = $InputFileName;\n", fp);
    fclose(fp);

    assert_eval_eq("$InputFileName", "\"\"", 0);
    snprintf(src, sizeof src, "Get[\"%s\"]; StringMatchQ[ifnSeen, \"/*%s\"]",
             path, strrchr(path, '/'));
    assert_eval_eq(src, "True", 0);
    /* Restored once the file has been read. */
    assert_eval_eq("$InputFileName", "\"\"", 0);
    unlink(path);
}

/* $TemporaryDirectory: a Protected System` string with no trailing separator,
 * naming an existing directory. */
static void test_temporary_directory(void) {
    assert_eval_eq("StringQ[$TemporaryDirectory] && StringLength[$TemporaryDirectory] > 0", "True", 0);
    assert_eval_eq("FileType[$TemporaryDirectory]", "Directory", 0);
    assert_eval_eq("StringTake[$TemporaryDirectory, -1] =!= \"/\" || $TemporaryDirectory === \"/\"", "True", 0);
    assert_eval_eq("{Attributes[$TemporaryDirectory], Context[$TemporaryDirectory]}",
                   "{{Protected}, \"System`\"}", 0);
    assert_eval_eq("tdp = FileNameJoin[{$TemporaryDirectory, \"mth_td_probe.txt\"}]; "
                   "WriteString[OpenWrite[tdp], \"x\"]; Close[tdp]; "
                   "{FileExistsQ[tdp], DeleteFile[tdp], FileExistsQ[tdp]}", "{True, Null, False}", 0);
}

/* $BaseDirectory / $UserBaseDirectory: Protected System` strings with no
 * trailing separator; the user base lives under $HomeDirectory (unless
 * MATHILDA_USERBASE overrides it). */
static void test_base_directories(void) {
    assert_eval_eq("{Attributes[$BaseDirectory], Attributes[$UserBaseDirectory], "
                   "Context[$BaseDirectory], Context[$UserBaseDirectory]}",
                   "{{Protected}, {Protected}, \"System`\", \"System`\"}", 0);
    assert_eval_eq("StringQ[$BaseDirectory] && StringQ[$UserBaseDirectory]", "True", 0);
    assert_eval_eq("StringEndsQ[$BaseDirectory, \"Mathilda\"]", "True", 0);
    if (!getenv("MATHILDA_USERBASE"))
        assert_eval_eq("StringStartsQ[$UserBaseDirectory, $HomeDirectory] && "
                       "StringEndsQ[$UserBaseDirectory, \"Mathilda\"]", "True", 0);
#if defined(__APPLE__)
    if (!getenv("MATHILDA_BASE"))
        assert_eval_eq("$BaseDirectory", "\"/Library/Mathilda\"", 0);
#endif
}

/* DirectoryQ: True only for an existing directory; String or File[...],
 * `~` expanded; any other argument (a list included) is False; a wrong
 * argument count is DirectoryQ::argx and unevaluated.  FileType and
 * FileExistsQ take File[...] and `~` the same way. */
static void test_directoryq(void) {
    assert_eval_eq("DirectoryQ[$TemporaryDirectory]", "True", 0);
    assert_eval_eq("DirectoryQ[File[$TemporaryDirectory]]", "True", 0);
    assert_eval_eq("DirectoryQ[\"~\"]", "True", 0);
    assert_eval_eq("DirectoryQ[\"~/\"]", "True", 0);
    assert_eval_eq("{DirectoryQ[\".\"], DirectoryQ[\"..\"], DirectoryQ[\"/\"]}", "{True, True, True}", 0);
    assert_eval_eq("DirectoryQ[\"/no/such/dir/mth\"]", "False", 0);
    assert_eval_eq("DirectoryQ[\"\"]", "False", 0);
    assert_eval_eq("{DirectoryQ[3], DirectoryQ[x], DirectoryQ[{\"/tmp\"}]}", "{False, False, False}", 0);
    assert_eval_eq("dqf = FileNameJoin[{$TemporaryDirectory, \"mth_dq_probe.txt\"}]; "
                   "WriteString[OpenWrite[dqf], \"x\"]; Close[dqf]; r = DirectoryQ[dqf]; DeleteFile[dqf]; r",
                   "False", 0);
    assert_eval_eq("Quiet[DirectoryQ[]]", "DirectoryQ[]", 0);
    assert_eval_eq("Check[DirectoryQ[\"/\", 2], \"msg\"]", "\"msg\"", 0);
    assert_eval_eq("Attributes[DirectoryQ]", "{Protected}", 0);
    assert_eval_eq("{FileType[\"~\"], FileType[File[$TemporaryDirectory]], FileExistsQ[\"~\"], "
                   "FileExistsQ[File[$TemporaryDirectory]]}", "{Directory, Directory, True, True}", 0);
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

    /* FileSize */
    TEST(test_filesize_matches_written_bytes);
    TEST(test_filesize_zero_for_empty_file);
    TEST(test_filesize_head_is_integer);
    TEST(test_filesize_missing_returns_failed);
    TEST(test_filesize_check_sees_message);
    TEST(test_filesize_unevaluated_on_bad_args);
    TEST(test_filesize_protected);

    /* FileExtension */
    TEST(test_fileextension_simple);
    TEST(test_fileextension_only_last);
    TEST(test_fileextension_none);
    TEST(test_fileextension_trailing_dot);
    TEST(test_fileextension_leading_dot_only);
    TEST(test_fileextension_ignores_directory);
    TEST(test_fileextension_directory_name);
    TEST(test_fileextension_unevaluated_on_bad_args);

    /* DirectoryName / FileNames / $InputFileName */
    TEST(test_directoryname_cases);
    TEST(test_directoryname_unevaluated_on_bad_args);
    TEST(test_filenames_tree);
    TEST(test_filenames_bad_args);
    TEST(test_inputfilename_get);
    TEST(test_directory_stack);
    TEST(test_environment);
    TEST(test_parentdirectory);
    TEST(test_filenametake);
    TEST(test_absolutefilename_filetype);
    TEST(test_tofilename);
    TEST(test_deletefile);
    TEST(test_syntaxq);
    TEST(test_import_text);
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
    TEST(test_temporary_directory);
    TEST(test_base_directories);
    TEST(test_directoryq);

    printf("All tests passed!\n");
    return 0;
}
