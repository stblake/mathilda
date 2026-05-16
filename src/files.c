/*
 * files.c — Filesystem-predicate and path-string builtins.
 *
 * Builtins:
 *   FileExistsQ["name"]   — touches the filesystem (lstat).
 *   FileExtension["name"] — pure string manipulation.
 *   FileBaseName["name"]  — pure string manipulation.
 *
 * Both string builtins follow the POSIX convention: the path
 * separator is `/`, the directory portion is everything up to and
 * including the final `/`, and the file-name component is everything
 * after.  The "extension" is the suffix following the LAST `.` in
 * the file-name component, but a leading `.` (e.g. ".bashrc") is
 * treated as part of the base name rather than an extension marker.
 *
 * No filesystem access happens for FileExtension / FileBaseName; they
 * operate purely on the given string.
 */

/* lstat is POSIX, not C99.  glibc hides it under -std=c99 unless we
 * declare the feature-test macro before any include.  Darwin exposes
 * it implicitly but accepts the macro. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "files.h"
#include "expr.h"
#include "symtab.h"
#include "attr.h"

#include <sys/stat.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

/* FileExistsQ["name"]
 *
 * Returns True if anything exists at the given path, False otherwise.
 * Uses lstat() rather than stat() so dangling symlinks (which are
 * themselves filesystem objects) are still reported as existing.
 *
 * Anything other than a single EXPR_STRING argument leaves the call
 * unevaluated (NULL return).  The evaluator owns `res`: on a
 * successful rewrite it frees the input itself after this builtin
 * returns, so we must NOT free it here. */
Expr* builtin_fileexistsq(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1) return NULL;

    Expr* arg = res->data.function.args[0];
    if (arg->type != EXPR_STRING) return NULL;

    struct stat st;
    int exists = (lstat(arg->data.string, &st) == 0);
    return expr_new_symbol(exists ? "True" : "False");
}

/* Split "path" into (directory-prefix-length, filename-component).
 * The filename component is everything after the last '/'; it may be
 * empty (when "path" is empty or ends with '/').  Pure string view —
 * no allocation. */
static const char* filename_component(const char* path) {
    const char* slash = strrchr(path, '/');
    return slash ? slash + 1 : path;
}

/* Find the extension start within `name` (the filename component
 * only — pass the output of filename_component, not a raw path).
 *
 * Returns the offset of the character AFTER the last '.' that
 * qualifies as an extension separator, or strlen(name) if there is
 * no extension.  By the Mathematica contract, a '.' qualifies only
 * when it is neither the first nor the last character of `name`:
 *
 *   "file.txt"  -> 5  ("txt")
 *   "a.b.c"     -> 4  ("c")
 *   "file"      -> 4  ("")
 *   "file."     -> 5  ("")  trailing dot is not an extension marker
 *   ".bashrc"   -> 7  ("")  leading dot is not an extension marker
 *   ""          -> 0  ("")
 */
static size_t extension_offset(const char* name) {
    size_t len = strlen(name);
    if (len == 0) return 0;
    /* Walk from the end looking for a qualifying '.'.  Stop before
     * index 0 so a lone leading dot is never picked up. */
    for (size_t i = len; i > 1; i--) {
        if (name[i - 1] == '.') {
            /* A '.' at the very end of `name` means no extension. */
            if (i == len) return len;
            return i;
        }
    }
    return len;
}

/* FileExtension["name"]
 *
 * Returns the substring after the last '.' in the filename component
 * of "name", excluding the dot.  Returns "" when no extension is
 * present.  Pure string manipulation — does not touch the
 * filesystem. */
Expr* builtin_fileextension(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1) return NULL;

    Expr* arg = res->data.function.args[0];
    if (arg->type != EXPR_STRING) return NULL;

    const char* name = filename_component(arg->data.string);
    size_t off = extension_offset(name);
    /* extension_offset returns len(name) when there is none; passing
     * name + len yields an empty string, which is the right answer. */
    return expr_new_string(name + off);
}

/* FileBaseName["name"]
 *
 * Returns the filename component of "name" with its trailing
 * extension removed.  When there is no extension the filename
 * component is returned verbatim.  Directory specifications are
 * always dropped.  Pure string manipulation — does not touch the
 * filesystem. */
Expr* builtin_filebasename(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1) return NULL;

    Expr* arg = res->data.function.args[0];
    if (arg->type != EXPR_STRING) return NULL;

    const char* name = filename_component(arg->data.string);
    size_t off = extension_offset(name);
    size_t base_len;
    if (off == strlen(name)) {
        /* No extension — keep the whole filename component, including
         * any trailing '.' (matches "file." -> "file."). */
        base_len = off;
    } else {
        /* off points just past the '.'; exclude the '.' from the
         * base name. */
        base_len = off - 1;
    }

    /* Build a NUL-terminated copy of the prefix so expr_new_string
     * (which strdup's its input) sees exactly the bytes we want.
     * Reusing arg->data.string is not safe — the evaluator will free
     * res (and with it arg) after we return. */
    char* buf = (char*)malloc(base_len + 1);
    if (!buf) return NULL;  /* leave res untouched on OOM */
    memcpy(buf, name, base_len);
    buf[base_len] = '\0';
    Expr* result = expr_new_string(buf);
    free(buf);
    return result;
}

void files_init(void) {
    /* Docstrings live in info.c alongside the other File-I/O entries.
     * All three builtins are Protected with no special evaluation
     * behavior — argument validation handles non-string inputs by
     * leaving the call unevaluated. */
    symtab_add_builtin("FileExistsQ", builtin_fileexistsq);
    symtab_get_def("FileExistsQ")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("FileExtension", builtin_fileextension);
    symtab_get_def("FileExtension")->attributes |= ATTR_PROTECTED;

    symtab_add_builtin("FileBaseName", builtin_filebasename);
    symtab_get_def("FileBaseName")->attributes |= ATTR_PROTECTED;
}
