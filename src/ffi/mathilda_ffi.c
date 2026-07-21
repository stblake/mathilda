/* mathilda_ffi.c — in-process C API for embedding the Mathilda kernel.
 *
 * See mathilda_ffi.h for the contract. This mirrors the parse -> evaluate ->
 * format pipeline that the sidecar's pipe mode (repl.c: pipe_process_input)
 * runs, but returns the formatted text to the caller instead of writing NDJSON
 * to stdout. No readline, no stdio loop, no process — suitable for iOS/Android
 * where the kernel is linked directly into the host app.
 */
/* setenv() is POSIX, hidden by glibc under -std=c99; request it (matches the
 * pattern in repl.c). Must precede any system header include. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "ffi/mathilda_ffi.h"

#include <stdlib.h>
#include <string.h>

#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "print.h"
#include "print_latex.h"
#include "symtab.h"
#include "core.h"
#include "loadmodule.h"
#include "version.h"

/* One-time init guard. The kernel's global symbol table is process-wide, so a
 * second symtab_init()/core_init() would leak or corrupt it. */
static int g_initialized = 0;

/* Optional home directory for the internal module tree, set before init. */
static char g_home[4096] = {0};

/* Return a heap copy of a C string literal so the caller can always free the
 * result of an eval with mathilda_ffi_free() regardless of the path taken. */
static char* ffi_strdup(const char* s) {
    size_t n = strlen(s) + 1;
    char* out = (char*)malloc(n);
    if (out) memcpy(out, s, n);
    return out;
}

void mathilda_ffi_set_home(const char* dir) {
    if (!dir || !*dir) return;
    /* Prefer the process environment so mathilda_resolve_internal() (which
     * reads $MATHILDA_HOME) picks it up as its first candidate. setenv is
     * available on iOS/Android/macOS/Linux; keep a copy for diagnostics. */
    strncpy(g_home, dir, sizeof(g_home) - 1);
    g_home[sizeof(g_home) - 1] = '\0';
    setenv("MATHILDA_HOME", g_home, 1);
}

void mathilda_ffi_init(void) {
    if (g_initialized) return;
    symtab_init();
    core_init();
    /* Load the internal bootstrap (init.m). If MATHILDA_HOME was set via
     * mathilda_ffi_set_home() it is the first search candidate; otherwise the
     * loader falls back to the exe-relative / installed / CWD ladder. A load
     * failure leaves a usable-but-reduced kernel (C builtins only). */
    mathilda_load_module("init.m");
    g_initialized = 1;
}

/* Shared front half: ensure init, parse, evaluate. Returns the evaluated
 * expression (caller frees) or NULL on parse/eval failure, with *parse_failed
 * distinguishing the two so callers can format an appropriate message. */
static Expr* ffi_parse_eval(const char* input, int* parse_failed) {
    *parse_failed = 0;
    if (!g_initialized) mathilda_ffi_init();
    if (!input) { *parse_failed = 1; return NULL; }

    Expr* parsed = parse_expression(input);
    if (!parsed) { *parse_failed = 1; return NULL; }

    Expr* evaluated = evaluate(parsed);
    expr_free(parsed);
    return evaluated; /* may be NULL: "evaluated to nothing" */
}

char* mathilda_ffi_eval(const char* input) {
    int parse_failed = 0;
    Expr* evaluated = ffi_parse_eval(input, &parse_failed);
    if (parse_failed) return ffi_strdup("$Failed (parse error)");
    if (!evaluated)   return ffi_strdup("");

    char* result = expr_to_string(evaluated);
    expr_free(evaluated);
    if (!result) return ffi_strdup("");
    return result; /* expr_to_string already returns caller-owned malloc'd mem */
}

char* mathilda_ffi_eval_latex(const char* input) {
    int parse_failed = 0;
    Expr* evaluated = ffi_parse_eval(input, &parse_failed);
    if (parse_failed) return ffi_strdup("\\text{\\$Failed (parse error)}");
    if (!evaluated)   return ffi_strdup("");

    char* latex = expr_to_latex(evaluated);
    expr_free(evaluated);
    if (!latex) return ffi_strdup("");
    return latex;
}

void mathilda_ffi_free(char* s) {
    free(s);
}

const char* mathilda_ffi_version(void) {
    return MATHILDA_VERSION_STRING;
}
