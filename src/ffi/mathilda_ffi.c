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
#include <stdio.h>

#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "print.h"
#include "print_latex.h"
#include "symtab.h"
#include "core.h"
#include "loadmodule.h"
#include "version.h"
#include "sym_names.h"
#include "graphics_json.h"

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

/* Minimal JSON string escaper — a self-contained copy of repl.c's json_escape
 * (that one is static and lives in repl.o, which is excluded from libmathilda.a).
 * Writes at most outlen-1 bytes plus a NUL. */
static void ffi_json_escape(const char* s, char* out, size_t outlen) {
    size_t i = 0;
    while (*s && i + 7 < outlen) {
        unsigned char c = (unsigned char)*s;
        if (c == '"')       { out[i++] = '\\'; out[i++] = '"'; }
        else if (c == '\\') { out[i++] = '\\'; out[i++] = '\\'; }
        else if (c == '\n') { out[i++] = '\\'; out[i++] = 'n'; }
        else if (c == '\r') { out[i++] = '\\'; out[i++] = 'r'; }
        else if (c == '\t') { out[i++] = '\\'; out[i++] = 't'; }
        else if (c < 0x20)  { i += (size_t)snprintf(out + i, outlen - i, "\\u%04x", (unsigned)c); }
        else                { out[i++] = (char)c; }
        s++;
    }
    out[i] = '\0';
}

/* Empty-expr sentinel: a lone result with no renderable payload. */
static char* ffi_empty_expr(void) {
    return ffi_strdup("{\"type\":\"expr\",\"payload\":\"\"}");
}

char* mathilda_ffi_eval_json(const char* input) {
    int parse_failed = 0;
    Expr* evaluated = ffi_parse_eval(input, &parse_failed);

    if (parse_failed) {
        /* Echo the (escaped) input so a stray char / bracket mismatch is
         * diagnosable rather than an opaque "parse error" (matches repl.c). */
        const char* in = input ? input : "";
        size_t esc_cap = strlen(in) * 6 + 8;
        char* esc = malloc(esc_cap);
        char* out = NULL;
        if (esc) {
            ffi_json_escape(in, esc, esc_cap);
            size_t bcap = esc_cap + 64;
            out = malloc(bcap);
            if (out)
                snprintf(out, bcap,
                    "{\"type\":\"error\",\"message\":\"Parse error: %s\"}", esc);
        }
        free(esc);
        return out ? out
                   : ffi_strdup("{\"type\":\"error\",\"message\":\"Parse error\"}");
    }
    if (!evaluated) return ffi_empty_expr();

    /* Graphics[...] / Graphics3D[...] → Plotly JSON payload for the notebook.
     * The front end auto-displays a top-level Graphics result, so Plot[...],
     * Show[...], Plot3D[...] and `g // Graphics` all land here. */
    if (evaluated->type == EXPR_FUNCTION
        && evaluated->data.function.head
        && evaluated->data.function.head->type == EXPR_SYMBOL) {
        const char* head_sym = evaluated->data.function.head->data.symbol.name;
        char* plotly = NULL;
        if (head_sym == SYM_Graphics)        plotly = graphics_to_plotly_json(evaluated);
        else if (head_sym == SYM_Graphics3D) plotly = graphics3d_to_plotly_json(evaluated);
        if (plotly) {
            expr_free(evaluated);
            size_t n = strlen(plotly) + 48;
            char* out = malloc(n);
            if (out) snprintf(out, n, "{\"type\":\"plot\",\"payload\":%s}", plotly);
            free(plotly);
            return out ? out : ffi_empty_expr();
        }
        /* Non-convertible Graphics (e.g. empty): emit nothing renderable rather
         * than dumping the raw Graphics[...] tree as red KaTeX. */
        if (head_sym == SYM_Graphics || head_sym == SYM_Graphics3D) {
            expr_free(evaluated);
            return ffi_empty_expr();
        }
    }

    char* text  = expr_to_string(evaluated);
    char* latex = expr_to_latex(evaluated);   /* must be produced before free */
    expr_free(evaluated);
    if (!text) { free(latex); return ffi_empty_expr(); }

    size_t tlen = strlen(text) * 6 + 4;
    char* tesc = malloc(tlen);
    if (tesc) ffi_json_escape(text, tesc, tlen);
    free(text);

    char* lesc = NULL;
    if (latex) {
        size_t llen = strlen(latex) * 6 + 4;
        lesc = malloc(llen);
        if (lesc) ffi_json_escape(latex, lesc, llen);
        free(latex);
    }

    char* out = NULL;
    if (tesc) {
        size_t n = strlen(tesc) + (lesc ? strlen(lesc) : 0) + 64;
        out = malloc(n);
        if (out) {
            if (lesc && lesc[0])
                snprintf(out, n,
                    "{\"type\":\"expr\",\"payload\":\"%s\",\"latex\":\"%s\"}", tesc, lesc);
            else
                snprintf(out, n, "{\"type\":\"expr\",\"payload\":\"%s\"}", tesc);
        }
    }
    free(tesc);
    free(lesc);
    return out ? out : ffi_empty_expr();
}

void mathilda_ffi_free(char* s) {
    free(s);
}

const char* mathilda_ffi_version(void) {
    return MATHILDA_VERSION_STRING;
}
