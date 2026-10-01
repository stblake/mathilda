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
#include <ctype.h>

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
#include "image.h"
#include "meminfo.h"

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
/* Worst-case escaped size of an n-byte string: each byte can become the 6-byte
 * "\uXXXX", plus room for the NUL and a little slack. Use everywhere a buffer
 * is sized for ffi_json_escape so the margin can't drift between call sites. */
#define FFI_JSON_ESC_BOUND(n) ((n) * 6 + 8)
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
        size_t esc_cap = FFI_JSON_ESC_BOUND(strlen(in));
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

    size_t tlen = FFI_JSON_ESC_BOUND(strlen(text));
    char* tesc = malloc(tlen);
    if (tesc) ffi_json_escape(text, tesc, tlen);
    free(text);

    char* lesc = NULL;
    if (latex) {
        size_t llen = FFI_JSON_ESC_BOUND(strlen(latex));
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

/* ====================================================================
 * Notebook-cell evaluation and editor services.
 *
 * In-process twin of repl.c's pipe mode (pipe_process_input / pipe_eval_
 * statement): it emits the SAME NDJSON event lines, so a front end maps one
 * event vocabulary whichever transport it uses. repl.c owns the canonical
 * implementation but lives in repl.o, which is excluded from libmathilda.a —
 * the same reason ffi_json_escape duplicates repl.c's json_escape. The two are
 * held in step by tools/check_pipe_protocol.py (exercising the sidecar) and by
 * the FFI cell tests. A later consolidation can point repl.c's loop at this.
 * ==================================================================== */

typedef struct { mathilda_ffi_sink sink; void* ctx; } FfiSink;

static void ffi_cell_emit(const FfiSink* s, const char* line) {
    if (s && s->sink) s->sink(s->ctx, line);
}

/* Session $Line for cell mode: one per statement, never reset — a fresh kernel
 * process starts a fresh session — so %, %%, In[n], Out[n] resolve as in the
 * REPL. Mirrors repl.c's g_pipe_line. */
static int g_ffi_line = 0;

/* Set $Line before evaluating (a `%` in the statement is Out[$Line - 1]). */
static void ffi_set_line(int line) {
    Expr* sym = expr_new_symbol(SYM_DollarLine);
    Expr* val = expr_new_integer(line);
    symtab_add_own_value("$Line", sym, val);
    expr_free(sym);
    expr_free(val);
}

/* In[line] = parsed / Out[line] = evaluated. symtab_add_down_value copies both
 * arguments; expr_new_function consumes the head and arg handed to it, so only
 * the pattern is ours to free. Both borrow their value argument. */
static void ffi_store_in(int line, Expr* parsed) {
    Expr* arg = expr_new_integer(line);
    Expr* pat = expr_new_function(expr_new_symbol(SYM_In), &arg, 1);
    symtab_add_down_value("In", pat, parsed);
    expr_free(pat);
}
static void ffi_store_out(int line, Expr* evaluated) {
    if (!evaluated) return;
    Expr* arg = expr_new_integer(line);
    Expr* pat = expr_new_function(expr_new_symbol(SYM_Out), &arg, 1);
    symtab_add_down_value("Out", pat, evaluated);
    expr_free(pat);
}

/* `?x` may sit at the tail of a CompoundExpression (`a = 5; ?Sin`); the head to
 * test is the final element. Mirrors repl.c's pipe_final_expr. */
static Expr* ffi_final_expr(Expr* e) {
    while (e && e->type == EXPR_FUNCTION && e->data.function.head
           && e->data.function.head->type == EXPR_SYMBOL
           && e->data.function.head->data.symbol.name == SYM_CompoundExpression
           && e->data.function.arg_count > 0)
        e = e->data.function.args[e->data.function.arg_count - 1];
    return e;
}
static const char* ffi_info_symbol(Expr* parsed) {
    Expr* e = ffi_final_expr(parsed);
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 1) return NULL;
    Expr* a = e->data.function.args[0];
    if (a && a->type == EXPR_SYMBOL) return a->data.symbol.name;
    if (a && a->type == EXPR_STRING) return a->data.string;
    return NULL;
}
static int ffi_is_info_query(Expr* parsed) {
    Expr* e = ffi_final_expr(parsed);
    return e && e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Information;
}

static uint64_t ffi_memory_bytes(void) {
    uint64_t b = 0;
    if (!meminfo_current(&b)) return 0;
    return b;
}

static void ffi_emit_done(const FfiSink* s, int id) {
    char buf[96];
    snprintf(buf, sizeof(buf), "{\"id\":%d,\"type\":\"done\",\"memory\":%llu}",
             id, (unsigned long long)ffi_memory_bytes());
    ffi_cell_emit(s, buf);
}
static void ffi_emit_line(const FfiSink* s, int id, int line) {
    char buf[96];
    snprintf(buf, sizeof(buf), "{\"id\":%d,\"type\":\"line\",\"line\":%d}", id, line);
    ffi_cell_emit(s, buf);
}

/* {"id":N,"type":<kind>,"<field>":"<escaped>"} for `len` bytes of `text`. */
static void ffi_emit_text(const FfiSink* s, int id, const char* kind,
                          const char* field, const char* text, size_t len) {
    char* raw = malloc(len + 1);
    if (!raw) return;
    memcpy(raw, text, len);
    raw[len] = '\0';
    size_t ecap = FFI_JSON_ESC_BOUND(len);
    char* esc = malloc(ecap);
    size_t lcap = ecap + strlen(kind) + strlen(field) + 64;
    char* line = esc ? malloc(lcap) : NULL;
    if (esc && line) {
        ffi_json_escape(raw, esc, ecap);
        snprintf(line, lcap, "{\"id\":%d,\"type\":\"%s\",\"%s\":\"%s\"}", id, kind, field, esc);
        ffi_cell_emit(s, line);
    }
    free(line);
    free(esc);
    free(raw);
}

static void ffi_emit_parse_error(const FfiSink* s, const char* input, int id) {
    size_t esc_cap = FFI_JSON_ESC_BOUND(strlen(input));
    char* esc = malloc(esc_cap);
    char* buf = NULL;
    if (esc) {
        ffi_json_escape(input, esc, esc_cap);
        size_t bcap = esc_cap + 128;
        buf = malloc(bcap);
        if (buf)
            snprintf(buf, bcap,
                "{\"id\":%d,\"type\":\"error\",\"message\":\"Parse error: %s\"}", id, esc);
    }
    if (buf) ffi_cell_emit(s, buf);
    else {
        char sbuf[128];
        snprintf(sbuf, sizeof(sbuf),
            "{\"id\":%d,\"type\":\"error\",\"message\":\"Parse error\"}", id);
        ffi_cell_emit(s, sbuf);
    }
    free(esc);
    free(buf);
}

/* ---- Print/message capture during one statement (POSIX open_memstream) ---- */
typedef struct {
    FILE*  saved_out; FILE* saved_err;
    FILE*  out; FILE* err;
    char*  out_buf; size_t out_len;
    char*  err_buf; size_t err_len;
} FfiCapture;

static void ffi_capture_begin(FfiCapture* c) {
    memset(c, 0, sizeof(*c));
#ifndef _WIN32
    fflush(stdout);
    fflush(stderr);
    c->out = open_memstream(&c->out_buf, &c->out_len);
    c->err = open_memstream(&c->err_buf, &c->err_len);
    if (!c->out || !c->err) {
        if (c->out) fclose(c->out);
        if (c->err) fclose(c->err);
        free(c->out_buf);
        free(c->err_buf);
        memset(c, 0, sizeof(*c));
        return;
    }
    c->saved_out = stdout;
    c->saved_err = stderr;
    stdout = c->out;
    stderr = c->err;
#endif
}

#ifndef _WIN32
/* One "message" line per funnelled Head::tag message; a line not beginning with
 * a non-blank "Head::tag" continues the previous one. Mirrors repl.c. */
static void ffi_emit_messages(const FfiSink* s, int id, const char* text) {
    const char* msg_start = NULL;
    const char* p = text;
    while (*p) {
        const char* eol = strchr(p, '\n');
        const char* next = eol ? eol + 1 : p + strlen(p);
        size_t len = (size_t)((eol ? eol : next) - p);
        int starts = len > 0 && !isspace((unsigned char)*p);
        if (starts) {
            const char* dc = strstr(p, "::");
            starts = dc && dc < p + len;
        }
        if (starts && msg_start) {
            size_t mlen = (size_t)(p - msg_start);
            while (mlen > 0 && (msg_start[mlen - 1] == '\n' || msg_start[mlen - 1] == '\r')) mlen--;
            if (mlen) ffi_emit_text(s, id, "message", "text", msg_start, mlen);
            msg_start = NULL;
        }
        if (!msg_start && len > 0) msg_start = p;
        p = next;
    }
    if (msg_start) {
        size_t mlen = strlen(msg_start);
        while (mlen > 0 && (msg_start[mlen - 1] == '\n' || msg_start[mlen - 1] == '\r')) mlen--;
        if (mlen) ffi_emit_text(s, id, "message", "text", msg_start, mlen);
    }
}
#endif

static void ffi_capture_end(const FfiSink* s, FfiCapture* c, int id) {
#ifndef _WIN32
    if (!c->out) return;
    stdout = c->saved_out;
    stderr = c->saved_err;
    fclose(c->out);
    fclose(c->err);
    if (c->out_buf && c->out_len > 0)
        ffi_emit_text(s, id, "stream", "text", c->out_buf, c->out_len);
    if (c->err_buf && c->err_len > 0)
        ffi_emit_messages(s, id, c->err_buf);
    free(c->out_buf);
    free(c->err_buf);
    memset(c, 0, sizeof(*c));
#else
    (void)s; (void)c; (void)id;
#endif
}

/* Evaluate ONE statement and emit its events. Borrows `parsed`; emits no
 * "done". Mirrors repl.c's pipe_eval_statement. */
static void ffi_eval_statement(const FfiSink* s, Expr* parsed, int id,
                               int show_result, int cell) {
    int info_query = ffi_is_info_query(parsed);
    char info_sym[128];
    info_sym[0] = '\0';
    if (info_query) {
        const char* isname = ffi_info_symbol(parsed);
        if (isname) {
            strncpy(info_sym, isname, sizeof(info_sym) - 1);
            info_sym[sizeof(info_sym) - 1] = '\0';
        }
    }

    int line = 0;
    if (cell) {
        line = ++g_ffi_line;
        ffi_set_line(line);
        ffi_store_in(line, parsed);
        ffi_emit_line(s, id, line);
    }

    FfiCapture cap;
    if (cell) ffi_capture_begin(&cap);
    Expr* evaluated = evaluate(parsed);
    if (cell) ffi_capture_end(s, &cap, id);

    if (cell) ffi_store_out(line, evaluated);
    if (!evaluated) return;

    if (cell && (!show_result
        || (evaluated->type == EXPR_SYMBOL && evaluated->data.symbol.name == SYM_Null))) {
        expr_free(evaluated);
        return;
    }

    /* ?Pat* -> the List of matching names, as a "names" event. */
    if (info_query && evaluated->type == EXPR_FUNCTION
        && evaluated->data.function.head
        && evaluated->data.function.head->type == EXPR_SYMBOL
        && evaluated->data.function.head->data.symbol.name == SYM_List) {
        size_t n = evaluated->data.function.arg_count;
        size_t cap_n = 64;
        for (size_t i = 0; i < n; i++) {
            Expr* e = evaluated->data.function.args[i];
            if (e->type == EXPR_STRING) cap_n += FFI_JSON_ESC_BOUND(strlen(e->data.string));
        }
        char* buf = malloc(cap_n);
        if (buf) {
            int off = snprintf(buf, cap_n, "{\"id\":%d,\"type\":\"names\",\"payload\":[", id);
            int first = 1;
            for (size_t i = 0; i < n && off > 0 && (size_t)off < cap_n; i++) {
                Expr* e = evaluated->data.function.args[i];
                if (e->type != EXPR_STRING) continue;
                size_t ecap = FFI_JSON_ESC_BOUND(strlen(e->data.string));
                char* esc = malloc(ecap);
                if (!esc) break;
                ffi_json_escape(e->data.string, esc, ecap);
                off += snprintf(buf + off, cap_n - (size_t)off, "%s\"%s\"", first ? "" : ",", esc);
                free(esc);
                first = 0;
            }
            if (off > 0 && (size_t)off < cap_n) snprintf(buf + off, cap_n - (size_t)off, "]}");
            ffi_cell_emit(s, buf);
            free(buf);
        }
        expr_free(evaluated);
        return;
    }

    /* ?sym -> the raw docstring as a "usage" event (not run through a renderer). */
    if (info_query && evaluated->type == EXPR_STRING) {
        const char* doc = evaluated->data.string;
        size_t dcap = FFI_JSON_ESC_BOUND(strlen(doc));
        char* esc = malloc(dcap);
        if (esc) {
            ffi_json_escape(doc, esc, dcap);
            size_t bcap = dcap + 256;
            char* buf = malloc(bcap);
            if (buf) {
                if (info_sym[0])
                    snprintf(buf, bcap,
                        "{\"id\":%d,\"type\":\"usage\",\"payload\":\"%s\",\"symbol\":\"%s\"}",
                        id, esc, info_sym);
                else
                    snprintf(buf, bcap,
                        "{\"id\":%d,\"type\":\"usage\",\"payload\":\"%s\"}", id, esc);
                ffi_cell_emit(s, buf);
                free(buf);
            }
            free(esc);
        }
        expr_free(evaluated);
        return;
    }

    /* Image[...] / Image3D[...] -> base64 RGBA "image" event. */
    if (evaluated->type == EXPR_FUNCTION
        && evaluated->data.function.head
        && evaluated->data.function.head->type == EXPR_SYMBOL) {
        const char* ih = evaluated->data.function.head->data.symbol.name;
        if (ih && (strcmp(ih, "Image") == 0 || strcmp(ih, "Image3D") == 0)) {
            char* ijson = image_to_json(evaluated);
            if (ijson) {
                expr_free(evaluated);
                size_t jl = strlen(ijson) + 64;
                char* jline = malloc(jl);
                if (jline) {
                    snprintf(jline, jl, "{\"id\":%d,\"type\":\"image\",\"payload\":%s}", id, ijson);
                    ffi_cell_emit(s, jline);
                    free(jline);
                }
                free(ijson);
                return;
            }
        }
    }

    /* Graphics[...] / Graphics3D[...] -> Plotly JSON "plot" event. */
    if (evaluated->type == EXPR_FUNCTION
        && evaluated->data.function.head
        && evaluated->data.function.head->type == EXPR_SYMBOL) {
        const char* head_sym = evaluated->data.function.head->data.symbol.name;
        char* plotly = NULL;
        if (head_sym == SYM_Graphics)        plotly = graphics_to_plotly_json(evaluated);
        else if (head_sym == SYM_Graphics3D) plotly = graphics3d_to_plotly_json(evaluated);
        if (plotly) {
            expr_free(evaluated);
            size_t jl = strlen(plotly) + 64;
            char* jline = malloc(jl);
            if (jline) {
                snprintf(jline, jl, "{\"id\":%d,\"type\":\"plot\",\"payload\":%s}", id, plotly);
                ffi_cell_emit(s, jline);
                free(jline);
            }
            free(plotly);
            return;
        }
        if (head_sym == SYM_Graphics || head_sym == SYM_Graphics3D) {
            expr_free(evaluated);
            return;
        }
    }

    /* Ordinary result: "expr" with plain payload and optional KaTeX latex. */
    char* result_str = expr_to_string(evaluated);
    char* latex_raw  = expr_to_latex(evaluated);   /* before expr_free */
    expr_free(evaluated);

    if (!result_str) {
        free(latex_raw);
        char buf[128];
        snprintf(buf, sizeof(buf),
                 "{\"id\":%d,\"type\":\"error\",\"message\":\"Out of memory\"}", id);
        ffi_cell_emit(s, buf);
        return;
    }
    size_t escaped_len = FFI_JSON_ESC_BOUND(strlen(result_str));
    char* escaped = malloc(escaped_len);
    if (!escaped) { free(result_str); free(latex_raw); return; }
    ffi_json_escape(result_str, escaped, escaped_len);
    free(result_str);

    char* latex_esc = NULL;
    if (latex_raw) {
        size_t llen = FFI_JSON_ESC_BOUND(strlen(latex_raw));
        latex_esc = malloc(llen);
        if (latex_esc) ffi_json_escape(latex_raw, latex_esc, llen);
        free(latex_raw);
    }

    size_t line_len = escaped_len + (latex_esc ? strlen(latex_esc) : 0) + 128;
    char* json_line = malloc(line_len);
    if (json_line) {
        if (latex_esc && latex_esc[0])
            snprintf(json_line, line_len,
                     "{\"id\":%d,\"type\":\"expr\",\"payload\":\"%s\",\"latex\":\"%s\"}",
                     id, escaped, latex_esc);
        else
            snprintf(json_line, line_len,
                     "{\"id\":%d,\"type\":\"expr\",\"payload\":\"%s\"}", id, escaped);
        ffi_cell_emit(s, json_line);
        free(json_line);
    }
    free(escaped);
    free(latex_esc);
}

/* True if `s` holds only whitespace and (* ... *) comments (so trailing such
 * after the last statement is not a syntax error). */
static int ffi_blanks_only(const char* s) {
    while (*s) {
        if (isspace((unsigned char)*s)) { s++; continue; }
        if (s[0] == '(' && s[1] == '*') {
            int depth = 1; s += 2;
            while (*s && depth > 0) {
                if (s[0] == '(' && s[1] == '*') { depth++; s += 2; }
                else if (s[0] == '*' && s[1] == ')') { depth--; s += 2; }
                else s++;
            }
            if (depth > 0) return 0;  /* unterminated comment is a real error */
        } else {
            return 0;
        }
    }
    return 1;
}

void mathilda_ffi_eval_cell(const char* code, int id, int cell,
                            mathilda_ffi_sink sink, void* ctx) {
    if (!g_initialized) mathilda_ffi_init();
    FfiSink s = { sink, ctx };
    if (!code) { ffi_emit_parse_error(&s, "", id); ffi_emit_done(&s, id); return; }

    if (!cell) {
        Expr* parsed = parse_expression(code);
        if (!parsed) ffi_emit_parse_error(&s, code, id);
        else { ffi_eval_statement(&s, parsed, id, 1, 0); expr_free(parsed); }
        ffi_emit_done(&s, id);
        return;
    }

    /* Parse the whole cell FIRST (a syntax error anywhere evaluates nothing),
     * splitting statements like a Mathematica input cell. */
    size_t n = 0, capn = 8;
    Expr** stmts = malloc(capn * sizeof(Expr*));
    int* shown = malloc(capn * sizeof(int));
    int failed = (stmts == NULL || shown == NULL);
    const char* p = code;
    while (!failed) {
        const char* start = p;
        Expr* e = parse_next_expression(&p);
        if (!e) {
            if (!ffi_blanks_only(start)) failed = 1;
            break;
        }
        if (n == capn) {
            size_t nc = capn * 2;
            Expr** ns = realloc(stmts, nc * sizeof(Expr*));
            if (ns) stmts = ns;
            int* nb = realloc(shown, nc * sizeof(int));
            if (nb) shown = nb;
            if (!ns || !nb) { expr_free(e); failed = 1; break; }
            capn = nc;
        }
        stmts[n] = e;
        shown[n] = !(p > start && p[-1] == ';');
        n++;
    }

    if (failed) ffi_emit_parse_error(&s, code, id);
    else for (size_t i = 0; i < n; i++) ffi_eval_statement(&s, stmts[i], id, shown[i], 1);
    for (size_t i = 0; i < n; i++) expr_free(stmts[i]);
    free(stmts);
    free(shown);
    ffi_emit_done(&s, id);
}

/* ---- Editor services: completeness check and name completion -------------- */

int mathilda_ffi_is_complete(const char* code) {
    if (!code) return 1;
    int depth = 0, comment = 0, in_string = 0;
    for (const char* p = code; *p; ) {
        if (in_string) {
            if (*p == '\\' && p[1]) { p += 2; continue; }
            if (*p == '"') in_string = 0;
            p++;
        } else if (comment > 0) {
            if (p[0] == '(' && p[1] == '*') { comment++; p += 2; }
            else if (p[0] == '*' && p[1] == ')') { comment--; p += 2; }
            else p++;
        } else if (p[0] == '(' && p[1] == '*') {
            comment++; p += 2;
        } else if (*p == '"') {
            in_string = 1; p++;
        } else {
            if (*p == '(' || *p == '[' || *p == '{') depth++;
            else if (*p == ')' || *p == ']' || *p == '}') depth--;
            p++;
        }
    }
    return !in_string && comment == 0 && depth <= 0;
}

/* Collector for symtab_for_each: names beginning with `prefix`. */
typedef struct { const char* prefix; size_t plen; char** names; size_t n, cap; } FfiComplete;
static void ffi_complete_visit(const char* name, SymbolDef* def, void* user) {
    (void)def;
    FfiComplete* c = (FfiComplete*)user;
    if (!name) return;
    if (c->plen && strncmp(name, c->prefix, c->plen) != 0) return;
    if (c->n == c->cap) {
        size_t nc = c->cap ? c->cap * 2 : 32;
        char** nn = realloc(c->names, nc * sizeof(char*));
        if (!nn) return;
        c->names = nn; c->cap = nc;
    }
    size_t len = strlen(name) + 1;
    char* dup = malloc(len);
    if (!dup) return;
    memcpy(dup, name, len);
    c->names[c->n++] = dup;
}
static int ffi_strcmp_ptr(const void* a, const void* b) {
    return strcmp(*(const char* const*)a, *(const char* const*)b);
}

char* mathilda_ffi_complete(const char* prefix) {
    if (!g_initialized) mathilda_ffi_init();
    FfiComplete c = { prefix ? prefix : "", prefix ? strlen(prefix) : 0, NULL, 0, 0 };
    symtab_for_each(ffi_complete_visit, &c);
    if (c.n > 1) qsort(c.names, c.n, sizeof(char*), ffi_strcmp_ptr);

    size_t cap = 4;
    for (size_t i = 0; i < c.n; i++) cap += FFI_JSON_ESC_BOUND(strlen(c.names[i])) + 3;
    char* out = malloc(cap);
    if (out) {
        size_t off = 0;
        out[off++] = '[';
        for (size_t i = 0; i < c.n; i++) {
            size_t ecap = FFI_JSON_ESC_BOUND(strlen(c.names[i]));
            char* esc = malloc(ecap);
            if (esc) {
                ffi_json_escape(c.names[i], esc, ecap);
                off += (size_t)snprintf(out + off, cap - off, "%s\"%s\"", i ? "," : "", esc);
                free(esc);
            }
        }
        if (off + 2 <= cap) { out[off++] = ']'; out[off] = '\0'; }
        else { free(out); out = NULL; }
    }
    for (size_t i = 0; i < c.n; i++) free(c.names[i]);
    free(c.names);
    return out ? out : ffi_strdup("[]");
}

void mathilda_ffi_free(char* s) {
    free(s);
}

const char* mathilda_ffi_version(void) {
    return MATHILDA_VERSION_STRING;
}
