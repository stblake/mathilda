/* fileno()/isatty() are POSIX, hidden by glibc under -std=c99; request them. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "regex_common.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "part.h"
#include "eval.h"
#include "message.h"
#include "symtab.h"
#include "repl_hooks.h"
#include "sym_names.h"
#include "show.h"
#include "render3d.h"
#include "graphics_json.h"
#include "image.h"
#include "meminfo.h"
#include "print_latex.h"
#include "version.h"
#include "quit.h"
#include "dialog.h"
#include "term/term_input.h"
#include "term/term_interrupt.h"
#include "term/term_readline.h"
#include "term/term_editor_tty.h"
#include "term/term_highlight.h"
#include "term/term_history.h"
#include "term/term_color.h"
#include "mtime.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
void mth_block_graveyard_flush(void);   /* modular.c: free Block values parked by a TimeConstrained unwind */

/* Portable isatty + fileno for pipe-mode detection */
#ifdef _WIN32
  #include <io.h>
  #ifndef isatty
    #define isatty _isatty
  #endif
  #ifndef fileno
    #define fileno _fileno
  #endif
#else
  #include <unistd.h>
#endif


#define MAX_INPUT_LEN 10240

/* Session history ($Line, In[n], Out[n]) and the blank-input helpers live in
 * dialog.c, shared with Dialog[]'s input loop. */

void process_input(const char* input, int line_number) {
    if (strlen(input) == 0) return;
    if (mth_is_blank_or_comment_only(input)) return;

    mth_session_set_line(line_number);

    /* $PreRead: text-level hook applied to the raw input string
     * before parsing. Pass-through (a strdup of `input`) when unset. */
    char* cooked_input = repl_apply_pre_read(input);
    if (!cooked_input) {
        printf("Parse error\n\n");
        return;
    }

    // Parse the (possibly hook-transformed) input
    Expr* parsed = parse_expression(cooked_input);
    free(cooked_input);
    if (!parsed) {
        printf("Parse error\n\n");
        return;
    }

    mth_session_store_in(line_number, parsed);

    /* $Pre: applied to the parsed expression. Consumes our reference
     * to `parsed`; we treat the result as the new input to evaluate. */
    Expr* pre_input = repl_apply_pre(expr_copy(parsed));
    double t0 = mth_mono_seconds();
    Expr* evaluated = evaluate(pre_input);
    double took = mth_mono_seconds() - t0;
    expr_free(pre_input);
    mth_block_graveyard_flush();   /* top-level boundary: nothing in flight */

    /* Quit[] was evaluated: the session is ending, so there is no Out[n] to
     * store or print. The caller sees mth_quit_requested() and shuts down. */
    if (mth_quit_requested()) {
        expr_free(evaluated);
        expr_free(parsed);
        return;
    }

    /* $Post: applied to the evaluator's result. Consumes our reference
     * to `evaluated`. */
    evaluated = repl_apply_post(evaluated);

    /* A NULL result means the evaluation produced nothing displayable
     * (e.g. a hook absorbed the value). Skip Out[n] storage and the
     * "Out[n]= " banner rather than crashing in expr_copy. */
    if (!evaluated) {
        expr_free(parsed);
        return;
    }

    /* Store Out[line_number] = evaluated (post-$Post, pre-$PrePrint:
     * Mathematica's documented ordering). The full `evaluated` is kept for the
     * formatted display below. */
    mth_session_store_out(line_number, evaluated);

    /* $PrePrint: applied only for display. Out[n] keeps the
     * pre-$PrePrint value above; here we render a possibly modified
     * copy. */
    Expr* to_print = repl_apply_pre_print(expr_copy(evaluated));

    /* `?sym` / Information[sym] yields the raw docstring as a String. Print it
     * as a formatted usage message (real newlines/tabs, no surrounding quotes
     * or InputForm escaping) — Mathematica's behavior — rather than as a quoted
     * string literal. Keyed on the *input* head so only help queries take this
     * path; an ordinary string result still prints quoted. */
    int is_info_query =
        parsed && parsed->type == EXPR_FUNCTION
        && parsed->data.function.head->type == EXPR_SYMBOL
        && parsed->data.function.head->data.symbol.name == SYM_Information
        && to_print && to_print->type == EXPR_STRING;

    /* Colour (and the elapsed-time note) only in the interactive session:
     * term_color_init runs nowhere else, so these are "" / false in -file. */
    bool c = term_color_stdout();
    printf("%sOut[%d]=%s ", term_style_on(TERM_STYLE_OUT, c), line_number, term_style_off(c));
    if (is_info_query) fputs(to_print->data.string, stdout);
    else if (!(c && term_editor_session_active() && term_hl_print_result(to_print)))
        expr_print(to_print);
    if (term_tty_stdout() && took >= 1.0) {
        /* A dim note on its own line, under the result. */
        char el[64];
        term_format_elapsed(el, sizeof(el), took);
        printf("%stook %s%s\n", term_style_on(TERM_STYLE_DIM, c), el, term_style_off(c));
    }
    printf("\n"); // extra blank line

    /* Mathematica's front end auto-displays a top-level Graphics[...] (or
     * Graphics3D[...], from Plot3D) result. This REPL is the sole "front
     * end", so it owns rendering: Show[]/Plot[]/Plot3D[] merely return such
     * an object and we render it here. Routing every display through one
     * path means `g // Graphics`, Show[...], Plot[...] and Plot3D[...] all
     * render identically, and a trailing `;` (which yields Null) correctly
     * suppresses the window. graphics_show/graphics3d_show borrow the expr
     * (no ownership transfer); on a non-graphics build their stubs print a
     * one-line "install raylib" hint instead. */
    /* MATHILDA_NO_WINDOW suppresses the display, for callers that evaluate expressions in bulk and do
     * not want a window per Graphics result. Documentation generation is the case that forced this:
     * site/generate.py re-verifies every documented example against this binary, and the ones calling
     * Plot or Manipulate opened a real Raylib window each -- dozens of them, over the user's work,
     * during what should be a silent batch job. The expression still evaluates and still prints; only
     * the window is withheld. */
    if (to_print && to_print->type == EXPR_FUNCTION
        && to_print->data.function.head->type == EXPR_SYMBOL
        && to_print->data.function.arg_count >= 1
        && getenv("MATHILDA_NO_WINDOW") == NULL) {
        if (to_print->data.function.head->data.symbol.name == SYM_Graphics) graphics_show(to_print);
        else if (to_print->data.function.head->data.symbol.name == SYM_Graphics3D) graphics3d_show(to_print);
    }

    expr_free(to_print);
    expr_free(parsed);
    expr_free(evaluated);
}

/* process_input for the interactive loop: the interrupt state is EVALUATING
 * for its duration (a Ctrl-C aborts the evaluation), and no interrupt or
 * top-level abort survives it in either direction -- a Ctrl-C that arrives
 * after the result was printed must not abort the next input. */
static void repl_process(const char* line, int line_number) {
    TermState prev = term_interrupt_set_state(TERM_STATE_EVALUATING);
    mth_interrupt_clear();
    (void)eval_take_top_abort();
    process_input(line, line_number);
    mth_interrupt_clear();
    term_interrupt_set_state(prev);
}

static void repl_banner(TermBackendKind kind) {
    printf("\nMathilda " MATHILDA_VERSION_STRING " - A small, open source computer algebra system.\n\n");
    if (kind != TERM_BACKEND_STDIO) {
        printf("This program is free, open source software and comes with ABSOLUTELY NO WARRANTY.\n\n");
        printf("Press Return to evaluate. An open bracket, string or comment continues\n");
        printf("on the next line; press Esc then Return to force evaluation.\n");
        if (kind == TERM_BACKEND_EDITOR)
            printf("Tab completes, Ctrl-R searches the history, Ctrl-_ undoes.\n");
        printf("Ctrl-C interrupts an evaluation (twice within a second force-quits).\n");
        printf("Exit by evaluating Quit[] or pressing Ctrl-D.\n\n");
    } else {
        printf("Exit by evaluating Quit[] or pressing Ctrl+Z (Windows) / Ctrl+D (Unix).\n\n");
    }
}

/* The interactive loop, over any line reader (term_input.h): readline with
 * its smart Return key, or plain fgets in a NO_READLINE build. With the smart
 * Return binding a single read returns the whole (possibly multi-line)
 * expression, so no accumulation buffer is needed. */
void repl_loop(TermReader* rd) {
    repl_banner(term_reader_kind(rd));
    int line_number = 1;
    char prompt[64];
    for (;;) {
        snprintf(prompt, sizeof(prompt), "In[%d]:= ", line_number);
        char* line = NULL;
        TermReadStatus st = term_read(rd, prompt, &line);
        if (st == TERM_READ_INTERRUPTED) continue;        /* same In[n] again */
        if (st == TERM_READ_EOF) {
            printf("\n");
            repl_apply_epilog();                          /* EOF: $Epilog first */
            break;
        }
        if (line[0] == '\0') { free(line); continue; }
        /* Before evaluating: a crash or a forced quit keeps the entry. */
        term_reader_history_add(rd, line);
        repl_process(line, line_number++);
        free(line);
        if (mth_quit_requested()) {
            /* Quit[] / Quit[n] evaluated anywhere in the input (a bare
             * Quit[], If[c, Quit[]], a Get'd file, ...): run $Epilog first. */
            repl_apply_epilog();
            break;
        }
    }
    printf("\n");
}

#include "core.h"
#include "loadmodule.h"
#include "context.h"

/* =====================================================================
 * Minimal NDJSON pipe-mode protocol
 *
 * When stdin is not a terminal (i.e. the frontend spawned us as a
 * sidecar), switch from the readline REPL to a simple line-based
 * protocol over stdio.
 *
 * Request  (one line on stdin, any length):
 *   {"id": N, "expr": "1+1"}               -- evaluate one expression
 *   {"id": N, "expr": "...", "cell": true} -- evaluate a notebook cell
 *   {"type": "ping"}                        -- readiness probe
 *   {"type": "quit"}                        -- graceful shutdown
 *
 * Response (one JSON object per line on stdout), for each request:
 *   {"id": N, "type": "stream",  "text": "hello\n"}          (cell only)
 *   {"id": N, "type": "message", "text": "Power::infy: ..."} (cell only)
 *   {"id": N, "type": "expr",  "payload": "2", "latex": "2", "slatex": "2"}
 *     latex  = TraditionalForm LaTeX, slatex = StandardForm LaTeX (the
 *     notebook's default display); either is absent when it has no typeset
 *     form (print_latex.h), and slatex alone for an outer TraditionalForm[...].
 *     slatex is the SELECTABLE form: its subexpressions carry
 *     \htmlData{mp=K}, and "sparts": [...] holds each one's InputForm text
 *     (sparts[K]; absent when the result is too large to tag).
 *   {"id": N, "type": "error", "message": "Parse error: ..."}
 *   ... (usage / names / image / plot, see pipe_eval_statement)
 *   {"id": N, "type": "done", "memory": BYTES}               -- always last
 *   {"type": "pong"}
 *
 * A plain request is one expression, and its behaviour is unchanged from
 * the protocol's first version: Print text goes to stdout raw, between
 * protocol lines; messages go to stderr; a Null result is sent as the
 * payload "Null". Tools depend on all three (site/generate.py keeps a
 * `x = ...;` setup line BECAUSE it sees "Null"; the audit tools and
 * book/tools/gen_compileprint.py read raw Print lines), so the notebook's
 * semantics are opt-in rather than a change under them.
 *
 * With "cell": true the request is a notebook input cell:
 *   - it may hold several statements, one per line or separated by ';',
 *     exactly as a Mathematica input cell does, each evaluated in turn;
 *   - each statement's Print output and messages are captured and sent as
 *     "stream" and "message" lines BEFORE its result, so they reach the
 *     cell instead of being dropped (stdout) or only logged (stderr);
 *   - a statement ending in ';', or evaluating to Null, sends no result.
 * An older kernel ignores the unknown "cell" key and answers as before,
 * and an older front end never sends it, so either side can be upgraded
 * alone.
 *
 * Protocol lines are written to the stdout the process started with
 * (g_pipe_out), never to whatever `stdout` names at the time: during an
 * evaluation `stdout` and `stderr` are pointed at in-memory streams to
 * capture Print output and messages (see PipeCapture).
 * ===================================================================*/

/* The real stdout, captured before any evaluation can redirect `stdout`. */
static FILE* g_pipe_out = NULL;

/* $Line for this session: the interactive loop's `line_number`, for a front end
 * that does not own a prompt. One per statement evaluated, never reset -- the
 * front end restarts the kernel to clear the session, which starts a fresh
 * process and with it a fresh counter. */
static int g_pipe_line = 0;

static void pipe_emit(const char* line) {
    FILE* out = g_pipe_out ? g_pipe_out : stdout;
    fputs(line, out);
    fputc('\n', out);
    fflush(out);
}

/* Append the UTF-8 encoding of code point `cp` to `buf` at `*i`. */
static void utf8_append(char* buf, size_t* i, unsigned long cp) {
    if (cp < 0x80) {
        buf[(*i)++] = (char)cp;
    } else if (cp < 0x800) {
        buf[(*i)++] = (char)(0xC0 | (cp >> 6));
        buf[(*i)++] = (char)(0x80 | (cp & 0x3F));
    } else if (cp < 0x10000) {
        buf[(*i)++] = (char)(0xE0 | (cp >> 12));
        buf[(*i)++] = (char)(0x80 | ((cp >> 6) & 0x3F));
        buf[(*i)++] = (char)(0x80 | (cp & 0x3F));
    } else {
        buf[(*i)++] = (char)(0xF0 | (cp >> 18));
        buf[(*i)++] = (char)(0x80 | ((cp >> 12) & 0x3F));
        buf[(*i)++] = (char)(0x80 | ((cp >> 6) & 0x3F));
        buf[(*i)++] = (char)(0x80 | (cp & 0x3F));
    }
}

/* Four hex digits at `p` as a number, or -1. */
static long hex4(const char* p) {
    long v = 0;
    for (int k = 0; k < 4; k++) {
        char c = p[k];
        v <<= 4;
        if (c >= '0' && c <= '9')      v |= c - '0';
        else if (c >= 'a' && c <= 'f') v |= c - 'a' + 10;
        else if (c >= 'A' && c <= 'F') v |= c - 'A' + 10;
        else return -1;
    }
    return v;
}

/* The decoded string value of `key` in the one-line JSON object `json`, as a
 * fresh malloc'd string (caller frees), or NULL if absent or not a string.
 *
 * Heap-allocated to the input's own length -- a decoded JSON string is never
 * longer than its encoding -- so a cell is not cut off at a fixed buffer size.
 * (It was: request and value both lived in 10 KB stack buffers.) Decodes every
 * JSON escape, including \b, \f and \uXXXX with surrogate pairs, since serde
 * writes control characters as \u00XX. */
static char* json_get_string_dup(const char* json, const char* key) {
    char search[256];
    snprintf(search, sizeof(search), "\"%s\"", key);
    const char* p = strstr(json, search);
    if (!p) return NULL;
    p += strlen(search);
    while (*p == ' ' || *p == '\t' || *p == ':') p++;
    if (*p != '"') return NULL;
    p++;
    char* buf = malloc(strlen(p) + 1);
    if (!buf) return NULL;
    size_t i = 0;
    while (*p && *p != '"') {
        if (*p == '\\' && *(p + 1)) {
            p++;
            switch (*p) {
                case '"':  buf[i++] = '"';  break;
                case '\\': buf[i++] = '\\'; break;
                case '/':  buf[i++] = '/';  break;
                case 'b':  buf[i++] = '\b'; break;
                case 'f':  buf[i++] = '\f'; break;
                case 'n':  buf[i++] = '\n'; break;
                case 'r':  buf[i++] = '\r'; break;
                case 't':  buf[i++] = '\t'; break;
                case 'u': {
                    long cp = hex4(p + 1);
                    if (cp < 0) { buf[i++] = 'u'; break; }
                    p += 4;
                    if (cp >= 0xD800 && cp <= 0xDBFF && p[1] == '\\' && p[2] == 'u') {
                        long lo = hex4(p + 3);
                        if (lo >= 0xDC00 && lo <= 0xDFFF) {
                            cp = 0x10000 + ((cp - 0xD800) << 10) + (lo - 0xDC00);
                            p += 6;
                        }
                    }
                    if (cp == 0) cp = 0xFFFD;  /* an embedded NUL would end the C string */
                    utf8_append(buf, &i, (unsigned long)cp);
                    break;
                }
                default:   buf[i++] = *p;   break;
            }
        } else {
            buf[i++] = *p;
        }
        p++;
    }
    buf[i] = '\0';
    return buf;
}

static int json_get_int(const char* json, const char* key, int* out) {
    char search[256];
    snprintf(search, sizeof(search), "\"%s\"", key);
    const char* p = strstr(json, search);
    if (!p) return 0;
    p += strlen(search);
    while (*p == ' ' || *p == '\t' || *p == ':') p++;
    if (!(*p == '-' || isdigit((unsigned char)*p))) return 0;
    *out = (int)strtol(p, NULL, 10);
    return 1;
}

/* True when `key` is present in the one-line JSON object `json` with the
 * literal value true. */
static bool json_get_true(const char* json, const char* key) {
    char search[256];
    snprintf(search, sizeof(search), "\"%s\"", key);
    const char* p = strstr(json, search);
    if (!p) return false;
    p += strlen(search);
    while (*p == ' ' || *p == '\t' || *p == ':') p++;
    return strncmp(p, "true", 4) == 0;
}

static void json_escape(const char* s, char* out, size_t outlen) {
    size_t i = 0;
    while (*s && i + 7 < outlen) {
        unsigned char c = (unsigned char)*s;
        if (c == '"') {
            out[i++] = '\\'; out[i++] = '"';
        } else if (c == '\\') {
            out[i++] = '\\'; out[i++] = '\\';
        } else if (c == '\n') {
            out[i++] = '\\'; out[i++] = 'n';
        } else if (c == '\r') {
            out[i++] = '\\'; out[i++] = 'r';
        } else if (c == '\t') {
            out[i++] = '\\'; out[i++] = 't';
        } else if (c < 0x20) {
            i += (size_t)snprintf(out + i, outlen - i, "\\u%04x", (unsigned)c);
        } else {
            out[i++] = (char)c;
        }
        s++;
    }
    out[i] = '\0';
}

/* `,"key":"<escaped value>"` for an optional JSON string field, or NULL when
 * the value is absent or empty (the field is then simply left out). Caller
 * frees. */
static char* pipe_json_field(const char* key, const char* value) {
    if (!value || !value[0]) return NULL;
    size_t vlen = strlen(value) * 6 + 4;
    size_t n = vlen + strlen(key) + 8;
    char* out = malloc(n);
    if (!out) return NULL;
    int k = snprintf(out, n, ",\"%s\":\"", key);
    json_escape(value, out + k, vlen);
    strcat(out, "\"");
    return out;
}

/* `,"key":["<escaped>", ...]` for an optional array of JSON strings, or NULL
 * when it is empty (the field is then left out). Caller frees. */
static char* pipe_json_str_array(const char* key, char* const* items, size_t n) {
    if (n == 0) return NULL;
    size_t cap = strlen(key) + 8;
    for (size_t i = 0; i < n; i++) cap += strlen(items[i]) * 6 + 8;
    char* out = malloc(cap);
    if (!out) return NULL;
    size_t k = (size_t)snprintf(out, cap, ",\"%s\":[", key);
    for (size_t i = 0; i < n; i++) {
        if (i) out[k++] = ',';
        out[k++] = '"';
        json_escape(items[i], out + k, cap - k);
        k += strlen(out + k);
        out[k++] = '"';
    }
    out[k++] = ']';
    out[k] = '\0';
    return out;
}

/* `?x` may sit at the end of a CompoundExpression -- `a = 5; ?Sin` -- whose
 * value IS that last element, so the head to test is the final one, not the
 * top-level CompoundExpression. Without unwrapping, `D[x,x]; ?Find*` fell
 * through to the ordinary expression path and the front end tried to typeset a
 * help result as mathematics. */
static Expr* pipe_final_expr(Expr* e) {
    while (e && e->type == EXPR_FUNCTION && e->data.function.head
           && e->data.function.head->type == EXPR_SYMBOL
           && e->data.function.head->data.symbol.name == SYM_CompoundExpression
           && e->data.function.arg_count > 0) {
        e = e->data.function.args[e->data.function.arg_count - 1];
    }
    return e;
}

/* The symbol `?name` asked about, or NULL. Borrowed -- do not free.
 *
 * The usage message carried only the text, so the notebook had a docstring and no way to know which
 * symbol it described, and therefore could not offer a link to that symbol's page. */
/* The kernel's resident bytes, for the notebook's status bar.
 *
 * Attached to the `done` message rather than exposed as a separate request: memory can only have
 * changed because something was evaluated, `done` is sent exactly then, and a poll would either
 * lag or add traffic for a number that was already available. Zero when the platform cannot report
 * it, which the front end shows as nothing rather than as "0 B". */
static uint64_t pipe_memory_bytes(void) {
    uint64_t b = 0;
    if (!meminfo_current(&b)) return 0;
    return b;
}

static const char* pipe_info_symbol(Expr* parsed) {
    Expr* e = pipe_final_expr(parsed);
    if (!e || e->type != EXPR_FUNCTION || e->data.function.arg_count != 1) return NULL;
    Expr* a = e->data.function.args[0];
    if (a && a->type == EXPR_SYMBOL) return a->data.symbol.name;
    if (a && a->type == EXPR_STRING) return a->data.string;
    return NULL;
}

static bool pipe_is_info_query(Expr* parsed) {
    Expr* e = pipe_final_expr(parsed);
    return e && e->type == EXPR_FUNCTION && e->data.function.head
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Information;
}

static void pipe_emit_done(int id) {
    char buf[96];
    snprintf(buf, sizeof(buf), "{\"id\":%d,\"type\":\"done\",\"memory\":%llu}",
             id, (unsigned long long)pipe_memory_bytes());
    pipe_emit(buf);
}

/* The session line this statement took, sent before its output.
 *
 * The notebook labels the cell In[n] with it. That label has to be the KERNEL's
 * $Line and not a count the front end keeps: `%3`, `In[3]` and `Out[3]` all
 * resolve against $Line, and one kernel serves every notebook on the canvas, so
 * a per-notebook counter would label a cell In[2] while `%2` addressed a line
 * from a different notebook entirely.
 *
 * Its own message rather than a field on the result, because a statement need
 * not produce one: `x = 5;` and a bare Print[] send no result at all, and their
 * lines count just the same. */
static void pipe_emit_line(int id, int line) {
    char buf[96];
    snprintf(buf, sizeof(buf), "{\"id\":%d,\"type\":\"line\",\"line\":%d}", id, line);
    pipe_emit(buf);
}

/* {"id":N,"type":<kind>,"<field>":"<escaped text>"} for `len` bytes of `text`.
 * Only reached from the fd-redirect capture path (pipe_emit_messages /
 * pipe_capture_end), which is itself #ifndef _WIN32; guarded to match so the
 * MinGW build does not see it as an unused static function (-Werror). */
#ifndef _WIN32
static void pipe_emit_text(int id, const char* kind, const char* field,
                           const char* text, size_t len) {
    char* raw = malloc(len + 1);
    if (!raw) return;
    memcpy(raw, text, len);
    raw[len] = '\0';
    size_t ecap = len * 6 + 8;
    char* esc = malloc(ecap);
    size_t lcap = ecap + strlen(kind) + strlen(field) + 64;
    char* line = esc ? malloc(lcap) : NULL;
    if (esc && line) {
        json_escape(raw, esc, ecap);
        snprintf(line, lcap, "{\"id\":%d,\"type\":\"%s\",\"%s\":\"%s\"}", id, kind, field, esc);
        pipe_emit(line);
    }
    free(line);
    free(esc);
    free(raw);
}
#endif /* !_WIN32 */

/* ---------------------------------------------------------------------
 * Capturing Print output and messages during one statement.
 *
 * Print writes to `stdout` and every message funnels through
 * mth_message_v to `stderr`. In pipe mode stdout IS the protocol channel,
 * so a Print used to land between protocol lines as raw text that the
 * front end discarded, and messages went to stderr, which it only logged.
 * For the length of each evaluation both the stdout and stderr file
 * DESCRIPTORS (1 and 2) are redirected to anonymous temp files; afterwards
 * the captured bytes are sent as "stream" and "message" lines ahead of the
 * statement's result.
 *
 * The redirect is at the descriptor level (dup2) rather than by reassigning
 * `stdout`/`stderr`: those are `FILE *const` on musl libc and cannot be
 * assigned (an earlier version did, and broke the aarch64-musl build).  The
 * descriptor approach is also a stricter catch-all — it captures writes from
 * any code or library that goes to fd 1/2, not only the FILE* handles.
 *
 * The cost is that output arrives when the statement finishes rather than
 * as it is printed, and that Print text and messages of ONE statement are
 * not interleaved with each other (Print text first). Streaming would need
 * a writer thread on a pipe, which is not worth it for a notebook that
 * shows the cell's output at the end anyway.
 *
 * Windows has no fork/dup2 plumbing here; there the streams are left alone.
 * --------------------------------------------------------------------- */
typedef struct {
    int    saved_out_fd;   /* dup() of fd 1, or -1 when capture is inactive */
    int    saved_err_fd;   /* dup() of fd 2 */
    FILE*  out_tmp;        /* anonymous temp file receiving redirected stdout */
    FILE*  err_tmp;        /* ... stderr */
    char*  out_buf;        /* captured bytes, filled on _end  */
    size_t out_len;
    char*  err_buf;
    size_t err_len;
} PipeCapture;

#ifndef _WIN32
/* Read a temp file written via a redirected descriptor into a fresh heap
 * buffer (NUL-terminated).  Returns NULL and *len==0 on empty/failure. */
static char* pipe_slurp_tmp(FILE* f, size_t* len) {
    *len = 0;
    if (!f) return NULL;
    if (fseek(f, 0, SEEK_END) != 0) return NULL;
    long sz = ftell(f);
    if (sz <= 0) { rewind(f); return NULL; }
    rewind(f);
    char* buf = malloc((size_t)sz + 1);
    if (!buf) return NULL;
    size_t got = fread(buf, 1, (size_t)sz, f);
    buf[got] = '\0';
    *len = got;
    return buf;
}
#endif

static void pipe_capture_begin(PipeCapture* c) {
    memset(c, 0, sizeof(*c));
    c->saved_out_fd = -1;
    c->saved_err_fd = -1;
#ifndef _WIN32
    fflush(stdout);
    fflush(stderr);
    c->out_tmp = tmpfile();
    c->err_tmp = tmpfile();
    if (!c->out_tmp || !c->err_tmp) goto fail;
    /* Save the real descriptors, then point fd 1/2 at the temp files. */
    c->saved_out_fd = dup(fileno(stdout));
    c->saved_err_fd = dup(fileno(stderr));
    if (c->saved_out_fd < 0 || c->saved_err_fd < 0) goto fail;
    if (dup2(fileno(c->out_tmp), fileno(stdout)) < 0 ||
        dup2(fileno(c->err_tmp), fileno(stderr)) < 0) goto fail;
    return;
fail:
    /* Could not redirect: leave the streams alone so output still reaches the
     * terminal, and mark the capture inactive. */
    if (c->saved_out_fd >= 0) { close(c->saved_out_fd); c->saved_out_fd = -1; }
    if (c->saved_err_fd >= 0) { close(c->saved_err_fd); c->saved_err_fd = -1; }
    if (c->out_tmp) { fclose(c->out_tmp); c->out_tmp = NULL; }
    if (c->err_tmp) { fclose(c->err_tmp); c->err_tmp = NULL; }
#endif
}

#ifndef _WIN32
/* Send captured messages, one "message" line each. A message starts on a line
 * that begins with a non-blank character and contains "::" (the Head::tag
 * form every funnelled message has); any other line continues the one before
 * it, which is how mth_message_cont writes a multi-line message. */
static void pipe_emit_messages(int id, const char* text) {
    const char* msg_start = NULL;
    const char* p = text;
    while (*p) {
        const char* eol = strchr(p, '\n');
        const char* next = eol ? eol + 1 : p + strlen(p);
        size_t len = (size_t)((eol ? eol : next) - p);
        bool starts = len > 0 && !isspace((unsigned char)*p);
        if (starts) {
            const char* dc = strstr(p, "::");
            starts = dc && dc < p + len;
        }
        if (starts && msg_start) {
            size_t mlen = (size_t)(p - msg_start);
            while (mlen > 0 && (msg_start[mlen - 1] == '\n' || msg_start[mlen - 1] == '\r')) mlen--;
            if (mlen) pipe_emit_text(id, "message", "text", msg_start, mlen);
            msg_start = NULL;
        }
        if (!msg_start && len > 0) msg_start = p;
        p = next;
    }
    if (msg_start) {
        size_t mlen = strlen(msg_start);
        while (mlen > 0 && (msg_start[mlen - 1] == '\n' || msg_start[mlen - 1] == '\r')) mlen--;
        if (mlen) pipe_emit_text(id, "message", "text", msg_start, mlen);
    }
}

#endif /* !_WIN32 */

/* Restore the descriptors and send what was captured: Print text, then
 * messages. */
static void pipe_capture_end(PipeCapture* c, int id) {
#ifndef _WIN32
    if (c->saved_out_fd < 0) return;   /* capture was inactive */
    /* Flush the redirected FILE* buffers down to the temp files, then put the
     * real descriptors back before reading/emitting (so emit writes to the
     * real stdout). */
    fflush(stdout);
    fflush(stderr);
    dup2(c->saved_out_fd, fileno(stdout));
    dup2(c->saved_err_fd, fileno(stderr));
    close(c->saved_out_fd);
    close(c->saved_err_fd);

    c->out_buf = pipe_slurp_tmp(c->out_tmp, &c->out_len);
    c->err_buf = pipe_slurp_tmp(c->err_tmp, &c->err_len);
    fclose(c->out_tmp);   /* removes the anonymous temp file */
    fclose(c->err_tmp);

    if (c->out_buf && c->out_len > 0)
        pipe_emit_text(id, "stream", "text", c->out_buf, c->out_len);
    if (c->err_buf && c->err_len > 0)
        pipe_emit_messages(id, c->err_buf);
    free(c->out_buf);
    free(c->err_buf);
    memset(c, 0, sizeof(*c));
    c->saved_out_fd = -1;
    c->saved_err_fd = -1;
#else
    (void)c; (void)id;
#endif
}

static void pipe_emit_parse_error(const char* input, int id) {
    /* Echo the exact received input in the error so a stray/invisible
     * character or bracket mismatch in the caller's text is diagnosable
     * rather than an opaque "Parse error". */
    size_t esc_cap = strlen(input) * 6 + 8;
    char* esc = malloc(esc_cap);
    char* buf = NULL;
    if (esc) {
        json_escape(input, esc, esc_cap);
        size_t bcap = esc_cap + 128;
        buf = malloc(bcap);
        if (buf)
            snprintf(buf, bcap,
                "{\"id\":%d,\"type\":\"error\",\"message\":\"Parse error: %s\"}",
                id, esc);
    }
    if (buf) {
        pipe_emit(buf);
    } else {
        char sbuf[128];
        snprintf(sbuf, sizeof(sbuf),
            "{\"id\":%d,\"type\":\"error\",\"message\":\"Parse error\"}", id);
        pipe_emit(sbuf);
    }
    free(esc);
    free(buf);
}

/* Evaluate ONE statement and send its output. Borrows `parsed`; sends no
 * "done" -- pipe_process_input does, once per request.
 *
 * In a notebook cell (`cell`), Print text and messages are captured and sent
 * first, and no result is sent when `show_result` is false (the statement
 * ended in ';') or the value is Null. A plain request keeps the original
 * behaviour: nothing captured, and Null sent as a payload like any value. */
static void pipe_eval_statement(Expr* parsed, int id, bool show_result, bool cell) {
    /* `?sym` / Information[sym] yields the raw docstring as a String, and a
     * usage message is not an expression: it must not be quoted, InputForm
     * escaped, or handed to the front end's math renderer. Captured from the
     * *input* head, exactly as the interactive REPL does above, so only help
     * queries take this path and an ordinary string result still comes back
     * quoted. */
    bool info_query = pipe_is_info_query(parsed);
    /* Copied into a buffer rather than kept as a borrowed pointer: `?"name"`
     * gives a STRING whose storage dies with the tree. */
    char info_sym[128];
    info_sym[0] = '\0';
    if (info_query) {
        const char* isname = pipe_info_symbol(parsed);
        if (isname) {
            strncpy(info_sym, isname, sizeof(info_sym) - 1);
            info_sym[sizeof(info_sym) - 1] = '\0';
        }
    }

    /* Session history, so `%` / `%%` / `%n` / In[n] / Out[n] work in a notebook
     * as they do in the interactive REPL.
     *
     * ONE LINE PER STATEMENT, which is how the REPL numbers them too: `%` then
     * means "the previous result" whether that came from the cell above or from
     * the line above inside this cell. Numbered before evaluating, because `%`
     * in this statement is Out[$Line - 1].
     *
     * Cells only, deliberately. A plain request is a one-shot evaluation from a
     * batch tool (the site generator, the audit sweeps, the book tools), which
     * has no numbered history to refer back to -- and holding every result alive
     * in Out[n] would grow a sweep of 10^7-element arrays without bound. `cell`
     * already draws exactly this line: notebook semantics on one side, one
     * expression in and one result out on the other. */
    int line = 0;
    if (cell) {
        line = ++g_pipe_line;
        mth_session_set_line(line);
        mth_session_store_in(line, parsed);
        pipe_emit_line(id, line);
    }

    PipeCapture cap;
    if (cell) pipe_capture_begin(&cap);
    Expr* evaluated = evaluate(parsed);
    mth_block_graveyard_flush();   /* top-level boundary: nothing in flight */
    if (cell) pipe_capture_end(&cap, id);

    /* Out[line] BEFORE every early return below: a result that is suppressed by
     * ';', is Null, or is a picture rather than an expression is still the
     * session's previous result, and `%` must find it. */
    if (cell) mth_session_store_out(line, evaluated);

    if (!evaluated) return;

    /* Mathematica shows no Out[] for a statement ending in ';' or for a Null
     * value -- `x = 5;`, `Print["hi"]` -- so in a cell neither sends a result. */
    if (cell && (!show_result
        || (evaluated->type == EXPR_SYMBOL && evaluated->data.symbol.name == SYM_Null))) {
        expr_free(evaluated);
        return;
    }

    /* `?Pat*` evaluates to the LIST of matching names. Emit it as its own
     * message so the notebook can lay it out as a grid; as a plain expression
     * it would be a single long braced line run through the math renderer. */
    if (info_query && evaluated->type == EXPR_FUNCTION
        && evaluated->data.function.head
        && evaluated->data.function.head->type == EXPR_SYMBOL
        && evaluated->data.function.head->data.symbol.name == SYM_List) {
        size_t n = evaluated->data.function.arg_count;
        size_t cap_n = 64;
        for (size_t i = 0; i < n; i++) {
            Expr* e = evaluated->data.function.args[i];
            if (e->type == EXPR_STRING) cap_n += strlen(e->data.string) * 6 + 8;
        }
        char* buf = malloc(cap_n);
        if (buf) {
            int off = snprintf(buf, cap_n, "{\"id\":%d,\"type\":\"names\",\"payload\":[", id);
            bool first = true;
            for (size_t i = 0; i < n && off > 0 && (size_t)off < cap_n; i++) {
                Expr* e = evaluated->data.function.args[i];
                if (e->type != EXPR_STRING) continue;
                size_t ecap = strlen(e->data.string) * 6 + 8;
                char* esc = malloc(ecap);
                if (!esc) break;
                json_escape(e->data.string, esc, ecap);
                off += snprintf(buf + off, cap_n - (size_t)off, "%s\"%s\"",
                                first ? "" : ",", esc);
                free(esc);
                first = false;
            }
            if (off > 0 && (size_t)off < cap_n) snprintf(buf + off, cap_n - (size_t)off, "]}");
            pipe_emit(buf);
            free(buf);
        }
        expr_free(evaluated);
        return;
    }

    if (info_query && evaluated->type == EXPR_STRING) {
        const char* doc = evaluated->data.string;
        size_t dcap = strlen(doc) * 6 + 8;
        char* esc = malloc(dcap);
        if (esc) {
            json_escape(doc, esc, dcap);
            size_t bcap = dcap + 256;
            char* buf = malloc(bcap);
            if (buf) {
                if (info_sym[0])
                    snprintf(buf, bcap,
                        "{\"id\":%d,\"type\":\"usage\",\"payload\":\"%s\","
                        "\"symbol\":\"%s\"}", id, esc, info_sym);
                else
                    snprintf(buf, bcap,
                        "{\"id\":%d,\"type\":\"usage\",\"payload\":\"%s\"}", id, esc);
                pipe_emit(buf);
                free(buf);
            }
            free(esc);
        }
        expr_free(evaluated);
        return;
    }

    /* Image[...] / Image3D[...] → RGBA for the notebook to draw on a canvas.
     *
     * Before this, an image came back as type "expr" and the notebook printed its expression -- which
     * after the storage canonicalisation reads `Image[NDArray[{{...}}], "Real"]`. Correct, and not a
     * picture. The front end had no image output kind at all. */
    if (evaluated->type == EXPR_FUNCTION
        && evaluated->data.function.head
        && evaluated->data.function.head->type == EXPR_SYMBOL) {
        const char* ih = evaluated->data.function.head->data.symbol.name;
        if (ih && (ih == SYM_Image || ih == SYM_Image3D)) {
            char* ijson = image_to_json(evaluated);
            if (ijson) {
                expr_free(evaluated);
                size_t jl = strlen(ijson) + 64;
                char* jline = malloc(jl);
                if (jline) {
                    snprintf(jline, jl, "{\"id\":%d,\"type\":\"image\",\"payload\":%s}",
                             id, ijson);
                    pipe_emit(jline);
                    free(jline);
                }
                free(ijson);
                return;
            }
            /* Not a well-formed image after all: fall through and print it as text. */
        }
    }

    /* Graphics[...] / Graphics3D[...] → Plotly JSON for the notebook. */
    if (evaluated->type == EXPR_FUNCTION
        && evaluated->data.function.head
        && evaluated->data.function.head->type == EXPR_SYMBOL) {
        const char* head_sym = evaluated->data.function.head->data.symbol.name;
        char* plotly = NULL;
        if (head_sym == SYM_Graphics)
            plotly = graphics_to_plotly_json(evaluated);
        else if (head_sym == SYM_Graphics3D)
            plotly = graphics3d_to_plotly_json(evaluated);
        if (plotly) {
            expr_free(evaluated);
            size_t json_len = strlen(plotly) + 64;
            char* json_line = malloc(json_len);
            if (json_line) {
                snprintf(json_line, json_len, "{\"id\":%d,\"type\":\"plot\",\"payload\":%s}",
                         id, plotly);
                pipe_emit(json_line);
                free(json_line);
            }
            free(plotly);
            return;
        }
        /* plotly == NULL (e.g. empty Graphics3D): nothing to draw. */
        if (head_sym == SYM_Graphics || head_sym == SYM_Graphics3D) {
            expr_free(evaluated);
            return;
        }
    }

    char* result_str = expr_to_string_display(evaluated);
    char* latex_raw  = expr_to_latex(evaluated);   /* must be before expr_free */
    LatexParts sparts;
    char* slatex_raw = expr_to_latex_selectable(evaluated, &sparts);
    expr_free(evaluated);

    if (!result_str) {
        free(latex_raw);
        free(slatex_raw);
        latex_parts_free(&sparts);
        char buf[128];
        snprintf(buf, sizeof(buf),
                 "{\"id\":%d,\"type\":\"error\",\"message\":\"Out of memory\"}", id);
        pipe_emit(buf);
        return;
    }

    size_t escaped_len = strlen(result_str) * 6 + 4;
    char* escaped = malloc(escaped_len);
    if (!escaped) {
        free(result_str);
        free(latex_raw);
        free(slatex_raw);
        latex_parts_free(&sparts);
        char buf[128];
        snprintf(buf, sizeof(buf),
                 "{\"id\":%d,\"type\":\"error\",\"message\":\"Out of memory\"}", id);
        pipe_emit(buf);
        return;
    }
    json_escape(result_str, escaped, escaped_len);
    free(result_str);

    /* latex_raw was produced above (before expr_free) — now escape it */
    char* latex_esc  = pipe_json_field("latex", latex_raw);
    char* slatex_esc = pipe_json_field("slatex", slatex_raw);
    char* sparts_esc = slatex_esc
        ? pipe_json_str_array("sparts", sparts.text, sparts.count) : NULL;
    free(latex_raw);
    free(slatex_raw);
    latex_parts_free(&sparts);

    size_t line_len = escaped_len + (latex_esc ? strlen(latex_esc) : 0)
                    + (slatex_esc ? strlen(slatex_esc) : 0)
                    + (sparts_esc ? strlen(sparts_esc) : 0) + 128;
    char* json_line = malloc(line_len);
    if (json_line) {
        snprintf(json_line, line_len,
                 "{\"id\":%d,\"type\":\"expr\",\"payload\":\"%s\"%s%s%s}",
                 id, escaped, latex_esc ? latex_esc : "", slatex_esc ? slatex_esc : "",
                 sparts_esc ? sparts_esc : "");
        pipe_emit(json_line);
        free(json_line);
    }
    free(escaped);
    free(latex_esc);
    free(slatex_esc);
    free(sparts_esc);
}

/* Evaluate one request, then send "done".
 *
 * A plain request is ONE expression, read with parse_expression as it always
 * was. A notebook cell (`cell`) may hold several statements: it is split with
 * parse_next_expression -- the statement reader
 * `-file` scripts use -- so statements on separate lines are separate
 * statements, as in a Mathematica input cell. parse_expression, used before,
 * reads exactly ONE expression and rejected `a = 1` NEWLINE `b = 2` as a parse
 * error, because the parser ends a statement at a top-level newline and the
 * second line was left over.
 *
 * The whole cell is parsed before anything is evaluated, so a syntax error
 * anywhere evaluates nothing (Mathematica checks a cell's syntax first too).
 * A ';' that ended a statement suppresses that statement's result; the reader
 * consumes it, so it is recognised as the character just before the cursor. */
static void pipe_process_input(const char* input, int id, bool cell) {
    if (!cell) {
        Expr* parsed = parse_expression(input);
        if (!parsed) {
            pipe_emit_parse_error(input, id);
        } else {
            pipe_eval_statement(parsed, id, true, false);
            expr_free(parsed);
        }
        pipe_emit_done(id);
        return;
    }

    size_t n = 0, cap = 8;
    Expr** stmts = malloc(cap * sizeof(Expr*));
    bool* shown = malloc(cap * sizeof(bool));
    bool failed = (stmts == NULL || shown == NULL);

    const char* p = input;
    while (!failed) {
        const char* start = p;
        Expr* e = parse_next_expression(&p);
        if (!e) {
            /* NULL is both end of input and a syntax error: anything but blanks
             * and comments left over is the error. */
            if (*mth_skip_blanks_and_comments(start) != '\0') failed = true;
            break;
        }
        if (n == cap) {
            size_t ncap = cap * 2;
            Expr** ns = realloc(stmts, ncap * sizeof(Expr*));
            if (ns) stmts = ns;
            bool* nb = realloc(shown, ncap * sizeof(bool));
            if (nb) shown = nb;
            if (!ns || !nb) { expr_free(e); failed = true; break; }
            cap = ncap;
        }
        stmts[n] = e;
        shown[n] = !(p > start && p[-1] == ';');
        n++;
    }

    if (failed) {
        pipe_emit_parse_error(input, id);
    } else {
        for (size_t i = 0; i < n; i++) pipe_eval_statement(stmts[i], id, shown[i], true);
    }
    for (size_t i = 0; i < n; i++) expr_free(stmts[i]);
    free(stmts);
    free(shown);
    pipe_emit_done(id);
}

/* Read one line of any length from `in` into a fresh malloc'd buffer, without
 * its line terminator. NULL at end of input. A request line used to be read
 * into a 10 KB stack buffer, and fgets silently split anything longer: the
 * first 10 KB went out as a truncated cell, the rest as garbage requests. */
static char* pipe_read_line(FILE* in) {
    size_t cap = 4096, len = 0;
    char* buf = malloc(cap);
    if (!buf) return NULL;
    int c;
    while ((c = getc(in)) != EOF && c != '\n') {
        if (len + 1 >= cap) {
            char* nb = realloc(buf, cap * 2);
            if (!nb) { free(buf); return NULL; }
            buf = nb;
            cap *= 2;
        }
        buf[len++] = (char)c;
    }
    if (c == EOF && len == 0) { free(buf); return NULL; }
    while (len > 0 && buf[len - 1] == '\r') len--;
    buf[len] = '\0';
    return buf;
}

/* Emit the structural-selection spans as a JSON array of [start, end] integer
 * pairs (byte offsets into the request's expr). Pure integers → no escaping. */
static void pipe_emit_spans(int id, const MthSpanSink* k) {
    size_t cap = 64 + (k ? k->n : 0) * 24;
    char* buf = malloc(cap);
    if (!buf) return;
    int off = snprintf(buf, cap, "{\"id\":%d,\"type\":\"spans\",\"payload\":[", id);
    if (k) {
        for (size_t i = 0; i < k->n && off > 0 && (size_t)off < cap; i++)
            off += snprintf(buf + off, cap - (size_t)off, "%s[%d,%d]",
                            i ? "," : "", k->v[i].start, k->v[i].end);
    }
    if (off > 0 && (size_t)off < cap) snprintf(buf + off, cap - (size_t)off, "]}");
    pipe_emit(buf);
    free(buf);
}

/* Handle a {"id":N,"expr":"...","spans":true} request: parse WITHOUT evaluating
 * and return the subexpression spans the notebook uses for structural selection,
 * then the terminating done. */
static void pipe_process_spans(const char* input, int id) {
    MthSpanSink* k = mth_parse_spans(input);
    pipe_emit_spans(id, k);
    mth_span_sink_free(k);
    pipe_emit_done(id);
}

static void pipe_mode_loop(void) {
    g_pipe_out = stdout;
    char* line;
    while ((line = pipe_read_line(stdin)) != NULL) {
        if (line[0] == '\0') { free(line); continue; }

        /* Control messages carry no "expr", which keeps an expression that
         * merely mentions "ping" or "quit" from being mistaken for one. */
        bool has_expr = strstr(line, "\"expr\"") != NULL;
        if (!has_expr && strstr(line, "\"ping\"")) {
            pipe_emit("{\"type\":\"pong\"}");
            free(line);
            continue;
        }
        if (!has_expr && strstr(line, "\"quit\"")) {
            free(line);
            fflush(stdout);
            break;
        }

        int id = 0;
        char* expr = NULL;
        if (json_get_int(line, "id", &id))
            expr = json_get_string_dup(line, "expr");
        if (expr) {
            /* A "spans" request parses for structural selection and never
             * evaluates; everything else is an ordinary evaluation. */
            if (json_get_true(line, "spans"))
                pipe_process_spans(expr, id);
            else
                pipe_process_input(expr, id, json_get_true(line, "cell"));
        }
        free(expr);
        free(line);
    }
}

/* =====================================================================
 * Script mode:  Mathilda -file script.m
 *
 * Runs a file the way `wolframscript -file` does: every expression in the
 * file is parsed and evaluated in order and nothing is echoed, so the
 * script's output is exactly what it Print[]s. This is the same evaluation
 * the Get["file"] builtin performs, but reached without a live session — the
 * loop below is spelled out here rather than delegating to
 * mathilda_run_file() because a script runner owes the caller a real
 * diagnostic and a nonzero exit status when the file has a syntax error,
 * where Get[] simply stops reading.
 * ===================================================================*/

/* Read all of `path` into a freshly malloc'd, NUL-terminated buffer.
 * Returns NULL if the file cannot be opened or read (caller reports). */
static char* read_whole_file(const char* path) {
    FILE* fp = fopen(path, "rb");
    if (!fp) return NULL;

    if (fseek(fp, 0, SEEK_END) != 0) { fclose(fp); return NULL; }
    long fsize = ftell(fp);
    if (fsize < 0) { fclose(fp); return NULL; }
    rewind(fp);

    char* buffer = malloc((size_t)fsize + 1);
    if (!buffer) { fclose(fp); return NULL; }

    size_t read_len = fread(buffer, 1, (size_t)fsize, fp);
    buffer[read_len] = '\0';          /* short read (text mode/CRLF) is fine */
    fclose(fp);
    return buffer;
}

/* Report an error at `pos` within the script text `buf`, in the
 * "file:line: message" form editors and CI logs already know how to read,
 * followed by the offending source line. */
static void report_error_at(const char* path, const char* buf, const char* pos,
                            const char* message) {
    int line = 1;
    const char* line_start = buf;
    for (const char* p = buf; p < pos && *p; p++) {
        if (*p == '\n') { line++; line_start = p + 1; }
    }
    const char* line_end = line_start;
    while (*line_end && *line_end != '\n' && *line_end != '\r') line_end++;

    fprintf(stderr, "%s:%d: %s\n", path, line, message);
    fprintf(stderr, "  %.*s\n", (int)(line_end - line_start), line_start);
}

/* Locate an unterminated string literal or (* ... *) comment in `buf`,
 * returning a pointer to its opener (and naming the construct in *what), or
 * NULL if the text is well formed.
 *
 * The lexer treats both as running to end of file, so in a script a single
 * stray `(*` would silently swallow every statement after it: the run would
 * exit 0 having quietly done half the work. Checking up front turns that into
 * a diagnostic before anything is evaluated. The scan alternates between the
 * two constructs deliberately — a `"` inside a comment and a `(*` inside a
 * string are both just text, and whichever opens first consumes the other. */
static const char* find_unterminated(const char* buf, const char** what) {
    const char* p = buf;
    while (*p) {
        if (*p == '"') {
            const char* open = p++;
            while (*p && *p != '"') p += (*p == '\\' && p[1]) ? 2 : 1;
            if (!*p) { *what = "string"; return open; }
            p++;
        } else if (p[0] == '(' && p[1] == '*') {
            const char* open = p;
            int depth = 1;
            p += 2;
            while (*p && depth > 0) {
                if (p[0] == '(' && p[1] == '*') { depth++; p += 2; }
                else if (p[0] == '*' && p[1] == ')') { depth--; p += 2; }
                else { p++; }
            }
            if (depth > 0) { *what = "comment"; return open; }
        } else {
            p++;
        }
    }
    return NULL;
}

/* Evaluate every expression in `path`. Returns the process exit status:
 * 0 on success, 1 if the file could not be read or contained a syntax
 * error. */
static int run_script_file(const char* path) {
    char* buffer = read_whole_file(path);
    if (!buffer) {
        fprintf(stderr, "Mathilda: cannot open file: %s\n", path);
        return 1;
    }

    /* Up-front lexical check: an unterminated string or comment would make
     * the parser swallow the rest of the file as if it were not there. */
    const char* what = NULL;
    const char* opener = find_unterminated(buffer, &what);
    if (opener) {
        char message[64];
        snprintf(message, sizeof(message), "unterminated %s", what);
        report_error_at(path, buffer, opener, message);
        free(buffer);
        return 1;
    }

    int status = 0;
    Expr* prev_input = mth_input_file_enter(path);   /* $InputFileName */
    const char* ptr = buffer;
    while (*ptr != '\0') {
        const char* stmt_start = ptr;
        Expr* parsed = parse_next_expression(&ptr);
        if (!parsed) {
            /* parse_next_expression returns NULL both at end of input and on
             * a syntax error. Trailing whitespace or comments are a normal
             * end of file; anything else left unconsumed is a real error,
             * and it starts at the first significant character — not at
             * stmt_start, which is still sitting on the newline that ended
             * the previous statement. */
            const char* err = mth_skip_blanks_and_comments(stmt_start);
            if (*err != '\0') {
                report_error_at(path, buffer, err, "syntax error");
                status = 1;
            }
            break;
        }
        /* evaluate() borrows its argument and returns a new tree; both are
         * ours to free. The result is discarded: a script speaks via
         * Print[], not via echoed values. */
        (void)eval_take_top_abort();   /* clear any stale flag */
        Expr* evaluated = evaluate(parsed);
        mth_block_graveyard_flush();   /* top-level boundary: nothing in flight */
        expr_free(evaluated);
        expr_free(parsed);
        if (eval_take_top_abort()) {
            /* Abort[] reached the top level: the script ends here, with
             * $Aborted shown as the result. */
            printf("$Aborted\n");
            break;
        }
        if (mth_quit_requested()) {
            /* Quit[n] ends the script here, with exit status n. */
            status = mth_quit_status();
            break;
        }
    }
    mth_input_file_leave(prev_input);

    free(buffer);

    /* A script run is a session; give $Epilog its one evaluation, exactly as
     * the interactive loop does on exit. */
    repl_apply_epilog();
    fflush(stdout);
    return status;
}

/* Dialog[] input readers (see main). Each shows the prompt and returns one
 * malloc'd line, or NULL at end of input. */
#ifndef NO_READLINE
/* A terminal: the session's line reader (ctx). A Ctrl-C at the prompt clears
 * the line and asks again; accepted lines go into the reader's history. */
static char* dialog_read_term(const char* prompt, void* ctx) {
    TermReader* rd = (TermReader*)ctx;
    for (;;) {
        char* line = NULL;
        TermReadStatus st = term_read(rd, prompt, &line);
        if (st == TERM_READ_INTERRUPTED) continue;
        if (st == TERM_READ_EOF) return NULL;
        term_reader_history_add(rd, line);
        return line;
    }
}
#endif

static char* dialog_read_stdin(const char* prompt, void* ctx) {
    (void)ctx;
    fputs(prompt, stdout);
    fflush(stdout);
    char* line = pipe_read_line(stdin);
    /* A terminal echoes the typed line and its newline; a pipe does not, so
     * end the prompt line here and the Out line starts on its own, as at a
     * terminal. */
    if (line && !isatty(fileno(stdin))) { fputc('\n', stdout); fflush(stdout); }
    return line;
}

static void print_usage(FILE* out, const char* prog) {
    fprintf(out,
        "Mathilda " MATHILDA_VERSION_STRING " - a small, open source computer algebra system.\n"
        "\n"
        "Usage: %s [options] [file]\n"
        "\n"
        "  -file <path>    evaluate every expression in <path>, then exit\n"
        "  -h, --help      show this message and exit\n"
        "  -v, --version   print version information and exit\n"
        "\n"
        "A bare <file> argument is equivalent to -file <file>. With no file,\n"
        "Mathilda starts the interactive REPL when stdin is a terminal, and\n"
        "otherwise speaks the NDJSON pipe protocol on stdio.\n",
        prog);
}

int main(int argc, char** argv) {
    const char* prog   = (argc > 0 && argv[0]) ? argv[0] : "Mathilda";
    const char* script = NULL;

    for (int i = 1; i < argc; i++) {
        const char* arg = argv[i];
        if (strcmp(arg, "-file") == 0 || strcmp(arg, "--file") == 0) {
            if (i + 1 >= argc) {
                fprintf(stderr, "%s: %s requires a path\n", prog, arg);
                return 2;
            }
            script = argv[++i];
        } else if (strcmp(arg, "-h") == 0 || strcmp(arg, "--help") == 0) {
            print_usage(stdout, prog);
            return 0;
        } else if (strcmp(arg, "-v") == 0 || strcmp(arg, "--version") == 0) {
            printf("%s\n", mathilda_version());
            return 0;
        } else if (arg[0] == '-' && arg[1] != '\0') {
            fprintf(stderr, "%s: unknown option: %s\n", prog, arg);
            print_usage(stderr, prog);
            return 2;
        } else if (!script) {
            script = arg;               /* bare path: Mathilda script.m */
        } else {
            fprintf(stderr, "%s: unexpected argument: %s\n", prog, arg);
            return 2;
        }
    }

    /* Detect pipe mode: when stdin is not a terminal the frontend has
     * spawned us as a sidecar and we communicate via NDJSON over stdio.
     * A -file run is a script, not a session, so it takes precedence: the
     * script's own output must not be wrapped in the pipe protocol just
     * because it was launched from a shell script with redirected stdin.
     * The interactive readline REPL is preserved when stdin is a tty. */
    int pipe_mode = !script && !isatty(fileno(stdin));

    if (pipe_mode) {
        /* Disable libc's stdout buffer so every response line is delivered
         * to the pipe immediately rather than accumulating. */
        setvbuf(stdout, NULL, _IONBF, 0);
    } else if (script && !isatty(fileno(stdout))) {
        /* Redirected script output is block-buffered by default, which holds
         * back a long benchmark's progress until it exits. Line-buffer it so
         * `Mathilda -file bench.m | tee log` streams as it runs. */
        setvbuf(stdout, NULL, _IOLBF, BUFSIZ);
    }

    symtab_init();
    core_init();

    /* Load the internal bootstrap (init.m). Path resolution is independent of
     * the current working directory (see mathilda_load_module), so a relocated
     * or installed binary still finds its bundled src/internal tree. If it
     * cannot be located the loader prints a LoadModule::nofile diagnostic —
     * far better than the previous silent load of a non-functional kernel. */
    mathilda_load_module("init.m");

    /* Where Dialog[] reads its lines: the terminal (readline) when there is
     * one, a script's stdin under -file, and nowhere in pipe mode, where stdin
     * is the protocol -- a dialog there sees end of input at once. */
    MthDialogIO dio = { NULL, NULL, NULL, false };
    TermReader* rd = NULL;
    TermHistory* hist = NULL;
    bool interactive_rl = false;   /* the full interactive readline session */
    if (!pipe_mode) {
        dio.interactive = isatty(fileno(stdin)) != 0;
        dio.read = dialog_read_stdin;
#ifndef NO_READLINE
        if (!script || dio.interactive) {
            /* The interactive session gets every terminal feature; a -file
             * run from a terminal gets the plain reader its Dialog always
             * had (readline and in-memory history, nothing else). */
            if (!script) {
                /* Colour is decided once, here, for this session only. */
                term_color_init(true);
                if (term_color_stderr())
                    mth_message_set_style(term_style_on(TERM_STYLE_MSG_HEAD, true),
                                          term_style_off(true));
                /* The history file: $MATHILDA_HISTFILE (empty disables it)
                 * or ~/.mathilda_history, $MATHILDA_HISTSIZE entries. With
                 * no file the store still filters the session's history. */
                size_t hmax = term_history_size_from_env(getenv("MATHILDA_HISTSIZE"));
                char* hp = term_history_resolve_path(getenv("MATHILDA_HISTFILE"), getenv("HOME"));
                hist = term_history_open(hp, hmax);
                free(hp);
                if (!hist) hist = term_history_open(NULL, hmax);
            }
            /* The in-tree editor on a capable terminal; readline when it is
             * declined (MATHILDA_EDITOR=readline, TERM=dumb) or cannot open. */
            if (!script && term_editor_select(getenv("MATHILDA_EDITOR"), getenv("TERM"),
                                              dio.interactive, isatty(fileno(stdout)) != 0)) {
                TermEditorOptions eo = { hist, term_color_stdout() };
                rd = term_reader_open_editor(&eo);
            }
            if (!rd) {
                TermReadlineOptions ro = { !script, hist };
                rd = term_reader_open_readline(&ro);
            }
            interactive_rl = rd && !script;
        }
#endif
        if (!rd && !script) rd = term_reader_open_stdio();
#ifndef NO_READLINE
        if (rd && dio.interactive) { dio.read = dialog_read_term; dio.ctx = rd; }
#endif
        mth_dialog_set_io(&dio);
    }

    int rc = 0;
    if (script) {
        rc = run_script_file(script);
    } else if (pipe_mode) {
        pipe_mode_loop();
    } else {
        /* Ctrl-C handling only for the readline session: the fgets loop of a
         * NO_READLINE build keeps the default (Ctrl-C ends the process). */
        if (interactive_rl) term_interrupt_install();
        repl_loop(rd);
        if (interactive_rl) term_interrupt_uninstall();
        if (mth_quit_requested()) rc = mth_quit_status();
    }
    if (rd) {
        mth_dialog_set_io(NULL);
        term_reader_close(rd);
    }
    term_history_close(hist);
    /* Tear down the context subsystem so its $Context / $ContextPath strings and
     * any open package frames are freed rather than lingering as leaks at exit.
     * context_init() is called from core_init(); this is its symmetric partner,
     * which was previously never invoked (dead code). Nothing runs after this. */
    regex_cache_clear();     /* compiled string-pattern cache (leak hygiene) */
    context_shutdown();
    return rc;
}
