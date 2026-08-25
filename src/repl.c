/* fileno()/isatty() are POSIX, hidden by glibc under -std=c99; request them. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "expr.h"
#include "parse.h"
#include "print.h"
#include "part.h"
#include "eval.h"
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

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

#ifndef NO_READLINE
  #include <readline/readline.h>
  #include <readline/history.h>
#endif

#define MAX_INPUT_LEN 10240

/* Advance past whitespace and (* ... *) comments (nested comments allowed),
 * returning the first character that is neither — the terminating NUL if the
 * text holds nothing else. An unterminated comment is deliberately NOT
 * skipped: the returned pointer is its opening '(', because that is a genuine
 * error and that is where it starts. Script mode uses the position to point
 * at the offending token; is_blank_or_comment_only() only asks whether
 * anything is left. */
static const char* skip_blanks_and_comments(const char* s) {
    while (*s) {
        if (isspace((unsigned char)*s)) {
            s++;
        } else if (s[0] == '(' && s[1] == '*') {
            const char* open = s;
            int depth = 1;
            s += 2;
            while (*s && depth > 0) {
                if (s[0] == '(' && s[1] == '*') { depth++; s += 2; }
                else if (s[0] == '*' && s[1] == ')') { depth--; s += 2; }
                else { s++; }
            }
            if (depth > 0) return open;  /* unterminated comment is a real error */
        } else {
            return s;
        }
    }
    return s;
}

/* True if `s` consists only of whitespace and (* ... *) comments. Used to
 * distinguish a no-op line from a genuine parse failure, so the REPL doesn't
 * shout "Parse error" at a stray comment. */
static int is_blank_or_comment_only(const char* s) {
    return *skip_blanks_and_comments(s) == '\0';
}

/* Mathematica strips a top-level NumberForm from the value stored in Out[n]/%,
 * so `%` (and arithmetic on it) sees the underlying number rather than the
 * inert format wrapper -- while the displayed Out[n]= line still uses the
 * wrapper's formatting. The stripping is TOP-LEVEL only: a nested NumberForm,
 * or one bound to a variable via Set, is preserved. Returns a BORROWED pointer
 * into `e` (or `e` itself when there is nothing to strip). */
static Expr* out_value_unwrapped(Expr* e) {
    if (e && e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_NumberForm
        && e->data.function.arg_count >= 1)
        return e->data.function.args[0];
    return e;
}

void process_input(const char* input, int line_number) {
    if (strlen(input) == 0) return;
    if (is_blank_or_comment_only(input)) return;

    // Update $Line
    Expr* line_sym = expr_new_symbol(SYM_DollarLine);
    Expr* line_val = expr_new_integer(line_number);
    symtab_add_own_value("$Line", line_sym, line_val);
    expr_free(line_sym);
    expr_free(line_val);

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

    /* Store In[line_number] = parsed BEFORE running $Pre/$Post so that
     * In[n] reflects what the user typed, not whatever a hook did. */
    Expr* in_sym = expr_new_symbol(SYM_In);
    Expr* in_arg = expr_new_integer(line_number);
    Expr* in_args[] = {in_arg};
    Expr* in_pattern = expr_new_function(in_sym, in_args, 1);
    symtab_add_down_value("In", in_pattern, parsed);
    expr_free(in_pattern);

    /* $Pre: applied to the parsed expression. Consumes our reference
     * to `parsed`; we treat the result as the new input to evaluate. */
    Expr* pre_input = repl_apply_pre(expr_copy(parsed));
    Expr* evaluated = evaluate(pre_input);
    expr_free(pre_input);

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
     * Mathematica's documented ordering). */
    Expr* out_sym = expr_new_symbol(SYM_Out);
    Expr* out_arg = expr_new_integer(line_number);
    Expr* out_args[] = {out_arg};
    Expr* out_pattern = expr_new_function(out_sym, out_args, 1);
    /* Store the unwrapped value (add_down_value copies its argument); the full
     * `evaluated` is kept for the formatted display below. */
    symtab_add_down_value("Out", out_pattern, out_value_unwrapped(evaluated));
    expr_free(out_pattern);

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

    printf("Out[%d]= ", line_number);
    if (is_info_query) fputs(to_print->data.string, stdout);
    else               expr_print(to_print);
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

#ifndef NO_READLINE

/* True when `s` is a syntactically complete Mathilda input: no unterminated
 * "..." string or (* ... *) comment, and every '(', '[' or '{' has a matching
 * closer. Newlines are whitespace to the lexer, so a balanced multi-line
 * buffer parses exactly as its single-line spelling. An *over*-closed buffer
 * (a stray ')') reports complete on purpose: appending text cannot repair it,
 * so it should submit now and let the parser flag the error rather than trap
 * the user on an endlessly growing line. The empty string is complete —
 * submitting it is a harmless no-op the caller skips. This is the predicate
 * behind the smart Return key: complete -> evaluate, incomplete -> open a
 * fresh continuation line.
 *
 * The string/comment lexing here mirrors find_unterminated() exactly (a '"'
 * inside a comment and a "(*" inside a string are both just text, and a
 * backslash escapes the next character inside a string) so the completeness
 * check never disagrees with the real lexer on where a token ends. */
static int mth_input_complete(const char* s) {
    int depth = 0;      /* net (), [], {} nesting outside strings/comments   */
    int comment = 0;    /* (* ... *) nesting depth                           */
    int in_string = 0;  /* inside a "..." literal                            */
    for (const char* p = s; *p; ) {
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

/* Return key handler: evaluate the buffer once it forms a complete
 * expression, otherwise open a new line so the user can keep typing. This is
 * the terminal-native substitute for a notebook's Shift+Enter — it needs no
 * enhanced-keyboard protocol and works on every terminal. */
static int mth_smart_return(int count, int key) {
    if (mth_input_complete(rl_line_buffer ? rl_line_buffer : ""))
        return rl_newline(count, key);   /* accept the whole buffer */
    rl_insert_text("\n");                /* still open: continue editing */
    return 0;
}

/* Esc-Return (and Alt/Meta-Return, where the terminal sends it) forces the
 * buffer to be evaluated even while mth_input_complete() still considers it
 * open. The escape hatch for a genuine syntax error the user wants to see,
 * so an unbalanced buffer can never trap them on a growing line. */
static int mth_force_return(int count, int key) {
    return rl_newline(count, key);
}

/* Install the completeness-driven Return bindings. Called once, when the
 * interactive loop starts. Bracketed paste is (re)enabled so a pasted block
 * with embedded newlines is inserted verbatim rather than firing the Return
 * handler on every line and submitting a fragment mid-paste. */
static void mth_setup_readline(void) {
    rl_variable_bind("enable-bracketed-paste", "on");
    rl_bind_key('\r', mth_smart_return);        /* RET / Ctrl-M */
    rl_bind_key('\n', mth_smart_return);        /* LFD / Ctrl-J */
    /* Esc-Return / Meta-Return force-submits an expression the completeness
     * check still considers open. rl_bind_keyseq() is a GNU Readline entry
     * point absent from Apple's libedit shim; RL_STATE_INITIALIZED is defined
     * only by GNU Readline, so it doubles as the "real readline" probe. */
#ifdef RL_STATE_INITIALIZED
    rl_bind_keyseq("\\e\\r", mth_force_return); /* Esc then Return */
    rl_bind_keyseq("\\e\\n", mth_force_return);
#else
    (void)mth_force_return;
#endif
}

void repl_loop() {
    printf("\nMathilda " MATHILDA_VERSION_STRING " - A small, open source computer algebra system.\n\n");
    printf("This program is free, open source software and comes with ABSOLUTELY NO WARRANTY.\n\n");
    printf("Press Return to evaluate. An open bracket, string or comment continues\n");
    printf("on the next line; press Esc then Return to force evaluation.\n");
    printf("Exit by evaluating Quit[] or CONTROL-C.\n\n");

    mth_setup_readline();

    int line_number = 1;
    char prompt[64];

    while (1) {
        snprintf(prompt, sizeof(prompt), "In[%d]:= ", line_number);

        /* With the smart Return binding a single readline() call returns the
         * whole (possibly multi-line) expression, so no accumulation buffer
         * is needed and readline owns the allocation. */
        char* line = readline(prompt);
        if (!line) {
            printf("\n");
            /* EOF: run $Epilog before tearing down. */
            repl_apply_epilog();
            break;
        }

        if (strlen(line) == 0) {
            free(line);
            continue;
        }

        add_history(line);

        if (strcmp(line, "Quit[]") == 0) {
            /* User-requested shutdown: run $Epilog first. */
            repl_apply_epilog();
            free(line);
            break;
        }

        process_input(line, line_number);
        line_number++;
        free(line);
    }

    printf("\n");
}
#else
/* Fallback interactive loop when readline is not available (e.g. Windows).
 * Uses fgets; no history or line-editing. Pipe mode bypasses this entirely. */
void repl_loop(void) {
    printf("\nMathilda " MATHILDA_VERSION_STRING " - A small, open source computer algebra system.\n\n");
    printf("Exit by evaluating Quit[] or pressing Ctrl+Z (Windows) / Ctrl+D (Unix).\n\n");

    char line[MAX_INPUT_LEN];
    int line_number = 1;

    while (1) {
        printf("In[%d]:= ", line_number);
        fflush(stdout);
        if (!fgets(line, sizeof(line), stdin)) {
            printf("\n");
            repl_apply_epilog();
            break;
        }
        /* Strip trailing newline / carriage-return. */
        size_t len = strlen(line);
        while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r'))
            line[--len] = '\0';
        if (len == 0) continue;
        if (strcmp(line, "Quit[]") == 0) {
            repl_apply_epilog();
            break;
        }
        process_input(line, line_number++);
    }
    printf("\n");
}
#endif

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
 * Request  (one line on stdin):
 *   {"id": N, "expr": "1+1"}    -- evaluate expression
 *   {"type": "ping"}             -- readiness probe
 *   {"type": "quit"}             -- graceful shutdown
 *
 * Response (one JSON object per line on stdout):
 *   {"id": N, "type": "expr",  "payload": "2"}
 *   {"id": N, "type": "error", "message": "Parse error"}
 *   {"id": N, "type": "done"}
 *   {"type": "pong"}
 *
 * stdout is set to unbuffered at startup so every response line is
 * delivered to the pipe immediately.
 * ===================================================================*/

static void pipe_emit(const char* line) {
    puts(line);
    fflush(stdout);
}

static int json_get_string(const char* json, const char* key,
                           char* buf, size_t buflen) {
    char search[256];
    snprintf(search, sizeof(search), "\"%s\"", key);
    const char* p = strstr(json, search);
    if (!p) return 0;
    p += strlen(search);
    while (*p == ' ' || *p == '\t' || *p == ':') p++;
    if (*p != '"') return 0;
    p++;
    size_t i = 0;
    while (*p && *p != '"' && i + 1 < buflen) {
        if (*p == '\\' && *(p + 1)) {
            p++;
            switch (*p) {
                case '"':  buf[i++] = '"';  break;
                case '\\': buf[i++] = '\\'; break;
                case '/':  buf[i++] = '/';  break;
                case 'n':  buf[i++] = '\n'; break;
                case 'r':  buf[i++] = '\r'; break;
                case 't':  buf[i++] = '\t'; break;
                default:   buf[i++] = *p;   break;
            }
        } else {
            buf[i++] = (char)*p;
        }
        p++;
    }
    buf[i] = '\0';
    return 1;
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

static void pipe_process_input(const char* input, int id) {
    Expr* parsed = parse_expression(input);
    if (!parsed) {
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
        char dbuf[64];
        snprintf(dbuf, sizeof(dbuf), "{\"id\":%d,\"type\":\"done\",\"memory\":%llu}",
                 id, (unsigned long long)pipe_memory_bytes());
        pipe_emit(dbuf);
        return;
    }

    /* `?sym` / Information[sym] yields the raw docstring as a String, and a
     * usage message is not an expression: it must not be quoted, InputForm
     * escaped, or handed to the front end's math renderer. Captured from the
     * *input* head before `parsed` is freed, exactly as the interactive REPL
     * does above, so only help queries take this path and an ordinary string
     * result still comes back quoted. */
    bool info_query = pipe_is_info_query(parsed);
    /* Captured BEFORE evaluate/expr_free, into a buffer rather than as a borrowed pointer.
     *
     * The first version read it after the free, which is undefined behaviour that happened to yield
     * nothing -- the notebook received no symbol and could not offer a documentation link. A symbol
     * name is interned and would have survived, but `?"name"` gives a STRING whose storage dies with
     * the tree, so the copy covers both. */
    char info_sym[128];
    info_sym[0] = '\0';
    if (info_query) {
        const char* isname = pipe_info_symbol(parsed);
        if (isname) {
            strncpy(info_sym, isname, sizeof(info_sym) - 1);
            info_sym[sizeof(info_sym) - 1] = '\0';
        }
    }

    Expr* evaluated = evaluate(parsed);
    expr_free(parsed);

    if (!evaluated) {
        char buf[64];
        snprintf(buf, sizeof(buf), "{\"id\":%d,\"type\":\"done\",\"memory\":%llu}",
                 id, (unsigned long long)pipe_memory_bytes());
        pipe_emit(buf);
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
        size_t cap = 64;
        for (size_t i = 0; i < n; i++) {
            Expr* e = evaluated->data.function.args[i];
            if (e->type == EXPR_STRING) cap += strlen(e->data.string) * 6 + 8;
        }
        char* buf = malloc(cap);
        if (buf) {
            int off = snprintf(buf, cap, "{\"id\":%d,\"type\":\"names\",\"payload\":[", id);
            bool first = true;
            for (size_t i = 0; i < n && off > 0 && (size_t)off < cap; i++) {
                Expr* e = evaluated->data.function.args[i];
                if (e->type != EXPR_STRING) continue;
                size_t ecap = strlen(e->data.string) * 6 + 8;
                char* esc = malloc(ecap);
                if (!esc) break;
                json_escape(e->data.string, esc, ecap);
                off += snprintf(buf + off, cap - (size_t)off, "%s\"%s\"",
                                first ? "" : ",", esc);
                free(esc);
                first = false;
            }
            if (off > 0 && (size_t)off < cap) snprintf(buf + off, cap - (size_t)off, "]}");
            pipe_emit(buf);
            free(buf);
        }
        expr_free(evaluated);
        char done[64];
        snprintf(done, sizeof(done), "{\"id\":%d,\"type\":\"done\",\"memory\":%llu}",
                 id, (unsigned long long)pipe_memory_bytes());
        pipe_emit(done);
        return;
    }

    if (info_query && evaluated->type == EXPR_STRING) {
        const char* doc = evaluated->data.string;
        size_t cap = strlen(doc) * 6 + 8;
        char* esc = malloc(cap);
        if (esc) {
            json_escape(doc, esc, cap);
            size_t bcap = cap + 64;
            char* buf = malloc(bcap);
            if (buf) {
                const char* isym = info_sym[0] ? info_sym : NULL;
                if (isym)
                    snprintf(buf, bcap,
                        "{\"id\":%d,\"type\":\"usage\",\"payload\":\"%s\","
                        "\"symbol\":\"%s\"}", id, esc, isym);
                else
                    snprintf(buf, bcap,
                        "{\"id\":%d,\"type\":\"usage\",\"payload\":\"%s\"}", id, esc);
                pipe_emit(buf);
                free(buf);
            }
            free(esc);
        }
        expr_free(evaluated);
        char done[64];
        snprintf(done, sizeof(done), "{\"id\":%d,\"type\":\"done\",\"memory\":%llu}",
                 id, (unsigned long long)pipe_memory_bytes());
        pipe_emit(done);
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
        if (ih && (strcmp(ih, "Image") == 0 || strcmp(ih, "Image3D") == 0)) {
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
                char idone[64];
                snprintf(idone, sizeof(idone), "{\"id\":%d,\"type\":\"done\",\"memory\":%llu}",
                 id, (unsigned long long)pipe_memory_bytes());
                pipe_emit(idone);
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
                strcpy(json_line, "{\"id\":");
                char id_buf[32]; snprintf(id_buf, sizeof(id_buf), "%d", id);
                strcat(json_line, id_buf);
                strcat(json_line, ",\"type\":\"plot\",\"payload\":");
                strcat(json_line, plotly);
                strcat(json_line, "}");
                pipe_emit(json_line);
                free(json_line);
            }
            free(plotly);
            char done[64];
            snprintf(done, sizeof(done), "{\"id\":%d,\"type\":\"done\",\"memory\":%llu}",
                 id, (unsigned long long)pipe_memory_bytes());
            pipe_emit(done);
            return;
        }
        /* plotly == NULL (e.g. empty Graphics3D): fall through to text. */
        if (head_sym == SYM_Graphics || head_sym == SYM_Graphics3D) {
            expr_free(evaluated);
            char done[64];
            snprintf(done, sizeof(done), "{\"id\":%d,\"type\":\"done\",\"memory\":%llu}",
                 id, (unsigned long long)pipe_memory_bytes());
            pipe_emit(done);
            return;
        }
    }

    char* result_str = expr_to_string(evaluated);
    char* latex_raw  = expr_to_latex(evaluated);   /* must be before expr_free */
    expr_free(evaluated);

    if (!result_str) {
        char buf[256];
        snprintf(buf, sizeof(buf),
                 "{\"id\":%d,\"type\":\"error\",\"message\":\"Out of memory\"}", id);
        pipe_emit(buf);
        snprintf(buf, sizeof(buf), "{\"id\":%d,\"type\":\"done\",\"memory\":%llu}",
                 id, (unsigned long long)pipe_memory_bytes());
        pipe_emit(buf);
        return;
    }

    size_t escaped_len = strlen(result_str) * 6 + 4;
    char* escaped = malloc(escaped_len);
    if (!escaped) {
        free(result_str);
        char buf[256];
        snprintf(buf, sizeof(buf),
                 "{\"id\":%d,\"type\":\"error\",\"message\":\"Out of memory\"}", id);
        pipe_emit(buf);
        snprintf(buf, sizeof(buf), "{\"id\":%d,\"type\":\"done\",\"memory\":%llu}",
                 id, (unsigned long long)pipe_memory_bytes());
        pipe_emit(buf);
        return;
    }
    json_escape(result_str, escaped, escaped_len);
    free(result_str);

    /* latex_raw was produced above (before expr_free) — now escape it */
    char* latex_esc  = NULL;
    if (latex_raw) {
        size_t llen = strlen(latex_raw) * 6 + 4;
        latex_esc = malloc(llen);
        if (latex_esc) json_escape(latex_raw, latex_esc, llen);
        free(latex_raw);
    }

    size_t line_len = escaped_len + (latex_esc ? strlen(latex_esc) : 0) + 128;
    char* json_line = malloc(line_len);
    if (json_line) {
        if (latex_esc && strlen(latex_esc) > 0) {
            snprintf(json_line, line_len,
                     "{\"id\":%d,\"type\":\"expr\",\"payload\":\"%s\",\"latex\":\"%s\"}",
                     id, escaped, latex_esc);
        } else {
            snprintf(json_line, line_len,
                     "{\"id\":%d,\"type\":\"expr\",\"payload\":\"%s\"}", id, escaped);
        }
        pipe_emit(json_line);
        free(json_line);
    }
    free(escaped);
    free(latex_esc);

    char done[64];
    snprintf(done, sizeof(done), "{\"id\":%d,\"type\":\"done\",\"memory\":%llu}",
                 id, (unsigned long long)pipe_memory_bytes());
    pipe_emit(done);
}

static void pipe_mode_loop(void) {
    char line[MAX_INPUT_LEN];
    while (fgets(line, sizeof(line), stdin)) {
        size_t len = strlen(line);
        while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r'))
            line[--len] = '\0';
        if (len == 0) continue;

        if (strstr(line, "\"ping\"")) {
            pipe_emit("{\"type\":\"pong\"}");
            continue;
        }
        if (strstr(line, "\"quit\"")) {
            fflush(stdout);
            break;
        }

        int id = 0;
        char expr_buf[MAX_INPUT_LEN];
        if (!json_get_int(line, "id", &id) ||
            !json_get_string(line, "expr", expr_buf, sizeof(expr_buf))) {
            continue;
        }
        pipe_process_input(expr_buf, id);
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
            const char* err = skip_blanks_and_comments(stmt_start);
            if (*err != '\0') {
                report_error_at(path, buffer, err, "syntax error");
                status = 1;
            }
            break;
        }
        /* evaluate() borrows its argument and returns a new tree; both are
         * ours to free. The result is discarded: a script speaks via
         * Print[], not via echoed values. */
        Expr* evaluated = evaluate(parsed);
        expr_free(evaluated);
        expr_free(parsed);
    }

    free(buffer);

    /* A script run is a session; give $Epilog its one evaluation, exactly as
     * the interactive loop does on exit. */
    repl_apply_epilog();
    fflush(stdout);
    return status;
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

    int rc = 0;
    if (script) {
        rc = run_script_file(script);
    } else if (pipe_mode) {
        pipe_mode_loop();
    } else {
        repl_loop();
    }
    /* Tear down the context subsystem so its $Context / $ContextPath strings and
     * any open package frames are freed rather than lingering as leaks at exit.
     * context_init() is called from core_init(); this is its symmetric partner,
     * which was previously never invoked (dead code). Nothing runs after this. */
    context_shutdown();
    return rc;
}
