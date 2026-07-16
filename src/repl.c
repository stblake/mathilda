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
#include "print_latex.h"
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

/* True if `s` consists only of whitespace and (* ... *) comments
 * (nested comments allowed). Used to distinguish a no-op line from a
 * genuine parse failure, so the REPL doesn't shout "Parse error" at a
 * stray comment. */
static int is_blank_or_comment_only(const char* s) {
    while (*s) {
        if (isspace((unsigned char)*s)) {
            s++;
        } else if (s[0] == '(' && s[1] == '*') {
            int depth = 1;
            s += 2;
            while (*s && depth > 0) {
                if (s[0] == '(' && s[1] == '*') { depth++; s += 2; }
                else if (s[0] == '*' && s[1] == ')') { depth--; s += 2; }
                else { s++; }
            }
            if (depth > 0) return 0;  /* unterminated comment is a real error */
        } else {
            return 0;
        }
    }
    return 1;
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
    symtab_add_down_value("Out", out_pattern, evaluated);
    expr_free(out_pattern);

    /* $PrePrint: applied only for display. Out[n] keeps the
     * pre-$PrePrint value above; here we render a possibly modified
     * copy. */
    Expr* to_print = repl_apply_pre_print(expr_copy(evaluated));

    printf("Out[%d]= ", line_number);
    expr_print(to_print);
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
    if (to_print && to_print->type == EXPR_FUNCTION
        && to_print->data.function.head->type == EXPR_SYMBOL
        && to_print->data.function.arg_count >= 1) {
        if (to_print->data.function.head->data.symbol.name == SYM_Graphics) graphics_show(to_print);
        else if (to_print->data.function.head->data.symbol.name == SYM_Graphics3D) graphics3d_show(to_print);
    }

    expr_free(to_print);
    expr_free(parsed);
    expr_free(evaluated);
}

#ifndef NO_READLINE
void repl_loop() {
    printf("\nMathilda - A tiny, LLM-generated, Mathematica-like computer algebra system.\n\n");
    printf("This program is free, open source software and comes with ABSOLUTELY NO WARRANTY.\n\n");
    printf("End a line with '\\' to enter a multiline expression. Press Return to evaluate.\n");
    printf("Exit by evaluating Quit[] or CONTROL-C.\n\n");

    int line_number = 1;
    char prompt[64];
    char full_input[MAX_INPUT_LEN] = {0};
    int in_multiline = 0;

    while (1) {
        int submit_now = 0;
        if (!in_multiline) {
            snprintf(prompt, sizeof(prompt), "In[%d]:= ", line_number);
        } else {
            prompt[0] = '\0';
        }

        char* line = readline(prompt);
        if (!line) {
            printf("\n");
            /* EOF: run $Epilog before tearing down. */
            repl_apply_epilog();
            break;
        }

        size_t line_len = strlen(line);

        /* Remove trailing whitespace to properly check for backslash. */
        size_t check_len = line_len;
        while (check_len > 0 && (line[check_len - 1] == ' ' || line[check_len - 1] == '\t')) {
            check_len--;
        }

        int has_backslash = 0;
        if (check_len > 0 && line[check_len - 1] == '\\') {
            has_backslash = 1;
            line[check_len - 1] = '\0';
            line_len = strlen(line);
        }

        /* Check buffer limits. */
        if (strlen(full_input) + line_len + 2 >= MAX_INPUT_LEN) {
            printf("Input too long!\n");
            full_input[0] = '\0';
            in_multiline = 0;
            free(line);
            continue;
        }

        if (in_multiline) {
            strcat(full_input, "\n");
        }
        strcat(full_input, line);

        if (has_backslash) {
            in_multiline = 1;
        } else {
            submit_now = 1;
        }

        if (submit_now) {
            if (strlen(full_input) == 0) {
                free(line);
                continue;
            }

            add_history(full_input);

            if (strcmp(full_input, "Quit[]") == 0) {
                /* User-requested shutdown: run $Epilog first. */
                repl_apply_epilog();
                free(line);
                break;
            }

            process_input(full_input, line_number);

            full_input[0] = '\0';
            in_multiline = 0;
            line_number++;
        }

        free(line);
    }

    printf("\n");
}
#else
/* Fallback interactive loop when readline is not available (e.g. Windows).
 * Uses fgets; no history or line-editing. Pipe mode bypasses this entirely. */
void repl_loop(void) {
    printf("\nMathilda - A tiny, Mathematica-like computer algebra system.\n\n");
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
        snprintf(dbuf, sizeof(dbuf), "{\"id\":%d,\"type\":\"done\"}", id);
        pipe_emit(dbuf);
        return;
    }

    Expr* evaluated = evaluate(parsed);
    expr_free(parsed);

    if (!evaluated) {
        char buf[64];
        snprintf(buf, sizeof(buf), "{\"id\":%d,\"type\":\"done\"}", id);
        pipe_emit(buf);
        return;
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
            snprintf(done, sizeof(done), "{\"id\":%d,\"type\":\"done\"}", id);
            pipe_emit(done);
            return;
        }
        /* plotly == NULL (e.g. empty Graphics3D): fall through to text. */
        if (head_sym == SYM_Graphics || head_sym == SYM_Graphics3D) {
            expr_free(evaluated);
            char done[64];
            snprintf(done, sizeof(done), "{\"id\":%d,\"type\":\"done\"}", id);
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
        snprintf(buf, sizeof(buf), "{\"id\":%d,\"type\":\"done\"}", id);
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
        snprintf(buf, sizeof(buf), "{\"id\":%d,\"type\":\"done\"}", id);
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
    snprintf(done, sizeof(done), "{\"id\":%d,\"type\":\"done\"}", id);
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

int main(void) {
    /* Detect pipe mode: when stdin is not a terminal the frontend has
     * spawned us as a sidecar and we communicate via NDJSON over stdio.
     * The interactive readline REPL is preserved when stdin is a tty. */
    int pipe_mode = !isatty(fileno(stdin));

    if (pipe_mode) {
        /* Disable libc's stdout buffer so every response line is delivered
         * to the pipe immediately rather than accumulating. */
        setvbuf(stdout, NULL, _IONBF, 0);
    }

    symtab_init();
    core_init();

    /* Load the internal bootstrap (init.m). Path resolution is independent of
     * the current working directory (see mathilda_load_module), so a relocated
     * or installed binary still finds its bundled src/internal tree. If it
     * cannot be located the loader prints a LoadModule::nofile diagnostic —
     * far better than the previous silent load of a non-functional kernel. */
    mathilda_load_module("init.m");

    if (pipe_mode) {
        pipe_mode_loop();
    } else {
        repl_loop();
    }
    return 0;
}
