/*
 * term_readline.c -- see term_readline.h. The only translation unit that
 * includes readline; under NO_READLINE it compiles to nothing.
 */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "term_readline.h"

typedef int term_readline_unused_;   /* keep the translation unit non-empty */

#ifndef NO_READLINE

#include "term_interrupt.h"
#include "term_lex.h"
#include "term_history.h"
#include "term_complete.h"
#include "term_color.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <readline/readline.h>
#include <readline/history.h>

/* GNU Readline >= 6.3 has rl_signal_event_hook and the line-state API; Apple's
 * libedit shim defines RL_READLINE_VERSION 0x0402 and no RL_STATE_INITIALIZED. */
#if defined(RL_READLINE_VERSION) && defined(RL_STATE_INITIALIZED) && RL_READLINE_VERSION >= 0x0603
#define TERM_RL_GNU 1
#endif

typedef struct {
    TermReadlineOptions opt;
} RlImpl;

/* Return key handler: evaluate the buffer once it forms a complete
 * expression, otherwise open a new line so the user can keep typing. This is
 * the terminal-native substitute for a notebook's Shift+Enter — it needs no
 * enhanced-keyboard protocol and works on every terminal. */
static int mth_smart_return(int count, int key) {
    if (term_lex_input_complete(rl_line_buffer ? rl_line_buffer : ""))
        return rl_newline(count, key);   /* accept the whole buffer */
    rl_insert_text("\n");                /* still open: continue editing */
    return 0;
}

/* Esc-Return (and Alt/Meta-Return, where the terminal sends it) forces the
 * buffer to be evaluated even while term_lex_input_complete() still considers it
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

#ifdef TERM_RL_GNU
/* Called by readline (not in signal context) after a signal interrupted its read.
 * readline's own SIGINT handler has already echoed ^C and re-raised to ours, which
 * set the reading flag. Clear the whole buffer, including a multi-line one, and
 * redraw a fresh prompt below it. */
static int term_rl_on_signal(void) {
    if (!term_interrupt_take_reading()) return 0;
    rl_point = rl_end; rl_redisplay();        /* cursor below the old buffer   */
    rl_replace_line("", 0); rl_point = rl_end = 0;
    rl_crlf(); rl_on_new_line(); rl_redisplay();
    return 0;
}
#endif

/* ---- Tab completion ------------------------------------------------------ */

/* Longest common prefix of the candidates, as readline's match[0]. */
static char* common_prefix(char** items, size_t n) {
    size_t len = strlen(items[0]);
    for (size_t i = 1; i < n; i++) {
        size_t k = 0;
        while (k < len && items[i][k] == items[0][k]) k++;
        len = k;
    }
    char* r = malloc(len + 1);
    if (!r) return NULL;
    memcpy(r, items[0], len);
    r[len] = '\0';
    return r;
}

static char** term_rl_complete(const char* text, int start, int end) {
    (void)text; (void)end;
    rl_attempted_completion_over = 1;              /* never fall back to file names */
    rl_completion_append_character = '\0';         /* GNU resets it per attempt      */
    TermCompletions c;
    term_complete(rl_line_buffer, (size_t)rl_point, &c);
    if (c.count == 0 || c.word_start != (size_t)start) { term_completions_free(&c); return NULL; }
    /* readline's contract: [0] = longest common prefix, [1..n] = matches, NULL. */
    char** m = malloc((c.count + 2) * sizeof(char*));
    char* lcp = m ? common_prefix(c.items, c.count) : NULL;
    if (!m || !lcp) { free(m); term_completions_free(&c); return NULL; }
    m[0] = lcp;
    for (size_t i = 0; i < c.count; i++) m[i + 1] = c.items[i];   /* ownership moves */
    m[c.count + 1] = NULL;
    free(c.items);                                  /* the array only */
    return m;
}

/* Word breaks: every ASCII character outside [A-Za-z0-9$`], so readline's
 * `start` is the term_lex word start (after a '"', which is a break). */
static void term_rl_setup_completion(void) {
    static char breaks[128];
    size_t k = 0;
    for (int ch = 1; ch < 128; ch++)
        if (!term_lex_is_word_char((char)ch)) breaks[k++] = (char)ch;
    breaks[k] = '\0';
    rl_completer_word_break_characters = breaks;
    rl_attempted_completion_function = term_rl_complete;
#ifndef TERM_RL_GNU
    rl_completion_append_character = '\0';
#endif
}

static TermReadStatus rl_read(TermReader* rd, const char* prompt, char** out) {
    RlImpl* im = (RlImpl*)rd->impl;
    char styled[256];
#ifdef TERM_RL_GNU
    /* The coloured prompt, with every SGR inside \001...\002 so readline
     * does not count it as width. libedit: plain (its support for the
     * markers is unverified, and a garbled redraw is worse than no colour). */
    if (im && im->opt.interactive_features && term_color_stdout()
        && term_style_wrap(styled, sizeof(styled), TERM_STYLE_IN, prompt, true, true)
           < (int)sizeof(styled))
        prompt = styled;
#else
    (void)im; (void)styled;
#endif
    char* line = readline(prompt);
    if (!line) {
        /* libedit has no signal hook: a Ctrl-C there can end readline() with
         * NULL. With the reading flag set that is an interrupted prompt, not
         * end of input -- the session must never exit on a Ctrl-C. */
        if (term_interrupt_take_reading()) return TERM_READ_INTERRUPTED;
        return TERM_READ_EOF;
    }
    *out = line;
    return TERM_READ_OK;
}

/* With a history store (the interactive session), the store filters the
 * entry (blank, comment-only, a repeat of the previous one) and appends it to
 * the file; readline's own list gets exactly what the store kept. Without
 * one, the plain -file Dialog reader's long-standing behaviour. */
static void rl_history_add(TermReader* rd, const char* entry) {
    RlImpl* im = (RlImpl*)rd->impl;
    if (im && im->opt.history) {
        if (term_history_add(im->opt.history, entry)) add_history(entry);
        return;
    }
    if (*entry) add_history(entry);
}

static void rl_close(TermReader* rd) {
    RlImpl* im = (RlImpl*)rd->impl;
#ifdef TERM_RL_GNU
    if (im && im->opt.interactive_features) rl_signal_event_hook = NULL;
#endif
    free(im);
}

static const TermReaderOps RL_OPS = { rl_read, rl_history_add, rl_close };

TermReader* term_reader_open_readline(const TermReadlineOptions* o) {
    TermReader* rd = malloc(sizeof(*rd));
    RlImpl* im = malloc(sizeof(*im));
    if (!rd || !im) { free(rd); free(im); return NULL; }
    memset(im, 0, sizeof(*im));
    if (o) im->opt = *o;
    rd->ops = &RL_OPS;
    rd->kind = TERM_BACKEND_READLINE;
    rd->impl = im;
    if (im->opt.interactive_features && im->opt.history) {
        TermHistory* h = im->opt.history;
        using_history();
        stifle_history((int)term_history_max(h));
        for (size_t i = 0; i < term_history_count(h); i++) add_history(term_history_get(h, i));
    }
    if (im->opt.interactive_features) {
        mth_setup_readline();
        term_rl_setup_completion();
#ifdef TERM_RL_GNU
        rl_signal_event_hook = term_rl_on_signal;
#endif
    }
    return rd;
}

#endif /* NO_READLINE */
