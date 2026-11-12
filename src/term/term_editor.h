/*
 * term_editor.h -- the in-tree multi-line line editor (the state machine).
 *
 * Pure C99 and terminal-free: bytes from the keyboard go in through
 * term_editor_feed, the drawing comes out through the write callback, one
 * call per fed batch. term_editor_tty.c puts it on a real terminal (raw mode,
 * resize, restore); the unit tests drive it with an in-memory sink.
 *
 * Editing model:
 *   - The buffer is UTF-8 text with '\n' line breaks. Up/Down move between
 *     lines (keeping a goal column) and recall history only from the first /
 *     last line. Return evaluates when term_lex_input_complete() says the
 *     buffer is complete, else opens a new line; Esc then Return (or
 *     Alt-Return) always evaluates. Continuation lines are drawn under the
 *     first one, indented by the prompt's width (visual only).
 *   - Emacs keys (Ctrl-A/E/B/F/D/H/K/U/W/Y/T/L, Alt-B/F/D/Backspace, Home,
 *     End, Delete, Ctrl/Alt-arrows), undo (Ctrl-_ or Ctrl-/), bracketed paste,
 *     Tab completion, Ctrl-R history search.
 *   - Ctrl-C clears the whole buffer (TERM_ED_INTERRUPT); Ctrl-D on an empty
 *     buffer is end of input (TERM_ED_EOF); Ctrl-Z asks the driver to suspend.
 *   - Live colouring through term_highlight, relexing only from the edited
 *     line; the bracket at or before the cursor and its partner are
 *     highlighted, a closer that matches no opener is shown in red.
 *
 * Performance contract (tests/test_term_editor_latency.c): one write per fed
 * batch, only changed rows redrawn, no allocation in steady-state editing
 * (every allocation goes through a counted hook).
 */
#ifndef MATHILDA_TERM_EDITOR_H
#define MATHILDA_TERM_EDITOR_H

#include <stdbool.h>
#include <stddef.h>

#include "term_complete.h"
#include "term_highlight.h"
#include "term_render.h"

struct TermHistory;

typedef enum {
    TERM_ED_CONTINUE = 0,   /* keep feeding                                 */
    TERM_ED_ACCEPT,         /* a line was accepted: term_editor_text()      */
    TERM_ED_EOF,            /* Ctrl-D on an empty buffer                    */
    TERM_ED_INTERRUPT,      /* Ctrl-C: the buffer was discarded             */
    TERM_ED_SUSPEND         /* Ctrl-Z: the driver suspends, then redraws    */
} TermEdStatus;

/* The one-line description of a completion candidate (e.g. its docstring's
 * first line); NULL for none. The returned text is only read before the next
 * call. */
typedef const char* (*TermDocFn)(const char* name, size_t len, void* ctx);

typedef struct {
    int                 cols, rows;    /* terminal size                        */
    bool                color;         /* SGR allowed                          */
    struct TermHistory* history;       /* entries for Up/Down and Ctrl-R, or NULL */
    TermHlClassify      classify;      /* symbol classes, or NULL              */
    void*               classify_ctx;
    TermDocFn           doc;           /* popup descriptions, or NULL          */
    void*               doc_ctx;
    TermWriteFn         write;         /* the terminal                         */
    void*               write_ctx;
    /* Completion candidates; NULL uses term_complete (the symbol table). */
    void (*complete)(const char* buf, size_t cursor, TermCompletions* out);
} TermEditorConfig;

typedef struct TermEditor TermEditor;

TermEditor*  term_editor_new(const TermEditorConfig* cfg);
void         term_editor_free(TermEditor* ed);

/* Start reading one entry: an empty buffer and the prompt (plain ASCII; it is
 * drawn in the prompt style when colour is on). cursor_abs_row is the
 * terminal row of the cursor (0-based) if known, else -1; start_col > 0 means
 * the cursor is not at the start of a line, so the prompt goes on a new one. */
void         term_editor_begin(TermEditor* ed, const char* prompt, int cursor_abs_row,
                               int start_col);

/* Feed input bytes (n may be 0). Keys are applied in order and the screen is
 * updated once, with one write, at the end. Processing stops at the first key
 * that ends the read (accept, EOF, interrupt, suspend); the bytes after it are
 * kept and processed by the first feed after the next term_editor_begin
 * (typeahead). An incomplete escape sequence at the end is kept until more
 * bytes arrive or term_editor_timeout is called. */
TermEdStatus term_editor_feed(TermEditor* ed, const char* bytes, size_t n);

/* True when a lone ESC (or a partial sequence) waits for its continuation;
 * the driver calls term_editor_timeout after ~50 ms without input. */
bool         term_editor_pending(const TermEditor* ed);
TermEdStatus term_editor_timeout(TermEditor* ed);

/* Set the terminal size between reads (draws nothing). */
void         term_editor_set_size(TermEditor* ed, int cols, int rows);
/* The terminal changed size during a read: redraw the area from its top. */
void         term_editor_resize(TermEditor* ed, int cols, int rows);
/* Redraw from scratch at the cursor's line (after a suspend), or after
 * clearing the screen (clear = true, as Ctrl-L does). */
void         term_editor_redraw(TermEditor* ed, bool clear);

/* The buffer (NUL-terminated) after TERM_ED_ACCEPT. */
const char*  term_editor_text(const TermEditor* ed, size_t* len);

/* Test hooks: allocations made by any editor since start-up, and sink writes
 * made by this editor. */
size_t       term_editor_alloc_count(void);
size_t       term_editor_write_count(const TermEditor* ed);
/* Test hook: the cursor's byte offset in the buffer. */
size_t       term_editor_cursor(const TermEditor* ed);

#endif /* MATHILDA_TERM_EDITOR_H */
