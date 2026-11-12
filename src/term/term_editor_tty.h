/*
 * term_editor_tty.h -- the in-tree editor on a real terminal: the
 * TERM_BACKEND_EDITOR line reader.
 *
 * POSIX only (termios, poll, SIGWINCH); on Windows term_reader_open_editor
 * returns NULL and the REPL keeps its other reader.
 *
 * While a line is being read the terminal is in raw mode with ISIG off, so
 * Ctrl-C, Ctrl-Z and Ctrl-\ arrive as keys. Before term_read returns, the
 * termios saved at the start of the read is restored, bracketed paste is
 * switched off and the editor's SIGWINCH and crash-signal handlers are
 * removed: from then on, until the next term_read, the terminal belongs to
 * whoever else draws on it (the evaluation's output, a status line). The
 * terminal is also restored by an atexit handler and, while reading, on
 * SIGTERM, SIGHUP, SIGQUIT, SIGSEGV, SIGBUS, SIGFPE, SIGILL and SIGABRT (the
 * previous disposition is then put back and the signal re-raised).
 */
#ifndef MATHILDA_TERM_EDITOR_TTY_H
#define MATHILDA_TERM_EDITOR_TTY_H

#include <stdbool.h>
#include <stddef.h>
#include "term_input.h"

struct TermHistory;

typedef struct {
    struct TermHistory* history;   /* the session's store, or NULL       */
    bool                color;     /* the session's stdout colour choice */
    bool                slash;     /* '/' opens the slash-command menu   */
} TermEditorOptions;

/* Whether the interactive session should use the in-tree editor:
 * both stdin and stdout are terminals, TERM is set and not "dumb", and
 * MATHILDA_EDITOR is unset, empty, or anything but "readline", "libedit",
 * "plain", "none", "off" or "0". Pure; the caller passes the environment. */
bool        term_editor_select(const char* env_editor, const char* env_term,
                               bool stdin_tty, bool stdout_tty);

/* The editor backend, or NULL (not a terminal, Windows, out of memory). */
TermReader* term_reader_open_editor(const TermEditorOptions* o);

/* Run a picker (term_editor.h) on the editor's terminal: raw mode for its
 * duration, the cursor hidden. TERM_READ_OK with the chosen index in *out,
 * TERM_READ_INTERRUPTED when cancelled, TERM_READ_EOF when the reader is not
 * the editor (the caller falls back to text) or the input went away. */
struct TermPick;
TermReadStatus term_reader_pick(TermReader* rd, const struct TermPick* p, size_t* out);

/* True while an editor reader is open: the session draws with the editor,
 * so its results may be coloured (term_result_print). */
bool        term_editor_session_active(void);

#endif /* MATHILDA_TERM_EDITOR_TTY_H */
