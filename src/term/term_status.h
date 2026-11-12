/*
 * term_status.h -- the REPL's live status line (spinner) during evaluation,
 * the "took" statistics after it, and the optional status line between
 * evaluations.
 *
 * Design (thoughts/shared/research/2026-11-11-repl-modernize-pr3.md): one
 * persistent forwarder thread, created on first use with every signal
 * blocked, owns the terminal while an evaluation runs. term_status_begin()
 * points fd 1 and fd 2 at one pipe (so stdout/stderr interleaving is kept
 * byte for byte) and clears the tty's ECHO bit (ISIG stays, so Ctrl-C still
 * raises SIGINT on the main thread). The thread copies everything from the
 * pipe to a saved descriptor of the real terminal and, between those writes,
 * draws one status line: after 200 ms, redrawn every 100 ms, only while the
 * cursor is at column 0, erased before any output and at the end.
 * term_status_end() puts the descriptors and ECHO back, then waits until the
 * thread has drained the pipe and erased the line.
 *
 * Handover with the line editor (PR 2): the editor owns the terminal from the
 * prompt until it returns a line; this module owns it from begin to end. At
 * both points the tty is in its cooked mode, fds 1/2 are the terminal, and the
 * cursor is at column 0. Dialog[] hands the terminal back to its reader with
 * term_status_suspend() / term_status_resume().
 *
 * Off (no thread, no redirection, nothing written; begin/end only measure)
 * unless term_status_init(true) decided otherwise -- see term_status_decide.
 * Pipe mode and -file never call term_status_init.
 */
#ifndef MATHILDA_TERM_STATUS_H
#define MATHILDA_TERM_STATUS_H

#include <stdbool.h>
#include <stdint.h>

struct Expr;

typedef struct {
    bool interactive;          /* the interactive readline session           */
    bool threads;              /* built with MATHILDA_THREADS, not Windows   */
    bool ttys;                 /* fd 0, 1, 2 are ttys; fd 1 and 2 the same one */
    const char* term;          /* $TERM                                       */
    const char* spinner;       /* $MATHILDA_SPINNER ("0" disables)            */
} TermStatusEnv;

/* The rule, as a pure function: interactive && threads && ttys && TERM set,
 * not "dumb" && MATHILDA_SPINNER is not "0". */
bool term_status_decide(const TermStatusEnv* env);

/* UTF-8 locale? The first non-empty of LC_ALL, LC_CTYPE, LANG names UTF-8. */
bool term_status_locale_utf8(const char* lc_all, const char* lc_ctype, const char* lang);

/* Decide once (reads the environment and the descriptors) and, when on, open
 * the pipe and install the exit / fatal-signal cleanup. Colour follows
 * term_color_stdout(), so call it after term_color_init(). Also reads
 * MATHILDA_STATUSLINE ("1" turns the between-evaluations line on) and
 * MATHILDA_SPINNER_DELAY_MS (default 200). */
void term_status_init(bool interactive);
bool term_status_enabled(void);        /* the spinner will run               */
void term_status_set_enabled(bool on); /* for a future /spinner command      */
bool term_status_unicode(void);        /* glyph set chosen at init           */

/* Around one top-level evaluation. `label` is shown next to the spinner (an
 * interned or static string; NULL: "Evaluating"). Nested calls only count. */
void term_status_begin(const char* label);

typedef struct {
    double elapsed;   /* wall seconds                                       */
    double cpu;       /* process CPU seconds, < 0 when unknown              */
    uint64_t peak;    /* peak resident bytes during the evaluation; 0 unknown */
} TermStatusStats;
void term_status_end(TermStatusStats* out);

/* Give the terminal back for a read inside an evaluation (Dialog[]), and take
 * it again afterwards. No-ops when the spinner is not running. */
void term_status_suspend(void);
void term_status_resume(void);

/* Stop the thread and close the descriptors (end of session). */
void term_status_shutdown(void);

/* The head to show for an input: Set/SetDelayed show their right-hand side,
 * CompoundExpression its last call. NULL for an atom. */
const char* term_status_label_of(const struct Expr* e);

/* The between-evaluations line ("last 3.21 s · mem 1.4 GB · Global`"). */
bool term_status_statusline(void);
void term_status_set_statusline(bool on);
void term_status_print_statusline(double last_elapsed);

/* ---- test seam ------------------------------------------------------------
 * Run the forwarder against arbitrary descriptors: tty_fd receives the
 * output and the status line, route[0..n_route-1] are pointed at the pipe
 * during an evaluation. echo_fd < 0 leaves termios alone. Replaces the
 * term_status_init decision. */
typedef struct {
    int tty_fd;
    int route[2];
    int n_route;
    int echo_fd;
    bool color, unicode;
    int delay_ms, period_ms, cols;
} TermStatusConfig;
bool term_status_open(const TermStatusConfig* cfg);

#endif /* MATHILDA_TERM_STATUS_H */
