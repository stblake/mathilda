/*
 * term_interrupt.h -- the user-interrupt flag (Ctrl-C) and its public API.
 *
 * Readline-free; compiled into the Mathilda binary AND libmathilda.a (only
 * src/term/term_readline.o is filtered out of LIB_OBJ), so an embedder that
 * never installs a handler still links the flag the evaluator polls.
 *
 * How an interrupt travels: something (the REPL's SIGINT handler, or a
 * subprocess builtin that saw its child die of SIGINT) calls
 * mth_interrupt_request(). The evaluator polls the flag at the top of every
 * fixed-point iteration and turns it into the ordinary Abort[] sentinel, so
 * every enclosing frame unwinds through its normal cleanup. The flag is
 * sticky: it stays set until the abort is consumed by the top level
 * (eval_report_uncaught_abort, also reached by every Dialog line) or by a
 * CheckAbort, so a speculative probe that drops the sentinel is re-aborted at
 * its next step. AbortProtect defers it: mth_interrupt_should_stop() is false
 * while an AbortProtect is active.
 *
 * Threads: the flag is a single volatile sig_atomic_t. Worker threads read it
 * (through tc_deadline_passed(), e.g. fi_stop_requested in FactorInteger's
 * parallel driver) without synchronisation. That read is racy by design and
 * benign: a worker sees the store a little late at worst, and the value only
 * ever moves between 0 and 1.
 */
#ifndef MATHILDA_TERM_INTERRUPT_H
#define MATHILDA_TERM_INTERRUPT_H

#include <signal.h>
#include <stdbool.h>

/* ---- public API (a fixed contract with other branches; do not change) ---- */

/* True while a user interrupt is in flight: from the request until the top level
 * (eval_report_uncaught_abort, also reached by every Dialog line) or a CheckAbort
 * consumes the resulting Abort[]. */
bool mth_interrupt_pending(void);

/* Raise an interrupt as if Ctrl-C had been pressed: sets the same flag the SIGINT
 * handler sets, so the evaluator's next poll raises the Abort[] sentinel.
 * Async-signal-safe (a single sig_atomic_t store); callable from any thread.
 *
 * A builtin that waits on a child process (system(), popen()/pclose(), any
 * blocking wait) must call this when the child died of SIGINT or exited with
 * status 130: system() ignores SIGINT in the parent while it waits, so the
 * REPL's handler never sees the keypress. */
void mth_interrupt_request(void);

/* ---- kernel-internal (stable, but not part of the cross-branch contract) ---- */
extern volatile sig_atomic_t g_mth_interrupt_pending;
#define MTH_INTERRUPT_PENDING() (g_mth_interrupt_pending != 0)   /* evaluator hot path */
void mth_interrupt_clear(void);          /* consumers only: top level, CheckAbort, REPL */
bool mth_interrupt_should_stop(void);    /* pending && abort_protect_depth() == 0     */

/* ---- REPL SIGINT handler (interactive readline loop only) -----------------
 *
 * The handler only stores sig_atomic_t values. What a Ctrl-C means depends on
 * what the REPL is doing:
 *   IDLE        ignored (start-up, shut-down, between states);
 *   READING     at a prompt: sets the reading flag, which the line reader
 *               turns into "clear the line and prompt again";
 *   EVALUATING  requests an interrupt (the evaluation ends as $Aborted). A
 *               second Ctrl-C within 1000 ms, still evaluating, writes a
 *               newline to stderr and exits with status 130 -- the escape
 *               hatch for a computation that never polls.
 * Nothing here is installed unless the REPL calls term_interrupt_install(): an
 * embedder, pipe mode and -file runs get no handler. */
typedef enum { TERM_STATE_IDLE, TERM_STATE_READING, TERM_STATE_EVALUATING } TermState;

/* POSIX: sigaction(SIGINT) with SA_RESTART (a Ctrl-C during Export or Get
 * must not make a stdio call fail with EINTR), saving the old disposition. */
void      term_interrupt_install(void);
/* Restore the disposition saved by term_interrupt_install. */
void      term_interrupt_uninstall(void);
/* Set the state; returns the previous one. Entering EVALUATING resets the
 * double-tap timer, so a Ctrl-C at the prompt never counts towards a quit. */
TermState term_interrupt_set_state(TermState s);
TermState term_interrupt_state(void);
/* Read and clear the prompt interrupt (the READING-state flag). */
bool      term_interrupt_take_reading(void);

/* Test seams: the handler body with an injected clock (milliseconds), and the
 * hard-exit function (default: C99 _Exit). */
void      term_interrupt_on_signal_at(long now_ms);
void      term_interrupt_set_exit_fn(void (*fn)(int status));

#endif /* MATHILDA_TERM_INTERRUPT_H */
