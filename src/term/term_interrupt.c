/*
 * term_interrupt.c -- see term_interrupt.h.
 *
 * Readline-free and in libmathilda.a: an embedder that never calls
 * term_interrupt_install() gets the flag and no handler.
 */
/* sigaction, clock_gettime and write are POSIX, hidden by glibc under -std=c99. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "term_interrupt.h"

#include "funcprog.h"   /* abort_protect_depth */

#include <stdlib.h>
#include <time.h>
#ifndef _WIN32
#include <unistd.h>
#endif

/* ---- the flag ------------------------------------------------------------ */

volatile sig_atomic_t g_mth_interrupt_pending = 0;

bool mth_interrupt_pending(void) {
    return g_mth_interrupt_pending != 0;
}

void mth_interrupt_request(void) {
    g_mth_interrupt_pending = 1;
}

void mth_interrupt_clear(void) {
    g_mth_interrupt_pending = 0;
}

/* The depth is read only when the flag is set, so the common (no interrupt)
 * case costs one load. */
bool mth_interrupt_should_stop(void) {
    return g_mth_interrupt_pending != 0 && abort_protect_depth() == 0;
}

/* ---- the REPL handler ---------------------------------------------------- */

#define TERM_DOUBLE_TAP_MS 1000L
#define TERM_MS_MASK ((1L << 30) - 1)    /* times are kept modulo 2^30 ms */

static volatile sig_atomic_t g_state = TERM_STATE_IDLE;
static volatile sig_atomic_t g_reading_interrupt = 0;
static volatile sig_atomic_t g_have_last = 0;
static volatile sig_atomic_t g_last_ms = 0;
/* stderr is a terminal: cached at install time because isatty() is not
 * async-signal-safe. Used to end the tty's echoed "^C" with a newline. */
static volatile sig_atomic_t g_stderr_tty = 0;

static void default_exit(int status) { _Exit(status); }
static void (*g_exit_fn)(int) = default_exit;

void term_interrupt_set_exit_fn(void (*fn)(int status)) {
    g_exit_fn = fn ? fn : default_exit;
}

TermState term_interrupt_set_state(TermState s) {
    TermState prev = (TermState)g_state;
    if (s == TERM_STATE_EVALUATING && prev != TERM_STATE_EVALUATING) g_have_last = 0;
    g_state = (sig_atomic_t)s;
    return prev;
}

TermState term_interrupt_state(void) {
    return (TermState)g_state;
}

bool term_interrupt_take_reading(void) {
    bool r = g_reading_interrupt != 0;
    g_reading_interrupt = 0;
    return r;
}

/* The handler body. Async-signal-safe: sig_atomic_t stores, write(2), _Exit. */
void term_interrupt_on_signal_at(long now_ms) {
    long now = now_ms & TERM_MS_MASK;
    switch ((TermState)g_state) {
    case TERM_STATE_READING:
        g_reading_interrupt = 1;
        break;
    case TERM_STATE_EVALUATING:
        if (g_have_last && ((now - (long)g_last_ms) & TERM_MS_MASK) < TERM_DOUBLE_TAP_MS) {
#ifndef _WIN32
            if (write(2, "\n", 1) < 0) { /* nothing to do: we are exiting */ }
#endif
            g_exit_fn(130);
            break;   /* reached only with a test exit function */
        }
        mth_interrupt_request();
#ifndef _WIN32
        /* The tty driver echoes "^C" with no newline; without this the
         * result prints on that line ("^COut[3]= $Aborted"). */
        if (g_stderr_tty && write(2, "\n", 1) < 0) { /* best effort */ }
#endif
        g_have_last = 1;
        g_last_ms = (sig_atomic_t)now;
        break;
    case TERM_STATE_IDLE:
    default:
        break;
    }
}

#ifndef _WIN32

static long monotonic_ms(void) {
    struct timespec ts;
    if (clock_gettime(CLOCK_MONOTONIC, &ts) != 0) return 0;
    return (long)((ts.tv_sec % 1000000L) * 1000L + ts.tv_nsec / 1000000L);
}

static void on_sigint(int sig) {
    (void)sig;
    term_interrupt_on_signal_at(monotonic_ms());
}

static struct sigaction g_saved_action;
static bool g_installed = false;

void term_interrupt_install(void) {
    if (g_installed) return;
    g_stderr_tty = isatty(2) ? 1 : 0;
    struct sigaction sa;
    sa.sa_handler = on_sigint;
    sigemptyset(&sa.sa_mask);
    sa.sa_flags = SA_RESTART;
    if (sigaction(SIGINT, &sa, &g_saved_action) == 0) g_installed = true;
}

void term_interrupt_uninstall(void) {
    if (!g_installed) return;
    sigaction(SIGINT, &g_saved_action, NULL);
    g_installed = false;
}

#else  /* _WIN32: never installed (NO_READLINE keeps the fgets loop), but it builds. */

static void (*g_saved_handler)(int) = SIG_DFL;
static bool g_installed = false;

static void on_sigint(int sig) {
    signal(SIGINT, on_sigint);   /* re-arm: Windows resets to SIG_DFL on entry */
    (void)sig;
    term_interrupt_on_signal_at((long)(clock() * 1000L / CLOCKS_PER_SEC));
}

void term_interrupt_install(void) {
    if (g_installed) return;
    g_saved_handler = signal(SIGINT, on_sigint);
    g_installed = true;
}

void term_interrupt_uninstall(void) {
    if (!g_installed) return;
    signal(SIGINT, g_saved_handler);
    g_installed = false;
}

#endif
