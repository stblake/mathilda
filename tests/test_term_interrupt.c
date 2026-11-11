/*
 * test_term_interrupt.c -- the REPL's SIGINT state machine
 * (src/term/term_interrupt.c) and the stdio line reader (src/term/term_input.c).
 *
 * The handler body is driven through term_interrupt_on_signal_at with an
 * injected clock, and the hard exit is replaced by a recorder, so the
 * double-tap rule is tested without killing the test.
 */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "test_utils.h"
#include "term/term_interrupt.h"
#include "term/term_input.h"

#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_exit_status = -1;
static void record_exit(int status) { g_exit_status = status; }

static void reset(void) {
    term_interrupt_set_state(TERM_STATE_IDLE);
    (void)term_interrupt_take_reading();
    mth_interrupt_clear();
    g_exit_status = -1;
}

static void test_reading(void) {
    /* 1. At a prompt: the reading flag, never the kernel flag. */
    reset();
    term_interrupt_set_state(TERM_STATE_READING);
    term_interrupt_on_signal_at(1000);
    ASSERT(!mth_interrupt_pending());
    ASSERT(term_interrupt_take_reading());
    ASSERT(!term_interrupt_take_reading());          /* read-and-clear */
    /* Two at a prompt only clear the line twice. */
    term_interrupt_on_signal_at(1100);
    term_interrupt_on_signal_at(1200);
    ASSERT(g_exit_status == -1);
    ASSERT(!mth_interrupt_pending());
}

static void test_evaluating(void) {
    /* 2. Evaluating: the kernel flag; a second within 1 s exits 130. */
    reset();
    term_interrupt_set_state(TERM_STATE_EVALUATING);
    term_interrupt_on_signal_at(1000);
    ASSERT(mth_interrupt_pending());
    ASSERT(!term_interrupt_take_reading());
    ASSERT(g_exit_status == -1);
    term_interrupt_on_signal_at(1800);
    ASSERT(g_exit_status == 130);

    /* 3. More than 1 s apart: no exit. */
    reset();
    term_interrupt_set_state(TERM_STATE_EVALUATING);
    term_interrupt_on_signal_at(1000);
    term_interrupt_on_signal_at(2500);
    ASSERT(g_exit_status == -1);
    ASSERT(mth_interrupt_pending());

    /* Across the 2^30 ms wrap. */
    reset();
    term_interrupt_set_state(TERM_STATE_EVALUATING);
    term_interrupt_on_signal_at((1L << 30) - 200);
    term_interrupt_on_signal_at((1L << 30) + 300);
    ASSERT(g_exit_status == 130);
}

static void test_idle_and_reset(void) {
    /* 4. Idle: ignored. */
    reset();
    term_interrupt_on_signal_at(1000);
    ASSERT(!mth_interrupt_pending());
    ASSERT(!term_interrupt_take_reading());

    /* 5. Re-entering EVALUATING resets the double-tap timer. */
    reset();
    term_interrupt_set_state(TERM_STATE_EVALUATING);
    term_interrupt_on_signal_at(1000);
    term_interrupt_set_state(TERM_STATE_READING);
    term_interrupt_set_state(TERM_STATE_EVALUATING);
    term_interrupt_on_signal_at(1500);
    ASSERT(g_exit_status == -1);

    /* set_state returns the previous state. */
    reset();
    ASSERT(term_interrupt_set_state(TERM_STATE_READING) == TERM_STATE_IDLE);
    ASSERT(term_interrupt_set_state(TERM_STATE_EVALUATING) == TERM_STATE_READING);
    ASSERT(term_interrupt_state() == TERM_STATE_EVALUATING);
}

static void test_real_signal(void) {
    /* 6. A real SIGINT through the installed handler. */
    reset();
    struct sigaction before, during, after;
    sigaction(SIGINT, NULL, &before);
    term_interrupt_install();
    sigaction(SIGINT, NULL, &during);
    ASSERT(during.sa_handler != before.sa_handler);
    ASSERT((during.sa_flags & SA_RESTART) != 0);
    term_interrupt_set_state(TERM_STATE_EVALUATING);
    raise(SIGINT);
    ASSERT(mth_interrupt_pending());
    term_interrupt_uninstall();
    sigaction(SIGINT, NULL, &after);
    ASSERT(after.sa_handler == before.sa_handler);
    ASSERT(after.sa_flags == before.sa_flags);
    reset();
}

static void test_stdio_reader(void) {
    /* 7. The stdio backend over a file-backed stdin: each line, then EOF. */
    char path[64];
    snprintf(path, sizeof(path), "/tmp/mth_term_input_%ld.txt", (long)getpid());
    FILE* f = fopen(path, "w");
    ASSERT(f != NULL);
    fputs("1+1\nf[x]\r\n\nlast", f);
    fclose(f);
    ASSERT(freopen(path, "r", stdin) != NULL);

    TermReader* rd = term_reader_open_stdio();
    ASSERT(rd && term_reader_kind(rd) == TERM_BACKEND_STDIO);
    const char* want[] = { "1+1", "f[x]", "", "last" };
    for (size_t i = 0; i < 4; i++) {
        char* line = NULL;
        term_interrupt_set_state(TERM_STATE_EVALUATING);
        TermReadStatus st = term_read(rd, "", &line);
        ASSERT(st == TERM_READ_OK);
        ASSERT_MSG(strcmp(line, want[i]) == 0, "line %zu is %s", i, line);
        /* The previous state is restored after the read. */
        ASSERT(term_interrupt_state() == TERM_STATE_EVALUATING);
        term_reader_history_add(rd, line);       /* a no-op for stdio */
        free(line);
    }
    char* line = (char*)"sentinel";
    ASSERT(term_read(rd, "", &line) == TERM_READ_EOF);
    ASSERT(line == NULL);
    term_reader_close(rd);
    term_reader_close(NULL);
    remove(path);
    reset();
}

int main(void) {
    term_interrupt_set_exit_fn(record_exit);
    TEST(test_reading);
    TEST(test_evaluating);
    TEST(test_idle_and_reset);
    TEST(test_real_signal);
    TEST(test_stdio_reader);
    printf("All term_interrupt tests passed.\n");
    return 0;
}
