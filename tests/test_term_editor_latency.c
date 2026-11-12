/*
 * test_term_editor_latency.c -- the editor's performance gate (REPL PR 2).
 *
 * A 200-line buffer of realistic input on a 120 x 40 colour terminal, then a
 * few thousand keystrokes fed one at a time, each timed from the feed to the
 * end of the one write it produces. Asserted:
 *   - the slowest keystroke takes under 2 ms (the budget is keystroke to
 *     screen; the write here goes to a sink, the terminal's own time is not
 *     the editor's);
 *   - no keystroke produces more than one write (one that changes nothing on
 *     screen, such as Ctrl-E at the end of a line, produces none);
 *   - steady-state editing allocates nothing;
 *   - the worst case, a '"' at the top that turns the whole buffer into a
 *     string and back, also stays under 2 ms.
 * The numbers are printed so a regression shows before it fails.
 */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L
#endif

#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include "term/term_editor.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static size_t g_bytes, g_writes;

static void sink(void* ctx, const char* s, size_t n) {
    (void)ctx; (void)s;
    g_bytes += n;
    g_writes++;
}

static double now_us(void) {
#if defined(CLOCK_MONOTONIC)
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1e6 + (double)ts.tv_nsec / 1e3;
#else
    return (double)clock() * 1e6 / (double)CLOCKS_PER_SEC;
#endif
}

static int cmp_double(const void* a, const void* b) {
    double x = *(const double*)a, y = *(const double*)b;
    return x < y ? -1 : x > y;
}

#define NKEYS 3000
static double g_t[NKEYS];

typedef struct { const char* name; const char* bytes; } KeyDef;

/* A typing-heavy mix: characters, brackets (match highlight), deletes,
 * moves within a line and between lines, line start / end. Net length 0. */
static const KeyDef CYCLE[] = {
    { "char", "a" }, { "char", "b" }, { "bracket", "[" }, { "char", "x" },
    { "bracket", "]" }, { "backspace", "\x7f" }, { "backspace", "\x7f" },
    { "backspace", "\x7f" }, { "backspace", "\x7f" }, { "backspace", "\x7f" },
    { "left", "\x1b[D" }, { "right", "\x1b[C" }, { "up", "\x1b[A" },
    { "down", "\x1b[B" }, { "ctrl-e", "\x05" }, { "ctrl-a", "\x01" },
    { "word", "\x1b" "f" }, { "char", "q" }, { "backspace", "\x7f" },
    { "ctrl-e", "\x05" },
};
#define NCYCLE (sizeof(CYCLE) / sizeof(CYCLE[0]))

static size_t g_silent;   /* keystrokes that changed nothing on screen */

static void run_keys(TermEditor* ed, size_t n, double* times, size_t* bad_writes) {
    for (size_t i = 0; i < n; i++) {
        const char* k = CYCLE[i % NCYCLE].bytes;
        size_t w0 = g_writes;
        double t0 = now_us();
        TermEdStatus st = term_editor_feed(ed, k, strlen(k));
        double t1 = now_us();
        ASSERT(st == TERM_ED_CONTINUE);
        if (times) times[i] = t1 - t0;
        if (bad_writes && g_writes - w0 > 1) (*bad_writes)++;
        if (bad_writes && g_writes == w0) g_silent++;
    }
}

int main(void) {
    symtab_init();
    core_init();

    TermEditorConfig cfg;
    memset(&cfg, 0, sizeof(cfg));
    cfg.cols = 120;
    cfg.rows = 40;
    cfg.color = true;
    cfg.classify = term_hl_classify_symtab;
    cfg.write = sink;
    TermEditor* ed = term_editor_new(&cfg);
    ASSERT(ed != NULL);
    term_editor_begin(ed, "In[1]:= ", 39, 0);

    /* 200 lines, pasted (one batch). */
    static char paste[64 * 1024];
    size_t pn = 0;
    pn += (size_t)snprintf(paste + pn, sizeof(paste) - pn, "\x1b[200~");
    for (int i = 0; i < 200; i++)
        pn += (size_t)snprintf(paste + pn, sizeof(paste) - pn,
                               "f%d[x_, y_] := Module[{t = Sin[x] + %d y^2}, t (* step %d *)]%s",
                               i, i, i, i < 199 ? ",\n" : "");
    pn += (size_t)snprintf(paste + pn, sizeof(paste) - pn, "\x1b[201~");
    ASSERT(term_editor_feed(ed, paste, pn) == TERM_ED_CONTINUE);
    for (int i = 0; i < 100; i++) term_editor_feed(ed, "\x1b[A", 3);   /* to the middle */

    /* Warm up until the undo log has reached its cap (then it is compacted in
     * place, never grown). */
    run_keys(ed, 3 * 2048, NULL, NULL);

    size_t a0 = term_editor_alloc_count();
    size_t bad_writes = 0;
    g_bytes = 0;
    run_keys(ed, NKEYS, g_t, &bad_writes);
    size_t allocs = term_editor_alloc_count() - a0;
    size_t bytes = g_bytes;

    qsort(g_t, NKEYS, sizeof(double), cmp_double);
    double sum = 0;
    for (size_t i = 0; i < NKEYS; i++) sum += g_t[i];
    double mean = sum / NKEYS, p50 = g_t[NKEYS / 2], p99 = g_t[NKEYS * 99 / 100], mx = g_t[NKEYS - 1];

    /* Worst case: a quote on the first line relexes all 200 lines (twice:
     * insert, then delete). */
    term_editor_feed(ed, "\x1b[5~", 4);                /* PageUp: buffer start */
    double worst = 0;
    for (int i = 0; i < 20; i++) {
        double t0 = now_us();
        term_editor_feed(ed, "\"", 1);
        double t1 = now_us();
        term_editor_feed(ed, "\x7f", 1);
        double t2 = now_us();
        if (t1 - t0 > worst) worst = t1 - t0;
        if (t2 - t1 > worst) worst = t2 - t1;
    }

    printf("latency (200-line buffer, 120x40, colour, %d keystrokes):\n", NKEYS);
    printf("  mean %.1f us  p50 %.1f us  p99 %.1f us  max %.1f us\n", mean, p50, p99, mx);
    printf("  whole-buffer relex (a quote at the top): max %.1f us\n", worst);
    printf("  keystrokes with more than one write: %zu  with none (no change): %zu\n",
           bad_writes, g_silent);
    printf("  bytes/keystroke %.1f  allocations %zu\n", (double)bytes / NKEYS, allocs);

    ASSERT_MSG(mx < 2000.0, "slowest keystroke %.1f us >= 2 ms", mx);
    ASSERT_MSG(worst < 2000.0, "whole-buffer relex %.1f us >= 2 ms", worst);
    ASSERT_MSG(bad_writes == 0, "%zu keystrokes made more than one write", bad_writes);
    ASSERT_MSG(allocs == 0, "%zu allocations in steady-state editing", allocs);

    term_editor_free(ed);
    printf("term_editor latency gate passed.\n");
    return 0;
}
