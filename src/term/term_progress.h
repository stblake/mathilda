/*
 * term_progress.h -- the loop-progress stack the REPL spinner reads.
 *
 * Iterators (Do, Table, Sum, Product, and the Nest/Fold family through
 * iter_run) push a slot when their interpreted loop starts, store the number
 * of completed steps once per iteration, and pop when the loop ends. The
 * interactive REPL's status thread (term_status.c) samples the stack about
 * ten times a second to draw "Table 3400/10000 [bar] 34% . eta 6s".
 *
 * Cost: a push or pop is a handful of stores; a step is ONE store of a long,
 * next to an evaluate() call that costs 50 ns to many microseconds. Nothing
 * here allocates, locks or branches on whether anyone is reading. Numeric
 * fast paths (numloop_*, packed Table) are not instrumented.
 *
 * Threads (C99 has no memory model, so this is the documented contract):
 *   - Only the evaluating (main) thread writes the stack.
 *   - Every slot field is a naturally aligned, word-sized volatile object
 *     (long or a pointer), so a reader sees either the old or the new value,
 *     never a torn one, on every supported target.
 *   - A reader may see a slot whose fields come from two different loops for
 *     one sample; it only formats them, so the worst case is one odd frame.
 *     `gen` (bumped by every push) lets it notice a new loop.
 *   - `head` always points at an interned symbol name or a string literal;
 *     both live for the whole process, so a stale pointer is still valid.
 * ThreadSanitizer reports these reads by design.
 *
 * Unwinding: Abort / Throw / Return travel as ordinary return values, so the
 * pop runs. A TimeConstrained timeout siglongjmps past it; the depth is a
 * plain int registered with tc_register_unwind_counter (term_progress_init),
 * so the timeout restores it like every other nesting counter. pop(token)
 * sets the depth to the token, so even a missed pop is repaired by the
 * enclosing loop's pop.
 */
#ifndef MATHILDA_TERM_PROGRESS_H
#define MATHILDA_TERM_PROGRESS_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#define MTH_PROGRESS_SLOTS 16

typedef struct {
    const char* volatile head;   /* "Table", "Do", ... (immortal string)        */
    volatile long done;          /* completed steps                             */
    volatile long total;         /* number of steps, or 0 when unknown/infinite */
    volatile unsigned long gen;  /* push generation, unique per loop            */
} MthProgressSlot;

extern MthProgressSlot g_mth_progress[MTH_PROGRESS_SLOTS];
extern int g_mth_progress_depth;          /* may exceed MTH_PROGRESS_SLOTS */
extern unsigned long g_mth_progress_gen;
extern bool g_mth_progress_registered;

/* Register the depth with TimeConstrained's unwind (idempotent). Called by the
 * REPL at start-up and by the first push; call it from an embedder's init to
 * cover a timeout around the very first loop. */
void term_progress_init(void);

/* Start a loop: returns the token to pass to step/pop. total 0 = unknown. */
static inline int mth_progress_push(const char* head, long total) {
    if (!g_mth_progress_registered) term_progress_init();
    int d = g_mth_progress_depth;
    if (d >= 0 && d < MTH_PROGRESS_SLOTS) {
        MthProgressSlot* s = &g_mth_progress[d];
        s->done = 0;
        s->total = total < 0 ? 0 : total;
        s->head = head;
        s->gen = ++g_mth_progress_gen;
    }
    g_mth_progress_depth = d + 1;
    return d;
}

/* Monitor[expr, mon] (below): set by the status thread when mon is due for
 * a refresh; mth_monitor_tick (frontend.c) evaluates it on this thread. */
extern volatile int g_mth_monitor_due;
void mth_monitor_tick(void);

/* `done` steps of the loop `tok` have completed. A loop step is where loop
 * variables change, so a due Monitor is refreshed here too (one load). */
static inline void mth_progress_step(int tok, long done) {
    if (tok >= 0 && tok < MTH_PROGRESS_SLOTS) g_mth_progress[tok].done = done;
    if (g_mth_monitor_due) mth_monitor_tick();
}

/* End the loop `tok` (and anything a missed pop left above it). */
static inline void mth_progress_pop(int tok) {
    g_mth_progress_depth = tok < 0 ? 0 : tok;
}

/* Steps of a numeric range {imin, imax, di}, for a progress total: 0 when it
 * is empty, not finite, or too long for a long. */
long mth_progress_range_count(double imin, double imax, double di);

/* ---- reader side (the status thread) ---------------------------------- */

typedef struct {
    const char* head;
    long done, total;
    unsigned long gen;
} MthProgressView;

/* Copy the outermost live slot and the one directly inside it (racy, see
 * above). Returns the live depth, capped at MTH_PROGRESS_SLOTS; 0 means no
 * loop is running, and `inner` is meaningful only when it is >= 2. */
int term_progress_sample(MthProgressView* outer, MthProgressView* inner);

/* ---- worker threads and the memory budget ------------------------------ */

/* A parallel driver is about to start (+n) or has joined (-n) worker threads.
 * Locked; called once per parallel region, never per element. */
void mth_progress_workers_add(int delta);
int  mth_progress_workers(void);

/* The memory budget the spinner reports against (0 = the machine's RAM). For a
 * future MemoryConstrained; nothing sets it today. */
void     mth_progress_set_memory_limit(uint64_t bytes);
uint64_t mth_progress_memory_limit(void);

/* ---- Monitor[expr, mon] ------------------------------------------------
 *
 * While a Monitor is active and the status line is drawing, the status
 * thread sets g_mth_monitor_due about ten times a second. The evaluator
 * looks at it on its periodic check (every 128 evaluation steps, next to the
 * TimeConstrained deadline check, so the common path pays nothing) and, when
 * it is set, evaluates `mon` on the main thread and publishes its text here;
 * the status thread copies the text into the line. Loop steps
 * (mth_progress_step) check it as well. Without a status line nothing sets
 * the flag and `mon` is never evaluated. */
#define MTH_MONITOR_TEXT 256

void mth_monitor_publish(const char* text);      /* main thread; NULL clears   */
bool mth_monitor_text(char* buf, size_t cap);    /* any thread; false if none  */
void mth_monitor_set_active(int depth);          /* Monitor nesting depth      */
bool mth_monitor_active(void);

#endif /* MATHILDA_TERM_PROGRESS_H */
