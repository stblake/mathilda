/*
 * term_progress.c -- see term_progress.h.
 *
 * In libmathilda.a (no readline, no terminal): the iterators that write the
 * stack are kernel code, so the stack must link into every build.
 */
#include "term_progress.h"

#include "core.h"   /* tc_register_unwind_counter */

#include <math.h>
#include <limits.h>
#include <string.h>
#ifdef MATHILDA_THREADS
#include <pthread.h>
#endif

MthProgressSlot g_mth_progress[MTH_PROGRESS_SLOTS];
int g_mth_progress_depth = 0;
unsigned long g_mth_progress_gen = 0;
bool g_mth_progress_registered = false;

void term_progress_init(void) {
    if (g_mth_progress_registered) return;
    g_mth_progress_registered = true;
    tc_register_unwind_counter(&g_mth_progress_depth);
}

long mth_progress_range_count(double imin, double imax, double di) {
    if (!(di != 0.0) || !isfinite(imin) || !isfinite(imax) || !isfinite(di)) return 0;
    double q = (imax - imin) / di;
    if (!(q >= 0.0)) return 0;
    q = floor(q + 1e-9) + 1.0;
    if (q >= (double)LONG_MAX) return 0;
    return (long)q;
}

int term_progress_sample(MthProgressView* outer, MthProgressView* inner) {
    int d = *(volatile int*)&g_mth_progress_depth;
    if (d <= 0) return 0;
    if (d > MTH_PROGRESS_SLOTS) d = MTH_PROGRESS_SLOTS;
    const MthProgressSlot* o = &g_mth_progress[0];
    const MthProgressSlot* i = &g_mth_progress[d >= 2 ? 1 : 0];
    if (outer) {
        outer->gen = o->gen; outer->head = o->head;
        outer->done = o->done; outer->total = o->total;
    }
    if (inner) {
        inner->gen = i->gen; inner->head = i->head;
        inner->done = i->done; inner->total = i->total;
    }
    return d;
}

/* ---- worker threads and the memory budget ------------------------------ */

static volatile int g_workers = 0;
static volatile uint64_t g_mem_limit = 0;
#ifdef MATHILDA_THREADS
static pthread_mutex_t g_workers_mu = PTHREAD_MUTEX_INITIALIZER;
#endif

void mth_progress_workers_add(int delta) {
#ifdef MATHILDA_THREADS
    pthread_mutex_lock(&g_workers_mu);
#endif
    int w = g_workers + delta;
    g_workers = w < 0 ? 0 : w;
#ifdef MATHILDA_THREADS
    pthread_mutex_unlock(&g_workers_mu);
#endif
}

int mth_progress_workers(void) { return g_workers; }

void mth_progress_set_memory_limit(uint64_t bytes) { g_mem_limit = bytes; }
uint64_t mth_progress_memory_limit(void) { return g_mem_limit; }

/* ---- Monitor[expr, mon] ---------------------------------------------------- */

volatile int g_mth_monitor_due = 0;
static volatile int g_monitor_active = 0;
static char g_monitor_text[MTH_MONITOR_TEXT];
#ifdef MATHILDA_THREADS
static pthread_mutex_t g_monitor_mu = PTHREAD_MUTEX_INITIALIZER;
#endif

void mth_monitor_publish(const char* text) {
#ifdef MATHILDA_THREADS
    pthread_mutex_lock(&g_monitor_mu);
#endif
    size_t n = 0;
    if (text) {
        n = strlen(text);
        if (n >= sizeof g_monitor_text) {
            n = sizeof g_monitor_text - 1;
            while (n > 0 && ((unsigned char)text[n] & 0xC0) == 0x80) n--;   /* whole UTF-8 */
        }
        memcpy(g_monitor_text, text, n);
    }
    g_monitor_text[n] = '\0';
#ifdef MATHILDA_THREADS
    pthread_mutex_unlock(&g_monitor_mu);
#endif
}

bool mth_monitor_text(char* buf, size_t cap) {
    if (!buf || cap == 0) return false;
    buf[0] = '\0';
    if (!g_monitor_active) return false;
#ifdef MATHILDA_THREADS
    pthread_mutex_lock(&g_monitor_mu);
#endif
    size_t n = strlen(g_monitor_text);
    if (n >= cap) n = cap - 1;
    memcpy(buf, g_monitor_text, n);
    buf[n] = '\0';
#ifdef MATHILDA_THREADS
    pthread_mutex_unlock(&g_monitor_mu);
#endif
    return n > 0;
}

void mth_monitor_set_active(int depth) {
    g_monitor_active = depth > 0;
    if (depth <= 0) { g_mth_monitor_due = 0; mth_monitor_publish(NULL); }
}

bool mth_monitor_active(void) { return g_monitor_active != 0; }
