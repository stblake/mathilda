/* src/linalg/mp_parallel.c -- the worker-thread split shared by the
 * arbitrary-precision dense kernels (Dot, LU).  See mp_parallel.h. */
#ifndef _POSIX_C_SOURCE
#define _POSIX_C_SOURCE 200809L        /* sysconf */
#endif
#if defined(__APPLE__) && !defined(_DARWIN_C_SOURCE)
#define _DARWIN_C_SOURCE 1             /* _SC_NPROCESSORS_ONLN under _POSIX_C_SOURCE */
#endif

#include "mp_parallel.h"
#include "core.h"

#include <stdbool.h>
#include <stddef.h>

#ifdef USE_MPFR
#include <gmp.h>
#include <mpfr.h>
#endif

#if defined(MATHILDA_THREADS) && defined(USE_MPFR)
#include <pthread.h>
#include <signal.h>
#include <unistd.h>
#ifdef USE_FLINT
#include <flint/flint.h>
#endif

typedef struct {
    MpRangeFn fn;
    void* ctx;
    int64_t lo, hi;
    int chunk;
} MpJob;

static void* mp_job_run(void* p)
{
    MpJob* j = (MpJob*)p;
    j->fn(j->ctx, j->lo, j->hi, j->chunk);
#ifdef USE_FLINT
    flint_cleanup();                   /* this thread's FLINT caches */
#endif
    return NULL;
}

static int mp_nthreads(int64_t n, double work)
{
    if (n < 2 || work < 4194304.0 || !mpfr_buildopt_tls_p()) return 1;
    long ncpu = sysconf(_SC_NPROCESSORS_ONLN);
    if (ncpu < 1) ncpu = 1;
    if (ncpu > MP_PAR_MAX) ncpu = MP_PAR_MAX;
    double by_work = work / 1048576.0;
    if ((double)ncpu > by_work) ncpu = (long)by_work;
    if (ncpu > n) ncpu = (long)n;
    return ncpu < 1 ? 1 : (int)ncpu;
}
#endif

int mp_parallel_for(MpRangeFn fn, void* ctx, int64_t n, double work)
{
    if (n <= 0) return 0;
#if defined(MATHILDA_THREADS) && defined(USE_MPFR)
    int nt = mp_nthreads(n, work);
    if (nt > 1) {
        MpJob job[MP_PAR_MAX];
        pthread_t th[MP_PAR_MAX];
        int started[MP_PAR_MAX];
        int64_t step = (n + nt - 1) / nt;
        int nj = 0;
        for (int t = 0; t < nt; t++) {
            int64_t lo = (int64_t)t * step;
            if (lo >= n) break;
            job[nj].fn = fn;
            job[nj].ctx = ctx;
            job[nj].lo = lo;
            job[nj].hi = lo + step < n ? lo + step : n;
            job[nj].chunk = nj;
            nj++;
        }
        sigset_t all, old;
        sigfillset(&all);
        pthread_sigmask(SIG_SETMASK, &all, &old);
        tc_async_defer_push();
        /* a worker's stack as large as the main thread's would be: FLINT's
         * kernels keep their temporaries on the stack (TMP_ALLOC, 8 KB a
         * block, many blocks live at once), and the platform default for a
         * secondary thread (512 KB on macOS) is overrun by a heap division
         * with a few hundred divisors -- the multi-modular Groebner
         * verification spun in the stack probe instead of finishing */
        pthread_attr_t attr;
        int have_attr = pthread_attr_init(&attr) == 0;
        if (have_attr) pthread_attr_setstacksize(&attr, (size_t)64 << 20);
        for (int i = 0; i < nj; i++)
            started[i] = pthread_create(&th[i], have_attr ? &attr : NULL, mp_job_run, &job[i]) == 0;
        if (have_attr) pthread_attr_destroy(&attr);
        pthread_sigmask(SIG_SETMASK, &old, NULL);
        for (int i = 0; i < nj; i++) if (started[i]) pthread_join(th[i], NULL);
        for (int i = 0; i < nj; i++) if (!started[i]) fn(ctx, job[i].lo, job[i].hi, i);
        tc_async_defer_pop();
        return nj;
    }
#else
    (void)work;
#endif
    fn(ctx, 0, n, 0);
    return 1;
}
