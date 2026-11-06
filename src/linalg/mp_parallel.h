#ifndef MATHILDA_LINALG_MP_PARALLEL_H
#define MATHILDA_LINALG_MP_PARALLEL_H

#include <stdint.h>

/* Most chunks mp_parallel_for splits a range into. */
#define MP_PAR_MAX 16

/* A chunk of work: the half-open index range [lo, hi), and the chunk's number
 * in [0, MP_PAR_MAX) so it can write a private slot of a per-chunk array. */
typedef void (*MpRangeFn)(void* ctx, int64_t lo, int64_t hi, int chunk);

/* Run fn over [0, n) in up to MP_PAR_MAX contiguous chunks on worker threads
 * when `work` (the job's cost in limb-products: multiply-adds times
 * (bits/64)^2) pays for them -- at least ~2^22 in all and ~2^20 per thread --
 * else serially as one chunk.  Returns the number of chunks used.
 *
 * For arbitrary-precision kernels (MPFR / FLINT arb numbers only: fn must
 * not touch an Expr or the evaluator).  The discipline of
 * src/numbertheory/factorinteger/fi_par.c: workers start with every signal
 * blocked, the region is bracketed by tc_async_defer_push/pop, the calling
 * thread only joins, each worker frees its FLINT caches on exit.  Serial when
 * MPFR is not built thread-safe or MATHILDA_THREADS is off. */
int mp_parallel_for(MpRangeFn fn, void* ctx, int64_t n, double work);

#endif
