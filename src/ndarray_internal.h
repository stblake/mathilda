#ifndef MATHILDA_NDARRAY_INTERNAL_H
#define MATHILDA_NDARRAY_INTERNAL_H

/* ---------------------------------------------------------------------------
 * Buffer-level helpers shared between the NDArray fast-path modules
 * (ndreduce.c, ndstruct.c). Not part of the public ndarray.h surface — these
 * expose raw double scratch buffers and are only meaningful to code that
 * already works directly on the flat machine-precision buffer.
 * -------------------------------------------------------------------------- */

#include "expr.h"     /* NDType */
#include <stddef.h>
#include <stdbool.h>

/* ---------------------------------------------------------------------------
 * Parallelism. The element-wise map path (ndarray.c) already splits large flat
 * ranges across worker threads; these expose the same machinery to the
 * reduction/structural modules. Guarded internally by MATHILDA_THREADS — when
 * threads are unavailable every call runs serially (one chunk). The size
 * threshold and thread cap live in ndarray.c.
 * -------------------------------------------------------------------------- */

#define NDARRAY_MAX_THREADS 16

/* A chunk callback processes the half-open flat range [lo, hi), writing to
 * disjoint output indices (no cross-chunk dependency). Returns false to signal
 * the kernel declined an element (caller degrades). */
typedef bool (*nd_chunk_fn)(void* ctx, size_t lo, size_t hi);

/* Run `fn` over [0, n) — threaded for large n, serial otherwise. Returns true
 * iff every chunk returned true. Use when outputs are index-disjoint (maps,
 * columnwise reductions where each thread owns whole output columns). */
bool nd_parallel_for(size_t n, nd_chunk_fn fn, void* ctx);

/* Parallel reduction: split [0, n) into up to NDARRAY_MAX_THREADS chunks, each
 * folding its range into a private slot of `ncomp` doubles (no shared writes,
 * so no races); returns the number of chunks used. The caller combines the
 * filled slots (sum: add the slots; min/max: reduce with min/max). `slots` must
 * hold at least NDARRAY_MAX_THREADS * ncomp doubles. */
typedef void (*nd_reduce_chunk_fn)(void* ctx, size_t lo, size_t hi, double* slot);
int nd_parallel_reduce(size_t n, nd_reduce_chunk_fn fn, void* ctx,
                       int ncomp, double* slots);

/* Copy the real parts of the `count` elements buf[base], buf[base+stride], ...
 * into out[0..count). Complex imaginary parts are dropped (callers use this
 * only on real dtypes or where the imaginary part is irrelevant). */
void nd_gather_real(const void* buf, NDType dt, size_t base, size_t stride,
                    size_t count, double* out);

/* Quickselect: partition scratch[0..n) in place so scratch[k] holds the k-th
 * smallest value (0-based) and every earlier slot is <= it. Returns scratch[k].
 * O(n) average; mutates scratch. */
double nd_select_kth(double* scratch, size_t n, size_t k);

/* Sort scratch[0..n) ascending in place (plain qsort by value). */
void nd_sort_ascending(double* scratch, size_t n);

#endif /* MATHILDA_NDARRAY_INTERNAL_H */
