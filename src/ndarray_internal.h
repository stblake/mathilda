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
