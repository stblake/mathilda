/*
 * src/linalg/numarray.h
 *
 * Shared, dtype-aware marshalling between Mathilda expressions and the flat
 * `double` buffers the BLAS/LAPACK bridges (blas_bridge.c, lapack_bridge.c)
 * hand to the numerical kernels.
 *
 * A "numeric array" argument is accepted in either of two shapes:
 *   - an `NDArray[...]` value (EXPR_NDARRAY, dense row-major real doubles), or
 *   - an ordinary nested `List` whose leaves are numeric (Integer, Real,
 *     BigInt, MPFR, `Rational[p,q]`, or `Complex[re,im]`).
 *
 * The routine being called fixes the element dtype (a `d*` routine wants real,
 * a `z*` routine wants complex), so the loaders take an explicit `want_complex`
 * flag rather than sniffing the data. Complex buffers store interleaved
 * (re, im) pairs per element, byte-identical to the LAPACK `complex*16` and
 * CBLAS `void`-pointer ABIs. A real-typed load of data that carries a nonzero
 * imaginary part fails (returns false) so the call is left unevaluated.
 *
 * Results come back as `NDArray[...]` for real data (consistent with the rest
 * of the numeric linear-algebra surface) and as a nested `List` of
 * `Complex[re,im]` for complex data, since NDArray is real-only. A complex
 * entry whose imaginary part is exactly zero is collapsed to a bare real,
 * matching the numericalize / machine-kernel convention.
 *
 * All loaders allocate the returned buffer with malloc; the caller owns it and
 * must free() it. On failure nothing is allocated and nothing leaks.
 */
#ifndef MATHILDA_LINALG_NUMARRAY_H
#define MATHILDA_LINALG_NUMARRAY_H

#include "expr.h"
#include <stdbool.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Read one numeric scalar leaf into (re, im). Recognises Integer, Real,
 * BigInt, MPFR, Rational[p,q], and Complex[re,im]. Returns false (leaving the
 * outputs untouched) for anything else. */
bool na_read_scalar(const Expr* e, double* re, double* im);

/* Load a rank-1 vector (an NDArray of rank 1, or a flat List of numeric
 * leaves) into a fresh double buffer. Real: n doubles. Complex: 2n interleaved
 * doubles. On success sets *n and *buf (caller frees *buf) and returns true;
 * returns false on wrong rank, an empty vector, a non-numeric leaf, or (when
 * want_complex is false) a leaf with a nonzero imaginary part. */
bool na_load_vector(const Expr* e, bool want_complex, int* n, double** buf);

/* Load a rank-2 matrix (an NDArray of rank 2, or a rectangular List-of-Lists)
 * into a fresh double buffer. `colmajor` selects the output layout:
 *   true  -> Fortran column-major, element (i,j) at index i + j*rows
 *   false -> C row-major,          element (i,j) at index i*cols + j
 * (each times 2, interleaved, when want_complex). On success sets *rows,
 * *cols, *buf (caller frees *buf) and returns true; false on wrong rank, an
 * empty/ragged matrix, a non-numeric leaf, or a real load of complex data. */
bool na_load_matrix(const Expr* e, bool want_complex, bool colmajor,
                    int* rows, int* cols, double** buf);

/* Build a scalar Expr: a bare Real when im == 0, else Complex[re, im]. */
Expr* na_scalar(double re, double im);

/* Build a length-n vector Expr from a flat buffer. Real -> NDArray of rank 1;
 * complex -> List of Complex[...] (im == 0 collapsed to Real). */
Expr* na_build_vector(const double* buf, int n, bool is_complex);

/* Build a rows x cols matrix Expr from a flat buffer laid out per `colmajor`
 * (as in na_load_matrix). Real -> NDArray of rank 2; complex -> nested List of
 * Complex[...] (im == 0 collapsed to Real). */
Expr* na_build_matrix(const double* buf, int rows, int cols, bool is_complex,
                      bool colmajor);

#ifdef __cplusplus
}
#endif

#endif /* MATHILDA_LINALG_NUMARRAY_H */
