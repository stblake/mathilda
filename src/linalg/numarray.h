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
 * Results come back PACKED for both real and complex data: a real result is a
 * float64 NDArray and a complex result is a complex64 NDArray (interleaved
 * re,im) -- NDArray is NOT real-only. `na_build_matrix` / `na_build_vector`
 * never box a per-element `Expr`; `na_scalar` still returns a bare
 * `Complex[re,im]` (im == 0 collapsed to real) for a rank-0 result.  Callers
 * that must match an input's representation stamp the result via
 * `na_build_*_as(..., na_result_presentation(input))`.
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
 * (as in na_load_matrix). Real -> float64 NDArray; complex -> complex64 NDArray
 * (both packed; NDArray is NOT real-only). im == 0 is NOT collapsed here. */
Expr* na_build_matrix(const double* buf, int rows, int cols, bool is_complex,
                      bool colmajor);

/* The presentation a result factor should carry so it stays packed (no
 * per-element Expr boxing) yet behaves correctly downstream: an NDArray input's
 * own present_as (a visible NDArray[...] input yields a visible result; a
 * transparent packed-list input stays transparent), else NDA_HEAD_LIST for a
 * boxed-List input -- a transparent packed-list that reads as a List and threads
 * correctly against any operand in reconstruction arithmetic. */
NDPresentation na_result_presentation(const Expr* input);

/* na_build_vector / na_build_matrix, then stamp present_as. Pass
 * na_result_presentation(input) so machine-numeric factors match the input's
 * representation without ever round-tripping through boxed Exprs. */
Expr* na_build_vector_as(const double* buf, int n, bool is_complex,
                         NDPresentation pres);
Expr* na_build_matrix_as(const double* buf, int rows, int cols, bool is_complex,
                         bool colmajor, NDPresentation pres);

#ifdef __cplusplus
}
#endif

#endif /* MATHILDA_LINALG_NUMARRAY_H */
