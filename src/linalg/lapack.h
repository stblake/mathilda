/*
 * src/linalg/lapack.h
 *
 * Platform-papering wrapper for BLAS and LAPACK.
 *
 * Mathilda uses BLAS/LAPACK for the machine-precision linear-algebra
 * kernels (QRDecomposition, and eventually Inverse, LinearSolve,
 * Det, Eigenvalues, LeastSquares, ...).  We deliberately speak the
 * Fortran ABI directly rather than depending on LAPACKE because:
 *
 *   - Apple's Accelerate framework ships the Fortran ABI plus CBLAS
 *     but does NOT ship LAPACKE row-major wrappers, and requiring a
 *     Homebrew lapack package would defeat the zero-install property
 *     of Accelerate on macOS.
 *   - Every LAPACK distribution (Accelerate, OpenBLAS, MKL, Netlib)
 *     exposes the same Fortran ABI symbols (`dgeqp3_`, etc.).
 *   - The row-major <-> column-major transpose is cheap and lives
 *     once inside our load/store helpers.
 *
 * Build-time policy: this header is safe to include unconditionally.
 * When USE_LAPACK is undefined the body collapses, and a stub
 * `mathilda_lapack_probe()` is provided so callers compile cleanly
 * on systems without BLAS/LAPACK (the makefile auto-degrades to
 * USE_LAPACK=0 when no library is found).
 */

#ifndef MATHILDA_LINALG_LAPACK_H
#define MATHILDA_LINALG_LAPACK_H

#ifdef __cplusplus
extern "C" {
#endif

#ifdef USE_LAPACK

#if defined(MATHILDA_USE_ACCELERATE)
/*
 * Apple Accelerate.  Opt in to the modern LAPACK interface
 * (still Fortran ABI, but with the `__LAPACK_int` typedef which
 * lets Accelerate use a wider int on Apple Silicon if it ever
 * wants to).  Defining ACCELERATE_NEW_LAPACK before including
 * <Accelerate/Accelerate.h> is the documented way to do this.
 * Available since macOS 13.3.
 */
#  ifndef ACCELERATE_NEW_LAPACK
#    define ACCELERATE_NEW_LAPACK
#  endif
#  include <Accelerate/Accelerate.h>
#else
/*
 * Linux / *BSD / Solaris with a system BLAS/LAPACK.  cblas.h is
 * universally provided; LAPACK symbols are declared `extern`
 * below in the Fortran ABI, so we do not require <lapacke.h>.
 */
#  include <cblas.h>
#endif

#include <stdint.h>

/*
 * Fortran ABI LAPACK prototypes.
 *
 * These are the routines we currently call from Phase 3 (machine
 * QR kernel).  On Apple Accelerate they are already declared in
 * <Accelerate/Accelerate.h> via <vecLib/lapack.h>, so the
 * declarations are guarded out on macOS to avoid redeclaration
 * conflicts.
 *
 * Calling convention (all platforms with USE_LAPACK):
 *   - all scalar arguments passed by pointer
 *   - matrices in column-major order
 *   - `info` parameter receives the LAPACK status code
 *     (0 = success, <0 = bad argument index, >0 = numerical issue)
 *
 * See https://netlib.org/lapack/ for full specification.
 */
#if !defined(MATHILDA_USE_ACCELERATE)

/* dgeqrf — compute QR factorization of a real M x N matrix
 *   A (M*N, column-major)         input/output
 *   tau (min(M,N))                output reflector scalars
 *   work (lwork)                  workspace
 */
void dgeqrf_(const int* m, const int* n, double* a, const int* lda,
             double* tau, double* work, const int* lwork, int* info);

/* dgeqp3 — QR factorization with column pivoting */
void dgeqp3_(const int* m, const int* n, double* a, const int* lda,
             int* jpvt, double* tau, double* work, const int* lwork,
             int* info);

/* dorgqr — form Q (M x K) from reflectors produced by dgeqrf/dgeqp3 */
void dorgqr_(const int* m, const int* n, const int* k, double* a,
             const int* lda, const double* tau, double* work,
             const int* lwork, int* info);

/* Complex variants (real and imaginary stored interleaved per element). */
void zgeqrf_(const int* m, const int* n, double* a, const int* lda,
             double* tau, double* work, const int* lwork, int* info);

void zgeqp3_(const int* m, const int* n, double* a, const int* lda,
             int* jpvt, double* tau, double* work, const int* lwork,
             double* rwork, int* info);

void zungqr_(const int* m, const int* n, const int* k, double* a,
             const int* lda, const double* tau, double* work,
             const int* lwork, int* info);

#endif /* !MATHILDA_USE_ACCELERATE */

#endif /* USE_LAPACK */

/*
 * Phase-1 probe.  Returns:
 *   1 when USE_LAPACK is defined AND the linked BLAS produces the
 *     correct answer for a trivial cblas_ddot call.
 *   0 when USE_LAPACK is undefined, OR the BLAS gave a wrong answer
 *     (which would indicate a misconfigured / corrupted distribution).
 *
 * Used by tests/test_lapack.c and by future linalg kernels that
 * want to gate LAPACK-only code paths at runtime.
 */
int mathilda_lapack_probe(void);

/* ---------------------------------------------------------------------
 * Thin C wrappers around the Fortran-ABI LAPACK QR routines.
 *
 *   - Hide the pointer-to-int Fortran calling convention.
 *   - Run a workspace query, allocate, dispatch, and free, so kernels
 *     only see a single function call.
 *   - Forward LAPACK's `info` code unchanged (0 success, <0 bad arg,
 *     >0 numerical issue) so the kernel can decide to fall back to
 *     the symbolic path.
 *
 * Matrices are passed in **column-major** layout (the native Fortran
 * convention).  Caller owns A, jpvt, tau buffers; the wrappers only
 * borrow them for the duration of the call.
 *
 * Complex variants (z*) take A and tau as flat double arrays with
 * interleaved (real, imag) pairs per element, matching the LAPACK
 * `complex*16` ABI.  A occupies `2 * lda * n` doubles; tau occupies
 * `2 * min(m, n)` doubles.
 *
 * All wrappers are stubs that return a negative info code when
 * USE_LAPACK is undefined, so call sites compile and link cleanly
 * without LAPACK.  The kernel detects the stub case and falls back to
 * the symbolic pipeline.
 * ------------------------------------------------------------------ */
int mat_lapack_dgeqrf(int m, int n, double* A, int lda, double* tau);
int mat_lapack_dgeqp3(int m, int n, double* A, int lda,
                      int* jpvt, double* tau);
int mat_lapack_dorgqr(int m, int n, int k, double* A, int lda,
                      const double* tau);

int mat_lapack_zgeqrf(int m, int n, double* A, int lda, double* tau);
int mat_lapack_zgeqp3(int m, int n, double* A, int lda,
                      int* jpvt, double* tau);
int mat_lapack_zungqr(int m, int n, int k, double* A, int lda,
                      const double* tau);

#ifdef __cplusplus
}
#endif

#endif /* MATHILDA_LINALG_LAPACK_H */
