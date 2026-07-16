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

/* dgetrf -- LU factorisation with partial row pivoting.
 *   A (M*N, column-major)         input/output (L unit-lower, U upper)
 *   ipiv (min(M,N), 1-indexed)    row pivots; row i was swapped with ipiv[i]
 */
void dgetrf_(const int* m, const int* n, double* a, const int* lda,
             int* ipiv, int* info);

/* dpotrf -- Cholesky factorisation of a real symmetric positive-definite
 * matrix A = U^T U (UPLO='U') or A = L L^T (UPLO='L').  Reads / writes
 * only the selected triangle; the opposite triangle is untouched.
 *   info: 0 success, <0 bad arg, >0 leading minor not positive definite. */
void dpotrf_(const char* uplo, const int* n, double* a, const int* lda,
             int* info);

/* dlange -- 1/inf/Frobenius norm of a real matrix.  `norm` is one of
 * '1', 'O' (1-norm), 'I' (inf-norm), 'F', 'E' (Frobenius), 'M' (max abs). */
double dlange_(const char* norm, const int* m, const int* n,
               const double* a, const int* lda, double* work);

/* dgecon -- reciprocal-condition estimate from a dgetrf factorisation. */
void dgecon_(const char* norm, const int* n, const double* a,
             const int* lda, const double* anorm, double* rcond,
             double* work, int* iwork, int* info);

/* Complex variants. */
void zgetrf_(const int* m, const int* n, double* a, const int* lda,
             int* ipiv, int* info);

/* zpotrf -- Cholesky factorisation of a complex Hermitian positive-definite
 * matrix A = U^H U (UPLO='U') or A = L L^H (UPLO='L'). */
void zpotrf_(const char* uplo, const int* n, double* a, const int* lda,
             int* info);

double zlange_(const char* norm, const int* m, const int* n,
               const double* a, const int* lda, double* work);

void zgecon_(const char* norm, const int* n, const double* a,
             const int* lda, const double* anorm, double* rcond,
             double* work, double* rwork, int* info);

/* dgesdd -- divide-and-conquer SVD of a real m x n matrix.
 *   jobz     'A' = all m cols of U and all n rows of V^T (square U, V^T)
 *            'S' = "thin": first min(m,n) cols of U / rows of V^T
 *            'O' = overwrite A with U or V^T
 *            'N' = no singular vectors
 *   A (m*n, column-major)         destroyed on exit
 *   s (min(m,n))                  singular values, descending
 *   U (m*m or m*min(m,n))         left singular vectors, column-major
 *   VT (n*n or min(m,n)*n)        V^H rows, column-major
 *   work / lwork / iwork          workspace */
void dgesdd_(const char* jobz, const int* m, const int* n,
             double* A, const int* lda,
             double* s,
             double* U, const int* ldu,
             double* VT, const int* ldvt,
             double* work, const int* lwork, int* iwork, int* info);

/* zgesdd -- divide-and-conquer SVD of a complex m x n matrix.
 * A, U, VT carry (re, im) pairs.  rwork is real, length depends on
 * (jobz, min(m,n)); LAPACK documents the formula. */
void zgesdd_(const char* jobz, const int* m, const int* n,
             double* A, const int* lda,
             double* s,
             double* U, const int* ldu,
             double* VT, const int* ldvt,
             double* work, const int* lwork,
             double* rwork, int* iwork, int* info);

/* dggsvd3 -- generalized SVD of two real matrices A (m x n) and B (p x n).
 * Returns matrices U (m x m), V (p x p), Q (n x n), and integer
 * partition (k, l) of the singular values into alpha (length n) and
 * beta (length n).  See LAPACK documentation for the precise output
 * structure. */
void dggsvd3_(const char* jobu, const char* jobv, const char* jobq,
              const int* m, const int* n, const int* p,
              int* k_out, int* l_out,
              double* A, const int* lda, double* B, const int* ldb,
              double* alpha, double* beta,
              double* U, const int* ldu,
              double* V, const int* ldv,
              double* Q, const int* ldq,
              double* work, const int* lwork, int* iwork, int* info);

/* zggsvd3 -- complex variant of dggsvd3. */
void zggsvd3_(const char* jobu, const char* jobv, const char* jobq,
              const int* m, const int* n, const int* p,
              int* k_out, int* l_out,
              double* A, const int* lda, double* B, const int* ldb,
              double* alpha, double* beta,
              double* U, const int* ldu,
              double* V, const int* ldv,
              double* Q, const int* ldq,
              double* work, const int* lwork,
              double* rwork, int* iwork, int* info);

/* dgesv / zgesv -- solve the square system A X = B by LU factorisation.
 *   A (n*n, col-major)  overwritten with its LU factors
 *   ipiv (n)            row pivots
 *   B (n*nrhs, col-major) overwritten with the solution X */
void dgesv_(const int* n, const int* nrhs, double* a, const int* lda,
            int* ipiv, double* b, const int* ldb, int* info);
void zgesv_(const int* n, const int* nrhs, double* a, const int* lda,
            int* ipiv, double* b, const int* ldb, int* info);

/* dgetri / zgetri -- inverse from an existing dgetrf/zgetrf LU factorisation. */
void dgetri_(const int* n, double* a, const int* lda, const int* ipiv,
             double* work, const int* lwork, int* info);
void zgetri_(const int* n, double* a, const int* lda, const int* ipiv,
             double* work, const int* lwork, int* info);

/* dgeev / zgeev -- eigenvalues and (right) eigenvectors of a general matrix.
 * Real: eigenvalues arrive split as (wr, wi); complex conjugate pairs share a
 * packed VR column (LAPACK convention). Complex: eigenvalues in w, VR complex.*/
void dgeev_(const char* jobvl, const char* jobvr, const int* n,
            double* a, const int* lda, double* wr, double* wi,
            double* vl, const int* ldvl, double* vr, const int* ldvr,
            double* work, const int* lwork, int* info);
void zgeev_(const char* jobvl, const char* jobvr, const int* n,
            double* a, const int* lda, double* w,
            double* vl, const int* ldvl, double* vr, const int* ldvr,
            double* work, const int* lwork, double* rwork, int* info);

/* dsyev / zheev -- eigenvalues (ascending) and eigenvectors of a
 * symmetric / Hermitian matrix; A is overwritten with the eigenvectors. */
void dsyev_(const char* jobz, const char* uplo, const int* n,
            double* a, const int* lda, double* w,
            double* work, const int* lwork, int* info);
void zheev_(const char* jobz, const char* uplo, const int* n,
            double* a, const int* lda, double* w,
            double* work, const int* lwork, double* rwork, int* info);

/* dgels / zgels -- least-squares / minimum-norm solve of A X = B (full rank).
 * B is max(m,n)*nrhs, col-major; the solution occupies its first n rows. */
void dgels_(const char* trans, const int* m, const int* n, const int* nrhs,
            double* a, const int* lda, double* b, const int* ldb,
            double* work, const int* lwork, int* info);
void zgels_(const char* trans, const int* m, const int* n, const int* nrhs,
            double* a, const int* lda, double* b, const int* ldb,
            double* work, const int* lwork, int* info);

/* dgesvd / zgesvd -- classic (QR-iteration) SVD.  jobu = jobvt = 'A' returns
 * full square U and V^H.  Complex variant needs rwork of length 5*min(m,n). */
void dgesvd_(const char* jobu, const char* jobvt, const int* m, const int* n,
             double* a, const int* lda, double* s,
             double* u, const int* ldu, double* vt, const int* ldvt,
             double* work, const int* lwork, int* info);
void zgesvd_(const char* jobu, const char* jobvt, const int* m, const int* n,
             double* a, const int* lda, double* s,
             double* u, const int* ldu, double* vt, const int* ldvt,
             double* work, const int* lwork, double* rwork, int* info);

/* dtrtrs / ztrtrs -- solve a triangular system A X = B (no workspace). */
void dtrtrs_(const char* uplo, const char* trans, const char* diag,
             const int* n, const int* nrhs, const double* a, const int* lda,
             double* b, const int* ldb, int* info);
void ztrtrs_(const char* uplo, const char* trans, const char* diag,
             const int* n, const int* nrhs, const double* a, const int* lda,
             double* b, const int* ldb, int* info);

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

/* LU factorisation with partial pivoting.  Caller owns A (column-major,
 * lda*n doubles for real, 2*lda*n for complex) and ipiv (min(m,n) ints).
 * On return A holds Doolittle's L (strict lower, unit diag) and U (upper)
 * superimposed; ipiv[i] is the 1-indexed row that was swapped to row i. */
int mat_lapack_dgetrf(int m, int n, double* A, int lda, int* ipiv);
int mat_lapack_zgetrf(int m, int n, double* A, int lda, int* ipiv);

/* Cholesky factorisation.  `uplo` is 'U' or 'L' selecting which triangle
 * of A is referenced.  On success A's selected triangle holds the Cholesky
 * factor; the opposite triangle is untouched.  Returns LAPACK info:
 *   0    success (A is positive definite),
 *   <0   bad argument,
 *   >0   k means the k-th leading minor is not positive definite.
 * When USE_LAPACK is undefined the stubs return -1 so call sites can
 * route to an in-house Cholesky fallback. */
int mat_lapack_dpotrf(char uplo, int n, double* A, int lda);
int mat_lapack_zpotrf(char uplo, int n, double* A, int lda);

/* L-infinity norm of an m x n matrix.  `norm_kind` selects '1', 'I',
 * 'F', 'M'.  Returns the norm value (always >= 0).  When USE_LAPACK
 * is off the stub returns -1.0 so callers can detect the no-LAPACK case. */
double mat_lapack_dlange(char norm_kind, int m, int n,
                         const double* A, int lda);
double mat_lapack_zlange(char norm_kind, int m, int n,
                         const double* A, int lda);

/* Reciprocal condition-number estimate using an LU factorisation.
 * Caller provides the original-matrix norm `anorm` (typically computed
 * via mat_lapack_dlange / mat_lapack_zlange with the matching kind).
 * On success `*rcond` is in (0, 1] and approximates 1 / cond(A); on
 * failure (non-zero info) `*rcond` is untouched.  Returns LAPACK info
 * (0 success, <0 bad arg, >0 numerical issue). */
int mat_lapack_dgecon(char norm_kind, int n, const double* A_LU,
                      int lda, double anorm, double* rcond);
int mat_lapack_zgecon(char norm_kind, int n, const double* A_LU,
                      int lda, double anorm, double* rcond);

/* Divide-and-conquer SVD wrappers.  `jobz` is 'A' (full square U/V^T),
 * 'S' (thin), 'N' (no vectors), or 'O' (overwrite A).  A is destroyed
 * on exit (callers must keep a copy if they need m later).  S has
 * length min(m, n); U has dimensions m x m ('A') or m x min(m,n) ('S');
 * VT has dimensions n x n ('A') or min(m,n) x n ('S').
 *
 * Complex (z) variants use interleaved (re, im) doubles in A, U, VT
 * and emit U / V^T entries in the same layout.
 *
 * Returns LAPACK info (0 success, <0 bad arg, >0 dbdsdc didn't converge).
 * Stubs return -1 when USE_LAPACK is undefined. */
int mat_lapack_dgesdd(char jobz, int m, int n, double* A, int lda,
                      double* S, double* U, int ldu,
                      double* VT, int ldvt);
int mat_lapack_zgesdd(char jobz, int m, int n, double* A, int lda,
                      double* S, double* U, int ldu,
                      double* VT, int ldvt);

/* Generalized SVD wrappers.  See LAPACK dggsvd3 / zggsvd3 for the
 * precise output layout (we forward k, l, alpha, beta, U, V, Q
 * unchanged).  A and B are destroyed on exit.  alpha and beta have
 * length n; together they encode the generalized singular values
 * (alpha[i] / beta[i] for those indices where beta > 0).
 *
 * `jobu`, `jobv`, `jobq` select which orthogonal/unitary factors are
 * computed (typically 'U', 'V', 'Q' for "yes" or 'N' for "no").
 *
 * Returns LAPACK info; stubs return -1 when USE_LAPACK is undefined. */
int mat_lapack_dggsvd3(char jobu, char jobv, char jobq,
                       int m, int n, int p,
                       int* k_out, int* l_out,
                       double* A, int lda, double* B, int ldb,
                       double* alpha, double* beta,
                       double* U, int ldu,
                       double* V, int ldv,
                       double* Q, int ldq);
int mat_lapack_zggsvd3(char jobu, char jobv, char jobq,
                       int m, int n, int p,
                       int* k_out, int* l_out,
                       double* A, int lda, double* B, int ldb,
                       double* alpha, double* beta,
                       double* U, int ldu,
                       double* V, int ldv,
                       double* Q, int ldq);

/* ---------------------------------------------------------------------
 * Higher-level solve / eigen / SVD wrappers used by lapack_bridge.c.
 * All matrices are column-major; complex variants carry interleaved
 * (re, im) doubles. Each returns the LAPACK info code (0 success, <0 bad
 * arg, >0 numerical issue); stubs return -1 without LAPACK.
 * ------------------------------------------------------------------ */

/* Solve A X = B (A n*n, B n*nrhs). On success A holds the LU factors, ipiv the
 * pivots, and B is overwritten with the solution X. */
int mat_lapack_dgesv(int n, int nrhs, double* A, int lda, int* ipiv,
                     double* B, int ldb);
int mat_lapack_zgesv(int n, int nrhs, double* A, int lda, int* ipiv,
                     double* B, int ldb);

/* Invert an already-factored matrix (dgetrf/zgetrf output in A + ipiv). */
int mat_lapack_dgetri(int n, double* A, int lda, const int* ipiv);
int mat_lapack_zgetri(int n, double* A, int lda, const int* ipiv);

/* General eigenproblem: eigenvalues and right eigenvectors (jobvl='N',
 * jobvr='V'). Real: (wr, wi) split, VR packed per LAPACK's real convention.
 * Complex: eigenvalues in w (2n interleaved), VR complex n*n. A is destroyed. */
int mat_lapack_dgeev(int n, double* A, int lda,
                     double* wr, double* wi, double* VR, int ldvr);
int mat_lapack_zgeev(int n, double* A, int lda,
                     double* w, double* VR, int ldvr);

/* Symmetric / Hermitian eigenproblem (jobz='V', uplo='U'): eigenvalues
 * ascending in w (length n, real), eigenvectors overwrite A (columns). */
int mat_lapack_dsyev(int n, double* A, int lda, double* w);
int mat_lapack_zheev(int n, double* A, int lda, double* w);

/* Least squares min ||A X - B|| (trans='N', full rank). A is m*n; B is
 * ldb*nrhs with ldb = max(m,n); the solution occupies its first n rows. */
int mat_lapack_dgels(int m, int n, int nrhs, double* A, int lda,
                     double* B, int ldb);
int mat_lapack_zgels(int m, int n, int nrhs, double* A, int lda,
                     double* B, int ldb);

/* Classic SVD (jobu=jobvt='A'): full square U (m*m) and VT (n*n), S length
 * min(m,n). A is destroyed. */
int mat_lapack_dgesvd(int m, int n, double* A, int lda,
                      double* S, double* U, int ldu, double* VT, int ldvt);
int mat_lapack_zgesvd(int m, int n, double* A, int lda,
                      double* S, double* U, int ldu, double* VT, int ldvt);

/* Triangular solve A X = B (uplo='U', trans='N', diag='N'); B overwritten. */
int mat_lapack_dtrtrs(int n, int nrhs, const double* A, int lda,
                      double* B, int ldb);
int mat_lapack_ztrtrs(int n, int nrhs, const double* A, int lda,
                      double* B, int ldb);

#ifdef __cplusplus
}
#endif

#endif /* MATHILDA_LINALG_LAPACK_H */
