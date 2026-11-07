/*
 * src/linalg/lapack.c
 *
 * BLAS/LAPACK helpers.
 *
 * Two concerns live here:
 *
 *   1. `mathilda_lapack_probe()` -- the Phase-1 linkage probe.
 *      Returns 1 when the configured BLAS responds correctly to a
 *      `cblas_ddot` call; 0 otherwise (including USE_LAPACK off).
 *
 *   2. Thin C wrappers around the Fortran-ABI QR routines used by
 *      `src/linalg/qrdecomp_machine.c` (Phase 3) -- and, in follow-up
 *      phases, by the Inverse / LinearSolve / Det / Eigenvalues /
 *      LeastSquares fast paths.  Sharing the workspace-query +
 *      pointer-arg convention here keeps the kernel files free of
 *      Fortran-flavoured boilerplate.
 *
 * The wrappers always exist; under USE_LAPACK=0 they return a
 * negative info code so call sites link cleanly and can route to the
 * symbolic fallback at runtime.
 */

#include "lapack.h"

#include <math.h>
#include <stdbool.h>   /* issue #128: <cblas.h> does not supply it on Linux */
#include <stdlib.h>
#include <string.h>

int mathilda_lapack_probe(void)
{
#ifdef USE_LAPACK
    /*
     * cblas_ddot is the simplest possible BLAS call: dot product of
     * two strided double vectors.  We use two parallel unit vectors;
     * the answer should be exactly 1.0 for any sensible BLAS.
     *
     * If BLAS is misconfigured (e.g. mixed 32/64-bit int convention,
     * mismatched calling convention, corrupted library) cblas_ddot
     * will either crash or return garbage, and the probe reports 0.
     */
    const double x[2] = { 1.0, 0.0 };
    const double y[2] = { 1.0, 0.0 };
    const double r = cblas_ddot(2, x, 1, y, 1);
    return (r == 1.0) ? 1 : 0;
#else
    return 0;
#endif
}

/* ---------------------------------------------------------------------
 * LAPACK QR wrappers.
 *
 * The Fortran ABI requires every scalar to be passed by pointer, and
 * (with the exception of `dorgqr`/`zungqr`) most routines need a
 * dynamically-sized workspace whose preferred length is reported by a
 * preliminary "query" call with `lwork = -1`.  We bundle the query +
 * allocate + dispatch + free into one C call per routine.
 *
 * Conventions across all wrappers:
 *   - `m`, `n`, `k`, `lda` are passed as C `int` and forwarded as
 *     pointers (LAPACK Fortran integers are `int` on every platform
 *     we target -- Accelerate's `__LAPACK_int` typedef resolves to
 *     `int` under ACCELERATE_NEW_LAPACK).
 *   - The matrix `A` is column-major.  For complex routines it stores
 *     interleaved (re, im) pairs as a flat double array; the LAPACK
 *     `complex*16` layout is byte-identical.
 *   - `jpvt` follows LAPACK's 1-indexed convention.  Pivot routines
 *     initialise it to zero on entry (the documented "no pivot
 *     preference" sentinel) before forwarding.
 *
 * Each wrapper returns the LAPACK `info` code unchanged.  Workspace
 * allocation failure surfaces as `info = -999` so the kernel can fall
 * back to the symbolic path.
 * ------------------------------------------------------------------ */

#ifdef USE_LAPACK

/* Compute lwork from a returned workspace-query scalar.  Floors at the
 * LAPACK-documented minimum so under-counted distributions still work. */
static int lapack_lwork(double query, int floor_min)
{
    int w = (int)query;
    if (w < floor_min) w = floor_min;
    return w;
}

/* Cast a `double*` (interleaved re/im pairs) to whatever pointer type
 * the platform LAPACK declares its `complex*16` arguments as.  On
 * Apple Accelerate this is `_Complex double*`, which C99 guarantees is
 * layout-compatible with `double[2]` (§6.2.5/13).  Elsewhere we declare
 * our own Fortran-ABI prototypes with plain `double*`, so the cast is
 * a no-op.  Routed through a tiny inline helper so each call site is
 * a single, self-documenting cast and the warning suppression lives in
 * exactly one place. */
#if defined(MATHILDA_USE_ACCELERATE)
typedef _Complex double mathilda_lapack_zptr_t;
#else
typedef double mathilda_lapack_zptr_t;
#endif
static inline mathilda_lapack_zptr_t* lapack_zptr(double* p)
{
    return (mathilda_lapack_zptr_t*)(void*)p;
}

int mat_lapack_dgeqrf(int m, int n, double* A, int lda, double* tau)
{
    int info = 0, lwork = -1;
    double query = 0.0;
    dgeqrf_(&m, &n, A, &lda, tau, &query, &lwork, &info);
    if (info != 0) return info;
    lwork = lapack_lwork(query, n > 1 ? n : 1);
    double* work = (double*)malloc((size_t)lwork * sizeof(double));
    if (!work) return -999;
    dgeqrf_(&m, &n, A, &lda, tau, work, &lwork, &info);
    free(work);
    return info;
}

int mat_lapack_dgeqp3(int m, int n, double* A, int lda,
                      int* jpvt, double* tau)
{
    for (int j = 0; j < n; j++) jpvt[j] = 0;  /* "no pivot preference" */
    int info = 0, lwork = -1;
    double query = 0.0;
    dgeqp3_(&m, &n, A, &lda, jpvt, tau, &query, &lwork, &info);
    if (info != 0) return info;
    /* Documented minimum is 3*n+1; we floor at that and prefer the
     * library-reported optimum. */
    lwork = lapack_lwork(query, 3 * n + 1);
    double* work = (double*)malloc((size_t)lwork * sizeof(double));
    if (!work) return -999;
    dgeqp3_(&m, &n, A, &lda, jpvt, tau, work, &lwork, &info);
    free(work);
    return info;
}

int mat_lapack_dorgqr(int m, int n, int k, double* A, int lda,
                      const double* tau)
{
    int info = 0, lwork = -1;
    double query = 0.0;
    /* tau is read-only in dorgqr; LAPACK declares it `const`-equivalent
     * but the Fortran prototype isn't const-qualified, so we cast. */
    double* tau_nc = (double*)tau;
    dorgqr_(&m, &n, &k, A, &lda, tau_nc, &query, &lwork, &info);
    if (info != 0) return info;
    lwork = lapack_lwork(query, n > 1 ? n : 1);
    double* work = (double*)malloc((size_t)lwork * sizeof(double));
    if (!work) return -999;
    dorgqr_(&m, &n, &k, A, &lda, tau_nc, work, &lwork, &info);
    free(work);
    return info;
}

int mat_lapack_zgeqrf(int m, int n, double* A, int lda, double* tau)
{
    int info = 0, lwork = -1;
    /* zgeqrf's `work` is complex*16 -- a pair of doubles.  The
     * workspace-query result lives in the real part of work[0]. */
    double query[2] = { 0.0, 0.0 };
    zgeqrf_(&m, &n, lapack_zptr(A), &lda, lapack_zptr(tau),
            lapack_zptr(query), &lwork, &info);
    if (info != 0) return info;
    lwork = lapack_lwork(query[0], n > 1 ? n : 1);
    double* work = (double*)malloc((size_t)lwork * 2 * sizeof(double));
    if (!work) return -999;
    zgeqrf_(&m, &n, lapack_zptr(A), &lda, lapack_zptr(tau),
            lapack_zptr(work), &lwork, &info);
    free(work);
    return info;
}

int mat_lapack_zgeqp3(int m, int n, double* A, int lda,
                      int* jpvt, double* tau)
{
    for (int j = 0; j < n; j++) jpvt[j] = 0;
    int info = 0, lwork = -1;
    double query[2] = { 0.0, 0.0 };
    /* rwork is real, length 2*n. */
    double* rwork = (double*)malloc((size_t)(2 * n) * sizeof(double));
    if (!rwork) return -999;
    zgeqp3_(&m, &n, lapack_zptr(A), &lda, jpvt, lapack_zptr(tau),
            lapack_zptr(query), &lwork, rwork, &info);
    if (info != 0) { free(rwork); return info; }
    lwork = lapack_lwork(query[0], n + 1);
    double* work = (double*)malloc((size_t)lwork * 2 * sizeof(double));
    if (!work) { free(rwork); return -999; }
    zgeqp3_(&m, &n, lapack_zptr(A), &lda, jpvt, lapack_zptr(tau),
            lapack_zptr(work), &lwork, rwork, &info);
    free(work);
    free(rwork);
    return info;
}

int mat_lapack_zungqr(int m, int n, int k, double* A, int lda,
                      const double* tau)
{
    int info = 0, lwork = -1;
    double query[2] = { 0.0, 0.0 };
    double* tau_nc = (double*)tau;
    zungqr_(&m, &n, &k, lapack_zptr(A), &lda, lapack_zptr(tau_nc),
            lapack_zptr(query), &lwork, &info);
    if (info != 0) return info;
    lwork = lapack_lwork(query[0], n > 1 ? n : 1);
    double* work = (double*)malloc((size_t)lwork * 2 * sizeof(double));
    if (!work) return -999;
    zungqr_(&m, &n, &k, lapack_zptr(A), &lda, lapack_zptr(tau_nc),
            lapack_zptr(work), &lwork, &info);
    free(work);
    return info;
}

int mat_lapack_dgetrf(int m, int n, double* A, int lda, int* ipiv)
{
    int info = 0;
    dgetrf_(&m, &n, A, &lda, ipiv, &info);
    return info;
}

int mat_lapack_zgetrf(int m, int n, double* A, int lda, int* ipiv)
{
    int info = 0;
    zgetrf_(&m, &n, lapack_zptr(A), &lda, ipiv, &info);
    return info;
}

int mat_lapack_dpotrf(char uplo, int n, double* A, int lda)
{
    int info = 0;
    char up[2] = { uplo, 0 };
    dpotrf_(up, &n, A, &lda, &info);
    return info;
}

int mat_lapack_zpotrf(char uplo, int n, double* A, int lda)
{
    int info = 0;
    char up[2] = { uplo, 0 };
    zpotrf_(up, &n, lapack_zptr(A), &lda, &info);
    return info;
}

double mat_lapack_dlange(char norm_kind, int m, int n,
                         const double* A, int lda)
{
    /* dlange's workspace requirement is m doubles for the 'I' norm,
     * zero otherwise.  Allocate unconditionally to keep the call site
     * simple. */
    double* work = NULL;
    if (norm_kind == 'I' || norm_kind == 'i') {
        work = (double*)malloc((size_t)(m > 0 ? m : 1) * sizeof(double));
        if (!work) return -1.0;
    }
    char nk[2] = { norm_kind, 0 };
    /* dlange reads A but the Fortran prototype isn't const-qualified.
     * The cast is the documented LAPACK calling pattern. */
    double r = dlange_(nk, &m, &n, (double*)A, &lda, work);
    if (work) free(work);
    return r;
}

double mat_lapack_zlange(char norm_kind, int m, int n,
                         const double* A, int lda)
{
    double* work = NULL;
    if (norm_kind == 'I' || norm_kind == 'i') {
        work = (double*)malloc((size_t)(m > 0 ? m : 1) * sizeof(double));
        if (!work) return -1.0;
    }
    char nk[2] = { norm_kind, 0 };
    double r = zlange_(nk, &m, &n, lapack_zptr((double*)A), &lda, work);
    if (work) free(work);
    return r;
}

int mat_lapack_dgecon(char norm_kind, int n, const double* A_LU,
                      int lda, double anorm, double* rcond)
{
    int info = 0;
    char nk[2] = { norm_kind, 0 };
    /* dgecon workspace: 4n doubles, n ints. */
    double* work = (double*)malloc((size_t)(4 * n > 1 ? 4 * n : 1)
                                    * sizeof(double));
    int* iwork = (int*)malloc((size_t)(n > 1 ? n : 1) * sizeof(int));
    if (!work || !iwork) {
        if (work) free(work);
        if (iwork) free(iwork);
        return -999;
    }
    dgecon_(nk, &n, (double*)A_LU, &lda, &anorm, rcond,
            work, iwork, &info);
    free(work);
    free(iwork);
    return info;
}

int mat_lapack_zgecon(char norm_kind, int n, const double* A_LU,
                      int lda, double anorm, double* rcond)
{
    int info = 0;
    char nk[2] = { norm_kind, 0 };
    /* zgecon workspace: 2n complex doubles (= 4n doubles), 2n real
     * doubles for rwork. */
    double* work  = (double*)malloc((size_t)(2 * n > 1 ? 2 * n : 1)
                                     * 2 * sizeof(double));
    double* rwork = (double*)malloc((size_t)(2 * n > 1 ? 2 * n : 1)
                                     * sizeof(double));
    if (!work || !rwork) {
        if (work) free(work);
        if (rwork) free(rwork);
        return -999;
    }
    zgecon_(nk, &n, lapack_zptr((double*)A_LU), &lda, &anorm, rcond,
            lapack_zptr(work), rwork, &info);
    free(work);
    free(rwork);
    return info;
}

/* ---------------------------------------------------------------------
 * SVD wrappers.
 *
 * dgesdd / zgesdd use the divide-and-conquer algorithm -- materially
 * faster than dgesvd / zgesvd for medium-large matrices and what
 * Mathematica uses internally.  The iwork buffer holds 8 * min(m, n)
 * Fortran integers (we use C `int`).
 *
 * For zgesdd, LAPACK documents the rwork size:
 *   jobz = 'N' :  max(1, 7 * min(m, n))
 *   else        :  max(1, 5*mn*mn + 7*mn) where mn = min(m, n)
 * We allocate the always-safe upper bound; for our typical usage
 * (jobz = 'A') this is the required size anyway.
 * ------------------------------------------------------------------ */
int mat_lapack_dgesdd(char jobz, int m, int n, double* A, int lda,
                      double* S, double* U, int ldu,
                      double* VT, int ldvt)
{
    if (m <= 0 || n <= 0) return -1;
    int info = 0, lwork = -1;
    int mn = (m < n) ? m : n;
    char jz[2] = { jobz, 0 };
    int* iwork = (int*)malloc((size_t)(8 * mn) * sizeof(int));
    if (!iwork) return -999;
    double query = 0.0;
    dgesdd_(jz, &m, &n, A, &lda, S, U, &ldu, VT, &ldvt,
            &query, &lwork, iwork, &info);
    if (info != 0) { free(iwork); return info; }
    lwork = lapack_lwork(query, 1);
    double* work = (double*)malloc((size_t)lwork * sizeof(double));
    if (!work) { free(iwork); return -999; }
    dgesdd_(jz, &m, &n, A, &lda, S, U, &ldu, VT, &ldvt,
            work, &lwork, iwork, &info);
    free(work);
    free(iwork);
    return info;
}

int mat_lapack_zgesdd(char jobz, int m, int n, double* A, int lda,
                      double* S, double* U, int ldu,
                      double* VT, int ldvt)
{
    if (m <= 0 || n <= 0) return -1;
    int info = 0, lwork = -1;
    int mn = (m < n) ? m : n;
    char jz[2] = { jobz, 0 };
    int* iwork = (int*)malloc((size_t)(8 * mn) * sizeof(int));
    if (!iwork) return -999;
    /* rwork: upper-bound formula for jobz != 'N'. */
    size_t rwork_n = (size_t)(5 * mn * mn + 7 * mn);
    if (rwork_n < 1) rwork_n = 1;
    double* rwork = (double*)malloc(rwork_n * sizeof(double));
    if (!rwork) { free(iwork); return -999; }
    double query[2] = { 0.0, 0.0 };
    zgesdd_(jz, &m, &n, lapack_zptr(A), &lda, S,
            lapack_zptr(U), &ldu, lapack_zptr(VT), &ldvt,
            lapack_zptr(query), &lwork, rwork, iwork, &info);
    if (info != 0) { free(rwork); free(iwork); return info; }
    lwork = lapack_lwork(query[0], 1);
    double* work = (double*)malloc((size_t)lwork * 2 * sizeof(double));
    if (!work) { free(rwork); free(iwork); return -999; }
    zgesdd_(jz, &m, &n, lapack_zptr(A), &lda, S,
            lapack_zptr(U), &ldu, lapack_zptr(VT), &ldvt,
            lapack_zptr(work), &lwork, rwork, iwork, &info);
    free(work);
    free(rwork);
    free(iwork);
    return info;
}

/* ---------------------------------------------------------------------
 * Generalized SVD wrappers.
 *
 * dggsvd3 / zggsvd3 compute the generalized SVD of (A, B) where A is
 * m x n and B is p x n.  The iwork buffer is of length n.  See the
 * LAPACK manual for the precise output structure of U, V, Q, alpha,
 * beta, k, l.
 * ------------------------------------------------------------------ */
int mat_lapack_dggsvd3(char jobu, char jobv, char jobq,
                       int m, int n, int p,
                       int* k_out, int* l_out,
                       double* A, int lda, double* B, int ldb,
                       double* alpha, double* beta,
                       double* U, int ldu,
                       double* V, int ldv,
                       double* Q, int ldq)
{
    if (m <= 0 || n <= 0 || p <= 0) return -1;
    int info = 0, lwork = -1;
    char ju[2] = { jobu, 0 };
    char jv[2] = { jobv, 0 };
    char jq[2] = { jobq, 0 };
    int* iwork = (int*)malloc((size_t)n * sizeof(int));
    if (!iwork) return -999;
    double query = 0.0;
    dggsvd3_(ju, jv, jq, &m, &n, &p, k_out, l_out,
             A, &lda, B, &ldb, alpha, beta,
             U, &ldu, V, &ldv, Q, &ldq,
             &query, &lwork, iwork, &info);
    if (info != 0) { free(iwork); return info; }
    lwork = lapack_lwork(query, 1);
    double* work = (double*)malloc((size_t)lwork * sizeof(double));
    if (!work) { free(iwork); return -999; }
    dggsvd3_(ju, jv, jq, &m, &n, &p, k_out, l_out,
             A, &lda, B, &ldb, alpha, beta,
             U, &ldu, V, &ldv, Q, &ldq,
             work, &lwork, iwork, &info);
    free(work);
    free(iwork);
    return info;
}

int mat_lapack_zggsvd3(char jobu, char jobv, char jobq,
                       int m, int n, int p,
                       int* k_out, int* l_out,
                       double* A, int lda, double* B, int ldb,
                       double* alpha, double* beta,
                       double* U, int ldu,
                       double* V, int ldv,
                       double* Q, int ldq)
{
    if (m <= 0 || n <= 0 || p <= 0) return -1;
    int info = 0, lwork = -1;
    char ju[2] = { jobu, 0 };
    char jv[2] = { jobv, 0 };
    char jq[2] = { jobq, 0 };
    int* iwork = (int*)malloc((size_t)n * sizeof(int));
    if (!iwork) return -999;
    /* rwork length is 2*n per the LAPACK docs. */
    double* rwork = (double*)malloc((size_t)(2 * n) * sizeof(double));
    if (!rwork) { free(iwork); return -999; }
    double query[2] = { 0.0, 0.0 };
    zggsvd3_(ju, jv, jq, &m, &n, &p, k_out, l_out,
             lapack_zptr(A), &lda, lapack_zptr(B), &ldb, alpha, beta,
             lapack_zptr(U), &ldu, lapack_zptr(V), &ldv,
             lapack_zptr(Q), &ldq,
             lapack_zptr(query), &lwork, rwork, iwork, &info);
    if (info != 0) { free(rwork); free(iwork); return info; }
    lwork = lapack_lwork(query[0], 1);
    double* work = (double*)malloc((size_t)lwork * 2 * sizeof(double));
    if (!work) { free(rwork); free(iwork); return -999; }
    zggsvd3_(ju, jv, jq, &m, &n, &p, k_out, l_out,
             lapack_zptr(A), &lda, lapack_zptr(B), &ldb, alpha, beta,
             lapack_zptr(U), &ldu, lapack_zptr(V), &ldv,
             lapack_zptr(Q), &ldq,
             lapack_zptr(work), &lwork, rwork, iwork, &info);
    free(work);
    free(rwork);
    free(iwork);
    return info;
}

/* ---------------------------------------------------------------------
 * Solve / eigen / SVD wrappers used by lapack_bridge.c.  Same
 * workspace-query + allocate + dispatch + free convention as above.
 * ------------------------------------------------------------------ */

int mat_lapack_dgesv(int n, int nrhs, double* A, int lda, int* ipiv,
                     double* B, int ldb)
{
    int info = 0;
    dgesv_(&n, &nrhs, A, &lda, ipiv, B, &ldb, &info);
    return info;
}

int mat_lapack_zgesv(int n, int nrhs, double* A, int lda, int* ipiv,
                     double* B, int ldb)
{
    int info = 0;
    zgesv_(&n, &nrhs, lapack_zptr(A), &lda, ipiv, lapack_zptr(B), &ldb, &info);
    return info;
}

int mat_lapack_dgetri(int n, double* A, int lda, const int* ipiv)
{
    int info = 0, lwork = -1;
    double query = 0.0;
    int* ip = (int*)ipiv;
    dgetri_(&n, A, &lda, ip, &query, &lwork, &info);
    if (info != 0) return info;
    lwork = lapack_lwork(query, n > 1 ? n : 1);
    double* work = (double*)malloc((size_t)lwork * sizeof(double));
    if (!work) return -999;
    dgetri_(&n, A, &lda, ip, work, &lwork, &info);
    free(work);
    return info;
}

int mat_lapack_zgetri(int n, double* A, int lda, const int* ipiv)
{
    int info = 0, lwork = -1;
    double query[2] = { 0.0, 0.0 };
    int* ip = (int*)ipiv;
    zgetri_(&n, lapack_zptr(A), &lda, ip, lapack_zptr(query), &lwork, &info);
    if (info != 0) return info;
    lwork = lapack_lwork(query[0], n > 1 ? n : 1);
    double* work = (double*)malloc((size_t)lwork * 2 * sizeof(double));
    if (!work) return -999;
    zgetri_(&n, lapack_zptr(A), &lda, ip, lapack_zptr(work), &lwork, &info);
    free(work);
    return info;
}

int mat_lapack_dgeev(int n, double* A, int lda,
                     double* wr, double* wi, double* VR, int ldvr)
{
    int info = 0, lwork = -1;
    double query = 0.0, vl = 0.0;
    int ldvl = 1;
    char jvl[2] = { 'N', 0 }, jvr[2] = { 'V', 0 };
    dgeev_(jvl, jvr, &n, A, &lda, wr, wi, &vl, &ldvl, VR, &ldvr,
           &query, &lwork, &info);
    if (info != 0) return info;
    lwork = lapack_lwork(query, 4 * n > 1 ? 4 * n : 1);
    double* work = (double*)malloc((size_t)lwork * sizeof(double));
    if (!work) return -999;
    dgeev_(jvl, jvr, &n, A, &lda, wr, wi, &vl, &ldvl, VR, &ldvr,
           work, &lwork, &info);
    free(work);
    return info;
}


/* ---------------------------------------------------------------------- *
 *  Complex Schur / eigenvector self-check and fallback.                   *
 *                                                                         *
 *  Reference LAPACK 3.12.1 (and the OpenBLAS builds that bundle it) has a *
 *  regression in the complex blocked multishift QR (Reference-LAPACK      *
 *  issue #1160, introduced by PR #1020): for n >= 75 -- the zlaqr0 path   *
 *  -- zgees returns a unitary Q and a triangular T with Q T Q^H != A, and *
 *  zgeev wrong eigenvectors (residual ~0.1 on a random 91 x 91), while    *
 *  the eigenvalues stay right. Found on the Ubuntu CI runner by a         *
 *  standalone zgees program. Every call therefore checks its answer with  *
 *  one random vector, O(n^2): A (Q x) against Q (T x), or A (V x) against *
 *  V (diag(w) x). A failing answer is recomputed by the unblocked chain   *
 *  zgehrd + zunghr + zlahqr (+ ztrevc) -- the path LAPACK itself takes    *
 *  below n = 75 -- on non-Accelerate builds (Accelerate's LAPACK predates *
 *  the regression).                                                       *
 * ---------------------------------------------------------------------- */

/* y = M x for a column-major complex n x n M (interleaved). */
static void zmv(int n, const double* M, const double* x, double* y)
{
    for (int i = 0; i < n; i++) { y[2 * i] = 0.0; y[2 * i + 1] = 0.0; }
    for (int j = 0; j < n; j++) {
        double xr = x[2 * j], xi = x[2 * j + 1];
        const double* c = M + 2 * (size_t)j * (size_t)n;
        for (int i = 0; i < n; i++) {
            y[2 * i]     += c[2 * i] * xr - c[2 * i + 1] * xi;
            y[2 * i + 1] += c[2 * i] * xi + c[2 * i + 1] * xr;
        }
    }
}

/* Does A0 X == X D hold, with D the triangular T (D_diag false) or diag(w)
 * (D_diag true)? Tested on one fixed pseudo-random x, relative to ||A0||. */
static bool zlapack_verify(int n, const double* A0, const double* X, const double* T,
                           const double* w, bool D_diag)
{
    if (n <= 0) return true;
    size_t vn = 2 * (size_t)n;
    double* x  = (double*)malloc(sizeof(double) * vn);
    double* u  = (double*)malloc(sizeof(double) * vn);
    double* y1 = (double*)malloc(sizeof(double) * vn);
    double* y2 = (double*)malloc(sizeof(double) * vn);
    if (!x || !u || !y1 || !y2) { free(x); free(u); free(y1); free(y2); return true; }
    unsigned s = 2463534242u;
    double xn = 0.0;
    for (size_t k = 0; k < vn; k++) {
        s ^= s << 13; s ^= s >> 17; s ^= s << 5;
        x[k] = (double)(s % 20001u) / 10000.0 - 1.0;
        xn += x[k] * x[k];
    }
    xn = sqrt(xn);
    /* y1 = A0 (X x) */
    zmv(n, X, x, u);
    zmv(n, A0, u, y1);
    /* y2 = X (D x) */
    if (D_diag) {
        for (int i = 0; i < n; i++) {
            double wr = w[2 * i], wi = w[2 * i + 1], xr = x[2 * i], xi = x[2 * i + 1];
            u[2 * i] = wr * xr - wi * xi;
            u[2 * i + 1] = wr * xi + wi * xr;
        }
    } else {
        /* T is upper triangular */
        for (int i = 0; i < n; i++) {
            double sr = 0.0, si = 0.0;
            for (int j = i; j < n; j++) {
                const double* t = T + 2 * ((size_t)j * (size_t)n + (size_t)i);
                sr += t[0] * x[2 * j] - t[1] * x[2 * j + 1];
                si += t[0] * x[2 * j + 1] + t[1] * x[2 * j];
            }
            u[2 * i] = sr;
            u[2 * i + 1] = si;
        }
    }
    zmv(n, X, u, y2);
    double an = 0.0, d = 0.0;
    for (size_t k = 0; k < 2 * (size_t)n * (size_t)n; k++) { double a = fabs(A0[k]); if (a > an) an = a; }
    for (size_t k = 0; k < vn; k++) { double e = fabs(y1[k] - y2[k]); if (e > d) d = e; }
    free(x); free(u); free(y1); free(y2);
    /* a correct factorisation is ~ n eps ||A|| ||x||; a broken one O(||A||) */
    return d <= 1e-8 * (an > 0.0 ? an : 1.0) * (xn > 0.0 ? xn : 1.0) * (double)n;
}

#if !defined(MATHILDA_USE_ACCELERATE)
/* Complex Schur form of A0 by the unblocked chain: A <- T, VS <- Q, w the
 * eigenvalues. 0 on success. */
static int zschur_unblocked(int n, const double* A0, double* A, double* w, double* VS)
{
    size_t nn = 2 * (size_t)n * (size_t)n;
    memcpy(A, A0, sizeof(double) * nn);
    double* tau = (double*)malloc(sizeof(double) * 2 * (size_t)(n > 1 ? n : 1));
    if (!tau) return -999;
    int ilo = 1, ihi = n, info;
    info = mat_lapack_zgehrd(n, ilo, ihi, A, n, tau);
    if (info != 0) { free(tau); return info; }
    memcpy(VS, A, sizeof(double) * nn);
    info = mat_lapack_zunghr(n, ilo, ihi, VS, n, tau);
    free(tau);
    if (info != 0) return info;
    for (int j = 0; j < n; j++)
        for (int i = j + 2; i < n; i++) {
            A[2 * ((size_t)j * (size_t)n + (size_t)i)] = 0.0;
            A[2 * ((size_t)j * (size_t)n + (size_t)i) + 1] = 0.0;
        }
    int wantt = 1, wantz = 1;
    zlahqr_(&wantt, &wantz, &n, &ilo, &ihi, A, &n, w, &ilo, &ihi, VS, &n, &info);
    return info;
}
#endif

int mat_lapack_zgeev(int n, double* A, int lda, double* w, double* VR, int ldvr)
{
    int info = 0, lwork = -1;
    double query[2] = { 0.0, 0.0 }, vl[2] = { 0.0, 0.0 };
    int ldvl = 1;
    char jvl[2] = { 'N', 0 }, jvr[2] = { 'V', 0 };
    double* rwork = (double*)malloc((size_t)(2 * n > 1 ? 2 * n : 1) * sizeof(double));
    if (!rwork) return -999;
    zgeev_(jvl, jvr, &n, lapack_zptr(A), &lda, lapack_zptr(w),
           lapack_zptr(vl), &ldvl, lapack_zptr(VR), &ldvr,
           lapack_zptr(query), &lwork, rwork, &info);
    if (info != 0) { free(rwork); return info; }
    lwork = lapack_lwork(query[0], 2 * n > 1 ? 2 * n : 1);
    double* work = (double*)malloc((size_t)lwork * 2 * sizeof(double));
    if (!work) { free(rwork); return -999; }
    double* A0 = (lda == n && ldvr == n) ? (double*)malloc(sizeof(double) * 2 * (size_t)n * (size_t)n) : NULL;
    if (A0) memcpy(A0, A, sizeof(double) * 2 * (size_t)n * (size_t)n);
    zgeev_(jvl, jvr, &n, lapack_zptr(A), &lda, lapack_zptr(w),
           lapack_zptr(vl), &ldvl, lapack_zptr(VR), &ldvr,
           lapack_zptr(work), &lwork, rwork, &info);
    free(work);
    free(rwork);
    if (info == 0 && A0 && !zlapack_verify(n, A0, VR, NULL, w, true)) {
#if !defined(MATHILDA_USE_ACCELERATE)
        /* the unblocked Schur form, then the eigenvectors of T back-
         * transformed by Q (ztrevc, HOWMNY = 'B'), normalised as zgeev does:
         * unit 2-norm, largest component real */
        info = zschur_unblocked(n, A0, A, w, VR);
        if (info == 0) {
            double* wk = (double*)malloc(sizeof(double) * 4 * (size_t)(n > 1 ? n : 1));
            double* rw = (double*)malloc(sizeof(double) * (size_t)(n > 1 ? n : 1));
            int mm = n, m = 0, one = 1;
            if (!wk || !rw) info = -999;
            else ztrevc_("R", "B", NULL, &n, A, &n, NULL, &one, VR, &n, &mm, &m, wk, rw, &info, 1, 1);
            free(wk);
            free(rw);
        }
        if (info == 0)
            for (int j = 0; j < n; j++) {
                double* c = VR + 2 * (size_t)j * (size_t)n;
                double nrm = 0.0, big = -1.0;
                int at = 0;
                for (int i = 0; i < n; i++) {
                    double a = c[2 * i] * c[2 * i] + c[2 * i + 1] * c[2 * i + 1];
                    nrm += a;
                    if (a > big) { big = a; at = i; }
                }
                nrm = sqrt(nrm);
                if (nrm == 0.0) continue;
                double pr = c[2 * at] / sqrt(big), pi = -c[2 * at + 1] / sqrt(big);   /* conj phase */
                for (int i = 0; i < n; i++) {
                    double r = c[2 * i], im = c[2 * i + 1];
                    c[2 * i] = (r * pr - im * pi) / nrm;
                    c[2 * i + 1] = (r * pi + im * pr) / nrm;
                }
                c[2 * at + 1] = 0.0;
            }
        if (info == 0 && !zlapack_verify(n, A0, VR, NULL, w, true)) info = -1001;
#else
        info = -1001;
#endif
    }
    free(A0);
    return info;
}

int mat_lapack_dgeev_values(int n, double* A, int lda, double* wr, double* wi)
{
    int info = 0, lwork = -1;
    double query = 0.0, vl = 0.0, vr = 0.0;
    int ldvl = 1, ldvr = 1;
    char jvl[2] = { 'N', 0 }, jvr[2] = { 'N', 0 };
    dgeev_(jvl, jvr, &n, A, &lda, wr, wi, &vl, &ldvl, &vr, &ldvr,
           &query, &lwork, &info);
    if (info != 0) return info;
    lwork = lapack_lwork(query, 3 * n > 1 ? 3 * n : 1);
    double* work = (double*)malloc((size_t)lwork * sizeof(double));
    if (!work) return -999;
    dgeev_(jvl, jvr, &n, A, &lda, wr, wi, &vl, &ldvl, &vr, &ldvr,
           work, &lwork, &info);
    free(work);
    return info;
}

int mat_lapack_zgeev_values(int n, double* A, int lda, double* w)
{
    int info = 0, lwork = -1;
    double query[2] = { 0.0, 0.0 }, vl[2] = { 0.0, 0.0 }, vr[2] = { 0.0, 0.0 };
    int ldvl = 1, ldvr = 1;
    char jvl[2] = { 'N', 0 }, jvr[2] = { 'N', 0 };
    double* rwork = (double*)malloc((size_t)(2 * n > 1 ? 2 * n : 1) * sizeof(double));
    if (!rwork) return -999;
    zgeev_(jvl, jvr, &n, lapack_zptr(A), &lda, lapack_zptr(w),
           lapack_zptr(vl), &ldvl, lapack_zptr(vr), &ldvr,
           lapack_zptr(query), &lwork, rwork, &info);
    if (info != 0) { free(rwork); return info; }
    lwork = lapack_lwork(query[0], 2 * n > 1 ? 2 * n : 1);
    double* work = (double*)malloc((size_t)lwork * 2 * sizeof(double));
    if (!work) { free(rwork); return -999; }
    zgeev_(jvl, jvr, &n, lapack_zptr(A), &lda, lapack_zptr(w),
           lapack_zptr(vl), &ldvl, lapack_zptr(vr), &ldvr,
           lapack_zptr(work), &lwork, rwork, &info);
    free(work);
    free(rwork);
    return info;
}

int mat_lapack_dsyev(int n, double* A, int lda, double* w)
{
    int info = 0, lwork = -1;
    double query = 0.0;
    char jz[2] = { 'V', 0 }, up[2] = { 'U', 0 };
    dsyev_(jz, up, &n, A, &lda, w, &query, &lwork, &info);
    if (info != 0) return info;
    lwork = lapack_lwork(query, 3 * n - 1 > 1 ? 3 * n - 1 : 1);
    double* work = (double*)malloc((size_t)lwork * sizeof(double));
    if (!work) return -999;
    dsyev_(jz, up, &n, A, &lda, w, work, &lwork, &info);
    free(work);
    return info;
}

int mat_lapack_dsyev_values(int n, double* A, int lda, double* w)
{
    int info = 0, lwork = -1;
    double query = 0.0;
    char jz[2] = { 'N', 0 }, up[2] = { 'U', 0 };
    dsyev_(jz, up, &n, A, &lda, w, &query, &lwork, &info);
    if (info != 0) return info;
    lwork = lapack_lwork(query, 3 * n - 1 > 1 ? 3 * n - 1 : 1);
    double* work = (double*)malloc((size_t)lwork * sizeof(double));
    if (!work) return -999;
    dsyev_(jz, up, &n, A, &lda, w, work, &lwork, &info);
    free(work);
    return info;
}

int mat_lapack_dsbevd(int n, int kd, const double* A, double* w, double* Z, int vectors)
{
    if (n <= 0 || kd < 0 || kd >= n) return -1;
    int ldab = kd + 1, info = 0, lwork = -1, liwork = -1, iq = 0, ldz = vectors ? n : 1;
    double query = 0.0, zdummy = 0.0;
    char jz[2] = { vectors ? 'V' : 'N', 0 }, up[2] = { 'U', 0 };
    double* AB = (double*)malloc((size_t)ldab * (size_t)n * sizeof(double));
    if (!AB) return -999;
    for (int j = 0; j < n; j++) {
        int i0 = j - kd > 0 ? j - kd : 0;
        for (int i = i0; i <= j; i++)
            AB[(size_t)(kd + i - j) + (size_t)j * (size_t)ldab] = A[(size_t)i * (size_t)n + (size_t)j];
    }
    double* zp = vectors ? Z : &zdummy;
    dsbevd_(jz, up, &n, &kd, AB, &ldab, w, zp, &ldz, &query, &lwork, &iq, &liwork, &info MTH_FSTRLEN2);
    if (info != 0) { free(AB); return info; }
    lwork = lapack_lwork(query, vectors ? 1 + 5 * n + 2 * n * n : 2 * n);
    liwork = iq > 0 ? iq : (vectors ? 3 + 5 * n : 1);
    double* work = (double*)malloc((size_t)lwork * sizeof(double));
    int* iwork = (int*)malloc((size_t)liwork * sizeof(int));
    if (!work || !iwork) { free(AB); free(work); free(iwork); return -999; }
    dsbevd_(jz, up, &n, &kd, AB, &ldab, w, zp, &ldz, work, &lwork, iwork, &liwork, &info MTH_FSTRLEN2);
    free(AB); free(work); free(iwork);
    return info;
}

/* Half-bandwidth of a symmetric n x n matrix: the largest |i - j| with a
 * nonzero entry (exact zeros only: band storage must reproduce A exactly).
 * Gives up (returns n) once the band is too wide to be worth using. */
static int sym_half_bandwidth(int n, const double* A, int give_up)
{
    int kd = 0;
    for (int i = 0; i < n; i++) {
        const double* row = A + (size_t)i * (size_t)n;
        for (int j = n - 1; j > i + kd; j--)
            if (row[j] != 0.0) { kd = j - i; break; }
        if (kd > give_up) return n;
    }
    return kd;
}

int mat_lapack_dsyev_auto(int n, double* A, double* w, int vectors)
{
    if (vectors) return mat_lapack_dsyev(n, A, n, w);
    /* Values only of a narrowly banded matrix: dsbevd's band reduction is
     * O(n^2 kd).  Measured against the dense routines below it wins while
     * kd < ~n/40 (n = 1000: 0.06 s at kd = 5 against 0.17 s; kd = 60 already
     * loses, 0.25 s), so the cut is kd <= n/64. */
    if (n >= 128) {
        int kd = sym_half_bandwidth(n, A, n / 64);
        if (kd <= n / 64 && mat_lapack_dsbevd(n, kd, A, w, NULL, 0) == 0) return 0;
    }
#if defined(MATHILDA_USE_ACCELERATE)
    /* Values only, large n: the two-stage reduction (dense -> band in BLAS-3
     * blocks, then band -> tridiagonal) does the O(n^3) part as matrix-matrix
     * products, where dsytrd spends half of it in matrix-vector ones: 1.8x on
     * a 1600 x 1600 symmetric matrix.  Below ~800 the band stage's overhead
     * loses (0.016 s against 0.011 s at n = 400).  Any failure falls back. */
    if (n >= MTH_DSYEV_2STAGE_MIN_N) {
        int info = 0, lwork = -1;
        double query = 0.0;
        char jz[2] = { 'N', 0 }, up[2] = { 'U', 0 };
        double* Ac = (double*)malloc((size_t)n * (size_t)n * sizeof(double));
        if (Ac) {
            memcpy(Ac, A, (size_t)n * (size_t)n * sizeof(double));
            dsyev_2stage_(jz, up, &n, Ac, &n, w, &query, &lwork, &info);
            if (info == 0) {
                lwork = lapack_lwork(query, 3 * n);
                double* work = (double*)malloc((size_t)lwork * sizeof(double));
                if (work) {
                    dsyev_2stage_(jz, up, &n, Ac, &n, w, work, &lwork, &info);
                    free(work);
                    if (info == 0) { free(Ac); return 0; }
                }
            }
            free(Ac);
        }
    }
#endif
    return mat_lapack_dsyev_values(n, A, n, w);
}

int mat_lapack_zheev(int n, double* A, int lda, double* w)
{
    int info = 0, lwork = -1;
    double query[2] = { 0.0, 0.0 };
    char jz[2] = { 'V', 0 }, up[2] = { 'U', 0 };
    size_t rwn = (size_t)(3 * n - 2 > 1 ? 3 * n - 2 : 1);
    double* rwork = (double*)malloc(rwn * sizeof(double));
    if (!rwork) return -999;
    zheev_(jz, up, &n, lapack_zptr(A), &lda, w,
           lapack_zptr(query), &lwork, rwork, &info);
    if (info != 0) { free(rwork); return info; }
    lwork = lapack_lwork(query[0], 2 * n - 1 > 1 ? 2 * n - 1 : 1);
    double* work = (double*)malloc((size_t)lwork * 2 * sizeof(double));
    if (!work) { free(rwork); return -999; }
    zheev_(jz, up, &n, lapack_zptr(A), &lda, w,
           lapack_zptr(work), &lwork, rwork, &info);
    free(work);
    free(rwork);
    return info;
}

int mat_lapack_zheevd(int n, double* A, int lda, double* w, int vectors)
{
    if (n <= 0) return 0;
    int info = 0, lwork = -1, lrwork = -1, liwork = -1;
    double query[2] = { 0.0, 0.0 }, rquery = 0.0;
    int iquery = 0;
    char jz[2] = { vectors ? 'V' : 'N', 0 }, up[2] = { 'U', 0 };
    size_t nn = (size_t)n * (size_t)lda;
    double* A0 = NULL;
    if (vectors) {                       /* the input, for the self-check */
        A0 = (double*)malloc(sizeof(double) * 2 * nn);
        if (A0) memcpy(A0, A, sizeof(double) * 2 * nn);
    }
    zheevd_(jz, up, &n, lapack_zptr(A), &lda, w, lapack_zptr(query), &lwork,
            &rquery, &lrwork, &iquery, &liwork, &info MTH_FSTRLEN2);
    if (info != 0) { free(A0); return info; }
    lwork  = lapack_lwork(query[0], 2 * n + n * n);
    lrwork = lapack_lwork(rquery, 1 + 5 * n + 2 * n * n);
    liwork = iquery > 3 + 5 * n ? iquery : 3 + 5 * n;
    double* work  = (double*)malloc((size_t)lwork * 2 * sizeof(double));
    double* rwork = (double*)malloc((size_t)lrwork * sizeof(double));
    int*    iwork = (int*)malloc((size_t)liwork * sizeof(int));
    if (!work || !rwork || !iwork) {
        free(work); free(rwork); free(iwork); free(A0);
        return -999;
    }
    zheevd_(jz, up, &n, lapack_zptr(A), &lda, w, lapack_zptr(work), &lwork,
            rwork, &lrwork, iwork, &liwork, &info MTH_FSTRLEN2);
    free(work); free(rwork); free(iwork);
    if (info == 0 && A0 && lda == n) {
        /* zheevd reads only the upper triangle; the self-check multiplies by
         * the full matrix, so complete the lower half from it. */
        for (int j = 0; j < n; j++)
            for (int i = j + 1; i < n; i++) {
                size_t lo = 2 * ((size_t)i + (size_t)j * (size_t)n);
                size_t hi = 2 * ((size_t)j + (size_t)i * (size_t)n);
                A0[lo] = A0[hi]; A0[lo + 1] = -A0[hi + 1];
            }
        double* wc = (double*)calloc(2 * (size_t)n, sizeof(double));
        if (wc) {
            for (int i = 0; i < n; i++) wc[2 * i] = w[i];
            if (!zlapack_verify(n, A0, A, NULL, wc, true)) info = -1001;
            free(wc);
        }
    }
    free(A0);
    return info;
}

int mat_lapack_dggev(int n, double* A, int lda, double* B, int ldb,
                     double* alphar, double* alphai, double* beta,
                     double* VR, int ldvr)
{
    int info = 0, lwork = -1;
    double query = 0.0, vl = 0.0;
    int ldvl = 1;
    char jvl[2] = { 'N', 0 }, jvr[2] = { 'V', 0 };
    dggev_(jvl, jvr, &n, A, &lda, B, &ldb, alphar, alphai, beta,
           &vl, &ldvl, VR, &ldvr, &query, &lwork, &info);
    if (info != 0) return info;
    lwork = lapack_lwork(query, 8 * n > 1 ? 8 * n : 1);
    double* work = (double*)malloc((size_t)lwork * sizeof(double));
    if (!work) return -999;
    dggev_(jvl, jvr, &n, A, &lda, B, &ldb, alphar, alphai, beta,
           &vl, &ldvl, VR, &ldvr, work, &lwork, &info);
    free(work);
    return info;
}

int mat_lapack_zggev(int n, double* A, int lda, double* B, int ldb,
                     double* alpha, double* beta, double* VR, int ldvr)
{
    int info = 0, lwork = -1;
    double query[2] = { 0.0, 0.0 }, vl[2] = { 0.0, 0.0 };
    int ldvl = 1;
    char jvl[2] = { 'N', 0 }, jvr[2] = { 'V', 0 };
    double* rwork = (double*)malloc((size_t)(8 * n > 1 ? 8 * n : 1) * sizeof(double));
    if (!rwork) return -999;
    zggev_(jvl, jvr, &n, lapack_zptr(A), &lda, lapack_zptr(B), &ldb,
           lapack_zptr(alpha), lapack_zptr(beta),
           lapack_zptr(vl), &ldvl, lapack_zptr(VR), &ldvr,
           lapack_zptr(query), &lwork, rwork, &info);
    if (info != 0) { free(rwork); return info; }
    lwork = lapack_lwork(query[0], 2 * n > 1 ? 2 * n : 1);
    double* work = (double*)malloc((size_t)lwork * 2 * sizeof(double));
    if (!work) { free(rwork); return -999; }
    zggev_(jvl, jvr, &n, lapack_zptr(A), &lda, lapack_zptr(B), &ldb,
           lapack_zptr(alpha), lapack_zptr(beta),
           lapack_zptr(vl), &ldvl, lapack_zptr(VR), &ldvr,
           lapack_zptr(work), &lwork, rwork, &info);
    free(work);
    free(rwork);
    return info;
}

int mat_lapack_dsygv(int n, double* A, int lda, double* B, int ldb, double* w)
{
    int info = 0, lwork = -1, itype = 1;
    double query = 0.0;
    char jz[2] = { 'V', 0 }, up[2] = { 'U', 0 };
    dsygv_(&itype, jz, up, &n, A, &lda, B, &ldb, w, &query, &lwork, &info);
    if (info != 0) return info;
    lwork = lapack_lwork(query, 3 * n - 1 > 1 ? 3 * n - 1 : 1);
    double* work = (double*)malloc((size_t)lwork * sizeof(double));
    if (!work) return -999;
    dsygv_(&itype, jz, up, &n, A, &lda, B, &ldb, w, work, &lwork, &info);
    free(work);
    return info;
}

int mat_lapack_zhegv(int n, double* A, int lda, double* B, int ldb, double* w)
{
    int info = 0, lwork = -1, itype = 1;
    double query[2] = { 0.0, 0.0 };
    char jz[2] = { 'V', 0 }, up[2] = { 'U', 0 };
    size_t rwn = (size_t)(3 * n - 2 > 1 ? 3 * n - 2 : 1);
    double* rwork = (double*)malloc(rwn * sizeof(double));
    if (!rwork) return -999;
    zhegv_(&itype, jz, up, &n, lapack_zptr(A), &lda, lapack_zptr(B), &ldb, w,
           lapack_zptr(query), &lwork, rwork, &info);
    if (info != 0) { free(rwork); return info; }
    lwork = lapack_lwork(query[0], 2 * n - 1 > 1 ? 2 * n - 1 : 1);
    double* work = (double*)malloc((size_t)lwork * 2 * sizeof(double));
    if (!work) { free(rwork); return -999; }
    zhegv_(&itype, jz, up, &n, lapack_zptr(A), &lda, lapack_zptr(B), &ldb, w,
           lapack_zptr(work), &lwork, rwork, &info);
    free(work);
    free(rwork);
    return info;
}

int mat_lapack_dgels(int m, int n, int nrhs, double* A, int lda,
                     double* B, int ldb)
{
    int info = 0, lwork = -1;
    double query = 0.0;
    char tr[2] = { 'N', 0 };
    int mn = (m < n) ? m : n;
    int mx = (mn > nrhs) ? mn : nrhs;
    int fl = mn + mx;
    dgels_(tr, &m, &n, &nrhs, A, &lda, B, &ldb, &query, &lwork, &info);
    if (info != 0) return info;
    lwork = lapack_lwork(query, fl > 1 ? fl : 1);
    double* work = (double*)malloc((size_t)lwork * sizeof(double));
    if (!work) return -999;
    dgels_(tr, &m, &n, &nrhs, A, &lda, B, &ldb, work, &lwork, &info);
    free(work);
    return info;
}

int mat_lapack_zgels(int m, int n, int nrhs, double* A, int lda,
                     double* B, int ldb)
{
    int info = 0, lwork = -1;
    double query[2] = { 0.0, 0.0 };
    char tr[2] = { 'N', 0 };
    int mn = (m < n) ? m : n;
    int mx = (mn > nrhs) ? mn : nrhs;
    int fl = mn + mx;
    zgels_(tr, &m, &n, &nrhs, lapack_zptr(A), &lda, lapack_zptr(B), &ldb,
           lapack_zptr(query), &lwork, &info);
    if (info != 0) return info;
    lwork = lapack_lwork(query[0], fl > 1 ? fl : 1);
    double* work = (double*)malloc((size_t)lwork * 2 * sizeof(double));
    if (!work) return -999;
    zgels_(tr, &m, &n, &nrhs, lapack_zptr(A), &lda, lapack_zptr(B), &ldb,
           lapack_zptr(work), &lwork, &info);
    free(work);
    return info;
}

int mat_lapack_dgesvd(int m, int n, double* A, int lda,
                      double* S, double* U, int ldu, double* VT, int ldvt)
{
    if (m <= 0 || n <= 0) return -1;
    int info = 0, lwork = -1;
    double query = 0.0;
    char ja[2] = { 'A', 0 }, jvt[2] = { 'A', 0 };
    int mn = (m < n) ? m : n;
    int mx = (m > n) ? m : n;
    int fl = 5 * mn;
    if (3 * mn + mx > fl) fl = 3 * mn + mx;
    dgesvd_(ja, jvt, &m, &n, A, &lda, S, U, &ldu, VT, &ldvt,
            &query, &lwork, &info);
    if (info != 0) return info;
    lwork = lapack_lwork(query, fl > 1 ? fl : 1);
    double* work = (double*)malloc((size_t)lwork * sizeof(double));
    if (!work) return -999;
    dgesvd_(ja, jvt, &m, &n, A, &lda, S, U, &ldu, VT, &ldvt,
            work, &lwork, &info);
    free(work);
    return info;
}

int mat_lapack_zgesvd(int m, int n, double* A, int lda,
                      double* S, double* U, int ldu, double* VT, int ldvt)
{
    if (m <= 0 || n <= 0) return -1;
    int info = 0, lwork = -1;
    double query[2] = { 0.0, 0.0 };
    char ja[2] = { 'A', 0 }, jvt[2] = { 'A', 0 };
    int mn = (m < n) ? m : n;
    int mx = (m > n) ? m : n;
    double* rwork = (double*)malloc((size_t)(5 * mn > 1 ? 5 * mn : 1) * sizeof(double));
    if (!rwork) return -999;
    zgesvd_(ja, jvt, &m, &n, lapack_zptr(A), &lda, S,
            lapack_zptr(U), &ldu, lapack_zptr(VT), &ldvt,
            lapack_zptr(query), &lwork, rwork, &info);
    if (info != 0) { free(rwork); return info; }
    int fl = 2 * mn + mx;
    lwork = lapack_lwork(query[0], fl > 1 ? fl : 1);
    double* work = (double*)malloc((size_t)lwork * 2 * sizeof(double));
    if (!work) { free(rwork); return -999; }
    zgesvd_(ja, jvt, &m, &n, lapack_zptr(A), &lda, S,
            lapack_zptr(U), &ldu, lapack_zptr(VT), &ldvt,
            lapack_zptr(work), &lwork, rwork, &info);
    free(work);
    free(rwork);
    return info;
}

int mat_lapack_dtrtrs(int n, int nrhs, const double* A, int lda,
                      double* B, int ldb)
{
    int info = 0;
    char up[2] = { 'U', 0 }, tr[2] = { 'N', 0 }, dg[2] = { 'N', 0 };
    dtrtrs_(up, tr, dg, &n, &nrhs, (double*)A, &lda, B, &ldb, &info);
    return info;
}

int mat_lapack_ztrtrs(int n, int nrhs, const double* A, int lda,
                      double* B, int ldb)
{
    int info = 0;
    char up[2] = { 'U', 0 }, tr[2] = { 'N', 0 }, dg[2] = { 'N', 0 };
    ztrtrs_(up, tr, dg, &n, &nrhs, lapack_zptr((double*)A), &lda,
            lapack_zptr(B), &ldb, &info);
    return info;
}

/* The literal 0 passed for the SELECT / SELCTG function-pointer arguments is a
 * null pointer constant (C11 6.3.2.3p3) that converts to any function-pointer
 * type; the NULL macro must NOT be used there because it may expand to
 * (void*)0, which does not implicitly convert to a function pointer and would
 * trip -Werror=incompatible-pointer-types.  SORT='N' means LAPACK never
 * references SELECT or BWORK, so a null for both is safe. */
/* ---- Hessenberg reduction (HessenbergDecomposition) ---------------------- */

int mat_lapack_dgehrd(int n, int ilo, int ihi, double* A, int lda, double* tau)
{
    int info = 0, lwork = -1;
    double query = 0.0;
    dgehrd_(&n, &ilo, &ihi, A, &lda, tau, &query, &lwork, &info);
    if (info != 0) return info;
    lwork = lapack_lwork(query, n > 1 ? n : 1);
    double* work = (double*)malloc((size_t)lwork * sizeof(double));
    if (!work) return -999;
    dgehrd_(&n, &ilo, &ihi, A, &lda, tau, work, &lwork, &info);
    free(work);
    return info;
}

int mat_lapack_zgehrd(int n, int ilo, int ihi, double* A, int lda, double* tau)
{
    int info = 0, lwork = -1;
    double query[2] = { 0.0, 0.0 };
    zgehrd_(&n, &ilo, &ihi, lapack_zptr(A), &lda, lapack_zptr(tau),
            lapack_zptr(query), &lwork, &info);
    if (info != 0) return info;
    lwork = lapack_lwork(query[0], n > 1 ? n : 1);
    double* work = (double*)malloc((size_t)lwork * 2 * sizeof(double));
    if (!work) return -999;
    zgehrd_(&n, &ilo, &ihi, lapack_zptr(A), &lda, lapack_zptr(tau),
            lapack_zptr(work), &lwork, &info);
    free(work);
    return info;
}

int mat_lapack_dorghr(int n, int ilo, int ihi, double* A, int lda, const double* tau)
{
    int info = 0, lwork = -1;
    double query = 0.0;
    dorghr_(&n, &ilo, &ihi, A, &lda, (double*)tau, &query, &lwork, &info);
    if (info != 0) return info;
    lwork = lapack_lwork(query, n > 1 ? n : 1);
    double* work = (double*)malloc((size_t)lwork * sizeof(double));
    if (!work) return -999;
    dorghr_(&n, &ilo, &ihi, A, &lda, (double*)tau, work, &lwork, &info);
    free(work);
    return info;
}

int mat_lapack_zunghr(int n, int ilo, int ihi, double* A, int lda, const double* tau)
{
    int info = 0, lwork = -1;
    double query[2] = { 0.0, 0.0 };
    zunghr_(&n, &ilo, &ihi, lapack_zptr(A), &lda, lapack_zptr((double*)tau),
            lapack_zptr(query), &lwork, &info);
    if (info != 0) return info;
    lwork = lapack_lwork(query[0], n > 1 ? n : 1);
    double* work = (double*)malloc((size_t)lwork * 2 * sizeof(double));
    if (!work) return -999;
    zunghr_(&n, &ilo, &ihi, lapack_zptr(A), &lda, lapack_zptr((double*)tau),
            lapack_zptr(work), &lwork, &info);
    free(work);
    return info;
}

int mat_lapack_dgees(int n, double* A, int lda,
                     double* wr, double* wi, double* VS, int ldvs)
{
    int info = 0, lwork = -1, sdim = 0;
    double query = 0.0;
    char jvs[2] = { 'V', 0 }, srt[2] = { 'N', 0 };
    dgees_(jvs, srt, 0, &n, A, &lda, &sdim, wr, wi, VS, &ldvs,
           &query, &lwork, NULL, &info MTH_FSTRLEN2);
    if (info != 0) return info;
    lwork = lapack_lwork(query, 3 * n > 1 ? 3 * n : 1);
    double* work = (double*)malloc((size_t)lwork * sizeof(double));
    if (!work) return -999;
    dgees_(jvs, srt, 0, &n, A, &lda, &sdim, wr, wi, VS, &ldvs,
           work, &lwork, NULL, &info MTH_FSTRLEN2);
    free(work);
    return info;
}

int mat_lapack_zgees(int n, double* A, int lda, double* w, double* VS, int ldvs)
{
    int info = 0, lwork = -1, sdim = 0;
    double query[2] = { 0.0, 0.0 };
    char jvs[2] = { 'V', 0 }, srt[2] = { 'N', 0 };
    double* rwork = (double*)malloc((size_t)(n > 1 ? n : 1) * sizeof(double));
    if (!rwork) return -999;
    zgees_(jvs, srt, 0, &n, lapack_zptr(A), &lda, &sdim, lapack_zptr(w),
           lapack_zptr(VS), &ldvs, lapack_zptr(query), &lwork, rwork, NULL, &info MTH_FSTRLEN2);
    if (info != 0) { free(rwork); return info; }
    lwork = lapack_lwork(query[0], 2 * n > 1 ? 2 * n : 1);
    double* work = (double*)malloc((size_t)lwork * 2 * sizeof(double));
    if (!work) { free(rwork); return -999; }
    double* A0 = (lda == n && ldvs == n) ? (double*)malloc(sizeof(double) * 2 * (size_t)n * (size_t)n) : NULL;
    if (A0) memcpy(A0, A, sizeof(double) * 2 * (size_t)n * (size_t)n);
    zgees_(jvs, srt, 0, &n, lapack_zptr(A), &lda, &sdim, lapack_zptr(w),
           lapack_zptr(VS), &ldvs, lapack_zptr(work), &lwork, rwork, NULL, &info MTH_FSTRLEN2);
    free(work);
    free(rwork);
    if (info == 0 && A0 && !zlapack_verify(n, A0, VS, A, w, false)) {
#if !defined(MATHILDA_USE_ACCELERATE)
        info = zschur_unblocked(n, A0, A, w, VS);
        if (info == 0 && !zlapack_verify(n, A0, VS, A, w, false)) info = -1001;
#else
        info = -1001;
#endif
    }
    free(A0);
    return info;
}

int mat_lapack_dgges(int n, double* A, int lda, double* B, int ldb,
                     double* alphar, double* alphai, double* beta,
                     double* VSL, int ldvsl, double* VSR, int ldvsr)
{
    int info = 0, lwork = -1, sdim = 0;
    double query = 0.0;
    char jl[2] = { 'V', 0 }, jr[2] = { 'V', 0 }, srt[2] = { 'N', 0 };
    dgges_(jl, jr, srt, 0, &n, A, &lda, B, &ldb, &sdim,
           alphar, alphai, beta, VSL, &ldvsl, VSR, &ldvsr,
           &query, &lwork, NULL, &info);
    if (info != 0) return info;
    int fl = 8 * n + 16;
    lwork = lapack_lwork(query, fl > 1 ? fl : 1);
    double* work = (double*)malloc((size_t)lwork * sizeof(double));
    if (!work) return -999;
    dgges_(jl, jr, srt, 0, &n, A, &lda, B, &ldb, &sdim,
           alphar, alphai, beta, VSL, &ldvsl, VSR, &ldvsr,
           work, &lwork, NULL, &info);
    free(work);
    return info;
}

int mat_lapack_zgges(int n, double* A, int lda, double* B, int ldb,
                     double* alpha, double* beta,
                     double* VSL, int ldvsl, double* VSR, int ldvsr)
{
    int info = 0, lwork = -1, sdim = 0;
    double query[2] = { 0.0, 0.0 };
    char jl[2] = { 'V', 0 }, jr[2] = { 'V', 0 }, srt[2] = { 'N', 0 };
    double* rwork = (double*)malloc((size_t)(8 * n > 1 ? 8 * n : 1) * sizeof(double));
    if (!rwork) return -999;
    zgges_(jl, jr, srt, 0, &n, lapack_zptr(A), &lda, lapack_zptr(B), &ldb, &sdim,
           lapack_zptr(alpha), lapack_zptr(beta),
           lapack_zptr(VSL), &ldvsl, lapack_zptr(VSR), &ldvsr,
           lapack_zptr(query), &lwork, rwork, NULL, &info);
    if (info != 0) { free(rwork); return info; }
    lwork = lapack_lwork(query[0], 2 * n > 1 ? 2 * n : 1);
    double* work = (double*)malloc((size_t)lwork * 2 * sizeof(double));
    if (!work) { free(rwork); return -999; }
    zgges_(jl, jr, srt, 0, &n, lapack_zptr(A), &lda, lapack_zptr(B), &ldb, &sdim,
           lapack_zptr(alpha), lapack_zptr(beta),
           lapack_zptr(VSL), &ldvsl, lapack_zptr(VSR), &ldvsr,
           lapack_zptr(work), &lwork, rwork, NULL, &info);
    free(work);
    free(rwork);
    return info;
}

int mat_lapack_dgebal(char job, int n, double* A, int lda,
                      int* ilo, int* ihi, double* scale)
{
    int info = 0;
    char jb[2] = { job, 0 };
    dgebal_(jb, &n, A, &lda, ilo, ihi, scale, &info);
    return info;
}

int mat_lapack_zgebal(char job, int n, double* A, int lda,
                      int* ilo, int* ihi, double* scale)
{
    int info = 0;
    char jb[2] = { job, 0 };
    zgebal_(jb, &n, lapack_zptr(A), &lda, ilo, ihi, scale, &info);
    return info;
}

int mat_lapack_dgebak(char job, char side, int n, int ilo, int ihi,
                      const double* scale, int m, double* V, int ldv)
{
    int info = 0;
    char jb[2] = { job, 0 }, sd[2] = { side, 0 };
    dgebak_(jb, sd, &n, &ilo, &ihi, scale, &m, V, &ldv, &info);
    return info;
}

int mat_lapack_zgebak(char job, char side, int n, int ilo, int ihi,
                      const double* scale, int m, double* V, int ldv)
{
    int info = 0;
    char jb[2] = { job, 0 }, sd[2] = { side, 0 };
    zgebak_(jb, sd, &n, &ilo, &ihi, scale, &m, lapack_zptr(V), &ldv, &info);
    return info;
}


/* ---------------------------------------------------------------------
 * Banded LU / Cholesky (see lapack.h).
 * ------------------------------------------------------------------- */
int mat_lapack_dgbtrf(int n, int kl, int ku, double* ab, int ldab, int* ipiv)
{
    int info = 0;
    dgbtrf_(&n, &n, &kl, &ku, ab, &ldab, ipiv, &info);
    return info;
}
int mat_lapack_zgbtrf(int n, int kl, int ku, double* ab, int ldab, int* ipiv)
{
    int info = 0;
    zgbtrf_(&n, &n, &kl, &ku, lapack_zptr(ab), &ldab, ipiv, &info);
    return info;
}
int mat_lapack_dgbtrs(char trans, int n, int kl, int ku, int nrhs, const double* ab,
                      int ldab, const int* ipiv, double* b, int ldb)
{
    int info = 0;
    char t[2] = { trans, 0 };
    dgbtrs_(t, &n, &kl, &ku, &nrhs, (double*)ab, &ldab, (int*)ipiv, b, &ldb, &info);
    return info;
}
int mat_lapack_zgbtrs(char trans, int n, int kl, int ku, int nrhs, const double* ab,
                      int ldab, const int* ipiv, double* b, int ldb)
{
    int info = 0;
    char t[2] = { trans, 0 };
    zgbtrs_(t, &n, &kl, &ku, &nrhs, lapack_zptr((double*)ab), &ldab, (int*)ipiv,
            lapack_zptr(b), &ldb, &info);
    return info;
}
int mat_lapack_dpbtrf(int n, int kd, double* ab, int ldab)
{
    int info = 0;
    dpbtrf_("L", &n, &kd, ab, &ldab, &info);
    return info;
}
int mat_lapack_zpbtrf(int n, int kd, double* ab, int ldab)
{
    int info = 0;
    zpbtrf_("L", &n, &kd, lapack_zptr(ab), &ldab, &info);
    return info;
}
int mat_lapack_dpbtrs(int n, int kd, int nrhs, const double* ab, int ldab, double* b, int ldb)
{
    int info = 0;
    dpbtrs_("L", &n, &kd, &nrhs, (double*)ab, &ldab, b, &ldb, &info);
    return info;
}
int mat_lapack_zpbtrs(int n, int kd, int nrhs, const double* ab, int ldab, double* b, int ldb)
{
    int info = 0;
    zpbtrs_("L", &n, &kd, &nrhs, lapack_zptr((double*)ab), &ldab, lapack_zptr(b), &ldb, &info);
    return info;
}

#else  /* !USE_LAPACK */
int mat_lapack_dgbtrf(int n, int kl, int ku, double* ab, int ldab, int* ipiv)
{ (void)n; (void)kl; (void)ku; (void)ab; (void)ldab; (void)ipiv; return -1; }
int mat_lapack_zgbtrf(int n, int kl, int ku, double* ab, int ldab, int* ipiv)
{ (void)n; (void)kl; (void)ku; (void)ab; (void)ldab; (void)ipiv; return -1; }
int mat_lapack_dgbtrs(char trans, int n, int kl, int ku, int nrhs, const double* ab,
                      int ldab, const int* ipiv, double* b, int ldb)
{ (void)trans; (void)n; (void)kl; (void)ku; (void)nrhs; (void)ab; (void)ldab; (void)ipiv; (void)b; (void)ldb; return -1; }
int mat_lapack_zgbtrs(char trans, int n, int kl, int ku, int nrhs, const double* ab,
                      int ldab, const int* ipiv, double* b, int ldb)
{ (void)trans; (void)n; (void)kl; (void)ku; (void)nrhs; (void)ab; (void)ldab; (void)ipiv; (void)b; (void)ldb; return -1; }
int mat_lapack_dpbtrf(int n, int kd, double* ab, int ldab)
{ (void)n; (void)kd; (void)ab; (void)ldab; return -1; }
int mat_lapack_zpbtrf(int n, int kd, double* ab, int ldab)
{ (void)n; (void)kd; (void)ab; (void)ldab; return -1; }
int mat_lapack_dpbtrs(int n, int kd, int nrhs, const double* ab, int ldab, double* b, int ldb)
{ (void)n; (void)kd; (void)nrhs; (void)ab; (void)ldab; (void)b; (void)ldb; return -1; }
int mat_lapack_zpbtrs(int n, int kd, int nrhs, const double* ab, int ldab, double* b, int ldb)
{ (void)n; (void)kd; (void)nrhs; (void)ab; (void)ldab; (void)b; (void)ldb; return -1; }


/* Stubs so call sites link cleanly when LAPACK isn't available.  Each
 * returns a negative info code; the kernel interprets any negative
 * info as "LAPACK unavailable, fall back to the symbolic pipeline". */
int mat_lapack_dgeqrf(int m, int n, double* A, int lda, double* tau)
{ (void)m; (void)n; (void)A; (void)lda; (void)tau; return -1; }
int mat_lapack_dgeqp3(int m, int n, double* A, int lda,
                      int* jpvt, double* tau)
{ (void)m; (void)n; (void)A; (void)lda; (void)jpvt; (void)tau; return -1; }
int mat_lapack_dorgqr(int m, int n, int k, double* A, int lda,
                      const double* tau)
{ (void)m; (void)n; (void)k; (void)A; (void)lda; (void)tau; return -1; }
int mat_lapack_zgeqrf(int m, int n, double* A, int lda, double* tau)
{ (void)m; (void)n; (void)A; (void)lda; (void)tau; return -1; }
int mat_lapack_zgeqp3(int m, int n, double* A, int lda,
                      int* jpvt, double* tau)
{ (void)m; (void)n; (void)A; (void)lda; (void)jpvt; (void)tau; return -1; }
int mat_lapack_zungqr(int m, int n, int k, double* A, int lda,
                      const double* tau)
{ (void)m; (void)n; (void)k; (void)A; (void)lda; (void)tau; return -1; }

int mat_lapack_dgetrf(int m, int n, double* A, int lda, int* ipiv)
{ (void)m; (void)n; (void)A; (void)lda; (void)ipiv; return -1; }
int mat_lapack_zgetrf(int m, int n, double* A, int lda, int* ipiv)
{ (void)m; (void)n; (void)A; (void)lda; (void)ipiv; return -1; }
int mat_lapack_dpotrf(char uplo, int n, double* A, int lda)
{ (void)uplo; (void)n; (void)A; (void)lda; return -1; }
int mat_lapack_zpotrf(char uplo, int n, double* A, int lda)
{ (void)uplo; (void)n; (void)A; (void)lda; return -1; }
double mat_lapack_dlange(char norm_kind, int m, int n,
                         const double* A, int lda)
{ (void)norm_kind; (void)m; (void)n; (void)A; (void)lda; return -1.0; }
double mat_lapack_zlange(char norm_kind, int m, int n,
                         const double* A, int lda)
{ (void)norm_kind; (void)m; (void)n; (void)A; (void)lda; return -1.0; }
int mat_lapack_dgecon(char norm_kind, int n, const double* A_LU,
                      int lda, double anorm, double* rcond)
{ (void)norm_kind; (void)n; (void)A_LU; (void)lda; (void)anorm; (void)rcond; return -1; }
int mat_lapack_zgecon(char norm_kind, int n, const double* A_LU,
                      int lda, double anorm, double* rcond)
{ (void)norm_kind; (void)n; (void)A_LU; (void)lda; (void)anorm; (void)rcond; return -1; }

int mat_lapack_dgesdd(char jobz, int m, int n, double* A, int lda,
                      double* S, double* U, int ldu,
                      double* VT, int ldvt)
{ (void)jobz; (void)m; (void)n; (void)A; (void)lda; (void)S;
  (void)U; (void)ldu; (void)VT; (void)ldvt; return -1; }
int mat_lapack_zgesdd(char jobz, int m, int n, double* A, int lda,
                      double* S, double* U, int ldu,
                      double* VT, int ldvt)
{ (void)jobz; (void)m; (void)n; (void)A; (void)lda; (void)S;
  (void)U; (void)ldu; (void)VT; (void)ldvt; return -1; }
int mat_lapack_dggsvd3(char jobu, char jobv, char jobq,
                       int m, int n, int p,
                       int* k_out, int* l_out,
                       double* A, int lda, double* B, int ldb,
                       double* alpha, double* beta,
                       double* U, int ldu,
                       double* V, int ldv,
                       double* Q, int ldq)
{ (void)jobu; (void)jobv; (void)jobq; (void)m; (void)n; (void)p;
  (void)k_out; (void)l_out; (void)A; (void)lda; (void)B; (void)ldb;
  (void)alpha; (void)beta; (void)U; (void)ldu; (void)V; (void)ldv;
  (void)Q; (void)ldq; return -1; }
int mat_lapack_zggsvd3(char jobu, char jobv, char jobq,
                       int m, int n, int p,
                       int* k_out, int* l_out,
                       double* A, int lda, double* B, int ldb,
                       double* alpha, double* beta,
                       double* U, int ldu,
                       double* V, int ldv,
                       double* Q, int ldq)
{ (void)jobu; (void)jobv; (void)jobq; (void)m; (void)n; (void)p;
  (void)k_out; (void)l_out; (void)A; (void)lda; (void)B; (void)ldb;
  (void)alpha; (void)beta; (void)U; (void)ldu; (void)V; (void)ldv;
  (void)Q; (void)ldq; return -1; }

int mat_lapack_dgesv(int n, int nrhs, double* A, int lda, int* ipiv,
                     double* B, int ldb)
{ (void)n; (void)nrhs; (void)A; (void)lda; (void)ipiv; (void)B; (void)ldb; return -1; }
int mat_lapack_zgesv(int n, int nrhs, double* A, int lda, int* ipiv,
                     double* B, int ldb)
{ (void)n; (void)nrhs; (void)A; (void)lda; (void)ipiv; (void)B; (void)ldb; return -1; }
int mat_lapack_dgetri(int n, double* A, int lda, const int* ipiv)
{ (void)n; (void)A; (void)lda; (void)ipiv; return -1; }
int mat_lapack_zgetri(int n, double* A, int lda, const int* ipiv)
{ (void)n; (void)A; (void)lda; (void)ipiv; return -1; }
int mat_lapack_dgeev(int n, double* A, int lda,
                     double* wr, double* wi, double* VR, int ldvr)
{ (void)n; (void)A; (void)lda; (void)wr; (void)wi; (void)VR; (void)ldvr; return -1; }
int mat_lapack_zgeev(int n, double* A, int lda, double* w, double* VR, int ldvr)
{ (void)n; (void)A; (void)lda; (void)w; (void)VR; (void)ldvr; return -1; }
int mat_lapack_dgeev_values(int n, double* A, int lda, double* wr, double* wi)
{ (void)n; (void)A; (void)lda; (void)wr; (void)wi; return -1; }
int mat_lapack_zgeev_values(int n, double* A, int lda, double* w)
{ (void)n; (void)A; (void)lda; (void)w; return -1; }
int mat_lapack_dsyev(int n, double* A, int lda, double* w)
{ (void)n; (void)A; (void)lda; (void)w; return -1; }
int mat_lapack_dsyev_values(int n, double* A, int lda, double* w)
{ (void)n; (void)A; (void)lda; (void)w; return -1; }
int mat_lapack_zheev(int n, double* A, int lda, double* w)
{ (void)n; (void)A; (void)lda; (void)w; return -1; }
int mat_lapack_zheevd(int n, double* A, int lda, double* w, int vectors)
{ (void)n; (void)A; (void)lda; (void)w; (void)vectors; return -1; }
int mat_lapack_dggev(int n, double* A, int lda, double* B, int ldb,
                     double* alphar, double* alphai, double* beta,
                     double* VR, int ldvr)
{ (void)n; (void)A; (void)lda; (void)B; (void)ldb; (void)alphar; (void)alphai;
  (void)beta; (void)VR; (void)ldvr; return -1; }
int mat_lapack_zggev(int n, double* A, int lda, double* B, int ldb,
                     double* alpha, double* beta, double* VR, int ldvr)
{ (void)n; (void)A; (void)lda; (void)B; (void)ldb; (void)alpha; (void)beta;
  (void)VR; (void)ldvr; return -1; }
int mat_lapack_dsygv(int n, double* A, int lda, double* B, int ldb, double* w)
{ (void)n; (void)A; (void)lda; (void)B; (void)ldb; (void)w; return -1; }
int mat_lapack_zhegv(int n, double* A, int lda, double* B, int ldb, double* w)
{ (void)n; (void)A; (void)lda; (void)B; (void)ldb; (void)w; return -1; }
int mat_lapack_dgels(int m, int n, int nrhs, double* A, int lda,
                     double* B, int ldb)
{ (void)m; (void)n; (void)nrhs; (void)A; (void)lda; (void)B; (void)ldb; return -1; }
int mat_lapack_zgels(int m, int n, int nrhs, double* A, int lda,
                     double* B, int ldb)
{ (void)m; (void)n; (void)nrhs; (void)A; (void)lda; (void)B; (void)ldb; return -1; }
int mat_lapack_dgesvd(int m, int n, double* A, int lda,
                      double* S, double* U, int ldu, double* VT, int ldvt)
{ (void)m; (void)n; (void)A; (void)lda; (void)S; (void)U; (void)ldu;
  (void)VT; (void)ldvt; return -1; }
int mat_lapack_zgesvd(int m, int n, double* A, int lda,
                      double* S, double* U, int ldu, double* VT, int ldvt)
{ (void)m; (void)n; (void)A; (void)lda; (void)S; (void)U; (void)ldu;
  (void)VT; (void)ldvt; return -1; }
int mat_lapack_dtrtrs(int n, int nrhs, const double* A, int lda,
                      double* B, int ldb)
{ (void)n; (void)nrhs; (void)A; (void)lda; (void)B; (void)ldb; return -1; }
int mat_lapack_ztrtrs(int n, int nrhs, const double* A, int lda,
                      double* B, int ldb)
{ (void)n; (void)nrhs; (void)A; (void)lda; (void)B; (void)ldb; return -1; }
int mat_lapack_dgehrd(int n, int ilo, int ihi, double* A, int lda, double* tau)
{ (void)n; (void)ilo; (void)ihi; (void)A; (void)lda; (void)tau; return -1; }
int mat_lapack_zgehrd(int n, int ilo, int ihi, double* A, int lda, double* tau)
{ (void)n; (void)ilo; (void)ihi; (void)A; (void)lda; (void)tau; return -1; }
int mat_lapack_dorghr(int n, int ilo, int ihi, double* A, int lda, const double* tau)
{ (void)n; (void)ilo; (void)ihi; (void)A; (void)lda; (void)tau; return -1; }
int mat_lapack_zunghr(int n, int ilo, int ihi, double* A, int lda, const double* tau)
{ (void)n; (void)ilo; (void)ihi; (void)A; (void)lda; (void)tau; return -1; }

int mat_lapack_dgees(int n, double* A, int lda,
                     double* wr, double* wi, double* VS, int ldvs)
{ (void)n; (void)A; (void)lda; (void)wr; (void)wi; (void)VS; (void)ldvs; return -1; }
int mat_lapack_zgees(int n, double* A, int lda, double* w, double* VS, int ldvs)
{ (void)n; (void)A; (void)lda; (void)w; (void)VS; (void)ldvs; return -1; }
int mat_lapack_dgges(int n, double* A, int lda, double* B, int ldb,
                     double* alphar, double* alphai, double* beta,
                     double* VSL, int ldvsl, double* VSR, int ldvsr)
{ (void)n; (void)A; (void)lda; (void)B; (void)ldb; (void)alphar; (void)alphai;
  (void)beta; (void)VSL; (void)ldvsl; (void)VSR; (void)ldvsr; return -1; }
int mat_lapack_zgges(int n, double* A, int lda, double* B, int ldb,
                     double* alpha, double* beta,
                     double* VSL, int ldvsl, double* VSR, int ldvsr)
{ (void)n; (void)A; (void)lda; (void)B; (void)ldb; (void)alpha; (void)beta;
  (void)VSL; (void)ldvsl; (void)VSR; (void)ldvsr; return -1; }
int mat_lapack_dgebal(char job, int n, double* A, int lda,
                      int* ilo, int* ihi, double* scale)
{ (void)job; (void)n; (void)A; (void)lda; (void)ilo; (void)ihi; (void)scale; return -1; }
int mat_lapack_zgebal(char job, int n, double* A, int lda,
                      int* ilo, int* ihi, double* scale)
{ (void)job; (void)n; (void)A; (void)lda; (void)ilo; (void)ihi; (void)scale; return -1; }
int mat_lapack_dgebak(char job, char side, int n, int ilo, int ihi,
                      const double* scale, int m, double* V, int ldv)
{ (void)job; (void)side; (void)n; (void)ilo; (void)ihi; (void)scale;
  (void)m; (void)V; (void)ldv; return -1; }
int mat_lapack_zgebak(char job, char side, int n, int ilo, int ihi,
                      const double* scale, int m, double* V, int ldv)
{ (void)job; (void)side; (void)n; (void)ilo; (void)ihi; (void)scale;
  (void)m; (void)V; (void)ldv; return -1; }

#endif /* USE_LAPACK */
