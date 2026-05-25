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

#include <stdlib.h>

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

#else  /* !USE_LAPACK */

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

#endif /* USE_LAPACK */
