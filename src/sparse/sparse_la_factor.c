/* Mathilda -- sparse linear algebra: direct factorisations (SpFactor).
 *
 * Three backends behind one interface (sparse_la.h):
 *
 *   UMFPACK  unsymmetric multifrontal LU (umfpack_dl_* / umfpack_zl_*, 64-bit
 *            indices, packed complex). Present when the build links SuiteSparse
 *            (USE_UMFPACK).
 *   CHOLMOD  supernodal Cholesky for Hermitian positive-definite matrices
 *            (cholmod_l_*). Present with USE_CHOLMOD.
 *   LAPACK   the matrix densified (column-major) and factored with getrf;
 *            used when SuiteSparse is absent, or asked for, and the dense
 *            form fits the densify cap.
 *
 * The CSR arrays of an SpMat, read as a compressed-COLUMN matrix, describe
 * M = A^T. Rather than transposing, the library factors M and every solve
 * picks the transposed system: A x = b is M^T x = b (UMFPACK_Aat, the plain
 * transpose, also for complex), and A^T x = b is M x = b (UMFPACK_A).
 * det(M) = det(A). For CHOLMOD a Hermitian A gives M = conj(A), so
 * A x = b is solved as M conj(x) = conj(b).
 */
#include "sparse_la.h"
#include "lapack.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef USE_UMFPACK
#include <umfpack.h>
#endif
#ifdef USE_CHOLMOD
#include <cholmod.h>
#endif

/* Densify cap of the LAPACK backend (elements), the same as the gate's. */
#define SPF_DENSE_MAX ((int64_t)1 << 27)

typedef enum { BK_UMFPACK, BK_CHOLMOD, BK_LAPACK } Backend;

struct SpFactor {
    Backend bk;
    int64_t n;
    bool cplx;
    double rcond;
    /* UMFPACK: the arrays it was factored from (needed again by solve) */
    int64_t* ap;
    int64_t* ai;
    double* ax;
    void* numeric;
    /* CHOLMOD */
#ifdef USE_CHOLMOD
    cholmod_common cc;
    cholmod_factor* L;
#endif
    /* LAPACK: column-major LU and its pivots */
    double* lu;
    int* ipiv;
};

int64_t spf_n(const SpFactor* f) { return f->n; }
bool spf_cplx(const SpFactor* f) { return f->cplx; }
double spf_rcond(const SpFactor* f) { return f->rcond; }

bool spf_badly_conditioned(const SpFactor* f) {
    if (f->rcond < 0.0) return false;
    return f->bk == BK_LAPACK ? f->rcond < 1e6 * 2.220446049250313e-16 : f->rcond < 1e-16;
}

const char* spf_backend(const SpFactor* f) {
    switch (f->bk) {
        case BK_UMFPACK: return "UMFPACK";
        case BK_CHOLMOD: return "CHOLMOD";
        default:         return "LAPACK";
    }
}

void spf_free(SpFactor* f) {
    if (!f) return;
#ifdef USE_UMFPACK
    if (f->bk == BK_UMFPACK && f->numeric) {
        if (f->cplx) umfpack_zl_free_numeric(&f->numeric);
        else umfpack_dl_free_numeric(&f->numeric);
    }
#endif
#ifdef USE_CHOLMOD
    if (f->bk == BK_CHOLMOD) {
        if (f->L) cholmod_l_free_factor(&f->L, &f->cc);
        cholmod_l_finish(&f->cc);
    }
#endif
    free(f->ap);
    free(f->ai);
    free(f->ax);
    free(f->lu);
    free(f->ipiv);
    free(f);
}

#if defined(USE_UMFPACK) || defined(USE_CHOLMOD)
/* Owned copies of the CSR arrays (UMFPACK keeps pointers into them). */
static bool copy_arrays(const SpMat* A, SpFactor* f) {
    int w = A->cplx ? 2 : 1;
    size_t nz = (size_t)(A->nnz > 0 ? A->nnz : 1);
    f->ap = malloc(sizeof(int64_t) * (size_t)(A->n + 1));
    f->ai = malloc(sizeof(int64_t) * nz);
    f->ax = malloc(sizeof(double) * (size_t)w * nz);
    if (!f->ap || !f->ai || !f->ax) return false;
    memcpy(f->ap, A->rp, sizeof(int64_t) * (size_t)(A->n + 1));
    memcpy(f->ai, A->ci, sizeof(int64_t) * (size_t)A->nnz);
    memcpy(f->ax, A->x, sizeof(double) * (size_t)w * (size_t)A->nnz);
    return true;
}
#endif

/* ---- UMFPACK -------------------------------------------------------------- */

#ifdef USE_UMFPACK
static SpfStatus factor_umfpack(const SpMat* A, SpFactor* f) {
    if (!copy_arrays(A, f)) return SPF_NOMEM;
    double control[UMFPACK_CONTROL], info[UMFPACK_INFO];
    void* symbolic = NULL;
    int64_t n = A->n;
    int st;
    if (f->cplx) {
        umfpack_zl_defaults(control);
        st = umfpack_zl_symbolic(n, n, f->ap, f->ai, f->ax, NULL, &symbolic, control, info);
        if (st == UMFPACK_OK)
            st = umfpack_zl_numeric(f->ap, f->ai, f->ax, NULL, symbolic, &f->numeric, control, info);
        if (symbolic) umfpack_zl_free_symbolic(&symbolic);
    } else {
        umfpack_dl_defaults(control);
        st = umfpack_dl_symbolic(n, n, f->ap, f->ai, f->ax, &symbolic, control, info);
        if (st == UMFPACK_OK)
            st = umfpack_dl_numeric(f->ap, f->ai, f->ax, symbolic, &f->numeric, control, info);
        if (symbolic) umfpack_dl_free_symbolic(&symbolic);
    }
    if (f->numeric) f->bk = BK_UMFPACK;      /* so spf_free releases it on every exit */
    if (st == UMFPACK_WARNING_singular_matrix) return SPF_SINGULAR;
    if (st == UMFPACK_ERROR_out_of_memory) return SPF_NOMEM;
    if (st != UMFPACK_OK) return SPF_FAILED;
    f->rcond = info[UMFPACK_RCOND];
    /* a pivot so small that the reciprocal condition underflows is
     * singular for every practical purpose */
    if (!(f->rcond > 0.0)) return SPF_SINGULAR;
    f->bk = BK_UMFPACK;
    return SPF_OK;
}

static bool solve_umfpack(const SpFactor* f, const double* B, double* X, int64_t nrhs, bool trans) {
    /* M = A^T as stored: A x = b is the transposed M system */
    int sys = trans ? UMFPACK_A : UMFPACK_Aat;
    double control[UMFPACK_CONTROL], info[UMFPACK_INFO];
    int w = f->cplx ? 2 : 1;
    size_t len = (size_t)w * (size_t)f->n;
    double* tmp = NULL;
    if (X == B) {
        tmp = malloc(sizeof(double) * len);
        if (!tmp) return false;
    }
    if (f->cplx) umfpack_zl_defaults(control);
    else umfpack_dl_defaults(control);
    for (int64_t r = 0; r < nrhs; r++) {
        const double* b = B + (size_t)r * len;
        double* x = tmp ? tmp : X + (size_t)r * len;
        int st = f->cplx
            ? umfpack_zl_solve(sys, f->ap, f->ai, f->ax, NULL, x, NULL, b, NULL, f->numeric, control, info)
            : umfpack_dl_solve(sys, f->ap, f->ai, f->ax, x, b, f->numeric, control, info);
        if (st != UMFPACK_OK && st != UMFPACK_WARNING_singular_matrix) { free(tmp); return false; }
        if (tmp) memcpy(X + (size_t)r * len, tmp, sizeof(double) * len);
    }
    free(tmp);
    return true;
}

static bool det_umfpack(const SpFactor* f, double* re, double* im) {
    double mx[2] = { 0.0, 0.0 }, ex = 0.0, info[UMFPACK_INFO];
    int st = f->cplx ? umfpack_zl_get_determinant(mx, NULL, &ex, f->numeric, info)
                     : umfpack_dl_get_determinant(mx, &ex, f->numeric, info);
    if (st != UMFPACK_OK && st != UMFPACK_WARNING_singular_matrix) return false;
    double s = pow(10.0, ex);
    *re = mx[0] * s;
    *im = f->cplx ? mx[1] * s : 0.0;
    return isfinite(*re) && isfinite(*im);
}
#endif

/* ---- CHOLMOD -------------------------------------------------------------- */

#ifdef USE_CHOLMOD
static void wrap_sparse(const SpFactor* f, cholmod_sparse* S) {
    memset(S, 0, sizeof *S);
    S->nrow = (size_t)f->n;
    S->ncol = (size_t)f->n;
    S->nzmax = (size_t)f->ap[f->n];
    S->p = f->ap;
    S->i = f->ai;
    S->x = f->ax;
    S->stype = -1;                 /* lower triangle of M; M is Hermitian */
    S->itype = CHOLMOD_LONG;
    S->xtype = f->cplx ? CHOLMOD_COMPLEX : CHOLMOD_REAL;
    S->dtype = CHOLMOD_DOUBLE;
    S->sorted = 1;
    S->packed = 1;
}

static SpfStatus factor_cholmod(const SpMat* A, SpFactor* f) {
    if (!copy_arrays(A, f)) return SPF_NOMEM;
    cholmod_l_start(&f->cc);
    f->bk = BK_CHOLMOD;                 /* so spf_free finishes the common */
    f->cc.print = 0;
    f->cc.error_handler = NULL;
    /* LL' always: the simplicial default is LDL', which succeeds on an
     * indefinite matrix, and positive definiteness is what this backend
     * promises (LinearSolve::npdef otherwise) */
    f->cc.supernodal = CHOLMOD_SUPERNODAL;
    f->cc.final_ll = 1;
    cholmod_sparse S;
    wrap_sparse(f, &S);
    f->L = cholmod_l_analyze(&S, &f->cc);
    if (!f->L) return f->cc.status == CHOLMOD_OUT_OF_MEMORY ? SPF_NOMEM : SPF_FAILED;
    cholmod_l_factorize(&S, f->L, &f->cc);
    if (f->cc.status == CHOLMOD_NOT_POSDEF || f->L->minor < f->L->n) return SPF_NOTPD;
    if (f->cc.status == CHOLMOD_OUT_OF_MEMORY) return SPF_NOMEM;
    if (f->cc.status < CHOLMOD_OK) return SPF_FAILED;
    f->rcond = cholmod_l_rcond(f->L, &f->cc);
    return SPF_OK;
}

static bool solve_cholmod(const SpFactor* f, const double* B, double* X, int64_t nrhs, bool trans) {
    int w = f->cplx ? 2 : 1;
    size_t len = (size_t)w * (size_t)f->n;
    /* conjugation: A x = b with M = conj(A) is M conj(x) = conj(b); the
     * transposed system A^T x = b is M x = b directly (real: no-ops) */
    bool conj = f->cplx && !trans;
    double* in = malloc(sizeof(double) * len * (size_t)nrhs);
    if (!in) return false;
    memcpy(in, B, sizeof(double) * len * (size_t)nrhs);
    if (conj)
        for (size_t k = 1; k < len * (size_t)nrhs; k += 2) in[k] = -in[k];
    cholmod_dense D;
    memset(&D, 0, sizeof D);
    D.nrow = (size_t)f->n;
    D.ncol = (size_t)nrhs;
    D.nzmax = (size_t)f->n * (size_t)nrhs;
    D.d = (size_t)f->n;
    D.x = in;
    D.xtype = f->cplx ? CHOLMOD_COMPLEX : CHOLMOD_REAL;
    D.dtype = CHOLMOD_DOUBLE;
    cholmod_common* cc = (cholmod_common*)&f->cc;
    cholmod_dense* Y = cholmod_l_solve(CHOLMOD_A, f->L, &D, cc);
    free(in);
    if (!Y) return false;
    memcpy(X, Y->x, sizeof(double) * len * (size_t)nrhs);
    cholmod_l_free_dense(&Y, cc);
    if (conj)
        for (size_t k = 1; k < len * (size_t)nrhs; k += 2) X[k] = -X[k];
    return true;
}
#endif

/* ---- LAPACK (dense) ------------------------------------------------------- */

static SpfStatus factor_lapack(const SpMat* A, SpFactor* f) {
    int64_t n = A->n;
    if (n > 0 && n > SPF_DENSE_MAX / n) return SPF_TOOBIG;
    int w = A->cplx ? 2 : 1;
    f->lu = calloc((size_t)w * (size_t)(n * n > 0 ? n * n : 1), sizeof(double));
    f->ipiv = malloc(sizeof(int) * (size_t)(n > 0 ? n : 1));
    if (!f->lu || !f->ipiv) return SPF_NOMEM;
    for (int64_t i = 0; i < n; i++)
        for (int64_t k = A->rp[i]; k < A->rp[i + 1]; k++) {
            size_t at = (size_t)w * ((size_t)A->ci[k] * (size_t)n + (size_t)i);  /* column-major */
            f->lu[at] = A->x[w * k];
            if (A->cplx) f->lu[at + 1] = A->x[w * k + 1];
        }
    double anorm = A->cplx ? mat_lapack_zlange('1', (int)n, (int)n, f->lu, (int)n)
                           : mat_lapack_dlange('1', (int)n, (int)n, f->lu, (int)n);
    int info = A->cplx ? mat_lapack_zgetrf((int)n, (int)n, f->lu, (int)n, f->ipiv)
                       : mat_lapack_dgetrf((int)n, (int)n, f->lu, (int)n, f->ipiv);
    if (info < 0) return SPF_FAILED;                  /* no LAPACK in this build */
    if (info > 0) return SPF_SINGULAR;
    f->bk = BK_LAPACK;
    double rc = -1.0;
    if (anorm >= 0.0) {
        int ci = A->cplx ? mat_lapack_zgecon('1', (int)n, f->lu, (int)n, anorm, &rc)
                         : mat_lapack_dgecon('1', (int)n, f->lu, (int)n, anorm, &rc);
        if (ci != 0) rc = -1.0;
    }
    f->rcond = rc;
    if (rc == 0.0) return SPF_SINGULAR;
    return SPF_OK;
}

/* One right-hand side through the LAPACK LU (P A = L U, column-major):
 * A x = b:   b <- P b, L y = b, U x = y
 * A^T x = b: U^T y = b, L^T z = y, x <- P^T z  */
static void lu_solve1(const SpFactor* f, double* x, bool trans) {
    int64_t n = f->n;
    const double* a = f->lu;
#define AR(i, j) a[(size_t)(j) * (size_t)n + (size_t)(i)]
#define ACR(i, j) a[2 * ((size_t)(j) * (size_t)n + (size_t)(i))]
#define ACI(i, j) a[2 * ((size_t)(j) * (size_t)n + (size_t)(i)) + 1]
    if (!f->cplx) {
        if (!trans) {
            for (int64_t i = 0; i < n; i++) {
                int64_t p = f->ipiv[i] - 1;
                if (p != i) { double t = x[i]; x[i] = x[p]; x[p] = t; }
            }
            for (int64_t i = 0; i < n; i++) {
                double s = x[i];
                for (int64_t j = 0; j < i; j++) s -= AR(i, j) * x[j];
                x[i] = s;
            }
            for (int64_t i = n; i-- > 0;) {
                double s = x[i];
                for (int64_t j = i + 1; j < n; j++) s -= AR(i, j) * x[j];
                x[i] = s / AR(i, i);
            }
        } else {
            for (int64_t i = 0; i < n; i++) {
                double s = x[i];
                for (int64_t j = 0; j < i; j++) s -= AR(j, i) * x[j];
                x[i] = s / AR(i, i);
            }
            for (int64_t i = n; i-- > 0;) {
                double s = x[i];
                for (int64_t j = i + 1; j < n; j++) s -= AR(j, i) * x[j];
                x[i] = s;
            }
            for (int64_t i = n; i-- > 0;) {
                int64_t p = f->ipiv[i] - 1;
                if (p != i) { double t = x[i]; x[i] = x[p]; x[p] = t; }
            }
        }
        return;
    }
    /* complex, interleaved */
    if (!trans) {
        for (int64_t i = 0; i < n; i++) {
            int64_t p = f->ipiv[i] - 1;
            if (p != i) {
                double tr = x[2 * i], ti = x[2 * i + 1];
                x[2 * i] = x[2 * p]; x[2 * i + 1] = x[2 * p + 1];
                x[2 * p] = tr; x[2 * p + 1] = ti;
            }
        }
        for (int64_t i = 0; i < n; i++) {
            double sr = x[2 * i], si = x[2 * i + 1];
            for (int64_t j = 0; j < i; j++) {
                double lr = ACR(i, j), li = ACI(i, j);
                sr -= lr * x[2 * j] - li * x[2 * j + 1];
                si -= lr * x[2 * j + 1] + li * x[2 * j];
            }
            x[2 * i] = sr; x[2 * i + 1] = si;
        }
        for (int64_t i = n; i-- > 0;) {
            double sr = x[2 * i], si = x[2 * i + 1];
            for (int64_t j = i + 1; j < n; j++) {
                double ur = ACR(i, j), ui = ACI(i, j);
                sr -= ur * x[2 * j] - ui * x[2 * j + 1];
                si -= ur * x[2 * j + 1] + ui * x[2 * j];
            }
            double dr = ACR(i, i), di = ACI(i, i), d = dr * dr + di * di;
            x[2 * i] = (sr * dr + si * di) / d;
            x[2 * i + 1] = (si * dr - sr * di) / d;
        }
    } else {
        for (int64_t i = 0; i < n; i++) {
            double sr = x[2 * i], si = x[2 * i + 1];
            for (int64_t j = 0; j < i; j++) {
                double ur = ACR(j, i), ui = ACI(j, i);
                sr -= ur * x[2 * j] - ui * x[2 * j + 1];
                si -= ur * x[2 * j + 1] + ui * x[2 * j];
            }
            double dr = ACR(i, i), di = ACI(i, i), d = dr * dr + di * di;
            x[2 * i] = (sr * dr + si * di) / d;
            x[2 * i + 1] = (si * dr - sr * di) / d;
        }
        for (int64_t i = n; i-- > 0;) {
            double sr = x[2 * i], si = x[2 * i + 1];
            for (int64_t j = i + 1; j < n; j++) {
                double lr = ACR(j, i), li = ACI(j, i);
                sr -= lr * x[2 * j] - li * x[2 * j + 1];
                si -= lr * x[2 * j + 1] + li * x[2 * j];
            }
            x[2 * i] = sr; x[2 * i + 1] = si;
        }
        for (int64_t i = n; i-- > 0;) {
            int64_t p = f->ipiv[i] - 1;
            if (p != i) {
                double tr = x[2 * i], ti = x[2 * i + 1];
                x[2 * i] = x[2 * p]; x[2 * i + 1] = x[2 * p + 1];
                x[2 * p] = tr; x[2 * p + 1] = ti;
            }
        }
    }
#undef AR
#undef ACR
#undef ACI
}

static bool det_lapack(const SpFactor* f, double* re, double* im) {
    double dr = 1.0, di = 0.0;
    int64_t n = f->n;
    for (int64_t i = 0; i < n; i++) {
        if (f->ipiv[i] - 1 != i) { dr = -dr; di = -di; }
        if (!f->cplx) {
            dr *= f->lu[(size_t)i * (size_t)n + (size_t)i];
        } else {
            size_t at = 2 * ((size_t)i * (size_t)n + (size_t)i);
            double ar = f->lu[at], ai = f->lu[at + 1];
            double t = dr * ar - di * ai;
            di = dr * ai + di * ar;
            dr = t;
        }
    }
    *re = dr;
    *im = di;
    return isfinite(dr) && isfinite(di);
}

/* ---- the interface -------------------------------------------------------- */

SpFactor* spf_factor(const SpMat* A, SpfWant want, SpfStatus* st) {
    *st = SPF_OK;
    if (A->m != A->n) { *st = SPF_NOTSQUARE; return NULL; }
    SpFactor* f = calloc(1, sizeof *f);
    if (!f) { *st = SPF_NOMEM; return NULL; }
    f->n = A->n;
    f->cplx = A->cplx;
    f->rcond = -1.0;
    f->bk = BK_LAPACK;
    SpfStatus s = SPF_FAILED;

    if (want == SPF_WANT_CHOL || (want == SPF_WANT_AUTO && spla_hermitian_posdiag(A))) {
#ifdef USE_CHOLMOD
        s = factor_cholmod(A, f);
        if (s == SPF_OK) return f;
        if (want == SPF_WANT_CHOL) { *st = s; spf_free(f); return NULL; }
        /* not positive definite after all: LU below, from a fresh struct */
        spf_free(f);
        f = calloc(1, sizeof *f);
        if (!f) { *st = SPF_NOMEM; return NULL; }
        f->n = A->n; f->cplx = A->cplx; f->rcond = -1.0; f->bk = BK_LAPACK;
#else
        if (want == SPF_WANT_CHOL) {
            /* Cholesky without CHOLMOD: positive definiteness via the dense
             * potrf test, then the LU (a Hermitian PD matrix needs no
             * pivoting, so the LU is as stable as the Cholesky) */
            int64_t n = A->n;
            if (n > 0 && n > SPF_DENSE_MAX / n) { *st = SPF_TOOBIG; spf_free(f); return NULL; }
            int w = A->cplx ? 2 : 1;
            double* d = calloc((size_t)w * (size_t)(n * n > 0 ? n * n : 1), sizeof(double));
            if (!d) { *st = SPF_NOMEM; spf_free(f); return NULL; }
            for (int64_t i = 0; i < n; i++)
                for (int64_t k = A->rp[i]; k < A->rp[i + 1]; k++) {
                    size_t at = (size_t)w * ((size_t)A->ci[k] * (size_t)n + (size_t)i);
                    d[at] = A->x[w * k];
                    if (A->cplx) d[at + 1] = A->x[w * k + 1];
                }
            int info = A->cplx ? mat_lapack_zpotrf('L', (int)n, d, (int)n)
                               : mat_lapack_dpotrf('L', (int)n, d, (int)n);
            free(d);
            if (info > 0 || !spla_hermitian_posdiag(A)) { *st = SPF_NOTPD; spf_free(f); return NULL; }
            if (info < 0) { *st = SPF_FAILED; spf_free(f); return NULL; }
        }
#endif
    }
#ifdef USE_UMFPACK
    if (want != SPF_WANT_DENSE) {
        s = factor_umfpack(A, f);
        if (s == SPF_OK) return f;
        *st = s;
        spf_free(f);
        return NULL;
    }
#endif
    s = factor_lapack(A, f);
    if (s == SPF_OK) return f;
    *st = s;
    spf_free(f);
    return NULL;
}

bool spf_solve(const SpFactor* f, const double* B, double* X, int64_t nrhs, bool trans) {
    switch (f->bk) {
#ifdef USE_UMFPACK
        case BK_UMFPACK: return solve_umfpack(f, B, X, nrhs, trans);
#endif
#ifdef USE_CHOLMOD
        case BK_CHOLMOD: return solve_cholmod(f, B, X, nrhs, trans);
#endif
        case BK_LAPACK: {
            int w = f->cplx ? 2 : 1;
            size_t len = (size_t)w * (size_t)f->n;
            if (X != B) memcpy(X, B, sizeof(double) * len * (size_t)nrhs);
            for (int64_t r = 0; r < nrhs; r++) lu_solve1(f, X + (size_t)r * len, trans);
            return true;
        }
        default:
            return false;
    }
}

bool spf_det(const SpFactor* f, double* re, double* im) {
    switch (f->bk) {
#ifdef USE_UMFPACK
        case BK_UMFPACK: return det_umfpack(f, re, im);
#endif
        case BK_LAPACK: return det_lapack(f, re, im);
        default: return false;
    }
}
