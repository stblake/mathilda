/* Mathilda -- sparse linear algebra: machine views, SpMV, structural tests.
 *
 * See sparse_la.h for the layer map. Everything in this file is O(nnz) and
 * allocation-light: a view borrows the SparseArray's own CSR buffers and only
 * copies what has the wrong width or type (int32 column indices, int64
 * values).
 */
#include "sparse_la.h"
#include "sparse_internal.h"
#include "ndarray_internal.h"
#include "numarray.h"
#include <math.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif
#include <stdlib.h>
#include <string.h>

bool spla_is_zero(const Expr* e) {
    if (sp_is_zero(e)) return true;
#ifdef USE_MPFR
    if (e && e->type == EXPR_MPFR) return mpfr_zero_p(e->data.mpfr) != 0;
#endif
    return false;
}

bool spla_view(const SparseArray* sa, bool promote_int, SpMat* A) {
    memset(A, 0, sizeof *A);
    if (!sa || sa->rank != 2) return false;
    if (!sp_is_zero(sa->bg)) return false;
    if ((sa->kind == SPV_EXPR || sa->kind == SPV_INT64) && !promote_int) return false;
    const SparsePattern* p = sa->pat;
    A->m = sa->dims[0];
    A->n = sa->dims[1];
    A->nnz = p->nnz;
    A->rp = p->rowptr;
    A->cplx = sa->kind == SPV_C128;
    if (p->idx32) {
        A->owned_ci = malloc(sizeof(int64_t) * (size_t)(p->nnz > 0 ? p->nnz : 1));
        if (!A->owned_ci) return false;
        const int32_t* c32 = p->colidx;
        for (int64_t k = 0; k < p->nnz; k++) A->owned_ci[k] = c32[k];
        A->ci = A->owned_ci;
    } else {
        A->ci = p->colidx;
    }
    if (sa->kind == SPV_EXPR) {
        /* boxed values (a mix of machine reals and complexes, which the
         * array keeps boxed as the reference does, or exact numbers in a
         * machine context): every one must read as a number */
        int64_t nz = p->nnz;
        double* re = malloc(sizeof(double) * 2 * (size_t)(nz > 0 ? nz : 1));
        if (!re) { free(A->owned_ci); memset(A, 0, sizeof *A); return false; }
        Expr* const* ev = sa->vals;
        bool cplx = false;
        for (int64_t k = 0; k < nz; k++) {
            double r, i;
            if (!na_read_scalar(ev[k], &r, &i) || !isfinite(r) || !isfinite(i)) {
                free(re); free(A->owned_ci); memset(A, 0, sizeof *A); return false;
            }
            re[2 * k] = r;
            re[2 * k + 1] = i;
            if (i != 0.0) cplx = true;
        }
        if (!cplx)
            for (int64_t k = 0; k < nz; k++) re[k] = re[2 * k];
        A->owned_x = re;
        A->x = re;
        A->cplx = cplx;
        return true;
    }
    if (sa->kind == SPV_INT64) {
        A->owned_x = malloc(sizeof(double) * (size_t)(p->nnz > 0 ? p->nnz : 1));
        if (!A->owned_x) { free(A->owned_ci); memset(A, 0, sizeof *A); return false; }
        const int64_t* iv = sa->vals;
        for (int64_t k = 0; k < p->nnz; k++) A->owned_x[k] = (double)iv[k];
        A->x = A->owned_x;
    } else {
        A->x = sa->vals;
    }
    return true;
}

void spla_view_free(SpMat* A) {
    if (!A) return;
    free(A->owned_ci);
    free(A->owned_x);
    /* a view built by spla_from_dense owns rp too; it marks that by setting
     * owned_ci (always) and keeping rp inside the same block -- see there */
    memset(A, 0, sizeof *A);
}

/* spla_from_dense keeps rp, ci and x in ONE allocation, recorded as owned_ci
 * (rp and x point into it), so spla_view_free releases everything. */
bool spla_from_dense(const Expr* e, SpMat* A) {
    memset(A, 0, sizeof *A);
    int rows = 0, cols = 0;
    double* buf = NULL;
    bool cplx = false;
    if (!na_load_matrix(e, false, false, &rows, &cols, &buf)) {
        if (!na_load_matrix(e, true, false, &rows, &cols, &buf)) return false;
        cplx = true;
    }
    int w = cplx ? 2 : 1;
    int64_t nnz = 0;
    for (int64_t k = 0; k < (int64_t)rows * cols; k++)
        if (buf[w * k] != 0.0 || (cplx && buf[w * k + 1] != 0.0)) nnz++;
    size_t bytes = sizeof(int64_t) * (size_t)(rows + 1) + sizeof(int64_t) * (size_t)(nnz ? nnz : 1)
                 + sizeof(double) * (size_t)w * (size_t)(nnz ? nnz : 1);
    int64_t* blk = malloc(bytes);
    if (!blk) { free(buf); return false; }
    int64_t* rp = blk;
    int64_t* ci = rp + rows + 1;
    double* x = (double*)(ci + (nnz ? nnz : 1));
    int64_t k = 0;
    rp[0] = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            const double* v = buf + (size_t)w * ((size_t)i * cols + j);
            if (v[0] != 0.0 || (cplx && v[1] != 0.0)) {
                ci[k] = j;
                x[w * k] = v[0];
                if (cplx) x[w * k + 1] = v[1];
                k++;
            }
        }
        rp[i + 1] = k;
    }
    free(buf);
    A->m = rows; A->n = cols; A->nnz = nnz;
    A->rp = rp; A->ci = ci; A->x = x; A->cplx = cplx;
    A->owned_ci = blk;
    return true;
}

/* ---- SpMV ---------------------------------------------------------------- */

typedef struct { const SpMat* A; const double* x; double* y; } MvCtx;

static bool mv_rows(void* c, size_t lo, size_t hi) {
    const MvCtx* t = c;
    const SpMat* A = t->A;
    if (!A->cplx) {
        for (size_t i = lo; i < hi; i++) {
            double acc = 0.0;
            for (int64_t k = A->rp[i]; k < A->rp[i + 1]; k++) acc += A->x[k] * t->x[A->ci[k]];
            t->y[i] = acc;
        }
    } else {
        for (size_t i = lo; i < hi; i++) {
            double re = 0.0, im = 0.0;
            for (int64_t k = A->rp[i]; k < A->rp[i + 1]; k++) {
                double ar = A->x[2 * k], ai = A->x[2 * k + 1];
                double xr = t->x[2 * A->ci[k]], xi = t->x[2 * A->ci[k] + 1];
                re += ar * xr - ai * xi;
                im += ar * xi + ai * xr;
            }
            t->y[2 * i] = re;
            t->y[2 * i + 1] = im;
        }
    }
    return true;
}

void spla_spmv(const SpMat* A, const double* x, double* y, bool trans) {
    if (!trans) {
        MvCtx c = { A, x, y };
        if (A->nnz > 200000) nd_parallel_for((size_t)A->m, mv_rows, &c);
        else mv_rows(&c, 0, (size_t)A->m);
        return;
    }
    /* y = A^T x: a scatter over the rows (serial; the transpose kernels that
     * need speed transpose once and use the forward product) */
    int w = A->cplx ? 2 : 1;
    memset(y, 0, sizeof(double) * (size_t)w * (size_t)A->n);
    for (int64_t i = 0; i < A->m; i++) {
        for (int64_t k = A->rp[i]; k < A->rp[i + 1]; k++) {
            int64_t j = A->ci[k];
            if (!A->cplx) {
                y[j] += A->x[k] * x[i];
            } else {
                double ar = A->x[2 * k], ai = A->x[2 * k + 1];
                double xr = x[2 * i], xi = x[2 * i + 1];
                y[2 * j] += ar * xr - ai * xi;
                y[2 * j + 1] += ar * xi + ai * xr;
            }
        }
    }
}

/* Error-free accumulation (Ogita, Rump & Oishi, "Accurate sum and dot
 * product", SIAM J. Sci. Comput. 26, 2005: Dot2): s + c carries the running
 * sum to about twice the working precision. */
static void acc_add(double* s, double* c, double v) {
    double t = *s + v;
    double z = t - *s;
    *c += (*s - (t - z)) + (v - z);
    *s = t;
}
static void acc_mul(double* s, double* c, double a, double b) {
    double p = a * b;
    acc_add(s, c, p);
    *c += fma(a, b, -p);
}

void spla_residual(const SpMat* A, const double* b, const double* x, double* r) {
    for (int64_t i = 0; i < A->m; i++) {
        if (!A->cplx) {
            double s = b[i], c = 0.0;
            for (int64_t k = A->rp[i]; k < A->rp[i + 1]; k++) acc_mul(&s, &c, -A->x[k], x[A->ci[k]]);
            r[i] = s + c;
        } else {
            double sr = b[2 * i], cr = 0.0, si = b[2 * i + 1], ci = 0.0;
            for (int64_t k = A->rp[i]; k < A->rp[i + 1]; k++) {
                double ar = A->x[2 * k], ai = A->x[2 * k + 1];
                const double* xv = x + 2 * A->ci[k];
                acc_mul(&sr, &cr, -ar, xv[0]);
                acc_mul(&sr, &cr, ai, xv[1]);
                acc_mul(&si, &ci, -ar, xv[1]);
                acc_mul(&si, &ci, -ai, xv[0]);
            }
            r[2 * i] = sr + cr;
            r[2 * i + 1] = si + ci;
        }
    }
}

bool spla_is_square(const SpMat* A) { return A->m == A->n; }

/* Binary search for column j in row i; index into ci/x or -1. */
static int64_t find_entry(const SpMat* A, int64_t i, int64_t j) {
    int64_t lo = A->rp[i], hi = A->rp[i + 1];
    while (lo < hi) {
        int64_t mid = lo + (hi - lo) / 2;
        if (A->ci[mid] < j) lo = mid + 1;
        else hi = mid;
    }
    return (lo < A->rp[i + 1] && A->ci[lo] == j) ? lo : -1;
}

bool spla_hermitian_posdiag(const SpMat* A) {
    if (A->m != A->n) return false;
    for (int64_t i = 0; i < A->m; i++) {
        bool diag = false;
        for (int64_t k = A->rp[i]; k < A->rp[i + 1]; k++) {
            int64_t j = A->ci[k];
            if (j == i) {
                if (!(A->x[A->cplx ? 2 * k : k] > 0.0)) return false;
                if (A->cplx && A->x[2 * k + 1] != 0.0) return false;
                diag = true;
                continue;
            }
            int64_t t = find_entry(A, j, i);
            if (t < 0) return false;
            if (!A->cplx) {
                if (A->x[t] != A->x[k]) return false;
            } else {
                if (A->x[2 * t] != A->x[2 * k] || A->x[2 * t + 1] != -A->x[2 * k + 1]) return false;
            }
        }
        if (!diag) return false;
    }
    return true;
}

bool spla_have_umfpack(void) {
#ifdef USE_UMFPACK
    return true;
#else
    return false;
#endif
}
bool spla_have_cholmod(void) {
#ifdef USE_CHOLMOD
    return true;
#else
    return false;
#endif
}
bool spla_have_arpack(void) {
#ifdef USE_ARPACK
    return true;
#else
    return false;
#endif
}
