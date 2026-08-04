/* ndcorrcov.c -- machine-precision buffer kernels for Covariance / Correlation.
 *
 * These are the single entry points nd_covariance / nd_correlation used by BOTH
 * the interpreter (src/stats/corrcov.c routes here when any argument is an
 * NDArray) and the Compile[] VM (the ND_FNS / ND_FN2S delegated-kernel tables in
 * src/compile/compile.c reference them), so the compiled answer is identical to
 * the REPL by construction — the faithful-degrade contract documented in
 * ndreduce.h.
 *
 * Shapes handled directly on the flat buffer:
 *   Covariance[v, w]  (two rank-1 real vectors)  -> a Real scalar
 *   Covariance[a]     (one rank-2 real matrix)   -> a p x p matrix (auto)
 *   Covariance[a, b]  (two rank-2 real matrices) -> a p x q matrix (cross)
 * and the same three for Correlation.  The vector reduction is a two-pass
 * threaded centered inner product; the matrix gram is A_c^T B_c via BLAS
 * (cblas_dsyrk for the symmetric auto form, cblas_dgemm for the cross form),
 * with a plain triple loop when USE_LAPACK is off.  p and q (the number of
 * variables) are small in practice; n (the number of observations) is large.
 *
 * Anything outside the real-machine case — an int64/bool buffer (whose exact
 * covariance is a Rational no float slot holds), a complex dtype, a shape
 * mismatch, a zero-variance column for Correlation — returns
 * ndarray_delist_and_reeval(res), so the exact/symbolic generic path answers and
 * the result is never wrong, only (in that case) not accelerated.
 *
 * `res` is the whole builtin call; it is BORROWED (never freed). The result is
 * freshly owned by the caller. */

#include <math.h>
#include <stdlib.h>

#include "expr.h"
#include "ndarray.h"
#include "ndarray_internal.h"   /* nd_parallel_reduce, NDARRAY_MAX_THREADS */
#include "numarray.h"           /* na_load_matrix */
#include "ndreduce.h"           /* our own prototypes */

#ifdef USE_LAPACK
#include "lapack.h"             /* CBLAS: cblas_dsyrk / cblas_dgemm */
#endif

/* ---------------------------------------------------------------- vector core */

/* A rank-1 float64 NDArray as a raw contiguous buffer, or NULL if not that. */
static const double* nd_real_vec(const Expr* a, size_t* n) {
    if (!a || a->type != EXPR_NDARRAY) return NULL;
    if (a->data.ndarray.rank != 1 || a->data.ndarray.dtype != NDT_FLOAT64) return NULL;
    *n = (size_t)a->data.ndarray.dims[0];
    return (const double*)a->data.ndarray.data;
}

/* Threaded Sum of a contiguous buffer (one private slot per chunk). */
typedef struct { const double* x; } nd_sum_ctx;
static void nd_sum_chunk(void* c, size_t lo, size_t hi, double* slot) {
    const double* x = ((const nd_sum_ctx*)c)->x;
    double s = 0.0;
    for (size_t i = lo; i < hi; i++) s += x[i];
    slot[0] = s;
}
static double nd_vec_mean(const double* x, size_t n) {
    nd_sum_ctx c = { x };
    double slots[NDARRAY_MAX_THREADS];
    int k = nd_parallel_reduce(n, nd_sum_chunk, &c, 1, slots);
    double acc = 0.0;
    for (int t = 0; t < k; t++) acc += slots[t];
    return acc / (double)n;
}

/* Threaded three-way centered reduction: Sum(dx*dy), Sum(dx*dx), Sum(dy*dy)
 * where dx = x-mx, dy = y-my. Covariance needs only the first; Correlation
 * needs all three (the (n-1) divisor cancels). */
typedef struct { const double* x; const double* y; double mx, my; } nd_cent_ctx;
static void nd_cent_chunk(void* c, size_t lo, size_t hi, double* slot) {
    const nd_cent_ctx* d = (const nd_cent_ctx*)c;
    double sxy = 0.0, sxx = 0.0, syy = 0.0;
    for (size_t i = lo; i < hi; i++) {
        double dx = d->x[i] - d->mx, dy = d->y[i] - d->my;
        sxy += dx * dy; sxx += dx * dx; syy += dy * dy;
    }
    slot[0] = sxy; slot[1] = sxx; slot[2] = syy;
}
/* Fills sums[0..2] = {Sxy, Sxx, Syy}. */
static void nd_centered_sums(const double* x, const double* y, size_t n,
                             double mx, double my, double sums[3]) {
    nd_cent_ctx c = { x, y, mx, my };
    double slots[NDARRAY_MAX_THREADS * 3];
    int k = nd_parallel_reduce(n, nd_cent_chunk, &c, 3, slots);
    sums[0] = sums[1] = sums[2] = 0.0;
    for (int t = 0; t < k; t++) {
        sums[0] += slots[t * 3 + 0];
        sums[1] += slots[t * 3 + 1];
        sums[2] += slots[t * 3 + 2];
    }
}

/* Covariance[v, w] on two rank-1 real vectors -> Real, or NULL to signal the
 * caller should degrade (wrong dtype / rank, unequal length, n < 2). */
static Expr* cov_vec(const Expr* av, const Expr* aw) {
    size_t nv, nw;
    const double* v = nd_real_vec(av, &nv);
    const double* w = nd_real_vec(aw, &nw);
    if (!v || !w || nv != nw || nv < 2) return NULL;
    double mv = nd_vec_mean(v, nv), mw = nd_vec_mean(w, nv);
    double s[3];
    nd_centered_sums(v, w, nv, mv, mw, s);
    return expr_new_real(s[0] / (double)(nv - 1));
}

/* Correlation[v, w] = Cov/(sigma_v sigma_w); the (n-1) cancels. NULL to degrade
 * (wrong shape, or a zero-variance vector whose correlation is Indeterminate). */
static Expr* corr_vec(const Expr* av, const Expr* aw) {
    size_t nv, nw;
    const double* v = nd_real_vec(av, &nv);
    const double* w = nd_real_vec(aw, &nw);
    if (!v || !w || nv != nw || nv < 2) return NULL;
    double mv = nd_vec_mean(v, nv), mw = nd_vec_mean(w, nv);
    double s[3];
    nd_centered_sums(v, w, nv, mv, mw, s);
    if (s[1] <= 0.0 || s[2] <= 0.0) return NULL;   /* Indeterminate: let generic path answer */
    return expr_new_real(s[0] / sqrt(s[1] * s[2]));
}

/* ---------------------------------------------------------------- matrix core */

/* Subtract each column's mean from that column (row-major n x cols, in place). */
static void center_columns(double* A, int n, int cols) {
    for (int j = 0; j < cols; j++) {
        double s = 0.0;
        for (int i = 0; i < n; i++) s += A[(size_t)i * cols + j];
        double m = s / (double)n;
        for (int i = 0; i < n; i++) A[(size_t)i * cols + j] -= m;
    }
}

/* C(p x q, row-major) = alpha * Ac^T Bc, Ac row-major n x p, Bc row-major n x q.
 * `symmetric` (auto-covariance, B == A, p == q) uses the symmetric-rank-k
 * kernel — half the flops and a guaranteed-symmetric result. */
static void gram(const double* A, const double* B, double* C,
                 int n, int p, int q, double alpha, int symmetric) {
#ifdef USE_LAPACK
    if (symmetric) {
        cblas_dsyrk(CblasRowMajor, CblasUpper, CblasTrans,
                    p, n, alpha, A, p, 0.0, C, p);
        for (int i = 0; i < p; i++)                  /* mirror upper -> lower */
            for (int j = 0; j < i; j++)
                C[(size_t)i * p + j] = C[(size_t)j * p + i];
        return;
    }
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                p, q, n, alpha, A, p, B, q, 0.0, C, q);
#else
    (void)symmetric;
    for (int i = 0; i < p; i++)
        for (int j = 0; j < q; j++) {
            double s = 0.0;
            for (int k = 0; k < n; k++)
                s += A[(size_t)k * p + i] * B[(size_t)k * q + j];
            C[(size_t)i * q + j] = alpha * s;
        }
#endif
}

/* Column sum-of-squares of a centered row-major n x cols matrix:
 * ss[j] = Sum_i A[i][j]^2 (the (n-1)*variance of column j). */
static void col_sumsq(const double* A, int n, int cols, double* ss) {
    for (int j = 0; j < cols; j++) {
        double s = 0.0;
        for (int i = 0; i < n; i++) {
            double x = A[(size_t)i * cols + j];
            s += x * x;
        }
        ss[j] = s;
    }
}

/* Presentation source: the first NDArray argument (there is always one on the
 * buffer path), so the result inherits packed-vs-visible transparency. */
static const Expr* first_ndarray(const Expr* res) {
    for (size_t i = 0; i < res->data.function.arg_count; i++) {
        const Expr* a = res->data.function.args[i];
        if (a && a->type == EXPR_NDARRAY) return a;
    }
    return NULL;
}

/* Shared matrix driver. `correlate` selects Correlation (normalize by the outer
 * product of column standard deviations) over Covariance. */
static Expr* matrix_form(Expr* res, const Expr* a, const Expr* b, int correlate) {
    int an, ap, bn, bq;
    double *A = NULL, *B = NULL;
    int same = (a == b);
    if (!na_load_matrix(a, false, false, &an, &ap, &A))
        return ndarray_delist_and_reeval(res);
    if (same) { B = A; bn = an; bq = ap; }
    else if (!na_load_matrix(b, false, false, &bn, &bq, &B)) {
        free(A);
        return ndarray_delist_and_reeval(res);
    }
    if (an != bn || an < 2) {
        free(A); if (!same) free(B);
        return ndarray_delist_and_reeval(res);
    }
    int n = an, p = ap, q = bq;

    center_columns(A, n, p);
    if (!same) center_columns(B, n, q);

    double* C = malloc(sizeof(double) * (size_t)p * (size_t)q);
    if (!C) {
        free(A);
        if (!same) free(B);
        return ndarray_delist_and_reeval(res);
    }

    if (!correlate) {
        gram(A, B, C, n, p, q, 1.0 / (double)(n - 1), same);
    } else {
        gram(A, B, C, n, p, q, 1.0, same);              /* raw cross-products */
        double* ssa = malloc(sizeof(double) * (size_t)p);
        double* ssb = same ? ssa : malloc(sizeof(double) * (size_t)q);
        int degenerate = (!ssa || !ssb);
        if (!degenerate) {
            col_sumsq(A, n, p, ssa);
            if (!same) col_sumsq(B, n, q, ssb);
            for (int i = 0; i < p && !degenerate; i++)
                if (ssa[i] <= 0.0) degenerate = 1;      /* zero variance -> Indeterminate */
            for (int j = 0; j < q && !degenerate; j++)
                if (ssb[j] <= 0.0) degenerate = 1;
        }
        if (degenerate) {                               /* degrade to the exact path */
            free(ssa);
            if (!same) free(ssb);
            free(C);
            free(A);
            if (!same) free(B);
            return ndarray_delist_and_reeval(res);
        }
        for (int i = 0; i < p; i++)
            for (int j = 0; j < q; j++) {
                if (same && i == j) { C[(size_t)i * q + j] = 1.0; continue; }
                C[(size_t)i * q + j] /= sqrt(ssa[i] * ssb[j]);
            }
        free(ssa);
        if (!same) free(ssb);
    }

    free(A);
    if (!same) free(B);

    const Expr* src = first_ndarray(res);
    int64_t dims[2] = { p, q };
    /* Inherit the input's presentation (packed List -> matrix, visible -> NDArray). */
    Expr* out = src ? expr_new_ndarray_like(src, 2, dims, C, NDT_FLOAT64)
                    : na_build_matrix(C, p, q, false, false);
    if (src) return out;          /* expr_new_ndarray_like adopted C */
    free(C);                      /* na_build_matrix copied C */
    return out;
}

/* --------------------------------------------------------------- entry points */

/* Any argument a real (float32/float64) NDArray? A machine Real makes the whole
 * result Real, so the buffer path is exact-for-real; an int64/bool-only set of
 * operands has an exact Rational answer and must take the generic List path. */
static int has_real_float(const Expr* res) {
    for (size_t i = 0; i < res->data.function.arg_count; i++) {
        const Expr* a = res->data.function.args[i];
        if (a && a->type == EXPR_NDARRAY &&
            (a->data.ndarray.dtype == NDT_FLOAT64 ||
             a->data.ndarray.dtype == NDT_FLOAT32))
            return 1;
    }
    return 0;
}

static int is_rank1_nd(const Expr* e) {
    return e && e->type == EXPR_NDARRAY && e->data.ndarray.rank == 1;
}

static Expr* corrcov_dispatch(Expr* res, int correlate) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (!has_real_float(res)) return ndarray_delist_and_reeval(res);

    if (argc == 2) {
        Expr* a = res->data.function.args[0];
        Expr* b = res->data.function.args[1];
        if (is_rank1_nd(a) && is_rank1_nd(b)) {         /* vector . vector -> scalar */
            Expr* r = correlate ? corr_vec(a, b) : cov_vec(a, b);
            return r ? r : ndarray_delist_and_reeval(res);
        }
        return matrix_form(res, a, b, correlate);       /* matrix x matrix -> matrix */
    }
    if (argc == 1) {
        Expr* a = res->data.function.args[0];
        return matrix_form(res, a, a, correlate);        /* auto-(co)variance */
    }
    return ndarray_delist_and_reeval(res);
}

Expr* nd_covariance(Expr* res)  { return corrcov_dispatch(res, 0); }
Expr* nd_correlation(Expr* res) { return corrcov_dispatch(res, 1); }
