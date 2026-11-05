/* Mathilda -- sparse linear algebra: matrix functions of a SparseArray.
 *
 *   MatrixExp[s]     a machine SparseArray answers a SparseArray (the
 *                    reference's head); the exponential is computed on the
 *                    dense form (it is dense in general) within the densify
 *                    cap. Exact / symbolic / MPFR input keeps the dense
 *                    path's List answer (the gate).
 *   MatrixExp[s, v]  the action exp(s) v for a machine SparseArray, by the
 *                    Al-Mohy--Higham truncated Taylor method on SpMV -- no
 *                    exp(s) is formed, so it scales to any size the SpMV
 *                    does. The answer is a dense vector, as in the reference.
 *
 * The theta_m table is shared with the dense action (src/linalg/
 * matfunc_machine.c, mf_expmv). The dense version falls back to the full
 * exponential when the Taylor cost exceeds a dense expm; here that fallback
 * is not available (the dense exponential is what we avoid), so the
 * iteration always runs: its cost is m* s SpMVs with s = ceil(||A - mu I||_1
 * / theta_m*), which the trace shift keeps small for the usual
 * diffusion-type matrices.
 */
#include "sparse_la.h"
#include "sparse_internal.h"
#include "matfunc_internal.h"
#include "numarray.h"
#include "common.h"
#include "sym_names.h"
#include "symtab.h"
#include "eval.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define MF_UNIT_ROUNDOFF 1.1102230246251565e-16

static double vec_inf(const double* x, int64_t n, int w) {
    double m = 0.0;
    for (int64_t i = 0; i < n; i++) {
        double a = w == 1 ? fabs(x[i]) : hypot(x[2 * i], x[2 * i + 1]);
        if (a > m) m = a;
    }
    return m;
}

/* exp(A) b for a square machine SpMat (b promoted to A's type by the
 * caller). out has n (complex: 2n) doubles. False on non-finite data. */
static bool sp_expmv(const SpMat* A, const double* b, double* out) {
    int64_t n = A->n;
    int w = A->cplx ? 2 : 1;
    size_t len = (size_t)w * (size_t)n;
    /* mu = trace / n, and the 1-norm of A - mu I (column sums) */
    double mur = 0.0, mui = 0.0;
    for (int64_t i = 0; i < n; i++)
        for (int64_t k = A->rp[i]; k < A->rp[i + 1]; k++)
            if (A->ci[k] == i) {
                mur += A->x[w * k];
                if (w == 2) mui += A->x[2 * k + 1];
            }
    mur /= (double)n;
    mui /= (double)n;
    double* colsum = calloc((size_t)(n > 0 ? n : 1), sizeof(double));
    double *x = malloc(sizeof(double) * len), *y = malloc(sizeof(double) * len);
    if (!colsum || !x || !y) { free(colsum); free(x); free(y); return false; }
    for (int64_t i = 0; i < n; i++) {
        bool diag = false;
        for (int64_t k = A->rp[i]; k < A->rp[i + 1]; k++) {
            int64_t j = A->ci[k];
            double re = A->x[w * k], im = w == 2 ? A->x[2 * k + 1] : 0.0;
            if (j == i) { re -= mur; im -= mui; diag = true; }
            colsum[j] += hypot(re, im);
        }
        if (!diag) colsum[i] += hypot(mur, mui);
    }
    double anorm = 0.0;
    for (int64_t j = 0; j < n; j++) if (colsum[j] > anorm) anorm = colsum[j];
    free(colsum);

    int mstar = 0;
    double sbest = 1.0;
    if (anorm > 0.0) {
        double best = INFINITY;
        for (size_t t = 0; t < EXPMV_NTHETA; t++) {
            double sc = ceil(anorm / EXPMV_THETA[t]);
            double c = (double)EXPMV_M[t] * sc;
            if (c < best) { best = c; mstar = EXPMV_M[t]; sbest = sc; }
        }
    }
    if (!isfinite(sbest) || sbest > 1e9) { free(x); free(y); return false; }
    int64_t s = (int64_t)sbest;
    /* eta = exp(mu / s) */
    double er = exp(mur / (double)s), etar = er * cos(mui / (double)s), etai = er * sin(mui / (double)s);
    memcpy(out, b, sizeof(double) * len);
    memcpy(x, b, sizeof(double) * len);
    for (int64_t i = 0; i < s; i++) {
        double c1 = vec_inf(x, n, w);
        for (int j = 1; j <= mstar; j++) {
            /* y = (A - mu I) x / (s j) */
            spla_spmv(A, x, y, false);
            double scale = 1.0 / ((double)s * (double)j);
            for (int64_t k = 0; k < n; k++) {
                if (w == 1) {
                    x[k] = (y[k] - mur * x[k]) * scale;
                } else {
                    double xr = x[2 * k], xi = x[2 * k + 1];
                    x[2 * k] = (y[2 * k] - (mur * xr - mui * xi)) * scale;
                    x[2 * k + 1] = (y[2 * k + 1] - (mur * xi + mui * xr)) * scale;
                }
            }
            double c2 = vec_inf(x, n, w);
            for (size_t k = 0; k < len; k++) out[k] += x[k];
            if (c1 + c2 <= MF_UNIT_ROUNDOFF * vec_inf(out, n, w)) break;
            c1 = c2;
        }
        for (int64_t k = 0; k < n; k++) {
            if (w == 1) {
                out[k] *= etar;
            } else {
                double fr = out[2 * k], fi = out[2 * k + 1];
                out[2 * k] = fr * etar - fi * etai;
                out[2 * k + 1] = fr * etai + fi * etar;
            }
        }
        memcpy(x, out, sizeof(double) * len);
    }
    free(x);
    free(y);
    for (size_t k = 0; k < len; k++) if (!isfinite(out[k])) return false;
    return true;
}

/* Is the SparseArray a machine one (typed machine values, or boxed values
 * holding a machine number)? */
static bool sparse_is_machine(const SparseArray* sa) {
    if (sa->kind == SPV_F64 || sa->kind == SPV_C128) return true;
    if (sa->kind != SPV_EXPR) return false;
    Expr* const* v = sa->vals;
    for (int64_t k = 0; k < sa->pat->nnz; k++) {
        CommonInexactInfo inf = common_scan_inexact(v[k]);
        if (inf.has_inexact) return inf.min_bits <= 53;
    }
    return false;
}

static Expr* eval_owned(Expr* e) {
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}

static Expr* native_matrixexp(const Expr* call) {
    size_t na = call->data.function.arg_count;
    Expr* const* a = call->data.function.args;
    if (na < 1 || na > 2 || a[0]->type != EXPR_SPARSE) return NULL;
    const SparseArray* sa = sp_of(a[0]);
    if (sa->rank != 2 || sa->dims[0] != sa->dims[1] || !sparse_is_machine(sa)) return NULL;

    if (na == 1) {
        /* SparseArray[MatrixExp[Normal[s]]]: the dense exponential, answered
         * with the reference's head */
        int64_t total = sp_total(2, sa->dims);
        if (total < 0 || total > ((int64_t)1 << 27)) return NULL;
        Expr* d = sp_normal(sa);
        if (!d) return NULL;
        Expr* ea[1] = { d };
        Expr* ex = eval_owned(expr_new_function(expr_new_symbol("MatrixExp"), ea, 1));
        if (!ex || (ex->type == EXPR_FUNCTION && ex->data.function.head->type == EXPR_SYMBOL &&
                    strcmp(ex->data.function.head->data.symbol.name, "MatrixExp") == 0)) {
            expr_free(ex);
            return NULL;
        }
        Expr* sa_args[1] = { ex };
        return eval_owned(expr_new_function(expr_new_symbol("SparseArray"), sa_args, 1));
    }

    /* MatrixExp[s, v] */
    SpMat A;
    if (!spla_view(sa, true, &A)) return NULL;
    const Expr* v = a[1];
    Expr* vd = NULL;
    if (v->type == EXPR_SPARSE) {
        if (sp_of(v)->rank != 1) { spla_view_free(&A); return NULL; }
        vd = sp_normal(sp_of(v));
        if (!vd) { spla_view_free(&A); return NULL; }
        v = vd;
    }
    int n = 0;
    double* b = NULL;
    bool cplx = A.cplx;
    if (!na_load_vector(v, cplx, &n, &b)) {
        if (cplx || !na_load_vector(v, true, &n, &b)) {
            expr_free(vd);
            spla_view_free(&A);
            return NULL;
        }
        cplx = true;
    }
    if ((int64_t)n != A.n) {
        free(b);
        expr_free(vd);
        spla_view_free(&A);
        return NULL;        /* the dense path reports the shape (::lslc) */
    }
    if (cplx && !A.cplx) {
        double* x = malloc(sizeof(double) * 2 * (size_t)(A.nnz > 0 ? A.nnz : 1));
        if (!x) { free(b); expr_free(vd); spla_view_free(&A); return NULL; }
        for (int64_t k = 0; k < A.nnz; k++) { x[2 * k] = A.x[k]; x[2 * k + 1] = 0.0; }
        free(A.owned_x);
        A.owned_x = x;
        A.x = x;
        A.cplx = true;
    }
    double* out = malloc(sizeof(double) * (cplx ? 2 : 1) * (size_t)n);
    bool ok = out && sp_expmv(&A, b, out);
    Expr* r = ok ? na_build_vector_as(out, n, cplx, na_result_presentation(a[1])) : NULL;
    free(out);
    free(b);
    expr_free(vd);
    spla_view_free(&A);
    return r;
}

void sparse_la_func_init(void) {
    symtab_set_sparse_native("MatrixExp", (const void*)native_matrixexp);
}
