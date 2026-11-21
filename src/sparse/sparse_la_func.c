/* Mathilda -- sparse linear algebra: matrix functions of a SparseArray.
 *
 *   MatrixExp[s]     a machine SparseArray answers a SparseArray (the
 *                    reference's head); the exponential is computed on the
 *                    dense form (it is dense in general) within the densify
 *                    cap. Exact / symbolic input, and most MPFR input
 *                    (mpfr_answers_sparse: unless ||s - mu I||_1 <= 1/2),
 *                    keeps the dense path's List answer (the gate).
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

/* The numeric class of a SparseArray's values: 0 exact (or symbolic), 1
 * machine, 2 arbitrary precision. The reference answers MatrixExp of a
 * machine SparseArray with a SparseArray and of an exact one with a List; of
 * an MPFR one see mpfr_answers_sparse. */
static int sparse_class(const SparseArray* sa) {
    /* no stored value: a machine background (N of an exact zero matrix)
     * makes it machine -- the reference answers the identity SparseArray */
    if (sa->pat->nnz == 0 && sa->bg && sa->bg->type == EXPR_REAL) return 1;
    if (sa->kind == SPV_F64 || sa->kind == SPV_C128) return 1;
    if (sa->kind != SPV_EXPR) return 0;
    Expr* const* v = sa->vals;
    long bits = 0;
    for (int64_t k = 0; k < sa->pat->nnz; k++) {
        CommonInexactInfo inf = common_scan_inexact(v[k]);
        if (inf.has_inexact && (bits == 0 || inf.min_bits < bits)) bits = inf.min_bits;
    }
    return bits == 0 ? 0 : bits <= 53 ? 1 : 2;
}

/* An arbitrary-precision zero (0``p, or a complex one with both parts so). */
static bool mpfr_zero_entry(const Expr* e) {
    if (e && e->type == EXPR_FUNCTION && e->data.function.arg_count == 2 &&
        e->data.function.head->type == EXPR_SYMBOL &&
        strcmp(e->data.function.head->data.symbol.name, "Complex") == 0)
        return e->data.function.args[0]->type != EXPR_INTEGER && spla_is_zero(e->data.function.args[0]) &&
               e->data.function.args[1]->type != EXPR_INTEGER && spla_is_zero(e->data.function.args[1]);
    return e && e->type != EXPR_INTEGER && spla_is_zero(e);
}

/* The arbitrary-precision exponential's zeros (the structurally zero blocks
 * of a reducible matrix come out as 0``p) are not stored: the reference's
 * SparseArray holds the numerically nonzero entries over an exact 0
 * background, so they become that background. */
static Expr* mpfr_zeros_to_background(Expr* m) {
    if (!m || m->type != EXPR_FUNCTION) return m;
    m = expr_unshare(m);                        /* copy-on-write, both levels */
    for (size_t i = 0; i < m->data.function.arg_count; i++) {
        Expr* row = m->data.function.args[i];
        if (!row || row->type != EXPR_FUNCTION) continue;
        bool any = false;
        for (size_t j = 0; j < row->data.function.arg_count && !any; j++)
            any = mpfr_zero_entry(row->data.function.args[j]);
        if (!any) continue;
        row = m->data.function.args[i] = expr_unshare(row);
        for (size_t j = 0; j < row->data.function.arg_count; j++)
            if (mpfr_zero_entry(row->data.function.args[j])) {
                expr_free(row->data.function.args[j]);
                row->data.function.args[j] = expr_new_integer(0);
            }
    }
    return m;
}

static Expr* eval_owned(Expr* e) {
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}

/* MatrixExp of an MPFR SparseArray: the reference answers a SparseArray
 * exactly when the trace-shifted matrix is small,
 *
 *     || A - mu I ||_1 <= 1/2,   mu = Tr[A] / n,
 *
 * the 1-norm (largest column sum of absolute values; |z| for a complex
 * entry) compared as a machine double -- presumably the branch on which its
 * exponential needs no squaring and keeps the sparse arithmetic. Every other
 * MPFR matrix answers the dense List. Established on ~560 reference probes
 * (orders 1..12; diagonal, triangular, nilpotent, banded, block, permutation
 * and random real/complex patterns; 17..50 digits) with no exception:
 * {{0,t},{0,0}} is sparse for t = 1/2 + 10^-17 and dense for 1/2 + 10^-16
 * (the double rounding of the norm), the 1-norm and not the infinity-,
 * 2- or Frobenius norm decides ({{0,.3,.3},{0,0,0},{0,0,0}} sparse, its
 * transpose dense), and the shift is the mean, not the midpoint, of the
 * diagonal (diag(0,0,1) dense, diag(0,1/2,1) sparse). Order 1 shifts to 0
 * and is always sparse. The norm is read from a machine copy, as the
 * reference compares it: it only chooses the head. */
static bool mpfr_answers_sparse(const Expr* s) {
    Expr* na[1] = { expr_copy((Expr*)s) };
    Expr* ms = eval_owned(expr_new_function(expr_new_symbol("N"), na, 1));
    bool sparse = false;
    SpMat A;
    if (ms && ms->type == EXPR_SPARSE && spla_view(sp_of(ms), true, &A)) {   /* square: checked by the caller */
        int64_t n = A.n;
        double* col = calloc((size_t)(n > 0 ? n : 1), sizeof(double));
        bool* hasd = calloc((size_t)(n > 0 ? n : 1), sizeof(bool));
        if (col && hasd) {
            double mr = 0.0, mi = 0.0;              /* mu = Tr[A] / n */
            for (int64_t i = 0; i < n; i++)
                for (int64_t k = A.rp[i]; k < A.rp[i + 1]; k++)
                    if (A.ci[k] == i) {
                        mr += A.cplx ? A.x[2 * k] : A.x[k];
                        mi += A.cplx ? A.x[2 * k + 1] : 0.0;
                    }
            mr /= (double)n;
            mi /= (double)n;
            for (int64_t i = 0; i < n; i++)
                for (int64_t k = A.rp[i]; k < A.rp[i + 1]; k++) {
                    int64_t j = A.ci[k];
                    double re = A.cplx ? A.x[2 * k] : A.x[k], im = A.cplx ? A.x[2 * k + 1] : 0.0;
                    if (j == i) { re -= mr; im -= mi; hasd[i] = true; }
                    col[j] += hypot(re, im);
                }
            double norm = 0.0;
            for (int64_t j = 0; j < n; j++) {
                if (!hasd[j]) col[j] += hypot(mr, mi);   /* an implicit 0 on the diagonal */
                if (col[j] > norm || isnan(col[j])) norm = col[j];
            }
            sparse = norm <= 0.5;
        }
        free(col);
        free(hasd);
        spla_view_free(&A);
    }
    expr_free(ms);
    return sparse;
}

static Expr* native_matrixexp(const Expr* call) {
    size_t na = call->data.function.arg_count;
    Expr* const* a = call->data.function.args;
    if (na < 1 || na > 2 || a[0]->type != EXPR_SPARSE) return NULL;
    const SparseArray* sa = sp_of(a[0]);
    int cls = sparse_class(sa);
    /* a nonzero background: the reference answers a List (the gate's dense
     * path) */
    if (sa->rank != 2 || sa->dims[0] != sa->dims[1] || cls == 0 || !spla_is_zero(sa->bg)) return NULL;
    if (na == 1 && cls == 2 && !mpfr_answers_sparse(a[0])) return NULL;   /* the dense List */

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
        if (cls == 2) ex = mpfr_zeros_to_background(ex);
        Expr* sa_args[1] = { ex };
        return eval_owned(expr_new_function(expr_new_symbol("SparseArray"), sa_args, 1));
    }

    /* MatrixExp[s, v]: the SpMV action is machine arithmetic; an
     * arbitrary-precision s keeps its precision on the dense path */
    if (cls != 1) return NULL;
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
