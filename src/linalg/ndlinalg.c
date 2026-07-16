/*
 * src/linalg/ndlinalg.c
 *
 * NDArray fast-path implementations for the linear-algebra builtins. See
 * ndlinalg.h for the contract. Every fast path:
 *   - loads the flat double buffer with numarray.c's na_load_* (which has its
 *     own EXPR_NDARRAY fast path: ndt_get straight into a col-/row-major buffer,
 *     widening float32 and marshalling complex),
 *   - computes with an in-house double kernel or a mat_lapack_* routine,
 *   - rebuilds an NDArray result via na_build_* / na_scalar (closed system),
 *   - or defers to linalg_delist_and_reeval(res) for shapes/dtypes it does not
 *     specialise (the routine's ordinary numeric-List path then runs).
 *
 * The real Det/Inverse/LinearSolve/MatrixRank paths use an in-house partial-
 * pivot double LU (nd_lu_real) so they stay O(n^3) in every build, with or
 * without LAPACK. Complex heavy paths use mat_lapack_z* when the runtime probe
 * reports LAPACK, else they defer.
 */

#include "ndlinalg.h"
#include "numarray.h"
#include "ndarray.h"
#include "lapack.h"
#include "eval.h"
#include "print.h"
#include "sym_names.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------------------------------------------------------------------ */
/*  Generic guard helpers                                              */
/* ------------------------------------------------------------------ */

bool linalg_call_has_ndarray(const Expr* res)
{
    if (!res || res->type != EXPR_FUNCTION) return false;
    for (size_t i = 0; i < res->data.function.arg_count; i++) {
        if (is_ndarray(res->data.function.args[i])) return true;
    }
    return false;
}

Expr* linalg_delist_and_reeval(Expr* res)
{
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;
    Expr** args = (n > 0) ? (Expr**)malloc(sizeof(Expr*) * n) : NULL;
    for (size_t i = 0; i < n; i++) {
        Expr* a = res->data.function.args[i];
        args[i] = is_ndarray(a) ? ndarray_to_nested_list(a) : expr_copy(a);
    }
    Expr* call = expr_new_function(expr_copy(res->data.function.head), args, n);
    if (args) free(args);
    /* The rebuilt call has no NDArray arguments, so evaluating it cannot
     * re-enter a linalg_call_has_ndarray guard: no recursion. */
    return eval_and_free(call);
}

/* Does the sole matrix/vector argument carry a complex dtype? */
static bool nd_arg_is_complex(const Expr* e)
{
    return is_ndarray(e) && ndt_is_complex(e->data.ndarray.dtype);
}

/* ------------------------------------------------------------------ */
/*  In-house real double LU (column-major, partial pivoting)           */
/* ------------------------------------------------------------------ */

/* Factor the n x n column-major matrix A in place. piv[k] is the row swapped
 * into position k. *sign accumulates the permutation sign. Returns the number
 * of nonzero pivots (== n for a nonsingular matrix); a zero pivot leaves the
 * remaining column untouched and does not increment the count. */
static int nd_lu_real(double* A, int n, int* piv, int* sign)
{
    *sign = 1;
    int rank = 0;
    for (int k = 0; k < n; k++) {
        int p = k;
        double mx = fabs(A[k + (size_t)k * n]);
        for (int i = k + 1; i < n; i++) {
            double v = fabs(A[i + (size_t)k * n]);
            if (v > mx) { mx = v; p = i; }
        }
        piv[k] = p;
        if (mx == 0.0) continue;              /* singular column */
        rank++;
        if (p != k) {
            for (int j = 0; j < n; j++) {
                double t = A[k + (size_t)j * n];
                A[k + (size_t)j * n] = A[p + (size_t)j * n];
                A[p + (size_t)j * n] = t;
            }
            *sign = -*sign;
        }
        double akk = A[k + (size_t)k * n];
        for (int i = k + 1; i < n; i++) A[i + (size_t)k * n] /= akk;
        for (int j = k + 1; j < n; j++) {
            double akj = A[k + (size_t)j * n];
            if (akj != 0.0)
                for (int i = k + 1; i < n; i++)
                    A[i + (size_t)j * n] -= A[i + (size_t)k * n] * akj;
        }
    }
    return rank;
}

/* Solve L U x = P b in place on one column-major rhs of length n, using the
 * factors produced by nd_lu_real. */
static void nd_lu_solve_col(const double* LU, int n, const int* piv, double* b)
{
    for (int k = 0; k < n; k++) {                 /* apply P */
        if (piv[k] != k) { double t = b[k]; b[k] = b[piv[k]]; b[piv[k]] = t; }
    }
    for (int i = 0; i < n; i++)                    /* forward: unit lower L */
        for (int k = 0; k < i; k++)
            b[i] -= LU[i + (size_t)k * n] * b[k];
    for (int i = n - 1; i >= 0; i--) {             /* back: upper U */
        for (int k = i + 1; k < n; k++)
            b[i] -= LU[i + (size_t)k * n] * b[k];
        b[i] /= LU[i + (size_t)i * n];
    }
}

/* ------------------------------------------------------------------ */
/*  Det                                                                */
/* ------------------------------------------------------------------ */

Expr* ndla_det(Expr* res)
{
    if (res->data.function.arg_count != 1) return linalg_delist_and_reeval(res);
    Expr* arg = res->data.function.args[0];
    if (!is_ndarray(arg) || arg->data.ndarray.rank != 2
        || arg->data.ndarray.dims[0] != arg->data.ndarray.dims[1]
        || arg->data.ndarray.dims[0] == 0) {
        /* non-square / wrong rank: let the real builtin emit Det::matsq */
        return linalg_delist_and_reeval(res);
    }
    bool cplx = nd_arg_is_complex(arg);

    if (!cplx) {
        int n, cc; double* A = NULL;
        if (!na_load_matrix(arg, false, true, &n, &cc, &A))
            return linalg_delist_and_reeval(res);
        int* piv = (int*)malloc(sizeof(int) * (size_t)n);
        int sign;
        nd_lu_real(A, n, piv, &sign);
        double det = (double)sign;
        for (int k = 0; k < n; k++) det *= A[k + (size_t)k * n];
        if (det == 0.0) det = 0.0;               /* normalise -0.0 to +0.0 */
        free(piv); free(A);
        return expr_new_real(det);
    }

    /* Complex: use LAPACK zgetrf when available, else defer. */
    if (!mathilda_lapack_probe()) return linalg_delist_and_reeval(res);
    int n, cc; double* A = NULL;
    if (!na_load_matrix(arg, true, true, &n, &cc, &A))
        return linalg_delist_and_reeval(res);
    int* piv = (int*)malloc(sizeof(int) * (size_t)n);
    int info = mat_lapack_zgetrf(n, n, A, n, piv);
    if (info < 0) { free(piv); free(A); return linalg_delist_and_reeval(res); }
    double dr = 1.0, di = 0.0;                      /* running complex product */
    for (int k = 0; k < n; k++) {
        double ur = A[2 * ((size_t)k + (size_t)k * n)];
        double ui = A[2 * ((size_t)k + (size_t)k * n) + 1];
        double nr = dr * ur - di * ui;
        double ni = dr * ui + di * ur;
        dr = nr; di = ni;
        if (piv[k] != k + 1) { dr = -dr; di = -di; }  /* ipiv is 1-indexed */
    }
    free(piv); free(A);
    return na_scalar(dr, di);
}

/* ------------------------------------------------------------------ */
/*  Inverse                                                            */
/* ------------------------------------------------------------------ */

Expr* ndla_inverse(Expr* res)
{
    if (res->data.function.arg_count != 1) return linalg_delist_and_reeval(res);
    Expr* arg = res->data.function.args[0];
    if (!is_ndarray(arg) || arg->data.ndarray.rank != 2
        || arg->data.ndarray.dims[0] != arg->data.ndarray.dims[1]
        || arg->data.ndarray.dims[0] == 0)
        return linalg_delist_and_reeval(res);
    bool cplx = nd_arg_is_complex(arg);

    if (!cplx) {
        int n, cc; double* A = NULL;
        if (!na_load_matrix(arg, false, true, &n, &cc, &A))
            return linalg_delist_and_reeval(res);
        int* piv = (int*)malloc(sizeof(int) * (size_t)n);
        int sign, rank = nd_lu_real(A, n, piv, &sign);
        if (rank < n) {                              /* singular: defer for Inverse::sing */
            free(piv); free(A);
            return linalg_delist_and_reeval(res);
        }
        /* Solve A X = I column by column (col-major result). */
        double* X = (double*)malloc(sizeof(double) * (size_t)n * (size_t)n);
        double* col = (double*)malloc(sizeof(double) * (size_t)n);
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) col[i] = (i == j) ? 1.0 : 0.0;
            nd_lu_solve_col(A, n, piv, col);
            for (int i = 0; i < n; i++) X[i + (size_t)j * n] = col[i];
        }
        Expr* out = na_build_matrix(X, n, n, false, true); /* col-major -> NDArray */
        free(col); free(X); free(piv); free(A);
        return out;
    }

    if (!mathilda_lapack_probe()) return linalg_delist_and_reeval(res);
    int n, cc; double* A = NULL;
    if (!na_load_matrix(arg, true, true, &n, &cc, &A))
        return linalg_delist_and_reeval(res);
    int* piv = (int*)malloc(sizeof(int) * (size_t)n);
    int info = mat_lapack_zgetrf(n, n, A, n, piv);
    if (info != 0) { free(piv); free(A); return linalg_delist_and_reeval(res); }
    info = mat_lapack_zgetri(n, A, n, piv);
    free(piv);
    if (info != 0) { free(A); return linalg_delist_and_reeval(res); }
    Expr* out = na_build_matrix(A, n, n, true, true);
    free(A);
    return out;
}

/* ------------------------------------------------------------------ */
/*  LinearSolve                                                        */
/* ------------------------------------------------------------------ */

/* Load the rhs (vector or matrix) of a solve into a col-major buffer with n
 * rows. Returns false (and defers via the caller) on shape/dtype mismatch.
 * *is_vec records whether the rhs was rank-1 (so the result is a vector). */
static bool nd_load_rhs(const Expr* b, bool want_complex, int n,
                        int* nrhs, bool* is_vec, double** buf)
{
    int bn = 0, bc = 0;
    if (is_ndarray(b) && b->data.ndarray.rank == 1) {
        if (!na_load_vector(b, want_complex, &bn, buf)) return false;
        if (bn != n) { free(*buf); *buf = NULL; return false; }
        *nrhs = 1; *is_vec = true; return true;
    }
    if (is_ndarray(b) && b->data.ndarray.rank == 2) {
        if (!na_load_matrix(b, want_complex, true, &bn, &bc, buf)) return false;
        if (bn != n) { free(*buf); *buf = NULL; return false; }
        *nrhs = bc; *is_vec = false; return true;
    }
    return false;   /* rhs is a plain List: defer to keep behaviour uniform */
}

Expr* ndla_linearsolve(Expr* res)
{
    if (res->data.function.arg_count != 2) return linalg_delist_and_reeval(res);
    Expr* m = res->data.function.args[0];
    Expr* b = res->data.function.args[1];
    if (!is_ndarray(m) || m->data.ndarray.rank != 2
        || m->data.ndarray.dims[0] != m->data.ndarray.dims[1]
        || m->data.ndarray.dims[0] == 0)
        return linalg_delist_and_reeval(res);
    bool cplx = nd_arg_is_complex(m) || nd_arg_is_complex(b);

    if (!cplx) {
        int n, cc; double* A = NULL;
        if (!na_load_matrix(m, false, true, &n, &cc, &A))
            return linalg_delist_and_reeval(res);
        int nrhs; bool is_vec; double* B = NULL;
        if (!nd_load_rhs(b, false, n, &nrhs, &is_vec, &B)) {
            free(A); return linalg_delist_and_reeval(res);
        }
        int* piv = (int*)malloc(sizeof(int) * (size_t)n);
        int sign, rank = nd_lu_real(A, n, piv, &sign);
        if (rank < n) {                              /* singular: defer */
            free(piv); free(A); free(B);
            return linalg_delist_and_reeval(res);
        }
        for (int j = 0; j < nrhs; j++) nd_lu_solve_col(A, n, piv, B + (size_t)j * n);
        Expr* out = is_vec ? na_build_vector(B, n, false)
                           : na_build_matrix(B, n, nrhs, false, true);
        free(piv); free(A); free(B);
        return out;
    }

    if (!mathilda_lapack_probe()) return linalg_delist_and_reeval(res);
    int n, cc; double* A = NULL;
    if (!na_load_matrix(m, true, true, &n, &cc, &A))
        return linalg_delist_and_reeval(res);
    int nrhs; bool is_vec; double* B = NULL;
    if (!nd_load_rhs(b, true, n, &nrhs, &is_vec, &B)) {
        free(A); return linalg_delist_and_reeval(res);
    }
    int* piv = (int*)malloc(sizeof(int) * (size_t)n);
    int info = mat_lapack_zgesv(n, nrhs, A, n, piv, B, n);
    free(piv); free(A);
    if (info != 0) { free(B); return linalg_delist_and_reeval(res); }
    Expr* out = is_vec ? na_build_vector(B, n, true)
                       : na_build_matrix(B, n, nrhs, true, true);
    free(B);
    return out;
}

/* ------------------------------------------------------------------ */
/*  MatrixRank (real; complex defers)                                  */
/* ------------------------------------------------------------------ */

Expr* ndla_matrixrank(Expr* res)
{
    if (res->data.function.arg_count != 1) return linalg_delist_and_reeval(res);
    Expr* arg = res->data.function.args[0];
    if (!is_ndarray(arg) || arg->data.ndarray.rank != 2 || nd_arg_is_complex(arg))
        return linalg_delist_and_reeval(res);

    int r, c; double* A = NULL;                     /* row-major r x c */
    if (!na_load_matrix(arg, false, false, &r, &c, &A))
        return linalg_delist_and_reeval(res);

    double maxabs = 0.0;
    for (size_t i = 0; i < (size_t)r * (size_t)c; i++) {
        double v = fabs(A[i]); if (v > maxabs) maxabs = v;
    }
    double tol = (double)((r > c) ? r : c) * 2.220446049250313e-16 * maxabs;

    int rank = 0;
    for (int col = 0, row = 0; col < c && row < r; col++) {
        int piv = row; double mx = fabs(A[(size_t)row * c + col]);
        for (int i = row + 1; i < r; i++) {
            double v = fabs(A[(size_t)i * c + col]);
            if (v > mx) { mx = v; piv = i; }
        }
        if (mx <= tol) continue;                     /* no pivot in this column */
        if (piv != row)
            for (int j = 0; j < c; j++) {
                double t = A[(size_t)row * c + j];
                A[(size_t)row * c + j] = A[(size_t)piv * c + j];
                A[(size_t)piv * c + j] = t;
            }
        double d = A[(size_t)row * c + col];
        for (int i = row + 1; i < r; i++) {
            double f = A[(size_t)i * c + col] / d;
            if (f != 0.0)
                for (int j = col; j < c; j++)
                    A[(size_t)i * c + j] -= f * A[(size_t)row * c + j];
        }
        rank++; row++;
    }
    free(A);
    return expr_new_integer(rank);
}

/* ------------------------------------------------------------------ */
/*  Tr (2-D, default Plus over the diagonal)                           */
/* ------------------------------------------------------------------ */

Expr* ndla_tr(Expr* res)
{
    if (res->data.function.arg_count != 1) return linalg_delist_and_reeval(res);
    Expr* arg = res->data.function.args[0];
    if (!is_ndarray(arg) || arg->data.ndarray.rank != 2)
        return linalg_delist_and_reeval(res);
    int r = (int)arg->data.ndarray.dims[0];
    int c = (int)arg->data.ndarray.dims[1];
    int d = (r < c) ? r : c;
    NDType dt = arg->data.ndarray.dtype;
    double sr = 0.0, si = 0.0;
    for (int i = 0; i < d; i++) {
        double re, im;
        ndt_get(arg->data.ndarray.data, (size_t)i * (size_t)c + (size_t)i, dt, &re, &im);
        sr += re; si += im;
    }
    return na_scalar(sr, si);
}

/* ------------------------------------------------------------------ */
/*  Norm (vectors: default 2-norm, integer/real p, or Infinity)        */
/* ------------------------------------------------------------------ */

Expr* ndla_norm(Expr* res)
{
    size_t argc = res->data.function.arg_count;
    if (argc < 1 || argc > 2) return linalg_delist_and_reeval(res);
    Expr* v = res->data.function.args[0];
    if (!is_ndarray(v) || v->data.ndarray.rank != 1)
        return linalg_delist_and_reeval(res);        /* matrix/spectral norm defers */

    /* Resolve p: default 2, positive Integer/Real, or Infinity. */
    double p = 2.0; bool p_inf = false;
    if (argc == 2) {
        Expr* pe = res->data.function.args[1];
        if (pe->type == EXPR_SYMBOL && pe->data.symbol.name == SYM_Infinity) p_inf = true;
        else if (pe->type == EXPR_INTEGER && pe->data.integer > 0) p = (double)pe->data.integer;
        else if (pe->type == EXPR_REAL && pe->data.real > 0.0) p = pe->data.real;
        else return linalg_delist_and_reeval(res);
    }

    int n; double* buf = NULL;
    bool cplx = nd_arg_is_complex(v);
    if (!na_load_vector(v, cplx, &n, &buf))
        return linalg_delist_and_reeval(res);

    double out;
    if (p_inf) {
        out = 0.0;
        for (int i = 0; i < n; i++) {
            double mag = cplx ? hypot(buf[2 * i], buf[2 * i + 1]) : fabs(buf[i]);
            if (mag > out) out = mag;
        }
    } else if (p == 2.0) {
        double s = 0.0;
        for (int i = 0; i < n; i++) {
            double mag = cplx ? hypot(buf[2 * i], buf[2 * i + 1]) : fabs(buf[i]);
            s += mag * mag;
        }
        out = sqrt(s);
    } else {
        double s = 0.0;
        for (int i = 0; i < n; i++) {
            double mag = cplx ? hypot(buf[2 * i], buf[2 * i + 1]) : fabs(buf[i]);
            s += pow(mag, p);
        }
        out = pow(s, 1.0 / p);
    }
    free(buf);
    return expr_new_real(out);
}

/* ------------------------------------------------------------------ */
/*  Normalize (vector / its 2-norm)                                    */
/* ------------------------------------------------------------------ */

Expr* ndla_normalize(Expr* res)
{
    if (res->data.function.arg_count != 1) return linalg_delist_and_reeval(res);
    Expr* v = res->data.function.args[0];
    if (!is_ndarray(v) || v->data.ndarray.rank != 1)
        return linalg_delist_and_reeval(res);
    bool cplx = nd_arg_is_complex(v);
    int n; double* buf = NULL;
    if (!na_load_vector(v, cplx, &n, &buf))
        return linalg_delist_and_reeval(res);

    double s = 0.0;
    for (int i = 0; i < n; i++) {
        double mag = cplx ? hypot(buf[2 * i], buf[2 * i + 1]) : fabs(buf[i]);
        s += mag * mag;
    }
    double nrm = sqrt(s);
    if (nrm > 0.0) {
        int comps = cplx ? 2 * n : n;
        for (int i = 0; i < comps; i++) buf[i] /= nrm;
    }
    Expr* out = na_build_vector(buf, n, cplx);
    free(buf);
    return out;
}

/* ------------------------------------------------------------------ */
/*  Cross (binary, length-3 vectors)                                   */
/* ------------------------------------------------------------------ */

Expr* ndla_cross(Expr* res)
{
    if (res->data.function.arg_count != 2) return linalg_delist_and_reeval(res);
    Expr* u = res->data.function.args[0];
    Expr* v = res->data.function.args[1];
    if (!is_ndarray(u) || !is_ndarray(v)
        || u->data.ndarray.rank != 1 || v->data.ndarray.rank != 1
        || u->data.ndarray.dims[0] != 3 || v->data.ndarray.dims[0] != 3)
        return linalg_delist_and_reeval(res);
    bool cplx = nd_arg_is_complex(u) || nd_arg_is_complex(v);
    int un, vn; double* U = NULL; double* V = NULL;
    if (!na_load_vector(u, cplx, &un, &U)) return linalg_delist_and_reeval(res);
    if (!na_load_vector(v, cplx, &vn, &V)) { free(U); return linalg_delist_and_reeval(res); }

    double out[6];
    if (!cplx) {
        out[0] = U[1] * V[2] - U[2] * V[1];
        out[1] = U[2] * V[0] - U[0] * V[2];
        out[2] = U[0] * V[1] - U[1] * V[0];
    } else {
        /* complex cross via interleaved (re,im); (a*b) is complex multiply */
        #define CR(i) U[2*(i)]
        #define CI(i) U[2*(i)+1]
        #define DR(i) V[2*(i)]
        #define DI(i) V[2*(i)+1]
        /* component k = U[a]*V[b] - U[b]*V[a] */
        int a[3] = {1, 2, 0}, bb[3] = {2, 0, 1};
        for (int k = 0; k < 3; k++) {
            int i = a[k], j = bb[k];
            double r1 = CR(i) * DR(j) - CI(i) * DI(j);
            double m1 = CR(i) * DI(j) + CI(i) * DR(j);
            double r2 = CR(j) * DR(i) - CI(j) * DI(i);
            double m2 = CR(j) * DI(i) + CI(j) * DR(i);
            out[2 * k]     = r1 - r2;
            out[2 * k + 1] = m1 - m2;
        }
        #undef CR
        #undef CI
        #undef DR
        #undef DI
    }
    free(U); free(V);
    return na_build_vector(out, 3, cplx);
}
