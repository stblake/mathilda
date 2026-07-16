/*
 * src/linalg/lapack_bridge.c
 *
 * `LAPACK`` context builtins — direct REPL access to the LAPACK driver
 * routines via the column-major Fortran wrappers in lapack.{c,h}. See
 * lapack_bridge.h for the design and numarray.h for marshalling.
 *
 * Every array argument is loaded column-major (colmajor=true) to match the
 * Fortran ABI; results are rebuilt from the column-major output buffers. Real
 * results come back as NDArray, complex results as nested List of Complex[...]
 * (numarray.h). Any shape mismatch, non-numeric leaf, or LAPACK failure
 * (singular system, non-convergence, not positive definite) yields NULL,
 * leaving the call unevaluated per the builtin contract.
 */

#include "lapack_bridge.h"
#include "numarray.h"
#include "symtab.h"
#include "attr.h"
#include "expr.h"

#include <stdlib.h>

#ifdef USE_LAPACK

#include "lapack.h"

/* Complex stride: 2 interleaved doubles per element, 1 for real. */
#define ST(cplx) ((cplx) ? (size_t)2 : (size_t)1)

static Expr* arg(Expr* res, size_t i) { return res->data.function.args[i]; }
static int argc_is(const Expr* res, size_t n)
{
    return res->type == EXPR_FUNCTION && res->data.function.arg_count == n;
}

static Expr* make_list(Expr** el, size_t n)
{
    return expr_new_function(expr_new_symbol("List"), el, n);
}

/* Load a right-hand side that is either a length-`n` vector or an n*nrhs
 * matrix, column-major. Sets *nrhs, *is_vector, *buf (caller frees). Returns
 * false on a row-count mismatch or non-numeric leaf. */
static int load_rhs(const Expr* e, int cplx, int n,
                    int* nrhs, int* is_vector, double** buf)
{
    int vn; double* vb = NULL;
    if (na_load_vector(e, cplx, &vn, &vb)) {
        if (vn != n) { free(vb); return 0; }
        *nrhs = 1; *is_vector = 1; *buf = vb;
        return 1;
    }
    int r, c; double* mb = NULL;
    if (na_load_matrix(e, cplx, 1 /*colmajor*/, &r, &c, &mb)) {
        if (r != n) { free(mb); return 0; }
        *nrhs = c; *is_vector = 0; *buf = mb;
        return 1;
    }
    return 0;
}

static Expr* build_rhs(const double* buf, int n, int nrhs, int is_vector, int cplx)
{
    return is_vector ? na_build_vector(buf, n, cplx)
                     : na_build_matrix(buf, n, nrhs, cplx, 1);
}

/* ================================================================== */
/*  Solve                                                              */
/* ================================================================== */

static Expr* lp_gesv(Expr* res, int cplx)
{
    if (!argc_is(res, 2)) return NULL;
    int n, na; double* A = NULL;
    if (!na_load_matrix(arg(res, 0), cplx, 1, &n, &na, &A)) return NULL;
    if (n != na) { free(A); return NULL; }
    int nrhs, isv; double* B = NULL;
    if (!load_rhs(arg(res, 1), cplx, n, &nrhs, &isv, &B)) { free(A); return NULL; }

    int* ipiv = (int*)malloc((size_t)n * sizeof(int));
    Expr* out = NULL;
    if (ipiv) {
        int info = cplx ? mat_lapack_zgesv(n, nrhs, A, n, ipiv, B, n)
                        : mat_lapack_dgesv(n, nrhs, A, n, ipiv, B, n);
        if (info == 0) out = build_rhs(B, n, nrhs, isv, cplx);
    }
    free(ipiv); free(A); free(B);
    return out;
}
static Expr* builtin_lapack_dgesv(Expr* res) { return lp_gesv(res, 0); }
static Expr* builtin_lapack_zgesv(Expr* res) { return lp_gesv(res, 1); }

static Expr* lp_trtrs(Expr* res, int cplx)
{
    if (!argc_is(res, 2)) return NULL;
    int n, na; double* A = NULL;
    if (!na_load_matrix(arg(res, 0), cplx, 1, &n, &na, &A)) return NULL;
    if (n != na) { free(A); return NULL; }
    int nrhs, isv; double* B = NULL;
    if (!load_rhs(arg(res, 1), cplx, n, &nrhs, &isv, &B)) { free(A); return NULL; }

    int info = cplx ? mat_lapack_ztrtrs(n, nrhs, A, n, B, n)
                    : mat_lapack_dtrtrs(n, nrhs, A, n, B, n);
    Expr* out = (info == 0) ? build_rhs(B, n, nrhs, isv, cplx) : NULL;
    free(A); free(B);
    return out;
}
static Expr* builtin_lapack_dtrtrs(Expr* res) { return lp_trtrs(res, 0); }
static Expr* builtin_lapack_ztrtrs(Expr* res) { return lp_trtrs(res, 1); }

static Expr* lp_gels(Expr* res, int cplx)
{
    if (!argc_is(res, 2)) return NULL;
    int m, n; double* A = NULL;
    if (!na_load_matrix(arg(res, 0), cplx, 1, &m, &n, &A)) return NULL;
    int nrhs, isv; double* B = NULL;
    if (!load_rhs(arg(res, 1), cplx, m, &nrhs, &isv, &B)) { free(A); return NULL; }

    size_t st = ST(cplx);
    int ldb = (m > n) ? m : n;
    double* Bbig = (double*)calloc(st * (size_t)ldb * (size_t)nrhs, sizeof(double));
    Expr* out = NULL;
    if (Bbig) {
        for (int j = 0; j < nrhs; j++)
            for (int i = 0; i < m; i++) {
                size_t s = st * ((size_t)i + (size_t)j * (size_t)m);
                size_t d = st * ((size_t)i + (size_t)j * (size_t)ldb);
                Bbig[d] = B[s];
                if (cplx) Bbig[d + 1] = B[s + 1];
            }
        int info = cplx ? mat_lapack_zgels(m, n, nrhs, A, m, Bbig, ldb)
                        : mat_lapack_dgels(m, n, nrhs, A, m, Bbig, ldb);
        if (info == 0) {
            double* X = (double*)malloc(st * (size_t)n * (size_t)nrhs * sizeof(double));
            if (X) {
                for (int j = 0; j < nrhs; j++)
                    for (int i = 0; i < n; i++) {
                        size_t s = st * ((size_t)i + (size_t)j * (size_t)ldb);
                        size_t d = st * ((size_t)i + (size_t)j * (size_t)n);
                        X[d] = Bbig[s];
                        if (cplx) X[d + 1] = Bbig[s + 1];
                    }
                out = build_rhs(X, n, nrhs, isv, cplx);
                free(X);
            }
        }
    }
    free(Bbig); free(A); free(B);
    return out;
}
static Expr* builtin_lapack_dgels(Expr* res) { return lp_gels(res, 0); }
static Expr* builtin_lapack_zgels(Expr* res) { return lp_gels(res, 1); }

/* ================================================================== */
/*  Factorizations                                                     */
/* ================================================================== */

static Expr* lp_getrf(Expr* res, int cplx)
{
    if (!argc_is(res, 1)) return NULL;
    int m, n; double* A = NULL;
    if (!na_load_matrix(arg(res, 0), cplx, 1, &m, &n, &A)) return NULL;
    int k = (m < n) ? m : n;
    int* ipiv = (int*)malloc((size_t)k * sizeof(int));
    Expr* out = NULL;
    if (ipiv) {
        int info = cplx ? mat_lapack_zgetrf(m, n, A, m, ipiv)
                        : mat_lapack_dgetrf(m, n, A, m, ipiv);
        if (info >= 0) {                    /* >0: singular U, factors still valid */
            Expr* LU = na_build_matrix(A, m, n, cplx, 1);
            Expr** pel = (Expr**)malloc((size_t)k * sizeof(Expr*));
            for (int i = 0; i < k; i++) pel[i] = expr_new_integer(ipiv[i]);
            Expr* piv = make_list(pel, (size_t)k);
            free(pel);
            Expr* el[2] = { LU, piv };
            out = make_list(el, 2);
        }
    }
    free(ipiv); free(A);
    return out;
}
static Expr* builtin_lapack_dgetrf(Expr* res) { return lp_getrf(res, 0); }
static Expr* builtin_lapack_zgetrf(Expr* res) { return lp_getrf(res, 1); }

static Expr* lp_getri(Expr* res, int cplx)
{
    if (!argc_is(res, 1)) return NULL;
    int n, na; double* A = NULL;
    if (!na_load_matrix(arg(res, 0), cplx, 1, &n, &na, &A)) return NULL;
    Expr* out = NULL;
    if (n == na) {
        int* ipiv = (int*)malloc((size_t)n * sizeof(int));
        if (ipiv) {
            int info = cplx ? mat_lapack_zgetrf(n, n, A, n, ipiv)
                            : mat_lapack_dgetrf(n, n, A, n, ipiv);
            if (info == 0) {
                info = cplx ? mat_lapack_zgetri(n, A, n, ipiv)
                            : mat_lapack_dgetri(n, A, n, ipiv);
                if (info == 0) out = na_build_matrix(A, n, n, cplx, 1);
            }
            free(ipiv);
        }
    }
    free(A);
    return out;
}
static Expr* builtin_lapack_dgetri(Expr* res) { return lp_getri(res, 0); }
static Expr* builtin_lapack_zgetri(Expr* res) { return lp_getri(res, 1); }

static Expr* lp_geqrf(Expr* res, int cplx)
{
    if (!argc_is(res, 1)) return NULL;
    int m, n; double* A = NULL;
    if (!na_load_matrix(arg(res, 0), cplx, 1, &m, &n, &A)) return NULL;
    int k = (m < n) ? m : n;
    size_t st = ST(cplx);
    double* tau = (double*)malloc(st * (size_t)k * sizeof(double));
    double* R = (double*)calloc(st * (size_t)k * (size_t)n, sizeof(double));
    Expr* out = NULL;
    if (tau && R) {
        int info = cplx ? mat_lapack_zgeqrf(m, n, A, m, tau)
                        : mat_lapack_dgeqrf(m, n, A, m, tau);
        if (info == 0) {
            /* Extract upper-trapezoidal R (k*n) before the reflectors below the
             * diagonal are consumed by orgqr/ungqr. */
            for (int j = 0; j < n; j++)
                for (int i = 0; i <= j && i < k; i++) {
                    size_t rd = st * ((size_t)i + (size_t)j * (size_t)k);
                    size_t ad = st * ((size_t)i + (size_t)j * (size_t)m);
                    R[rd] = A[ad];
                    if (cplx) R[rd + 1] = A[ad + 1];
                }
            int info2 = cplx ? mat_lapack_zungqr(m, k, k, A, m, tau)
                             : mat_lapack_dorgqr(m, k, k, A, m, tau);
            if (info2 == 0) {
                Expr* Q  = na_build_matrix(A, m, k, cplx, 1);   /* first k cols */
                Expr* Rm = na_build_matrix(R, k, n, cplx, 1);
                Expr* el[2] = { Q, Rm };
                out = make_list(el, 2);
            }
        }
    }
    free(tau); free(R); free(A);
    return out;
}
static Expr* builtin_lapack_dgeqrf(Expr* res) { return lp_geqrf(res, 0); }
static Expr* builtin_lapack_zgeqrf(Expr* res) { return lp_geqrf(res, 1); }

static Expr* lp_potrf(Expr* res, int cplx)
{
    if (!argc_is(res, 1)) return NULL;
    int n, na; double* A = NULL;
    if (!na_load_matrix(arg(res, 0), cplx, 1, &n, &na, &A)) return NULL;
    if (n != na) { free(A); return NULL; }
    Expr* out = NULL;
    int info = cplx ? mat_lapack_zpotrf('L', n, A, n)
                    : mat_lapack_dpotrf('L', n, A, n);
    if (info == 0) {
        size_t st = ST(cplx);
        /* Zero the untouched strict-upper triangle (col-major entry (i,j), i<j). */
        for (int j = 0; j < n; j++)
            for (int i = 0; i < j; i++) {
                size_t o = st * ((size_t)i + (size_t)j * (size_t)n);
                A[o] = 0.0;
                if (cplx) A[o + 1] = 0.0;
            }
        out = na_build_matrix(A, n, n, cplx, 1);
    }
    free(A);
    return out;
}
static Expr* builtin_lapack_dpotrf(Expr* res) { return lp_potrf(res, 0); }
static Expr* builtin_lapack_zpotrf(Expr* res) { return lp_potrf(res, 1); }

/* ================================================================== */
/*  SVD                                                                */
/* ================================================================== */

static Expr* lp_svd(Expr* res, int cplx, int divide_conquer)
{
    if (!argc_is(res, 1)) return NULL;
    int m, n; double* A = NULL;
    if (!na_load_matrix(arg(res, 0), cplx, 1, &m, &n, &A)) return NULL;
    int mn = (m < n) ? m : n;
    size_t st = ST(cplx);
    double* S  = (double*)malloc((size_t)mn * sizeof(double));
    double* U  = (double*)malloc(st * (size_t)m * (size_t)m * sizeof(double));
    double* VT = (double*)malloc(st * (size_t)n * (size_t)n * sizeof(double));
    Expr* out = NULL;
    if (S && U && VT) {
        int info;
        if (divide_conquer)
            info = cplx ? mat_lapack_zgesdd('A', m, n, A, m, S, U, m, VT, n)
                        : mat_lapack_dgesdd('A', m, n, A, m, S, U, m, VT, n);
        else
            info = cplx ? mat_lapack_zgesvd(m, n, A, m, S, U, m, VT, n)
                        : mat_lapack_dgesvd(m, n, A, m, S, U, m, VT, n);
        if (info == 0) {
            Expr* Ue = na_build_matrix(U, m, m, cplx, 1);
            Expr* Se = na_build_vector(S, mn, 0);
            Expr* Ve = na_build_matrix(VT, n, n, cplx, 1);
            Expr* el[3] = { Ue, Se, Ve };
            out = make_list(el, 3);
        }
    }
    free(S); free(U); free(VT); free(A);
    return out;
}
static Expr* builtin_lapack_dgesdd(Expr* res) { return lp_svd(res, 0, 1); }
static Expr* builtin_lapack_zgesdd(Expr* res) { return lp_svd(res, 1, 1); }
static Expr* builtin_lapack_dgesvd(Expr* res) { return lp_svd(res, 0, 0); }
static Expr* builtin_lapack_zgesvd(Expr* res) { return lp_svd(res, 1, 0); }

/* ================================================================== */
/*  Eigenproblems                                                      */
/* ================================================================== */

static Expr* lp_dgeev(Expr* res)
{
    if (!argc_is(res, 1)) return NULL;
    int n, na; double* A = NULL;
    if (!na_load_matrix(arg(res, 0), 0, 1, &n, &na, &A)) return NULL;
    if (n != na) { free(A); return NULL; }

    double* wr = (double*)malloc((size_t)n * sizeof(double));
    double* wi = (double*)malloc((size_t)n * sizeof(double));
    double* VR = (double*)malloc((size_t)n * (size_t)n * sizeof(double));
    double* zb = (double*)malloc(2 * (size_t)n * sizeof(double));
    Expr* out = NULL;
    if (wr && wi && VR && zb) {
        int info = mat_lapack_dgeev(n, A, n, wr, wi, VR, n);
        if (info == 0) {
            Expr** vals = (Expr**)malloc((size_t)n * sizeof(Expr*));
            Expr** vecs = (Expr**)malloc((size_t)n * sizeof(Expr*));
            for (int j = 0; j < n; j++) {
                vals[j] = na_scalar(wr[j], wi[j]);
                for (int r = 0; r < n; r++) {
                    double re, im;
                    if (wi[j] == 0.0) {         /* real eigenvector */
                        re = VR[(size_t)r + (size_t)j * n]; im = 0.0;
                    } else if (wi[j] > 0.0) {   /* first of a conjugate pair */
                        re = VR[(size_t)r + (size_t)j * n];
                        im = VR[(size_t)r + (size_t)(j + 1) * n];
                    } else {                    /* second of the pair (conjugate) */
                        re = VR[(size_t)r + (size_t)(j - 1) * n];
                        im = -VR[(size_t)r + (size_t)j * n];
                    }
                    zb[2 * r] = re; zb[2 * r + 1] = im;
                }
                vecs[j] = na_build_vector(zb, n, 1);
            }
            Expr* el[2] = { make_list(vals, (size_t)n), make_list(vecs, (size_t)n) };
            free(vals); free(vecs);
            out = make_list(el, 2);
        }
    }
    free(wr); free(wi); free(VR); free(zb); free(A);
    return out;
}

static Expr* lp_zgeev(Expr* res)
{
    if (!argc_is(res, 1)) return NULL;
    int n, na; double* A = NULL;
    if (!na_load_matrix(arg(res, 0), 1, 1, &n, &na, &A)) return NULL;
    if (n != na) { free(A); return NULL; }

    double* w  = (double*)malloc(2 * (size_t)n * sizeof(double));
    double* VR = (double*)malloc(2 * (size_t)n * (size_t)n * sizeof(double));
    double* zb = (double*)malloc(2 * (size_t)n * sizeof(double));
    Expr* out = NULL;
    if (w && VR && zb) {
        int info = mat_lapack_zgeev(n, A, n, w, VR, n);
        if (info == 0) {
            Expr** vals = (Expr**)malloc((size_t)n * sizeof(Expr*));
            Expr** vecs = (Expr**)malloc((size_t)n * sizeof(Expr*));
            for (int j = 0; j < n; j++) {
                vals[j] = na_scalar(w[2 * j], w[2 * j + 1]);
                for (int r = 0; r < n; r++) {
                    size_t off = 2 * ((size_t)r + (size_t)j * n);
                    zb[2 * r] = VR[off]; zb[2 * r + 1] = VR[off + 1];
                }
                vecs[j] = na_build_vector(zb, n, 1);
            }
            Expr* el[2] = { make_list(vals, (size_t)n), make_list(vecs, (size_t)n) };
            free(vals); free(vecs);
            out = make_list(el, 2);
        }
    }
    free(w); free(VR); free(zb); free(A);
    return out;
}
static Expr* builtin_lapack_dgeev(Expr* res) { return lp_dgeev(res); }
static Expr* builtin_lapack_zgeev(Expr* res) { return lp_zgeev(res); }

static Expr* lp_syev(Expr* res, int cplx)
{
    if (!argc_is(res, 1)) return NULL;
    int n, na; double* A = NULL;
    if (!na_load_matrix(arg(res, 0), cplx, 1, &n, &na, &A)) return NULL;
    if (n != na) { free(A); return NULL; }

    double* w = (double*)malloc((size_t)n * sizeof(double));
    size_t st = ST(cplx);
    double* col = (double*)malloc(st * (size_t)n * sizeof(double));
    Expr* out = NULL;
    if (w && col) {
        int info = cplx ? mat_lapack_zheev(n, A, n, w)
                        : mat_lapack_dsyev(n, A, n, w);
        if (info == 0) {
            Expr** vals = (Expr**)malloc((size_t)n * sizeof(Expr*));
            Expr** vecs = (Expr**)malloc((size_t)n * sizeof(Expr*));
            for (int j = 0; j < n; j++) {
                vals[j] = expr_new_real(w[j]);
                for (int r = 0; r < n; r++) {
                    size_t o = st * ((size_t)r + (size_t)j * (size_t)n);
                    col[st * r] = A[o];
                    if (cplx) col[st * r + 1] = A[o + 1];
                }
                vecs[j] = na_build_vector(col, n, cplx);
            }
            Expr* el[2] = { make_list(vals, (size_t)n), make_list(vecs, (size_t)n) };
            free(vals); free(vecs);
            out = make_list(el, 2);
        }
    }
    free(w); free(col); free(A);
    return out;
}
static Expr* builtin_lapack_dsyev(Expr* res) { return lp_syev(res, 0); }
static Expr* builtin_lapack_zheev(Expr* res) { return lp_syev(res, 1); }

/* ================================================================== */
/*  Norm                                                               */
/* ================================================================== */

static Expr* lp_lange(Expr* res, int cplx)
{
    if (!argc_is(res, 1)) return NULL;
    int m, n; double* A = NULL;
    if (!na_load_matrix(arg(res, 0), cplx, 1, &m, &n, &A)) return NULL;
    double r = cplx ? mat_lapack_zlange('F', m, n, A, m)
                    : mat_lapack_dlange('F', m, n, A, m);
    free(A);
    return (r < 0.0) ? NULL : expr_new_real(r);
}
static Expr* builtin_lapack_dlange(Expr* res) { return lp_lange(res, 0); }
static Expr* builtin_lapack_zlange(Expr* res) { return lp_lange(res, 1); }

/* ================================================================== */
/*  Registration                                                       */
/* ================================================================== */

static void reg(const char* name, Expr* (*fn)(Expr*), const char* doc)
{
    symtab_add_builtin(name, fn);
    symtab_get_def(name)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(name, doc);
}

void lapack_bridge_init(void)
{
    reg("LAPACK`dgesv", builtin_lapack_dgesv,
        "LAPACK`dgesv[A, B] solves the square real system A.X = B by LU "
        "factorisation (LAPACK dgesv). B is a length-n vector or an n*k "
        "matrix; the result matches its shape. Unevaluated on a non-square, "
        "singular, or mismatched system.");
    reg("LAPACK`zgesv", builtin_lapack_zgesv,
        "LAPACK`zgesv[A, B] solves the square complex system A.X = B by LU "
        "factorisation (LAPACK zgesv).");
    reg("LAPACK`dtrtrs", builtin_lapack_dtrtrs,
        "LAPACK`dtrtrs[A, B] solves the triangular real system A.X = B where A "
        "is the upper-triangular part of the square matrix A (non-unit "
        "diagonal), via LAPACK dtrtrs.");
    reg("LAPACK`ztrtrs", builtin_lapack_ztrtrs,
        "LAPACK`ztrtrs[A, B] solves the complex upper-triangular system "
        "A.X = B via LAPACK ztrtrs.");
    reg("LAPACK`dgels", builtin_lapack_dgels,
        "LAPACK`dgels[A, B] gives the least-squares (m>=n) or minimum-norm "
        "(m<n) full-rank solution of A.X = B for a real m*n matrix A, via "
        "LAPACK dgels. B is a length-m vector or m*k matrix.");
    reg("LAPACK`zgels", builtin_lapack_zgels,
        "LAPACK`zgels[A, B] gives the least-squares / minimum-norm solution of "
        "the complex system A.X = B via LAPACK zgels.");

    reg("LAPACK`dgetrf", builtin_lapack_dgetrf,
        "LAPACK`dgetrf[A] gives {LU, pivots} for the real m*n matrix A: the "
        "combined unit-lower/upper LU factors and the 1-based row pivots, via "
        "LAPACK dgetrf.");
    reg("LAPACK`zgetrf", builtin_lapack_zgetrf,
        "LAPACK`zgetrf[A] gives {LU, pivots} for the complex m*n matrix A via "
        "LAPACK zgetrf.");
    reg("LAPACK`dgetri", builtin_lapack_dgetri,
        "LAPACK`dgetri[A] gives the inverse of the square real matrix A "
        "(factored internally with dgetrf, then inverted with LAPACK dgetri). "
        "Unevaluated for a singular or non-square matrix.");
    reg("LAPACK`zgetri", builtin_lapack_zgetri,
        "LAPACK`zgetri[A] gives the inverse of the square complex matrix A via "
        "dgetrf + LAPACK zgetri.");
    reg("LAPACK`dgeqrf", builtin_lapack_dgeqrf,
        "LAPACK`dgeqrf[A] gives the economy QR factorisation {Q, R} of the "
        "real m*n matrix A (Q is m*k, R is k*n upper-trapezoidal, k=Min[m,n], "
        "A == Q.R), via LAPACK dgeqrf + dorgqr.");
    reg("LAPACK`zgeqrf", builtin_lapack_zgeqrf,
        "LAPACK`zgeqrf[A] gives the economy QR factorisation {Q, R} of the "
        "complex m*n matrix A via LAPACK zgeqrf + zungqr.");
    reg("LAPACK`dpotrf", builtin_lapack_dpotrf,
        "LAPACK`dpotrf[A] gives the lower Cholesky factor L (A == L.Transpose[L]) "
        "of the real symmetric positive-definite matrix A, via LAPACK dpotrf. "
        "Unevaluated when A is not positive definite.");
    reg("LAPACK`zpotrf", builtin_lapack_zpotrf,
        "LAPACK`zpotrf[A] gives the lower Cholesky factor L "
        "(A == L.ConjugateTranspose[L]) of the Hermitian positive-definite "
        "matrix A, via LAPACK zpotrf.");

    reg("LAPACK`dgesdd", builtin_lapack_dgesdd,
        "LAPACK`dgesdd[A] gives the full SVD {U, S, VT} of the real m*n matrix "
        "A (A == U.DiagonalMatrix[S].VT) via the divide-and-conquer LAPACK "
        "dgesdd; U is m*m, VT is n*n, S has length Min[m,n].");
    reg("LAPACK`zgesdd", builtin_lapack_zgesdd,
        "LAPACK`zgesdd[A] gives the full SVD {U, S, VT} of the complex m*n "
        "matrix A via the divide-and-conquer LAPACK zgesdd.");
    reg("LAPACK`dgesvd", builtin_lapack_dgesvd,
        "LAPACK`dgesvd[A] gives the full SVD {U, S, VT} of the real m*n matrix "
        "A via the QR-iteration LAPACK dgesvd.");
    reg("LAPACK`zgesvd", builtin_lapack_zgesvd,
        "LAPACK`zgesvd[A] gives the full SVD {U, S, VT} of the complex m*n "
        "matrix A via the QR-iteration LAPACK zgesvd.");

    reg("LAPACK`dgeev", builtin_lapack_dgeev,
        "LAPACK`dgeev[A] gives {values, vectors} for the general real square "
        "matrix A: eigenvalues (real or complex) and the matching right "
        "eigenvectors, via LAPACK dgeev.");
    reg("LAPACK`zgeev", builtin_lapack_zgeev,
        "LAPACK`zgeev[A] gives {values, vectors} for the general complex "
        "square matrix A via LAPACK zgeev.");
    reg("LAPACK`dsyev", builtin_lapack_dsyev,
        "LAPACK`dsyev[A] gives {values, vectors} for the real symmetric matrix "
        "A: ascending real eigenvalues and orthonormal eigenvectors (upper "
        "triangle referenced), via LAPACK dsyev.");
    reg("LAPACK`zheev", builtin_lapack_zheev,
        "LAPACK`zheev[A] gives {values, vectors} for the complex Hermitian "
        "matrix A: ascending real eigenvalues and orthonormal eigenvectors, "
        "via LAPACK zheev.");

    reg("LAPACK`dlange", builtin_lapack_dlange,
        "LAPACK`dlange[A] gives the Frobenius norm of the real matrix A via "
        "LAPACK dlange.");
    reg("LAPACK`zlange", builtin_lapack_zlange,
        "LAPACK`zlange[A] gives the Frobenius norm of the complex matrix A via "
        "LAPACK zlange.");
}

#else /* !USE_LAPACK */

void lapack_bridge_init(void) { /* no LAPACK: nothing to register */ }

#endif /* USE_LAPACK */
