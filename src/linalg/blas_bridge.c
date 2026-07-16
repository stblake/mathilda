/*
 * src/linalg/blas_bridge.c
 *
 * `BLAS`` context builtins — direct REPL access to the BLAS kernels via the
 * platform CBLAS interface (Apple Accelerate or a system cblas). See
 * blas_bridge.h for the design and numarray.h for the marshalling contract.
 *
 * CBLAS is used in row-major mode (CblasRowMajor), which matches NDArray's
 * native layout, so no transpose is needed on the BLAS side (unlike the
 * column-major Fortran LAPACK path in lapack_bridge.c).
 *
 * Each builtin loads its array arguments into fresh double buffers, calls the
 * kernel (which may overwrite one of those owned buffers — never the caller's
 * Expr), builds the result, frees the buffers, and returns. Any shape mismatch
 * or non-numeric leaf yields NULL, leaving the call unevaluated per the builtin
 * contract.
 */

#include "blas_bridge.h"
#include "numarray.h"
#include "symtab.h"
#include "attr.h"
#include "expr.h"

#include <stdlib.h>

#ifdef USE_LAPACK

#include "lapack.h"   /* pulls in Accelerate.h / cblas.h (CBLAS declarations) */

/* CBLAS's complex routines take their complex operands as void-pointers on a
 * system cblas, but as `_Complex double*` on Apple Accelerate
 * (ACCELERATE_NEW_LAPACK). An interleaved (re, im) double buffer is
 * layout-compatible with either; routing it through a `void*` lets the
 * implicit void*<->object-pointer conversion keep both toolchains free of
 * -Wincompatible-pointer-types warnings. */
static void* vp(double* p) { return (void*)p; }

/* Read a scalar argument as a real (rejecting a nonzero imaginary part). */
static int rd_real(const Expr* e, double* v)
{
    double re, im;
    if (!na_read_scalar(e, &re, &im) || im != 0.0) return 0;
    *v = re;
    return 1;
}

/* Read a scalar argument as a complex (re, im) pair. */
static int rd_cplx(const Expr* e, double out[2])
{
    double re, im;
    if (!na_read_scalar(e, &re, &im)) return 0;
    out[0] = re; out[1] = im;
    return 1;
}

static Expr* arg(Expr* res, size_t i) { return res->data.function.args[i]; }
static int argc_is(const Expr* res, size_t n)
{
    return res->type == EXPR_FUNCTION && res->data.function.arg_count == n;
}

/* ================================================================== */
/*  Level 1 — vector/vector                                            */
/* ================================================================== */

static Expr* builtin_blas_ddot(Expr* res)
{
    if (!argc_is(res, 2)) return NULL;
    int nx, ny; double *x = NULL, *y = NULL;
    if (!na_load_vector(arg(res, 0), false, &nx, &x)) return NULL;
    if (!na_load_vector(arg(res, 1), false, &ny, &y)) { free(x); return NULL; }
    Expr* out = (nx == ny) ? expr_new_real(cblas_ddot(nx, x, 1, y, 1)) : NULL;
    free(x); free(y);
    return out;
}

static Expr* builtin_blas_dnrm2(Expr* res)
{
    if (!argc_is(res, 1)) return NULL;
    int n; double* x = NULL;
    if (!na_load_vector(arg(res, 0), false, &n, &x)) return NULL;
    Expr* out = expr_new_real(cblas_dnrm2(n, x, 1));
    free(x);
    return out;
}

static Expr* builtin_blas_dasum(Expr* res)
{
    if (!argc_is(res, 1)) return NULL;
    int n; double* x = NULL;
    if (!na_load_vector(arg(res, 0), false, &n, &x)) return NULL;
    Expr* out = expr_new_real(cblas_dasum(n, x, 1));
    free(x);
    return out;
}

static Expr* builtin_blas_idamax(Expr* res)
{
    if (!argc_is(res, 1)) return NULL;
    int n; double* x = NULL;
    if (!na_load_vector(arg(res, 0), false, &n, &x)) return NULL;
    /* CBLAS returns a 0-based index; report the 1-based position. */
    Expr* out = expr_new_integer((int64_t)cblas_idamax(n, x, 1) + 1);
    free(x);
    return out;
}

static Expr* builtin_blas_daxpy(Expr* res)
{
    if (!argc_is(res, 3)) return NULL;
    double alpha;
    if (!rd_real(arg(res, 0), &alpha)) return NULL;
    int nx, ny; double *x = NULL, *y = NULL;
    if (!na_load_vector(arg(res, 1), false, &nx, &x)) return NULL;
    if (!na_load_vector(arg(res, 2), false, &ny, &y)) { free(x); return NULL; }
    Expr* out = NULL;
    if (nx == ny) {
        cblas_daxpy(nx, alpha, x, 1, y, 1);
        out = na_build_vector(y, ny, false);
    }
    free(x); free(y);
    return out;
}

static Expr* builtin_blas_dscal(Expr* res)
{
    if (!argc_is(res, 2)) return NULL;
    double alpha;
    if (!rd_real(arg(res, 0), &alpha)) return NULL;
    int n; double* x = NULL;
    if (!na_load_vector(arg(res, 1), false, &n, &x)) return NULL;
    cblas_dscal(n, alpha, x, 1);
    Expr* out = na_build_vector(x, n, false);
    free(x);
    return out;
}

static Expr* blas_zdot(Expr* res, int conj)
{
    if (!argc_is(res, 2)) return NULL;
    int nx, ny; double *x = NULL, *y = NULL;
    if (!na_load_vector(arg(res, 0), true, &nx, &x)) return NULL;
    if (!na_load_vector(arg(res, 1), true, &ny, &y)) { free(x); return NULL; }
    Expr* out = NULL;
    if (nx == ny) {
        double dot[2];
        if (conj) cblas_zdotc_sub(nx, vp(x), 1, vp(y), 1, vp(dot));
        else      cblas_zdotu_sub(nx, vp(x), 1, vp(y), 1, vp(dot));
        out = na_scalar(dot[0], dot[1]);
    }
    free(x); free(y);
    return out;
}
static Expr* builtin_blas_zdotu(Expr* res) { return blas_zdot(res, 0); }
static Expr* builtin_blas_zdotc(Expr* res) { return blas_zdot(res, 1); }

static Expr* builtin_blas_dznrm2(Expr* res)
{
    if (!argc_is(res, 1)) return NULL;
    int n; double* x = NULL;
    if (!na_load_vector(arg(res, 0), true, &n, &x)) return NULL;
    Expr* out = expr_new_real(cblas_dznrm2(n, vp(x), 1));
    free(x);
    return out;
}

static Expr* builtin_blas_zaxpy(Expr* res)
{
    if (!argc_is(res, 3)) return NULL;
    double alpha[2];
    if (!rd_cplx(arg(res, 0), alpha)) return NULL;
    int nx, ny; double *x = NULL, *y = NULL;
    if (!na_load_vector(arg(res, 1), true, &nx, &x)) return NULL;
    if (!na_load_vector(arg(res, 2), true, &ny, &y)) { free(x); return NULL; }
    Expr* out = NULL;
    if (nx == ny) {
        cblas_zaxpy(nx, vp(alpha), vp(x), 1, vp(y), 1);
        out = na_build_vector(y, ny, true);
    }
    free(x); free(y);
    return out;
}

static Expr* builtin_blas_zscal(Expr* res)
{
    if (!argc_is(res, 2)) return NULL;
    double alpha[2];
    if (!rd_cplx(arg(res, 0), alpha)) return NULL;
    int n; double* x = NULL;
    if (!na_load_vector(arg(res, 1), true, &n, &x)) return NULL;
    cblas_zscal(n, vp(alpha), vp(x), 1);
    Expr* out = na_build_vector(x, n, true);
    free(x);
    return out;
}

/* ================================================================== */
/*  Level 2 — matrix/vector                                            */
/* ================================================================== */

static Expr* builtin_blas_dgemv(Expr* res)
{
    if (!argc_is(res, 5)) return NULL;
    double alpha, beta;
    if (!rd_real(arg(res, 0), &alpha)) return NULL;
    int m, n; double* A = NULL;
    if (!na_load_matrix(arg(res, 1), false, false, &m, &n, &A)) return NULL;
    int nx; double* x = NULL;
    if (!na_load_vector(arg(res, 2), false, &nx, &x)) { free(A); return NULL; }
    int ny; double* y = NULL;
    if (!rd_real(arg(res, 3), &beta)
        || !na_load_vector(arg(res, 4), false, &ny, &y)) { free(A); free(x); return NULL; }
    Expr* out = NULL;
    if (nx == n && ny == m) {
        cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, alpha, A, n, x, 1, beta, y, 1);
        out = na_build_vector(y, m, false);
    }
    free(A); free(x); free(y);
    return out;
}

static Expr* builtin_blas_zgemv(Expr* res)
{
    if (!argc_is(res, 5)) return NULL;
    double alpha[2], beta[2];
    if (!rd_cplx(arg(res, 0), alpha)) return NULL;
    int m, n; double* A = NULL;
    if (!na_load_matrix(arg(res, 1), true, false, &m, &n, &A)) return NULL;
    int nx; double* x = NULL;
    if (!na_load_vector(arg(res, 2), true, &nx, &x)) { free(A); return NULL; }
    int ny; double* y = NULL;
    if (!rd_cplx(arg(res, 3), beta)
        || !na_load_vector(arg(res, 4), true, &ny, &y)) { free(A); free(x); return NULL; }
    Expr* out = NULL;
    if (nx == n && ny == m) {
        cblas_zgemv(CblasRowMajor, CblasNoTrans, m, n, vp(alpha), vp(A), n,
                    vp(x), 1, vp(beta), vp(y), 1);
        out = na_build_vector(y, m, true);
    }
    free(A); free(x); free(y);
    return out;
}

static Expr* builtin_blas_dger(Expr* res)
{
    if (!argc_is(res, 4)) return NULL;
    double alpha;
    if (!rd_real(arg(res, 0), &alpha)) return NULL;
    int mx, ny; double *x = NULL, *y = NULL;
    if (!na_load_vector(arg(res, 1), false, &mx, &x)) return NULL;
    if (!na_load_vector(arg(res, 2), false, &ny, &y)) { free(x); return NULL; }
    int m, n; double* A = NULL;
    if (!na_load_matrix(arg(res, 3), false, false, &m, &n, &A)) { free(x); free(y); return NULL; }
    Expr* out = NULL;
    if (m == mx && n == ny) {
        cblas_dger(CblasRowMajor, m, n, alpha, x, 1, y, 1, A, n);
        out = na_build_matrix(A, m, n, false, false);
    }
    free(x); free(y); free(A);
    return out;
}

static Expr* builtin_blas_dtrmv(Expr* res)
{
    if (!argc_is(res, 2)) return NULL;
    int n, nc; double* A = NULL;
    if (!na_load_matrix(arg(res, 0), false, false, &n, &nc, &A)) return NULL;
    int nx; double* x = NULL;
    if (n != nc || !na_load_vector(arg(res, 1), false, &nx, &x)) { free(A); return NULL; }
    Expr* out = NULL;
    if (nx == n) {
        cblas_dtrmv(CblasRowMajor, CblasUpper, CblasNoTrans, CblasNonUnit,
                    n, A, n, x, 1);
        out = na_build_vector(x, n, false);
    }
    free(A); free(x);
    return out;
}

/* ================================================================== */
/*  Level 3 — matrix/matrix                                            */
/* ================================================================== */

static Expr* builtin_blas_dgemm(Expr* res)
{
    if (!argc_is(res, 5)) return NULL;
    double alpha, beta;
    if (!rd_real(arg(res, 0), &alpha)) return NULL;
    int m, k; double* A = NULL;
    if (!na_load_matrix(arg(res, 1), false, false, &m, &k, &A)) return NULL;
    int kb, n; double* B = NULL;
    if (!na_load_matrix(arg(res, 2), false, false, &kb, &n, &B)) { free(A); return NULL; }
    int mc, nc; double* C = NULL;
    if (!rd_real(arg(res, 3), &beta)
        || !na_load_matrix(arg(res, 4), false, false, &mc, &nc, &C)) { free(A); free(B); return NULL; }
    Expr* out = NULL;
    if (kb == k && mc == m && nc == n) {
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    m, n, k, alpha, A, k, B, n, beta, C, n);
        out = na_build_matrix(C, m, n, false, false);
    }
    free(A); free(B); free(C);
    return out;
}

static Expr* builtin_blas_zgemm(Expr* res)
{
    if (!argc_is(res, 5)) return NULL;
    double alpha[2], beta[2];
    if (!rd_cplx(arg(res, 0), alpha)) return NULL;
    int m, k; double* A = NULL;
    if (!na_load_matrix(arg(res, 1), true, false, &m, &k, &A)) return NULL;
    int kb, n; double* B = NULL;
    if (!na_load_matrix(arg(res, 2), true, false, &kb, &n, &B)) { free(A); return NULL; }
    int mc, nc; double* C = NULL;
    if (!rd_cplx(arg(res, 3), beta)
        || !na_load_matrix(arg(res, 4), true, false, &mc, &nc, &C)) { free(A); free(B); return NULL; }
    Expr* out = NULL;
    if (kb == k && mc == m && nc == n) {
        cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                    m, n, k, vp(alpha), vp(A), k, vp(B), n, vp(beta), vp(C), n);
        out = na_build_matrix(C, m, n, true, false);
    }
    free(A); free(B); free(C);
    return out;
}

static Expr* builtin_blas_dsymm(Expr* res)
{
    if (!argc_is(res, 5)) return NULL;
    double alpha, beta;
    if (!rd_real(arg(res, 0), &alpha)) return NULL;
    int m, ma; double* A = NULL;
    if (!na_load_matrix(arg(res, 1), false, false, &m, &ma, &A)) return NULL;
    int mb, n; double* B = NULL;
    if (m != ma || !na_load_matrix(arg(res, 2), false, false, &mb, &n, &B)) { free(A); return NULL; }
    int mc, nc; double* C = NULL;
    if (!rd_real(arg(res, 3), &beta)
        || !na_load_matrix(arg(res, 4), false, false, &mc, &nc, &C)) { free(A); free(B); return NULL; }
    Expr* out = NULL;
    if (mb == m && mc == m && nc == n) {
        cblas_dsymm(CblasRowMajor, CblasLeft, CblasUpper,
                    m, n, alpha, A, m, B, n, beta, C, n);
        out = na_build_matrix(C, m, n, false, false);
    }
    free(A); free(B); free(C);
    return out;
}

static Expr* builtin_blas_dtrsm(Expr* res)
{
    if (!argc_is(res, 3)) return NULL;
    double alpha;
    if (!rd_real(arg(res, 0), &alpha)) return NULL;
    int m, ma; double* A = NULL;
    if (!na_load_matrix(arg(res, 1), false, false, &m, &ma, &A)) return NULL;
    int mb, n; double* B = NULL;
    if (m != ma || !na_load_matrix(arg(res, 2), false, false, &mb, &n, &B)) { free(A); return NULL; }
    Expr* out = NULL;
    if (mb == m) {
        cblas_dtrsm(CblasRowMajor, CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit,
                    m, n, alpha, A, m, B, n);
        out = na_build_matrix(B, m, n, false, false);
    }
    free(A); free(B);
    return out;
}

static Expr* builtin_blas_dsyrk(Expr* res)
{
    if (!argc_is(res, 4)) return NULL;
    double alpha, beta;
    if (!rd_real(arg(res, 0), &alpha)) return NULL;
    int n, k; double* A = NULL;
    if (!na_load_matrix(arg(res, 1), false, false, &n, &k, &A)) return NULL;
    int nc, nc2; double* C = NULL;
    if (!rd_real(arg(res, 2), &beta)
        || !na_load_matrix(arg(res, 3), false, false, &nc, &nc2, &C)) { free(A); return NULL; }
    Expr* out = NULL;
    if (nc == n && nc2 == n) {
        cblas_dsyrk(CblasRowMajor, CblasUpper, CblasNoTrans,
                    n, k, alpha, A, k, beta, C, n);
        /* dsyrk writes only the upper triangle; mirror it to the lower. */
        for (int i = 0; i < n; i++)
            for (int j = 0; j < i; j++)
                C[(size_t)i * n + j] = C[(size_t)j * n + i];
        out = na_build_matrix(C, n, n, false, false);
    }
    free(A); free(C);
    return out;
}

static Expr* builtin_blas_zherk(Expr* res)
{
    if (!argc_is(res, 4)) return NULL;
    double alpha, beta;   /* zherk's alpha/beta are real */
    if (!rd_real(arg(res, 0), &alpha)) return NULL;
    int n, k; double* A = NULL;
    if (!na_load_matrix(arg(res, 1), true, false, &n, &k, &A)) return NULL;
    int nc, nc2; double* C = NULL;
    if (!rd_real(arg(res, 2), &beta)
        || !na_load_matrix(arg(res, 3), true, false, &nc, &nc2, &C)) { free(A); return NULL; }
    Expr* out = NULL;
    if (nc == n && nc2 == n) {
        cblas_zherk(CblasRowMajor, CblasUpper, CblasNoTrans,
                    n, k, alpha, vp(A), k, beta, vp(C), n);
        /* Mirror the computed upper triangle to the lower as the Hermitian
         * conjugate; force a real diagonal (zherk guarantees it). */
        for (int i = 0; i < n; i++) {
            C[2 * ((size_t)i * n + i) + 1] = 0.0;
            for (int j = 0; j < i; j++) {
                size_t lo = 2 * ((size_t)i * n + j);
                size_t up = 2 * ((size_t)j * n + i);
                C[lo]     =  C[up];
                C[lo + 1] = -C[up + 1];
            }
        }
        out = na_build_matrix(C, n, n, true, false);
    }
    free(A); free(C);
    return out;
}

/* ================================================================== */
/*  Registration                                                       */
/* ================================================================== */

static void reg(const char* name, Expr* (*fn)(Expr*), const char* doc)
{
    symtab_add_builtin(name, fn);
    symtab_get_def(name)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(name, doc);
}

void blas_bridge_init(void)
{
    /* Level 1 */
    reg("BLAS`ddot", builtin_blas_ddot,
        "BLAS`ddot[x, y] gives the dot product x.y of two equal-length real "
        "vectors via the BLAS ddot kernel. Unevaluated on a length mismatch or "
        "non-real entry.");
    reg("BLAS`dnrm2", builtin_blas_dnrm2,
        "BLAS`dnrm2[x] gives the Euclidean (2-)norm of the real vector x via "
        "the BLAS dnrm2 kernel.");
    reg("BLAS`dasum", builtin_blas_dasum,
        "BLAS`dasum[x] gives the sum of absolute values of the real vector x "
        "via the BLAS dasum kernel.");
    reg("BLAS`idamax", builtin_blas_idamax,
        "BLAS`idamax[x] gives the 1-based index of the entry of largest "
        "absolute value in the real vector x via the BLAS idamax kernel.");
    reg("BLAS`daxpy", builtin_blas_daxpy,
        "BLAS`daxpy[alpha, x, y] gives alpha x + y for a real scalar alpha and "
        "equal-length real vectors x, y via the BLAS daxpy kernel.");
    reg("BLAS`dscal", builtin_blas_dscal,
        "BLAS`dscal[alpha, x] gives alpha x for a real scalar alpha and real "
        "vector x via the BLAS dscal kernel.");
    reg("BLAS`zdotu", builtin_blas_zdotu,
        "BLAS`zdotu[x, y] gives the unconjugated dot product x.y of two "
        "equal-length complex vectors via the BLAS zdotu kernel.");
    reg("BLAS`zdotc", builtin_blas_zdotc,
        "BLAS`zdotc[x, y] gives the conjugated dot product Conjugate[x].y of "
        "two equal-length complex vectors via the BLAS zdotc kernel.");
    reg("BLAS`dznrm2", builtin_blas_dznrm2,
        "BLAS`dznrm2[x] gives the Euclidean (2-)norm of the complex vector x "
        "via the BLAS dznrm2 kernel.");
    reg("BLAS`zaxpy", builtin_blas_zaxpy,
        "BLAS`zaxpy[alpha, x, y] gives alpha x + y for a complex scalar alpha "
        "and equal-length complex vectors x, y via the BLAS zaxpy kernel.");
    reg("BLAS`zscal", builtin_blas_zscal,
        "BLAS`zscal[alpha, x] gives alpha x for a complex scalar alpha and "
        "complex vector x via the BLAS zscal kernel.");

    /* Level 2 */
    reg("BLAS`dgemv", builtin_blas_dgemv,
        "BLAS`dgemv[alpha, A, x, beta, y] gives alpha A.x + beta y for a real "
        "m*n matrix A, length-n vector x, and length-m vector y via the BLAS "
        "dgemv kernel. Unevaluated on a shape mismatch.");
    reg("BLAS`zgemv", builtin_blas_zgemv,
        "BLAS`zgemv[alpha, A, x, beta, y] gives alpha A.x + beta y for a "
        "complex m*n matrix A, length-n vector x, and length-m vector y via "
        "the BLAS zgemv kernel.");
    reg("BLAS`dger", builtin_blas_dger,
        "BLAS`dger[alpha, x, y, A] gives the rank-1 update alpha x.y^T + A for "
        "a real scalar alpha, length-m vector x, length-n vector y, and m*n "
        "matrix A via the BLAS dger kernel.");
    reg("BLAS`dtrmv", builtin_blas_dtrmv,
        "BLAS`dtrmv[A, x] gives A.x where A is treated as the upper-triangular "
        "part of the square real matrix A (non-unit diagonal) via the BLAS "
        "dtrmv kernel.");

    /* Level 3 */
    reg("BLAS`dgemm", builtin_blas_dgemm,
        "BLAS`dgemm[alpha, A, B, beta, C] gives alpha A.B + beta C for real "
        "matrices with A m*k, B k*n, C m*n via the BLAS dgemm kernel. "
        "Unevaluated on a shape mismatch or non-real entry.");
    reg("BLAS`zgemm", builtin_blas_zgemm,
        "BLAS`zgemm[alpha, A, B, beta, C] gives alpha A.B + beta C for complex "
        "matrices with A m*k, B k*n, C m*n via the BLAS zgemm kernel.");
    reg("BLAS`dsymm", builtin_blas_dsymm,
        "BLAS`dsymm[alpha, A, B, beta, C] gives alpha A.B + beta C where A is "
        "the m*m symmetric matrix given by its upper triangle, B and C are "
        "m*n, via the BLAS dsymm kernel (left side).");
    reg("BLAS`dtrsm", builtin_blas_dtrsm,
        "BLAS`dtrsm[alpha, A, B] solves A.X = alpha B for X, where A is the "
        "m*m upper-triangular (non-unit) matrix and B is m*n, via the BLAS "
        "dtrsm kernel (left side).");
    reg("BLAS`dsyrk", builtin_blas_dsyrk,
        "BLAS`dsyrk[alpha, A, beta, C] gives the symmetric rank-k update "
        "alpha A.A^T + beta C where A is n*k and C is the n*n symmetric matrix "
        "given by its upper triangle, via the BLAS dsyrk kernel; the full "
        "symmetric result is returned.");
    reg("BLAS`zherk", builtin_blas_zherk,
        "BLAS`zherk[alpha, A, beta, C] gives the Hermitian rank-k update "
        "alpha A.A^H + beta C for real scalars alpha, beta, complex n*k matrix "
        "A, and Hermitian n*n C (upper triangle), via the BLAS zherk kernel; "
        "the full Hermitian result is returned.");
}

#else /* !USE_LAPACK */

void blas_bridge_init(void) { /* no BLAS: nothing to register */ }

#endif /* USE_LAPACK */
