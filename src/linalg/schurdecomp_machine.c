/*
 * schurdecomp_machine.c -- machine-precision SchurDecomposition via LAPACK.
 *
 *   standard, real     -> dgees   (real quasi-triangular Schur form)
 *   standard, complex  -> zgees   (complex upper-triangular Schur form; also the
 *                                   RealBlockDiagonalForm -> False path for a
 *                                   real matrix, whose real buffer is widened)
 *   generalized, real  -> dgges   (QZ, real quasi-triangular S / T)
 *   generalized, cplx  -> zgges   (QZ, complex upper-triangular S / T)
 *   Pivoting -> True   -> dgebal balances the matrix, dgebak reconstructs the
 *                         scaling/permutation d = P*D by back-transforming the
 *                         identity, and the Schur form is computed on the
 *                         balanced matrix, so m . d == d . q . t . q^H.
 *
 * The matrix is loaded REAL-FIRST (schur_load_cm with want_complex=false), which
 * hits na_load_matrix's float64 memcpy/transpose fast path for a packed/NDArray
 * real input; a genuinely complex matrix fails the real load and is re-loaded
 * complex.  Factors are built with na_build_matrix_as, which keeps them packed
 * (float64 / complex64 NDArray, no per-element Expr boxing) and stamps the
 * presentation from na_result_presentation(m): a visible NDArray input yields
 * visible factors, a transparent packed-list or boxed-List input yields
 * transparent packed-list factors that thread correctly in reconstruction.
 *
 * Memory contract: standard builtin ownership (SPEC.md §4).  Never frees the
 * input; every malloc'd buffer is released on every exit path.
 */

#include "schurdecomp.h"
#include "schurdecomp_internal.h"
#include "numarray.h"
#include "lapack.h"
#include "sym_names.h"
#include "eval.h"
#include "expr.h"

#include <stdlib.h>
#include <stdbool.h>
#include <stddef.h>

/* Load an n x n numeric matrix (List or NDArray) into a column-major double
 * buffer (real -> n*n, complex -> 2*n*n interleaved).  On a first-pass miss the
 * matrix is numericalised with N[] and retried, so symbolic-but-numeric
 * constants (Pi, E, Sqrt[2], ...) load; a genuinely symbolic matrix still
 * fails.  Returns true with *buf owned by the caller, or false with *buf NULL. */
static bool schur_load_cm(const Expr* m, bool want_complex, int n, double** buf) {
    *buf = NULL;
    int rows = 0, cols = 0;
    if (na_load_matrix((Expr*)m, want_complex, true, &rows, &cols, buf)
        && rows == n && cols == n)
        return true;
    if (*buf) { free(*buf); *buf = NULL; }

    Expr* nm = eval_and_free(expr_new_function(expr_new_symbol("N"),
                             (Expr*[]){ expr_copy((Expr*)m) }, 1));
    rows = cols = 0;
    bool ok = na_load_matrix(nm, want_complex, true, &rows, &cols, buf)
              && rows == n && cols == n;
    if (!ok && *buf) { free(*buf); *buf = NULL; }
    expr_free(nm);
    return ok;
}

/* Widen a real column-major buffer to interleaved complex (im = 0). */
static double* schur_widen(const double* re, size_t nn) {
    double* z = (double*)malloc(2 * nn * sizeof(double));
    if (!z) return NULL;
    for (size_t k = 0; k < nn; k++) { z[2 * k] = re[k]; z[2 * k + 1] = 0.0; }
    return z;
}

/* Build a factor: a packed float64/complex64 NDArray carrying `pres`. */
static Expr* schur_build(const double* buf, int n, bool is_complex,
                         bool colmajor, NDPresentation pres) {
    return na_build_matrix_as(buf, n, n, is_complex, colmajor, pres);
}

/* Assemble a List result from an already-owned array of factor Exprs. */
static Expr* schur_list(Expr** items, size_t k) {
    return expr_new_function(expr_new_symbol(SYM_List), items, k);
}

/* Build the balancing d = P*D by back-transforming the identity (Pivoting). */
static double* schur_pivot_matrix(int n, int ilo, int ihi, const double* scale,
                                  bool is_complex) {
    size_t nn = (size_t)n * (size_t)n;
    size_t comps = is_complex ? 2 : 1;
    double* D = (double*)malloc(comps * nn * sizeof(double));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            size_t off = comps * ((size_t)i + (size_t)j * (size_t)n);
            D[off] = (i == j) ? 1.0 : 0.0;
            if (is_complex) D[off + 1] = 0.0;
        }
    if (is_complex) mat_lapack_zgebak('B', 'R', n, ilo, ihi, scale, n, D, n);
    else            mat_lapack_dgebak('B', 'R', n, ilo, ihi, scale, n, D, n);
    return D;
}

/* ------------------------------------------------------------------------ *
 *  Standard Schur.                                                          *
 * ------------------------------------------------------------------------ */

Expr* schur_machine_standard(const Expr* m, int n, const SchurOpts* opts) {
    if (!mathilda_lapack_probe()) return NULL;

    NDPresentation pres = na_result_presentation(m);
    size_t nn = (size_t)n * (size_t)n;

    /* Real-first load (fast path); fall back to complex for a complex matrix. */
    double* Ar = NULL;  /* real col-major (dgees overwrites it with T)        */
    double* Az = NULL;  /* complex col-major (zgees overwrites it with T)     */
    bool is_real = schur_load_cm(m, /*want_complex=*/false, n, &Ar);
    if (!is_real && !schur_load_cm(m, /*want_complex=*/true, n, &Az))
        return NULL;    /* non-numeric leaf -> leave the call unevaluated */

    bool use_real = is_real && opts->real_block_diagonal_form;
    Expr* result = NULL;

    if (use_real) {
        double* wr = (double*)malloc(sizeof(double) * (size_t)n);
        double* wi = (double*)malloc(sizeof(double) * (size_t)n);
        double* VS = (double*)malloc(sizeof(double) * nn);
        double* D  = NULL;

        int ok = 1;
        if (opts->pivoting) {
            int ilo = 1, ihi = n;
            double* scale = (double*)malloc(sizeof(double) * (size_t)n);
            if (mat_lapack_dgebal('B', n, Ar, n, &ilo, &ihi, scale) == 0)
                D = schur_pivot_matrix(n, ilo, ihi, scale, false);
            else
                ok = 0;
            free(scale);
        }

        if (ok && mat_lapack_dgees(n, Ar, n, wr, wi, VS, n) == 0) {
            Expr* q = schur_build(VS, n, false, true, pres);
            Expr* t = schur_build(Ar, n, false, true, pres);
            if (opts->pivoting && D) {
                Expr* d = schur_build(D, n, false, true, pres);
                Expr* items[3] = { q, t, d };
                result = schur_list(items, 3);
            } else {
                Expr* items[2] = { q, t };
                result = schur_list(items, 2);
            }
        }

        free(wr); free(wi); free(VS);
        if (D) free(D);
    } else {
        /* Complex driver: genuinely complex, or RealBlockDiagonalForm -> False
         * on a real matrix (widen the real buffer to interleaved complex). */
        if (!Az) Az = schur_widen(Ar, nn);
        double* w  = (double*)malloc(sizeof(double) * 2 * (size_t)n);
        double* VS = (double*)malloc(sizeof(double) * 2 * nn);
        double* D  = NULL;

        int ok = (Az != NULL);
        if (ok && opts->pivoting) {
            int ilo = 1, ihi = n;
            double* scale = (double*)malloc(sizeof(double) * (size_t)n);
            if (mat_lapack_zgebal('B', n, Az, n, &ilo, &ihi, scale) == 0)
                D = schur_pivot_matrix(n, ilo, ihi, scale, true);
            else
                ok = 0;
            free(scale);
        }

        if (ok && mat_lapack_zgees(n, Az, n, w, VS, n) == 0) {
            Expr* q = schur_build(VS, n, true, true, pres);
            Expr* t = schur_build(Az, n, true, true, pres);
            if (opts->pivoting && D) {
                Expr* d = schur_build(D, n, true, true, pres);
                Expr* items[3] = { q, t, d };
                result = schur_list(items, 3);
            } else {
                Expr* items[2] = { q, t };
                result = schur_list(items, 2);
            }
        }

        free(w); free(VS);
        if (D) free(D);
    }

    free(Ar); free(Az);
    return result;
}

/* ------------------------------------------------------------------------ *
 *  Generalized (QZ) Schur.                                                  *
 * ------------------------------------------------------------------------ */

Expr* schur_machine_generalized(const Expr* m, const Expr* a, int n,
                                const SchurOpts* opts) {
    if (!mathilda_lapack_probe()) return NULL;

    NDPresentation pres = na_result_presentation(m);
    size_t nn = (size_t)n * (size_t)n;

    double *Ar = NULL, *Br = NULL;   /* real col-major (dgges overwrites)   */
    double *Az = NULL, *Bz = NULL;   /* complex col-major (zgges overwrites) */
    bool m_real = schur_load_cm(m, false, n, &Ar);
    if (!m_real && !schur_load_cm(m, true, n, &Az)) return NULL;
    bool a_real = schur_load_cm(a, false, n, &Br);
    if (!a_real && !schur_load_cm(a, true, n, &Bz)) { free(Ar); free(Az); return NULL; }

    bool use_real = m_real && a_real && opts->real_block_diagonal_form;
    Expr* result = NULL;

    if (use_real) {
        double* alphar = (double*)malloc(sizeof(double) * (size_t)n);
        double* alphai = (double*)malloc(sizeof(double) * (size_t)n);
        double* beta   = (double*)malloc(sizeof(double) * (size_t)n);
        double* VSL = (double*)malloc(sizeof(double) * nn);
        double* VSR = (double*)malloc(sizeof(double) * nn);

        if (mat_lapack_dgges(n, Ar, n, Br, n, alphar, alphai, beta,
                             VSL, n, VSR, n) == 0) {
            Expr* q = schur_build(VSL, n, false, true, pres);
            Expr* s = schur_build(Ar,  n, false, true, pres);
            Expr* p = schur_build(VSR, n, false, true, pres);
            Expr* t = schur_build(Br,  n, false, true, pres);
            Expr* items[4] = { q, s, p, t };
            result = schur_list(items, 4);
        }

        free(alphar); free(alphai); free(beta); free(VSL); free(VSR);
    } else {
        if (!Az) Az = schur_widen(Ar, nn);
        if (!Bz) Bz = schur_widen(Br, nn);
        double* alpha = (double*)malloc(sizeof(double) * 2 * (size_t)n);
        double* beta  = (double*)malloc(sizeof(double) * 2 * (size_t)n);
        double* VSL = (double*)malloc(sizeof(double) * 2 * nn);
        double* VSR = (double*)malloc(sizeof(double) * 2 * nn);

        if (Az && Bz &&
            mat_lapack_zgges(n, Az, n, Bz, n, alpha, beta, VSL, n, VSR, n) == 0) {
            Expr* q = schur_build(VSL, n, true, true, pres);
            Expr* s = schur_build(Az,  n, true, true, pres);
            Expr* p = schur_build(VSR, n, true, true, pres);
            Expr* t = schur_build(Bz,  n, true, true, pres);
            Expr* items[4] = { q, s, p, t };
            result = schur_list(items, 4);
        }

        free(alpha); free(beta); free(VSL); free(VSR);
    }

    free(Ar); free(Br); free(Az); free(Bz);
    return result;
}
