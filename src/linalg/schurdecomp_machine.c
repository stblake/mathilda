/*
 * schurdecomp_machine.c -- machine-precision SchurDecomposition via LAPACK.
 *
 *   standard, real     -> dgees   (real quasi-triangular Schur form)
 *   standard, complex  -> zgees   (complex upper-triangular Schur form; also the
 *                                   RealBlockDiagonalForm -> False path for a
 *                                   real matrix, whose input is loaded complex)
 *   generalized, real  -> dgges   (QZ, real quasi-triangular S / T)
 *   generalized, cplx  -> zgges   (QZ, complex upper-triangular S / T)
 *   Pivoting -> True   -> dgebal balances the matrix, dgebak reconstructs the
 *                         scaling/permutation d = P*D by back-transforming the
 *                         identity, and the Schur form is computed on the
 *                         balanced matrix, so m . d == d . q . t . q^H.
 *
 * Every matrix argument (NDArray or boxed List-of-Lists) is marshalled through
 * numarray.c's na_load_matrix into a Fortran column-major double buffer;
 * schur_load_cm additionally numericalises symbolic constants (Pi, E, ...) via
 * N[] so a numeric-but-not-yet-evaluated matrix still loads.  Every factor is
 * rebuilt with schur_build, which matches the result representation to the
 * input: a boxed-List input yields boxed nested Lists (so List arithmetic in a
 * reconstruction threads element-wise), an NDArray input yields an NDArray.
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
#include "pack.h"
#include "ndarray.h"
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

/* Rebuild an n x n factor.  na_build_matrix hands back an NDArray for real data;
 * when the input was a plain boxed List we unpack it so the whole result stays
 * boxed and threads correctly against the user's boxed matrix in a
 * reconstruction (a plain List minus an atomic NDArray mis-threads). */
static Expr* schur_build(const double* buf, int n, bool is_complex,
                         bool colmajor, bool packed) {
    Expr* mm = na_build_matrix(buf, n, n, is_complex, colmajor);
    if (!packed && mm && is_ndarray(mm)) {
        Expr* boxed = pack_unpack(mm);
        expr_free(mm);
        return boxed;
    }
    return mm;
}

/* Assemble a List result from an already-owned array of factor Exprs. */
static Expr* schur_list(Expr** items, size_t k) {
    return expr_new_function(expr_new_symbol(SYM_List), items, k);
}

/* ------------------------------------------------------------------------ *
 *  Standard Schur.                                                          *
 * ------------------------------------------------------------------------ */

Expr* schur_machine_standard(const Expr* m, int n, const SchurOpts* opts) {
    if (!mathilda_lapack_probe()) return NULL;

    double* Az = NULL;   /* complex column-major, 2*n*n interleaved */
    if (!schur_load_cm(m, /*want_complex=*/true, n, &Az))
        return NULL;     /* non-numeric leaf -> leave the call unevaluated */

    bool packed = is_ndarray(m);
    size_t nn = (size_t)n * (size_t)n;

    /* Real iff every imaginary component is exactly zero. */
    bool all_imag_zero = true;
    for (size_t k = 0; k < nn; k++)
        if (Az[2 * k + 1] != 0.0) { all_imag_zero = false; break; }

    bool use_real = all_imag_zero && opts->real_block_diagonal_form;
    Expr* result = NULL;

    if (use_real) {
        double* A  = (double*)malloc(sizeof(double) * nn);        /* -> Schur T */
        double* wr = (double*)malloc(sizeof(double) * (size_t)n);
        double* wi = (double*)malloc(sizeof(double) * (size_t)n);
        double* VS = (double*)malloc(sizeof(double) * nn);        /* Schur q    */
        double* D  = NULL;                                        /* pivoting d */
        for (size_t k = 0; k < nn; k++) A[k] = Az[2 * k];

        int ok = 1;
        if (opts->pivoting) {
            int ilo = 1, ihi = n;
            double* scale = (double*)malloc(sizeof(double) * (size_t)n);
            if (mat_lapack_dgebal('B', n, A, n, &ilo, &ihi, scale) == 0) {
                D = (double*)malloc(sizeof(double) * nn);
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                        D[(size_t)i + (size_t)j * (size_t)n] = (i == j) ? 1.0 : 0.0;
                if (mat_lapack_dgebak('B', 'R', n, ilo, ihi, scale, n, D, n) != 0)
                    ok = 0;
            } else {
                ok = 0;
            }
            free(scale);
        }

        if (ok && mat_lapack_dgees(n, A, n, wr, wi, VS, n) == 0) {
            Expr* q = schur_build(VS, n, false, true, packed);
            Expr* t = schur_build(A,  n, false, true, packed);
            if (opts->pivoting && D) {
                Expr* d = schur_build(D, n, false, true, packed);
                Expr* items[3] = { q, t, d };
                result = schur_list(items, 3);
            } else {
                Expr* items[2] = { q, t };
                result = schur_list(items, 2);
            }
        }

        free(A); free(wr); free(wi); free(VS);
        if (D) free(D);
    } else {
        /* Complex path (genuinely complex input, or RealBlockDiagonalForm ->
         * False on a real matrix, which asks for complex upper-triangular t). */
        double* w  = (double*)malloc(sizeof(double) * 2 * (size_t)n);
        double* VS = (double*)malloc(sizeof(double) * 2 * nn);
        double* D  = NULL;

        int ok = 1;
        if (opts->pivoting) {
            int ilo = 1, ihi = n;
            double* scale = (double*)malloc(sizeof(double) * (size_t)n);
            if (mat_lapack_zgebal('B', n, Az, n, &ilo, &ihi, scale) == 0) {
                D = (double*)malloc(sizeof(double) * 2 * nn);
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++) {
                        size_t off = (size_t)i + (size_t)j * (size_t)n;
                        D[2 * off]     = (i == j) ? 1.0 : 0.0;
                        D[2 * off + 1] = 0.0;
                    }
                if (mat_lapack_zgebak('B', 'R', n, ilo, ihi, scale, n, D, n) != 0)
                    ok = 0;
            } else {
                ok = 0;
            }
            free(scale);
        }

        if (ok && mat_lapack_zgees(n, Az, n, w, VS, n) == 0) {
            Expr* q = schur_build(VS, n, true, true, packed);
            Expr* t = schur_build(Az, n, true, true, packed);
            if (opts->pivoting && D) {
                Expr* d = schur_build(D, n, true, true, packed);
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

    free(Az);
    return result;
}

/* ------------------------------------------------------------------------ *
 *  Generalized (QZ) Schur.                                                  *
 * ------------------------------------------------------------------------ */

Expr* schur_machine_generalized(const Expr* m, const Expr* a, int n,
                                const SchurOpts* opts) {
    if (!mathilda_lapack_probe()) return NULL;

    double* Az = NULL;   /* m, complex column-major */
    double* Bz = NULL;   /* a, complex column-major */
    if (!schur_load_cm(m, true, n, &Az)) return NULL;
    if (!schur_load_cm(a, true, n, &Bz)) { free(Az); return NULL; }

    bool packed = is_ndarray(m) || is_ndarray(a);
    size_t nn = (size_t)n * (size_t)n;

    bool all_imag_zero = true;
    for (size_t k = 0; k < nn; k++)
        if (Az[2 * k + 1] != 0.0 || Bz[2 * k + 1] != 0.0) { all_imag_zero = false; break; }

    bool use_real = all_imag_zero && opts->real_block_diagonal_form;
    Expr* result = NULL;

    if (use_real) {
        double* A = (double*)malloc(sizeof(double) * nn);  /* -> S */
        double* B = (double*)malloc(sizeof(double) * nn);  /* -> T */
        double* alphar = (double*)malloc(sizeof(double) * (size_t)n);
        double* alphai = (double*)malloc(sizeof(double) * (size_t)n);
        double* beta   = (double*)malloc(sizeof(double) * (size_t)n);
        double* VSL = (double*)malloc(sizeof(double) * nn);  /* q */
        double* VSR = (double*)malloc(sizeof(double) * nn);  /* p */
        for (size_t k = 0; k < nn; k++) { A[k] = Az[2 * k]; B[k] = Bz[2 * k]; }

        if (mat_lapack_dgges(n, A, n, B, n, alphar, alphai, beta,
                             VSL, n, VSR, n) == 0) {
            Expr* q = schur_build(VSL, n, false, true, packed);
            Expr* s = schur_build(A,   n, false, true, packed);
            Expr* p = schur_build(VSR, n, false, true, packed);
            Expr* t = schur_build(B,   n, false, true, packed);
            Expr* items[4] = { q, s, p, t };
            result = schur_list(items, 4);
        }

        free(A); free(B); free(alphar); free(alphai); free(beta);
        free(VSL); free(VSR);
    } else {
        double* alpha = (double*)malloc(sizeof(double) * 2 * (size_t)n);
        double* beta  = (double*)malloc(sizeof(double) * 2 * (size_t)n);
        double* VSL = (double*)malloc(sizeof(double) * 2 * nn);
        double* VSR = (double*)malloc(sizeof(double) * 2 * nn);

        if (mat_lapack_zgges(n, Az, n, Bz, n, alpha, beta, VSL, n, VSR, n) == 0) {
            Expr* q = schur_build(VSL, n, true, true, packed);
            Expr* s = schur_build(Az,  n, true, true, packed);
            Expr* p = schur_build(VSR, n, true, true, packed);
            Expr* t = schur_build(Bz,  n, true, true, packed);
            Expr* items[4] = { q, s, p, t };
            result = schur_list(items, 4);
        }

        free(alpha); free(beta); free(VSL); free(VSR);
    }

    free(Az); free(Bz);
    return result;
}
