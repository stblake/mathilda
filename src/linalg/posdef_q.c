/*
 * src/linalg/posdef_q.c
 *
 * PositiveDefiniteMatrixQ -- the explicit positive-definiteness predicate.
 *
 * A matrix m is positive definite iff Re[Conjugate[x] . m . x] > 0 for
 * every nonzero vector x.  Equivalently, the Hermitian part
 * H = (m + ConjugateTranspose[m]) / 2 has only positive eigenvalues, i.e.
 * admits a Cholesky factorisation H = U^H U with a real positive diagonal.
 *
 * For numeric input we load the matrix into a column-major double buffer,
 * form the Hermitian part in-place into the upper triangle, and dispatch
 * to LAPACK's `dpotrf` / `zpotrf`.  Cholesky returns info == 0 iff the
 * matrix is positive definite, so the predicate is exactly that test.
 * When USE_LAPACK is unavailable we fall back to an in-house Cholesky.
 *
 * For symbolic / non-numeric input we return False -- "explicitly
 * positive definite" follows the Mathematica convention that the
 * predicate refuses to make claims it cannot prove.
 *
 * Diagnostics:
 *   - argc != 1 -> `PositiveDefiniteMatrixQ::argx` to stderr, the call
 *     is left unevaluated (mirrors SquareMatrixQ).
 * Shape rejections that return False: non-list input, empty list,
 * non-square (including ragged), and 3-D tensors.
 */

#include "linalg.h"
#include "ndlinalg.h"
#include "lapack.h"
#include "expr.h"
#include "sym_names.h"

#include <gmp.h>
#ifdef USE_MPFR
#  include <mpfr.h>
#endif

#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/* True iff `e` is the head-symbol-equals-List function form. */
static bool pdq_is_list(const Expr* e) {
    return e && e->type == EXPR_FUNCTION
           && e->data.function.head->type == EXPR_SYMBOL
           && e->data.function.head->data.symbol == SYM_List;
}

/* Coerce a leaf to a (re, im) double pair.  Recognises Integer, BigInt,
 * Real, Rational[p, q] (exact), Complex[re, im], and (under USE_MPFR)
 * MPFR.  Returns false for any other leaf, which causes the caller to
 * fall back to the symbolic / "False" path.
 *
 * Mirrors lu_mach_leaf_to_double in ludecomp_machine.c; kept local so the
 * two modules stay independent translation units and posdef_q.c can be
 * dropped without disturbing the LU pipeline. */
static bool pdq_leaf_to_double(Expr* e, double* out_re, double* out_im) {
    *out_im = 0.0;
    if (!e) return false;
    switch (e->type) {
        case EXPR_REAL:    *out_re = e->data.real;                    return true;
        case EXPR_INTEGER: *out_re = (double)e->data.integer;         return true;
        case EXPR_BIGINT:  *out_re = mpz_get_d(e->data.bigint);       return true;
#ifdef USE_MPFR
        case EXPR_MPFR:
            *out_re = mpfr_get_d(e->data.mpfr, MPFR_RNDN);
            return true;
#endif
        case EXPR_FUNCTION:
            if (e->data.function.head->type == EXPR_SYMBOL) {
                const char* h = e->data.function.head->data.symbol;
                if (h == SYM_Rational && e->data.function.arg_count == 2) {
                    double p, q, dummy;
                    if (pdq_leaf_to_double(e->data.function.args[0], &p, &dummy)
                        && pdq_leaf_to_double(e->data.function.args[1], &q, &dummy)
                        && q != 0.0) {
                        *out_re = p / q;
                        return true;
                    }
                    return false;
                }
                if (h == SYM_Complex && e->data.function.arg_count == 2) {
                    double r, i, dummy;
                    if (pdq_leaf_to_double(e->data.function.args[0], &r, &dummy)
                        && pdq_leaf_to_double(e->data.function.args[1], &i, &dummy)) {
                        *out_re = r;
                        *out_im = i;
                        return true;
                    }
                    return false;
                }
            }
            return false;
        default:
            return false;
    }
}

/* In-house Cholesky for a real symmetric n x n matrix laid out
 * column-major in A.  Reads only the upper triangle (UPLO='U').  Returns
 * 0 on success, k > 0 if the k-th leading minor is not positive definite.
 *
 * Used only when USE_LAPACK is unavailable -- otherwise dpotrf is faster
 * and benefits from BLAS3 blocking on large matrices. */
static int pdq_chol_real_inplace(double* A, int n) {
    for (int i = 0; i < n; i++) {
        double diag = A[i + i * n];
        for (int k = 0; k < i; k++) {
            double u = A[k + i * n];
            diag -= u * u;
        }
        /* NaN-safe via the negated comparison: !(>) catches NaN and -0. */
        if (!(diag > 0.0)) return i + 1;
        double u_ii = sqrt(diag);
        A[i + i * n] = u_ii;
        for (int j = i + 1; j < n; j++) {
            double s = A[i + j * n];
            for (int k = 0; k < i; k++) {
                s -= A[k + i * n] * A[k + j * n];
            }
            A[i + j * n] = s / u_ii;
        }
    }
    return 0;
}

/* In-house Cholesky for a complex Hermitian n x n matrix, A interleaved
 * (re, im) per element and column-major.  Reads only the upper triangle.
 * Same return convention as pdq_chol_real_inplace. */
static int pdq_chol_complex_inplace(double* A, int n) {
    for (int i = 0; i < n; i++) {
        /* Diagonal entry of a Hermitian matrix is real; we extract only
         * the real component and ignore any imaginary noise. */
        double diag = A[2 * (i + i * n)];
        for (int k = 0; k < i; k++) {
            double r = A[2 * (k + i * n)];
            double s = A[2 * (k + i * n) + 1];
            diag -= r * r + s * s;
        }
        if (!(diag > 0.0)) return i + 1;
        double u_ii = sqrt(diag);
        A[2 * (i + i * n)]     = u_ii;
        A[2 * (i + i * n) + 1] = 0.0;
        for (int j = i + 1; j < n; j++) {
            double sr = A[2 * (i + j * n)];
            double si = A[2 * (i + j * n) + 1];
            for (int k = 0; k < i; k++) {
                double a_r = A[2 * (k + i * n)];
                double a_i = A[2 * (k + i * n) + 1];
                double b_r = A[2 * (k + j * n)];
                double b_i = A[2 * (k + j * n) + 1];
                /* subtract conj(U[k,i]) * U[k,j]:
                 *   conj(a + i b) * (c + i d) = (a c + b d) + i(a d - b c) */
                sr -= a_r * b_r + a_i * b_i;
                si -= a_r * b_i - a_i * b_r;
            }
            A[2 * (i + j * n)]     = sr / u_ii;
            A[2 * (i + j * n) + 1] = si / u_ii;
        }
    }
    return 0;
}

Expr* builtin_positive_definite_matrix_q(Expr* res) {
    if (linalg_call_has_ndarray(res)) return linalg_delist_and_reeval(res);
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 1) {
        fprintf(stderr,
                "PositiveDefiniteMatrixQ::argx: PositiveDefiniteMatrixQ "
                "called with %zu argument%s; 1 argument is expected.\n",
                argc, argc == 1 ? "" : "s");
        return NULL;
    }

    Expr* m = res->data.function.args[0];

    /* Shape gate: must be a non-empty square List of Lists with no
     * deeper nesting.  Identical to the gate used by SquareMatrixQ. */
    if (!pdq_is_list(m)) return expr_new_symbol(SYM_False);
    size_t n_sz = m->data.function.arg_count;
    if (n_sz == 0) return expr_new_symbol(SYM_False);
    for (size_t i = 0; i < n_sz; i++) {
        Expr* row = m->data.function.args[i];
        if (!pdq_is_list(row)) return expr_new_symbol(SYM_False);
        if (row->data.function.arg_count != n_sz)
            return expr_new_symbol(SYM_False);
        for (size_t j = 0; j < n_sz; j++) {
            if (pdq_is_list(row->data.function.args[j]))
                return expr_new_symbol(SYM_False);
        }
    }
    if (n_sz > (size_t)INT_MAX) return expr_new_symbol(SYM_False);
    int n = (int)n_sz;

    /* Probe numericity and detect whether any entry is non-real.  Any
     * non-coercible leaf flips us into the "symbolic -> False" branch. */
    bool is_complex = false;
    for (int i = 0; i < n; i++) {
        Expr* row = m->data.function.args[i];
        for (int j = 0; j < n; j++) {
            double re, im;
            if (!pdq_leaf_to_double(row->data.function.args[j], &re, &im)) {
                return expr_new_symbol(SYM_False);
            }
            if (im != 0.0) is_complex = true;
        }
    }

    /* Build the Hermitian part H = (m + ConjugateTranspose[m]) / 2 into
     * a column-major double buffer.  Only the upper triangle is needed
     * by potrf (UPLO='U'); we leave the strict lower triangle filled
     * with the (untouched) input for simpler bookkeeping.
     *
     * Note: building the Hermitian part is essential because
     * PositiveDefiniteMatrixQ does NOT require the input to be Hermitian
     * -- only the Hermitian part determines Re[Conjugate[x].m.x]. */
    size_t stride = is_complex ? 2 : 1;
    double* A = (double*)malloc(stride * (size_t)n * (size_t)n
                                * sizeof(double));
    if (!A) return expr_new_symbol(SYM_False);

    for (int i = 0; i < n; i++) {
        Expr* row = m->data.function.args[i];
        for (int j = 0; j < n; j++) {
            double re, im;
            pdq_leaf_to_double(row->data.function.args[j], &re, &im);
            size_t off = stride * ((size_t)i + (size_t)j * (size_t)n);
            A[off] = re;
            if (is_complex) A[off + 1] = im;
        }
    }

    /* Symmetrise / Hermitise into the upper triangle in place: for i < j
     *   H[i,j] = (A[i,j] + Conjugate(A[j,i])) / 2
     * For the diagonal (i == j) of a Hermitian matrix we keep the real
     * part only; the imaginary component of (m + m^H)/2 vanishes by
     * construction.  potrf won't read the imaginary part of the diagonal
     * anyway, but we zero it to keep the in-house fallback deterministic. */
    if (!is_complex) {
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                double aij = A[i + j * n];
                double aji = A[j + i * n];
                A[i + j * n] = 0.5 * (aij + aji);
            }
        }
    } else {
        for (int i = 0; i < n; i++) {
            A[2 * (i + i * n) + 1] = 0.0;
            for (int j = i + 1; j < n; j++) {
                double aij_r = A[2 * (i + j * n)];
                double aij_i = A[2 * (i + j * n) + 1];
                double aji_r = A[2 * (j + i * n)];
                double aji_i = A[2 * (j + i * n) + 1];
                /* (a + Conjugate(b)) / 2 with a = A[i,j], b = A[j,i]. */
                A[2 * (i + j * n)]     = 0.5 * (aij_r + aji_r);
                A[2 * (i + j * n) + 1] = 0.5 * (aij_i - aji_i);
            }
        }
    }

    /* Necessary condition: a Hermitian PD matrix has strictly positive
     * diagonal.  Catches obvious rejections without invoking Cholesky. */
    bool diag_ok = true;
    for (int i = 0; i < n; i++) {
        double d = is_complex ? A[2 * (i + i * n)] : A[i + i * n];
        if (!(d > 0.0)) { diag_ok = false; break; }
    }
    if (!diag_ok) {
        free(A);
        return expr_new_symbol(SYM_False);
    }

    /* Dispatch to LAPACK's blocked Cholesky when available; otherwise
     * fall back to the in-house implementation.  Either returns 0 iff
     * the Hermitian part of the matrix is positive definite. */
    int info;
    if (is_complex) {
        info = mat_lapack_zpotrf('U', n, A, n);
        if (info < 0) info = pdq_chol_complex_inplace(A, n);
    } else {
        info = mat_lapack_dpotrf('U', n, A, n);
        if (info < 0) info = pdq_chol_real_inplace(A, n);
    }
    free(A);
    return expr_new_symbol(info == 0 ? "True" : "False");
}
