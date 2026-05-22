/* src/linalg/ludecomp_mpfr.c
 *
 * Arbitrary-precision LU kernel: Doolittle's elimination with partial
 * row pivoting over MPFR arrays.
 *
 * Invoked by lu_dispatch when common_scan_inexact reports inexact
 * input at min_bits > 53.  On any failure (USE_MPFR off, a leaf the
 * loader can't reduce to an MPFR value, a zero pivot during
 * elimination) the kernel returns NULL and lu_dispatch falls back to
 * the symbolic Doolittle pipeline.
 *
 * Algorithm.  For each step k in [0, n):
 *
 *   1. Pick the row in [k, n) with the largest |A[i, k]| as the
 *      pivot; swap into row k.
 *   2. For i in [k+1, n):
 *         A[i, k] /= A[k, k]                    (L entry)
 *         For j in [k+1, n):
 *             A[i, j] -= A[i, k] * A[k, j]      (Schur update)
 *
 * Condition-number estimate.  After the factorisation we form
 * A^{-1} = U^{-1} * L^{-1} * P via two back-substitutions then a
 * matrix product, and report ||A||_inf * ||A^{-1}||_inf.  This is
 * O(n^3) which matches the cost of the factorisation itself; far from
 * the tight LAPACK estimator but adequate for the typical
 * arbitrary-precision use case (a few dozen rows at thousands of bits).
 *
 * Memory contract.  Standard builtin contract: this file does NOT
 * call expr_free on the input `m` -- the evaluator owns it.  Every
 * mpfr_t initialised here is cleared along every exit path.
 */

#include "ludecomp_internal.h"
#include "linalg.h"
#include "expr.h"
#include "print.h"
#include "sym_names.h"
#include "common.h"
#include "numeric.h"

#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef USE_MPFR

Expr* lu_mpfr_dispatch(Expr* m, int n)
{
    (void)m; (void)n;
    return NULL;
}

#else /* USE_MPFR */

/* ------------------------------------------------------------------ *
 *  Allocation helpers.  Duplicated from qrdecomp_mpfr.c so this        *
 *  translation unit stays self-contained.                              *
 * ------------------------------------------------------------------ */
static mpfr_t* lum_array_alloc(size_t count, mpfr_prec_t bits)
{
    if (count == 0) return NULL;
    mpfr_t* a = (mpfr_t*)malloc(sizeof(mpfr_t) * count);
    for (size_t i = 0; i < count; i++) mpfr_init2(a[i], bits);
    return a;
}

static void lum_array_free(mpfr_t* a, size_t count)
{
    if (!a) return;
    for (size_t i = 0; i < count; i++) mpfr_clear(a[i]);
    free(a);
}

static bool lum_leaf_is_complex(Expr* e)
{
    if (!e) return false;
    mpfr_t r, im;
    mpfr_init2(r, 53);
    mpfr_init2(im, 53);
    bool is_inexact = false;
    bool ok = get_approx_mpfr(e, r, im, &is_inexact);
    bool result = ok && (mpfr_zero_p(im) == 0);
    mpfr_clear(r);
    mpfr_clear(im);
    return result;
}

/* Load the n x n input matrix into freshly-allocated row-major MPFR
 * arrays.  Row-major (rather than column-major) is convenient because
 * the Doolittle kernel always accesses whole rows for both pivoting
 * swaps and the elimination's outer loop. */
static bool lum_load_matrix(Expr* m, int n, mpfr_prec_t bits,
                            mpfr_t** out_re, mpfr_t** out_im,
                            bool* out_is_complex)
{
    if (m->type != EXPR_FUNCTION) return false;

    bool is_complex = false;
    for (int i = 0; i < n && !is_complex; i++) {
        Expr* row = m->data.function.args[i];
        if (!row || row->type != EXPR_FUNCTION) return false;
        for (int j = 0; j < n && !is_complex; j++) {
            if (lum_leaf_is_complex(row->data.function.args[j])) {
                is_complex = true;
            }
        }
    }

    size_t total = (size_t)n * (size_t)n;
    mpfr_t* A_re = lum_array_alloc(total, bits);
    mpfr_t* A_im = is_complex ? lum_array_alloc(total, bits) : NULL;

    mpfr_t tmp_im;
    mpfr_init2(tmp_im, bits);
    for (int i = 0; i < n; i++) {
        Expr* row = m->data.function.args[i];
        for (int j = 0; j < n; j++) {
            size_t off = (size_t)i * (size_t)n + (size_t)j;
            bool is_inexact = false;
            if (!get_approx_mpfr(row->data.function.args[j],
                                 A_re[off], tmp_im, &is_inexact)) {
                mpfr_clear(tmp_im);
                lum_array_free(A_re, total);
                lum_array_free(A_im, total);
                return false;
            }
            if (is_complex) mpfr_set(A_im[off], tmp_im, MPFR_RNDN);
        }
    }
    mpfr_clear(tmp_im);

    *out_re = A_re;
    *out_im = A_im;
    *out_is_complex = is_complex;
    return true;
}

/* Magnitude of a (re, im) cell.  Real case: |re|. */
static void lum_abs(const mpfr_t* re, const mpfr_t* im, size_t off,
                    bool is_complex, mpfr_t out)
{
    if (is_complex) mpfr_hypot(out, re[off], im[off], MPFR_RNDN);
    else            mpfr_abs(out,   re[off], MPFR_RNDN);
}

/* Swap whole rows in a row-major n x n matrix. */
static void lum_swap_rows(mpfr_t* A_re, mpfr_t* A_im, int n,
                          int a, int b, bool is_complex)
{
    if (a == b) return;
    for (int j = 0; j < n; j++) {
        mpfr_swap(A_re[(size_t)a * n + j], A_re[(size_t)b * n + j]);
        if (is_complex) mpfr_swap(A_im[(size_t)a * n + j],
                                   A_im[(size_t)b * n + j]);
    }
}

/* a /= b for a complex scalar a and b (a, b stored as (re, im) refs).
 * Uses 1 / (br + I bi) = (br - I bi) / (br^2 + bi^2). */
static void lum_cdiv(mpfr_t* ar, mpfr_t* ai,
                     const mpfr_t br, const mpfr_t bi,
                     mpfr_prec_t bits)
{
    mpfr_t den, t1, ar_in, ai_in;
    mpfr_init2(den,  bits);
    mpfr_init2(t1,   bits);
    mpfr_init2(ar_in, bits);
    mpfr_init2(ai_in, bits);
    mpfr_sqr(den, br, MPFR_RNDN);
    mpfr_sqr(t1,  bi, MPFR_RNDN);
    mpfr_add(den, den, t1, MPFR_RNDN);
    mpfr_set(ar_in, *ar, MPFR_RNDN);
    mpfr_set(ai_in, *ai, MPFR_RNDN);
    /* new ar = (ar * br + ai * bi) / den */
    mpfr_mul(t1, ar_in, br, MPFR_RNDN);
    mpfr_mul(*ar, ai_in, bi, MPFR_RNDN);
    mpfr_add(*ar, *ar, t1, MPFR_RNDN);
    mpfr_div(*ar, *ar, den, MPFR_RNDN);
    /* new ai = (ai * br - ar * bi) / den */
    mpfr_mul(t1, ai_in, br, MPFR_RNDN);
    mpfr_mul(*ai, ar_in, bi, MPFR_RNDN);
    mpfr_sub(*ai, t1, *ai, MPFR_RNDN);
    mpfr_div(*ai, *ai, den, MPFR_RNDN);
    mpfr_clear(den);
    mpfr_clear(t1);
    mpfr_clear(ar_in);
    mpfr_clear(ai_in);
}

/* a -= b * c for complex scalars.  All operands stored as (re, im)
 * refs; b and c are read-only. */
static void lum_csub_mul(mpfr_t* ar, mpfr_t* ai,
                          const mpfr_t br, const mpfr_t bi,
                          const mpfr_t cr, const mpfr_t ci,
                          mpfr_prec_t bits)
{
    mpfr_t pr, pi, t;
    mpfr_init2(pr, bits);
    mpfr_init2(pi, bits);
    mpfr_init2(t,  bits);
    /* (br + I bi)(cr + I ci) = (br cr - bi ci) + I (br ci + bi cr) */
    mpfr_mul(pr, br, cr, MPFR_RNDN);
    mpfr_mul(t,  bi, ci, MPFR_RNDN);
    mpfr_sub(pr, pr, t,  MPFR_RNDN);
    mpfr_mul(pi, br, ci, MPFR_RNDN);
    mpfr_mul(t,  bi, cr, MPFR_RNDN);
    mpfr_add(pi, pi, t,  MPFR_RNDN);
    mpfr_sub(*ar, *ar, pr, MPFR_RNDN);
    mpfr_sub(*ai, *ai, pi, MPFR_RNDN);
    mpfr_clear(pr);
    mpfr_clear(pi);
    mpfr_clear(t);
}

/* ------------------------------------------------------------------ *
 *  Doolittle elimination.  Overwrites A with the combined LU; perm    *
 *  starts as [1..n] and is updated to the pivoting permutation.       *
 *  Returns true on success.  *out_singular is set to true iff a zero  *
 *  pivot was encountered (factorisation continues; matches            *
 *  Mathematica behaviour).                                             *
 * ------------------------------------------------------------------ */
static bool lum_factor(mpfr_t* A_re, mpfr_t* A_im, int n,
                        mpfr_prec_t bits, bool is_complex,
                        int* perm, bool* out_singular)
{
    mpfr_t mag, best;
    mpfr_init2(mag,  bits);
    mpfr_init2(best, bits);
    *out_singular = false;

    for (int k = 0; k < n; k++) {
        /* Pivot: largest |A[i, k]| over i in [k, n). */
        int pivot_row = k;
        lum_abs(A_re, A_im, (size_t)k * n + k, is_complex, best);
        for (int i = k + 1; i < n; i++) {
            lum_abs(A_re, A_im, (size_t)i * n + k, is_complex, mag);
            if (mpfr_cmp(mag, best) > 0) {
                pivot_row = i;
                mpfr_set(best, mag, MPFR_RNDN);
            }
        }
        if (mpfr_zero_p(best)) {
            *out_singular = true;
            /* Skip the elimination for this step; U[k, k] stays zero. */
            continue;
        }
        if (pivot_row != k) {
            lum_swap_rows(A_re, A_im, n, k, pivot_row, is_complex);
            int tp = perm[k]; perm[k] = perm[pivot_row]; perm[pivot_row] = tp;
        }

        /* Eliminate. */
        for (int i = k + 1; i < n; i++) {
            size_t ik = (size_t)i * n + k;
            size_t kk = (size_t)k * n + k;
            if (is_complex) {
                lum_cdiv(&A_re[ik], &A_im[ik], A_re[kk], A_im[kk], bits);
            } else {
                mpfr_div(A_re[ik], A_re[ik], A_re[kk], MPFR_RNDN);
            }
            for (int j = k + 1; j < n; j++) {
                size_t ij = (size_t)i * n + j;
                size_t kj = (size_t)k * n + j;
                if (is_complex) {
                    lum_csub_mul(&A_re[ij], &A_im[ij],
                                  A_re[ik], A_im[ik],
                                  A_re[kj], A_im[kj], bits);
                } else {
                    mpfr_t prod;
                    mpfr_init2(prod, bits);
                    mpfr_mul(prod, A_re[ik], A_re[kj], MPFR_RNDN);
                    mpfr_sub(A_re[ij], A_re[ij], prod, MPFR_RNDN);
                    mpfr_clear(prod);
                }
            }
        }
    }
    mpfr_clear(mag);
    mpfr_clear(best);
    return true;
}

/* ------------------------------------------------------------------ *
 *  L-infinity norm of a row-major n x n matrix.                         *
 *      max_i sum_j |A[i, j]|                                            *
 * ------------------------------------------------------------------ */
static void lum_norm_inf(const mpfr_t* A_re, const mpfr_t* A_im, int n,
                          bool is_complex, mpfr_prec_t bits, mpfr_t out)
{
    mpfr_t row_sum, mag;
    mpfr_init2(row_sum, bits);
    mpfr_init2(mag,     bits);
    mpfr_set_zero(out, 1);
    for (int i = 0; i < n; i++) {
        mpfr_set_zero(row_sum, 1);
        for (int j = 0; j < n; j++) {
            lum_abs(A_re, A_im, (size_t)i * n + j, is_complex, mag);
            mpfr_add(row_sum, row_sum, mag, MPFR_RNDN);
        }
        if (mpfr_cmp(row_sum, out) > 0) mpfr_set(out, row_sum, MPFR_RNDN);
    }
    mpfr_clear(row_sum);
    mpfr_clear(mag);
}

/* Solve L y = e_col where L is the strictly-lower triangle of LU (with
 * unit diagonal implicit), then U x = y in place.  Writes the
 * solution vector into out_re / out_im (length n).  Used during the
 * condition-number computation. */
static void lum_solve_LU_canonical(const mpfr_t* LU_re, const mpfr_t* LU_im,
                                    int n, bool is_complex, int col,
                                    mpfr_t* out_re, mpfr_t* out_im,
                                    mpfr_prec_t bits)
{
    /* Right-hand side e_col, length n. */
    for (int i = 0; i < n; i++) {
        mpfr_set_zero(out_re[i], 1);
        if (is_complex) mpfr_set_zero(out_im[i], 1);
    }
    mpfr_set_si(out_re[col], 1, MPFR_RNDN);

    /* Forward substitution: L y = b.  L[i, j] = LU[i, j] for j < i;
     * L[i, i] = 1. */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            size_t ij = (size_t)i * n + j;
            if (is_complex) {
                lum_csub_mul(&out_re[i], &out_im[i],
                              LU_re[ij], LU_im[ij],
                              out_re[j], out_im[j], bits);
            } else {
                mpfr_t prod;
                mpfr_init2(prod, bits);
                mpfr_mul(prod, LU_re[ij], out_re[j], MPFR_RNDN);
                mpfr_sub(out_re[i], out_re[i], prod, MPFR_RNDN);
                mpfr_clear(prod);
            }
        }
    }
    /* Back substitution: U x = y.  U[i, j] = LU[i, j] for j >= i. */
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i + 1; j < n; j++) {
            size_t ij = (size_t)i * n + j;
            if (is_complex) {
                lum_csub_mul(&out_re[i], &out_im[i],
                              LU_re[ij], LU_im[ij],
                              out_re[j], out_im[j], bits);
            } else {
                mpfr_t prod;
                mpfr_init2(prod, bits);
                mpfr_mul(prod, LU_re[ij], out_re[j], MPFR_RNDN);
                mpfr_sub(out_re[i], out_re[i], prod, MPFR_RNDN);
                mpfr_clear(prod);
            }
        }
        size_t ii = (size_t)i * n + i;
        if (is_complex) {
            lum_cdiv(&out_re[i], &out_im[i], LU_re[ii], LU_im[ii], bits);
        } else {
            mpfr_div(out_re[i], out_re[i], LU_re[ii], MPFR_RNDN);
        }
    }
}

/* Compute ||A^{-1}||_inf via explicit inverse.  A^{-1}'s column k is
 * the solution of LU * x = P * e_k where P is the row permutation
 * encoded by perm (perm[i] is the source row that ended up at row i).
 * Equivalently: solve LU * x = e_{perm^{-1}[k]}, then sum row magnitudes
 * across all columns.
 *
 * Returns false if any diagonal of U is zero (singular). */
static bool lum_inv_norm_inf(const mpfr_t* LU_re, const mpfr_t* LU_im,
                              const int* perm, int n, bool is_complex,
                              mpfr_prec_t bits, mpfr_t out)
{
    /* Check for zero U diagonals -- singular if any. */
    for (int i = 0; i < n; i++) {
        size_t ii = (size_t)i * n + i;
        if (mpfr_zero_p(LU_re[ii])
            && (!is_complex || mpfr_zero_p(LU_im[ii]))) {
            return false;
        }
    }

    /* perm[i] is the original row at row i of LU.  So x = A^{-1} * b
     * is computed by setting b' = P * b (which permutes b's entries
     * by perm), then solving LU x = b'.  For b = e_k, b'[i] = 1 iff
     * perm[i] == k+1 (i.e. i = perm^{-1}[k+1]).  Equivalently, the
     * "permuted col" of the identity is the row index whose perm[i]
     * equals (k+1). */

    int* pinv = (int*)malloc((size_t)n * sizeof(int));
    for (int i = 0; i < n; i++) pinv[perm[i] - 1] = i;

    /* Row magnitude accumulators for A^{-1}. */
    mpfr_t* row_sums = lum_array_alloc((size_t)n, bits);
    for (int i = 0; i < n; i++) mpfr_set_zero(row_sums[i], 1);

    mpfr_t* col_re = lum_array_alloc((size_t)n, bits);
    mpfr_t* col_im = is_complex ? lum_array_alloc((size_t)n, bits) : NULL;
    mpfr_t mag;
    mpfr_init2(mag, bits);

    for (int k = 0; k < n; k++) {
        lum_solve_LU_canonical(LU_re, LU_im, n, is_complex, pinv[k],
                                col_re, col_im, bits);
        for (int i = 0; i < n; i++) {
            if (is_complex) mpfr_hypot(mag, col_re[i], col_im[i], MPFR_RNDN);
            else            mpfr_abs(mag,   col_re[i], MPFR_RNDN);
            mpfr_add(row_sums[i], row_sums[i], mag, MPFR_RNDN);
        }
    }

    mpfr_set_zero(out, 1);
    for (int i = 0; i < n; i++) {
        if (mpfr_cmp(row_sums[i], out) > 0)
            mpfr_set(out, row_sums[i], MPFR_RNDN);
    }

    lum_array_free(row_sums, (size_t)n);
    lum_array_free(col_re,   (size_t)n);
    lum_array_free(col_im,   is_complex ? (size_t)n : 0);
    mpfr_clear(mag);
    free(pinv);
    return true;
}

/* ------------------------------------------------------------------ *
 *  Output builders.                                                     *
 * ------------------------------------------------------------------ */
static Expr* lum_make_scalar(const mpfr_t re, const mpfr_t im,
                              bool is_complex)
{
    if (!is_complex || mpfr_zero_p(im)) return expr_new_mpfr_copy(re);
    Expr** args = (Expr**)malloc(sizeof(Expr*) * 2);
    args[0] = expr_new_mpfr_copy(re);
    args[1] = expr_new_mpfr_copy(im);
    Expr* z = expr_new_function(expr_new_symbol("Complex"), args, 2);
    free(args);
    return z;
}

static Expr* lum_build_lu(const mpfr_t* LU_re, const mpfr_t* LU_im,
                           int n, bool is_complex, mpfr_prec_t bits)
{
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    mpfr_t zero_im;
    mpfr_init2(zero_im, bits);
    mpfr_set_zero(zero_im, 1);
    for (int i = 0; i < n; i++) {
        Expr** elems = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
        for (int j = 0; j < n; j++) {
            size_t off = (size_t)i * n + j;
            if (is_complex) {
                elems[j] = lum_make_scalar(LU_re[off], LU_im[off], true);
            } else {
                elems[j] = lum_make_scalar(LU_re[off], zero_im, false);
            }
        }
        rows[i] = expr_new_function(
            expr_new_symbol("List"), elems, (size_t)n);
        free(elems);
    }
    mpfr_clear(zero_im);
    Expr* lu = expr_new_function(
        expr_new_symbol("List"), rows, (size_t)n);
    free(rows);
    return lu;
}

static Expr* lum_build_perm(const int* perm, int n)
{
    Expr** elems = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    for (int k = 0; k < n; k++) elems[k] = expr_new_integer(perm[k]);
    Expr* p = expr_new_function(
        expr_new_symbol("List"), elems, (size_t)n);
    free(elems);
    return p;
}

/* ------------------------------------------------------------------ *
 *  Kernel.                                                              *
 * ------------------------------------------------------------------ */
Expr* lu_mpfr_dispatch(Expr* m, int n)
{
    static uint64_t sing_warn_counter = 0;

    CommonInexactInfo info = common_scan_inexact(m);
    if (!info.has_inexact) return NULL;
    mpfr_prec_t bits = (mpfr_prec_t)info.min_bits;
    if (bits < 53) bits = 53;

    mpfr_t* A_re = NULL;
    mpfr_t* A_im = NULL;
    bool is_complex = false;
    if (!lum_load_matrix(m, n, bits, &A_re, &A_im, &is_complex)) {
        return NULL;
    }

    /* Snapshot ||A||_inf before the factorisation overwrites A. */
    mpfr_t anorm;
    mpfr_init2(anorm, bits);
    lum_norm_inf(A_re, A_im, n, is_complex, bits, anorm);

    int* perm = (int*)malloc((size_t)n * sizeof(int));
    for (int j = 0; j < n; j++) perm[j] = j + 1;

    bool singular = false;
    bool ok = lum_factor(A_re, A_im, n, bits, is_complex, perm, &singular);
    if (!ok) {
        mpfr_clear(anorm);
        lum_array_free(A_re, (size_t)n * n);
        lum_array_free(A_im, is_complex ? (size_t)n * n : 0);
        free(perm);
        return NULL;
    }

    if (singular) {
        if (!sing_warn_counter) {
            sing_warn_counter = 1;
            char* s = expr_to_string(m);
            fprintf(stderr,
                "LUDecomposition::sing: Matrix %s is singular.\n", s);
            free(s);
        }
    }

    /* Condition number = ||A||_inf * ||A^{-1}||_inf. */
    mpfr_t inv_norm, cond;
    mpfr_init2(inv_norm, bits);
    mpfr_init2(cond,     bits);
    bool nonsing = lum_inv_norm_inf(A_re, A_im, perm, n, is_complex,
                                     bits, inv_norm);
    if (nonsing) mpfr_mul(cond, anorm, inv_norm, MPFR_RNDN);
    else         mpfr_set_inf(cond, 1);

    Expr* lu = lum_build_lu(A_re, A_im, n, is_complex, bits);
    Expr* p  = lum_build_perm(perm, n);
    Expr* c  = expr_new_mpfr_copy(cond);

    lum_array_free(A_re, (size_t)n * n);
    lum_array_free(A_im, is_complex ? (size_t)n * n : 0);
    free(perm);
    mpfr_clear(anorm);
    mpfr_clear(inv_norm);
    mpfr_clear(cond);

    Expr** items = (Expr**)malloc(sizeof(Expr*) * 3);
    items[0] = lu; items[1] = p; items[2] = c;
    Expr* result = expr_new_function(expr_new_symbol("List"), items, 3);
    free(items);
    return result;
}

#endif /* USE_MPFR */
