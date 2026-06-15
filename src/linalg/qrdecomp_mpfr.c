/* src/linalg/qrdecomp_mpfr.c
 *
 * Phase-4 arbitrary-precision QR kernel: Householder reflections over
 * MPFR arrays.
 *
 * Invoked by qr_dispatch (qrdecomp.c) when common_scan_inexact reports
 * an inexact input at min_bits > 53 (i.e. the matrix carries at least
 * one MPFR leaf above IEEE double precision).  On any failure --
 * USE_MPFR off, a leaf the loader can't reduce to an MPFR value,
 * rank-deficient input without pivoting -- the kernel returns NULL and
 * qr_dispatch falls back to the symbolic Modified-Gram-Schmidt
 * pipeline (which is precision-agnostic via rationalise / exact-MGS /
 * numericalise).
 *
 * Algorithm.  Standard Householder QR, exactly as a textbook would
 * lay it out, with column pivoting (Businger-Golub) when Pivoting ->
 * True is requested.  For each reflector step k = 0 .. min(n, p) - 1:
 *
 *     1. (if pivoting) pick the column in [k, p) with the largest
 *        residual squared norm and swap it into position k.
 *     2. Compute
 *           norm    = || A[k:n, k] ||_2                     (real)
 *           |x_0|   = | A[k, k] |
 *           beta    = -(x_0 / |x_0|) * norm                 (real or complex)
 *           v[0]    = x_0 - beta
 *           v[i]    = A[k + i, k]                           i = 1 .. n-k-1
 *           tau     = 1 / (norm * (norm + |x_0|))
 *        (When x_0 == 0 we pick beta = norm and tau = 1 / norm^2,
 *        which is the limit of the formula above.)
 *     3. Apply H = I - tau v v^H to A[k:n, k:p] from the left.
 *        After this A[k, k] = beta and A[k+1:n, k] = 0 in exact
 *        arithmetic; in MPFR they are within rounding.
 *     4. Apply H to Q[:, k:n] from the right so Q := Q * H accumulates
 *        the orthonormal factor.
 *     5. Record R[k, k] = beta and R[k, c] = A[k, c] for c > k.
 *
 * Rank detection.  After the loop we count R diagonals exceeding
 * 2^(-bits/2) * |R[0, 0]|.  With pivoting the diagonals are
 * non-increasing in magnitude, so the first below-threshold diagonal
 * terminates the rank count.  Without pivoting we still apply the
 * same cutoff after the fact and truncate the output -- but if any
 * reflector along the way fired against a zero-norm column (genuine
 * mid-stream rank deficiency that pivoting would have moved to the
 * end) we bail with NULL and let the symbolic dispatcher handle it.
 *
 * Output.  Q is column-major n*n; we emit its first `rank` columns
 * (conjugated-transposed) as the rank x n list `q`.  R is row-major
 * rank_cap*p; we emit its first `rank` rows as `r`.  When pivoting is
 * on we additionally emit the p x p permutation matrix P with
 * P[perm[j] - 1, j] = 1 (matching the LAPACK jpvt convention used by
 * the machine kernel).
 *
 * Memory contract.  Standard builtin contract: this file does NOT
 * call expr_free on the input `m` -- the evaluator owns it
 * (MEMORY.md / SPEC.md S4.1).  Every mpfr_t initialised here is
 * cleared along every exit path.
 */

#include "qrdecomp_internal.h"
#include "linalg.h"
#include "expr.h"
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

/* MPFR unavailable: short-circuit; qr_dispatch falls back to symbolic. */
Expr* qr_mpfr_dispatch(Expr* m, int n, int p, const QrOpts* opts)
{
    (void)m; (void)n; (void)p; (void)opts;
    return NULL;
}

#else /* USE_MPFR */

/* ------------------------------------------------------------------ *
 *  One-shot warning helper.  Mirrors matsol_warn_once /                *
 *  qr_machine_warn_once: a 64-bit counter passed by pointer suppresses *
 *  repeated reports of the same condition across a REPL session.       *
 * ------------------------------------------------------------------ */
static void qrm_warn_once(uint64_t* counter, const char* msg)
{
    if (*counter) return;
    *counter = 1;
    fprintf(stderr, "%s", msg);
}

/* ------------------------------------------------------------------ *
 *  Complex-detection probe.  Used during the first matrix scan so we   *
 *  know whether to allocate the im[] plane.  Performs an MPFR-level     *
 *  approximation at a token precision (53 bits is enough -- we only    *
 *  care whether the imaginary component is zero).                      *
 * ------------------------------------------------------------------ */
static bool qrm_leaf_is_complex(Expr* e)
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

/* ------------------------------------------------------------------ *
 *  MPFR array alloc / free helpers.  Duplicated here (rather than       *
 *  re-using eigen_common.c's mpfr_array_alloc) so this file is          *
 *  self-contained and does not pull in eigen_internal.h.  mpfr_init2    *
 *  leaves cells at NaN, so all callers must explicitly assign before    *
 *  reading.                                                             *
 * ------------------------------------------------------------------ */
static mpfr_t* qrm_array_alloc(size_t count, mpfr_prec_t bits)
{
    if (count == 0) return NULL;
    mpfr_t* a = (mpfr_t*)malloc(sizeof(mpfr_t) * count);
    for (size_t i = 0; i < count; i++) mpfr_init2(a[i], bits);
    return a;
}

static void qrm_array_free(mpfr_t* a, size_t count)
{
    if (!a) return;
    for (size_t i = 0; i < count; i++) mpfr_clear(a[i]);
    free(a);
}

/* ------------------------------------------------------------------ *
 *  Load the n x p input matrix into freshly-allocated column-major     *
 *  MPFR arrays at `bits` precision.                                     *
 *                                                                       *
 *  *out_im is non-NULL iff the matrix has at least one entry with a    *
 *  non-zero imaginary part.  Returns false (with nothing leaked and    *
 *  *out_re / *out_im undefined) on any cell that cannot be reduced     *
 *  to a numeric MPFR value.                                             *
 * ------------------------------------------------------------------ */
static bool qrm_load_matrix(Expr* m, int n, int p, mpfr_prec_t bits,
                            mpfr_t** out_re, mpfr_t** out_im,
                            bool* out_is_complex)
{
    if (m->type != EXPR_FUNCTION) return false;

    /* First pass: detect any complex entries.  We only need to allocate
     * the im[] plane when at least one cell carries a non-zero imag
     * component. */
    bool is_complex = false;
    for (int i = 0; i < n && !is_complex; i++) {
        Expr* row = m->data.function.args[i];
        if (!row || row->type != EXPR_FUNCTION) return false;
        for (int j = 0; j < p && !is_complex; j++) {
            if (qrm_leaf_is_complex(row->data.function.args[j])) {
                is_complex = true;
            }
        }
    }

    size_t total = (size_t)n * (size_t)p;
    mpfr_t* A_re = qrm_array_alloc(total, bits);
    mpfr_t* A_im = is_complex ? qrm_array_alloc(total, bits) : NULL;

    /* Second pass: column-major copy.  Cell (i, j) lives at i + j*n.
     * get_approx_mpfr writes both re and im even for purely real inputs
     * (im set to 0 in that case), so we always pass a tmp_im. */
    mpfr_t tmp_im;
    mpfr_init2(tmp_im, bits);
    for (int i = 0; i < n; i++) {
        Expr* row = m->data.function.args[i];
        for (int j = 0; j < p; j++) {
            size_t off = (size_t)i + (size_t)j * (size_t)n;
            bool is_inexact = false;
            if (!get_approx_mpfr(row->data.function.args[j],
                                 A_re[off], tmp_im, &is_inexact)) {
                mpfr_clear(tmp_im);
                qrm_array_free(A_re, total);
                qrm_array_free(A_im, total);
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

/* ------------------------------------------------------------------ *
 *  Initialise Q to the n x n identity (column-major).                  *
 * ------------------------------------------------------------------ */
static void qrm_init_identity(mpfr_t* Q_re, mpfr_t* Q_im, int n,
                              bool is_complex)
{
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            size_t off = (size_t)i + (size_t)j * (size_t)n;
            mpfr_set_si(Q_re[off], i == j ? 1 : 0, MPFR_RNDN);
            if (is_complex) mpfr_set_zero(Q_im[off], 1);
        }
    }
}

/* ------------------------------------------------------------------ *
 *  Squared norm of sub-column A[row_lo:n, col].                         *
 *      real:    sum_i      A_re[i + col*n]^2                           *
 *      complex: sum_i      A_re[..]^2 + A_im[..]^2                     *
 * ------------------------------------------------------------------ */
static void qrm_col_norm_sq(const mpfr_t* A_re, const mpfr_t* A_im,
                            int n, int col, int row_lo,
                            bool is_complex, mpfr_t out)
{
    mpfr_set_zero(out, 1);
    mpfr_t tmp;
    mpfr_init2(tmp, mpfr_get_prec(out));
    for (int i = row_lo; i < n; i++) {
        size_t off = (size_t)i + (size_t)col * (size_t)n;
        mpfr_sqr(tmp, A_re[off], MPFR_RNDN);
        mpfr_add(out, out, tmp, MPFR_RNDN);
        if (is_complex) {
            mpfr_sqr(tmp, A_im[off], MPFR_RNDN);
            mpfr_add(out, out, tmp, MPFR_RNDN);
        }
    }
    mpfr_clear(tmp);
}

/* Swap two column-major columns of an n x ? matrix. */
static void qrm_swap_columns(mpfr_t* A_re, mpfr_t* A_im, int n,
                             int col_a, int col_b, bool is_complex)
{
    if (col_a == col_b) return;
    for (int i = 0; i < n; i++) {
        size_t off_a = (size_t)i + (size_t)col_a * (size_t)n;
        size_t off_b = (size_t)i + (size_t)col_b * (size_t)n;
        mpfr_swap(A_re[off_a], A_re[off_b]);
        if (is_complex) mpfr_swap(A_im[off_a], A_im[off_b]);
    }
}

/* ------------------------------------------------------------------ *
 *  Householder QR core.                                                 *
 *                                                                       *
 *  Inputs:                                                              *
 *    A_re, A_im     column-major n*p, OVERWRITTEN during the run.       *
 *    n, p           dimensions.                                         *
 *    bits           working precision.                                  *
 *    is_complex     true iff A_im is allocated.                          *
 *    with_pivoting  do column pivoting (Businger-Golub).                *
 *                                                                       *
 *  Outputs:                                                             *
 *    Q_re, Q_im     column-major n*n; caller pre-allocates and inits    *
 *                   to identity.  Updated to Q = H_0 H_1 ... H_{k-1}.   *
 *    R_re, R_im     row-major rank_cap*p; caller pre-allocates and       *
 *                   zeroes.  Filled row by row as steps complete.       *
 *    perm           length p, 1-indexed (matches LAPACK jpvt).          *
 *                   Caller initialises to [1, 2, ..., p].                *
 *    rank_out       numerical rank (<= min(n, p)).                       *
 *                                                                       *
 *  Returns true on success.  Returns false on rank deficiency without   *
 *  pivoting (the caller falls back to the symbolic dispatcher, which    *
 *  handles non-pivoted rank-deficient input correctly).                  *
 * ------------------------------------------------------------------ */
static bool qrm_householder(
    mpfr_t* A_re, mpfr_t* A_im,
    int n, int p, mpfr_prec_t bits,
    bool is_complex, bool with_pivoting,
    mpfr_t* Q_re, mpfr_t* Q_im,
    mpfr_t* R_re, mpfr_t* R_im,
    int* perm, int* rank_out)
{
    int rank_cap = (n < p) ? n : p;
    *rank_out = 0;
    if (rank_cap == 0) return true;

    /* Scratch arrays.  v lives in n-sized space (technically max n-step
     * is needed but allocating n covers every step), w in p-sized, z in
     * n-sized.  Pre-allocating once outside the loop saves k allocations. */
    mpfr_t* v_re = qrm_array_alloc((size_t)n, bits);
    mpfr_t* v_im = is_complex ? qrm_array_alloc((size_t)n, bits) : NULL;
    mpfr_t* w_re = qrm_array_alloc((size_t)p, bits);
    mpfr_t* w_im = is_complex ? qrm_array_alloc((size_t)p, bits) : NULL;
    mpfr_t* z_re = qrm_array_alloc((size_t)n, bits);
    mpfr_t* z_im = is_complex ? qrm_array_alloc((size_t)n, bits) : NULL;

    /* Pivot cache: residual squared norms for every column. */
    mpfr_t* col_norms_sq = NULL;
    if (with_pivoting) {
        col_norms_sq = qrm_array_alloc((size_t)p, bits);
        for (int j = 0; j < p; j++) {
            qrm_col_norm_sq(A_re, A_im, n, j, 0, is_complex, col_norms_sq[j]);
        }
    }

    mpfr_t norm, norm_sq, abs_x0, beta_re, beta_im, tau;
    mpfr_t t1, threshold, max_diag, th_sq;
    mpfr_init2(norm,      bits);
    mpfr_init2(norm_sq,   bits);
    mpfr_init2(abs_x0,    bits);
    mpfr_init2(beta_re,   bits);
    mpfr_init2(beta_im,   bits);
    mpfr_init2(tau,       bits);
    mpfr_init2(t1,        bits);
    mpfr_init2(threshold, bits);
    mpfr_init2(max_diag,  bits);
    mpfr_init2(th_sq,     bits);
    mpfr_set_zero(max_diag,  1);
    mpfr_set_zero(threshold, 1);
    mpfr_set_zero(th_sq,     1);

    bool ok = true;
    for (int step = 0; step < rank_cap; step++) {

        /* ----- pivot selection ----- */
        if (with_pivoting) {
            int best = step;
            for (int j = step + 1; j < p; j++) {
                if (mpfr_cmp(col_norms_sq[j], col_norms_sq[best]) > 0) {
                    best = j;
                }
            }
            if (best != step) {
                qrm_swap_columns(A_re, A_im, n, step, best, is_complex);
                /* Also swap the entries in already-stored R rows.
                 * LAPACK's dgeqp3 gets this for free (R lives in the
                 * upper triangle of A, so swapping A columns swaps R
                 * columns too); since we keep R in a separate buffer,
                 * we have to mirror the rearrangement manually. */
                for (int i = 0; i < step; i++) {
                    size_t off_a = (size_t)i * (size_t)p + (size_t)step;
                    size_t off_b = (size_t)i * (size_t)p + (size_t)best;
                    mpfr_swap(R_re[off_a], R_re[off_b]);
                    if (is_complex) mpfr_swap(R_im[off_a], R_im[off_b]);
                }
                int tmp_pi = perm[step];
                perm[step] = perm[best];
                perm[best] = tmp_pi;
                mpfr_swap(col_norms_sq[step], col_norms_sq[best]);
            }

            /* Rank-deficient stop: pivot residual norm^2 below threshold.
             * threshold and th_sq are zero until the first beta has been
             * recorded (i.e. step > 0), so this is dormant on step 0. */
            if (step > 0 && mpfr_cmp(col_norms_sq[step], th_sq) <= 0) {
                break;
            }
        }

        /* ----- reflector ----- */
        qrm_col_norm_sq(A_re, A_im, n, step, step, is_complex, norm_sq);
        mpfr_sqrt(norm, norm_sq, MPFR_RNDN);

        /* Zero residual at step k: with pivoting we already screened
         * (so this would mean true rank-deficiency surfacing); without
         * pivoting a zero norm is the rank-deficient mid-stream case
         * that the symbolic kernel handles cleanly. */
        if (mpfr_zero_p(norm)) {
            if (with_pivoting) break;
            ok = false;
            goto cleanup;
        }

        size_t off_diag = (size_t)step + (size_t)step * (size_t)n;
        if (is_complex) {
            mpfr_hypot(abs_x0, A_re[off_diag], A_im[off_diag], MPFR_RNDN);
        } else {
            mpfr_abs(abs_x0, A_re[off_diag], MPFR_RNDN);
        }

        /* beta = -(x_0 / |x_0|) * norm.  When x_0 == 0 the phase is
         * arbitrary; pick beta = +norm (real). */
        if (mpfr_zero_p(abs_x0)) {
            mpfr_set(beta_re, norm, MPFR_RNDN);
            if (is_complex) mpfr_set_zero(beta_im, 1);
        } else {
            mpfr_div(t1, A_re[off_diag], abs_x0, MPFR_RNDN);
            mpfr_mul(beta_re, t1, norm, MPFR_RNDN);
            mpfr_neg(beta_re, beta_re, MPFR_RNDN);
            if (is_complex) {
                mpfr_div(t1, A_im[off_diag], abs_x0, MPFR_RNDN);
                mpfr_mul(beta_im, t1, norm, MPFR_RNDN);
                mpfr_neg(beta_im, beta_im, MPFR_RNDN);
            }
        }

        /* tau = 1 / (norm * (norm + |x_0|)) -- real positive. */
        mpfr_add(t1, norm, abs_x0, MPFR_RNDN);
        mpfr_mul(t1, t1, norm, MPFR_RNDN);
        mpfr_ui_div(tau, 1, t1, MPFR_RNDN);

        /* v[0] = x_0 - beta;  v[i] = x_i for i > 0.  v lives in
         * positions 0..n-step-1 (length n - step). */
        int m_eff = n - step;
        mpfr_sub(v_re[0], A_re[off_diag], beta_re, MPFR_RNDN);
        if (is_complex) mpfr_sub(v_im[0], A_im[off_diag], beta_im, MPFR_RNDN);
        for (int i = 1; i < m_eff; i++) {
            size_t off = (size_t)(step + i) + (size_t)step * (size_t)n;
            mpfr_set(v_re[i], A_re[off], MPFR_RNDN);
            if (is_complex) mpfr_set(v_im[i], A_im[off], MPFR_RNDN);
        }

        /* ----- apply H = I - tau v v^H to A[step:n, step:p] -----
         *  For each column c in [step, p):
         *      w_c = sum_i conj(v[i]) * A[step+i, c]                (scalar)
         *      A[step+i, c] -= tau * v[i] * w_c                     (vector)
         *  Real case collapses the conj() to identity. */
        for (int c = step; c < p; c++) {
            int wi = c - step;
            mpfr_set_zero(w_re[wi], 1);
            if (is_complex) mpfr_set_zero(w_im[wi], 1);
            for (int i = 0; i < m_eff; i++) {
                size_t off = (size_t)(step + i) + (size_t)c * (size_t)n;
                /* w_re += v_re * A_re  (always) */
                mpfr_mul(t1, v_re[i], A_re[off], MPFR_RNDN);
                mpfr_add(w_re[wi], w_re[wi], t1, MPFR_RNDN);
                if (is_complex) {
                    /* (v_re - I v_im)(A_re + I A_im)
                     *  = (v_re A_re + v_im A_im) + I (v_re A_im - v_im A_re) */
                    mpfr_mul(t1, v_im[i], A_im[off], MPFR_RNDN);
                    mpfr_add(w_re[wi], w_re[wi], t1, MPFR_RNDN);
                    mpfr_mul(t1, v_re[i], A_im[off], MPFR_RNDN);
                    mpfr_add(w_im[wi], w_im[wi], t1, MPFR_RNDN);
                    mpfr_mul(t1, v_im[i], A_re[off], MPFR_RNDN);
                    mpfr_sub(w_im[wi], w_im[wi], t1, MPFR_RNDN);
                }
            }
            /* w *= tau */
            mpfr_mul(w_re[wi], w_re[wi], tau, MPFR_RNDN);
            if (is_complex) mpfr_mul(w_im[wi], w_im[wi], tau, MPFR_RNDN);

            for (int i = 0; i < m_eff; i++) {
                size_t off = (size_t)(step + i) + (size_t)c * (size_t)n;
                /* A_re -= v_re * w_re  (always) */
                mpfr_mul(t1, v_re[i], w_re[wi], MPFR_RNDN);
                mpfr_sub(A_re[off], A_re[off], t1, MPFR_RNDN);
                if (is_complex) {
                    /* (v_re + I v_im)(w_re + I w_im)
                     *  = (v_re w_re - v_im w_im) + I (v_re w_im + v_im w_re) */
                    mpfr_mul(t1, v_im[i], w_im[wi], MPFR_RNDN);
                    mpfr_add(A_re[off], A_re[off], t1, MPFR_RNDN);
                    mpfr_mul(t1, v_re[i], w_im[wi], MPFR_RNDN);
                    mpfr_sub(A_im[off], A_im[off], t1, MPFR_RNDN);
                    mpfr_mul(t1, v_im[i], w_re[wi], MPFR_RNDN);
                    mpfr_sub(A_im[off], A_im[off], t1, MPFR_RNDN);
                }
            }
        }

        /* ----- update cached column norms for the remaining
         *       columns (Businger-Golub eq. 4) -----
         *  After the reflector A[step, c] is the "popped-off" row of R;
         *  the residual norm of column c in [step+1, n) is reduced by
         *  |A[step, c]|^2 from its pre-reflector value. */
        if (with_pivoting) {
            for (int c = step + 1; c < p; c++) {
                size_t off = (size_t)step + (size_t)c * (size_t)n;
                mpfr_sqr(t1, A_re[off], MPFR_RNDN);
                mpfr_sub(col_norms_sq[c], col_norms_sq[c], t1, MPFR_RNDN);
                if (is_complex) {
                    mpfr_sqr(t1, A_im[off], MPFR_RNDN);
                    mpfr_sub(col_norms_sq[c], col_norms_sq[c], t1, MPFR_RNDN);
                }
                /* Floor at zero to absorb subtractive rounding noise. */
                if (mpfr_sgn(col_norms_sq[c]) < 0)
                    mpfr_set_zero(col_norms_sq[c], 1);
            }
        }

        /* ----- apply H to Q[:, step:n] from the right ------------------
         *  Q := Q * H = Q - tau (Q * v) v^H.
         *  For each row r:
         *    z_r = sum_i Q[r, step+i] * v[i]
         *    Q[r, step+i] -= tau * z_r * conj(v[i]) */
        for (int r = 0; r < n; r++) {
            mpfr_set_zero(z_re[r], 1);
            if (is_complex) mpfr_set_zero(z_im[r], 1);
        }
        for (int i = 0; i < m_eff; i++) {
            for (int r = 0; r < n; r++) {
                size_t off = (size_t)r + (size_t)(step + i) * (size_t)n;
                mpfr_mul(t1, Q_re[off], v_re[i], MPFR_RNDN);
                mpfr_add(z_re[r], z_re[r], t1, MPFR_RNDN);
                if (is_complex) {
                    /* (Q_re + I Q_im)(v_re + I v_im)
                     *  = (Q_re v_re - Q_im v_im) + I (Q_re v_im + Q_im v_re) */
                    mpfr_mul(t1, Q_im[off], v_im[i], MPFR_RNDN);
                    mpfr_sub(z_re[r], z_re[r], t1, MPFR_RNDN);
                    mpfr_mul(t1, Q_re[off], v_im[i], MPFR_RNDN);
                    mpfr_add(z_im[r], z_im[r], t1, MPFR_RNDN);
                    mpfr_mul(t1, Q_im[off], v_re[i], MPFR_RNDN);
                    mpfr_add(z_im[r], z_im[r], t1, MPFR_RNDN);
                }
            }
        }
        for (int r = 0; r < n; r++) {
            mpfr_mul(z_re[r], z_re[r], tau, MPFR_RNDN);
            if (is_complex) mpfr_mul(z_im[r], z_im[r], tau, MPFR_RNDN);
        }
        for (int i = 0; i < m_eff; i++) {
            for (int r = 0; r < n; r++) {
                size_t off = (size_t)r + (size_t)(step + i) * (size_t)n;
                mpfr_mul(t1, z_re[r], v_re[i], MPFR_RNDN);
                mpfr_sub(Q_re[off], Q_re[off], t1, MPFR_RNDN);
                if (is_complex) {
                    /* (z_re + I z_im)(v_re - I v_im)
                     *  = (z_re v_re + z_im v_im) + I (z_im v_re - z_re v_im) */
                    mpfr_mul(t1, z_im[r], v_im[i], MPFR_RNDN);
                    mpfr_sub(Q_re[off], Q_re[off], t1, MPFR_RNDN);
                    mpfr_mul(t1, z_im[r], v_re[i], MPFR_RNDN);
                    mpfr_sub(Q_im[off], Q_im[off], t1, MPFR_RNDN);
                    mpfr_mul(t1, z_re[r], v_im[i], MPFR_RNDN);
                    mpfr_add(Q_im[off], Q_im[off], t1, MPFR_RNDN);
                }
            }
        }

        /* ----- record R[step, step:p] = A[step, step:p] -----
         *  A[step, step] is now beta in exact arithmetic; we set it
         *  explicitly from the stored beta_re/im to keep the row clean
         *  of accumulated rounding noise.  R[step, c] for c > step
         *  comes from A (i.e. (Q^H · column c) at the step row). */
        size_t r_diag = (size_t)step * (size_t)p + (size_t)step;
        mpfr_set(R_re[r_diag], beta_re, MPFR_RNDN);
        if (is_complex) mpfr_set(R_im[r_diag], beta_im, MPFR_RNDN);
        for (int c = step + 1; c < p; c++) {
            size_t src = (size_t)step + (size_t)c * (size_t)n;
            size_t dst = (size_t)step * (size_t)p + (size_t)c;
            mpfr_set(R_re[dst], A_re[src], MPFR_RNDN);
            if (is_complex) mpfr_set(R_im[dst], A_im[src], MPFR_RNDN);
        }

        /* Track |beta| and update the rank-revealing threshold. */
        if (is_complex) mpfr_hypot(t1, beta_re, beta_im, MPFR_RNDN);
        else            mpfr_abs(t1, beta_re, MPFR_RNDN);
        if (mpfr_cmp(t1, max_diag) > 0) mpfr_set(max_diag, t1, MPFR_RNDN);

        /* threshold = max_diag * 2^(-bits/2);  th_sq = threshold^2. */
        mpfr_set(threshold, max_diag, MPFR_RNDN);
        mpfr_div_2si(threshold, threshold, (long)(bits / 2), MPFR_RNDN);
        mpfr_sqr(th_sq, threshold, MPFR_RNDN);

        (*rank_out) = step + 1;
    }

cleanup:
    /* Numerical rank refinement: count R diagonals above the threshold,
     * starting from i = 0 and stopping at the first below-threshold
     * entry.  With pivoting the R diagonals are non-increasing so this
     * is exact; without pivoting it's a best-effort guard. */
    if (ok && *rank_out > 0) {
        int true_rank = 0;
        mpfr_t mag;
        mpfr_init2(mag, bits);
        for (int i = 0; i < *rank_out; i++) {
            size_t off = (size_t)i * (size_t)p + (size_t)i;
            if (is_complex) mpfr_hypot(mag, R_re[off], R_im[off], MPFR_RNDN);
            else            mpfr_abs(mag,   R_re[off], MPFR_RNDN);
            if (mpfr_cmp(mag, threshold) > 0) true_rank = i + 1;
            else break;
        }
        *rank_out = true_rank;
        mpfr_clear(mag);
    }

    qrm_array_free(v_re, (size_t)n);
    qrm_array_free(v_im, is_complex ? (size_t)n : 0);
    qrm_array_free(w_re, (size_t)p);
    qrm_array_free(w_im, is_complex ? (size_t)p : 0);
    qrm_array_free(z_re, (size_t)n);
    qrm_array_free(z_im, is_complex ? (size_t)n : 0);
    if (col_norms_sq) qrm_array_free(col_norms_sq, (size_t)p);

    mpfr_clear(norm);
    mpfr_clear(norm_sq);
    mpfr_clear(abs_x0);
    mpfr_clear(beta_re);
    mpfr_clear(beta_im);
    mpfr_clear(tau);
    mpfr_clear(t1);
    mpfr_clear(threshold);
    mpfr_clear(max_diag);
    mpfr_clear(th_sq);
    return ok;
}

/* ------------------------------------------------------------------ *
 *  Output: build a scalar Expr from an (re, im) MPFR pair.              *
 *                                                                       *
 *  When is_complex is false (or im is identically zero) we emit a bare  *
 *  EXPR_MPFR; otherwise Complex[EXPR_MPFR, EXPR_MPFR].  The collapse    *
 *  to-real-when-im-is-zero matches the machine kernel.                  *
 * ------------------------------------------------------------------ */
static Expr* qrm_make_scalar(const mpfr_t re, const mpfr_t im,
                             bool is_complex)
{
    if (!is_complex || mpfr_zero_p(im)) return expr_new_mpfr_copy(re);
    Expr** args = (Expr**)malloc(sizeof(Expr*) * 2);
    args[0] = expr_new_mpfr_copy(re);
    args[1] = expr_new_mpfr_copy(im);
    Expr* z = expr_new_function(expr_new_symbol(SYM_Complex), args, 2);
    free(args);
    return z;
}

/* Build the q output (rank x n) from the column-major Q buffer (n x n).
 *   Real input:    q[j, i] = Q[i, j]            (Transpose)
 *   Complex input: q[j, i] = Conjugate(Q[i, j]) (ConjugateTranspose)
 *
 * Convention matches the symbolic and machine kernels -- preserves the
 * public identity m == ConjugateTranspose[q] . r for both real and
 * complex matrices without ever emitting Conjugate[...] heads. */
static Expr* qrm_build_q(const mpfr_t* Q_re, const mpfr_t* Q_im,
                         int n, int rank, bool is_complex,
                         mpfr_prec_t bits)
{
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)rank);
    mpfr_t neg_im;
    mpfr_init2(neg_im, bits);
    mpfr_set_zero(neg_im, 1);
    for (int j = 0; j < rank; j++) {
        Expr** elems = NULL;
        if (n > 0) elems = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
        for (int i = 0; i < n; i++) {
            size_t off = (size_t)i + (size_t)j * (size_t)n;
            if (is_complex) {
                mpfr_neg(neg_im, Q_im[off], MPFR_RNDN);
                elems[i] = qrm_make_scalar(Q_re[off], neg_im, true);
            } else {
                elems[i] = qrm_make_scalar(Q_re[off], neg_im, false);
            }
        }
        rows[j] = expr_new_function(expr_new_symbol(SYM_List),
                                    elems, (size_t)n);
        if (elems) free(elems);
    }
    mpfr_clear(neg_im);
    Expr* q = expr_new_function(expr_new_symbol(SYM_List),
                                rows, (size_t)rank);
    free(rows);
    return q;
}

/* Build the r output (rank x p) from the row-major R buffer
 * (rank_cap x p, only the first `rank` rows are read).  Entries
 * strictly below the leading echelon are guaranteed zero by qrm_householder
 * (R_im / R_re were zeroed by the caller and only the upper trapezoid
 * was ever written). */
static Expr* qrm_build_r(const mpfr_t* R_re, const mpfr_t* R_im,
                         int rank, int p, bool is_complex,
                         mpfr_prec_t bits)
{
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)rank);
    mpfr_t zero_im;
    mpfr_init2(zero_im, bits);
    mpfr_set_zero(zero_im, 1);
    for (int j = 0; j < rank; j++) {
        Expr** elems = NULL;
        if (p > 0) elems = (Expr**)malloc(sizeof(Expr*) * (size_t)p);
        for (int k = 0; k < p; k++) {
            size_t off = (size_t)j * (size_t)p + (size_t)k;
            if (is_complex) {
                elems[k] = qrm_make_scalar(R_re[off], R_im[off], true);
            } else {
                elems[k] = qrm_make_scalar(R_re[off], zero_im, false);
            }
        }
        rows[j] = expr_new_function(expr_new_symbol(SYM_List),
                                    elems, (size_t)p);
        if (elems) free(elems);
    }
    mpfr_clear(zero_im);
    Expr* r = expr_new_function(expr_new_symbol(SYM_List),
                                rows, (size_t)rank);
    free(rows);
    return r;
}

/* Build a p x p Integer permutation matrix from a 1-indexed perm array.
 *   P[perm[j] - 1, j] = 1 -- exactly matches the machine kernel's
 *   mach_build_perm convention so unit tests can compare against the
 *   same Integer 0 / 1 entries. */
static Expr* qrm_build_perm(const int* perm, int p)
{
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)p);
    for (int i = 0; i < p; i++) {
        Expr** elems = NULL;
        if (p > 0) elems = (Expr**)malloc(sizeof(Expr*) * (size_t)p);
        for (int j = 0; j < p; j++) {
            int src_col = perm[j] - 1;
            elems[j] = expr_new_integer(src_col == i ? 1 : 0);
        }
        rows[i] = expr_new_function(expr_new_symbol(SYM_List),
                                    elems, (size_t)p);
        if (elems) free(elems);
    }
    Expr* P = expr_new_function(expr_new_symbol(SYM_List),
                                rows, (size_t)p);
    free(rows);
    return P;
}

/* ------------------------------------------------------------------ *
 *  The kernel.                                                          *
 *                                                                       *
 *  Lives behind the qr_mpfr_dispatch entry point declared in            *
 *  qrdecomp_internal.h.  Returns NULL on any failure path; the caller   *
 *  (qr_dispatch) treats NULL as "fall back to symbolic" and never       *
 *  frees the input.                                                     *
 * ------------------------------------------------------------------ */
Expr* qr_mpfr_dispatch(Expr* m, int n, int p, const QrOpts* opts)
{
    static uint64_t mpfr_warn_counter = 0;

    /* Working precision = min input precision (matches the
     * inexact-in / inexact-out convention used by the symbolic and
     * machine kernels). */
    CommonInexactInfo info = common_scan_inexact(m);
    if (!info.has_inexact) return NULL;
    mpfr_prec_t bits = (mpfr_prec_t)info.min_bits;
    if (bits < 53) bits = 53;

    /* Load matrix.  Non-numeric leaves bail to symbolic. */
    mpfr_t* A_re = NULL;
    mpfr_t* A_im = NULL;
    bool is_complex = false;
    if (!qrm_load_matrix(m, n, p, bits, &A_re, &A_im, &is_complex)) {
        return NULL;
    }

    int rank_cap = (n < p) ? n : p;
    size_t total_A = (size_t)n * (size_t)p;
    size_t total_Q = (size_t)n * (size_t)n;
    size_t total_R = (size_t)rank_cap * (size_t)p;

    /* Q identity, R zero. */
    mpfr_t* Q_re = qrm_array_alloc(total_Q, bits);
    mpfr_t* Q_im = is_complex ? qrm_array_alloc(total_Q, bits) : NULL;
    qrm_init_identity(Q_re, Q_im, n, is_complex);

    mpfr_t* R_re = qrm_array_alloc(total_R, bits);
    mpfr_t* R_im = is_complex ? qrm_array_alloc(total_R, bits) : NULL;
    for (size_t t = 0; t < total_R; t++) {
        mpfr_set_zero(R_re[t], 1);
        if (R_im) mpfr_set_zero(R_im[t], 1);
    }

    int* perm = (int*)malloc((size_t)p * sizeof(int));
    for (int j = 0; j < p; j++) perm[j] = j + 1;

    int rank = 0;
    bool ok = qrm_householder(A_re, A_im, n, p, bits, is_complex,
                              opts->pivoting,
                              Q_re, Q_im, R_re, R_im, perm, &rank);

    /* A is no longer needed (overwritten during reflectors). */
    qrm_array_free(A_re, total_A);
    qrm_array_free(A_im, is_complex ? total_A : 0);

    if (!ok) {
        qrm_warn_once(&mpfr_warn_counter,
            "QRDecomposition: MPFR fast path hit rank-deficient input "
            "without pivoting; falling back to symbolic kernel.\n");
        qrm_array_free(Q_re, total_Q);
        qrm_array_free(Q_im, is_complex ? total_Q : 0);
        qrm_array_free(R_re, total_R);
        qrm_array_free(R_im, is_complex ? total_R : 0);
        free(perm);
        return NULL;
    }

    /* Build outputs.  Empty-rank case -> empty `{}` matrices (same as
     * symbolic / machine kernels). */
    Expr* q;
    Expr* r;
    if (rank == 0) {
        q = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
        r = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    } else {
        q = qrm_build_q(Q_re, Q_im, n, rank, is_complex, bits);
        r = qrm_build_r(R_re, R_im, rank, p, is_complex, bits);
    }

    qrm_array_free(Q_re, total_Q);
    qrm_array_free(Q_im, is_complex ? total_Q : 0);
    qrm_array_free(R_re, total_R);
    qrm_array_free(R_im, is_complex ? total_R : 0);

    Expr* result;
    if (opts->pivoting) {
        Expr* P = qrm_build_perm(perm, p);
        free(perm);
        Expr** items = (Expr**)malloc(sizeof(Expr*) * 3);
        items[0] = q; items[1] = r; items[2] = P;
        result = expr_new_function(expr_new_symbol(SYM_List), items, 3);
        free(items);
    } else {
        free(perm);
        Expr** items = (Expr**)malloc(sizeof(Expr*) * 2);
        items[0] = q; items[1] = r;
        result = expr_new_function(expr_new_symbol(SYM_List), items, 2);
        free(items);
    }
    return result;
}

#endif /* USE_MPFR */
