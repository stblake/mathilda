/* svdecomp_mpfr.c
 *
 * Arbitrary-precision SingularValueDecomposition kernel.
 *
 * Algorithm: one-sided Jacobi SVD over column-major MPFR arrays at
 * the input's working precision.  One-sided Jacobi is the standard
 * choice for arbitrary-precision SVD: it converges quadratically, it
 * needs no auxiliary tridiagonal step, and it generalises cleanly to
 * complex matrices (paired re/im arrays, complex Jacobi rotations).
 *
 * Picked by svd_dispatch (svdecomp.c) when the input contains at
 * least one inexact leaf and the minimum precision is > 53 bits.
 *
 * Both real and complex inputs are supported.  Complex matrices use
 * paired (re, im) MPFR buffers throughout: A_re/A_im, V_re/V_im,
 * U_re/U_im are allocated as separate column-major arrays.  The 2x2
 * complex Jacobi rotation factors as a phase rotation that makes the
 * column inner product real-positive, followed by a real Jacobi
 * rotation in the same plane.  See svdm_jacobi_sweeps_complex.
 *
 * Limitations (return NULL so the call falls back to the symbolic
 * dispatcher, which numericalises and re-dispatches when possible):
 *
 *   - Generalized SVD ({m, a} input).  Lands in Phase E.
 *
 * High-level flow (real, n x p):
 *
 *   1. Load m into A_re (column-major n*p) at `bits` precision.
 *   2. Initialise V_re to the p x p identity.
 *   3. Sweep:  for each pair (i, j), 0 <= i < j < p, compute
 *        alpha = || A[:, i] ||^2,
 *        gamma = || A[:, j] ||^2,
 *        beta  = < A[:, i], A[:, j] >.
 *      If 2*|beta| > tol * (alpha + gamma) the columns aren't
 *      sufficiently orthogonal -- compute the Jacobi rotation (c, s)
 *      that diagonalises [alpha beta; beta gamma] and apply it to
 *      A's columns i, j and to V's columns i, j.
 *      Repeat sweeps until the off-diagonal-norm bound falls under
 *      tol * Frobenius(A) (typically 2-3 sweeps for well-conditioned
 *      matrices, more for poorly-scaled ones).
 *   4. Singular values: sigma_i = || A[:, i] || after convergence.
 *   5. Left singular vectors: U[:, i] = A[:, i] / sigma_i for sigma_i > 0;
 *      orthonormal-complete the remaining n - rank columns via plain
 *      Gram-Schmidt against the identity at MPFR precision.
 *   6. Sort sigma descending and permute U / V columns to match.
 *   7. Wrap U (n x n), Sigma (n x p), V (p x p) as Mathilda lists with
 *      EXPR_MPFR entries at the input precision.
 *   8. Hand to svd_apply_postprocess for tolerance, truncation, and
 *      TargetStructure.
 *
 * Memory contract.  Standard builtin contract (see svdecomp.h): this
 * file does NOT call expr_free on the input matrix.  Every mpfr_t
 * initialised here is cleared along every exit path.
 */

#include "svdecomp.h"
#include "svdecomp_internal.h"
#include "linalg.h"
#include "common.h"
#include "sym_names.h"
#include "numeric.h"
#include "message.h"   /* mth_message: Quiet/Check funnel */
#include "eval.h"
#include "mp_parallel.h"

#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#ifdef USE_FLINT
#include <flint/flint.h>
#include <flint/arb.h>
#endif
#endif

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef USE_MPFR

Expr* svd_mpfr_dispatch(const SvdArgs* args, int n, int p, int n_a)
{
    (void)args; (void)n; (void)p; (void)n_a;
    return NULL;
}

#else  /* USE_MPFR */

/* ---------------------------------------------------------------------
 * MPFR array helpers (duplicated from qrdecomp_mpfr.c for self-
 * containment -- same style and contract).  mpfr_init2 leaves cells at
 * NaN; every caller must explicitly assign before reading. */

static mpfr_t* svdm_array_alloc(size_t count, mpfr_prec_t bits)
{
    if (count == 0) return NULL;
    mpfr_t* a = (mpfr_t*)malloc(sizeof(mpfr_t) * count);
    for (size_t i = 0; i < count; i++) mpfr_init2(a[i], bits);
    return a;
}

static void svdm_array_free(mpfr_t* a, size_t count)
{
    if (!a) return;
    for (size_t i = 0; i < count; i++) mpfr_clear(a[i]);
    free(a);
}

/* True iff the matrix has at least one entry with a non-zero imaginary
 * component (probed at token 53-bit precision -- we only care about
 * zero-vs-non-zero, not magnitude). */
static bool svdm_leaf_is_complex(Expr* e)
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

/* Walk the matrix and return true iff any leaf has a non-zero
 * imaginary component (probed at 53 bits -- we only care about
 * zero-vs-non-zero, not magnitude).  Used by the MPFR dispatcher to
 * pick the real or complex Jacobi kernel. */
static bool svdm_matrix_has_complex(Expr* m, int n, int p)
{
    if (m->type != EXPR_FUNCTION) return false;
    for (int i = 0; i < n; i++) {
        Expr* row = m->data.function.args[i];
        if (!row || row->type != EXPR_FUNCTION) return false;
        for (int j = 0; j < p; j++) {
            if (svdm_leaf_is_complex(row->data.function.args[j])) return true;
        }
    }
    return false;
}

/* Load the n x p matrix into a column-major MPFR buffer.  Returns
 * false on any cell that doesn't reduce to a numeric MPFR value or if
 * the matrix carries imaginary content (the caller routes complex
 * input through svdm_load_complex instead). */
static bool svdm_load_real(Expr* m, int n, int p, mpfr_prec_t bits,
                           mpfr_t** out_A)
{
    if (m->type != EXPR_FUNCTION) return false;
    if (svdm_matrix_has_complex(m, n, p)) return false;

    size_t total = (size_t)n * (size_t)p;
    mpfr_t* A = svdm_array_alloc(total, bits);

    mpfr_t tmp_im;
    mpfr_init2(tmp_im, bits);
    for (int i = 0; i < n; i++) {
        Expr* row = m->data.function.args[i];
        for (int j = 0; j < p; j++) {
            size_t off = (size_t)i + (size_t)j * (size_t)n;
            bool is_inexact = false;
            if (!get_approx_mpfr(row->data.function.args[j],
                                 A[off], tmp_im, &is_inexact)) {
                mpfr_clear(tmp_im);
                svdm_array_free(A, total);
                return false;
            }
        }
    }
    mpfr_clear(tmp_im);

    *out_A = A;
    return true;
}

/* Load the n x p matrix into paired column-major (re, im) MPFR
 * buffers.  Both A_re and A_im are allocated as separate `total`-sized
 * arrays.  Returns false on any cell that doesn't reduce to a numeric
 * value (Integer / Rational / Real / Complex / MPFR).  On failure
 * nothing leaks: both buffers are released before return. */
static bool svdm_load_complex(Expr* m, int n, int p, mpfr_prec_t bits,
                              mpfr_t** out_re, mpfr_t** out_im)
{
    if (m->type != EXPR_FUNCTION) return false;

    size_t total = (size_t)n * (size_t)p;
    mpfr_t* A_re = svdm_array_alloc(total, bits);
    mpfr_t* A_im = svdm_array_alloc(total, bits);

    for (int i = 0; i < n; i++) {
        Expr* row = m->data.function.args[i];
        if (!row || row->type != EXPR_FUNCTION) {
            svdm_array_free(A_re, total);
            svdm_array_free(A_im, total);
            return false;
        }
        for (int j = 0; j < p; j++) {
            size_t off = (size_t)i + (size_t)j * (size_t)n;
            bool is_inexact = false;
            if (!get_approx_mpfr(row->data.function.args[j],
                                 A_re[off], A_im[off], &is_inexact)) {
                svdm_array_free(A_re, total);
                svdm_array_free(A_im, total);
                return false;
            }
        }
    }

    *out_re = A_re;
    *out_im = A_im;
    return true;
}

/* ---------------------------------------------------------------------
 * Jacobi-rotation primitive.
 *
 *   Given a 2x2 symmetric block [alpha beta; beta gamma], compute the
 *   rotation (c, s) such that
 *     [ c  s][alpha beta][c -s] = [lam1  0 ]
 *     [-s  c][beta gamma][s  c]   [ 0  lam2]
 *
 *   Standard formula:
 *     if beta == 0:   c = 1, s = 0
 *     else:
 *       tau = (gamma - alpha) / (2 * beta)
 *       t   = sign(tau) / (|tau| + sqrt(1 + tau^2))   (or 1/(2 tau) for large |tau|)
 *       c   = 1 / sqrt(1 + t^2)
 *       s   = t * c
 *
 * Apply to columns i and j of an n x ? matrix C (column-major):
 *     col_i' =  c * col_i + s * col_j
 *     col_j' = -s * col_i + c * col_j
 *
 * Used to update both A and V during each rotation. */
static void svdm_apply_rotation(mpfr_t* C, int n, int col_i, int col_j,
                                mpfr_t c, mpfr_t s, mpfr_prec_t bits)
{
    mpfr_t ai, aj, t1, t2;
    mpfr_init2(ai, bits); mpfr_init2(aj, bits);
    mpfr_init2(t1, bits); mpfr_init2(t2, bits);
    for (int r = 0; r < n; r++) {
        size_t oi = (size_t)r + (size_t)col_i * (size_t)n;
        size_t oj = (size_t)r + (size_t)col_j * (size_t)n;
        mpfr_set(ai, C[oi], MPFR_RNDN);
        mpfr_set(aj, C[oj], MPFR_RNDN);
        mpfr_mul(t1, c, ai, MPFR_RNDN);
        mpfr_mul(t2, s, aj, MPFR_RNDN);
        mpfr_add(C[oi], t1, t2, MPFR_RNDN);
        mpfr_mul(t1, s, ai, MPFR_RNDN);
        mpfr_mul(t2, c, aj, MPFR_RNDN);
        mpfr_sub(C[oj], t2, t1, MPFR_RNDN);
    }
    mpfr_clear(ai); mpfr_clear(aj); mpfr_clear(t1); mpfr_clear(t2);
}

/* Inner product <A[:, i], A[:, j]> for two columns of an n x ? matrix. */
static void svdm_inner(const mpfr_t* A, int n, int i, int j,
                       mpfr_t out, mpfr_prec_t bits)
{
    mpfr_set_zero(out, 1);
    mpfr_t tmp; mpfr_init2(tmp, bits);
    for (int r = 0; r < n; r++) {
        size_t oi = (size_t)r + (size_t)i * (size_t)n;
        size_t oj = (size_t)r + (size_t)j * (size_t)n;
        mpfr_mul(tmp, A[oi], A[oj], MPFR_RNDN);
        mpfr_add(out, out, tmp, MPFR_RNDN);
    }
    mpfr_clear(tmp);
}

/* ---------------------------------------------------------------------
 * One-sided Jacobi SVD.
 *
 * Inputs:
 *   A          n x p column-major MPFR matrix.  Overwritten with
 *              [sigma_1 u_1 | sigma_2 u_2 | ... | sigma_p u_p] on exit.
 *   V          p x p column-major MPFR matrix.  Caller initialises to
 *              the identity; updated to the accumulated right-side
 *              orthogonal transform.
 *   bits       working precision.
 *
 * Returns the number of sweeps performed.  Convergence threshold is
 * 2^(-bits/2): each sweep should reduce the off-diagonal norm
 * quadratically, so a handful of sweeps suffice in practice.
 * ------------------------------------------------------------------ */
static int svdm_jacobi_sweeps(mpfr_t* A, int n, int p, mpfr_t* V,
                              mpfr_prec_t bits)
{
    mpfr_t alpha, gamma, beta, tau, t, c, s, abs_t, tmp;
    mpfr_init2(alpha, bits); mpfr_init2(gamma, bits);
    mpfr_init2(beta,  bits); mpfr_init2(tau,   bits);
    mpfr_init2(t,     bits); mpfr_init2(c,     bits);
    mpfr_init2(s,     bits); mpfr_init2(abs_t, bits);
    mpfr_init2(tmp,   bits);

    /* Convergence tolerance n * 2^-bits of the column norms (the one-sided
     * Jacobi criterion of Demmel-Veselic, LAPACK's m * eps).  Stopping at
     * half the precision, 2^(-bits/2), as this used to, leaves the columns
     * orthogonal -- and so U and V accurate -- to half the digits the result
     * displays; quadratic convergence makes the full criterion cost only a
     * sweep or two more.  n * 2^-bits stays above the rounding floor of an
     * n-term inner product, so it cannot stall. */
    mpfr_t tol; mpfr_init2(tol, bits);
    mpfr_set_si(tol, n, MPFR_RNDN);
    mpfr_mul_2si(tol, tol, -(long)bits, MPFR_RNDN);

    int max_sweeps = 30;     /* safety upper bound; usually <= 5 needed */
    int sweep = 0;
    for (sweep = 0; sweep < max_sweeps; sweep++) {
        bool rotated = false;
        for (int i = 0; i < p - 1; i++) {
            for (int j = i + 1; j < p; j++) {
                svdm_inner(A, n, i, i, alpha, bits);
                svdm_inner(A, n, j, j, gamma, bits);
                svdm_inner(A, n, i, j, beta,  bits);
                /* off-diagonal test: 2|beta| > tol * sqrt(alpha * gamma) */
                mpfr_abs(tmp, beta, MPFR_RNDN);
                mpfr_mul_2si(tmp, tmp, 1, MPFR_RNDN);   /* 2 |beta| */
                mpfr_mul(t, alpha, gamma, MPFR_RNDN);
                mpfr_sqrt(t, t, MPFR_RNDN);
                mpfr_mul(t, t, tol, MPFR_RNDN);
                if (mpfr_cmp(tmp, t) <= 0) continue;    /* converged enough */
                rotated = true;
                /* tau = (alpha - gamma) / (2 * beta).
                 *
                 * Sign matches the rotation convention used in
                 * svdm_apply_rotation:
                 *   [new_i, new_j] = [a_i, a_j] . [c -s; s c]
                 * which makes the off-diagonal of the rotated 2x2
                 * block (c^2 - s^2) beta + cs (gamma - alpha).
                 * Zeroing it gives cot(2 theta) = (alpha - gamma)/(2 beta). */
                mpfr_sub(tau, alpha, gamma, MPFR_RNDN);
                mpfr_mul_2si(beta, beta, 1, MPFR_RNDN); /* 2 beta */
                mpfr_div(tau, tau, beta, MPFR_RNDN);
                /* t = sign(tau) / (|tau| + sqrt(1 + tau^2)) */
                mpfr_abs(abs_t, tau, MPFR_RNDN);
                mpfr_sqr(tmp, tau, MPFR_RNDN);
                mpfr_add_ui(tmp, tmp, 1, MPFR_RNDN);
                mpfr_sqrt(tmp, tmp, MPFR_RNDN);
                mpfr_add(tmp, abs_t, tmp, MPFR_RNDN);
                mpfr_set_si(t, mpfr_sgn(tau) >= 0 ? 1 : -1, MPFR_RNDN);
                mpfr_div(t, t, tmp, MPFR_RNDN);
                /* c = 1 / sqrt(1 + t^2) */
                mpfr_sqr(tmp, t, MPFR_RNDN);
                mpfr_add_ui(tmp, tmp, 1, MPFR_RNDN);
                mpfr_sqrt(tmp, tmp, MPFR_RNDN);
                mpfr_si_div(c, 1, tmp, MPFR_RNDN);
                /* s = t * c */
                mpfr_mul(s, t, c, MPFR_RNDN);
                /* Apply to A's columns i, j and V's columns i, j. */
                svdm_apply_rotation(A, n, i, j, c, s, bits);
                if (V) svdm_apply_rotation(V, p, i, j, c, s, bits);  /* NULL: values only */
            }
        }
        if (!rotated) break;
    }

    mpfr_clear(alpha); mpfr_clear(gamma); mpfr_clear(beta);
    mpfr_clear(tau);   mpfr_clear(t);     mpfr_clear(c);
    mpfr_clear(s);     mpfr_clear(abs_t); mpfr_clear(tmp);
    mpfr_clear(tol);
    return sweep;
}

/* Forward declaration: defined alongside the real wrapper.  Used by
 * svdm_wrap_cm_complex below. */
static Expr* svdm_mpfr_to_expr(const mpfr_t v, mpfr_prec_t bits);

/* ---------------------------------------------------------------------
 * Complex one-sided Jacobi SVD support.
 *
 * The complex 2x2 Hermitian Gram block looks like
 *   G = [alpha   beta  ]    alpha, gamma in R+,   beta in C
 *       [beta^*  gamma ]
 * with alpha = ||A[:,i]||^2, gamma = ||A[:,j]||^2,
 * beta = <A[:,i], A[:,j]>_H = sum_r conj(A[r,i]) * A[r,j].
 *
 * Diagonalise via a phase rotation that makes beta real positive,
 * followed by a real Jacobi rotation in the same column pair.  Writing
 * e^(-i phi) = (beta_re - i beta_im) / |beta| the combined 2x2 unitary
 * is
 *   U = [ c          -s          ]
 *       [ s e^(-iphi)  c e^(-iphi)]
 * which is verified unitary: c^2+s^2 = 1 and the cross diagonals
 * cancel.  Applied to columns (i, j) of A (and the same to V) it
 * zeroes the (i, j) entry of the Gram matrix to working precision in
 * one sweep step. */

/* Complex Hermitian inner product <A[:,i], A[:,j]>_H = sum_r
 *   conj(A[r,i]) * A[r,j] = (Ai_re Aj_re + Ai_im Aj_im)
 *                          + i (Ai_re Aj_im - Ai_im Aj_re).
 * Returns the real and imaginary parts in out_re, out_im. */
static void svdm_inner_complex(const mpfr_t* A_re, const mpfr_t* A_im,
                                int n, int i, int j,
                                mpfr_t out_re, mpfr_t out_im,
                                mpfr_prec_t bits)
{
    mpfr_set_zero(out_re, 1);
    mpfr_set_zero(out_im, 1);
    mpfr_t tmp; mpfr_init2(tmp, bits);
    for (int r = 0; r < n; r++) {
        size_t oi = (size_t)r + (size_t)i * (size_t)n;
        size_t oj = (size_t)r + (size_t)j * (size_t)n;
        mpfr_mul(tmp, A_re[oi], A_re[oj], MPFR_RNDN);
        mpfr_add(out_re, out_re, tmp, MPFR_RNDN);
        mpfr_mul(tmp, A_im[oi], A_im[oj], MPFR_RNDN);
        mpfr_add(out_re, out_re, tmp, MPFR_RNDN);
        mpfr_mul(tmp, A_re[oi], A_im[oj], MPFR_RNDN);
        mpfr_add(out_im, out_im, tmp, MPFR_RNDN);
        mpfr_mul(tmp, A_im[oi], A_re[oj], MPFR_RNDN);
        mpfr_sub(out_im, out_im, tmp, MPFR_RNDN);
    }
    mpfr_clear(tmp);
}

/* Squared column norm ||A[:, i]||^2 = sum_r (A_re[r,i]^2 + A_im[r,i]^2). */
static void svdm_norm_sq_complex(const mpfr_t* A_re, const mpfr_t* A_im,
                                  int n, int i,
                                  mpfr_t out, mpfr_prec_t bits)
{
    mpfr_set_zero(out, 1);
    mpfr_t tmp; mpfr_init2(tmp, bits);
    for (int r = 0; r < n; r++) {
        size_t off = (size_t)r + (size_t)i * (size_t)n;
        mpfr_sqr(tmp, A_re[off], MPFR_RNDN);
        mpfr_add(out, out, tmp, MPFR_RNDN);
        mpfr_sqr(tmp, A_im[off], MPFR_RNDN);
        mpfr_add(out, out, tmp, MPFR_RNDN);
    }
    mpfr_clear(tmp);
}

/* Apply the complex column rotation
 *   v_i_new = c v_i + s   (ec aj_re + es aj_im,  ec aj_im - es aj_re)
 *   v_j_new = -s v_i + c (ec aj_re + es aj_im,  ec aj_im - es aj_re)
 * where ec = beta_re / |beta|, es = beta_im / |beta| (so
 * e^(-i phi) = ec - i es and (e^(-iphi)) * aj is the bracketed pair).
 * Updates columns i and j of an n x ? paired (re, im) matrix C. */
static void svdm_apply_rotation_complex(mpfr_t* C_re, mpfr_t* C_im,
                                          int n, int col_i, int col_j,
                                          mpfr_t c, mpfr_t s,
                                          mpfr_t ec, mpfr_t es,
                                          mpfr_prec_t bits)
{
    mpfr_t ai_re, ai_im, aj_re, aj_im, ezj_re, ezj_im, t1, t2;
    mpfr_init2(ai_re, bits); mpfr_init2(ai_im, bits);
    mpfr_init2(aj_re, bits); mpfr_init2(aj_im, bits);
    mpfr_init2(ezj_re, bits); mpfr_init2(ezj_im, bits);
    mpfr_init2(t1, bits); mpfr_init2(t2, bits);
    for (int r = 0; r < n; r++) {
        size_t oi = (size_t)r + (size_t)col_i * (size_t)n;
        size_t oj = (size_t)r + (size_t)col_j * (size_t)n;
        mpfr_set(ai_re, C_re[oi], MPFR_RNDN);
        mpfr_set(ai_im, C_im[oi], MPFR_RNDN);
        mpfr_set(aj_re, C_re[oj], MPFR_RNDN);
        mpfr_set(aj_im, C_im[oj], MPFR_RNDN);
        /* ezj = e^(-i phi) * aj = (ec aj_re + es aj_im) + i (ec aj_im - es aj_re) */
        mpfr_mul(t1, ec, aj_re, MPFR_RNDN);
        mpfr_mul(t2, es, aj_im, MPFR_RNDN);
        mpfr_add(ezj_re, t1, t2, MPFR_RNDN);
        mpfr_mul(t1, ec, aj_im, MPFR_RNDN);
        mpfr_mul(t2, es, aj_re, MPFR_RNDN);
        mpfr_sub(ezj_im, t1, t2, MPFR_RNDN);
        /* v_i_new = c * ai + s * ezj */
        mpfr_mul(t1, c, ai_re, MPFR_RNDN);
        mpfr_mul(t2, s, ezj_re, MPFR_RNDN);
        mpfr_add(C_re[oi], t1, t2, MPFR_RNDN);
        mpfr_mul(t1, c, ai_im, MPFR_RNDN);
        mpfr_mul(t2, s, ezj_im, MPFR_RNDN);
        mpfr_add(C_im[oi], t1, t2, MPFR_RNDN);
        /* v_j_new = -s * ai + c * ezj */
        mpfr_mul(t1, c, ezj_re, MPFR_RNDN);
        mpfr_mul(t2, s, ai_re, MPFR_RNDN);
        mpfr_sub(C_re[oj], t1, t2, MPFR_RNDN);
        mpfr_mul(t1, c, ezj_im, MPFR_RNDN);
        mpfr_mul(t2, s, ai_im, MPFR_RNDN);
        mpfr_sub(C_im[oj], t1, t2, MPFR_RNDN);
    }
    mpfr_clear(ai_re); mpfr_clear(ai_im);
    mpfr_clear(aj_re); mpfr_clear(aj_im);
    mpfr_clear(ezj_re); mpfr_clear(ezj_im);
    mpfr_clear(t1); mpfr_clear(t2);
}

/* Complex one-sided Jacobi sweeps.  A is n x p paired (re, im),
 * V is p x p paired (re, im) and starts as the identity.  Returns the
 * number of sweeps actually performed. */
static int svdm_jacobi_sweeps_complex(mpfr_t* A_re, mpfr_t* A_im,
                                       int n, int p,
                                       mpfr_t* V_re, mpfr_t* V_im,
                                       mpfr_prec_t bits)
{
    mpfr_t alpha, gamma, beta_re, beta_im, abs_beta, abs_beta2;
    mpfr_t ec, es, tau, t, c, s, abs_t, tmp;
    mpfr_init2(alpha, bits);    mpfr_init2(gamma, bits);
    mpfr_init2(beta_re, bits);  mpfr_init2(beta_im, bits);
    mpfr_init2(abs_beta, bits); mpfr_init2(abs_beta2, bits);
    mpfr_init2(ec, bits);       mpfr_init2(es, bits);
    mpfr_init2(tau, bits);      mpfr_init2(t, bits);
    mpfr_init2(c, bits);        mpfr_init2(s, bits);
    mpfr_init2(abs_t, bits);    mpfr_init2(tmp, bits);

    mpfr_t tol; mpfr_init2(tol, bits);       /* n * 2^-bits, as the real kernel */
    mpfr_set_si(tol, n, MPFR_RNDN);
    mpfr_mul_2si(tol, tol, -(long)bits, MPFR_RNDN);

    int max_sweeps = 30;
    int sweep = 0;
    for (sweep = 0; sweep < max_sweeps; sweep++) {
        bool rotated = false;
        for (int i = 0; i < p - 1; i++) {
            for (int j = i + 1; j < p; j++) {
                svdm_norm_sq_complex(A_re, A_im, n, i, alpha, bits);
                svdm_norm_sq_complex(A_re, A_im, n, j, gamma, bits);
                svdm_inner_complex(A_re, A_im, n, i, j,
                                     beta_re, beta_im, bits);
                /* |beta|^2 = beta_re^2 + beta_im^2 */
                mpfr_sqr(abs_beta2, beta_re, MPFR_RNDN);
                mpfr_sqr(tmp,       beta_im, MPFR_RNDN);
                mpfr_add(abs_beta2, abs_beta2, tmp, MPFR_RNDN);
                mpfr_sqrt(abs_beta, abs_beta2, MPFR_RNDN);
                /* Convergence: 2 |beta| <= tol * sqrt(alpha * gamma). */
                mpfr_mul_2si(tmp, abs_beta, 1, MPFR_RNDN);    /* 2|beta| */
                mpfr_mul(t, alpha, gamma, MPFR_RNDN);
                mpfr_sqrt(t, t, MPFR_RNDN);
                mpfr_mul(t, t, tol, MPFR_RNDN);
                if (mpfr_cmp(tmp, t) <= 0) continue;
                rotated = true;
                /* Phase rotation: ec = beta_re/|beta|, es = beta_im/|beta|. */
                mpfr_div(ec, beta_re, abs_beta, MPFR_RNDN);
                mpfr_div(es, beta_im, abs_beta, MPFR_RNDN);
                /* After making beta real-positive (value = |beta|) the
                 * real Jacobi rotation diagonalises [alpha |beta|; |beta| gamma]:
                 *   tau = (alpha - gamma) / (2 |beta|)
                 *   t   = sign(tau) / (|tau| + sqrt(1 + tau^2))
                 *   c   = 1/sqrt(1 + t^2), s = t c. */
                mpfr_sub(tau, alpha, gamma, MPFR_RNDN);
                mpfr_mul_2si(tmp, abs_beta, 1, MPFR_RNDN);
                mpfr_div(tau, tau, tmp, MPFR_RNDN);
                mpfr_abs(abs_t, tau, MPFR_RNDN);
                mpfr_sqr(tmp, tau, MPFR_RNDN);
                mpfr_add_ui(tmp, tmp, 1, MPFR_RNDN);
                mpfr_sqrt(tmp, tmp, MPFR_RNDN);
                mpfr_add(tmp, abs_t, tmp, MPFR_RNDN);
                mpfr_set_si(t, mpfr_sgn(tau) >= 0 ? 1 : -1, MPFR_RNDN);
                mpfr_div(t, t, tmp, MPFR_RNDN);
                mpfr_sqr(tmp, t, MPFR_RNDN);
                mpfr_add_ui(tmp, tmp, 1, MPFR_RNDN);
                mpfr_sqrt(tmp, tmp, MPFR_RNDN);
                mpfr_si_div(c, 1, tmp, MPFR_RNDN);
                mpfr_mul(s, t, c, MPFR_RNDN);
                /* Apply complex rotation to A's and V's columns i, j. */
                svdm_apply_rotation_complex(A_re, A_im, n, i, j,
                                              c, s, ec, es, bits);
                if (V_re)                                  /* NULL: values only */
                    svdm_apply_rotation_complex(V_re, V_im, p, i, j,
                                                  c, s, ec, es, bits);
            }
        }
        if (!rotated) break;
    }

    mpfr_clear(alpha);    mpfr_clear(gamma);
    mpfr_clear(beta_re);  mpfr_clear(beta_im);
    mpfr_clear(abs_beta); mpfr_clear(abs_beta2);
    mpfr_clear(ec); mpfr_clear(es);
    mpfr_clear(tau); mpfr_clear(t); mpfr_clear(c); mpfr_clear(s);
    mpfr_clear(abs_t); mpfr_clear(tmp); mpfr_clear(tol);
    return sweep;
}

/* Extract sigma and U columns from post-Jacobi complex A.
 *   sigma[j] = ||A[:, j]||  (real, non-negative)
 *   U[:, j]  = A[:, j] / sigma[j] when sigma[j] > 0
 *              left as zero        otherwise. */
static void svdm_extract_U_sigma_complex(const mpfr_t* A_re,
                                           const mpfr_t* A_im,
                                           int n, int p,
                                           mpfr_t* U_re, mpfr_t* U_im,
                                           mpfr_t* sigma,
                                           mpfr_prec_t bits)
{
    mpfr_t nsq, norm, inv;
    mpfr_init2(nsq, bits); mpfr_init2(norm, bits); mpfr_init2(inv, bits);
    for (int j = 0; j < p; j++) {
        svdm_norm_sq_complex(A_re, A_im, n, j, nsq, bits);
        mpfr_sqrt(norm, nsq, MPFR_RNDN);
        mpfr_set(sigma[j], norm, MPFR_RNDN);
        if (mpfr_zero_p(norm)) {
            for (int r = 0; r < n; r++) {
                size_t off = (size_t)r + (size_t)j * (size_t)n;
                mpfr_set_zero(U_re[off], 1);
                mpfr_set_zero(U_im[off], 1);
            }
            continue;
        }
        mpfr_si_div(inv, 1, norm, MPFR_RNDN);
        for (int r = 0; r < n; r++) {
            size_t off = (size_t)r + (size_t)j * (size_t)n;
            mpfr_mul(U_re[off], A_re[off], inv, MPFR_RNDN);
            mpfr_mul(U_im[off], A_im[off], inv, MPFR_RNDN);
        }
    }
    mpfr_clear(nsq); mpfr_clear(norm); mpfr_clear(inv);
}

/* Sort sigma in descending order, swapping U and V paired columns. */
static void svdm_sort_descending_complex(mpfr_t* sigma,
                                          mpfr_t* U_re, mpfr_t* U_im,
                                          int n,
                                          mpfr_t* V_re, mpfr_t* V_im,
                                          int p)
{
    for (int i = 1; i < p; i++) {
        for (int j = i; j > 0; j--) {
            if (mpfr_cmp(sigma[j], sigma[j - 1]) <= 0) break;
            mpfr_swap(sigma[j], sigma[j - 1]);
            for (int r = 0; r < n; r++) {
                size_t a = (size_t)r + (size_t)(j - 1) * (size_t)n;
                size_t b = (size_t)r + (size_t)j * (size_t)n;
                mpfr_swap(U_re[a], U_re[b]);
                mpfr_swap(U_im[a], U_im[b]);
            }
            for (int r = 0; r < p; r++) {
                size_t a = (size_t)r + (size_t)(j - 1) * (size_t)p;
                size_t b = (size_t)r + (size_t)j * (size_t)p;
                mpfr_swap(V_re[a], V_re[b]);
                mpfr_swap(V_im[a], V_im[b]);
            }
        }
    }
}

/* Gram-Schmidt orthonormal completion for complex U (paired re/im).
 * Identical to the real path but with Hermitian inner products and
 * complex scaling.  Standard basis seeds are real (im = 0). */
static int svdm_complete_U_complex(mpfr_t* U_re, mpfr_t* U_im,
                                    mpfr_t* sigma, int n, int p,
                                    mpfr_prec_t bits)
{
    int rank = 0;
    for (int j = 0; j < p; j++) {
        if (!mpfr_zero_p(sigma[j])) rank++;
    }
    int filled = rank;
    if (filled >= n) return rank;

    mpfr_t coeff_re, coeff_im, tmp1, tmp2, nsq, norm, inv;
    mpfr_init2(coeff_re, bits); mpfr_init2(coeff_im, bits);
    mpfr_init2(tmp1, bits);     mpfr_init2(tmp2, bits);
    mpfr_init2(nsq, bits);      mpfr_init2(norm, bits);
    mpfr_init2(inv, bits);

    mpfr_t* v_re = svdm_array_alloc((size_t)n, bits);
    mpfr_t* v_im = svdm_array_alloc((size_t)n, bits);

    for (int seed = 0; seed < n && filled < n; seed++) {
        for (int r = 0; r < n; r++) {
            mpfr_set_si(v_re[r], r == seed ? 1 : 0, MPFR_RNDN);
            mpfr_set_zero(v_im[r], 1);
        }
        /* v -= sum_j <U[:, j], v>_H * U[:, j] for j in [0, filled).
         * <U[:, j], v>_H = sum_r conj(U[r,j]) * v[r]
         *               = sum_r (U_re v_re + U_im v_im)
         *                  + i (U_re v_im - U_im v_re). */
        for (int j = 0; j < filled; j++) {
            mpfr_set_zero(coeff_re, 1);
            mpfr_set_zero(coeff_im, 1);
            for (int r = 0; r < n; r++) {
                size_t off = (size_t)r + (size_t)j * (size_t)n;
                mpfr_mul(tmp1, U_re[off], v_re[r], MPFR_RNDN);
                mpfr_add(coeff_re, coeff_re, tmp1, MPFR_RNDN);
                mpfr_mul(tmp1, U_im[off], v_im[r], MPFR_RNDN);
                mpfr_add(coeff_re, coeff_re, tmp1, MPFR_RNDN);
                mpfr_mul(tmp1, U_re[off], v_im[r], MPFR_RNDN);
                mpfr_add(coeff_im, coeff_im, tmp1, MPFR_RNDN);
                mpfr_mul(tmp1, U_im[off], v_re[r], MPFR_RNDN);
                mpfr_sub(coeff_im, coeff_im, tmp1, MPFR_RNDN);
            }
            /* v -= coeff * U[:, j], complex multiplication. */
            for (int r = 0; r < n; r++) {
                size_t off = (size_t)r + (size_t)j * (size_t)n;
                mpfr_mul(tmp1, coeff_re, U_re[off], MPFR_RNDN);
                mpfr_mul(tmp2, coeff_im, U_im[off], MPFR_RNDN);
                mpfr_sub(tmp1, tmp1, tmp2, MPFR_RNDN);
                mpfr_sub(v_re[r], v_re[r], tmp1, MPFR_RNDN);
                mpfr_mul(tmp1, coeff_re, U_im[off], MPFR_RNDN);
                mpfr_mul(tmp2, coeff_im, U_re[off], MPFR_RNDN);
                mpfr_add(tmp1, tmp1, tmp2, MPFR_RNDN);
                mpfr_sub(v_im[r], v_im[r], tmp1, MPFR_RNDN);
            }
        }
        mpfr_set_zero(nsq, 1);
        for (int r = 0; r < n; r++) {
            mpfr_sqr(tmp1, v_re[r], MPFR_RNDN);
            mpfr_add(nsq, nsq, tmp1, MPFR_RNDN);
            mpfr_sqr(tmp1, v_im[r], MPFR_RNDN);
            mpfr_add(nsq, nsq, tmp1, MPFR_RNDN);
        }
        mpfr_set_si(norm, 1, MPFR_RNDN);
        mpfr_mul_2si(norm, norm, -(long)bits, MPFR_RNDN);
        if (mpfr_cmp(nsq, norm) <= 0) continue;
        mpfr_sqrt(norm, nsq, MPFR_RNDN);
        mpfr_si_div(inv, 1, norm, MPFR_RNDN);
        for (int r = 0; r < n; r++) {
            size_t off = (size_t)r + (size_t)filled * (size_t)n;
            mpfr_mul(U_re[off], v_re[r], inv, MPFR_RNDN);
            mpfr_mul(U_im[off], v_im[r], inv, MPFR_RNDN);
        }
        filled++;
    }

    svdm_array_free(v_re, (size_t)n);
    svdm_array_free(v_im, (size_t)n);
    mpfr_clear(coeff_re); mpfr_clear(coeff_im);
    mpfr_clear(tmp1); mpfr_clear(tmp2);
    mpfr_clear(nsq); mpfr_clear(norm); mpfr_clear(inv);
    return rank;
}

/* Wrap a paired (re, im) column-major MPFR buffer as a row-major
 * Mathilda matrix.  Cells with a non-zero imaginary part emit
 * Complex[re, im]; pure-real cells stay as MPFR scalars (so downstream
 * Dot/Norm doesn't drag Complex through real-only computations). */
static Expr* svdm_wrap_cm_complex(const mpfr_t* buf_re, const mpfr_t* buf_im,
                                    int rows, int cols, int lda,
                                    mpfr_prec_t bits)
{
    Expr** row_exprs = (Expr**)malloc(sizeof(Expr*) * (size_t)rows);
    for (int i = 0; i < rows; i++) {
        Expr** elems = (Expr**)malloc(sizeof(Expr*) * (size_t)cols);
        for (int j = 0; j < cols; j++) {
            size_t off = (size_t)i + (size_t)j * (size_t)lda;
            if (mpfr_zero_p(buf_im[off])) {
                elems[j] = svdm_mpfr_to_expr(buf_re[off], bits);
            } else {
                Expr** ca = (Expr**)malloc(sizeof(Expr*) * 2);
                ca[0] = svdm_mpfr_to_expr(buf_re[off], bits);
                ca[1] = svdm_mpfr_to_expr(buf_im[off], bits);
                elems[j] = expr_new_function(expr_new_symbol(SYM_Complex),
                                               ca, 2);
                free(ca);
            }
        }
        row_exprs[i] = expr_new_function(expr_new_symbol(SYM_List),
                                          elems, (size_t)cols);
        free(elems);
    }
    Expr* m = expr_new_function(expr_new_symbol(SYM_List),
                                  row_exprs, (size_t)rows);
    free(row_exprs);
    return m;
}

/* ---------------------------------------------------------------------
 * Extract singular values and left singular vectors from the post-
 * Jacobi A.  On entry A's columns are sigma_i u_i; on exit:
 *   sigma[i] = || A[:, i] ||
 *   U[:, i]  = A[:, i] / sigma[i]  when sigma[i] > 0;
 *              left as zero          when sigma[i] = 0 (we fix those
 *              up later via the orthonormal-completion pass on U).
 * ------------------------------------------------------------------ */
static void svdm_extract_U_sigma(const mpfr_t* A, int n, int p,
                                 mpfr_t* U, mpfr_t* sigma,
                                 mpfr_prec_t bits)
{
    mpfr_t nsq, norm, inv;
    mpfr_init2(nsq, bits);
    mpfr_init2(norm, bits);
    mpfr_init2(inv, bits);
    for (int j = 0; j < p; j++) {
        svdm_inner(A, n, j, j, nsq, bits);
        mpfr_sqrt(norm, nsq, MPFR_RNDN);
        mpfr_set(sigma[j], norm, MPFR_RNDN);
        if (mpfr_zero_p(norm)) {
            for (int r = 0; r < n; r++) {
                size_t off = (size_t)r + (size_t)j * (size_t)n;
                mpfr_set_zero(U[off], 1);
            }
            continue;
        }
        mpfr_si_div(inv, 1, norm, MPFR_RNDN);
        for (int r = 0; r < n; r++) {
            size_t oa = (size_t)r + (size_t)j * (size_t)n;
            size_t ou = (size_t)r + (size_t)j * (size_t)n;
            mpfr_mul(U[ou], A[oa], inv, MPFR_RNDN);
        }
    }
    mpfr_clear(nsq); mpfr_clear(norm); mpfr_clear(inv);
}

/* ---------------------------------------------------------------------
 * Sort singular values in descending order; carry along U and V
 * columns.  Simple insertion sort (p is typically small for MPFR
 * workloads). */
static void svdm_sort_descending(mpfr_t* sigma, mpfr_t* U, int n,
                                 mpfr_t* V, int p)
{
    for (int i = 1; i < p; i++) {
        for (int j = i; j > 0; j--) {
            if (mpfr_cmp(sigma[j], sigma[j - 1]) <= 0) break;
            mpfr_swap(sigma[j], sigma[j - 1]);
            for (int r = 0; r < n; r++) {
                size_t a = (size_t)r + (size_t)(j - 1) * (size_t)n;
                size_t b = (size_t)r + (size_t)j * (size_t)n;
                mpfr_swap(U[a], U[b]);
            }
            for (int r = 0; r < p; r++) {
                size_t a = (size_t)r + (size_t)(j - 1) * (size_t)p;
                size_t b = (size_t)r + (size_t)j * (size_t)p;
                mpfr_swap(V[a], V[b]);
            }
        }
    }
}

/* ---------------------------------------------------------------------
 * Orthonormal completion of U.
 *
 * After Jacobi + extract: U has p columns (or min(n, p) -- we operate
 * uniformly on the first `built` columns where built == p in the
 * non-truncated case).  When the matrix has rank < n we need
 * n - rank more orthonormal columns to fill U out to its full n x n
 * shape; if the existing U columns are already a complete basis of
 * R^n (rank == n) the completion is a no-op.
 *
 * Algorithm: classical Gram-Schmidt against standard basis vectors,
 * skipping vectors with a too-small post-projection norm.  At MPFR
 * precision this is fast and numerically clean (much simpler than the
 * symbolic case which has to fight Together blow-up). */
static int svdm_complete_U(mpfr_t* U, mpfr_t* sigma, int n, int p,
                           mpfr_prec_t bits)
{
    /* Count the non-zero singular values. */
    int rank = 0;
    for (int j = 0; j < p; j++) {
        if (!mpfr_zero_p(sigma[j])) rank++;
    }
    /* Compact U so its first `rank` columns are the non-zero ones.
     * Carry the sigmas along to keep alignment.  After Jacobi + sort,
     * non-zero sigmas already come first, so this is a no-op in the
     * typical case -- but explicit to be safe. */

    int filled = rank;
    if (filled >= n) return rank;

    /* Try each standard basis vector e_i; project out components along
     * the existing orthonormal columns of U; if the residual norm is
     * non-zero, normalise and add as the next column. */
    mpfr_t coeff, tmp, nsq, norm, inv;
    mpfr_init2(coeff, bits); mpfr_init2(tmp,  bits);
    mpfr_init2(nsq,   bits); mpfr_init2(norm, bits);
    mpfr_init2(inv,   bits);

    mpfr_t* v = svdm_array_alloc((size_t)n, bits);

    /* U has space for n columns even though only the first `p` were
     * initialised.  We grow into columns [rank, n). */
    for (int seed = 0; seed < n && filled < n; seed++) {
        for (int r = 0; r < n; r++) {
            mpfr_set_si(v[r], r == seed ? 1 : 0, MPFR_RNDN);
        }
        /* v -= sum_j <U[:, j], v> * U[:, j] for j in [0, filled). */
        for (int j = 0; j < filled; j++) {
            mpfr_set_zero(coeff, 1);
            for (int r = 0; r < n; r++) {
                size_t off = (size_t)r + (size_t)j * (size_t)n;
                mpfr_mul(tmp, U[off], v[r], MPFR_RNDN);
                mpfr_add(coeff, coeff, tmp, MPFR_RNDN);
            }
            for (int r = 0; r < n; r++) {
                size_t off = (size_t)r + (size_t)j * (size_t)n;
                mpfr_mul(tmp, coeff, U[off], MPFR_RNDN);
                mpfr_sub(v[r], v[r], tmp, MPFR_RNDN);
            }
        }
        /* Norm test: |v|^2 > 2^(-bits) (avoid spurious tiny vectors). */
        mpfr_set_zero(nsq, 1);
        for (int r = 0; r < n; r++) {
            mpfr_sqr(tmp, v[r], MPFR_RNDN);
            mpfr_add(nsq, nsq, tmp, MPFR_RNDN);
        }
        mpfr_set_si(norm, 1, MPFR_RNDN);
        mpfr_mul_2si(norm, norm, -(long)bits, MPFR_RNDN);
        if (mpfr_cmp(nsq, norm) <= 0) continue;
        mpfr_sqrt(norm, nsq, MPFR_RNDN);
        mpfr_si_div(inv, 1, norm, MPFR_RNDN);
        for (int r = 0; r < n; r++) {
            size_t off = (size_t)r + (size_t)filled * (size_t)n;
            mpfr_mul(U[off], v[r], inv, MPFR_RNDN);
        }
        filled++;
    }

    svdm_array_free(v, (size_t)n);
    mpfr_clear(coeff); mpfr_clear(tmp);   mpfr_clear(nsq);
    mpfr_clear(norm);  mpfr_clear(inv);
    return rank;
}

/* Convert an MPFR cell to an Expr*: EXPR_MPFR at the working precision
 * for non-zero values; an explicit MPFR zero (also at the working
 * precision) for hard zeros.  Using expr_new_mpfr_copy preserves the
 * input precision and follows the same construction convention as the
 * other linalg MPFR kernels. */
static Expr* svdm_mpfr_to_expr(const mpfr_t v, mpfr_prec_t bits)
{
    (void)bits;
    return expr_new_mpfr_copy(v);
}

/* Wrap an n x p column-major MPFR buffer as a row-major Mathilda
 * matrix.  Steals nothing -- mpfr_set copies cell-by-cell. */
static Expr* svdm_wrap_cm(const mpfr_t* buf, int rows, int cols, int lda,
                          mpfr_prec_t bits)
{
    Expr** row_exprs = (Expr**)malloc(sizeof(Expr*) * (size_t)rows);
    for (int i = 0; i < rows; i++) {
        Expr** elems = (Expr**)malloc(sizeof(Expr*) * (size_t)cols);
        for (int j = 0; j < cols; j++) {
            size_t off = (size_t)i + (size_t)j * (size_t)lda;
            elems[j] = svdm_mpfr_to_expr(buf[off], bits);
        }
        row_exprs[i] = expr_new_function(expr_new_symbol(SYM_List),
                                          elems, (size_t)cols);
        free(elems);
    }
    Expr* m = expr_new_function(expr_new_symbol(SYM_List),
                                  row_exprs, (size_t)rows);
    free(row_exprs);
    return m;
}

/* Build the rectangular n x p Sigma matrix from the length-p sigma
 * vector.  Sigma[i][j] = sigma_i when i == j < min(n, p); MPFR zero at
 * the working precision elsewhere.  Crucially we use MPFR zeros (not
 * Real 0.0) so that downstream Dot products stay at the input
 * precision -- mixing MPFR and machine Real demotes the result to
 * MachinePrecision in Mathilda's arithmetic. */
static Expr* svdm_sigma_rect(const mpfr_t* sigma, int n, int p,
                             mpfr_prec_t bits)
{
    int mn = (n < p) ? n : p;
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    for (int i = 0; i < n; i++) {
        Expr** elems = (Expr**)malloc(sizeof(Expr*) * (size_t)p);
        for (int j = 0; j < p; j++) {
            if (i == j && i < mn) {
                elems[j] = svdm_mpfr_to_expr(sigma[i], bits);
            } else {
                elems[j] = expr_new_mpfr_from_si(0, bits);
            }
        }
        rows[i] = expr_new_function(expr_new_symbol(SYM_List),
                                      elems, (size_t)p);
        free(elems);
    }
    Expr* m = expr_new_function(expr_new_symbol(SYM_List),
                                  rows, (size_t)n);
    free(rows);
    return m;
}

/* ------------------------------------------------------------------ *
 *  Thin-or-square Jacobi core (n >= p).                               *
 *                                                                     *
 *  One-sided Jacobi orthogonalises the p columns of the n x p matrix  *
 *  m, so it needs p <= n: U is n x n and receives one column per      *
 *  column of A.  (A wide matrix is handed in conjugate-transposed by  *
 *  the dispatcher -- writing p > n columns into an n x n U overran    *
 *  the buffer.)  On success *U_out is the n x n U, *V_out the p x p V *
 *  and *sigma_out the p singular values in descending order (caller   *
 *  frees with svdm_array_free(.., p)).                                *
 * ------------------------------------------------------------------ */
static bool svdm_core_complex(Expr* m, int n, int p, mpfr_prec_t bits,
                              Expr** U_out, Expr** V_out,
                              mpfr_t** sigma_out)
{
    mpfr_t* A_re = NULL;
    mpfr_t* A_im = NULL;
    if (!svdm_load_complex(m, n, p, bits, &A_re, &A_im)) {
        return false;
    }

    size_t Vsz = (size_t)p * (size_t)p;
    size_t Usz = (size_t)n * (size_t)n;
    mpfr_t* V_re = svdm_array_alloc(Vsz, bits);
    mpfr_t* V_im = svdm_array_alloc(Vsz, bits);
    for (int j = 0; j < p; j++) {
        for (int i = 0; i < p; i++) {
            size_t off = (size_t)i + (size_t)j * (size_t)p;
            mpfr_set_si(V_re[off], i == j ? 1 : 0, MPFR_RNDN);
            mpfr_set_zero(V_im[off], 1);
        }
    }
    mpfr_t* U_re = svdm_array_alloc(Usz, bits);
    mpfr_t* U_im = svdm_array_alloc(Usz, bits);
    for (size_t k = 0; k < Usz; k++) {
        mpfr_set_zero(U_re[k], 1);
        mpfr_set_zero(U_im[k], 1);
    }
    mpfr_t* sigma = svdm_array_alloc((size_t)p, bits);

    (void)svdm_jacobi_sweeps_complex(A_re, A_im, n, p, V_re, V_im, bits);
    svdm_extract_U_sigma_complex(A_re, A_im, n, p, U_re, U_im, sigma, bits);
    svdm_sort_descending_complex(sigma, U_re, U_im, n, V_re, V_im, p);
    (void)svdm_complete_U_complex(U_re, U_im, sigma, n, p, bits);

    *U_out = svdm_wrap_cm_complex(U_re, U_im, n, n, n, bits);
    *V_out = svdm_wrap_cm_complex(V_re, V_im, p, p, p, bits);
    *sigma_out = sigma;

    svdm_array_free(A_re, (size_t)n * (size_t)p);
    svdm_array_free(A_im, (size_t)n * (size_t)p);
    svdm_array_free(V_re, Vsz); svdm_array_free(V_im, Vsz);
    svdm_array_free(U_re, Usz); svdm_array_free(U_im, Usz);
    return true;
}

static bool svdm_core_real(Expr* m, int n, int p, mpfr_prec_t bits,
                           Expr** U_out, Expr** V_out, mpfr_t** sigma_out)
{
    mpfr_t* A = NULL;
    if (!svdm_load_real(m, n, p, bits, &A)) {
        return false;
    }

    /* V starts as the p x p identity. */
    mpfr_t* V = svdm_array_alloc((size_t)p * (size_t)p, bits);
    for (int j = 0; j < p; j++) {
        for (int i = 0; i < p; i++) {
            size_t off = (size_t)i + (size_t)j * (size_t)p;
            mpfr_set_si(V[off], i == j ? 1 : 0, MPFR_RNDN);
        }
    }

    /* U is all zeros so the unused n - p columns (when n > p) read as
     * zero before the orthonormal-completion pass fills them. */
    mpfr_t* U = svdm_array_alloc((size_t)n * (size_t)n, bits);
    for (size_t k = 0; k < (size_t)n * (size_t)n; k++) {
        mpfr_set_zero(U[k], 1);
    }
    mpfr_t* sigma = svdm_array_alloc((size_t)p, bits);

    (void)svdm_jacobi_sweeps(A, n, p, V, bits);
    svdm_extract_U_sigma(A, n, p, U, sigma, bits);
    svdm_sort_descending(sigma, U, n, V, p);
    (void)svdm_complete_U(U, sigma, n, p, bits);

    *U_out = svdm_wrap_cm(U, n, n, n, bits);
    *V_out = svdm_wrap_cm(V, p, p, p, bits);
    *sigma_out = sigma;

    svdm_array_free(A, (size_t)n * (size_t)p);
    svdm_array_free(V, (size_t)p * (size_t)p);
    svdm_array_free(U, (size_t)n * (size_t)n);
    return true;
}

/* ------------------------------------------------------------------ *
 *  Public entry: MPFR Jacobi SVD dispatcher.                          *
 * ------------------------------------------------------------------ */
Expr* svd_mpfr_dispatch(const SvdArgs* args, int n, int p, int n_a)
{
    if (args->generalized) {
        /* Generalized SVD at MPFR precision: Mathilda doesn't yet have
         * a native Paige/Van Loan QR-reduce + 2-pair Jacobi loop at
         * MPFR precision (the algorithm is ~400 LOC of dense linear
         * algebra and warrants its own change).  Until then we fall
         * back to symbolic dispatch, which numericalises to 53-bit
         * Reals and routes through LAPACK.  Emit a one-shot warning so
         * the precision downgrade is visible. */
        static uint64_t gmpfr_warn = 0;
        if (!gmpfr_warn) {
            gmpfr_warn = 1;
            mth_message("SingularValueDecomposition", "gmpdwn", "Generalized SVD of high-precision MPFR input is computed at machine precision via LAPACK -- a native MPFR generalized kernel is not yet implemented.");
        }
        (void)n_a;
        return NULL;   /* falls through to svd_symbolic_dispatch */
    }

    (void)n_a;
    CommonInexactInfo info = common_scan_inexact(args->m);
    if (!info.has_inexact || info.min_bits <= 53) {
        return NULL;
    }
    mpfr_prec_t bits = (mpfr_prec_t)info.min_bits;

    /* The Jacobi core needs n >= p.  A wide matrix is decomposed through
     * its conjugate transpose: A^H = U' S' V'^H  gives  A = V' S'^T U'^H,
     * so U = V' (n x n), V = U' (p x p), and the min(n, p) = n singular
     * values are shared. */
    bool wide = p > n;
    Expr* work = args->m;
    Expr* owned = NULL;
    if (wide) {
        Expr** ca = (Expr**)malloc(sizeof(Expr*));
        ca[0] = expr_copy(args->m);
        Expr* call = expr_new_function(expr_new_symbol(SYM_ConjugateTranspose),
                                       ca, 1);
        free(ca);
        owned = evaluate(call);
        expr_free(call);
        work = owned;
    }
    int wn = wide ? p : n;     /* rows of the matrix handed to the core */
    int wp = wide ? n : p;     /* its columns: wp <= wn                 */

    Expr* U = NULL;
    Expr* V = NULL;
    mpfr_t* sigma = NULL;
    bool ok = svdm_matrix_has_complex(work, wn, wp)
                ? svdm_core_complex(work, wn, wp, bits, &U, &V, &sigma)
                : svdm_core_real(work, wn, wp, bits, &U, &V, &sigma);
    if (owned) expr_free(owned);
    if (!ok) return NULL;

    Expr* sigma_mat = svdm_sigma_rect(sigma, n, p, bits);
    svdm_array_free(sigma, (size_t)wp);

    Expr** items = (Expr**)malloc(sizeof(Expr*) * 3);
    items[0] = wide ? V : U;
    items[1] = sigma_mat;
    items[2] = wide ? U : V;
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), items, 3);
    free(items);

    int mn = (n < p) ? n : p;
    return svd_apply_postprocess(result, args, n, p, mn);
}

#ifdef USE_FLINT
/* ------------------------------------------------------------------ *
 *  Singular values by bidiagonalisation (real input, values only).    *
 *                                                                     *
 *  One-sided Jacobi costs ~ 3 n p^2 / 2 multiply-adds per sweep and   *
 *  an ill-conditioned matrix takes 10+ sweeps (a 120 x 120 Kahan      *
 *  matrix at 60 digits: 4 s).  The standard route instead:            *
 *                                                                     *
 *   1. Golub-Kahan bidiagonalisation A = U B V^T by Householder       *
 *      reflectors from the left and the right (4 n p^2 - 4 p^3 / 3    *
 *      multiply-adds, once).  Every inner product is one fused        *
 *      arb_approx_dot (accumulated in fixed point, rounded once) and  *
 *      every entry update one fused multiply-add; the columns (left   *
 *      reflector) and rows (right reflector) are independent and are *
 *      split over threads (mp_parallel_for).                          *
 *   2. The singular values of the upper bidiagonal B by the implicit- *
 *      shift Golub-Kahan QR iteration (Golub & Van Loan, Alg. 8.6.2): *
 *      Wilkinson shift from the trailing 2 x 2 of B^T B, a column /   *
 *      row Givens bulge chase, deflation once |e_i| is negligible     *
 *      against |d_i| + |d_i+1|, and a zero diagonal entry chased out  *
 *      along its row (or up its column at the bottom of a block).     *
 *      O(p) per sweep, a few sweeps per value.                        *
 *                                                                     *
 *  Accuracy: both stages are backward stable, so each value is        *
 *  within ~ p 2^-w ||A||_2 of the exact one -- the guarantee of the   *
 *  machine kernel (LAPACK gesdd) at w bits.  The work runs at the     *
 *  caller's precision plus SVDB_GUARD bits and is rounded to it at    *
 *  the end.                                                           *
 * ------------------------------------------------------------------ */
#define SVDB_GUARD 32

typedef struct {
    arb_ptr A;            /* n x p row-major midpoints */
    arb_ptr u;            /* reflector: rows k.. (left) or cols k+1.. (right) */
    arb_ptr wv;           /* left reflector: 2 u . A[k.., j] per column j */
    size_t n, p, k;
    slong prec;
} SvdBi;

/* Left reflector on column k, rows k..n-1, first half: for each column
 * j > k, wv_j = 2 u . A[k.., j]. */
static void svdb_left_dot(void* c0, int64_t lo, int64_t hi, int chunk)
{
    SvdBi* c = (SvdBi*)c0;
    (void)chunk;
    size_t n = c->n, p = c->p, k = c->k;
    for (int64_t t = lo; t < hi; t++) {
        size_t j = k + 1 + (size_t)t;
        arb_approx_dot(c->wv + j, NULL, 0, c->u + k, 1, c->A + k * p + j, (slong)p,
                       (slong)(n - k), c->prec);
        arf_mul_2exp_si(arb_midref(c->wv + j), arb_midref(c->wv + j), 1);
    }
}

/* ... second half, by rows (contiguous): A[i, j] -= u_i wv_j for i >= k,
 * j > k. */
static void svdb_left_upd(void* c0, int64_t lo, int64_t hi, int chunk)
{
    SvdBi* c = (SvdBi*)c0;
    (void)chunk;
    size_t p = c->p, k = c->k;
    for (int64_t t = lo; t < hi; t++) {
        size_t i = k + (size_t)t;
        arb_ptr row = c->A + i * p;
        arf_srcptr ui = arb_midref(c->u + i);
        if (arf_is_zero(ui)) continue;
        for (size_t j = k + 1; j < p; j++)
            arf_submul(arb_midref(row + j), ui, arb_midref(c->wv + j), c->prec, ARF_RND_NEAR);
    }
}

/* Right reflector on row k, columns k+1..p-1: for each row i > k,
 * w = A[i, k+1..] . u,  A[i, k+1..] -= (2 w) u. */
static void svdb_right(void* c0, int64_t lo, int64_t hi, int chunk)
{
    SvdBi* c = (SvdBi*)c0;
    (void)chunk;
    size_t p = c->p, k = c->k;
    arb_t w;
    arb_init(w);
    for (int64_t t = lo; t < hi; t++) {
        size_t i = k + 1 + (size_t)t;
        arb_ptr row = c->A + i * p;
        arb_approx_dot(w, NULL, 0, row + k + 1, 1, c->u + k + 1, 1, (slong)(p - k - 1),
                       c->prec);
        arf_mul_2exp_si(arb_midref(w), arb_midref(w), 1);
        for (size_t j = k + 1; j < p; j++)
            arf_submul(arb_midref(row + j), arb_midref(c->u + j), arb_midref(w), c->prec, ARF_RND_NEAR);
    }
    arb_clear(w);
}

/* The unit Householder vector u (contiguous, len cells) with
 * (I - 2 u u^T) x = alpha e_1 for x = len entries at stride st.  False
 * (alpha = 0, no reflection) when x is zero. */
static bool svdb_house(arb_srcptr x, slong st, size_t len, arb_ptr u, arf_t alpha,
                       slong prec)
{
    arb_t s;
    arf_t nrm;
    arb_init(s); arf_init(nrm);
    arb_approx_dot(s, NULL, 0, x, st, x, st, (slong)len, prec);
    bool ok = !arf_is_zero(arb_midref(s));
    if (ok) {
        arf_sqrt(nrm, arb_midref(s), prec, ARF_RND_NEAR);
        if (arf_sgn(arb_midref(x)) >= 0) arf_neg(alpha, nrm); else arf_set(alpha, nrm);
        arf_sub(arb_midref(u), arb_midref(x), alpha, prec, ARF_RND_NEAR);
        for (size_t i = 1; i < len; i++)
            arf_set(arb_midref(u + i), arb_midref(x + (slong)i * st));
        arb_approx_dot(s, NULL, 0, u, 1, u, 1, (slong)len, prec);
        if (arf_is_zero(arb_midref(s))) ok = false;
        else {
            arf_sqrt(nrm, arb_midref(s), prec, ARF_RND_NEAR);
            for (size_t i = 0; i < len; i++)
                arf_div(arb_midref(u + i), arb_midref(u + i), nrm, prec, ARF_RND_NEAR);
        }
    }
    if (!ok) arf_zero(alpha);
    arb_clear(s); arf_clear(nrm);
    return ok;
}

/* Givens: r = hypot(f, g), c = f / r, s = g / r (c = 1, s = 0 for g = 0). */
static void svdb_givens(mpfr_t c, mpfr_t s, mpfr_t r, const mpfr_t f, const mpfr_t g)
{
    if (mpfr_zero_p(g)) {
        mpfr_set_ui(c, 1, MPFR_RNDN); mpfr_set_zero(s, 1); mpfr_set(r, f, MPFR_RNDN);
        return;
    }
    mpfr_hypot(r, f, g, MPFR_RNDN);
    mpfr_div(c, f, r, MPFR_RNDN);
    mpfr_div(s, g, r, MPFR_RNDN);
}

/* x' = c x + s y,  y' = c y - s x  (t1, t2 scratch at x's precision). */
static void svdb_rot(mpfr_t x, mpfr_t y, const mpfr_t c, const mpfr_t s,
                     mpfr_t t1, mpfr_t t2)
{
    mpfr_mul(t1, s, y, MPFR_RNDN);
    mpfr_fma(t1, c, x, t1, MPFR_RNDN);
    mpfr_mul(t2, s, x, MPFR_RNDN);
    mpfr_fms(y, c, y, t2, MPFR_RNDN);
    mpfr_swap(x, t1);
}

/* Singular values of the upper bidiagonal (d[0..p-1], e[0..p-2]), left in
 * d as |values|.  False if the iteration does not converge. */
static bool svdb_bidiag_qr(mpfr_t* d, mpfr_t* e, int p, mpfr_prec_t w)
{
    mpfr_t c, s, r, f, g, t1, t2, tol, abstol, mu, a, b, x;
    mpfr_init2(c, w); mpfr_init2(s, w); mpfr_init2(r, w); mpfr_init2(f, w);
    mpfr_init2(g, w); mpfr_init2(t1, w); mpfr_init2(t2, w); mpfr_init2(tol, w);
    mpfr_init2(abstol, w); mpfr_init2(mu, w); mpfr_init2(a, w); mpfr_init2(b, w);
    mpfr_init2(x, w);

    /* Relative deflation threshold 4 2^-w; the absolute one, below which
     * an entry is noise of the backward-stable reduction anyway,
     * p 2^-w ||B|| with ||B|| <= max_i |d_i| + |e_i|. */
    mpfr_set_ui_2exp(tol, 1, 2 - (long)w, MPFR_RNDN);
    mpfr_set_zero(abstol, 1);
    for (int i = 0; i < p; i++) {
        mpfr_abs(t1, d[i], MPFR_RNDN);
        if (i + 1 < p) { mpfr_abs(t2, e[i], MPFR_RNDN); mpfr_add(t1, t1, t2, MPFR_RNDN); }
        if (mpfr_cmp(t1, abstol) > 0) mpfr_set(abstol, t1, MPFR_RNDN);
    }
    mpfr_mul_ui(abstol, abstol, (unsigned long)p, MPFR_RNDN);
    mpfr_mul_2si(abstol, abstol, -(long)w, MPFR_RNDN);

    bool ok = true;
    long iter = 0, maxit = 60L * p + 100;
    int hi = p - 1;
    while (hi > 0) {
        for (int i = 0; i < hi; i++) {
            if (mpfr_zero_p(e[i])) continue;
            mpfr_abs(t1, d[i], MPFR_RNDN);
            mpfr_abs(t2, d[i + 1], MPFR_RNDN);
            mpfr_add(t1, t1, t2, MPFR_RNDN);
            mpfr_mul(t1, t1, tol, MPFR_RNDN);
            if (mpfr_cmpabs(e[i], t1) <= 0 || mpfr_cmpabs(e[i], abstol) <= 0)
                mpfr_set_zero(e[i], 1);
        }
        while (hi > 0 && mpfr_zero_p(e[hi - 1])) hi--;
        if (hi == 0) break;
        int lo = hi - 1;
        while (lo > 0 && !mpfr_zero_p(e[lo - 1])) lo--;
        if (++iter > maxit) { ok = false; break; }

        /* A negligible diagonal entry splits the block. */
        int z = -1;
        for (int i = lo; i <= hi; i++)
            if (mpfr_cmpabs(d[i], abstol) <= 0) { z = i; break; }
        if (z >= 0) {
            mpfr_set_zero(d[z], 1);
            if (z < hi) {
                /* Rotate rows j and z (j = z+1..hi) to chase x = B(z, j)
                 * off the right end: (d_j, x) -> (r, 0), and the row-z
                 * entry of column j+1 becomes -s e_j. */
                mpfr_set(x, e[z], MPFR_RNDN);
                mpfr_set_zero(e[z], 1);
                for (int j = z + 1; j <= hi; j++) {
                    svdb_givens(c, s, r, d[j], x);
                    mpfr_set(d[j], r, MPFR_RNDN);
                    if (j < hi) {
                        mpfr_mul(x, s, e[j], MPFR_RNDN);
                        mpfr_neg(x, x, MPFR_RNDN);
                        mpfr_mul(e[j], c, e[j], MPFR_RNDN);
                    }
                }
            } else {
                /* Rotate columns j and hi (j = hi-1..lo) to chase
                 * x = B(j, hi) up: (d_j, x) -> (r, 0), and the column-hi
                 * entry of row j-1 becomes -s e_(j-1). */
                mpfr_set(x, e[hi - 1], MPFR_RNDN);
                mpfr_set_zero(e[hi - 1], 1);
                for (int j = hi - 1; j >= lo; j--) {
                    svdb_givens(c, s, r, d[j], x);
                    mpfr_set(d[j], r, MPFR_RNDN);
                    if (j > lo) {
                        mpfr_mul(x, s, e[j - 1], MPFR_RNDN);
                        mpfr_neg(x, x, MPFR_RNDN);
                        mpfr_mul(e[j - 1], c, e[j - 1], MPFR_RNDN);
                    }
                }
            }
            continue;
        }

        /* Wilkinson shift: the eigenvalue of the trailing 2 x 2
         * [[a, f], [f, b]] of B^T B nearer b. */
        int m = hi - 1;
        mpfr_sqr(a, d[m], MPFR_RNDN);
        if (m > lo) { mpfr_sqr(t1, e[m - 1], MPFR_RNDN); mpfr_add(a, a, t1, MPFR_RNDN); }
        mpfr_sqr(b, d[hi], MPFR_RNDN);
        mpfr_sqr(t1, e[m], MPFR_RNDN);
        mpfr_add(b, b, t1, MPFR_RNDN);
        mpfr_mul(f, d[m], e[m], MPFR_RNDN);
        mpfr_sub(t1, a, b, MPFR_RNDN);
        mpfr_div_2ui(t1, t1, 1, MPFR_RNDN);              /* delta */
        mpfr_hypot(t2, t1, f, MPFR_RNDN);
        if (mpfr_sgn(t1) < 0) mpfr_neg(t2, t2, MPFR_RNDN);
        mpfr_add(t1, t1, t2, MPFR_RNDN);
        if (mpfr_zero_p(t1)) mpfr_set(mu, b, MPFR_RNDN);
        else {
            mpfr_sqr(t2, f, MPFR_RNDN);
            mpfr_div(t2, t2, t1, MPFR_RNDN);
            mpfr_sub(mu, b, t2, MPFR_RNDN);
        }

        /* The bulge chase.  (f, g) is the pair the next rotation zeroes. */
        mpfr_sqr(f, d[lo], MPFR_RNDN);
        mpfr_sub(f, f, mu, MPFR_RNDN);
        mpfr_mul(g, d[lo], e[lo], MPFR_RNDN);
        for (int k = lo; k < hi; k++) {
            /* columns k, k+1 */
            svdb_givens(c, s, r, f, g);
            if (k > lo) mpfr_set(e[k - 1], r, MPFR_RNDN);
            svdb_rot(d[k], e[k], c, s, t1, t2);
            mpfr_mul(g, s, d[k + 1], MPFR_RNDN);          /* bulge B(k+1, k) */
            mpfr_mul(d[k + 1], c, d[k + 1], MPFR_RNDN);
            /* rows k, k+1 */
            svdb_givens(c, s, r, d[k], g);
            mpfr_set(d[k], r, MPFR_RNDN);
            svdb_rot(e[k], d[k + 1], c, s, t1, t2);
            mpfr_set(f, e[k], MPFR_RNDN);
            if (k + 1 < hi) {
                mpfr_mul(g, s, e[k + 1], MPFR_RNDN);      /* bulge B(k, k+2) */
                mpfr_mul(e[k + 1], c, e[k + 1], MPFR_RNDN);
            }
        }
    }
    for (int i = 0; i < p; i++) mpfr_abs(d[i], d[i], MPFR_RNDN);

    mpfr_clear(c); mpfr_clear(s); mpfr_clear(r); mpfr_clear(f); mpfr_clear(g);
    mpfr_clear(t1); mpfr_clear(t2); mpfr_clear(tol); mpfr_clear(abstol);
    mpfr_clear(mu); mpfr_clear(a); mpfr_clear(b); mpfr_clear(x);
    return ok;
}

/* One pass at working precision w: the p singular values (unsorted) of the
 * column-major n x p real matrix Acm, p <= n, into sigma.  False on
 * non-convergence. */
static bool svdm_values_bidiag_at(const mpfr_t* Acm, int n, int p, mpfr_prec_t w,
                                  mpfr_t* sigma)
{
    slong prec = (slong)w;
    size_t N = (size_t)n, P = (size_t)p;
    arb_ptr A = _arb_vec_init((slong)(N * P));
    arb_ptr u = _arb_vec_init((slong)N);
    arb_ptr wv = _arb_vec_init((slong)P);
    for (size_t i = 0; i < N; i++)
        for (size_t j = 0; j < P; j++)
            arf_set_mpfr(arb_midref(A + i * P + j), Acm[i + j * N]);
    mpfr_t* d = svdm_array_alloc(P, w);
    mpfr_t* e = svdm_array_alloc(P, w);
    arf_t alpha;
    arf_init(alpha);
    SvdBi c;
    c.A = A; c.u = u; c.wv = wv; c.n = N; c.p = P; c.prec = prec;
    double limbs = (double)w / 64.0;
    for (size_t k = 0; k < P; k++) {
        c.k = k;
        /* left reflector: zero A[k+1.., k] */
        if (svdb_house(A + k * P + k, (slong)P, N - k, u + k, alpha, prec) && k + 1 < P) {
            double work = (double)(N - k) * (double)(P - k - 1) * limbs * limbs;
            mp_parallel_for(svdb_left_dot, &c, (int64_t)(P - k - 1), work);
            mp_parallel_for(svdb_left_upd, &c, (int64_t)(N - k), work);
        }
        arf_get_mpfr(d[k], alpha, MPFR_RNDN);
        if (k + 1 >= P) break;
        /* right reflector: zero A[k, k+2..] */
        if (k + 2 < P) {
            if (svdb_house(A + k * P + k + 1, 1, P - k - 1, u + k + 1, alpha, prec))
                mp_parallel_for(svdb_right, &c, (int64_t)(N - k - 1),
                                2.0 * (double)(N - k - 1) * (double)(P - k - 1) * limbs * limbs);
            arf_get_mpfr(e[k], alpha, MPFR_RNDN);
        } else {
            arf_get_mpfr(e[k], arb_midref(A + k * P + k + 1), MPFR_RNDN);
        }
    }
    mpfr_set_zero(e[P - 1], 1);
    arf_clear(alpha);
    _arb_vec_clear(A, (slong)(N * P));
    _arb_vec_clear(u, (slong)N);
    _arb_vec_clear(wv, (slong)P);
    bool ok = svdb_bidiag_qr(d, e, p, w);
    if (ok)
        for (int i = 0; i < p; i++) mpfr_set(sigma[i], d[i], MPFR_RNDN);
    svdm_array_free(d, P);
    svdm_array_free(e, P);
    return ok;
}

/* The singular values of the exact (binary) input to `bits` relative
 * accuracy.  A backward-stable pass at w bits leaves sigma_i with an
 * absolute error ~ 2^-w sigma_max, i.e. log2(sigma_max / sigma_i) bits of
 * its relative accuracy are lost -- the 1.2e-42 smallest singular value of
 * a 50 x 50 equispaced Vandermonde matrix at 100 digits came out with 57.
 * So when the spread of the values that survive the default tolerance
 * (8 n 2^-bits sigma_max, the caller's floor) exceeds the guard, the pass
 * is repeated once with the spread added to the working precision. */
static bool svdm_values_bidiag(const mpfr_t* Acm, int n, int p, mpfr_prec_t bits,
                               mpfr_t* sigma)
{
    mpfr_prec_t w = bits + SVDB_GUARD;
    if (!svdm_values_bidiag_at(Acm, n, p, w, sigma)) return false;
    mpfr_t smax, floor_, smin;
    mpfr_init2(smax, 64); mpfr_init2(floor_, 64); mpfr_init2(smin, 64);
    mpfr_set_zero(smax, 1);
    for (int i = 0; i < p; i++)
        if (mpfr_cmp(sigma[i], smax) > 0) mpfr_set(smax, sigma[i], MPFR_RNDN);
    mpfr_mul_ui(floor_, smax, 8UL * (unsigned long)n, MPFR_RNDN);
    mpfr_mul_2si(floor_, floor_, -(long)bits, MPFR_RNDN);
    mpfr_set(smin, smax, MPFR_RNDN);
    for (int i = 0; i < p; i++)
        if (mpfr_cmp(sigma[i], floor_) > 0 && mpfr_cmp(sigma[i], smin) < 0)
            mpfr_set(smin, sigma[i], MPFR_RNDN);
    long spread = 0;
    if (!mpfr_zero_p(smin))
        spread = (long)mpfr_get_exp(smax) - (long)mpfr_get_exp(smin) + 1;
    mpfr_clear(smax); mpfr_clear(floor_); mpfr_clear(smin);
    if (spread > SVDB_GUARD / 2)
        return svdm_values_bidiag_at(Acm, n, p, w + (mpfr_prec_t)spread, sigma);
    return true;
}
#endif /* USE_FLINT */

/* ------------------------------------------------------------------ *
 *  Singular values only (SingularValueList, Norm[m, 2]).              *
 *                                                                     *
 *  Real input with FLINT: bidiagonalisation + bidiagonal QR           *
 *  (svdm_values_bidiag).  Otherwise, and as the fallback should that  *
 *  iteration not converge, the same one-sided Jacobi as the full      *
 *  decomposition, with V not accumulated and U never formed: after    *
 *  the sweeps the columns of A are mutually orthogonal and their      *
 *  2-norms are the singular values.                                   *
 *  Requires p <= n (the caller hands a wide matrix in transposed, the *
 *  singular values being those of the conjugate transpose).           *
 * ------------------------------------------------------------------ */
mpfr_t* svd_mpfr_values(Expr* m, int n, int p, mpfr_prec_t bits)
{
    if (p > n || p < 1) return NULL;
    mpfr_t* sigma = NULL;
    size_t total = (size_t)n * (size_t)p;
    mpfr_t acc, t2;
    if (svdm_matrix_has_complex(m, n, p)) {
        mpfr_t *A_re = NULL, *A_im = NULL;
        if (!svdm_load_complex(m, n, p, bits, &A_re, &A_im)) return NULL;
        (void)svdm_jacobi_sweeps_complex(A_re, A_im, n, p, NULL, NULL, bits);
        sigma = svdm_array_alloc((size_t)p, bits);
        mpfr_init2(acc, bits); mpfr_init2(t2, bits);
        for (int j = 0; j < p; j++) {
            mpfr_set_zero(acc, 1);
            for (int i = 0; i < n; i++) {
                size_t off = (size_t)i + (size_t)j * (size_t)n;
                mpfr_sqr(t2, A_re[off], MPFR_RNDN); mpfr_add(acc, acc, t2, MPFR_RNDN);
                mpfr_sqr(t2, A_im[off], MPFR_RNDN); mpfr_add(acc, acc, t2, MPFR_RNDN);
            }
            mpfr_sqrt(sigma[j], acc, MPFR_RNDN);
        }
        svdm_array_free(A_re, total);
        svdm_array_free(A_im, total);
    } else {
        mpfr_t* A = NULL;
        if (!svdm_load_real(m, n, p, bits, &A)) return NULL;
        sigma = svdm_array_alloc((size_t)p, bits);
        mpfr_init2(acc, bits); mpfr_init2(t2, bits);
#ifdef USE_FLINT
        if (svdm_values_bidiag((const mpfr_t*)A, n, p, bits, sigma)) {
            svdm_array_free(A, total);
            goto sorted;
        }
#endif
        (void)svdm_jacobi_sweeps(A, n, p, NULL, bits);
        for (int j = 0; j < p; j++) {
            mpfr_set_zero(acc, 1);
            for (int i = 0; i < n; i++) {
                mpfr_sqr(t2, A[(size_t)i + (size_t)j * (size_t)n], MPFR_RNDN);
                mpfr_add(acc, acc, t2, MPFR_RNDN);
            }
            mpfr_sqrt(sigma[j], acc, MPFR_RNDN);
        }
        svdm_array_free(A, total);
    }
#ifdef USE_FLINT
sorted:
#endif
    mpfr_clear(acc); mpfr_clear(t2);
    /* descending (insertion sort: p is small) */
    for (int i = 1; i < p; i++)
        for (int j = i; j > 0 && mpfr_cmp(sigma[j - 1], sigma[j]) < 0; j--)
            mpfr_swap(sigma[j - 1], sigma[j]);
    return sigma;
}

void svd_mpfr_values_free(mpfr_t* sigma, int p)
{
    svdm_array_free(sigma, (size_t)p);
}

#endif /* USE_MPFR */
