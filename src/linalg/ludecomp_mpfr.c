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
#include "message.h"   /* mth_message: Quiet/Check funnel */
#include "core.h"      /* tc_deadline_passed / tc_check_deadline */
#include "term/term_interrupt.h"   /* once-per-session warnings */

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

Expr* lu_mpfr_dispatch(Expr* m, int rows, int cols)
{
    (void)m; (void)rows; (void)cols;
    return NULL;
}

Expr* mpfr_det_dispatch(Expr* m, int n)
{
    (void)m; (void)n;
    return NULL;
}

MpfrLU* mpfr_lu_factor(const Expr* m, int n, long bits, bool* singular)
{
    (void)m; (void)n; (void)bits; *singular = false;
    return NULL;
}

MpfrLU* mpfr_lu_from_parts(const Expr* lu, const int* perm, int n, long bits)
{
    (void)lu; (void)perm; (void)n; (void)bits;
    return NULL;
}

Expr* mpfr_lu_solve(const MpfrLU* f, const Expr* b) { (void)f; (void)b; return NULL; }
long  mpfr_lu_bits(const MpfrLU* f) { (void)f; return 0; }
void  mpfr_lu_free(MpfrLU* f) { (void)f; }

#else /* USE_MPFR */

/* ------------------------------------------------------------------ *
 *  Allocation helpers.  Duplicated from qrdecomp_mpfr.c so this        *
 *  translation unit stays self-contained.                              *
 * ------------------------------------------------------------------ */
/* NULL also when a TimeConstrained deadline passes while the entries are
 * initialised (polled every 2^16): an n = 10^4 matrix is 10^8 limb blocks,
 * tens of seconds of allocation with no evaluator step in between. */
static mpfr_t* lum_array_alloc(size_t count, mpfr_prec_t bits)
{
    if (count == 0) return NULL;
    mpfr_t* a = (mpfr_t*)malloc(sizeof(mpfr_t) * count);
    if (!a) return NULL;
    for (size_t i = 0; i < count; i++) {
        if ((i & 0xffff) == 0 && i > 0 && tc_deadline_passed()) {
            while (i-- > 0) mpfr_clear(a[i]);
            free(a);
            return NULL;
        }
        mpfr_init2(a[i], bits);
    }
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

/* Load the rows x cols input matrix into freshly-allocated row-major
 * MPFR arrays.  Row-major (rather than column-major) is convenient
 * because the Doolittle kernel always accesses whole rows for both
 * pivoting swaps and the elimination's outer loop. */
static bool lum_load_matrix(Expr* m, int rows, int cols, mpfr_prec_t bits,
                            mpfr_t** out_re, mpfr_t** out_im,
                            bool* out_is_complex)
{
    if (m->type != EXPR_FUNCTION) return false;

    bool is_complex = false;
    for (int i = 0; i < rows && !is_complex; i++) {
        Expr* row = m->data.function.args[i];
        if (!row || row->type != EXPR_FUNCTION) return false;
        for (int j = 0; j < cols && !is_complex; j++) {
            if (lum_leaf_is_complex(row->data.function.args[j])) {
                is_complex = true;
            }
        }
    }

    size_t total = (size_t)rows * (size_t)cols;
    mpfr_t* A_re = lum_array_alloc(total, bits);
    mpfr_t* A_im = is_complex ? lum_array_alloc(total, bits) : NULL;
    if (!A_re || (is_complex && !A_im)) {        /* no memory, or a timeout */
        lum_array_free(A_re, A_re ? total : 0);
        lum_array_free(A_im, A_im ? total : 0);
        return false;
    }

    mpfr_t tmp_im;
    mpfr_init2(tmp_im, bits);
    for (int i = 0; i < rows; i++) {
        Expr* row = m->data.function.args[i];
        if (tc_deadline_passed()) {
            mpfr_clear(tmp_im);
            lum_array_free(A_re, total);
            lum_array_free(A_im, A_im ? total : 0);
            return false;
        }
        for (int j = 0; j < cols; j++) {
            size_t off = (size_t)i * (size_t)cols + (size_t)j;
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

/* Swap whole rows in a row-major matrix of width `cols`. */
static void lum_swap_rows(mpfr_t* A_re, mpfr_t* A_im, int cols,
                          int a, int b, bool is_complex)
{
    if (a == b) return;
    for (int j = 0; j < cols; j++) {
        mpfr_swap(A_re[(size_t)a * cols + j],
                  A_re[(size_t)b * cols + j]);
        if (is_complex) mpfr_swap(A_im[(size_t)a * cols + j],
                                   A_im[(size_t)b * cols + j]);
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
 *  Doolittle elimination.  Overwrites A (rows x cols, row-major) with *
 *  the combined LU; perm starts as [1..min(rows, cols)] and is        *
 *  updated to the pivoting permutation.  Returns true on success.     *
 *  *out_singular is set to true iff a zero pivot was encountered      *
 *  (factorisation continues; matches Mathematica behaviour).          *
 * ------------------------------------------------------------------ */
#ifdef USE_FLINT
/* ------------------------------------------------------------------ *
 *  Real matrices: Crout's ordering of the same elimination, every     *
 *  entry one fused dot product.                                       *
 *                                                                       *
 *  Step k first finishes column k of the active rows,                  *
 *      a_ik -= sum_{j<k} l_ij u_jk          (i >= k),                   *
 *  picks the pivot among them (largest |a_ik|, as above), swaps it    *
 *  into row k, finishes row k of U,                                    *
 *      a_kj -= sum_{i<k} l_ki u_ij          (j > k),                    *
 *  and scales column k of L by 1/u_kk.  In exact arithmetic these are *
 *  the values, pivots and permutation of the right-looking loop; here *
 *  each sum is one arb_approx_dot rounded once (no per-term rounding, *
 *  no per-term temporaries), and the independent dot products of each *
 *  phase are split over threads (mp_parallel_for).  At 500 digits a   *
 *  256 x 256 factorisation went from 3.5 s to ~0.5 s.                  *
 * ------------------------------------------------------------------ */
typedef struct {
    arb_ptr F;          /* rows x cols, row-major */
    int cols;
    int k;
    slong prec;
    arb_ptr tmp[MP_PAR_MAX];
} LumCrout;

/* a_ik -= L[i, :k] . U[:k, k] for i in [k + lo, k + hi) */
static void lum_crout_col(void* p, int64_t lo, int64_t hi, int chunk)
{
    LumCrout* c = (LumCrout*)p;
    arb_ptr t = c->tmp[chunk];
    for (int64_t i = c->k + lo; i < c->k + hi; i++) {
        arb_ptr x = c->F + i * c->cols + c->k;
        arb_approx_dot(t, x, 1, c->F + i * c->cols, 1, c->F + c->k, c->cols, c->k, c->prec);
        arb_swap(x, t);
    }
}

/* a_kj -= L[k, :k] . U[:k, j] for j in [k + 1 + lo, k + 1 + hi) */
static void lum_crout_row(void* p, int64_t lo, int64_t hi, int chunk)
{
    LumCrout* c = (LumCrout*)p;
    arb_ptr t = c->tmp[chunk];
    for (int64_t j = c->k + 1 + lo; j < c->k + 1 + hi; j++) {
        arb_ptr x = c->F + (int64_t)c->k * c->cols + j;
        arb_approx_dot(t, x, 1, c->F + (int64_t)c->k * c->cols, 1, c->F + j, c->cols, c->k, c->prec);
        arb_swap(x, t);
    }
}

static bool lum_factor_crout(mpfr_t* A, int rows, int cols, mpfr_prec_t bits,
                             int* perm, bool* out_singular)
{
    int steps = rows < cols ? rows : cols;
    size_t total = (size_t)rows * (size_t)cols;
    LumCrout c;
    c.F = _arb_vec_init((slong)total);
    c.cols = cols;
    c.prec = (slong)bits;
    for (int i = 0; i < MP_PAR_MAX; i++) { c.tmp[i] = (arb_ptr)flint_malloc(sizeof(arb_struct)); arb_init(c.tmp[i]); }
    for (size_t i = 0; i < total; i++) arf_set_mpfr(arb_midref(c.F + i), A[i]);
    double limbs = (double)bits / 64.0;
    *out_singular = false;
    bool timed_out = false;
    for (int k = 0; k < steps; k++) {
        /* the TimeConstrained poll (see lum_factor) */
        if (tc_deadline_passed()) { timed_out = true; break; }
        c.k = k;
        if (k > 0)
            mp_parallel_for(lum_crout_col, &c, rows - k,
                            (double)(rows - k) * (double)k * limbs * limbs);
        int piv = k;
        for (int i = k + 1; i < rows; i++)
            if (arf_cmpabs(arb_midref(c.F + (int64_t)i * cols + k),
                           arb_midref(c.F + (int64_t)piv * cols + k)) > 0)
                piv = i;
        if (arf_is_zero(arb_midref(c.F + (int64_t)piv * cols + k))) {
            /* a zero column: U[k, k] = 0, the L entries below are 0 too */
            *out_singular = true;
            if (k > 0 && k + 1 < cols)
                mp_parallel_for(lum_crout_row, &c, cols - k - 1,
                                (double)(cols - k - 1) * (double)k * limbs * limbs);
            continue;
        }
        if (piv != k) {
            _arb_vec_swap(c.F + (int64_t)k * cols, c.F + (int64_t)piv * cols, cols);
            int tp = perm[k]; perm[k] = perm[piv]; perm[piv] = tp;
        }
        if (k > 0 && k + 1 < cols)
            mp_parallel_for(lum_crout_row, &c, cols - k - 1,
                            (double)(cols - k - 1) * (double)k * limbs * limbs);
        arb_ptr d = c.F + (int64_t)k * cols + k;
        for (int i = k + 1; i < rows; i++) {
            arb_ptr x = c.F + (int64_t)i * cols + k;
            arb_div(x, x, d, c.prec);
            mag_zero(arb_radref(x));
        }
    }
    if (!timed_out)
        for (size_t i = 0; i < total; i++) arf_get_mpfr(A[i], arb_midref(c.F + i), MPFR_RNDN);
    for (int i = 0; i < MP_PAR_MAX; i++) { arb_clear(c.tmp[i]); flint_free(c.tmp[i]); }
    _arb_vec_clear(c.F, (slong)total);
    return !timed_out;
}
#endif /* USE_FLINT */

/* Returns false -- the matrix half-factored, to be discarded -- when an
 * active TimeConstrained deadline passes.  The factorisation is the one long
 * C loop of an arbitrary-precision LinearSolve / LUDecomposition / Det with
 * no evaluator step inside it, and most of its time is spent in worker
 * threads (mp_parallel_for), inside a tc_async_defer region where the
 * SIGPROF timeout may not unwind; without this poll a 10000 x 10000 factor
 * at 60 digits ran > 12 min past a 60 s cap.  The step loop polls the
 * deadline (non-jumping) on the main thread; callers free their arrays and
 * then take the timeout jump with tc_check_deadline(). */
static bool lum_factor(mpfr_t* A_re, mpfr_t* A_im, int rows, int cols,
                        mpfr_prec_t bits, bool is_complex,
                        int* perm, bool* out_singular)
{
#ifdef USE_FLINT
    if (!is_complex)
        return lum_factor_crout(A_re, rows, cols, bits, perm, out_singular);
#endif
    int steps = (rows < cols) ? rows : cols;
    mpfr_t mag, best;
    mpfr_init2(mag,  bits);
    mpfr_init2(best, bits);
    mpfr_t prod;
    mpfr_init2(prod, bits);
    *out_singular = false;
    bool timed_out = false;

    for (int k = 0; k < steps; k++) {
        if (tc_deadline_passed()) { timed_out = true; break; }
        /* Pivot: largest |A[i, k]| over i in [k, rows). */
        int pivot_row = k;
        lum_abs(A_re, A_im, (size_t)k * cols + k, is_complex, best);
        for (int i = k + 1; i < rows; i++) {
            lum_abs(A_re, A_im, (size_t)i * cols + k, is_complex, mag);
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
            lum_swap_rows(A_re, A_im, cols, k, pivot_row, is_complex);
            int tp = perm[k]; perm[k] = perm[pivot_row]; perm[pivot_row] = tp;
        }

        /* Eliminate. */
        for (int i = k + 1; i < rows; i++) {
            size_t ik = (size_t)i * cols + k;
            size_t kk = (size_t)k * cols + k;
            if (is_complex) {
                lum_cdiv(&A_re[ik], &A_im[ik], A_re[kk], A_im[kk], bits);
            } else {
                mpfr_div(A_re[ik], A_re[ik], A_re[kk], MPFR_RNDN);
            }
            for (int j = k + 1; j < cols; j++) {
                size_t ij = (size_t)i * cols + j;
                size_t kj = (size_t)k * cols + j;
                if (is_complex) {
                    lum_csub_mul(&A_re[ij], &A_im[ij],
                                  A_re[ik], A_im[ik],
                                  A_re[kj], A_im[kj], bits);
                } else {
                    mpfr_mul(prod, A_re[ik], A_re[kj], MPFR_RNDN);
                    mpfr_sub(A_re[ij], A_re[ij], prod, MPFR_RNDN);
                }
            }
        }
    }
    mpfr_clear(prod);
    mpfr_clear(mag);
    mpfr_clear(best);
    return !timed_out;
}

/* ------------------------------------------------------------------ *
 *  L-infinity norm of a row-major rows x cols matrix.                  *
 *      max_i sum_j |A[i, j]|                                            *
 * ------------------------------------------------------------------ */
static void lum_norm_inf(const mpfr_t* A_re, const mpfr_t* A_im,
                          int rows, int cols,
                          bool is_complex, mpfr_prec_t bits, mpfr_t out)
{
    mpfr_t row_sum, mag;
    mpfr_init2(row_sum, bits);
    mpfr_init2(mag,     bits);
    mpfr_set_zero(out, 1);
    for (int i = 0; i < rows; i++) {
        mpfr_set_zero(row_sum, 1);
        for (int j = 0; j < cols; j++) {
            lum_abs(A_re, A_im, (size_t)i * cols + j, is_complex, mag);
            mpfr_add(row_sum, row_sum, mag, MPFR_RNDN);
        }
        if (mpfr_cmp(row_sum, out) > 0) mpfr_set(out, row_sum, MPFR_RNDN);
    }
    mpfr_clear(row_sum);
    mpfr_clear(mag);
}

/* In-place solve of LU * x = rhs, where LU packs the Doolittle L
 * (strict-lower, unit diag implicit) and U (upper).  rhs is read,
 * then overwritten with the solution x.  Used by both the explicit
 * inverse and the Hager-Higham one-norm estimator. */
static void lum_solve_LU_inplace(const mpfr_t* LU_re, const mpfr_t* LU_im,
                                  int n, bool is_complex,
                                  mpfr_t* rhs_re, mpfr_t* rhs_im,
                                  mpfr_prec_t bits)
{
    mpfr_t prod;
    mpfr_init2(prod, bits);
    /* Forward substitution: L y = b. */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            size_t ij = (size_t)i * n + j;
            if (is_complex) {
                lum_csub_mul(&rhs_re[i], &rhs_im[i],
                              LU_re[ij], LU_im[ij],
                              rhs_re[j], rhs_im[j], bits);
            } else {
                mpfr_mul(prod, LU_re[ij], rhs_re[j], MPFR_RNDN);
                mpfr_sub(rhs_re[i], rhs_re[i], prod, MPFR_RNDN);
            }
        }
    }
    /* Back substitution: U x = y. */
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i + 1; j < n; j++) {
            size_t ij = (size_t)i * n + j;
            if (is_complex) {
                lum_csub_mul(&rhs_re[i], &rhs_im[i],
                              LU_re[ij], LU_im[ij],
                              rhs_re[j], rhs_im[j], bits);
            } else {
                mpfr_mul(prod, LU_re[ij], rhs_re[j], MPFR_RNDN);
                mpfr_sub(rhs_re[i], rhs_re[i], prod, MPFR_RNDN);
            }
        }
        size_t ii = (size_t)i * n + i;
        if (is_complex) {
            lum_cdiv(&rhs_re[i], &rhs_im[i], LU_re[ii], LU_im[ii], bits);
        } else {
            mpfr_div(rhs_re[i], rhs_re[i], LU_re[ii], MPFR_RNDN);
        }
    }
    mpfr_clear(prod);
}

/* In-place solve of (LU)^T * x = rhs.  (LU)^T = U^T * L^T, so we
 * forward-substitute through U^T (lower triangular when transposed)
 * then back-substitute through L^T (upper triangular when transposed,
 * with unit diagonal).  Used by Hager-Higham to apply A^{-T}. */
static void lum_solve_LU_T_inplace(const mpfr_t* LU_re, const mpfr_t* LU_im,
                                    int n, bool is_complex,
                                    mpfr_t* rhs_re, mpfr_t* rhs_im,
                                    mpfr_prec_t bits)
{
    mpfr_t prod;
    mpfr_init2(prod, bits);
    /* Forward sub on U^T: U^T[i, j] = LU[j, i] for j <= i; diagonal
     * is U[i, i] = LU[i, i] (NOT unit). */
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < i; j++) {
            /* U^T[i, j] = U[j, i] = LU[j, i]. */
            size_t ji = (size_t)j * n + i;
            if (is_complex) {
                lum_csub_mul(&rhs_re[i], &rhs_im[i],
                              LU_re[ji], LU_im[ji],
                              rhs_re[j], rhs_im[j], bits);
            } else {
                mpfr_mul(prod, LU_re[ji], rhs_re[j], MPFR_RNDN);
                mpfr_sub(rhs_re[i], rhs_re[i], prod, MPFR_RNDN);
            }
        }
        size_t ii = (size_t)i * n + i;
        if (is_complex) {
            lum_cdiv(&rhs_re[i], &rhs_im[i], LU_re[ii], LU_im[ii], bits);
        } else {
            mpfr_div(rhs_re[i], rhs_re[i], LU_re[ii], MPFR_RNDN);
        }
    }
    /* Back sub on L^T: L^T[i, j] = L[j, i] = LU[j, i] for j > i;
     * diagonal is L[i, i] = 1 (unit). */
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i + 1; j < n; j++) {
            /* L^T[i, j] = L[j, i] = LU[j, i]. */
            size_t ji = (size_t)j * n + i;
            if (is_complex) {
                lum_csub_mul(&rhs_re[i], &rhs_im[i],
                              LU_re[ji], LU_im[ji],
                              rhs_re[j], rhs_im[j], bits);
            } else {
                mpfr_mul(prod, LU_re[ji], rhs_re[j], MPFR_RNDN);
                mpfr_sub(rhs_re[i], rhs_re[i], prod, MPFR_RNDN);
            }
        }
        /* L diagonal is 1; no division. */
    }
    mpfr_clear(prod);
}

/* Convenience wrapper: solve LU * x = e_col.  Used by the legacy
 * explicit-inverse condition estimator (kept for the singular check). */
static void lum_solve_LU_canonical(const mpfr_t* LU_re, const mpfr_t* LU_im,
                                    int n, bool is_complex, int col,
                                    mpfr_t* out_re, mpfr_t* out_im,
                                    mpfr_prec_t bits)
{
    for (int i = 0; i < n; i++) {
        mpfr_set_zero(out_re[i], 1);
        if (is_complex) mpfr_set_zero(out_im[i], 1);
    }
    mpfr_set_si(out_re[col], 1, MPFR_RNDN);
    lum_solve_LU_inplace(LU_re, LU_im, n, is_complex,
                          out_re, out_im, bits);
}

/* Apply the row permutation P to a vector in place:
 *    v'[i] := v[perm[i] - 1]
 * Used between LU-solve stages in the Hager-Higham estimator. */
static void lum_apply_perm(const int* perm, int n, mpfr_t* v_re,
                            mpfr_t* v_im, bool is_complex)
{
    mpfr_t* tmp_re = lum_array_alloc((size_t)n, mpfr_get_prec(v_re[0]));
    mpfr_t* tmp_im = is_complex
                   ? lum_array_alloc((size_t)n, mpfr_get_prec(v_re[0]))
                   : NULL;
    for (int i = 0; i < n; i++) {
        int src = perm[i] - 1;
        mpfr_set(tmp_re[i], v_re[src], MPFR_RNDN);
        if (is_complex) mpfr_set(tmp_im[i], v_im[src], MPFR_RNDN);
    }
    for (int i = 0; i < n; i++) {
        mpfr_set(v_re[i], tmp_re[i], MPFR_RNDN);
        if (is_complex) mpfr_set(v_im[i], tmp_im[i], MPFR_RNDN);
    }
    lum_array_free(tmp_re, (size_t)n);
    if (is_complex) lum_array_free(tmp_im, (size_t)n);
}

/* Apply the inverse row permutation P^T:
 *    v'[perm[i] - 1] := v[i]
 * Used to finish A^{-T} v after the L^T / U^T solves. */
static void lum_apply_perm_T(const int* perm, int n, mpfr_t* v_re,
                              mpfr_t* v_im, bool is_complex)
{
    mpfr_t* tmp_re = lum_array_alloc((size_t)n, mpfr_get_prec(v_re[0]));
    mpfr_t* tmp_im = is_complex
                   ? lum_array_alloc((size_t)n, mpfr_get_prec(v_re[0]))
                   : NULL;
    for (int i = 0; i < n; i++) {
        int dst = perm[i] - 1;
        mpfr_set(tmp_re[dst], v_re[i], MPFR_RNDN);
        if (is_complex) mpfr_set(tmp_im[dst], v_im[i], MPFR_RNDN);
    }
    for (int i = 0; i < n; i++) {
        mpfr_set(v_re[i], tmp_re[i], MPFR_RNDN);
        if (is_complex) mpfr_set(v_im[i], tmp_im[i], MPFR_RNDN);
    }
    lum_array_free(tmp_re, (size_t)n);
    if (is_complex) lum_array_free(tmp_im, (size_t)n);
}

/* ------------------------------------------------------------------ *
 *  Hager-Higham one-norm estimator (real case).                        *
 *                                                                       *
 *  Estimates ||A^{-1}||_inf by computing ||A^{-T}||_1 via the           *
 *  Hager-Higham iteration (Higham, "FORTRAN codes for estimating the   *
 *  one-norm of a real or complex matrix...", 1988).  This is the       *
 *  algorithm LAPACK's `dlacn2` uses inside `dgecon`.                    *
 *                                                                       *
 *  Cost: each iteration is two triangular solves (each O(n^2)).         *
 *  Convergence typically in 2-5 iterations.  Total work O(n^2) vs the  *
 *  O(n^3) of the explicit-inverse estimator.                            *
 *                                                                       *
 *  Algorithm (for B = A^{-T}, estimating ||B||_1):                      *
 *    1. x = (1/n) * e                                                   *
 *    2. for iter = 1..max_iter:                                         *
 *         y = B * x        (= A^{-T} x: solve A^T y = x)                *
 *         gamma = ||y||_1                                               *
 *         xi[k] = sign(y[k])                                            *
 *         z = B^T * xi     (= A^{-1} xi: solve A z = xi)                *
 *         j = argmax_k |z[k]|                                           *
 *         if iter >= 2 and |z[j]| <= dot(z, x): return gamma            *
 *         x = e_j                                                       *
 *    3. return gamma                                                    *
 *                                                                       *
 *  Returns false on singular U.  Output written to *out.                 *
 * ------------------------------------------------------------------ */
static bool lum_inv_norm_hager(const mpfr_t* LU_re, const mpfr_t* LU_im,
                                const int* perm, int n,
                                mpfr_prec_t bits, mpfr_t out)
{
    /* Check for zero U diagonals -- singular if any. */
    for (int i = 0; i < n; i++) {
        size_t ii = (size_t)i * n + i;
        if (mpfr_zero_p(LU_re[ii])) return false;
    }

    enum { MAX_ITER = 5 };
    mpfr_t* x  = lum_array_alloc((size_t)n, bits);
    mpfr_t* xi = lum_array_alloc((size_t)n, bits);
    mpfr_t gamma_cur, gamma_best, dot_zx, abs_zj, t;
    mpfr_init2(gamma_cur,  bits);
    mpfr_init2(gamma_best, bits);
    mpfr_init2(dot_zx,     bits);
    mpfr_init2(abs_zj,     bits);
    mpfr_init2(t,          bits);
    mpfr_set_zero(gamma_best, 1);

    /* x = (1/n) * e. */
    for (int i = 0; i < n; i++) {
        mpfr_set_ui(x[i], 1, MPFR_RNDN);
        mpfr_div_ui(x[i], x[i], (unsigned long)n, MPFR_RNDN);
    }

    int prev_j = -1;
    bool ok_estimate = false;

    for (int iter = 0; iter < MAX_ITER; iter++) {
        /* y = A^{-T} x: solve A^T y = x.  A^T = U^T L^T P^{-T}.
         *    -> w = (U^T)^{-1} x          [solve U^T w = x]
         *    -> z = (L^T)^{-1} w          [solve L^T z = w]
         *    -> y = P^T z                 [unpermute] */
        for (int i = 0; i < n; i++) mpfr_set(xi[i], x[i], MPFR_RNDN);
        lum_solve_LU_T_inplace(LU_re, LU_im, n, false, xi, NULL, bits);
        lum_apply_perm_T(perm, n, xi, NULL, false);
        /* xi now holds y = A^{-T} x. */

        /* gamma = ||y||_1. */
        mpfr_set_zero(gamma_cur, 1);
        for (int i = 0; i < n; i++) {
            mpfr_abs(t, xi[i], MPFR_RNDN);
            mpfr_add(gamma_cur, gamma_cur, t, MPFR_RNDN);
        }
        if (mpfr_cmp(gamma_cur, gamma_best) > 0) {
            mpfr_set(gamma_best, gamma_cur, MPFR_RNDN);
            ok_estimate = true;
        }

        /* xi := sign(y) (sign(0) = 1). */
        for (int i = 0; i < n; i++) {
            if (mpfr_sgn(xi[i]) < 0) mpfr_set_si(xi[i], -1, MPFR_RNDN);
            else                     mpfr_set_si(xi[i],  1, MPFR_RNDN);
        }

        /* z = A^{-1} xi: solve A z = xi.  A = P^{-1} L U.
         *    -> v = P xi             [permute]
         *    -> z = U^{-1} L^{-1} v  [LU solve] */
        lum_apply_perm(perm, n, xi, NULL, false);
        lum_solve_LU_inplace(LU_re, LU_im, n, false, xi, NULL, bits);
        /* xi now holds z = A^{-1} sign(y). */

        /* j = argmax_i |z[i]|. */
        int j = 0;
        mpfr_abs(abs_zj, xi[0], MPFR_RNDN);
        for (int i = 1; i < n; i++) {
            mpfr_abs(t, xi[i], MPFR_RNDN);
            if (mpfr_cmp(t, abs_zj) > 0) {
                j = i;
                mpfr_set(abs_zj, t, MPFR_RNDN);
            }
        }

        /* Convergence: after at least one update, if |z[j]| <= z^T x,
         * iteration cannot improve gamma -- stop. */
        if (iter >= 1) {
            mpfr_set_zero(dot_zx, 1);
            for (int i = 0; i < n; i++) {
                mpfr_mul(t, xi[i], x[i], MPFR_RNDN);
                mpfr_add(dot_zx, dot_zx, t, MPFR_RNDN);
            }
            if (mpfr_cmp(abs_zj, dot_zx) <= 0) break;
        }

        /* Defensive guard: if we cycle on the same basis index, stop. */
        if (j == prev_j) break;
        prev_j = j;

        /* x = e_j. */
        for (int i = 0; i < n; i++) mpfr_set_zero(x[i], 1);
        mpfr_set_ui(x[j], 1, MPFR_RNDN);
    }

    mpfr_set(out, gamma_best, MPFR_RNDN);

    lum_array_free(x,  (size_t)n);
    lum_array_free(xi, (size_t)n);
    mpfr_clear(gamma_cur);
    mpfr_clear(gamma_best);
    mpfr_clear(dot_zx);
    mpfr_clear(abs_zj);
    mpfr_clear(t);
    return ok_estimate;
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
    Expr* z = expr_new_function(expr_new_symbol(SYM_Complex), args, 2);
    free(args);
    return z;
}

static Expr* lum_build_lu(const mpfr_t* LU_re, const mpfr_t* LU_im,
                           int rows, int cols, bool is_complex,
                           mpfr_prec_t bits)
{
    Expr** row_exprs = (Expr**)malloc(sizeof(Expr*) * (size_t)rows);
    mpfr_t zero_im;
    mpfr_init2(zero_im, bits);
    mpfr_set_zero(zero_im, 1);
    for (int i = 0; i < rows; i++) {
        Expr** elems = (Expr**)malloc(sizeof(Expr*) * (size_t)cols);
        for (int j = 0; j < cols; j++) {
            size_t off = (size_t)i * cols + j;
            if (is_complex) {
                elems[j] = lum_make_scalar(LU_re[off], LU_im[off], true);
            } else {
                elems[j] = lum_make_scalar(LU_re[off], zero_im, false);
            }
        }
        row_exprs[i] = expr_new_function(
            expr_new_symbol(SYM_List), elems, (size_t)cols);
        free(elems);
    }
    mpfr_clear(zero_im);
    Expr* lu = expr_new_function(
        expr_new_symbol(SYM_List), row_exprs, (size_t)rows);
    free(row_exprs);
    return lu;
}

static Expr* lum_build_perm(const int* perm, int n)
{
    Expr** elems = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    for (int k = 0; k < n; k++) elems[k] = expr_new_integer(perm[k]);
    Expr* p = expr_new_function(
        expr_new_symbol(SYM_List), elems, (size_t)n);
    free(elems);
    return p;
}

/* ------------------------------------------------------------------ *
 *  Kernel.                                                              *
 *                                                                       *
 *  Supports rectangular input.  The explicit-inverse condition         *
 *  estimator (lum_inv_norm_inf) requires a square matrix; for          *
 *  non-square input we set the c slot to exact Integer 0 (the          *
 *  condition number is undefined for non-square A).                     *
 * ------------------------------------------------------------------ */
Expr* lu_mpfr_dispatch(Expr* m, int rows, int cols)
{
    static uint64_t sing_warn_counter = 0;

    CommonInexactInfo info = common_scan_inexact(m);
    if (!info.has_inexact) return NULL;
    mpfr_prec_t bits = (mpfr_prec_t)info.min_bits;
    if (bits < 53) bits = 53;

    bool square = (rows == cols);
    size_t total = (size_t)rows * (size_t)cols;

    mpfr_t* A_re = NULL;
    mpfr_t* A_im = NULL;
    bool is_complex = false;
    if (!lum_load_matrix(m, rows, cols, bits, &A_re, &A_im, &is_complex)) {
        tc_check_deadline();          /* a timeout while loading: nothing held */
        return NULL;
    }

    /* Snapshot ||A||_inf before the factorisation overwrites A.  Only
     * needed for the condition number, which is square-only. */
    mpfr_t anorm;
    mpfr_init2(anorm, bits);
    if (square) {
        lum_norm_inf(A_re, A_im, rows, cols, is_complex, bits, anorm);
    } else {
        mpfr_set_zero(anorm, 1);
    }

    /* perm: full row permutation, length `rows`.  Doolittle only
     * swaps the first `steps` entries; rows past steps stay at
     * identity.  Matches Mathematica's convention. */
    int* perm = (int*)malloc((size_t)rows * sizeof(int));
    for (int i = 0; i < rows; i++) perm[i] = i + 1;

    bool singular = false;
    bool ok = lum_factor(A_re, A_im, rows, cols, bits, is_complex,
                          perm, &singular);
    if (!ok) {
        mpfr_clear(anorm);
        lum_array_free(A_re, total);
        lum_array_free(A_im, is_complex ? total : 0);
        free(perm);
        tc_check_deadline();          /* a timeout: everything is freed */
        return NULL;
    }

    if (singular) {
        /* While an interrupt unwinds, the message is muted (message.c), so
         * the once-per-session flag must not be spent on it. */
        if (!sing_warn_counter && !mth_interrupt_pending()) {
            sing_warn_counter = 1;
            char* s = expr_to_string(m);
            mth_message("LUDecomposition", "sing", "Matrix %s is singular.", s);
            free(s);
        }
    }

    /* Condition number = ||A||_inf * ||A^{-1}||_inf.  Square only;
     * non-square gets exact Integer 0.
     *
     * Real matrices: Hager-Higham one-norm estimator on A^{-T}.  Cost
     * O(n^2) per iteration with 2-5 iterations typical -- matches
     * LAPACK's dgecon strategy.  Returns an estimate that is a lower
     * bound but typically within 10% of the true value.
     *
     * Complex matrices: keep the explicit-inverse approach (the HH
     * sign-vector iteration in the complex case requires a separate
     * implementation).  Cost O(n^3) but exact.
     *
     * Badly-conditioned warning (::luc): if cond > 2^(bits) the
     * matrix is too ill-conditioned for the working precision to
     * yield a meaningful factorisation.  Mirrors Mathematica's
     * LUDecomposition::luc behaviour; suppressed if ::sing already
     * fired (the kernel only computes condition for nonsing, so the
     * branch ordering naturally handles this). */
    static uint64_t luc_warn_counter = 0;
    Expr* c;
    if (!square) {
        c = expr_new_integer(0);
    } else {
        mpfr_t inv_norm, cond;
        mpfr_init2(inv_norm, bits);
        mpfr_init2(cond,     bits);
        bool nonsing;
        if (is_complex) {
            nonsing = lum_inv_norm_inf(A_re, A_im, perm, rows, is_complex,
                                        bits, inv_norm);
        } else {
            nonsing = lum_inv_norm_hager(A_re, A_im, perm, rows,
                                          bits, inv_norm);
        }
        if (nonsing) mpfr_mul(cond, anorm, inv_norm, MPFR_RNDN);
        else         mpfr_set_inf(cond, 1);

        if (nonsing && !luc_warn_counter && !mth_interrupt_pending()) {
            /* threshold = 2^(bits).  Compare cond > threshold. */
            mpfr_t threshold;
            mpfr_init2(threshold, bits);
            mpfr_set_ui_2exp(threshold, 1, (mpfr_exp_t)bits, MPFR_RNDN);
            if (mpfr_cmp(cond, threshold) > 0) {
                luc_warn_counter = 1;
                char* s = expr_to_string(m);
                mth_message("LUDecomposition", "luc", "Result for LUDecomposition of badly conditioned matrix %s may contain significant numerical errors.", s);
                free(s);
            }
            mpfr_clear(threshold);
        }
        c = expr_new_mpfr_copy(cond);
        mpfr_clear(inv_norm);
        mpfr_clear(cond);
    }

    Expr* lu = lum_build_lu(A_re, A_im, rows, cols, is_complex, bits);
    Expr* p  = lum_build_perm(perm, rows);

    lum_array_free(A_re, total);
    lum_array_free(A_im, is_complex ? total : 0);
    free(perm);
    mpfr_clear(anorm);

    Expr** items = (Expr**)malloc(sizeof(Expr*) * 3);
    items[0] = lu; items[1] = p; items[2] = c;
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), items, 3);
    free(items);
    return result;
}

/* ------------------------------------------------------------------ *
 *  Arbitrary-precision determinant.                                    *
 *                                                                       *
 *  Det[m] for a genuine MPFR (or otherwise inexact-numeric) n x n       *
 *  matrix.  builtin_det reaches this only after the exact FLINT path    *
 *  declines, and it exists to replace the O(n!) Laplace fallback that   *
 *  otherwise hangs at n >= ~12 for inexact input.                       *
 *                                                                       *
 *  det(A) = sign(P) * prod(U[k,k]), where PA = LU and sign(P) is the    *
 *  parity of the pivoting permutation.  Reuses lum_load_matrix /        *
 *  lum_factor.  Returns an EXPR_MPFR (real) or Complex[mpfr, mpfr]      *
 *  (complex), NULL on a non-numeric leaf / USE_MPFR off.  The wide MPFR *
 *  exponent means the product never spuriously overflows.               *
 * ------------------------------------------------------------------ */
static int lum_perm_parity(const int* perm, int n)
{
    /* perm holds a permutation of 1..n (1-indexed).  Parity = (-1)^(number
     * of transpositions) computed by cycle decomposition. */
    char* seen = (char*)calloc((size_t)n, 1);
    int transpositions = 0;
    for (int i = 0; i < n; i++) {
        if (seen[i]) continue;
        int j = i, len = 0;
        while (!seen[j]) { seen[j] = 1; j = perm[j] - 1; len++; }
        transpositions += len - 1;
    }
    free(seen);
    return (transpositions & 1) ? -1 : 1;
}

Expr* mpfr_det_dispatch(Expr* m, int n)
{
    if (n <= 0) return NULL;

    CommonInexactInfo info = common_scan_inexact(m);
    if (!info.has_inexact) return NULL;         /* exact input: not our job */
    mpfr_prec_t bits = (mpfr_prec_t)info.min_bits;
    if (bits < 53) bits = 53;

    size_t total = (size_t)n * (size_t)n;
    mpfr_t* A_re = NULL;
    mpfr_t* A_im = NULL;
    bool is_complex = false;
    if (!lum_load_matrix(m, n, n, bits, &A_re, &A_im, &is_complex)) {
        tc_check_deadline();                    /* a timeout while loading */
        return NULL;                            /* symbolic leaf: caller falls back */
    }

    int* perm = (int*)malloc((size_t)n * sizeof(int));
    for (int i = 0; i < n; i++) perm[i] = i + 1;

    bool singular = false;
    if (!lum_factor(A_re, A_im, n, n, bits, is_complex, perm, &singular)) {
        lum_array_free(A_re, total);
        lum_array_free(A_im, is_complex ? total : 0);
        free(perm);
        tc_check_deadline();          /* a timeout: everything is freed */
        return NULL;
    }

    int sign = lum_perm_parity(perm, n);
    free(perm);

    /* det = sign * prod(diag).  A zero pivot leaves a zero on the diagonal,
     * so a singular matrix yields det = 0 for free. */
    Expr* out;
    if (is_complex) {
        mpfr_t pr, pi, tr, ti, t;
        mpfr_inits2(bits, pr, pi, tr, ti, t, (mpfr_ptr)0);
        mpfr_set_si(pr, sign, MPFR_RNDN);
        mpfr_set_zero(pi, 1);
        for (int k = 0; k < n; k++) {
            size_t kk = (size_t)k * n + k;
            /* (pr + I pi)(A_re + I A_im) */
            mpfr_mul(tr, pr, A_re[kk], MPFR_RNDN); mpfr_mul(t, pi, A_im[kk], MPFR_RNDN);
            mpfr_sub(tr, tr, t, MPFR_RNDN);
            mpfr_mul(ti, pr, A_im[kk], MPFR_RNDN); mpfr_mul(t, pi, A_re[kk], MPFR_RNDN);
            mpfr_add(ti, ti, t, MPFR_RNDN);
            mpfr_set(pr, tr, MPFR_RNDN); mpfr_set(pi, ti, MPFR_RNDN);
        }
        /* Normalise signed zeros so a singular matrix reads 0, not -0. */
        if (mpfr_zero_p(pr)) mpfr_setsign(pr, pr, 0, MPFR_RNDN);
        if (mpfr_zero_p(pi)) mpfr_setsign(pi, pi, 0, MPFR_RNDN);
        out = lum_make_scalar(pr, pi, true);
        mpfr_clears(pr, pi, tr, ti, t, (mpfr_ptr)0);
    } else {
        mpfr_t prod;
        mpfr_init2(prod, bits);
        mpfr_set_si(prod, sign, MPFR_RNDN);
        for (int k = 0; k < n; k++)
            mpfr_mul(prod, prod, A_re[(size_t)k * n + k], MPFR_RNDN);
        if (mpfr_zero_p(prod)) mpfr_setsign(prod, prod, 0, MPFR_RNDN);
        out = expr_new_mpfr_move(prod);
    }

    lum_array_free(A_re, total);
    lum_array_free(A_im, is_complex ? total : 0);
    return out;
}

/* ------------------------------------------------------------------ *
 *  A reusable MPFR LU factor.                                          *
 *                                                                       *
 *  LinearSolve[m] for an arbitrary-precision m caches one of these (a  *
 *  LinearSolveFunction then costs two triangular solves per right-hand *
 *  side instead of a fresh elimination), and LUBackSubstitution builds *
 *  one from the {lu, p, c} triple LUDecomposition returns.             *
 * ------------------------------------------------------------------ */
struct MpfrLU {
    int n;
    mpfr_prec_t bits;
    bool cplx;
    mpfr_t* re;     /* n*n row-major combined L\U (unit L implicit)      */
    mpfr_t* im;     /* NULL for a real factor                            */
    int* perm;      /* 1-indexed: row i of LU is row perm[i] of the matrix */
};

void mpfr_lu_free(MpfrLU* f)
{
    if (!f) return;
    size_t t = (size_t)f->n * (size_t)f->n;
    lum_array_free(f->re, t);
    lum_array_free(f->im, t);
    free(f->perm);
    free(f);
}

long mpfr_lu_bits(const MpfrLU* f) { return f ? (long)f->bits : 0; }

MpfrLU* mpfr_lu_factor(const Expr* m, int n, long bits, bool* singular)
{
    *singular = false;
    if (n <= 0 || bits < 2) return NULL;
    MpfrLU* f = (MpfrLU*)calloc(1, sizeof *f);
    if (!f) return NULL;
    f->n = n;
    f->bits = (mpfr_prec_t)bits;
    if (!lum_load_matrix((Expr*)m, n, n, f->bits, &f->re, &f->im, &f->cplx)) {
        free(f);
        tc_check_deadline();          /* a timeout while loading: nothing held */
        return NULL;
    }
    f->perm = (int*)malloc((size_t)n * sizeof(int));
    for (int i = 0; i < n; i++) f->perm[i] = i + 1;

    mpfr_t maxabs, mag;
    mpfr_init2(maxabs, f->bits);
    mpfr_init2(mag, f->bits);
    mpfr_set_zero(maxabs, 1);
    for (size_t k = 0; k < (size_t)n * (size_t)n; k++) {
        lum_abs(f->re, f->im, k, f->cplx, mag);
        if (mpfr_cmp(mag, maxabs) > 0) mpfr_set(maxabs, mag, MPFR_RNDN);
    }
    bool zero_pivot = false;
    if (!lum_factor(f->re, f->im, n, n, f->bits, f->cplx, f->perm, &zero_pivot)) {
        mpfr_clear(maxabs);
        mpfr_clear(mag);
        mpfr_lu_free(f);
        tc_check_deadline();          /* a timeout: everything is freed */
        return NULL;
    }
    /* A pivot at or below 64 n 2^(1 - bits) max|m| is the roundoff left by
     * eliminating a singular matrix, not a value: the same zero tolerance
     * the numeric row reduction (rref_numeric.c) applies.  It is relative to
     * the working precision, never to machine epsilon, so a matrix whose
     * condition number is far beyond 2^53 but below 2^bits is regular. */
    mpfr_mul_ui(maxabs, maxabs, (unsigned long)(64 * n), MPFR_RNDN);
    mpfr_mul_2si(maxabs, maxabs, 1 - (long)f->bits, MPFR_RNDN);
    *singular = zero_pivot;
    for (int k = 0; k < n && !*singular; k++) {
        lum_abs(f->re, f->im, (size_t)k * n + k, f->cplx, mag);
        if (mpfr_cmp(mag, maxabs) <= 0) *singular = true;
    }
    mpfr_clear(maxabs);
    mpfr_clear(mag);
    return f;
}

MpfrLU* mpfr_lu_from_parts(const Expr* lu, const int* perm, int n, long bits)
{
    if (n <= 0 || bits < 2) return NULL;
    MpfrLU* f = (MpfrLU*)calloc(1, sizeof *f);
    if (!f) return NULL;
    f->n = n;
    f->bits = (mpfr_prec_t)bits;
    if (!lum_load_matrix((Expr*)lu, n, n, f->bits, &f->re, &f->im, &f->cplx)) {
        free(f);
        tc_check_deadline();          /* a timeout while loading: nothing held */
        return NULL;
    }
    f->perm = (int*)malloc((size_t)n * sizeof(int));
    memcpy(f->perm, perm, (size_t)n * sizeof(int));
    return f;
}

/* x -= a * y over (re, im) pairs; an absent imaginary part (NULL) is 0.
 * t is caller-owned scratch at the working precision. */
static void mlu_submul(mpfr_ptr xr, mpfr_ptr xi,
                       mpfr_srcptr ar, mpfr_srcptr ai,
                       mpfr_srcptr yr, mpfr_srcptr yi, mpfr_ptr t)
{
    mpfr_mul(t, ar, yr, MPFR_RNDN);
    mpfr_sub(xr, xr, t, MPFR_RNDN);
    if (!xi) return;
    if (ai && yi) {
        mpfr_mul(t, ai, yi, MPFR_RNDN);
        mpfr_add(xr, xr, t, MPFR_RNDN);
    }
    if (yi) {
        mpfr_mul(t, ar, yi, MPFR_RNDN);
        mpfr_sub(xi, xi, t, MPFR_RNDN);
    }
    if (ai) {
        mpfr_mul(t, ai, yr, MPFR_RNDN);
        mpfr_sub(xi, xi, t, MPFR_RNDN);
    }
}

/* The shape of a right-hand side: a vector of n scalars (*k = 0) or an
 * n x k matrix.  False for anything else. */
static bool mlu_rhs_shape(const Expr* b, int n, int* k)
{
    if (!b || b->type != EXPR_FUNCTION || b->data.function.head->type != EXPR_SYMBOL ||
        b->data.function.head->data.symbol.name != SYM_List ||
        (int)b->data.function.arg_count != n)
        return false;
    const Expr* r0 = b->data.function.args[0];
    bool is_mat = r0->type == EXPR_FUNCTION && r0->data.function.head->type == EXPR_SYMBOL &&
                  r0->data.function.head->data.symbol.name == SYM_List;
    if (!is_mat) { *k = 0; return true; }
    int kk = (int)r0->data.function.arg_count;
    if (kk == 0) return false;
    for (int i = 0; i < n; i++) {
        const Expr* ri = b->data.function.args[i];
        if (ri->type != EXPR_FUNCTION || ri->data.function.head->type != EXPR_SYMBOL ||
            ri->data.function.head->data.symbol.name != SYM_List ||
            (int)ri->data.function.arg_count != kk)
            return false;
    }
    *k = kk;
    return true;
}

Expr* mpfr_lu_solve(const MpfrLU* f, const Expr* b)
{
    if (!f) return NULL;
    int n = f->n, k = 0;
    if (!mlu_rhs_shape(b, n, &k)) return NULL;
    int w = k ? k : 1;
    size_t tot = (size_t)n * (size_t)w;
    mpfr_prec_t bits = f->bits;

    /* Load b[[perm]] (row i of the work array is row perm[i] of b). */
    bool cplx = f->cplx;
    for (int i = 0; i < n && !cplx; i++) {
        const Expr* ri = b->data.function.args[i];
        for (int c = 0; c < w && !cplx; c++)
            if (lum_leaf_is_complex(k ? ri->data.function.args[c] : (Expr*)ri)) cplx = true;
    }
    mpfr_t* xr = lum_array_alloc(tot, bits);
    mpfr_t* xi = cplx ? lum_array_alloc(tot, bits) : NULL;
    mpfr_t tim;
    mpfr_init2(tim, bits);
    bool ok = true;
    for (int i = 0; i < n && ok; i++) {
        const Expr* ri = b->data.function.args[f->perm[i] - 1];
        for (int c = 0; c < w && ok; c++) {
            const Expr* e = k ? ri->data.function.args[c] : ri;
            size_t o = (size_t)i * w + c;
            ok = get_approx_mpfr(e, xr[o], tim, NULL);
            if (ok && cplx) mpfr_set(xi[o], tim, MPFR_RNDN);
        }
    }
    if (!ok) {
        mpfr_clear(tim);
        lum_array_free(xr, tot);
        lum_array_free(xi, cplx ? tot : 0);
        return NULL;
    }

    mpfr_t t;
    mpfr_init2(t, bits);
    const mpfr_t* Lr = f->re;
    const mpfr_t* Li = f->im;
    /* Forward substitution through the unit-diagonal L. */
    for (int i = 0; i < n; i++)
        for (int j = 0; j < i; j++) {
            size_t ij = (size_t)i * n + j;
            if (mpfr_zero_p(Lr[ij]) && (!Li || mpfr_zero_p(Li[ij]))) continue;
            for (int c = 0; c < w; c++) {
                size_t xo = (size_t)i * w + c, yo = (size_t)j * w + c;
                mlu_submul(xr[xo], cplx ? xi[xo] : NULL, Lr[ij], Li ? Li[ij] : NULL,
                           xr[yo], cplx ? xi[yo] : NULL, t);
            }
        }
    /* Back substitution through U. */
    mpfr_t dr, di, den, ar, ai;
    mpfr_inits2(bits, dr, di, den, ar, ai, (mpfr_ptr)0);
    for (int i = n - 1; i >= 0; i--) {
        for (int j = i + 1; j < n; j++) {
            size_t ij = (size_t)i * n + j;
            if (mpfr_zero_p(Lr[ij]) && (!Li || mpfr_zero_p(Li[ij]))) continue;
            for (int c = 0; c < w; c++) {
                size_t xo = (size_t)i * w + c, yo = (size_t)j * w + c;
                mlu_submul(xr[xo], cplx ? xi[xo] : NULL, Lr[ij], Li ? Li[ij] : NULL,
                           xr[yo], cplx ? xi[yo] : NULL, t);
            }
        }
        size_t ii = (size_t)i * n + i;
        if (!cplx) {
            for (int c = 0; c < w; c++)
                mpfr_div(xr[(size_t)i * w + c], xr[(size_t)i * w + c], Lr[ii], MPFR_RNDN);
            continue;
        }
        /* 1/(dr + i di) = (dr - i di)/(dr^2 + di^2) */
        mpfr_set(dr, Lr[ii], MPFR_RNDN);
        if (Li) mpfr_set(di, Li[ii], MPFR_RNDN); else mpfr_set_zero(di, 1);
        mpfr_sqr(den, dr, MPFR_RNDN);
        mpfr_sqr(t, di, MPFR_RNDN);
        mpfr_add(den, den, t, MPFR_RNDN);
        for (int c = 0; c < w; c++) {
            size_t o = (size_t)i * w + c;
            mpfr_set(ar, xr[o], MPFR_RNDN);
            mpfr_set(ai, xi[o], MPFR_RNDN);
            mpfr_mul(xr[o], ar, dr, MPFR_RNDN);
            mpfr_mul(t, ai, di, MPFR_RNDN);
            mpfr_add(xr[o], xr[o], t, MPFR_RNDN);
            mpfr_div(xr[o], xr[o], den, MPFR_RNDN);
            mpfr_mul(xi[o], ai, dr, MPFR_RNDN);
            mpfr_mul(t, ar, di, MPFR_RNDN);
            mpfr_sub(xi[o], xi[o], t, MPFR_RNDN);
            mpfr_div(xi[o], xi[o], den, MPFR_RNDN);
        }
    }
    mpfr_clears(dr, di, den, ar, ai, (mpfr_ptr)0);
    mpfr_clear(t);

    /* Build the answer: a vector, or an n x k matrix. */
    mpfr_set_zero(tim, 1);
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    for (int i = 0; i < n; i++) {
        if (!k) {
            rows[i] = lum_make_scalar(xr[i], cplx ? xi[i] : tim, cplx);
            continue;
        }
        Expr** el = (Expr**)malloc(sizeof(Expr*) * (size_t)k);
        for (int c = 0; c < k; c++) {
            size_t o = (size_t)i * k + c;
            el[c] = lum_make_scalar(xr[o], cplx ? xi[o] : tim, cplx);
        }
        rows[i] = expr_new_function(expr_new_symbol(SYM_List), el, (size_t)k);
        free(el);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)n);
    free(rows);
    mpfr_clear(tim);
    lum_array_free(xr, tot);
    lum_array_free(xi, cplx ? tot : 0);
    return out;
}

#endif /* USE_MPFR */
