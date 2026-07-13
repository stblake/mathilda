/* svdecomp_machine.c
 *
 * Machine-precision LAPACK kernel for SingularValueDecomposition.
 *
 *   - Standard SVD:     dgesdd / zgesdd (divide-and-conquer).
 *   - Generalized SVD:  dggsvd3 / zggsvd3.
 *
 * Selected by svd_dispatch (svdecomp.c) when the input contains at
 * least one inexact leaf and the minimum precision is <= 53 bits.
 * On any soft failure (USE_LAPACK off, non-numeric leaf, LAPACK
 * info != 0) we return NULL without freeing the input; svd_dispatch
 * then falls through to the symbolic dispatcher.
 *
 * High-level flow for the standard SVD form (the only one wired up
 * here in Phase 3; the generalized form lands in Phase 6):
 *
 *   1. Load m into a column-major double buffer (interleaved (re, im)
 *      for complex).
 *   2. Call dgesdd('A', ...) -- "all m cols of U, all n rows of V^H"
 *      so we get the full square U (n x n), Sigma vector (length
 *      min(n, p)), and V^H (p x p) in one call.
 *   3. Wrap U row-major; build the rectangular n x p Sigma diagonal;
 *      take Conjugate-Transpose of V^H to get V.
 *   4. Apply Tolerance / Truncation / TargetStructure via
 *      svd_apply_postprocess.
 *
 * Memory contract.  Standard builtin contract (see svdecomp.h).
 */

#include "svdecomp.h"
#include "svdecomp_internal.h"
#include "linalg.h"
#include "lapack.h"
#include "sym_names.h"
#include "common.h"

#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef USE_LAPACK
/* The static helpers below are only reachable from the USE_LAPACK branch
 * of svd_machine_dispatch.  Wrapping them removes -Wunused-function noise
 * on builds compiled without LAPACK. */

/* ---------------------------------------------------------------------
 * One-shot warning helper for the LAPACK fast path.  Same pattern as
 * qr_machine_warn_once in qrdecomp_machine.c. */
static void svd_machine_warn_once(uint64_t* counter, const char* msg)
{
    if (*counter) return;
    *counter = 1;
    fprintf(stderr, "%s", msg);
}

/* ---------------------------------------------------------------------
 * Numeric-leaf -> double conversion.  Same recognised forms as the
 * other machine kernels: Integer, BigInt, Real, MPFR, Rational[p, q],
 * Complex[re, im].  Anything else triggers a fall-back to symbolic. */
static bool svdm_leaf_to_double(Expr* e, double* out_re, double* out_im)
{
    *out_im = 0.0;
    if (!e) return false;
    switch (e->type) {
        case EXPR_REAL:    *out_re = e->data.real;            return true;
        case EXPR_INTEGER: *out_re = (double)e->data.integer; return true;
        case EXPR_BIGINT:  *out_re = mpz_get_d(e->data.bigint); return true;
#ifdef USE_MPFR
        case EXPR_MPFR:
            *out_re = mpfr_get_d(e->data.mpfr, MPFR_RNDN);
            return true;
#endif
        case EXPR_FUNCTION:
            if (e->data.function.head->type == EXPR_SYMBOL) {
                const char* h = e->data.function.head->data.symbol.name;
                if (h == SYM_Rational && e->data.function.arg_count == 2) {
                    double p, q, dummy;
                    if (svdm_leaf_to_double(e->data.function.args[0], &p, &dummy)
                        && svdm_leaf_to_double(e->data.function.args[1], &q, &dummy)
                        && q != 0.0) {
                        *out_re = p / q;
                        return true;
                    }
                    return false;
                }
                if (h == SYM_Complex && e->data.function.arg_count == 2) {
                    double r, i, dummy;
                    if (svdm_leaf_to_double(e->data.function.args[0], &r, &dummy)
                        && svdm_leaf_to_double(e->data.function.args[1], &i, &dummy)) {
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

static bool svdm_leaf_is_complex(Expr* e)
{
    double r, i;
    return svdm_leaf_to_double(e, &r, &i) && i != 0.0;
}

/* Load an n x p Mathilda matrix into a freshly-allocated column-major
 * double buffer.  Returns true on success.  On failure (any leaf isn't
 * a numeric form) returns false with *out_A and *out_is_complex
 * untouched and nothing leaked.
 *
 * `force_complex` lets the generalized-SVD caller force the complex
 * code path when one of the two matrices has imaginary content -- both
 * must be loaded in the same layout so LAPACK sees matching types. */
static bool svdm_load_matrix(Expr* m, int n, int p, bool force_complex,
                             double** out_A, bool* out_is_complex)
{
    if (m->type != EXPR_FUNCTION) return false;

    bool is_complex = force_complex;
    if (!is_complex) {
        for (int i = 0; i < n && !is_complex; i++) {
            Expr* row = m->data.function.args[i];
            if (row->type != EXPR_FUNCTION) return false;
            for (int j = 0; j < p && !is_complex; j++) {
                if (svdm_leaf_is_complex(row->data.function.args[j])) {
                    is_complex = true;
                }
            }
        }
    }

    size_t stride = is_complex ? 2 : 1;
    double* A = (double*)malloc(stride * (size_t)n * (size_t)p * sizeof(double));
    if (!A) return false;

    for (int i = 0; i < n; i++) {
        Expr* row = m->data.function.args[i];
        if (row->type != EXPR_FUNCTION
            || (int)row->data.function.arg_count != p) {
            free(A);
            return false;
        }
        for (int j = 0; j < p; j++) {
            double re, im;
            if (!svdm_leaf_to_double(row->data.function.args[j], &re, &im)) {
                free(A);
                return false;
            }
            size_t off = stride * ((size_t)i + (size_t)j * (size_t)n);
            A[off] = re;
            if (is_complex) A[off + 1] = im;
        }
    }

    *out_A = A;
    *out_is_complex = is_complex;
    return true;
}

/* Build a Mathilda Complex[re, im] scalar, collapsing to a bare Real
 * when im == 0.0 (matching the existing numericalize convention). */
static Expr* svdm_scalar(double re, double im, bool is_complex)
{
    if (!is_complex || im == 0.0) return expr_new_real(re);
    Expr** args = (Expr**)malloc(sizeof(Expr*) * 2);
    args[0] = expr_new_real(re);
    args[1] = expr_new_real(im);
    Expr* z = expr_new_function(expr_new_symbol(SYM_Complex), args, 2);
    free(args);
    return z;
}

/* Wrap a column-major double buffer as a row-major Mathilda matrix.
 *
 *   U_cm[i + j * lda]   ->   matrix[i][j]
 *
 * For complex inputs each entry is two consecutive doubles; the
 * `transpose_conj` flag, when true, emits matrix[i][j] = conj(buf[j + i * lda])
 * instead -- used to convert LAPACK's V^H output into V. */
static Expr* svdm_wrap_cm(const double* buf, int rows, int cols,
                          int lda, bool is_complex, bool transpose_conj)
{
    size_t stride = is_complex ? 2 : 1;
    Expr** row_exprs = (Expr**)malloc(sizeof(Expr*) * (size_t)rows);
    for (int i = 0; i < rows; i++) {
        Expr** elems = (Expr**)malloc(sizeof(Expr*) * (size_t)cols);
        for (int j = 0; j < cols; j++) {
            size_t off;
            if (transpose_conj) {
                off = stride * ((size_t)j + (size_t)i * (size_t)lda);
            } else {
                off = stride * ((size_t)i + (size_t)j * (size_t)lda);
            }
            double re = buf[off];
            double im = is_complex ? buf[off + 1] : 0.0;
            if (is_complex && transpose_conj) im = -im;
            elems[j] = svdm_scalar(re, im, is_complex);
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

/* Build the rectangular n x p Sigma matrix from the length-mn singular
 * value vector (mn = min(n, p)).  Diagonal entries are sigma_i; off-
 * diagonal entries are the exact Integer 0 (matches Mathematica's
 * formatting, e.g. {{12.4778, 0.}, {0., 5.65202}, {0., 0.}}). */
static Expr* svdm_sigma_rect(const double* S, int mn, int n, int p)
{
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
    for (int i = 0; i < n; i++) {
        Expr** elems = (Expr**)malloc(sizeof(Expr*) * (size_t)p);
        for (int j = 0; j < p; j++) {
            if (i == j && i < mn) elems[j] = expr_new_real(S[i]);
            else                   elems[j] = expr_new_real(0.0);
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

/* ---------------------------------------------------------------------
 * Apply LAPACK's default rank-revealing cutoff to the singular value
 * vector.  Entries below the threshold are set to literal 0.0; this
 * matches dgesdd's own internal behaviour and lets the post-processing
 * Tolerance pass operate on a cleaner vector.  We do NOT zero anything
 * when the user supplied an explicit Tolerance option -- that case is
 * handled later in svd_apply_postprocess so the user's threshold takes
 * priority over the LAPACK default. */
static void svdm_default_tolerance(double* S, int mn, int n, int p)
{
    if (mn == 0) return;
    double max_n_p = (n > p) ? (double)n : (double)p;
    double tol = max_n_p * DBL_EPSILON * S[0];
    for (int i = 0; i < mn; i++) {
        if (S[i] < tol) S[i] = 0.0;
    }
}
#endif /* USE_LAPACK */

/* ------------------------------------------------------------------ *
 *  Generalized SVD (LAPACK dggsvd3 / zggsvd3).                        *
 *                                                                     *
 *  Builds {{u, ua}, {sigma, sigma_a}, v} from LAPACK's GSVD output.   *
 *  Mathematica's convention reconstructs                              *
 *    m == u . sigma . ConjugateTranspose[v]                           *
 *    a == ua . sigma_a . ConjugateTranspose[v]                        *
 *  and LAPACK returns                                                 *
 *    A = U . D1 . (0 R) . Q^H                                         *
 *    B = V . D2 . (0 R) . Q^H                                         *
 *  so we take u = U, ua = V, v = Q, and materialise sigma = D1.(0 R), *
 *  sigma_a = D2.(0 R).  Those matrices are rectangular (Lm x Ln and   *
 *  Lp x Ln) and pseudo-diagonal -- the alpha[i] / beta[i] generalised *
 *  singular values land at positions tied to the (K, L) block         *
 *  structure of D1 / D2.  We materialise them as dense Mathilda       *
 *  lists for transparency.                                            *
 *                                                                     *
 *  Variable names in this block:                                      *
 *    Lm = LAPACK m = rows of A = Mathilda n                           *
 *    Ln = LAPACK n = shared cols of A and B = Mathilda p              *
 *    Lp = LAPACK p = rows of B = Mathilda n_a                         *
 * ------------------------------------------------------------------ */
#ifdef USE_LAPACK
/* Scan both matrices of an {m, a} pair and detect imaginary content
 * anywhere.  Used so both matrices share the same complex/real layout
 * when handed to LAPACK. */
static bool svdm_gen_any_complex(Expr* m, int n, int p)
{
    if (m->type != EXPR_FUNCTION) return false;
    for (int i = 0; i < n; i++) {
        Expr* row = m->data.function.args[i];
        if (row->type != EXPR_FUNCTION) return false;
        for (int j = 0; j < p; j++) {
            if (svdm_leaf_is_complex(row->data.function.args[j])) return true;
        }
    }
    return false;
}

/* Extract R from LAPACK's destroyed A / B buffers into a packed
 * (kpl x kpl) column-major buffer.  Handles both layouts:
 *   - Lm >= kpl: R sits in A[0..kpl-1, Ln-kpl..Ln-1] in one block.
 *   - Lm <  kpl: R splits as R11 in A and R22 in B (LAPACK
 *                rearranges the upper-triangular factor across the
 *                two matrices when A doesn't have enough rows).
 * stride = 1 for real, 2 for complex (each (r, c) holds (re, im) pair).
 * Returns NULL on allocation failure. */
static double* svdm_extract_R(const double* A, int lda,
                              const double* B, int ldb,
                              int Lm, int Ln, int K, int L,
                              bool is_complex)
{
    int kpl = K + L;
    size_t stride = is_complex ? 2 : 1;
    if (kpl == 0) {
        double* R = (double*)calloc(stride, sizeof(double));
        return R;  /* degenerate but valid */
    }
    double* R = (double*)calloc(stride * (size_t)kpl * (size_t)kpl,
                                  sizeof(double));
    if (!R) return NULL;
    int zero_cols = Ln - kpl;
    if (Lm >= kpl) {
        for (int c = 0; c < kpl; c++) {
            for (int r = 0; r <= c; r++) {  /* upper triangular */
                size_t dst = stride * ((size_t)r + (size_t)c * (size_t)kpl);
                size_t src = stride * ((size_t)r
                                       + (size_t)(zero_cols + c) * (size_t)lda);
                R[dst] = A[src];
                if (is_complex) R[dst + 1] = A[src + 1];
            }
        }
    } else {
        /* R11 is at A[0..Lm-1, zero_cols..Ln-1]. */
        for (int c = 0; c < kpl; c++) {
            int r_lim = (c < Lm - 1) ? c : Lm - 1;
            for (int r = 0; r <= r_lim; r++) {
                size_t dst = stride * ((size_t)r + (size_t)c * (size_t)kpl);
                size_t src = stride * ((size_t)r
                                       + (size_t)(zero_cols + c) * (size_t)lda);
                R[dst] = A[src];
                if (is_complex) R[dst + 1] = A[src + 1];
            }
        }
        /* R22 is at B[Lm-K..L-1, Ln+Lm-K-L..Ln-1] -> R[Lm..kpl-1, Lm..kpl-1].
         * Derivation: B row index = (Lm-K) + (r-Lm) = r-K; col index =
         * (zero_cols+Lm) + (c-Lm) = zero_cols + c. */
        for (int c = Lm; c < kpl; c++) {
            for (int r = Lm; r <= c; r++) {
                size_t dst = stride * ((size_t)r + (size_t)c * (size_t)kpl);
                size_t src = stride * ((size_t)(r - K)
                                       + (size_t)(zero_cols + c) * (size_t)ldb);
                R[dst] = B[src];
                if (is_complex) R[dst + 1] = B[src + 1];
            }
        }
    }
    return R;
}

/* Build sigma_m = D1 . (0 R) as an Lm x Ln row-major Mathilda matrix.
 * D1 has structure:
 *   row i in [0, min(K, Lm)):       identity (picks up R[i, jp])
 *   row i in [K, min(K+L, Lm)):     scales row i of R by alpha[i]
 *   row i in [K+L, Lm):             all zero
 * Columns 0..(Ln-K-L-1) of sigma_m are always zero (the leading-zero
 * block of (0 R)).  alpha values are always real (LAPACK convention);
 * for complex R, the scaling is real x complex. */
static Expr* svdm_build_sigma_m(const double* R, int kpl,
                                 const double* alpha,
                                 int Lm, int Ln, int K, int L,
                                 bool is_complex)
{
    int zero_cols = Ln - kpl;
    int row_K  = (K     < Lm) ? K     : Lm;
    int row_KL = (K + L < Lm) ? K + L : Lm;
    size_t stride = is_complex ? 2 : 1;

    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)Lm);
    for (int i = 0; i < Lm; i++) {
        Expr** elems = (Expr**)malloc(sizeof(Expr*) * (size_t)Ln);
        for (int j = 0; j < Ln; j++) {
            double val_re = 0.0, val_im = 0.0;
            if (j >= zero_cols) {
                int jp = j - zero_cols;
                size_t r_off = stride * ((size_t)i + (size_t)jp * (size_t)kpl);
                if (i < row_K) {
                    val_re = R[r_off];
                    if (is_complex) val_im = R[r_off + 1];
                } else if (i < row_KL) {
                    val_re = alpha[i] * R[r_off];
                    if (is_complex) val_im = alpha[i] * R[r_off + 1];
                }
            }
            elems[j] = svdm_scalar(val_re, val_im, is_complex);
        }
        rows[i] = expr_new_function(expr_new_symbol(SYM_List),
                                      elems, (size_t)Ln);
        free(elems);
    }
    Expr* m = expr_new_function(expr_new_symbol(SYM_List),
                                  rows, (size_t)Lm);
    free(rows);
    return m;
}

/* Build sigma_a = D2 . (0 R) as an Lp x Ln row-major Mathilda matrix.
 * D2 has only one nonzero per row: D2[i, K+i] = beta[K+i] for i in
 * [0, L) and zero elsewhere (works uniformly across both M >= K+L and
 * M < K+L cases because LAPACK fills beta[M..K+L-1] = 1 in the
 * second case, which is exactly the identity-block value).  beta is
 * always real; R may be complex. */
static Expr* svdm_build_sigma_a(const double* R, int kpl,
                                 const double* beta,
                                 int Lp, int Ln, int K, int L,
                                 bool is_complex)
{
    int zero_cols = Ln - kpl;
    int row_L = (L < Lp) ? L : Lp;
    size_t stride = is_complex ? 2 : 1;

    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)Lp);
    for (int i = 0; i < Lp; i++) {
        Expr** elems = (Expr**)malloc(sizeof(Expr*) * (size_t)Ln);
        for (int j = 0; j < Ln; j++) {
            double val_re = 0.0, val_im = 0.0;
            if (j >= zero_cols && i < row_L) {
                int jp = j - zero_cols;
                size_t r_off = stride * ((size_t)(K + i)
                                          + (size_t)jp * (size_t)kpl);
                val_re = beta[K + i] * R[r_off];
                if (is_complex) val_im = beta[K + i] * R[r_off + 1];
            }
            elems[j] = svdm_scalar(val_re, val_im, is_complex);
        }
        rows[i] = expr_new_function(expr_new_symbol(SYM_List),
                                      elems, (size_t)Ln);
        free(elems);
    }
    Expr* m = expr_new_function(expr_new_symbol(SYM_List),
                                  rows, (size_t)Lp);
    free(rows);
    return m;
}

/* Generalized-SVD entry.  Loads {m, a}, calls dggsvd3 or zggsvd3 based
 * on whether either matrix has imaginary content, assembles the
 * {{u, ua}, {sigma, sigma_a}, v} result, and hands off to the shared
 * post-processing.  Returns NULL on any soft failure (LAPACK
 * unavailable, non-numeric leaf, LAPACK info != 0). */
static Expr* svdm_dispatch_generalized(const SvdArgs* args,
                                         int n, int p, int n_a)
{
    static uint64_t lapack_warn_counter = 0;

    int Lm = n, Ln = p, Lp = n_a;

    /* Force a uniform complex layout for both A and B if either has any
     * imaginary content -- LAPACK needs them in matching types. */
    bool is_complex = svdm_gen_any_complex(args->m, Lm, Ln)
                    || svdm_gen_any_complex(args->a, Lp, Ln);

    double* A_cm = NULL;
    double* B_cm = NULL;
    bool dummy = false;
    if (!svdm_load_matrix(args->m, Lm, Ln, is_complex, &A_cm, &dummy))
        return NULL;
    if (!svdm_load_matrix(args->a, Lp, Ln, is_complex, &B_cm, &dummy)) {
        free(A_cm);
        return NULL;
    }

    size_t stride = is_complex ? 2 : 1;
    int lda = Lm, ldb = Lp, ldu = Lm, ldv = Lp, ldq = Ln;
    int K = 0, L = 0;
    double* alpha = (double*)malloc((size_t)Ln * sizeof(double));
    double* beta  = (double*)malloc((size_t)Ln * sizeof(double));
    double* U     = (double*)malloc(stride * (size_t)Lm * (size_t)Lm
                                     * sizeof(double));
    double* V     = (double*)malloc(stride * (size_t)Lp * (size_t)Lp
                                     * sizeof(double));
    double* Q     = (double*)malloc(stride * (size_t)Ln * (size_t)Ln
                                     * sizeof(double));
    if (!alpha || !beta || !U || !V || !Q) {
        free(A_cm); free(B_cm);
        free(alpha); free(beta); free(U); free(V); free(Q);
        return NULL;
    }

    int info = is_complex
        ? mat_lapack_zggsvd3('U', 'V', 'Q', Lm, Ln, Lp, &K, &L,
                              A_cm, lda, B_cm, ldb,
                              alpha, beta,
                              U, ldu, V, ldv, Q, ldq)
        : mat_lapack_dggsvd3('U', 'V', 'Q', Lm, Ln, Lp, &K, &L,
                              A_cm, lda, B_cm, ldb,
                              alpha, beta,
                              U, ldu, V, ldv, Q, ldq);
    if (info != 0) {
        if (info < 0) {
            svd_machine_warn_once(&lapack_warn_counter,
                "SingularValueDecomposition: LAPACK fast path unavailable "
                "for generalized SVD; leaving call unevaluated.\n");
        } else {
            svd_machine_warn_once(&lapack_warn_counter,
                "SingularValueDecomposition: LAPACK generalized SVD "
                "reported a numerical issue; leaving call unevaluated.\n");
        }
        free(A_cm); free(B_cm);
        free(alpha); free(beta); free(U); free(V); free(Q);
        return NULL;
    }

    double* R = svdm_extract_R(A_cm, lda, B_cm, ldb, Lm, Ln, K, L,
                                 is_complex);
    free(A_cm); free(B_cm);
    if (!R) {
        free(alpha); free(beta); free(U); free(V); free(Q);
        return NULL;
    }
    int kpl = K + L;

    Expr* u_mat   = svdm_wrap_cm(U, Lm, Lm, ldu, is_complex, false);
    Expr* ua_mat  = svdm_wrap_cm(V, Lp, Lp, ldv, is_complex, false);
    Expr* v_mat   = svdm_wrap_cm(Q, Ln, Ln, ldq, is_complex, false);
    Expr* sigma_m = svdm_build_sigma_m(R, kpl, alpha, Lm, Ln, K, L,
                                         is_complex);
    Expr* sigma_a = svdm_build_sigma_a(R, kpl, beta,  Lp, Ln, K, L,
                                         is_complex);

    free(U); free(V); free(Q); free(R); free(alpha); free(beta);

    Expr** u_pair_args = (Expr**)malloc(sizeof(Expr*) * 2);
    u_pair_args[0] = u_mat;  u_pair_args[1] = ua_mat;
    Expr* u_pair = expr_new_function(expr_new_symbol(SYM_List),
                                       u_pair_args, 2);
    free(u_pair_args);

    Expr** s_pair_args = (Expr**)malloc(sizeof(Expr*) * 2);
    s_pair_args[0] = sigma_m; s_pair_args[1] = sigma_a;
    Expr* s_pair = expr_new_function(expr_new_symbol(SYM_List),
                                       s_pair_args, 2);
    free(s_pair_args);

    Expr** items = (Expr**)malloc(sizeof(Expr*) * 3);
    items[0] = u_pair; items[1] = s_pair; items[2] = v_mat;
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), items, 3);
    free(items);

    return svd_apply_postprocess(result, args, n, p, kpl);
}
#endif /* USE_LAPACK */

/* ------------------------------------------------------------------ *
 *  Machine-precision SVD kernel entry.                                *
 * ------------------------------------------------------------------ */
Expr* svd_machine_dispatch(const SvdArgs* args, int n, int p, int n_a)
{
#ifndef USE_LAPACK
    (void)args; (void)n; (void)p; (void)n_a;
    return NULL;
#else
    static uint64_t lapack_warn_counter = 0;
    if (args->generalized) {
        return svdm_dispatch_generalized(args, n, p, n_a);
    }
    (void)n_a;

    /* Load m into a column-major double buffer. */
    double* A_cm = NULL;
    bool is_complex = false;
    if (!svdm_load_matrix(args->m, n, p, false, &A_cm, &is_complex)) {
        return NULL;
    }

    int mn = (n < p) ? n : p;
    int lda  = n;
    int ldu  = n;
    int ldvt = p;
    size_t stride = is_complex ? 2 : 1;

    double* S  = (double*)malloc((size_t)mn * sizeof(double));
    double* U  = (double*)malloc(stride * (size_t)n * (size_t)n * sizeof(double));
    double* VT = (double*)malloc(stride * (size_t)p * (size_t)p * sizeof(double));
    if (!S || !U || !VT) {
        if (S) free(S);
        if (U) free(U);
        if (VT) free(VT);
        free(A_cm);
        return NULL;
    }

    int info = is_complex
        ? mat_lapack_zgesdd('A', n, p, A_cm, lda, S, U, ldu, VT, ldvt)
        : mat_lapack_dgesdd('A', n, p, A_cm, lda, S, U, ldu, VT, ldvt);
    free(A_cm);

    if (info != 0) {
        if (info < 0) {
            svd_machine_warn_once(&lapack_warn_counter,
                "SingularValueDecomposition: LAPACK fast path "
                "unavailable; falling back to symbolic kernel.\n");
        } else {
            svd_machine_warn_once(&lapack_warn_counter,
                "SingularValueDecomposition: LAPACK divide-and-conquer "
                "SVD reported a numerical issue; falling back to "
                "symbolic kernel.\n");
        }
        free(S); free(U); free(VT);
        return NULL;
    }

    /* Apply the default rank-revealing cutoff only when the user did
     * NOT supply an explicit Tolerance (otherwise we'd double-cut and
     * the user's threshold would be irrelevant). */
    if (!args->tolerance) {
        svdm_default_tolerance(S, mn, n, p);
    }

    /* Build u (n x n), sigma (n x p), v (p x p, = ConjugateTranspose[VT]). */
    Expr* u_mat     = svdm_wrap_cm(U, n, n, ldu, is_complex, false);
    Expr* sigma_mat = svdm_sigma_rect(S, mn, n, p);
    Expr* v_mat     = svdm_wrap_cm(VT, p, p, ldvt, is_complex, true);
    free(S); free(U); free(VT);

    /* Build {u, sigma, v}. */
    Expr** items = (Expr**)malloc(sizeof(Expr*) * 3);
    items[0] = u_mat; items[1] = sigma_mat; items[2] = v_mat;
    Expr* result = expr_new_function(expr_new_symbol(SYM_List), items, 3);
    free(items);

    /* Hand the result to the shared post-processing for truncation,
     * tolerance, and TargetStructure. */
    return svd_apply_postprocess(result, args, n, p, mn);
#endif /* USE_LAPACK */
}
