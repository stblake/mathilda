#include "eigen.h"
#include "eigen_internal.h"
#include "linalg.h"
#include "eval.h"
#include "symtab.h"
#include "attr.h"
#include "poly.h"
#include "sym_names.h"
#include "sym_intern.h"
#include "common.h"
#include "numeric.h"
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <gmp.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif


/* ============================================================ *
 *  Phase 4: numerical "Banded" method (Hermitian-only).          *
 *                                                                 *
 *  For a real symmetric (or complex Hermitian) matrix with        *
 *  half-bandwidth b << n, we reduce to symmetric tridiagonal      *
 *  via two-sided Givens rotations that zero one off-band entry    *
 *  at a time and chase the resulting bulge down the matrix        *
 *  (Schwarz 1968 / Rutishauser).  The resulting symmetric         *
 *  tridiagonal eigenproblem is then handed off to the existing    *
 *  Phase 2 symmetric tridiag QR (direct_symtridiag_qr).           *
 *                                                                 *
 *  Dense storage (n*n) is used so we share workspace and printers *
 *  with the Direct kernels.  Each Givens applies as O(n) row/col  *
 *  updates; with b initial rotations per outer step and ~n/b      *
 *  chases each, the total cost is O(n^3) operations with the      *
 *  dense form.  A future sparse-aware Givens variant (touching    *
 *  only the band plus the in-flight bulge) would bring this down  *
 *  to the asymptotic O(n^2 b) the algorithm promises.  Even with  *
 *  dense storage, the cache footprint and absolute flop count are *
 *  smaller than dense Householder when b << n (no rank-2 trailing *
 *  block updates), so for narrowly-banded Hermitian matrices the  *
 *  Banded kernel still wins in practice.                          *
 *                                                                 *
 *  Complex Hermitian banded uses paired re/im Givens and reuses   *
 *  the Phase 2 phase-correction step (direct_phase_correct_tridiag)
 *  so the complex tridiagonal collapses to a real symmetric one   *
 *  before the final QR sweep.                                     *
 *                                                                 *
 *  LAPACK-HOOK: this whole block maps to dsbtrd / zhbtrd (band -> *
 *  tridiagonal) followed by dstedc / dsteqr (tridiagonal QR), or  *
 *  the combined wrappers dsbevd / zhbevd, when USE_LAPACK is set. *
 * ============================================================ */

/* Banded sub-options (presently a stub -- the Method -> "Banded"
 * form takes no sub-options in this initial implementation, but the
 * machinery is here so future "Bandwidth" overrides or QR-tolerance
 * knobs can be added without touching the call sites). */
typedef struct {
    int placeholder;
} BandedOpts;

static void banded_set_defaults(BandedOpts* o) { (void)o; }

static void banded_parse_subopts(Expr* method_value, BandedOpts* opts) {
    (void)method_value;
    banded_set_defaults(opts);
}

/* Heuristic: prefer Banded when the matrix is Hermitian and the
 * half-bandwidth is small enough that the reduction's b/n constant
 * advantage pays off.  Threshold: b <= max(8, n/4).  Anything denser
 * just goes through Direct, which has fewer Givens to chase. */
static bool banded_automatic_prefers_real(const MatD* A, size_t b,
                                            double sym_tol) {
    (void)sym_tol;
    size_t n = A->n;
    if (n <= 8) return false;       /* tiny -- Direct is faster + simpler */
    if (b >= n - 1) return false;   /* fully dense band */
    size_t thresh = (n / 4 > 8) ? n / 4 : 8;
    return b <= thresh;
}

/* Apply two-sided real Givens rotation G^T A G symmetrically.
 * G = [c s; -s c] in plane (p, q).  Full-width row / col updates so
 * the routine is correct regardless of where the in-flight bulge sits
 * during a chase.  When `want_Q`, post-multiplies Q by G in cols p, q.
 *
 * The two passes (row then col) are mathematically equivalent to
 * computing (G^T A) G in one shot but keep the indexing simple at the
 * cost of O(n) extra arithmetic per rotation -- acceptable for the
 * Banded kernel since the dominant cost is still O(n^2 b) total. */
static void banded_givens_two_sided_real(double* A, size_t n,
                                          size_t p, size_t q,
                                          double c, double s,
                                          double* Q, bool want_Q) {
    for (size_t j = 0; j < n; j++) {
        double ap = A[p * n + j];
        double aq = A[q * n + j];
        A[p * n + j] =  c * ap + s * aq;
        A[q * n + j] = -s * ap + c * aq;
    }
    for (size_t i = 0; i < n; i++) {
        double ip = A[i * n + p];
        double iq = A[i * n + q];
        A[i * n + p] =  c * ip + s * iq;
        A[i * n + q] = -s * ip + c * iq;
    }
    if (want_Q) {
        for (size_t i = 0; i < n; i++) {
            double qp = Q[i * n + p];
            double qq = Q[i * n + q];
            Q[i * n + p] =  c * qp + s * qq;
            Q[i * n + q] = -s * qp + c * qq;
        }
    }
}

/* Reduce a real symmetric banded matrix A (n*n, half-bandwidth b) to
 * symmetric tridiagonal form in place via Givens rotations with
 * bulge chasing.  On output:
 *   diag[i]   = A[i, i]
 *   sub[i]    = A[i + 1, i]   (i in [0, n-2])
 *   A         = symmetric tridiagonal (entries beyond the tridiagonal
 *               are numerically zero).
 *   Q         = orthogonal accumulator (n*n), updated when want_Q.
 *
 * Algorithm (one annihilation per inner step):
 *   For each "anchor column" k = 0..n-3:
 *     For d = b, b-1, ..., 2:
 *       Apply Givens in plane (k+d-1, k+d) to zero A[k+d, k] using the
 *       in-band element A[k+d-1, k].  By symmetry, A[k, k+d] is also
 *       zeroed.  This introduces a bulge at A[k+d-1, k+d+b] (and at
 *       the symmetric (k+d+b, k+d-1)).
 *
 *       Chase: as long as the bulge column is within range, apply
 *       another Givens in plane (cq + b - 1, cq + b) on the bulge
 *       column to eliminate the bulge entry -- which then propagates
 *       a new bulge b columns further along.  The chain terminates
 *       when the new bulge would fall past column n - 1.
 *
 *   After all outer steps the matrix is tridiagonal; we copy the
 *   diag / sub-diag out for direct_symtridiag_qr.
 */
static void banded_to_tridiag_real_sym(double* A, size_t n, size_t b,
                                         double* diag, double* sub,
                                         double* Q, bool want_Q) {
    if (want_Q) {
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++) Q[i * n + j] = (i == j) ? 1.0 : 0.0;
    }

    if (n >= 3 && b >= 2) {
        for (size_t k = 0; k + 2 < n; k++) {
            size_t d_max = b;
            if (k + d_max >= n) d_max = n - 1 - k;
            for (size_t d = d_max; d >= 2; d--) {
                size_t p = k + d - 1;
                size_t q = k + d;

                double a = A[p * n + k];
                double bb = A[q * n + k];
                if (bb == 0.0) continue;
                double r = hypot(a, bb);
                if (r == 0.0) continue;
                double c = a / r;
                double sn = bb / r;

                banded_givens_two_sided_real(A, n, p, q, c, sn, Q, want_Q);
                A[q * n + k] = 0.0;
                A[k * n + q] = 0.0;

                /* Chase the bulge introduced at (p, q + b) / (q + b, p). */
                size_t cp = p;
                size_t cq = q;
                while (cq + b < n) {
                    size_t bcol = cp;          /* bulge sits at row cq+b, col cp */
                    size_t new_p = cq + b - 1;
                    size_t new_q = cq + b;

                    double aa  = A[new_p * n + bcol];
                    double bb2 = A[new_q * n + bcol];
                    if (bb2 == 0.0) break;
                    double rr  = hypot(aa, bb2);
                    if (rr == 0.0) break;
                    double cc  = aa  / rr;
                    double ss  = bb2 / rr;

                    banded_givens_two_sided_real(A, n, new_p, new_q,
                                                   cc, ss, Q, want_Q);
                    A[new_q * n + bcol] = 0.0;
                    A[bcol * n + new_q] = 0.0;

                    cp = new_p;
                    cq = new_q;
                }
            }
        }
    }

    for (size_t i = 0; i < n; i++) diag[i] = A[i * n + i];
    for (size_t i = 0; i + 1 < n; i++) sub[i] = A[(i + 1) * n + i];
}

/* Top-level "Banded" kernel for a real symmetric matrix at machine
 * precision.  Mirrors direct_real_sym_machine; the only difference is
 * the tridiagonalisation step uses Givens band reduction. */
static Expr* banded_real_sym_machine(const MatD* A, size_t b,
                                       MateigenWant want, Expr* k_spec) {
    size_t n = A->n;
    if (n == 0) return NULL;

    double* W = (double*)malloc(sizeof(double) * n * n);
    memcpy(W, A->re, sizeof(double) * n * n);

    double* diag = (double*)malloc(sizeof(double) * n);
    double* sub  = (double*)calloc(n, sizeof(double));   /* n-1 + slack */
    bool want_Q  = (want & MATEIGEN_WANT_VECTORS) != 0;
    double* Q    = want_Q ? (double*)malloc(sizeof(double) * n * n) : NULL;

    banded_to_tridiag_real_sym(W, n, b, diag, sub, Q, want_Q);

    int qr_status = direct_symtridiag_qr(diag, sub, n, Q, want_Q);
    free(W);

    if (qr_status != 0) {
        free(diag); free(sub);
        if (Q) free(Q);
        return NULL;
    }

    size_t* perm = (size_t*)malloc(sizeof(size_t) * n);
    direct_sort_perm_desc_abs(diag, n, perm);

    Expr* out = (want_Q)
        ? direct_build_real_eigenvector_list(Q, n, perm)
        : direct_build_real_eigenvalue_list(diag, n, perm);

    free(diag); free(sub); free(perm);
    if (Q) free(Q);

    return direct_apply_k_spec_list(out, k_spec);
}

/* Apply two-sided complex Hermitian Givens U^H A U symmetrically.
 *   U = [c  -s_conj; s  c]  (c real, |c|^2 + |s|^2 = 1).
 *
 * Row update (left multiply by U^H):
 *   row p' =  c row_p + s_conj row_q
 *   row q' = -s row_p + c row_q
 *
 * Col update (right multiply by U):
 *   col p' =  c col_p + s col_q
 *   col q' = -s_conj col_p + c col_q
 *
 * Each "row" and "col" entry is a complex number stored as paired
 * (re, im) doubles; the formulas expand to standard 2x2 complex
 * arithmetic. */
static void banded_givens_two_sided_complex(double* A_re, double* A_im,
                                              size_t n,
                                              size_t p, size_t q,
                                              double c,
                                              double s_re, double s_im,
                                              double* Q_re, double* Q_im,
                                              bool want_Q) {
    /* Row update. */
    for (size_t j = 0; j < n; j++) {
        double apr = A_re[p * n + j];
        double api = A_im[p * n + j];
        double aqr = A_re[q * n + j];
        double aqi = A_im[q * n + j];
        /* p' = c*p + s_conj*q;   s_conj = s_re - i s_im
         *      (s_re - i s_im)(aqr + i aqi)
         *    = (s_re*aqr + s_im*aqi) + i (s_re*aqi - s_im*aqr) */
        A_re[p * n + j] = c * apr + s_re * aqr + s_im * aqi;
        A_im[p * n + j] = c * api + s_re * aqi - s_im * aqr;
        /* q' = -s*p + c*q;   -s = -(s_re + i s_im)
         *       -(s_re + i s_im)(apr + i api)
         *     = -(s_re*apr - s_im*api) - i (s_re*api + s_im*apr) */
        A_re[q * n + j] = -(s_re * apr - s_im * api) + c * aqr;
        A_im[q * n + j] = -(s_re * api + s_im * apr) + c * aqi;
    }
    /* Col update. */
    for (size_t i = 0; i < n; i++) {
        double ipr = A_re[i * n + p];
        double ipi = A_im[i * n + p];
        double iqr = A_re[i * n + q];
        double iqi = A_im[i * n + q];
        /* col p' = c*col_p + s*col_q;   s = s_re + i s_im
         *         (s_re + i s_im)(iqr + i iqi)
         *       = (s_re*iqr - s_im*iqi) + i (s_re*iqi + s_im*iqr) */
        A_re[i * n + p] = c * ipr + s_re * iqr - s_im * iqi;
        A_im[i * n + p] = c * ipi + s_re * iqi + s_im * iqr;
        /* col q' = -s_conj*col_p + c*col_q;   -s_conj = -(s_re - i s_im)
         *           -(s_re - i s_im)(ipr + i ipi)
         *         = -(s_re*ipr + s_im*ipi) + i (s_im*ipr - s_re*ipi) */
        A_re[i * n + q] = -(s_re * ipr + s_im * ipi) + c * iqr;
        A_im[i * n + q] =  (s_im * ipr - s_re * ipi) + c * iqi;
    }
    if (want_Q) {
        for (size_t i = 0; i < n; i++) {
            double qpr = Q_re[i * n + p];
            double qpi = Q_im[i * n + p];
            double qqr = Q_re[i * n + q];
            double qqi = Q_im[i * n + q];
            /* Same right-multiply pattern as the col update above. */
            Q_re[i * n + p] = c * qpr + s_re * qqr - s_im * qqi;
            Q_im[i * n + p] = c * qpi + s_re * qqi + s_im * qqr;
            Q_re[i * n + q] = -(s_re * qpr + s_im * qpi) + c * qqr;
            Q_im[i * n + q] =  (s_im * qpr - s_re * qpi) + c * qqi;
        }
    }
}

/* Reduce a complex Hermitian banded matrix A (n*n, half-bandwidth b)
 * to a complex Hermitian tridiagonal form with real diagonal and
 * generally complex sub-diagonal.  Caller must then run
 * direct_phase_correct_tridiag to make the sub-diagonal real-positive
 * before handing off to direct_symtridiag_qr.
 *
 * Givens parameters (Wilkinson's "real-c complex-s" form):
 *   r       = hypot(|a|, |b|)
 *   c       = |a| / r                                (real)
 *   s       = b * conj(a) / (|a| * r)                (complex)
 *   if |a| == 0:    c = 0, s = b / |b|               (pure swap)
 * After U^H [a; b] = [phase(a)*r; 0], so the chosen entry becomes a
 * complex of magnitude r and phase carried from the upper component.
 * That residual phase is mopped up later by the phase-correction step. */
static void banded_to_tridiag_complex_hermitian(double* A_re, double* A_im,
                                                  size_t n, size_t b,
                                                  double* diag,
                                                  double* sub_re,
                                                  double* sub_im,
                                                  double* Q_re, double* Q_im,
                                                  bool want_Q) {
    if (want_Q) {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                Q_re[i * n + j] = (i == j) ? 1.0 : 0.0;
                Q_im[i * n + j] = 0.0;
            }
        }
    }

    if (n >= 3 && b >= 2) {
        for (size_t k = 0; k + 2 < n; k++) {
            size_t d_max = b;
            if (k + d_max >= n) d_max = n - 1 - k;
            for (size_t d = d_max; d >= 2; d--) {
                size_t p = k + d - 1;
                size_t q = k + d;

                double ar = A_re[p * n + k];
                double ai = A_im[p * n + k];
                double br = A_re[q * n + k];
                double bi = A_im[q * n + k];
                double mag_a = hypot(ar, ai);
                double mag_b = hypot(br, bi);
                if (mag_b == 0.0) continue;
                double r = hypot(mag_a, mag_b);
                if (r == 0.0) continue;

                double c, s_re, s_im;
                if (mag_a == 0.0) {
                    /* Pure swap. */
                    c = 0.0;
                    s_re = br / mag_b;
                    s_im = bi / mag_b;
                } else {
                    c = mag_a / r;
                    /* s = b * conj(a) / (|a| * r)
                     *   = ((br + i bi)(ar - i ai)) / (|a| * r)
                     *   = (br*ar + bi*ai) + i (bi*ar - br*ai)   / (|a| * r) */
                    double den = mag_a * r;
                    s_re = (br * ar + bi * ai) / den;
                    s_im = (bi * ar - br * ai) / den;
                }

                banded_givens_two_sided_complex(A_re, A_im, n, p, q,
                                                  c, s_re, s_im,
                                                  Q_re, Q_im, want_Q);
                A_re[q * n + k] = 0.0; A_im[q * n + k] = 0.0;
                A_re[k * n + q] = 0.0; A_im[k * n + q] = 0.0;

                size_t cp = p;
                size_t cq = q;
                while (cq + b < n) {
                    size_t bcol  = cp;
                    size_t new_p = cq + b - 1;
                    size_t new_q = cq + b;

                    double aar = A_re[new_p * n + bcol];
                    double aai = A_im[new_p * n + bcol];
                    double bb_r = A_re[new_q * n + bcol];
                    double bb_i = A_im[new_q * n + bcol];
                    double ma  = hypot(aar, aai);
                    double mb  = hypot(bb_r, bb_i);
                    if (mb == 0.0) break;
                    double rr  = hypot(ma, mb);
                    if (rr == 0.0) break;

                    double cc, ssr, ssi;
                    if (ma == 0.0) {
                        cc = 0.0;
                        ssr = bb_r / mb;
                        ssi = bb_i / mb;
                    } else {
                        cc = ma / rr;
                        double den = ma * rr;
                        ssr = (bb_r * aar + bb_i * aai) / den;
                        ssi = (bb_i * aar - bb_r * aai) / den;
                    }

                    banded_givens_two_sided_complex(A_re, A_im, n, new_p, new_q,
                                                      cc, ssr, ssi,
                                                      Q_re, Q_im, want_Q);
                    A_re[new_q * n + bcol] = 0.0; A_im[new_q * n + bcol] = 0.0;
                    A_re[bcol * n + new_q] = 0.0; A_im[bcol * n + new_q] = 0.0;

                    cp = new_p;
                    cq = new_q;
                }
            }
        }
    }

    /* The diagonal is real for a Hermitian matrix; the off-diagonal
     * carries phase that direct_phase_correct_tridiag mops up. */
    for (size_t i = 0; i < n; i++) diag[i] = A_re[i * n + i];
    for (size_t i = 0; i + 1 < n; i++) {
        sub_re[i] = A_re[(i + 1) * n + i];
        sub_im[i] = A_im[(i + 1) * n + i];
    }
}

/* Top-level "Banded" kernel for a complex Hermitian matrix at machine
 * precision.  Mirrors direct_complex_hermitian_machine, swapping the
 * dense Householder tridiagonalisation for the band-Givens reduction. */
static Expr* banded_complex_hermitian_machine(const MatD* A, size_t b,
                                                 MateigenWant want,
                                                 Expr* k_spec) {
    size_t n = A->n;
    if (n == 0) return NULL;

    double* W_re = (double*)malloc(sizeof(double) * n * n);
    double* W_im = (double*)malloc(sizeof(double) * n * n);
    memcpy(W_re, A->re, sizeof(double) * n * n);
    memcpy(W_im, A->im, sizeof(double) * n * n);

    double* diag   = (double*)malloc(sizeof(double) * n);
    double* sub_re = (double*)calloc(n, sizeof(double));
    double* sub_im = (double*)calloc(n, sizeof(double));

    bool want_Q   = (want & MATEIGEN_WANT_VECTORS) != 0;
    double* Q_re  = want_Q ? (double*)malloc(sizeof(double) * n * n) : NULL;
    double* Q_im  = want_Q ? (double*)malloc(sizeof(double) * n * n) : NULL;

    banded_to_tridiag_complex_hermitian(W_re, W_im, n, b,
                                          diag, sub_re, sub_im,
                                          Q_re, Q_im, want_Q);
    free(W_re); free(W_im);

    direct_phase_correct_tridiag(sub_re, sub_im, n, Q_re, Q_im, want_Q);

    double* Z = want_Q ? (double*)malloc(sizeof(double) * n * n) : NULL;
    if (want_Q) {
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++) Z[i * n + j] = (i == j) ? 1.0 : 0.0;
    }
    int qr_status = direct_symtridiag_qr(diag, sub_re, n, Z, want_Q);

    free(sub_re); free(sub_im);

    if (qr_status != 0) {
        free(diag);
        if (Q_re) { free(Q_re); free(Q_im); }
        if (Z) free(Z);
        return NULL;
    }

    size_t* perm = (size_t*)malloc(sizeof(size_t) * n);
    direct_sort_perm_desc_abs(diag, n, perm);

    Expr* out;
    if (want_Q) {
        double* V_re = (double*)malloc(sizeof(double) * n * n);
        double* V_im = (double*)malloc(sizeof(double) * n * n);
        direct_compose_complex_Q_real_Z(Q_re, Q_im, Z, n, V_re, V_im);
        out = direct_build_complex_hermitian_eigvec_list(V_re, V_im, n, perm);
        free(V_re); free(V_im);
    } else {
        out = direct_build_real_eigenvalue_list(diag, n, perm);
    }

    free(diag); free(perm);
    if (Q_re) { free(Q_re); free(Q_im); }
    if (Z) free(Z);

    return direct_apply_k_spec_list(out, k_spec);
}

/* Top-level "Banded" machine dispatcher.  Returns NULL when the input
 * isn't Hermitian (or symmetric), or when the matrix is already dense
 * (b == n - 1) so the caller falls back to Direct without wasting a
 * Givens sweep.  `opts` is unused at present -- see BandedOpts. */
static Expr* banded_dispatch_machine(Expr* m, Expr* a, int64_t n,
                                       MateigenWant want, Expr* k_spec,
                                       const BandedOpts* opts) {
    (void)opts;
    if (a != NULL) return NULL;
    if (n <= 0)    return NULL;

    MatD A;
    if (!matD_load(m, (size_t)n, &A)) return NULL;

    Expr* out = NULL;
    if (A.is_complex) {
        double norm = matD_norm_inf_complex(&A);
        double herm_tol = 1e-12 * (norm == 0.0 ? 1.0 : norm) * (double)A.n;
        if (!matD_is_hermitian(&A, herm_tol)) { matD_free(&A); return NULL; }
        size_t b = matD_bandwidth_complex(&A, herm_tol);
        if (A.n >= 2 && b >= A.n - 1) {
            /* Already dense -- Banded would just churn Givens. */
            matD_free(&A);
            return NULL;
        }
        out = banded_complex_hermitian_machine(&A, b, want, k_spec);
    } else {
        double norm = matD_norm_inf_real(A.re, A.n);
        double sym_tol = 1e-12 * (norm == 0.0 ? 1.0 : norm) * (double)A.n;
        if (!matD_is_real_symmetric(&A, sym_tol)) { matD_free(&A); return NULL; }
        size_t b = matD_bandwidth_real(A.re, A.n, sym_tol);
        if (A.n >= 2 && b >= A.n - 1) {
            matD_free(&A);
            return NULL;
        }
        out = banded_real_sym_machine(&A, b, want, k_spec);
    }
    matD_free(&A);
    return out;
}

/* Automatic dispatch heuristic: prefer Banded when the matrix is
 * Hermitian and the half-bandwidth is small.  Reads the matrix once
 * via matD_load (which is also what banded_dispatch_machine does -- a
 * shared cache would shave one parse but the load is cheap relative
 * to a single Givens sweep, so we don't bother). */
bool banded_automatic_prefers(Expr* m, int64_t n) {
    if (n <= 8) return false;
    MatD A;
    if (!matD_load(m, (size_t)n, &A)) return false;
    bool out = false;
    if (A.is_complex) {
        double norm = matD_norm_inf_complex(&A);
        double herm_tol = 1e-12 * (norm == 0.0 ? 1.0 : norm) * (double)A.n;
        if (matD_is_hermitian(&A, herm_tol)) {
            size_t b = matD_bandwidth_complex(&A, herm_tol);
            out = banded_automatic_prefers_real(&A, b, herm_tol);
        }
    } else {
        double norm = matD_norm_inf_real(A.re, A.n);
        double sym_tol = 1e-12 * (norm == 0.0 ? 1.0 : norm) * (double)A.n;
        if (matD_is_real_symmetric(&A, sym_tol)) {
            size_t b = matD_bandwidth_real(A.re, A.n, sym_tol);
            out = banded_automatic_prefers_real(&A, b, sym_tol);
        }
    }
    matD_free(&A);
    return out;
}
#ifdef USE_MPFR

/* Apply two-sided real Givens G^T A G in plane (p, q) at MPFR
 * precision.  Scratch cells (`s1`, `s2`) are caller-supplied. */
static void banded_givens_two_sided_real_M(mpfr_t* A, size_t n,
                                             size_t p, size_t q,
                                             const mpfr_t c, const mpfr_t s,
                                             mpfr_t* Q, bool want_Q,
                                             mpfr_t s1, mpfr_t s2) {
    /* Row update. */
    for (size_t j = 0; j < n; j++) {
        mpfr_mul(s1, c, A[p * n + j], MPFR_RNDN);
        mpfr_mul(s2, s, A[q * n + j], MPFR_RNDN);
        mpfr_add(s1, s1, s2, MPFR_RNDN);                /* s1 = new A[p, j] */
        mpfr_mul(s2, s, A[p * n + j], MPFR_RNDN);
        mpfr_neg(s2, s2, MPFR_RNDN);                    /* -s * A[p, j] */
        mpfr_set(A[p * n + j], s1, MPFR_RNDN);
        mpfr_mul(s1, c, A[q * n + j], MPFR_RNDN);
        mpfr_add(A[q * n + j], s1, s2, MPFR_RNDN);      /* new A[q, j] */
    }
    /* Col update. */
    for (size_t i = 0; i < n; i++) {
        mpfr_mul(s1, c, A[i * n + p], MPFR_RNDN);
        mpfr_mul(s2, s, A[i * n + q], MPFR_RNDN);
        mpfr_add(s1, s1, s2, MPFR_RNDN);                /* new A[i, p] */
        mpfr_mul(s2, s, A[i * n + p], MPFR_RNDN);
        mpfr_neg(s2, s2, MPFR_RNDN);
        mpfr_set(A[i * n + p], s1, MPFR_RNDN);
        mpfr_mul(s1, c, A[i * n + q], MPFR_RNDN);
        mpfr_add(A[i * n + q], s1, s2, MPFR_RNDN);
    }
    if (want_Q) {
        for (size_t i = 0; i < n; i++) {
            mpfr_mul(s1, c, Q[i * n + p], MPFR_RNDN);
            mpfr_mul(s2, s, Q[i * n + q], MPFR_RNDN);
            mpfr_add(s1, s1, s2, MPFR_RNDN);
            mpfr_mul(s2, s, Q[i * n + p], MPFR_RNDN);
            mpfr_neg(s2, s2, MPFR_RNDN);
            mpfr_set(Q[i * n + p], s1, MPFR_RNDN);
            mpfr_mul(s1, c, Q[i * n + q], MPFR_RNDN);
            mpfr_add(Q[i * n + q], s1, s2, MPFR_RNDN);
        }
    }
}

/* MPFR variant of banded_to_tridiag_real_sym.  `tmp` holds >= 5
 * pre-initialised scratch cells (c, s, r, s1, s2). */
static void banded_to_tridiag_real_sym_M(mpfr_t* A, size_t n, size_t b,
                                           mpfr_t* diag, mpfr_t* sub,
                                           mpfr_t* Q, bool want_Q,
                                           mpfr_t* tmp /* >= 5 cells */) {
    if (want_Q) {
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                mpfr_set_ui(Q[i * n + j], (i == j) ? 1 : 0, MPFR_RNDN);
    }

    mpfr_t* c  = &tmp[0];
    mpfr_t* s  = &tmp[1];
    mpfr_t* r  = &tmp[2];
    mpfr_t* s1 = &tmp[3];
    mpfr_t* s2 = &tmp[4];

    if (n >= 3 && b >= 2) {
        for (size_t k = 0; k + 2 < n; k++) {
            size_t d_max = b;
            if (k + d_max >= n) d_max = n - 1 - k;
            for (size_t d = d_max; d >= 2; d--) {
                size_t p = k + d - 1;
                size_t q = k + d;
                if (mpfr_zero_p(A[q * n + k])) continue;
                mpfr_hypot(*r, A[p * n + k], A[q * n + k], MPFR_RNDN);
                if (mpfr_zero_p(*r)) continue;
                mpfr_div(*c, A[p * n + k], *r, MPFR_RNDN);
                mpfr_div(*s, A[q * n + k], *r, MPFR_RNDN);

                banded_givens_two_sided_real_M(A, n, p, q, *c, *s,
                                                  Q, want_Q, *s1, *s2);
                mpfr_set_zero(A[q * n + k], 1);
                mpfr_set_zero(A[k * n + q], 1);

                size_t cp = p;
                size_t cq = q;
                while (cq + b < n) {
                    size_t bcol  = cp;
                    size_t new_p = cq + b - 1;
                    size_t new_q = cq + b;
                    if (mpfr_zero_p(A[new_q * n + bcol])) break;
                    mpfr_hypot(*r, A[new_p * n + bcol],
                                   A[new_q * n + bcol], MPFR_RNDN);
                    if (mpfr_zero_p(*r)) break;
                    mpfr_div(*c, A[new_p * n + bcol], *r, MPFR_RNDN);
                    mpfr_div(*s, A[new_q * n + bcol], *r, MPFR_RNDN);

                    banded_givens_two_sided_real_M(A, n, new_p, new_q,
                                                      *c, *s, Q, want_Q,
                                                      *s1, *s2);
                    mpfr_set_zero(A[new_q * n + bcol], 1);
                    mpfr_set_zero(A[bcol * n + new_q], 1);

                    cp = new_p;
                    cq = new_q;
                }
            }
        }
    }

    for (size_t i = 0; i < n; i++) mpfr_set(diag[i], A[i * n + i], MPFR_RNDN);
    for (size_t i = 0; i + 1 < n; i++)
        mpfr_set(sub[i], A[(i + 1) * n + i], MPFR_RNDN);
}

/* MPFR orchestrator for real symmetric banded -- mirrors
 * direct_real_sym_mpfr with the band-Givens reduction substituted. */
static Expr* banded_real_sym_mpfr(const MatM* A, size_t b,
                                    MateigenWant want, Expr* k_spec) {
    size_t n = A->n;
    mpfr_prec_t bits = A->bits;
    bool want_Q = (want & MATEIGEN_WANT_VECTORS) != 0;

    if (n == 0) {
        Expr* empty = expr_new_function(expr_new_symbol("List"), NULL, 0);
        return direct_apply_k_spec_list(empty, k_spec);
    }

    mpfr_t* W = mpfr_array_alloc(n * n, bits);
    for (size_t i = 0; i < n * n; i++) mpfr_set(W[i], A->re[i], MPFR_RNDN);

    mpfr_t* diag = mpfr_array_alloc(n, bits);
    mpfr_t* sub  = mpfr_array_alloc(n > 0 ? n : 1, bits);   /* n-1 + slack */
    mpfr_t* Q    = want_Q ? mpfr_array_alloc(n * n, bits) : NULL;
    mpfr_t* tmp_band = mpfr_array_alloc(5,  bits);
    mpfr_t* tmp_qr   = mpfr_array_alloc(12, bits);

    if (n == 1) {
        mpfr_set(diag[0], W[0], MPFR_RNDN);
        if (want_Q) mpfr_set_ui(Q[0], 1, MPFR_RNDN);
    } else {
        banded_to_tridiag_real_sym_M(W, n, b, diag, sub, Q, want_Q, tmp_band);
        direct_symtridiag_qr_M(diag, sub, n, bits, Q, want_Q, tmp_qr);
    }

    size_t* perm = (size_t*)malloc(sizeof(size_t) * n);
    direct_sort_perm_desc_abs_M(diag, n, perm);

    Expr* result = (want_Q)
        ? direct_build_real_eigenvector_list_M(Q, n, perm)
        : direct_build_real_eigenvalue_list_M(diag, n, perm);
    free(perm);

    mpfr_array_free(W,    n * n);
    mpfr_array_free(diag, n);
    mpfr_array_free(sub,  n > 0 ? n : 1);
    if (Q) mpfr_array_free(Q, n * n);
    mpfr_array_free(tmp_band, 5);
    mpfr_array_free(tmp_qr,   12);

    return direct_apply_k_spec_list(result, k_spec);
}

/* Apply two-sided complex Hermitian Givens U^H A U in plane (p, q) at
 * MPFR precision.  Same algebraic formulas as
 * banded_givens_two_sided_complex; each scalar op uses MPFR scratch. */
static void banded_givens_two_sided_complex_M(mpfr_t* A_re, mpfr_t* A_im,
                                                 size_t n, size_t p, size_t q,
                                                 const mpfr_t c,
                                                 const mpfr_t s_re,
                                                 const mpfr_t s_im,
                                                 mpfr_t* Q_re, mpfr_t* Q_im,
                                                 bool want_Q,
                                                 mpfr_t t0, mpfr_t t1,
                                                 mpfr_t t2, mpfr_t t3,
                                                 mpfr_t scratch) {
    /* Row update. */
    for (size_t j = 0; j < n; j++) {
        /* (apr, api) = A[p, j]; (aqr, aqi) = A[q, j] -- snapshot. */
        mpfr_set(t0, A_re[p * n + j], MPFR_RNDN);
        mpfr_set(t1, A_im[p * n + j], MPFR_RNDN);
        mpfr_set(t2, A_re[q * n + j], MPFR_RNDN);
        mpfr_set(t3, A_im[q * n + j], MPFR_RNDN);

        /* new A[p, j]_re = c*apr + s_re*aqr + s_im*aqi */
        mpfr_mul(scratch, c, t0, MPFR_RNDN);
        mpfr_mul(A_re[p * n + j], s_re, t2, MPFR_RNDN);
        mpfr_add(A_re[p * n + j], A_re[p * n + j], scratch, MPFR_RNDN);
        mpfr_mul(scratch, s_im, t3, MPFR_RNDN);
        mpfr_add(A_re[p * n + j], A_re[p * n + j], scratch, MPFR_RNDN);
        /* new A[p, j]_im = c*api + s_re*aqi - s_im*aqr */
        mpfr_mul(scratch, c, t1, MPFR_RNDN);
        mpfr_mul(A_im[p * n + j], s_re, t3, MPFR_RNDN);
        mpfr_add(A_im[p * n + j], A_im[p * n + j], scratch, MPFR_RNDN);
        mpfr_mul(scratch, s_im, t2, MPFR_RNDN);
        mpfr_sub(A_im[p * n + j], A_im[p * n + j], scratch, MPFR_RNDN);
        /* new A[q, j]_re = -(s_re*apr - s_im*api) + c*aqr */
        mpfr_mul(A_re[q * n + j], s_re, t0, MPFR_RNDN);
        mpfr_mul(scratch, s_im, t1, MPFR_RNDN);
        mpfr_sub(A_re[q * n + j], A_re[q * n + j], scratch, MPFR_RNDN);
        mpfr_neg(A_re[q * n + j], A_re[q * n + j], MPFR_RNDN);
        mpfr_mul(scratch, c, t2, MPFR_RNDN);
        mpfr_add(A_re[q * n + j], A_re[q * n + j], scratch, MPFR_RNDN);
        /* new A[q, j]_im = -(s_re*api + s_im*apr) + c*aqi */
        mpfr_mul(A_im[q * n + j], s_re, t1, MPFR_RNDN);
        mpfr_mul(scratch, s_im, t0, MPFR_RNDN);
        mpfr_add(A_im[q * n + j], A_im[q * n + j], scratch, MPFR_RNDN);
        mpfr_neg(A_im[q * n + j], A_im[q * n + j], MPFR_RNDN);
        mpfr_mul(scratch, c, t3, MPFR_RNDN);
        mpfr_add(A_im[q * n + j], A_im[q * n + j], scratch, MPFR_RNDN);
    }
    /* Col update. */
    for (size_t i = 0; i < n; i++) {
        mpfr_set(t0, A_re[i * n + p], MPFR_RNDN);
        mpfr_set(t1, A_im[i * n + p], MPFR_RNDN);
        mpfr_set(t2, A_re[i * n + q], MPFR_RNDN);
        mpfr_set(t3, A_im[i * n + q], MPFR_RNDN);

        /* new A[i, p]_re = c*ipr + s_re*iqr - s_im*iqi */
        mpfr_mul(scratch, c, t0, MPFR_RNDN);
        mpfr_mul(A_re[i * n + p], s_re, t2, MPFR_RNDN);
        mpfr_add(A_re[i * n + p], A_re[i * n + p], scratch, MPFR_RNDN);
        mpfr_mul(scratch, s_im, t3, MPFR_RNDN);
        mpfr_sub(A_re[i * n + p], A_re[i * n + p], scratch, MPFR_RNDN);
        /* new A[i, p]_im = c*ipi + s_re*iqi + s_im*iqr */
        mpfr_mul(scratch, c, t1, MPFR_RNDN);
        mpfr_mul(A_im[i * n + p], s_re, t3, MPFR_RNDN);
        mpfr_add(A_im[i * n + p], A_im[i * n + p], scratch, MPFR_RNDN);
        mpfr_mul(scratch, s_im, t2, MPFR_RNDN);
        mpfr_add(A_im[i * n + p], A_im[i * n + p], scratch, MPFR_RNDN);
        /* new A[i, q]_re = -(s_re*ipr + s_im*ipi) + c*iqr */
        mpfr_mul(A_re[i * n + q], s_re, t0, MPFR_RNDN);
        mpfr_mul(scratch, s_im, t1, MPFR_RNDN);
        mpfr_add(A_re[i * n + q], A_re[i * n + q], scratch, MPFR_RNDN);
        mpfr_neg(A_re[i * n + q], A_re[i * n + q], MPFR_RNDN);
        mpfr_mul(scratch, c, t2, MPFR_RNDN);
        mpfr_add(A_re[i * n + q], A_re[i * n + q], scratch, MPFR_RNDN);
        /* new A[i, q]_im = (s_im*ipr - s_re*ipi) + c*iqi */
        mpfr_mul(A_im[i * n + q], s_im, t0, MPFR_RNDN);
        mpfr_mul(scratch, s_re, t1, MPFR_RNDN);
        mpfr_sub(A_im[i * n + q], A_im[i * n + q], scratch, MPFR_RNDN);
        mpfr_mul(scratch, c, t3, MPFR_RNDN);
        mpfr_add(A_im[i * n + q], A_im[i * n + q], scratch, MPFR_RNDN);
    }
    if (want_Q) {
        for (size_t i = 0; i < n; i++) {
            mpfr_set(t0, Q_re[i * n + p], MPFR_RNDN);
            mpfr_set(t1, Q_im[i * n + p], MPFR_RNDN);
            mpfr_set(t2, Q_re[i * n + q], MPFR_RNDN);
            mpfr_set(t3, Q_im[i * n + q], MPFR_RNDN);

            /* Q col p' = c*Q_p + s*Q_q (same as col-p update above). */
            mpfr_mul(scratch, c, t0, MPFR_RNDN);
            mpfr_mul(Q_re[i * n + p], s_re, t2, MPFR_RNDN);
            mpfr_add(Q_re[i * n + p], Q_re[i * n + p], scratch, MPFR_RNDN);
            mpfr_mul(scratch, s_im, t3, MPFR_RNDN);
            mpfr_sub(Q_re[i * n + p], Q_re[i * n + p], scratch, MPFR_RNDN);
            mpfr_mul(scratch, c, t1, MPFR_RNDN);
            mpfr_mul(Q_im[i * n + p], s_re, t3, MPFR_RNDN);
            mpfr_add(Q_im[i * n + p], Q_im[i * n + p], scratch, MPFR_RNDN);
            mpfr_mul(scratch, s_im, t2, MPFR_RNDN);
            mpfr_add(Q_im[i * n + p], Q_im[i * n + p], scratch, MPFR_RNDN);

            /* Q col q' = -s_conj * Q_p + c*Q_q. */
            mpfr_mul(Q_re[i * n + q], s_re, t0, MPFR_RNDN);
            mpfr_mul(scratch, s_im, t1, MPFR_RNDN);
            mpfr_add(Q_re[i * n + q], Q_re[i * n + q], scratch, MPFR_RNDN);
            mpfr_neg(Q_re[i * n + q], Q_re[i * n + q], MPFR_RNDN);
            mpfr_mul(scratch, c, t2, MPFR_RNDN);
            mpfr_add(Q_re[i * n + q], Q_re[i * n + q], scratch, MPFR_RNDN);
            mpfr_mul(Q_im[i * n + q], s_im, t0, MPFR_RNDN);
            mpfr_mul(scratch, s_re, t1, MPFR_RNDN);
            mpfr_sub(Q_im[i * n + q], Q_im[i * n + q], scratch, MPFR_RNDN);
            mpfr_mul(scratch, c, t3, MPFR_RNDN);
            mpfr_add(Q_im[i * n + q], Q_im[i * n + q], scratch, MPFR_RNDN);
        }
    }
}

/* MPFR variant of banded_to_tridiag_complex_hermitian.  Uses 12
 * scratch cells: 0=c, 1=s_re, 2=s_im, 3=mag_a, 4=mag_b, 5=r,
 * 6..9 = (t0, t1, t2, t3) snapshots, 10=scratch, 11=den. */
static void banded_to_tridiag_complex_hermitian_M(mpfr_t* A_re, mpfr_t* A_im,
                                                     size_t n, size_t b,
                                                     mpfr_t* diag,
                                                     mpfr_t* sub_re,
                                                     mpfr_t* sub_im,
                                                     mpfr_t* Q_re, mpfr_t* Q_im,
                                                     bool want_Q,
                                                     mpfr_t* tmp /* >= 12 */) {
    if (want_Q) {
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++) {
                mpfr_set_ui(Q_re[i * n + j], (i == j) ? 1 : 0, MPFR_RNDN);
                mpfr_set_zero(Q_im[i * n + j], 1);
            }
    }

    mpfr_t* c     = &tmp[0];
    mpfr_t* s_re  = &tmp[1];
    mpfr_t* s_im  = &tmp[2];
    mpfr_t* ma    = &tmp[3];
    mpfr_t* mb    = &tmp[4];
    mpfr_t* rr    = &tmp[5];
    mpfr_t* t0    = &tmp[6];
    mpfr_t* t1    = &tmp[7];
    mpfr_t* t2    = &tmp[8];
    mpfr_t* t3    = &tmp[9];
    mpfr_t* scrat = &tmp[10];
    mpfr_t* den   = &tmp[11];

    if (n >= 3 && b >= 2) {
        for (size_t k = 0; k + 2 < n; k++) {
            size_t d_max = b;
            if (k + d_max >= n) d_max = n - 1 - k;
            for (size_t d = d_max; d >= 2; d--) {
                size_t p = k + d - 1;
                size_t q = k + d;
                mpfr_hypot(*ma, A_re[p * n + k], A_im[p * n + k], MPFR_RNDN);
                mpfr_hypot(*mb, A_re[q * n + k], A_im[q * n + k], MPFR_RNDN);
                if (mpfr_zero_p(*mb)) continue;
                mpfr_hypot(*rr, *ma, *mb, MPFR_RNDN);
                if (mpfr_zero_p(*rr)) continue;

                if (mpfr_zero_p(*ma)) {
                    mpfr_set_zero(*c, 1);
                    mpfr_div(*s_re, A_re[q * n + k], *mb, MPFR_RNDN);
                    mpfr_div(*s_im, A_im[q * n + k], *mb, MPFR_RNDN);
                } else {
                    mpfr_div(*c, *ma, *rr, MPFR_RNDN);
                    /* den = |a| * r */
                    mpfr_mul(*den, *ma, *rr, MPFR_RNDN);
                    /* s_re = (br*ar + bi*ai) / den */
                    mpfr_mul(*s_re, A_re[q * n + k], A_re[p * n + k], MPFR_RNDN);
                    mpfr_mul(*scrat, A_im[q * n + k], A_im[p * n + k], MPFR_RNDN);
                    mpfr_add(*s_re, *s_re, *scrat, MPFR_RNDN);
                    mpfr_div(*s_re, *s_re, *den, MPFR_RNDN);
                    /* s_im = (bi*ar - br*ai) / den */
                    mpfr_mul(*s_im, A_im[q * n + k], A_re[p * n + k], MPFR_RNDN);
                    mpfr_mul(*scrat, A_re[q * n + k], A_im[p * n + k], MPFR_RNDN);
                    mpfr_sub(*s_im, *s_im, *scrat, MPFR_RNDN);
                    mpfr_div(*s_im, *s_im, *den, MPFR_RNDN);
                }

                banded_givens_two_sided_complex_M(A_re, A_im, n, p, q,
                                                     *c, *s_re, *s_im,
                                                     Q_re, Q_im, want_Q,
                                                     *t0, *t1, *t2, *t3, *scrat);
                mpfr_set_zero(A_re[q * n + k], 1); mpfr_set_zero(A_im[q * n + k], 1);
                mpfr_set_zero(A_re[k * n + q], 1); mpfr_set_zero(A_im[k * n + q], 1);

                size_t cp = p;
                size_t cq = q;
                while (cq + b < n) {
                    size_t bcol  = cp;
                    size_t new_p = cq + b - 1;
                    size_t new_q = cq + b;
                    mpfr_hypot(*ma, A_re[new_p * n + bcol],
                                    A_im[new_p * n + bcol], MPFR_RNDN);
                    mpfr_hypot(*mb, A_re[new_q * n + bcol],
                                    A_im[new_q * n + bcol], MPFR_RNDN);
                    if (mpfr_zero_p(*mb)) break;
                    mpfr_hypot(*rr, *ma, *mb, MPFR_RNDN);
                    if (mpfr_zero_p(*rr)) break;

                    if (mpfr_zero_p(*ma)) {
                        mpfr_set_zero(*c, 1);
                        mpfr_div(*s_re, A_re[new_q * n + bcol], *mb, MPFR_RNDN);
                        mpfr_div(*s_im, A_im[new_q * n + bcol], *mb, MPFR_RNDN);
                    } else {
                        mpfr_div(*c, *ma, *rr, MPFR_RNDN);
                        mpfr_mul(*den, *ma, *rr, MPFR_RNDN);
                        mpfr_mul(*s_re, A_re[new_q * n + bcol],
                                          A_re[new_p * n + bcol], MPFR_RNDN);
                        mpfr_mul(*scrat, A_im[new_q * n + bcol],
                                            A_im[new_p * n + bcol], MPFR_RNDN);
                        mpfr_add(*s_re, *s_re, *scrat, MPFR_RNDN);
                        mpfr_div(*s_re, *s_re, *den, MPFR_RNDN);
                        mpfr_mul(*s_im, A_im[new_q * n + bcol],
                                          A_re[new_p * n + bcol], MPFR_RNDN);
                        mpfr_mul(*scrat, A_re[new_q * n + bcol],
                                            A_im[new_p * n + bcol], MPFR_RNDN);
                        mpfr_sub(*s_im, *s_im, *scrat, MPFR_RNDN);
                        mpfr_div(*s_im, *s_im, *den, MPFR_RNDN);
                    }

                    banded_givens_two_sided_complex_M(A_re, A_im, n, new_p, new_q,
                                                         *c, *s_re, *s_im,
                                                         Q_re, Q_im, want_Q,
                                                         *t0, *t1, *t2, *t3, *scrat);
                    mpfr_set_zero(A_re[new_q * n + bcol], 1);
                    mpfr_set_zero(A_im[new_q * n + bcol], 1);
                    mpfr_set_zero(A_re[bcol * n + new_q], 1);
                    mpfr_set_zero(A_im[bcol * n + new_q], 1);

                    cp = new_p;
                    cq = new_q;
                }
            }
        }
    }

    for (size_t i = 0; i < n; i++) mpfr_set(diag[i], A_re[i * n + i], MPFR_RNDN);
    for (size_t i = 0; i + 1 < n; i++) {
        mpfr_set(sub_re[i], A_re[(i + 1) * n + i], MPFR_RNDN);
        mpfr_set(sub_im[i], A_im[(i + 1) * n + i], MPFR_RNDN);
    }
}

/* MPFR orchestrator for complex Hermitian banded. */
static Expr* banded_complex_hermitian_mpfr(const MatM* A, size_t b,
                                              MateigenWant want, Expr* k_spec) {
    size_t n = A->n;
    mpfr_prec_t bits = A->bits;
    if (n == 0) return NULL;

    bool want_Q = (want & MATEIGEN_WANT_VECTORS) != 0;

    mpfr_t* W_re = mpfr_array_alloc(n * n, bits);
    mpfr_t* W_im = mpfr_array_alloc(n * n, bits);
    for (size_t i = 0; i < n * n; i++) {
        mpfr_set(W_re[i], A->re[i], MPFR_RNDN);
        mpfr_set(W_im[i], A->im[i], MPFR_RNDN);
    }

    mpfr_t* diag   = mpfr_array_alloc(n, bits);
    mpfr_t* sub_re = mpfr_array_alloc(n, bits);
    mpfr_t* sub_im = mpfr_array_alloc(n, bits);
    mpfr_t* Q_re   = want_Q ? mpfr_array_alloc(n * n, bits) : NULL;
    mpfr_t* Q_im   = want_Q ? mpfr_array_alloc(n * n, bits) : NULL;
    mpfr_t* tmp    = mpfr_array_alloc(12, bits);

    if (n == 1) {
        mpfr_set(diag[0], W_re[0], MPFR_RNDN);
        if (want_Q) {
            mpfr_set_ui(Q_re[0], 1, MPFR_RNDN);
            mpfr_set_zero(Q_im[0], 1);
        }
    } else {
        banded_to_tridiag_complex_hermitian_M(W_re, W_im, n, b,
                                                 diag, sub_re, sub_im,
                                                 Q_re, Q_im, want_Q, tmp);
        direct_phase_correct_tridiag_M(sub_re, sub_im, n, bits,
                                         Q_re, Q_im, want_Q);
    }

    mpfr_array_free(W_re, n * n);
    mpfr_array_free(W_im, n * n);
    mpfr_array_free(tmp, 12);

    /* Real symmetric tridiag QR (sub_re is now real positive). */
    mpfr_t* Z = want_Q ? mpfr_array_alloc(n * n, bits) : NULL;
    if (want_Q) {
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                mpfr_set_si(Z[i * n + j], (i == j) ? 1 : 0, MPFR_RNDN);
    }
    mpfr_t* tmp_qr = mpfr_array_alloc(12, bits);
    int qr_status = (n >= 2)
        ? direct_symtridiag_qr_M(diag, sub_re, n, bits, Z, want_Q, tmp_qr)
        : 0;
    mpfr_array_free(tmp_qr, 12);
    mpfr_array_free(sub_re, n);
    mpfr_array_free(sub_im, n);

    if (qr_status != 0) {
        mpfr_array_free(diag, n);
        if (Q_re) { mpfr_array_free(Q_re, n * n); mpfr_array_free(Q_im, n * n); }
        if (Z)    mpfr_array_free(Z, n * n);
        return NULL;
    }

    size_t* perm = (size_t*)malloc(sizeof(size_t) * n);
    direct_sort_perm_desc_abs_M(diag, n, perm);

    Expr* result;
    if (want_Q) {
        mpfr_t* V_re = mpfr_array_alloc(n * n, bits);
        mpfr_t* V_im = mpfr_array_alloc(n * n, bits);
        direct_compose_complex_Q_real_Z_M(Q_re, Q_im, Z, n, bits, V_re, V_im);
        result = direct_build_complex_hermitian_eigvec_list_M(V_re, V_im, n, perm);
        mpfr_array_free(V_re, n * n);
        mpfr_array_free(V_im, n * n);
    } else {
        result = direct_build_real_eigenvalue_list_M(diag, n, perm);
    }

    mpfr_array_free(diag, n);
    if (Q_re) { mpfr_array_free(Q_re, n * n); mpfr_array_free(Q_im, n * n); }
    if (Z)    mpfr_array_free(Z, n * n);
    free(perm);

    return direct_apply_k_spec_list(result, k_spec);
}

/* MPFR Banded dispatcher.  Returns NULL when the matrix isn't
 * Hermitian or is already fully dense (b >= n - 1). */
static Expr* banded_dispatch_mpfr(Expr* m, Expr* a, int64_t n,
                                    mpfr_prec_t bits,
                                    MateigenWant want, Expr* k_spec,
                                    const BandedOpts* opts) {
    (void)opts;
    if (a != NULL) return NULL;
    if (n <= 0)    return NULL;

    MatM A;
    if (!matM_load(m, (size_t)n, bits, &A)) return NULL;

    /* Tolerance: n * 2^{-bits+4} * ||A||_inf. */
    mpfr_t norm, tol, factor;
    mpfr_init2(norm,   bits);
    mpfr_init2(tol,    bits);
    mpfr_init2(factor, bits);
    if (A.is_complex) matM_norm_inf_complex(&A, norm);
    else              matM_norm_inf_real(A.re, A.n, bits, norm);
    if (mpfr_zero_p(norm)) mpfr_set_ui(norm, 1, MPFR_RNDN);
    mpfr_set_ui(factor, 1, MPFR_RNDN);
    mpfr_div_2si(factor, factor, (long)bits - 4, MPFR_RNDN);
    mpfr_mul_ui(factor, factor, (unsigned long)A.n, MPFR_RNDN);
    mpfr_mul(tol, norm, factor, MPFR_RNDN);

    Expr* out = NULL;
    if (A.is_complex) {
        if (matM_is_hermitian(&A, tol)) {
            size_t b = matM_bandwidth_complex(&A, tol);
            if (!(A.n >= 2 && b >= A.n - 1)) {
                out = banded_complex_hermitian_mpfr(&A, b, want, k_spec);
            }
        }
    } else {
        if (matM_is_real_symmetric(&A, tol)) {
            size_t b = matM_bandwidth_real(A.re, A.n, tol);
            if (!(A.n >= 2 && b >= A.n - 1)) {
                out = banded_real_sym_mpfr(&A, b, want, k_spec);
            }
        }
    }

    mpfr_clear(norm); mpfr_clear(tol); mpfr_clear(factor);
    matM_free(&A);
    return out;
}
#endif /* USE_MPFR */

/* Top-level "Banded" dispatcher.  Returns NULL when the input matrix
 * is non-Hermitian or already dense (no band structure to exploit),
 * so the caller falls back to Direct or the symbolic path. */
Expr* banded_dispatch(Expr* m, Expr* a, int64_t n,
                               MateigenWant want, Expr* k_spec,
                               Expr* method_value) {
    if (a != NULL) return NULL;
    if (n <= 0)    return NULL;

    BandedOpts opts;
    banded_parse_subopts(method_value, &opts);

    CommonInexactInfo info = common_scan_inexact(m);
    if (info.has_inexact && info.min_bits > 53) {
        Expr* out = banded_dispatch_mpfr(m, a, n,
                                           (mpfr_prec_t)info.min_bits,
                                           want, k_spec, &opts);
        if (out) return out;
    }
    return banded_dispatch_machine(m, a, n, want, k_spec, &opts);
}
