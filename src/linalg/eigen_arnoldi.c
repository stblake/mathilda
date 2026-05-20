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
 *  Phase 3: numerical "Arnoldi" method (machine precision).      *
 *                                                                 *
 *  Builds an m-step Krylov subspace K_m(A, v0) via classical      *
 *  Gram-Schmidt with one re-orthogonalisation pass ("twice is     *
 *  enough"), reducing A to an m x m upper Hessenberg H_m.         *
 *  Eigenvalues of H_m (Ritz values) approximate the extreme       *
 *  eigenvalues of A; Ritz vectors V_m y_i (where y_i diagonalises *
 *  H_m) approximate the corresponding eigenvectors.               *
 *                                                                 *
 *  The H_m diagonalisation is delegated to Phase 2's Francis QR   *
 *  pipeline -- this is the central reuse that lets Arnoldi share  *
 *  all of the deflation, conjugate-pair handling, and Schur back- *
 *  substitution machinery without duplication.                    *
 *                                                                 *
 *  Default basis size m = max(2k, 20) capped at n, where k is the *
 *  requested eigenvalue count (k_spec).  When k_spec is absent,   *
 *  m = n and Arnoldi is mathematically equivalent to a full       *
 *  Hessenberg reduction.                                          *
 *                                                                 *
 *  Sub-options (Method -> {"Arnoldi", ...}):                      *
 *    - "BasisSize"     -> integer m                                *
 *    - "MaxIterations" -> integer (forwarded to inner QR sweeps)   *
 *    - "Tolerance"     -> Real (breakdown tolerance for ||w||)     *
 *                                                                 *
 *  Real and complex inputs go through separate kernels:           *
 *    - arnoldi_real_general_machine: pure-real Arnoldi.            *
 *    - arnoldi_complex_general_machine: paired re/im Arnoldi with  *
 *      a 2mu x 2mu real block-embedding for the H_m diagonal-     *
 *      isation (mirrors direct_complex_general_machine).          *
 *                                                                 *
 *  Hermitian / symmetric inputs flow through these same kernels:  *
 *  H_m comes out (nearly) symmetric / Hermitian for free; we do   *
 *  not yet special-case Lanczos because the test corpus matrices  *
 *  are small (<= 50) so the per-step savings are not worth a      *
 *  second implementation path.                                    *
 * ============================================================ */

/* User-tunable Arnoldi options, parsed from
 *   Method -> {"Arnoldi", "BasisSize" -> m, "MaxIterations" -> N,
 *              "Tolerance" -> t}.
 * Defaults are applied by arnoldi_set_defaults; unspecified fields
 * keep their default values. */
typedef struct {
    int64_t basis_size;       /* m; 0 = auto from k */
    int64_t max_iterations;   /* inner QR sweep cap (forwarded; unused for now) */
    double  tolerance;        /* breakdown tolerance for ||w|| */
    bool    tolerance_given;
} ArnoldiOpts;

static void arnoldi_set_defaults(ArnoldiOpts* o) {
    o->basis_size = 0;
    o->max_iterations = 1000;
    o->tolerance = 0.0;
    o->tolerance_given = false;
}

/* Coerce an option value to double; returns NaN on failure. */
static double arnoldi_coerce_double(Expr* e) {
    if (!e) return NAN;
    if (e->type == EXPR_INTEGER) return (double)e->data.integer;
    if (e->type == EXPR_REAL)    return e->data.real;
    if (e->type == EXPR_MPFR)    return mpfr_get_d(e->data.mpfr, MPFR_RNDN);
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Rational
        && e->data.function.arg_count == 2) {
        Expr* p = e->data.function.args[0];
        Expr* q = e->data.function.args[1];
        if (p->type == EXPR_INTEGER && q->type == EXPR_INTEGER
            && q->data.integer != 0)
            return (double)p->data.integer / (double)q->data.integer;
    }
    return NAN;
}

/* Match either an interned-symbol key or a string key against canonical
 * Arnoldi sub-option names.  Returns 1 if `key` denotes `name`. */
static bool arnoldi_subopt_key_eq(Expr* key, const char* name,
                                    const char* sym_intern) {
    if (!key) return false;
    if (key->type == EXPR_STRING && strcmp(key->data.string, name) == 0)
        return true;
    if (key->type == EXPR_SYMBOL && sym_intern
        && key->data.symbol == sym_intern)
        return true;
    return false;
}

/* Parse Method -> {"Arnoldi", sub-rules...} into opts.  When the RHS is
 * the bare string/symbol "Arnoldi", defaults apply.  Unknown sub-rules
 * are silently ignored. */
static void arnoldi_parse_subopts(Expr* method_value, ArnoldiOpts* opts) {
    arnoldi_set_defaults(opts);
    if (!method_value) return;
    if (method_value->type != EXPR_FUNCTION) return;
    Expr* head = method_value->data.function.head;
    if (head->type != EXPR_SYMBOL || head->data.symbol != SYM_List) return;
    for (size_t i = 1; i < method_value->data.function.arg_count; i++) {
        Expr* rule = method_value->data.function.args[i];
        if (rule->type != EXPR_FUNCTION) continue;
        if (rule->data.function.arg_count != 2) continue;
        Expr* rh = rule->data.function.head;
        if (rh->type != EXPR_SYMBOL) continue;
        if (rh->data.symbol != SYM_Rule
            && rh->data.symbol != SYM_RuleDelayed) continue;
        Expr* key = rule->data.function.args[0];
        Expr* val = rule->data.function.args[1];
        if (arnoldi_subopt_key_eq(key, "BasisSize", SYM_BasisSize)
            && val->type == EXPR_INTEGER) {
            opts->basis_size = val->data.integer;
        } else if (arnoldi_subopt_key_eq(key, "MaxIterations",
                                          SYM_MaxIterations)
                   && val->type == EXPR_INTEGER) {
            opts->max_iterations = val->data.integer;
        } else if (arnoldi_subopt_key_eq(key, "Tolerance", SYM_Tolerance)) {
            double d = arnoldi_coerce_double(val);
            if (!isnan(d)) { opts->tolerance = d; opts->tolerance_given = true; }
        }
    }
}

/* Number of eigenvalues requested by k_spec (0 = "all").  Accepts
 * Integer k, Integer -k, or UpTo[k]; ignores anything else. */
static size_t arnoldi_target_k(Expr* k_spec, size_t n) {
    if (!k_spec) return 0;
    if (k_spec->type == EXPR_INTEGER) {
        int64_t k = k_spec->data.integer;
        if (k < 0) k = -k;
        return ((size_t)k > n) ? n : (size_t)k;
    }
    if (k_spec->type == EXPR_FUNCTION
        && k_spec->data.function.head->type == EXPR_SYMBOL
        && k_spec->data.function.head->data.symbol == SYM_UpTo
        && k_spec->data.function.arg_count == 1
        && k_spec->data.function.args[0]->type == EXPR_INTEGER) {
        int64_t k = k_spec->data.function.args[0]->data.integer;
        if (k < 0) k = 0;
        return ((size_t)k > n) ? n : (size_t)k;
    }
    return 0;
}

/* Default basis size m = max(2k, 20), capped at n.  k = 0 -> m = n. */
static size_t arnoldi_default_basis_size(size_t target_k, size_t n) {
    size_t m;
    if (target_k == 0) {
        m = n;
    } else {
        m = 2 * target_k;
        if (m < 20) m = 20;
    }
    if (m > n) m = n;
    if (m < 1) m = 1;
    return m;
}

/* Heuristic for the Automatic dispatcher: route to Arnoldi when the
 * caller asked for a small k (k <= max(20, n/10)) AND the matrix is
 * large enough for Arnoldi's asymptotic savings to matter.  Returns
 * false when k_spec is absent so Automatic stays on the full-spectrum
 * Direct path.
 *
 * The n > 32 floor is essential: for small matrices, Direct's full
 * Hessenberg + Francis QR is both faster (no Krylov-loss roundoff
 * penalty) and more robust (Arnoldi's deterministic starting vector
 * can be orthogonal to specific eigenspaces -- e.g. anti-symmetric
 * eigenvectors of small symmetric Toeplitz matrices).  Explicit
 * Method -> "Arnoldi" bypasses this guard. */
bool arnoldi_automatic_prefers(Expr* k_spec, size_t n) {
    if (n <= 32) return false;
    size_t k = arnoldi_target_k(k_spec, n);
    if (k == 0) return false;
    size_t thresh = (n / 10 > 20) ? n / 10 : 20;
    return k <= thresh;
}

/* Real general Arnoldi at machine precision.
 *
 * Returns a freshly-allocated List of mu eigenvalues / eigenvectors
 * sorted by descending |lambda|, then trimmed by k_spec.  mu <= m is
 * the achieved basis size (smaller than m only on lucky breakdown).
 *
 * LAPACK-HOOK: the inner matrix-vector products and CGS dots/axpys are
 * the hottest loops here; under USE_LAPACK they would become dgemv +
 * ddot + daxpy calls.  The H_m diagonalisation is the same hook point
 * as direct_real_general_machine -- dhseqr / dorghr / dtrevc3. */
static Expr* arnoldi_real_general_machine(const MatD* A,
                                          MateigenWant want, Expr* k_spec,
                                          const ArnoldiOpts* opts) {
    size_t n = A->n;
    if (n == 0) return NULL;
    bool want_Q = (want & MATEIGEN_WANT_VECTORS) != 0;

    size_t target_k = arnoldi_target_k(k_spec, n);
    size_t m = (opts->basis_size > 0)
        ? ((size_t)opts->basis_size > n ? n : (size_t)opts->basis_size)
        : arnoldi_default_basis_size(target_k, n);
    if (m < 1) m = 1;

    double normA = matD_norm_inf_real(A->re, n);
    double scale = (normA == 0.0) ? 1.0 : normA;
    double tol = opts->tolerance_given ? opts->tolerance
                 : 1e-12 * (double)n * scale;
    double hard_floor = 1e-16 * scale;
    if (tol < hard_floor) tol = hard_floor;

    /* V: (m+1) columns of length n in column-major layout.
     *   V[j*n + i] = (v_j)_i. */
    double* V = (double*)calloc(n * (m + 1), sizeof(double));
    /* H: m x m row-major, upper Hessenberg, zero-initialised. */
    double* H = (double*)calloc(m * m, sizeof(double));
    if (!V || !H) { free(V); free(H); return NULL; }

    /* Seed v_0 = (1, 2, ..., n) / ||(1, 2, ..., n)||.  Deterministic;
     * non-trivial overlap with most eigenvectors.  Linear-ramp is
     * preferred over all-ones because it breaks the symmetric /
     * antisymmetric decomposition typical of symmetric tridiagonal
     * Toeplitz matrices (where ones is orthogonal to every other
     * eigenvector). */
    {
        double s2 = 0.0;
        for (size_t i = 0; i < n; i++) {
            double x = (double)(i + 1);
            s2 += x * x;
        }
        double inv = 1.0 / sqrt(s2);
        for (size_t i = 0; i < n; i++) V[i] = (double)(i + 1) * inv;
    }

    size_t mu = m;
    for (size_t j = 0; j < m; j++) {
        const double* vj = V + j * n;
        double* w = V + (j + 1) * n;

        /* w = A * vj.  A is row-major: A[i*n + k] = A_ik. */
        for (size_t i = 0; i < n; i++) {
            double s = 0.0;
            for (size_t k = 0; k < n; k++) s += A->re[i * n + k] * vj[k];
            w[i] = s;
        }

        /* Classical Gram-Schmidt: H[i,j] = <w, v_i>; w -= H[i,j] v_i. */
        for (size_t i = 0; i <= j; i++) {
            const double* vi = V + i * n;
            double dot = 0.0;
            for (size_t k = 0; k < n; k++) dot += w[k] * vi[k];
            H[i * m + j] = dot;
            for (size_t k = 0; k < n; k++) w[k] -= dot * vi[k];
        }
        /* One re-orthogonalisation pass -- "twice is enough" (Daniel et
         * al. 1976).  CGS loses orthogonality after one pass; the
         * correction restores it to machine precision in practice. */
        for (size_t i = 0; i <= j; i++) {
            const double* vi = V + i * n;
            double r = 0.0;
            for (size_t k = 0; k < n; k++) r += w[k] * vi[k];
            H[i * m + j] += r;
            for (size_t k = 0; k < n; k++) w[k] -= r * vi[k];
        }
        double nrm2 = 0.0;
        for (size_t k = 0; k < n; k++) nrm2 += w[k] * w[k];
        double nrm = sqrt(nrm2);
        if (nrm < tol) {
            /* Lucky breakdown: A has an invariant subspace of dim j+1
             * containing v_0.  The j+1 Ritz values of H_{j+1} are exact
             * eigenvalues of A (in floating point). */
            mu = j + 1;
            for (size_t k = 0; k < n; k++) w[k] = 0.0;
            break;
        }
        for (size_t k = 0; k < n; k++) w[k] /= nrm;
        if (j + 1 < m) H[(j + 1) * m + j] = nrm;
    }

    /* Extract mu x mu H sub-matrix for QR.  direct_qr_real_general
     * modifies H in place; using a copy keeps the Arnoldi factorisation
     * consistent for the (optional) eigenvector lift. */
    double* Hk = (double*)malloc(sizeof(double) * mu * mu);
    for (size_t i = 0; i < mu; i++)
        for (size_t j = 0; j < mu; j++)
            Hk[i * mu + j] = H[i * m + j];

    double* Qh = NULL;
    if (want_Q) {
        Qh = (double*)malloc(sizeof(double) * mu * mu);
        for (size_t i = 0; i < mu; i++)
            for (size_t j = 0; j < mu; j++)
                Qh[i * mu + j] = (i == j) ? 1.0 : 0.0;
    }

    double* eval_re = (double*)calloc(mu, sizeof(double));
    double* eval_im = (double*)calloc(mu, sizeof(double));
    int qr_status = direct_qr_real_general(Hk, mu, eval_re, eval_im, Qh);
    if (qr_status != 0) {
        free(V); free(H); free(Hk); free(eval_re); free(eval_im);
        if (Qh) free(Qh);
        return NULL;
    }

    size_t* perm = (size_t*)malloc(sizeof(size_t) * mu);
    direct_sort_perm_desc_abs_complex(eval_re, eval_im, mu, perm);

    Expr* out;
    if (want_Q) {
        /* Eigenvectors of H_m in the original basis -- schur_compute_eigvecs
         * gives Vh[sp * mu + j] = j-th component of sp-th (sorted) eigvec. */
        double* Vh_re = (double*)malloc(sizeof(double) * mu * mu);
        double* Vh_im = (double*)malloc(sizeof(double) * mu * mu);
        schur_compute_eigvecs(Hk, Qh, mu, eval_re, eval_im, perm,
                                Vh_re, Vh_im);

        /* Lift to A-eigenvectors via the Arnoldi basis:
         *   VA[k][i] = sum_{j=0}^{mu-1} V_m[j][i] * Vh[k][j]
         * with column-major V_m (V[j*n + i]) and row-eigvec Vh. */
        double* VA_re = (double*)calloc(mu * n, sizeof(double));
        double* VA_im = (double*)calloc(mu * n, sizeof(double));
        for (size_t k = 0; k < mu; k++) {
            for (size_t i = 0; i < n; i++) {
                double sr = 0.0, si = 0.0;
                for (size_t j = 0; j < mu; j++) {
                    sr += V[j * n + i] * Vh_re[k * mu + j];
                    si += V[j * n + i] * Vh_im[k * mu + j];
                }
                VA_re[k * n + i] = sr;
                VA_im[k * n + i] = si;
            }
        }
        /* Re-normalise: V_m is orthonormal and Vh is unitary, so the
         * product is theoretically unit-norm; re-normalise to absorb
         * accumulated CGS roundoff. */
        for (size_t k = 0; k < mu; k++) {
            double s = 0.0;
            for (size_t i = 0; i < n; i++) {
                s += VA_re[k * n + i] * VA_re[k * n + i]
                   + VA_im[k * n + i] * VA_im[k * n + i];
            }
            if (s > 0.0) {
                double inv = 1.0 / sqrt(s);
                for (size_t i = 0; i < n; i++) {
                    VA_re[k * n + i] *= inv;
                    VA_im[k * n + i] *= inv;
                }
            }
        }
        out = mateigen_build_complex_eigvec_list_rect(VA_re, VA_im, mu, n);
        free(VA_re); free(VA_im); free(Vh_re); free(Vh_im);
    } else {
        out = direct_build_complex_eigenvalue_list(eval_re, eval_im, mu, perm);
    }

    free(V); free(H); free(Hk); free(eval_re); free(eval_im); free(perm);
    if (Qh) free(Qh);

    return direct_apply_k_spec_list(out, k_spec);
}

/* Complex general Arnoldi at machine precision.
 *
 * Builds a complex Arnoldi factorisation using paired re/im storage
 * (no libmpc dependency).  H_m is a small mu x mu complex upper
 * Hessenberg matrix; we route it through a 2mu x 2mu real block
 * embedding (same trick as direct_complex_general_machine) so the
 * existing real Francis QR + Schur back-substitution pipeline can
 * extract eigenvalues and eigenvectors.
 *
 * Algorithm:
 *   1. Seed v_0 = (1+0i, ..., 1+0i) / sqrt(n).
 *   2. For j = 0..m-1:
 *      a. w = A * v_j (complex matvec).
 *      b. CGS: H[i,j] = <v_i, w>_hermitian = sum conj(v_i[k]) w[k].
 *         w -= H[i,j] v_i.
 *      c. Re-orth: r = <v_i, w>; H[i,j] += r; w -= r v_i.
 *      d. h_{j+1,j} = ||w||_2.
 *      e. Lucky breakdown if h_{j+1,j} < tol.
 *      f. v_{j+1} = w / h_{j+1,j}.
 *   3. Diagonalise H_m via the 2mu x 2mu real block embedding M.
 *   4. Recover H_m eigenpairs via candidate x = (a - d) + i(b + c) and
 *      grouped Gram-Schmidt for repeated eigenvalues (same protocol as
 *      direct_complex_general_machine).
 *   5. Lift to A-eigenvectors via V_m.
 *
 * LAPACK-HOOK: each loop body would become a zgemv + zdotc + zaxpy
 * sequence under USE_LAPACK; H_m diagonalisation maps to zhseqr +
 * ztrevc3. */
static Expr* arnoldi_complex_general_machine(const MatD* A,
                                              MateigenWant want, Expr* k_spec,
                                              const ArnoldiOpts* opts) {
    size_t n = A->n;
    if (n == 0) return NULL;
    bool want_Q = (want & MATEIGEN_WANT_VECTORS) != 0;

    size_t target_k = arnoldi_target_k(k_spec, n);
    size_t m = (opts->basis_size > 0)
        ? ((size_t)opts->basis_size > n ? n : (size_t)opts->basis_size)
        : arnoldi_default_basis_size(target_k, n);
    if (m < 1) m = 1;

    double normA = matD_norm_inf_complex(A);
    double scale = (normA == 0.0) ? 1.0 : normA;
    double tol = opts->tolerance_given ? opts->tolerance
                 : 1e-12 * (double)n * scale;
    double hard_floor = 1e-16 * scale;
    if (tol < hard_floor) tol = hard_floor;

    /* V: (m+1) columns of complex length n.  Paired arrays. */
    double* V_re = (double*)calloc(n * (m + 1), sizeof(double));
    double* V_im = (double*)calloc(n * (m + 1), sizeof(double));
    /* H: m x m complex upper Hessenberg, paired arrays. */
    double* H_re = (double*)calloc(m * m, sizeof(double));
    double* H_im = (double*)calloc(m * m, sizeof(double));
    if (!V_re || !V_im || !H_re || !H_im) {
        free(V_re); free(V_im); free(H_re); free(H_im);
        return NULL;
    }

    /* Seed v_0 = (1, 2, ..., n) / ||(1, 2, ..., n)||  (real, unit norm).
     * Linear ramp -- see arnoldi_real_general_machine for the rationale. */
    {
        double s2 = 0.0;
        for (size_t i = 0; i < n; i++) {
            double x = (double)(i + 1);
            s2 += x * x;
        }
        double inv = 1.0 / sqrt(s2);
        for (size_t i = 0; i < n; i++) V_re[i] = (double)(i + 1) * inv;
    }

    size_t mu = m;
    for (size_t j = 0; j < m; j++) {
        const double* vj_re = V_re + j * n;
        const double* vj_im = V_im + j * n;
        double* w_re = V_re + (j + 1) * n;
        double* w_im = V_im + (j + 1) * n;

        /* w = A * v_j (complex matvec).  A is row-major complex via re/im. */
        for (size_t i = 0; i < n; i++) {
            double sr = 0.0, si = 0.0;
            for (size_t k = 0; k < n; k++) {
                double ar = A->re[i * n + k];
                double ai = A->im[i * n + k];
                double xr = vj_re[k];
                double xi = vj_im[k];
                sr += ar * xr - ai * xi;
                si += ar * xi + ai * xr;
            }
            w_re[i] = sr;
            w_im[i] = si;
        }

        /* CGS: H[i,j] = <v_i, w>_hermitian = sum conj(v_i[k]) * w[k]. */
        for (size_t i = 0; i <= j; i++) {
            const double* vi_re = V_re + i * n;
            const double* vi_im = V_im + i * n;
            double dr = 0.0, di = 0.0;
            for (size_t k = 0; k < n; k++) {
                /* conj(v_i)[k] = (vi_re - i vi_im) */
                dr += vi_re[k] * w_re[k] + vi_im[k] * w_im[k];
                di += vi_re[k] * w_im[k] - vi_im[k] * w_re[k];
            }
            H_re[i * m + j] = dr;
            H_im[i * m + j] = di;
            /* w -= (dr + i di) v_i */
            for (size_t k = 0; k < n; k++) {
                double pr = dr * vi_re[k] - di * vi_im[k];
                double pi = dr * vi_im[k] + di * vi_re[k];
                w_re[k] -= pr;
                w_im[k] -= pi;
            }
        }
        /* Re-orth (twice is enough). */
        for (size_t i = 0; i <= j; i++) {
            const double* vi_re = V_re + i * n;
            const double* vi_im = V_im + i * n;
            double dr = 0.0, di = 0.0;
            for (size_t k = 0; k < n; k++) {
                dr += vi_re[k] * w_re[k] + vi_im[k] * w_im[k];
                di += vi_re[k] * w_im[k] - vi_im[k] * w_re[k];
            }
            H_re[i * m + j] += dr;
            H_im[i * m + j] += di;
            for (size_t k = 0; k < n; k++) {
                double pr = dr * vi_re[k] - di * vi_im[k];
                double pi = dr * vi_im[k] + di * vi_re[k];
                w_re[k] -= pr;
                w_im[k] -= pi;
            }
        }
        /* h_{j+1,j} = ||w||_2 (real, non-negative). */
        double nrm2 = 0.0;
        for (size_t k = 0; k < n; k++) {
            nrm2 += w_re[k] * w_re[k] + w_im[k] * w_im[k];
        }
        double nrm = sqrt(nrm2);
        if (nrm < tol) {
            mu = j + 1;
            for (size_t k = 0; k < n; k++) { w_re[k] = 0.0; w_im[k] = 0.0; }
            break;
        }
        double inv = 1.0 / nrm;
        for (size_t k = 0; k < n; k++) {
            w_re[k] *= inv;
            w_im[k] *= inv;
        }
        if (j + 1 < m) H_re[(j + 1) * m + j] = nrm;
        /* H_im[(j+1)*m + j] stays zero -- sub-diagonal is real non-negative. */
    }

    /* H_m is now a mu x mu complex matrix.  Build 2mu x 2mu real block
     * embedding M = [[Re(H_m), -Im(H_m)], [Im(H_m), Re(H_m)]] and run
     * the real Francis QR pipeline (same trick as direct_complex_general
     * _machine). */
    size_t N = 2 * mu;
    double* Mb = (double*)malloc(sizeof(double) * N * N);
    for (size_t i = 0; i < mu; i++) {
        for (size_t j = 0; j < mu; j++) {
            double r = H_re[i * m + j];
            double s = H_im[i * m + j];
            Mb[i * N + j]                   =  r;
            Mb[i * N + (j + mu)]            = -s;
            Mb[(i + mu) * N + j]            =  s;
            Mb[(i + mu) * N + (j + mu)]     =  r;
        }
    }

    double* Qm = (double*)malloc(sizeof(double) * N * N);
    for (size_t i = 0; i < N; i++)
        for (size_t j = 0; j < N; j++) Qm[i * N + j] = (i == j) ? 1.0 : 0.0;
    double* ub = (double*)malloc(sizeof(double) * N);
    direct_hessenberg_real(Mb, N, ub, Qm);
    free(ub);

    double* M_eval_re = (double*)calloc(N, sizeof(double));
    double* M_eval_im = (double*)calloc(N, sizeof(double));
    int qr_status = direct_qr_real_general(Mb, N, M_eval_re, M_eval_im, Qm);
    if (qr_status != 0) {
        free(V_re); free(V_im); free(H_re); free(H_im);
        free(Mb); free(Qm); free(M_eval_re); free(M_eval_im);
        return NULL;
    }

    /* Eigenvectors of M (in the original basis), one row per eigenvalue. */
    double* M_evec_re = (double*)malloc(sizeof(double) * N * N);
    double* M_evec_im = (double*)malloc(sizeof(double) * N * N);
    size_t* id_perm = (size_t*)malloc(sizeof(size_t) * N);
    for (size_t i = 0; i < N; i++) id_perm[i] = i;
    schur_compute_eigvecs(Mb, Qm, N, M_eval_re, M_eval_im, id_perm,
                            M_evec_re, M_evec_im);
    free(id_perm); free(Mb); free(Qm);

    /* Recover H_m's eigenpairs from M's spectrum (split spec by +J / -J).
     * For each M-eigenvector w (in the embedded basis) of eigenvalue mu,
     * the H_m candidate is x = (a - d) + i (b + c) with
     *   w_re = [a; b], w_im = [c; d].
     * Grouped Gram-Schmidt deduplicates per-eigenvalue (mu, conj mu)
     * pairs and yields exactly m_A(mu) eigenvectors per distinct mu --
     * summing to mu. */
    double spec_norm = 0.0;
    for (size_t i = 0; i < N; i++) {
        double a = hypot(M_eval_re[i], M_eval_im[i]);
        if (a > spec_norm) spec_norm = a;
    }
    double group_tol = 1e-8 * (spec_norm == 0.0 ? 1.0 : spec_norm) * (double)N;
    double extract_threshold = sqrt((double)mu) * 1e-9;

    int* used = (int*)calloc(N, sizeof(int));
    double* HA_eval_re = (double*)malloc(sizeof(double) * mu);
    double* HA_eval_im = (double*)malloc(sizeof(double) * mu);
    double* HA_evec_re = (double*)calloc(mu * mu, sizeof(double));
    double* HA_evec_im = (double*)calloc(mu * mu, sizeof(double));
    double* cand_re = (double*)malloc(sizeof(double) * mu);
    double* cand_im = (double*)malloc(sizeof(double) * mu);
    size_t produced = 0;

    for (size_t i = 0; i < N && produced < mu; i++) {
        if (used[i]) continue;
        used[i] = 1;
        size_t group_start = produced;
        for (size_t jj = i; jj < N && produced < mu; jj++) {
            if (jj != i) {
                if (used[jj]) continue;
                double dr = M_eval_re[jj] - M_eval_re[i];
                double di = M_eval_im[jj] - M_eval_im[i];
                if (hypot(dr, di) > group_tol) continue;
                used[jj] = 1;
            }
            for (size_t l = 0; l < mu; l++) {
                double a = M_evec_re[jj * N + l];
                double b = M_evec_re[jj * N + (l + mu)];
                double c = M_evec_im[jj * N + l];
                double d_im = M_evec_im[jj * N + (l + mu)];
                cand_re[l] = a - d_im;
                cand_im[l] = b + c;
            }
            /* Orthogonalise against already-emitted vectors in group. */
            for (int pass = 0; pass < 2; pass++) {
                for (size_t f = group_start; f < produced; f++) {
                    double pr = 0.0, pi = 0.0;
                    for (size_t l = 0; l < mu; l++) {
                        double vr = HA_evec_re[f * mu + l];
                        double vi = HA_evec_im[f * mu + l];
                        pr += vr * cand_re[l] + vi * cand_im[l];
                        pi += vr * cand_im[l] - vi * cand_re[l];
                    }
                    for (size_t l = 0; l < mu; l++) {
                        double vr = HA_evec_re[f * mu + l];
                        double vi = HA_evec_im[f * mu + l];
                        double pvr = pr * vr - pi * vi;
                        double pvi = pr * vi + pi * vr;
                        cand_re[l] -= pvr;
                        cand_im[l] -= pvi;
                    }
                }
            }
            double nrm2 = 0.0;
            for (size_t l = 0; l < mu; l++) {
                nrm2 += cand_re[l] * cand_re[l] + cand_im[l] * cand_im[l];
            }
            if (nrm2 < extract_threshold * extract_threshold) continue;
            double inv = 1.0 / sqrt(nrm2);
            for (size_t l = 0; l < mu; l++) {
                HA_evec_re[produced * mu + l] = cand_re[l] * inv;
                HA_evec_im[produced * mu + l] = cand_im[l] * inv;
            }
            HA_eval_re[produced] = M_eval_re[i];
            HA_eval_im[produced] = M_eval_im[i];
            produced++;
        }
    }
    free(used); free(M_eval_re); free(M_eval_im);
    free(M_evec_re); free(M_evec_im);
    free(cand_re); free(cand_im);

    if (produced != mu) {
        free(V_re); free(V_im); free(H_re); free(H_im);
        free(HA_eval_re); free(HA_eval_im);
        free(HA_evec_re); free(HA_evec_im);
        return NULL;
    }

    /* Sort H_m eigenvalues by descending |lambda|. */
    size_t* perm = (size_t*)malloc(sizeof(size_t) * mu);
    direct_sort_perm_desc_abs_complex(HA_eval_re, HA_eval_im, mu, perm);

    Expr* out;
    if (want_Q) {
        /* Lift H_m eigvecs to A eigvecs:
         *   VA[k][i] = sum_{j=0}^{mu-1} V_m[j][i] * y_k[j]
         * where y_k is the perm[k]-th H_m eigenvector. */
        double* VA_re = (double*)calloc(mu * n, sizeof(double));
        double* VA_im = (double*)calloc(mu * n, sizeof(double));
        for (size_t k = 0; k < mu; k++) {
            size_t src = perm[k];
            for (size_t i = 0; i < n; i++) {
                double sr = 0.0, si = 0.0;
                for (size_t j = 0; j < mu; j++) {
                    double vr_basis_re = V_re[j * n + i];
                    double vr_basis_im = V_im[j * n + i];
                    double yr = HA_evec_re[src * mu + j];
                    double yi = HA_evec_im[src * mu + j];
                    /* (vr + i vi)(yr + i yi) */
                    sr += vr_basis_re * yr - vr_basis_im * yi;
                    si += vr_basis_re * yi + vr_basis_im * yr;
                }
                VA_re[k * n + i] = sr;
                VA_im[k * n + i] = si;
            }
        }
        /* Re-normalise. */
        for (size_t k = 0; k < mu; k++) {
            double s = 0.0;
            for (size_t i = 0; i < n; i++) {
                s += VA_re[k * n + i] * VA_re[k * n + i]
                   + VA_im[k * n + i] * VA_im[k * n + i];
            }
            if (s > 0.0) {
                double inv = 1.0 / sqrt(s);
                for (size_t i = 0; i < n; i++) {
                    VA_re[k * n + i] *= inv;
                    VA_im[k * n + i] *= inv;
                }
            }
        }
        out = mateigen_build_complex_eigvec_list_rect(VA_re, VA_im, mu, n);
        free(VA_re); free(VA_im);
    } else {
        out = direct_build_complex_eigenvalue_list(HA_eval_re, HA_eval_im,
                                                     mu, perm);
    }

    free(V_re); free(V_im); free(H_re); free(H_im);
    free(HA_eval_re); free(HA_eval_im);
    free(HA_evec_re); free(HA_evec_im);
    free(perm);

    return direct_apply_k_spec_list(out, k_spec);
}

/* Top-level machine-precision Arnoldi dispatcher. */
static Expr* arnoldi_dispatch_machine(Expr* m, Expr* a, int64_t n,
                                       MateigenWant want, Expr* k_spec,
                                       const ArnoldiOpts* opts) {
    if (a != NULL) return NULL;          /* generalised: symbolic only */
    if (n <= 0)    return NULL;

    MatD A;
    if (!matD_load(m, (size_t)n, &A)) return NULL;

    Expr* out;
    if (A.is_complex) {
        out = arnoldi_complex_general_machine(&A, want, k_spec, opts);
    } else {
        out = arnoldi_real_general_machine(&A, want, k_spec, opts);
    }
    matD_free(&A);
    return out;
}
#ifdef USE_MPFR

/* Real general Arnoldi at MPFR precision. */
static Expr* arnoldi_real_general_mpfr(const MatM* A,
                                        MateigenWant want, Expr* k_spec,
                                        const ArnoldiOpts* opts) {
    size_t n = A->n;
    mpfr_prec_t bits = A->bits;
    if (n == 0) return NULL;
    bool want_Q = (want & MATEIGEN_WANT_VECTORS) != 0;

    size_t target_k = arnoldi_target_k(k_spec, n);
    size_t m = (opts->basis_size > 0)
        ? ((size_t)opts->basis_size > n ? n : (size_t)opts->basis_size)
        : arnoldi_default_basis_size(target_k, n);
    if (m < 1) m = 1;

    /* Tolerance.  Default ~ n * ||A|| / 2^(bits-10): sqrt-eps-like floor
     * because Arnoldi's CGS accumulates more roundoff than Householder. */
    mpfr_t normA, tol;
    mpfr_init2(normA, bits);
    mpfr_init2(tol, bits);
    matM_norm_inf_real(A->re, n, bits, normA);
    if (mpfr_zero_p(normA)) mpfr_set_ui(normA, 1, MPFR_RNDN);
    if (opts->tolerance_given) {
        mpfr_set_d(tol, opts->tolerance, MPFR_RNDN);
    } else {
        mpfr_set(tol, normA, MPFR_RNDN);
        mpfr_mul_ui(tol, tol, (unsigned long)n, MPFR_RNDN);
        mpfr_div_2si(tol, tol, (long)bits - 10, MPFR_RNDN);
    }

    /* V (column-major n*(m+1)) and H (row-major m*m). */
    mpfr_t* V = mpfr_array_alloc(n * (m + 1), bits);
    mpfr_t* H = mpfr_array_alloc(m * m, bits);
    for (size_t i = 0; i < n * (m + 1); i++) mpfr_set_zero(V[i], 1);
    for (size_t i = 0; i < m * m; i++)       mpfr_set_zero(H[i], 1);

    /* Seed v_0 = (1, 2, ..., n) / ||(1, 2, ..., n)||. */
    {
        mpfr_t s2, x, inv;
        mpfr_init2(s2, bits); mpfr_init2(x, bits); mpfr_init2(inv, bits);
        mpfr_set_zero(s2, 1);
        for (size_t i = 0; i < n; i++) {
            mpfr_set_ui(x, (unsigned long)(i + 1), MPFR_RNDN);
            mpfr_mul(x, x, x, MPFR_RNDN);
            mpfr_add(s2, s2, x, MPFR_RNDN);
        }
        mpfr_sqrt(s2, s2, MPFR_RNDN);
        mpfr_ui_div(inv, 1, s2, MPFR_RNDN);
        for (size_t i = 0; i < n; i++)
            mpfr_mul_ui(V[i], inv, (unsigned long)(i + 1), MPFR_RNDN);
        mpfr_clear(s2); mpfr_clear(x); mpfr_clear(inv);
    }

    /* Pre-initialised scratch. */
    mpfr_t sum, dot, r, nrm2, nrm, prod;
    mpfr_init2(sum, bits);  mpfr_init2(dot, bits);
    mpfr_init2(r, bits);    mpfr_init2(nrm2, bits);
    mpfr_init2(nrm, bits);  mpfr_init2(prod, bits);

    size_t mu = m;
    for (size_t j = 0; j < m; j++) {
        mpfr_t* vj = V + j * n;
        mpfr_t* w  = V + (j + 1) * n;

        /* w = A * vj */
        for (size_t i = 0; i < n; i++) {
            mpfr_set_zero(sum, 1);
            for (size_t k = 0; k < n; k++) {
                mpfr_mul(prod, A->re[i * n + k], vj[k], MPFR_RNDN);
                mpfr_add(sum, sum, prod, MPFR_RNDN);
            }
            mpfr_set(w[i], sum, MPFR_RNDN);
        }
        /* CGS */
        for (size_t i = 0; i <= j; i++) {
            mpfr_t* vi = V + i * n;
            mpfr_set_zero(dot, 1);
            for (size_t k = 0; k < n; k++) {
                mpfr_mul(prod, w[k], vi[k], MPFR_RNDN);
                mpfr_add(dot, dot, prod, MPFR_RNDN);
            }
            mpfr_set(H[i * m + j], dot, MPFR_RNDN);
            for (size_t k = 0; k < n; k++) {
                mpfr_mul(prod, dot, vi[k], MPFR_RNDN);
                mpfr_sub(w[k], w[k], prod, MPFR_RNDN);
            }
        }
        /* Re-orth */
        for (size_t i = 0; i <= j; i++) {
            mpfr_t* vi = V + i * n;
            mpfr_set_zero(r, 1);
            for (size_t k = 0; k < n; k++) {
                mpfr_mul(prod, w[k], vi[k], MPFR_RNDN);
                mpfr_add(r, r, prod, MPFR_RNDN);
            }
            mpfr_add(H[i * m + j], H[i * m + j], r, MPFR_RNDN);
            for (size_t k = 0; k < n; k++) {
                mpfr_mul(prod, r, vi[k], MPFR_RNDN);
                mpfr_sub(w[k], w[k], prod, MPFR_RNDN);
            }
        }
        /* h_{j+1,j} = ||w|| */
        mpfr_set_zero(nrm2, 1);
        for (size_t k = 0; k < n; k++) {
            mpfr_mul(prod, w[k], w[k], MPFR_RNDN);
            mpfr_add(nrm2, nrm2, prod, MPFR_RNDN);
        }
        mpfr_sqrt(nrm, nrm2, MPFR_RNDN);
        if (mpfr_cmp(nrm, tol) < 0) {
            mu = j + 1;
            for (size_t k = 0; k < n; k++) mpfr_set_zero(w[k], 1);
            break;
        }
        for (size_t k = 0; k < n; k++) mpfr_div(w[k], w[k], nrm, MPFR_RNDN);
        if (j + 1 < m) mpfr_set(H[(j + 1) * m + j], nrm, MPFR_RNDN);
    }

    /* Diagonalise H_m via the Phase 2d-B Francis QR at MPFR precision. */
    mpfr_t* Hk = mpfr_array_alloc(mu * mu, bits);
    for (size_t i = 0; i < mu; i++)
        for (size_t j = 0; j < mu; j++)
            mpfr_set(Hk[i * mu + j], H[i * m + j], MPFR_RNDN);
    mpfr_t* Qh = NULL;
    if (want_Q) {
        Qh = mpfr_array_alloc(mu * mu, bits);
        for (size_t i = 0; i < mu; i++)
            for (size_t j = 0; j < mu; j++)
                mpfr_set_si(Qh[i * mu + j], (i == j) ? 1 : 0, MPFR_RNDN);
    }
    mpfr_t* eval_re = mpfr_array_alloc(mu, bits);
    mpfr_t* eval_im = mpfr_array_alloc(mu, bits);
    for (size_t i = 0; i < mu; i++) {
        mpfr_set_zero(eval_re[i], 1);
        mpfr_set_zero(eval_im[i], 1);
    }
    mpfr_t* tmp = mpfr_array_alloc(14, bits);
    int qr_status = direct_qr_real_general_M(Hk, mu, bits,
                                              eval_re, eval_im, Qh, tmp);
    if (qr_status != 0) {
        mpfr_clear(sum); mpfr_clear(dot); mpfr_clear(r);
        mpfr_clear(nrm2); mpfr_clear(nrm); mpfr_clear(prod);
        mpfr_clear(normA); mpfr_clear(tol);
        mpfr_array_free(V, n * (m + 1));
        mpfr_array_free(H, m * m);
        mpfr_array_free(Hk, mu * mu);
        if (Qh) mpfr_array_free(Qh, mu * mu);
        mpfr_array_free(eval_re, mu);
        mpfr_array_free(eval_im, mu);
        mpfr_array_free(tmp, 14);
        return NULL;
    }

    size_t* perm = (size_t*)malloc(sizeof(size_t) * mu);
    direct_sort_perm_desc_abs_complex_M(eval_re, eval_im, mu, perm);

    Expr* out;
    if (want_Q) {
        mpfr_t* Vh_re = mpfr_array_alloc(mu * mu, bits);
        mpfr_t* Vh_im = mpfr_array_alloc(mu * mu, bits);
        schur_compute_eigvecs_M(Hk, Qh, mu, bits,
                                  eval_re, eval_im, perm, Vh_re, Vh_im);

        /* Lift: VA[k*n + i] = sum_j V[j*n + i] * Vh[k*mu + j]. */
        mpfr_t* VA_re = mpfr_array_alloc(mu * n, bits);
        mpfr_t* VA_im = mpfr_array_alloc(mu * n, bits);
        for (size_t k = 0; k < mu; k++) {
            for (size_t i = 0; i < n; i++) {
                mpfr_set_zero(sum, 1);
                mpfr_set_zero(dot, 1);
                for (size_t j = 0; j < mu; j++) {
                    mpfr_mul(prod, V[j * n + i], Vh_re[k * mu + j], MPFR_RNDN);
                    mpfr_add(sum, sum, prod, MPFR_RNDN);
                    mpfr_mul(prod, V[j * n + i], Vh_im[k * mu + j], MPFR_RNDN);
                    mpfr_add(dot, dot, prod, MPFR_RNDN);
                }
                mpfr_set(VA_re[k * n + i], sum, MPFR_RNDN);
                mpfr_set(VA_im[k * n + i], dot, MPFR_RNDN);
            }
        }
        /* Re-normalise. */
        for (size_t k = 0; k < mu; k++) {
            mpfr_set_zero(nrm2, 1);
            for (size_t i = 0; i < n; i++) {
                mpfr_mul(prod, VA_re[k * n + i], VA_re[k * n + i], MPFR_RNDN);
                mpfr_add(nrm2, nrm2, prod, MPFR_RNDN);
                mpfr_mul(prod, VA_im[k * n + i], VA_im[k * n + i], MPFR_RNDN);
                mpfr_add(nrm2, nrm2, prod, MPFR_RNDN);
            }
            if (mpfr_sgn(nrm2) > 0) {
                mpfr_sqrt(nrm, nrm2, MPFR_RNDN);
                for (size_t i = 0; i < n; i++) {
                    mpfr_div(VA_re[k * n + i], VA_re[k * n + i], nrm, MPFR_RNDN);
                    mpfr_div(VA_im[k * n + i], VA_im[k * n + i], nrm, MPFR_RNDN);
                }
            }
        }
        out = mateigen_build_complex_eigvec_list_rect_M(VA_re, VA_im, mu, n);
        mpfr_array_free(VA_re, mu * n);
        mpfr_array_free(VA_im, mu * n);
        mpfr_array_free(Vh_re, mu * mu);
        mpfr_array_free(Vh_im, mu * mu);
    } else {
        out = direct_build_complex_eigenvalue_list_M(eval_re, eval_im, mu, perm);
    }

    free(perm);
    mpfr_clear(sum); mpfr_clear(dot); mpfr_clear(r);
    mpfr_clear(nrm2); mpfr_clear(nrm); mpfr_clear(prod);
    mpfr_clear(normA); mpfr_clear(tol);
    mpfr_array_free(V, n * (m + 1));
    mpfr_array_free(H, m * m);
    mpfr_array_free(Hk, mu * mu);
    if (Qh) mpfr_array_free(Qh, mu * mu);
    mpfr_array_free(eval_re, mu);
    mpfr_array_free(eval_im, mu);
    mpfr_array_free(tmp, 14);

    return direct_apply_k_spec_list(out, k_spec);
}

/* Complex general Arnoldi at MPFR precision.  Mirrors the machine
 * kernel: paired re/im V_m and H_m via mpfr_t; H_m diagonalised
 * through a 2mu x 2mu real block embedding routed through
 * direct_qr_real_general_M; grouped complex Gram-Schmidt at MPFR
 * precision deduplicates the doubled spectrum. */
static Expr* arnoldi_complex_general_mpfr(const MatM* A,
                                            MateigenWant want, Expr* k_spec,
                                            const ArnoldiOpts* opts) {
    size_t n = A->n;
    mpfr_prec_t bits = A->bits;
    if (n == 0) return NULL;
    bool want_Q = (want & MATEIGEN_WANT_VECTORS) != 0;

    size_t target_k = arnoldi_target_k(k_spec, n);
    size_t m = (opts->basis_size > 0)
        ? ((size_t)opts->basis_size > n ? n : (size_t)opts->basis_size)
        : arnoldi_default_basis_size(target_k, n);
    if (m < 1) m = 1;

    mpfr_t normA, tol;
    mpfr_init2(normA, bits);
    mpfr_init2(tol, bits);
    matM_norm_inf_complex(A, normA);
    if (mpfr_zero_p(normA)) mpfr_set_ui(normA, 1, MPFR_RNDN);
    if (opts->tolerance_given) {
        mpfr_set_d(tol, opts->tolerance, MPFR_RNDN);
    } else {
        mpfr_set(tol, normA, MPFR_RNDN);
        mpfr_mul_ui(tol, tol, (unsigned long)n, MPFR_RNDN);
        mpfr_div_2si(tol, tol, (long)bits - 10, MPFR_RNDN);
    }

    /* Paired V (column-major) and H (row-major). */
    mpfr_t* V_re = mpfr_array_alloc(n * (m + 1), bits);
    mpfr_t* V_im = mpfr_array_alloc(n * (m + 1), bits);
    mpfr_t* H_re = mpfr_array_alloc(m * m, bits);
    mpfr_t* H_im = mpfr_array_alloc(m * m, bits);
    for (size_t i = 0; i < n * (m + 1); i++) {
        mpfr_set_zero(V_re[i], 1);
        mpfr_set_zero(V_im[i], 1);
    }
    for (size_t i = 0; i < m * m; i++) {
        mpfr_set_zero(H_re[i], 1);
        mpfr_set_zero(H_im[i], 1);
    }

    /* Seed v_0 = (1, 2, ..., n) / ||(1, 2, ..., n)||. */
    {
        mpfr_t s2, x, inv;
        mpfr_init2(s2, bits); mpfr_init2(x, bits); mpfr_init2(inv, bits);
        mpfr_set_zero(s2, 1);
        for (size_t i = 0; i < n; i++) {
            mpfr_set_ui(x, (unsigned long)(i + 1), MPFR_RNDN);
            mpfr_mul(x, x, x, MPFR_RNDN);
            mpfr_add(s2, s2, x, MPFR_RNDN);
        }
        mpfr_sqrt(s2, s2, MPFR_RNDN);
        mpfr_ui_div(inv, 1, s2, MPFR_RNDN);
        for (size_t i = 0; i < n; i++)
            mpfr_mul_ui(V_re[i], inv, (unsigned long)(i + 1), MPFR_RNDN);
        mpfr_clear(s2); mpfr_clear(x); mpfr_clear(inv);
    }

    mpfr_t sum_r, sum_i, dot_r, dot_i, pr, pi, prod, q;
    mpfr_t nrm2, nrm;
    mpfr_init2(sum_r, bits); mpfr_init2(sum_i, bits);
    mpfr_init2(dot_r, bits); mpfr_init2(dot_i, bits);
    mpfr_init2(pr, bits);    mpfr_init2(pi, bits);
    mpfr_init2(prod, bits);  mpfr_init2(q, bits);
    mpfr_init2(nrm2, bits);  mpfr_init2(nrm, bits);

    size_t mu = m;
    for (size_t j = 0; j < m; j++) {
        mpfr_t* vj_r = V_re + j * n;
        mpfr_t* vj_i = V_im + j * n;
        mpfr_t* w_r  = V_re + (j + 1) * n;
        mpfr_t* w_i  = V_im + (j + 1) * n;

        /* w = A * vj  (complex matvec) */
        for (size_t i = 0; i < n; i++) {
            mpfr_set_zero(sum_r, 1);
            mpfr_set_zero(sum_i, 1);
            for (size_t k = 0; k < n; k++) {
                /* (ar + i ai)(xr + i xi) */
                mpfr_mul(prod, A->re[i * n + k], vj_r[k], MPFR_RNDN);
                mpfr_add(sum_r, sum_r, prod, MPFR_RNDN);
                mpfr_mul(prod, A->im[i * n + k], vj_i[k], MPFR_RNDN);
                mpfr_sub(sum_r, sum_r, prod, MPFR_RNDN);
                mpfr_mul(prod, A->re[i * n + k], vj_i[k], MPFR_RNDN);
                mpfr_add(sum_i, sum_i, prod, MPFR_RNDN);
                mpfr_mul(prod, A->im[i * n + k], vj_r[k], MPFR_RNDN);
                mpfr_add(sum_i, sum_i, prod, MPFR_RNDN);
            }
            mpfr_set(w_r[i], sum_r, MPFR_RNDN);
            mpfr_set(w_i[i], sum_i, MPFR_RNDN);
        }

        /* CGS: H[i,j] = sum conj(v_i) * w */
        for (size_t i = 0; i <= j; i++) {
            mpfr_t* vi_r = V_re + i * n;
            mpfr_t* vi_i = V_im + i * n;
            mpfr_set_zero(dot_r, 1); mpfr_set_zero(dot_i, 1);
            for (size_t k = 0; k < n; k++) {
                /* conj(v_i)[k] * w[k] = (vi_r - i vi_i)(w_r + i w_i) */
                mpfr_mul(prod, vi_r[k], w_r[k], MPFR_RNDN);
                mpfr_add(dot_r, dot_r, prod, MPFR_RNDN);
                mpfr_mul(prod, vi_i[k], w_i[k], MPFR_RNDN);
                mpfr_add(dot_r, dot_r, prod, MPFR_RNDN);
                mpfr_mul(prod, vi_r[k], w_i[k], MPFR_RNDN);
                mpfr_add(dot_i, dot_i, prod, MPFR_RNDN);
                mpfr_mul(prod, vi_i[k], w_r[k], MPFR_RNDN);
                mpfr_sub(dot_i, dot_i, prod, MPFR_RNDN);
            }
            mpfr_set(H_re[i * m + j], dot_r, MPFR_RNDN);
            mpfr_set(H_im[i * m + j], dot_i, MPFR_RNDN);
            /* w -= H[i,j] * v_i */
            for (size_t k = 0; k < n; k++) {
                /* (dot_r + i dot_i)(vi_r + i vi_i) */
                mpfr_mul(pr, dot_r, vi_r[k], MPFR_RNDN);
                mpfr_mul(prod, dot_i, vi_i[k], MPFR_RNDN);
                mpfr_sub(pr, pr, prod, MPFR_RNDN);
                mpfr_mul(pi, dot_r, vi_i[k], MPFR_RNDN);
                mpfr_mul(prod, dot_i, vi_r[k], MPFR_RNDN);
                mpfr_add(pi, pi, prod, MPFR_RNDN);
                mpfr_sub(w_r[k], w_r[k], pr, MPFR_RNDN);
                mpfr_sub(w_i[k], w_i[k], pi, MPFR_RNDN);
            }
        }
        /* Re-orth */
        for (size_t i = 0; i <= j; i++) {
            mpfr_t* vi_r = V_re + i * n;
            mpfr_t* vi_i = V_im + i * n;
            mpfr_set_zero(dot_r, 1); mpfr_set_zero(dot_i, 1);
            for (size_t k = 0; k < n; k++) {
                mpfr_mul(prod, vi_r[k], w_r[k], MPFR_RNDN);
                mpfr_add(dot_r, dot_r, prod, MPFR_RNDN);
                mpfr_mul(prod, vi_i[k], w_i[k], MPFR_RNDN);
                mpfr_add(dot_r, dot_r, prod, MPFR_RNDN);
                mpfr_mul(prod, vi_r[k], w_i[k], MPFR_RNDN);
                mpfr_add(dot_i, dot_i, prod, MPFR_RNDN);
                mpfr_mul(prod, vi_i[k], w_r[k], MPFR_RNDN);
                mpfr_sub(dot_i, dot_i, prod, MPFR_RNDN);
            }
            mpfr_add(H_re[i * m + j], H_re[i * m + j], dot_r, MPFR_RNDN);
            mpfr_add(H_im[i * m + j], H_im[i * m + j], dot_i, MPFR_RNDN);
            for (size_t k = 0; k < n; k++) {
                mpfr_mul(pr, dot_r, vi_r[k], MPFR_RNDN);
                mpfr_mul(prod, dot_i, vi_i[k], MPFR_RNDN);
                mpfr_sub(pr, pr, prod, MPFR_RNDN);
                mpfr_mul(pi, dot_r, vi_i[k], MPFR_RNDN);
                mpfr_mul(prod, dot_i, vi_r[k], MPFR_RNDN);
                mpfr_add(pi, pi, prod, MPFR_RNDN);
                mpfr_sub(w_r[k], w_r[k], pr, MPFR_RNDN);
                mpfr_sub(w_i[k], w_i[k], pi, MPFR_RNDN);
            }
        }

        /* h = ||w||_2 (real, non-negative) */
        mpfr_set_zero(nrm2, 1);
        for (size_t k = 0; k < n; k++) {
            mpfr_mul(prod, w_r[k], w_r[k], MPFR_RNDN);
            mpfr_add(nrm2, nrm2, prod, MPFR_RNDN);
            mpfr_mul(prod, w_i[k], w_i[k], MPFR_RNDN);
            mpfr_add(nrm2, nrm2, prod, MPFR_RNDN);
        }
        mpfr_sqrt(nrm, nrm2, MPFR_RNDN);
        if (mpfr_cmp(nrm, tol) < 0) {
            mu = j + 1;
            for (size_t k = 0; k < n; k++) {
                mpfr_set_zero(w_r[k], 1);
                mpfr_set_zero(w_i[k], 1);
            }
            break;
        }
        for (size_t k = 0; k < n; k++) {
            mpfr_div(w_r[k], w_r[k], nrm, MPFR_RNDN);
            mpfr_div(w_i[k], w_i[k], nrm, MPFR_RNDN);
        }
        if (j + 1 < m) mpfr_set(H_re[(j + 1) * m + j], nrm, MPFR_RNDN);
    }

    /* Build 2mu x 2mu real block embedding M and run Phase 2d-B QR. */
    size_t N = 2 * mu;
    mpfr_t* Mb = mpfr_array_alloc(N * N, bits);
    for (size_t i = 0; i < mu; i++) {
        for (size_t j = 0; j < mu; j++) {
            mpfr_set(Mb[i * N + j], H_re[i * m + j], MPFR_RNDN);
            mpfr_neg(Mb[i * N + (j + mu)], H_im[i * m + j], MPFR_RNDN);
            mpfr_set(Mb[(i + mu) * N + j], H_im[i * m + j], MPFR_RNDN);
            mpfr_set(Mb[(i + mu) * N + (j + mu)], H_re[i * m + j], MPFR_RNDN);
        }
    }

    mpfr_t* Qm = mpfr_array_alloc(N * N, bits);
    for (size_t i = 0; i < N; i++)
        for (size_t j = 0; j < N; j++)
            mpfr_set_si(Qm[i * N + j], (i == j) ? 1 : 0, MPFR_RNDN);
    mpfr_t* ub = mpfr_array_alloc(N, bits);
    mpfr_t* tmpM = mpfr_array_alloc(14, bits);
    if (N >= 3) direct_hessenberg_real_M(Mb, N, bits, ub, Qm, tmpM);

    mpfr_t* M_eval_re = mpfr_array_alloc(N, bits);
    mpfr_t* M_eval_im = mpfr_array_alloc(N, bits);
    for (size_t i = 0; i < N; i++) {
        mpfr_set_zero(M_eval_re[i], 1);
        mpfr_set_zero(M_eval_im[i], 1);
    }
    int qr_status = direct_qr_real_general_M(Mb, N, bits,
                                              M_eval_re, M_eval_im, Qm, tmpM);
    if (qr_status != 0) {
        mpfr_clear(sum_r); mpfr_clear(sum_i);
        mpfr_clear(dot_r); mpfr_clear(dot_i);
        mpfr_clear(pr); mpfr_clear(pi);
        mpfr_clear(prod); mpfr_clear(q);
        mpfr_clear(nrm2); mpfr_clear(nrm);
        mpfr_clear(normA); mpfr_clear(tol);
        mpfr_array_free(V_re, n * (m + 1));
        mpfr_array_free(V_im, n * (m + 1));
        mpfr_array_free(H_re, m * m);
        mpfr_array_free(H_im, m * m);
        mpfr_array_free(Mb, N * N);
        mpfr_array_free(Qm, N * N);
        mpfr_array_free(ub, N);
        mpfr_array_free(tmpM, 14);
        mpfr_array_free(M_eval_re, N);
        mpfr_array_free(M_eval_im, N);
        return NULL;
    }

    /* Eigenvectors of M (in original basis). */
    mpfr_t* M_evec_re = mpfr_array_alloc(N * N, bits);
    mpfr_t* M_evec_im = mpfr_array_alloc(N * N, bits);
    size_t* id_perm = (size_t*)malloc(sizeof(size_t) * N);
    for (size_t i = 0; i < N; i++) id_perm[i] = i;
    schur_compute_eigvecs_M(Mb, Qm, N, bits, M_eval_re, M_eval_im,
                              id_perm, M_evec_re, M_evec_im);
    free(id_perm);
    mpfr_array_free(Mb, N * N);
    mpfr_array_free(Qm, N * N);
    mpfr_array_free(ub, N);
    mpfr_array_free(tmpM, 14);

    /* Recover H_m eigenpairs via x = (a - d) + i(b + c) + grouped GS. */
    mpfr_t spec_norm, mag, dr_diff, di_diff;
    mpfr_t group_tol, threshold_sq;
    mpfr_init2(spec_norm, bits); mpfr_init2(mag, bits);
    mpfr_init2(dr_diff, bits); mpfr_init2(di_diff, bits);
    mpfr_init2(group_tol, bits); mpfr_init2(threshold_sq, bits);

    mpfr_set_zero(spec_norm, 1);
    for (size_t i = 0; i < N; i++) {
        mpfr_hypot(mag, M_eval_re[i], M_eval_im[i], MPFR_RNDN);
        if (mpfr_cmp(mag, spec_norm) > 0) mpfr_set(spec_norm, mag, MPFR_RNDN);
    }
    if (mpfr_zero_p(spec_norm)) mpfr_set_ui(spec_norm, 1, MPFR_RNDN);
    /* group_tol ~ N * ||spec|| / 2^(bits/2 - 6) */
    mpfr_set(group_tol, spec_norm, MPFR_RNDN);
    mpfr_mul_ui(group_tol, group_tol, (unsigned long)N, MPFR_RNDN);
    mpfr_div_2si(group_tol, group_tol, (long)(bits / 2) - 6, MPFR_RNDN);
    /* threshold ~ sqrt(mu) / 2^(bits/2 - 6); threshold_sq is that squared. */
    {
        mpfr_t t;
        mpfr_init2(t, bits);
        mpfr_set_ui(t, (unsigned long)mu, MPFR_RNDN);
        mpfr_sqrt(t, t, MPFR_RNDN);
        mpfr_div_2si(t, t, (long)(bits / 2) - 6, MPFR_RNDN);
        mpfr_mul(threshold_sq, t, t, MPFR_RNDN);
        mpfr_clear(t);
    }

    int* used = (int*)calloc(N, sizeof(int));
    mpfr_t* HA_eval_re = mpfr_array_alloc(mu, bits);
    mpfr_t* HA_eval_im = mpfr_array_alloc(mu, bits);
    mpfr_t* HA_evec_re = mpfr_array_alloc(mu * mu, bits);
    mpfr_t* HA_evec_im = mpfr_array_alloc(mu * mu, bits);
    for (size_t i = 0; i < mu; i++) {
        mpfr_set_zero(HA_eval_re[i], 1);
        mpfr_set_zero(HA_eval_im[i], 1);
    }
    for (size_t i = 0; i < mu * mu; i++) {
        mpfr_set_zero(HA_evec_re[i], 1);
        mpfr_set_zero(HA_evec_im[i], 1);
    }
    mpfr_t* cand_re = mpfr_array_alloc(mu, bits);
    mpfr_t* cand_im = mpfr_array_alloc(mu, bits);
    mpfr_t inv;
    mpfr_init2(inv, bits);

    size_t produced = 0;
    for (size_t i = 0; i < N && produced < mu; i++) {
        if (used[i]) continue;
        used[i] = 1;
        size_t group_start = produced;
        for (size_t jj = i; jj < N && produced < mu; jj++) {
            if (jj != i) {
                if (used[jj]) continue;
                mpfr_sub(dr_diff, M_eval_re[jj], M_eval_re[i], MPFR_RNDN);
                mpfr_sub(di_diff, M_eval_im[jj], M_eval_im[i], MPFR_RNDN);
                mpfr_hypot(mag, dr_diff, di_diff, MPFR_RNDN);
                if (mpfr_cmp(mag, group_tol) > 0) continue;
                used[jj] = 1;
            }
            /* Candidate x = (a - d) + i (b + c). */
            for (size_t l = 0; l < mu; l++) {
                mpfr_sub(cand_re[l],
                         M_evec_re[jj * N + l],
                         M_evec_im[jj * N + (l + mu)], MPFR_RNDN);
                mpfr_add(cand_im[l],
                         M_evec_re[jj * N + (l + mu)],
                         M_evec_im[jj * N + l],        MPFR_RNDN);
            }
            /* Complex CGS (twice) against the group's emitted vectors. */
            for (int pass = 0; pass < 2; pass++) {
                for (size_t f = group_start; f < produced; f++) {
                    mpfr_set_zero(pr, 1); mpfr_set_zero(pi, 1);
                    for (size_t l = 0; l < mu; l++) {
                        /* conj(V_f) . cand = (vr - i vi)(cand_re + i cand_im) */
                        mpfr_mul(prod, HA_evec_re[f * mu + l], cand_re[l], MPFR_RNDN);
                        mpfr_add(pr, pr, prod, MPFR_RNDN);
                        mpfr_mul(prod, HA_evec_im[f * mu + l], cand_im[l], MPFR_RNDN);
                        mpfr_add(pr, pr, prod, MPFR_RNDN);
                        mpfr_mul(prod, HA_evec_re[f * mu + l], cand_im[l], MPFR_RNDN);
                        mpfr_add(pi, pi, prod, MPFR_RNDN);
                        mpfr_mul(prod, HA_evec_im[f * mu + l], cand_re[l], MPFR_RNDN);
                        mpfr_sub(pi, pi, prod, MPFR_RNDN);
                    }
                    for (size_t l = 0; l < mu; l++) {
                        /* (pr + i pi)(vr + i vi) */
                        mpfr_mul(dr_diff, pr, HA_evec_re[f * mu + l], MPFR_RNDN);
                        mpfr_mul(prod,    pi, HA_evec_im[f * mu + l], MPFR_RNDN);
                        mpfr_sub(dr_diff, dr_diff, prod, MPFR_RNDN);
                        mpfr_mul(di_diff, pr, HA_evec_im[f * mu + l], MPFR_RNDN);
                        mpfr_mul(prod,    pi, HA_evec_re[f * mu + l], MPFR_RNDN);
                        mpfr_add(di_diff, di_diff, prod, MPFR_RNDN);
                        mpfr_sub(cand_re[l], cand_re[l], dr_diff, MPFR_RNDN);
                        mpfr_sub(cand_im[l], cand_im[l], di_diff, MPFR_RNDN);
                    }
                }
            }
            /* Check norm against threshold. */
            mpfr_set_zero(nrm2, 1);
            for (size_t l = 0; l < mu; l++) {
                mpfr_mul(prod, cand_re[l], cand_re[l], MPFR_RNDN);
                mpfr_add(nrm2, nrm2, prod, MPFR_RNDN);
                mpfr_mul(prod, cand_im[l], cand_im[l], MPFR_RNDN);
                mpfr_add(nrm2, nrm2, prod, MPFR_RNDN);
            }
            if (mpfr_cmp(nrm2, threshold_sq) < 0) continue;
            mpfr_sqrt(nrm, nrm2, MPFR_RNDN);
            mpfr_ui_div(inv, 1, nrm, MPFR_RNDN);
            for (size_t l = 0; l < mu; l++) {
                mpfr_mul(HA_evec_re[produced * mu + l], cand_re[l], inv, MPFR_RNDN);
                mpfr_mul(HA_evec_im[produced * mu + l], cand_im[l], inv, MPFR_RNDN);
            }
            mpfr_set(HA_eval_re[produced], M_eval_re[i], MPFR_RNDN);
            mpfr_set(HA_eval_im[produced], M_eval_im[i], MPFR_RNDN);
            produced++;
        }
    }
    free(used);
    mpfr_array_free(M_eval_re, N); mpfr_array_free(M_eval_im, N);
    mpfr_array_free(M_evec_re, N * N); mpfr_array_free(M_evec_im, N * N);
    mpfr_array_free(cand_re, mu); mpfr_array_free(cand_im, mu);
    mpfr_clear(spec_norm); mpfr_clear(mag);
    mpfr_clear(dr_diff); mpfr_clear(di_diff);
    mpfr_clear(group_tol); mpfr_clear(threshold_sq);
    mpfr_clear(inv);

    if (produced != mu) {
        mpfr_clear(sum_r); mpfr_clear(sum_i);
        mpfr_clear(dot_r); mpfr_clear(dot_i);
        mpfr_clear(pr); mpfr_clear(pi);
        mpfr_clear(prod); mpfr_clear(q);
        mpfr_clear(nrm2); mpfr_clear(nrm);
        mpfr_clear(normA); mpfr_clear(tol);
        mpfr_array_free(V_re, n * (m + 1));
        mpfr_array_free(V_im, n * (m + 1));
        mpfr_array_free(H_re, m * m);
        mpfr_array_free(H_im, m * m);
        mpfr_array_free(HA_eval_re, mu);
        mpfr_array_free(HA_eval_im, mu);
        mpfr_array_free(HA_evec_re, mu * mu);
        mpfr_array_free(HA_evec_im, mu * mu);
        return NULL;
    }

    /* Sort H_m eigenvalues by descending |lambda|. */
    size_t* perm = (size_t*)malloc(sizeof(size_t) * mu);
    direct_sort_perm_desc_abs_complex_M(HA_eval_re, HA_eval_im, mu, perm);

    Expr* out;
    if (want_Q) {
        /* Lift to A-eigenvectors: VA[k][i] = sum_j V_m[j][i] * y_k[j]
         * where y_k is the perm[k]-th H_m eigenvector (complex). */
        mpfr_t* VA_re = mpfr_array_alloc(mu * n, bits);
        mpfr_t* VA_im = mpfr_array_alloc(mu * n, bits);
        for (size_t k = 0; k < mu; k++) {
            size_t src = perm[k];
            for (size_t i = 0; i < n; i++) {
                mpfr_set_zero(sum_r, 1);
                mpfr_set_zero(sum_i, 1);
                for (size_t j = 0; j < mu; j++) {
                    /* (vr + i vi)(yr + i yi) */
                    mpfr_mul(prod, V_re[j * n + i], HA_evec_re[src * mu + j], MPFR_RNDN);
                    mpfr_add(sum_r, sum_r, prod, MPFR_RNDN);
                    mpfr_mul(prod, V_im[j * n + i], HA_evec_im[src * mu + j], MPFR_RNDN);
                    mpfr_sub(sum_r, sum_r, prod, MPFR_RNDN);
                    mpfr_mul(prod, V_re[j * n + i], HA_evec_im[src * mu + j], MPFR_RNDN);
                    mpfr_add(sum_i, sum_i, prod, MPFR_RNDN);
                    mpfr_mul(prod, V_im[j * n + i], HA_evec_re[src * mu + j], MPFR_RNDN);
                    mpfr_add(sum_i, sum_i, prod, MPFR_RNDN);
                }
                mpfr_set(VA_re[k * n + i], sum_r, MPFR_RNDN);
                mpfr_set(VA_im[k * n + i], sum_i, MPFR_RNDN);
            }
        }
        /* Re-normalise. */
        for (size_t k = 0; k < mu; k++) {
            mpfr_set_zero(nrm2, 1);
            for (size_t i = 0; i < n; i++) {
                mpfr_mul(prod, VA_re[k * n + i], VA_re[k * n + i], MPFR_RNDN);
                mpfr_add(nrm2, nrm2, prod, MPFR_RNDN);
                mpfr_mul(prod, VA_im[k * n + i], VA_im[k * n + i], MPFR_RNDN);
                mpfr_add(nrm2, nrm2, prod, MPFR_RNDN);
            }
            if (mpfr_sgn(nrm2) > 0) {
                mpfr_sqrt(nrm, nrm2, MPFR_RNDN);
                for (size_t i = 0; i < n; i++) {
                    mpfr_div(VA_re[k * n + i], VA_re[k * n + i], nrm, MPFR_RNDN);
                    mpfr_div(VA_im[k * n + i], VA_im[k * n + i], nrm, MPFR_RNDN);
                }
            }
        }
        out = mateigen_build_complex_eigvec_list_rect_M(VA_re, VA_im, mu, n);
        mpfr_array_free(VA_re, mu * n);
        mpfr_array_free(VA_im, mu * n);
    } else {
        out = direct_build_complex_eigenvalue_list_M(HA_eval_re, HA_eval_im,
                                                       mu, perm);
    }
    free(perm);

    mpfr_clear(sum_r); mpfr_clear(sum_i);
    mpfr_clear(dot_r); mpfr_clear(dot_i);
    mpfr_clear(pr); mpfr_clear(pi);
    mpfr_clear(prod); mpfr_clear(q);
    mpfr_clear(nrm2); mpfr_clear(nrm);
    mpfr_clear(normA); mpfr_clear(tol);
    mpfr_array_free(V_re, n * (m + 1));
    mpfr_array_free(V_im, n * (m + 1));
    mpfr_array_free(H_re, m * m);
    mpfr_array_free(H_im, m * m);
    mpfr_array_free(HA_eval_re, mu);
    mpfr_array_free(HA_eval_im, mu);
    mpfr_array_free(HA_evec_re, mu * mu);
    mpfr_array_free(HA_evec_im, mu * mu);

    return direct_apply_k_spec_list(out, k_spec);
}

/* MPFR dispatcher for the Arnoldi method. */
static Expr* arnoldi_dispatch_mpfr(Expr* m, Expr* a, int64_t n,
                                    mpfr_prec_t bits,
                                    MateigenWant want, Expr* k_spec,
                                    const ArnoldiOpts* opts) {
    if (a != NULL) return NULL;
    if (n <= 0)    return NULL;

    MatM A;
    if (!matM_load(m, (size_t)n, bits, &A)) return NULL;

    Expr* out;
    if (A.is_complex) out = arnoldi_complex_general_mpfr(&A, want, k_spec, opts);
    else              out = arnoldi_real_general_mpfr(&A, want, k_spec, opts);

    matM_free(&A);
    return out;
}
#endif /* USE_MPFR */

/* Top-level "Arnoldi" dispatcher.  Parses sub-options from method_value
 * (the right-hand side of the Method -> ... rule) and routes through
 * the MPFR kernel when the input carries MPFR precision; otherwise
 * uses the machine kernel.  Returns NULL when the input shape isn't
 * supported so the caller can fall back to Direct or the symbolic
 * path. */
Expr* arnoldi_dispatch(Expr* m, Expr* a, int64_t n,
                                MateigenWant want, Expr* k_spec,
                                Expr* method_value) {
    if (a != NULL) return NULL;          /* generalised: symbolic only */
    if (n <= 0)    return NULL;

    ArnoldiOpts opts;
    arnoldi_parse_subopts(method_value, &opts);

    CommonInexactInfo info = common_scan_inexact(m);
    if (info.has_inexact && info.min_bits > 53) {
        Expr* out = arnoldi_dispatch_mpfr(m, a, n,
                                           (mpfr_prec_t)info.min_bits,
                                           want, k_spec, &opts);
        if (out) return out;
        /* MPFR Arnoldi failed -- fall through to the machine kernel,
         * which will coerce MPFR cells to doubles via matD_load. */
    }
    return arnoldi_dispatch_machine(m, a, n, want, k_spec, &opts);
}
