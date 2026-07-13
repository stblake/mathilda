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

/* M_PI is a POSIX/GNU extension; glibc hides it under -std=c99.  The
 * FEAST kernel maps Gauss-Legendre nodes on [-1, 1] to the upper-half
 * elliptic contour via t_e = (pi/2)(1 + xi_e), so we need a guaranteed
 * pi constant.  Matches the pattern in src/numeric.c / src/trig.c. */
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


/* ============================================================ *
 *  Phase 5: numerical "FEAST" method (skeleton).                  *
 *                                                                 *
 *  FEAST (Polizzi 2009) is a Hermitian-only spectral-projector    *
 *  method that returns every eigenpair whose eigenvalue lies in   *
 *  a user-supplied real interval [a, b].  The full algorithm is   *
 *
 *    1. Pick a subspace size m0 (>= count of eigenvalues in       *
 *       [a, b]) and initialise Y in C^{n x m0}.                   *
 *    2. Approximate the spectral projector P_[a,b](A) Y via       *
 *       Gauss-Legendre quadrature on the upper half of the        *
 *       elliptic contour through (a, 0) and (b, 0); Schwarz       *
 *       symmetry halves the number of complex linear solves.     *
 *    3. Form A_q = Q^* A Q, B_q = Q^* Q and solve the small       *
 *       generalised Hermitian-definite problem A_q y = lambda     *
 *       B_q y by Cholesky reduction + Direct.                     *
 *    4. Filter Ritz pairs against [a, b] and iterate until the    *
 *       residual converges or MaxIterations is reached.           *
 *                                                                 *
 *  Phase 1 (this commit) ships the option grammar, dispatcher     *
 *  skeleton, and the `feast_automatic_prefers` heuristic.  Both   *
 *  feast_dispatch_machine and feast_dispatch_mpfr return NULL,    *
 *  so Method -> "FEAST" silently falls through to Direct -- the   *
 *  observable behaviour is unchanged from before Phase 1.  The    *
 *  per-precision kernels arrive in Phases 2 (real symmetric, mac- *
 *  hine), 3 (complex Hermitian, machine), and 4 (MPFR analogues). *
 * ============================================================ */

/* FEAST sub-options.  All fields default to a sentinel meaning
 * "kernel chooses": interval_given == false signals a missing
 * Interval, basis_size == 0 means "auto-size from k_spec and n",
 * and tolerance == 0.0 means "use the precision-aware default".
 *
 * The interval bounds are stored as doubles since FEAST is a
 * floating-point algorithm; MPFR inputs widen these to the working
 * precision inside the kernel rather than carrying mpfr_t through
 * the option struct. */
typedef struct {
    double  interval_low;
    double  interval_high;
    bool    interval_given;
    int64_t contour_points;   /* N_e; supported by kernels: 2, 4, 8, 16 */
    int64_t subspace_size;    /* m0; 0 = auto */
    int64_t max_iterations;
    double  tolerance;        /* 0 = precision-aware default */
    bool    tolerance_given;
} FeastOpts;

static void feast_set_defaults(FeastOpts* o) {
    o->interval_low     = 0.0;
    o->interval_high    = 0.0;
    o->interval_given   = false;
    o->contour_points   = 8;
    o->subspace_size    = 0;
    o->max_iterations   = 20;
    o->tolerance        = 0.0;
    o->tolerance_given  = false;
}

/* Match an interned-symbol or string sub-option key, mirroring
 * arnoldi_subopt_key_eq.  Returns true when `key` denotes `name`. */
static bool feast_subopt_key_eq(Expr* key, const char* name,
                                  const char* sym_intern) {
    if (!key) return false;
    if (key->type == EXPR_STRING && strcmp(key->data.string, name) == 0)
        return true;
    if (key->type == EXPR_SYMBOL && sym_intern
        && key->data.symbol.name == sym_intern)
        return true;
    return false;
}

/* Coerce a numeric Expr leaf to double.  Mirrors arnoldi_coerce_double
 * (kept independent so the FEAST section stays self-contained and to
 * leave room for FEAST-specific value coercions later). */
static double feast_coerce_double(Expr* e) {
    if (!e) return NAN;
    if (e->type == EXPR_INTEGER) return (double)e->data.integer;
    if (e->type == EXPR_REAL)    return e->data.real;
    if (e->type == EXPR_MPFR)    return mpfr_get_d(e->data.mpfr, MPFR_RNDN);
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Rational
        && e->data.function.arg_count == 2) {
        Expr* p = e->data.function.args[0];
        Expr* q = e->data.function.args[1];
        if (p->type == EXPR_INTEGER && q->type == EXPR_INTEGER
            && q->data.integer != 0)
            return (double)p->data.integer / (double)q->data.integer;
    }
    return NAN;
}

/* Parse the right-hand side of a Method rule into FeastOpts.  The RHS
 * may be the bare string/symbol "FEAST" (all defaults; interval missing)
 * or a List of the form {"FEAST", "Key" -> value, ...}.  Unknown keys
 * are silently ignored, matching the Arnoldi / Banded conventions.
 *
 * Recognised keys:
 *   "Interval"      -> {a, b}        required at the call site
 *   "ContourPoints" -> integer       N_e
 *   "SubspaceSize"  -> integer       m_0
 *   "MaxIterations" -> integer       outer iteration cap
 *   "Tolerance"     -> numeric       residual stopping criterion
 */
static void feast_parse_subopts(Expr* method_value, FeastOpts* opts) {
    feast_set_defaults(opts);
    if (!method_value) return;
    if (method_value->type != EXPR_FUNCTION) return;
    Expr* head = method_value->data.function.head;
    if (head->type != EXPR_SYMBOL || head->data.symbol.name != SYM_List) return;

    for (size_t i = 1; i < method_value->data.function.arg_count; i++) {
        Expr* rule = method_value->data.function.args[i];
        if (rule->type != EXPR_FUNCTION) continue;
        if (rule->data.function.arg_count != 2) continue;
        Expr* rh = rule->data.function.head;
        if (rh->type != EXPR_SYMBOL) continue;
        if (rh->data.symbol.name != SYM_Rule
            && rh->data.symbol.name != SYM_RuleDelayed) continue;

        Expr* key = rule->data.function.args[0];
        Expr* val = rule->data.function.args[1];

        if (feast_subopt_key_eq(key, "Interval", SYM_Interval)) {
            /* Expect a List of exactly two numeric entries.  Anything
             * else leaves interval_given == false so the dispatcher
             * falls through to Direct. */
            if (val->type != EXPR_FUNCTION) continue;
            Expr* vh = val->data.function.head;
            if (vh->type != EXPR_SYMBOL || vh->data.symbol.name != SYM_List) continue;
            if (val->data.function.arg_count != 2) continue;
            double a = feast_coerce_double(val->data.function.args[0]);
            double b = feast_coerce_double(val->data.function.args[1]);
            if (isnan(a) || isnan(b)) continue;
            if (a > b) { double t = a; a = b; b = t; }
            opts->interval_low   = a;
            opts->interval_high  = b;
            opts->interval_given = true;
        } else if (feast_subopt_key_eq(key, "ContourPoints",
                                          SYM_ContourPoints)
                   && val->type == EXPR_INTEGER) {
            opts->contour_points = val->data.integer;
        } else if (feast_subopt_key_eq(key, "SubspaceSize",
                                          SYM_SubspaceSize)
                   && val->type == EXPR_INTEGER) {
            opts->subspace_size = val->data.integer;
        } else if (feast_subopt_key_eq(key, "MaxIterations",
                                          SYM_MaxIterations)
                   && val->type == EXPR_INTEGER) {
            opts->max_iterations = val->data.integer;
        } else if (feast_subopt_key_eq(key, "Tolerance", SYM_Tolerance)) {
            double d = feast_coerce_double(val);
            if (!isnan(d)) { opts->tolerance = d; opts->tolerance_given = true; }
        }
    }
}

/* Heuristic for the Automatic dispatcher.  FEAST is never auto-selected:
 * it requires a user-supplied interval and silently picking a contour
 * would be surprising.  Direct / Arnoldi / Banded already cover the
 * Automatic policy. */
bool feast_automatic_prefers(Expr* m, int64_t n, Expr* method_value) {
    (void)m;
    (void)n;
    (void)method_value;
    return false;
}

/* Emit a once-per-process stderr warning when a FEAST run had to fall
 * back to Direct (sub-space too small, LU singular at a quadrature
 * node, Cholesky failure on B_q, non-convergence within
 * MaxIterations).  Tagging the reason helps the user pick a smarter
 * SubspaceSize / Interval / ContourPoints when they retry.  Kept
 * silent on the first call to avoid noisy output when the kernel
 * cleanly converges. */
static void feast_warn_fallback(const char* reason) {
    static bool warned = false;
    if (warned) return;
    warned = true;
    fprintf(stderr,
        "Eigenvalues::feast: FEAST machine kernel could not converge "
        "(%s); falling back to the Direct method.\n",
        reason ? reason : "unknown");
}

/* --------- Complex LU factor / solve on paired re/im arrays --------- *
 *                                                                       *
 *  Used by FEAST to solve the shifted system (z_e I - A) X = Y at each  *
 *  Gauss-Legendre node z_e on the contour.  A is row-major n*n with     *
 *  separate re[] / im[] storage, matching the convention used by MatD   *
 *  elsewhere in this file.  Partial pivoting on |z| keeps the           *
 *  factorisation stable when (z_e I - A) gets close to a real           *
 *  eigenvalue (rare, since FEAST chooses contours away from the         *
 *  spectrum).                                                           *
 *                                                                       *
 *  Returns 0 on success, -1 on numerical singularity at a pivot.        *
 * --------------------------------------------------------------------- */
static int feast_complex_lu_factor(double* A_re, double* A_im,
                                     size_t n, int* piv) {
    for (size_t k = 0; k < n; k++) {
        size_t pivot = k;
        double max_abs = hypot(A_re[k * n + k], A_im[k * n + k]);
        for (size_t i = k + 1; i < n; i++) {
            double a = hypot(A_re[i * n + k], A_im[i * n + k]);
            if (a > max_abs) { max_abs = a; pivot = i; }
        }
        if (max_abs == 0.0) return -1;
        piv[k] = (int)pivot;
        if (pivot != k) {
            for (size_t j = 0; j < n; j++) {
                double t = A_re[k * n + j];
                A_re[k * n + j] = A_re[pivot * n + j];
                A_re[pivot * n + j] = t;
                t = A_im[k * n + j];
                A_im[k * n + j] = A_im[pivot * n + j];
                A_im[pivot * n + j] = t;
            }
        }
        double pkr = A_re[k * n + k];
        double pki = A_im[k * n + k];
        double pk_denom = pkr * pkr + pki * pki;
        for (size_t i = k + 1; i < n; i++) {
            double ar = A_re[i * n + k];
            double ai = A_im[i * n + k];
            /* factor = a_ik / a_kk = (ar+i ai)*(pkr-i pki)/|pk|^2. */
            double fr = (ar * pkr + ai * pki) / pk_denom;
            double fi = (ai * pkr - ar * pki) / pk_denom;
            A_re[i * n + k] = fr;
            A_im[i * n + k] = fi;
            for (size_t j = k + 1; j < n; j++) {
                double xr = A_re[k * n + j];
                double xi = A_im[k * n + j];
                A_re[i * n + j] -= fr * xr - fi * xi;
                A_im[i * n + j] -= fr * xi + fi * xr;
            }
        }
    }
    return 0;
}

/* Solve LUx = Pb in place over a single complex RHS b (length n),
 * using the factorisation produced by feast_complex_lu_factor.  Unit-
 * diagonal L occupies the strict lower triangle of LU; U occupies the
 * upper triangle including the diagonal. */
static void feast_complex_lu_solve(const double* LU_re,
                                     const double* LU_im,
                                     const int* piv, size_t n,
                                     double* b_re, double* b_im) {
    for (size_t k = 0; k < n; k++) {
        int p = piv[k];
        if ((size_t)p != k) {
            double t = b_re[k]; b_re[k] = b_re[p]; b_re[p] = t;
                   t = b_im[k]; b_im[k] = b_im[p]; b_im[p] = t;
        }
    }
    /* Forward solve L y = Pb (L unit-lower). */
    for (size_t i = 1; i < n; i++) {
        double sr = b_re[i], si = b_im[i];
        for (size_t j = 0; j < i; j++) {
            double lr = LU_re[i * n + j];
            double li = LU_im[i * n + j];
            double yr = b_re[j];
            double yi = b_im[j];
            sr -= lr * yr - li * yi;
            si -= lr * yi + li * yr;
        }
        b_re[i] = sr;
        b_im[i] = si;
    }
    /* Backward solve U x = y. */
    for (size_t ii = n; ii-- > 0; ) {
        double sr = b_re[ii], si = b_im[ii];
        for (size_t j = ii + 1; j < n; j++) {
            double ur = LU_re[ii * n + j];
            double ui = LU_im[ii * n + j];
            double xr = b_re[j];
            double xi = b_im[j];
            sr -= ur * xr - ui * xi;
            si -= ur * xi + ui * xr;
        }
        double dr = LU_re[ii * n + ii];
        double di = LU_im[ii * n + ii];
        double denom = dr * dr + di * di;
        b_re[ii] = (sr * dr + si * di) / denom;
        b_im[ii] = (si * dr - sr * di) / denom;
    }
}

/* --------- Gauss-Legendre nodes / weights on [-1, 1] -----------------*
 *                                                                       *
 *  Standard tabulated values (Abramowitz & Stegun 25.4.30, or any       *
 *  numerical-analysis textbook).  Nodes are ordered ascending; weights  *
 *  mirror the symmetric pairing.  FEAST maps these to the upper-half    *
 *  elliptic contour through (a, 0) and (b, 0) via t_e = (pi/2)(1+xi_e). *
 *                                                                       *
 *  Convention: Ne == number of upper-half quadrature nodes == number    *
 *  of complex linear solves per FEAST iteration.  This matches          *
 *  Polizzi's FEAST 4.0 reference implementation where Ne=8 gives 8      *
 *  solves and convergence to ~1e-12 on cleanly-separated spectra.       *
 * --------------------------------------------------------------------- */
static const double feast_gl_xi_2[2] = {
    -0.5773502691896257645091488,
     0.5773502691896257645091488
};
static const double feast_gl_w_2[2] = { 1.0, 1.0 };
static const double feast_gl_xi_4[4] = {
    -0.8611363115940525752239465,
    -0.3399810435848562648026658,
     0.3399810435848562648026658,
     0.8611363115940525752239465
};
static const double feast_gl_w_4[4] = {
    0.3478548451374538573730639,
    0.6521451548625461426269361,
    0.6521451548625461426269361,
    0.3478548451374538573730639
};
static const double feast_gl_xi_8[8] = {
    -0.9602898564975362316835609,
    -0.7966664774136267395915539,
    -0.5255324099163289858177390,
    -0.1834346424956498049394761,
     0.1834346424956498049394761,
     0.5255324099163289858177390,
     0.7966664774136267395915539,
     0.9602898564975362316835609
};
static const double feast_gl_w_8[8] = {
    0.1012285362903762591525314,
    0.2223810344533744705443560,
    0.3137066458778872873379622,
    0.3626837833783619829651504,
    0.3626837833783619829651504,
    0.3137066458778872873379622,
    0.2223810344533744705443560,
    0.1012285362903762591525314
};
static const double feast_gl_xi_16[16] = {
    -0.9894009349916499325961542,
    -0.9445750230732325760779884,
    -0.8656312023878317438804679,
    -0.7554044083550030338951012,
    -0.6178762444026437484466718,
    -0.4580167776572273863424194,
    -0.2816035507792589132304605,
    -0.0950125098376374401853193,
     0.0950125098376374401853193,
     0.2816035507792589132304605,
     0.4580167776572273863424194,
     0.6178762444026437484466718,
     0.7554044083550030338951012,
     0.8656312023878317438804679,
     0.9445750230732325760779884,
     0.9894009349916499325961542
};
static const double feast_gl_w_16[16] = {
    0.0271524594117540948517806,
    0.0622535239386478928628438,
    0.0951585116824927848099251,
    0.1246289712555338720524763,
    0.1495959888165767320815017,
    0.1691565193950025381893121,
    0.1826034150449235888667637,
    0.1894506104550684962853967,
    0.1894506104550684962853967,
    0.1826034150449235888667637,
    0.1691565193950025381893121,
    0.1495959888165767320815017,
    0.1246289712555338720524763,
    0.0951585116824927848099251,
    0.0622535239386478928628438,
    0.0271524594117540948517806
};

/* Return the GL order matching the user's ContourPoints request,
 * populating *xi and *w with the matching node/weight tables.  Any
 * value outside {2, 4, 8, 16} silently falls back to 8 (matching
 * FEAST's default). */
static int64_t feast_gl_lookup(int64_t ne, const double** xi,
                                 const double** w) {
    switch (ne) {
        case 2:  *xi = feast_gl_xi_2;  *w = feast_gl_w_2;  return 2;
        case 4:  *xi = feast_gl_xi_4;  *w = feast_gl_w_4;  return 4;
        case 16: *xi = feast_gl_xi_16; *w = feast_gl_w_16; return 16;
        case 8:
        default: *xi = feast_gl_xi_8;  *w = feast_gl_w_8;  return 8;
    }
}

/* Deterministic LCG-seeded subspace for FEAST initialisation.  Fills Y
 * (n*m0 row-major) with values in (-1, 1).  Determinism keeps CI
 * reproducible at the cost of a fixed-seed dependency -- Mathematica's
 * "ResetEigenvalues" sub-option, which re-randomises on each call, is
 * a future addition we can layer on top.  Constants from Knuth (TAOCP
 * vol 2 sec 3.3.4 table 1, line 25). */
static void feast_seed_subspace(double* Y, size_t n, size_t m0,
                                  uint64_t seed) {
    uint64_t state = seed ? seed : 0xC0FFEE05u;
    for (size_t i = 0; i < n * m0; i++) {
        state = state * 6364136223846793005ull + 1442695040888963407ull;
        uint32_t v = (uint32_t)(state >> 32);
        Y[i] = ((double)v * (2.0 / 4294967295.0)) - 1.0;
    }
}

/* Cholesky factorisation of an n x n SPD matrix B (row-major) into a
 * lower-triangular L (row-major, strict upper triangle zeroed).
 * Returns 0 on success, -1 on a non-positive pivot (rank-deficient or
 * indefinite B).  Allocated in fresh storage so the caller can keep B
 * intact for the back-transform of eigenvectors. */
static int feast_cholesky(const double* B, size_t n, double* L) {
    memset(L, 0, sizeof(double) * n * n);
    for (size_t j = 0; j < n; j++) {
        double s = B[j * n + j];
        for (size_t k = 0; k < j; k++) {
            double v = L[j * n + k];
            s -= v * v;
        }
        if (s <= 0.0) return -1;
        double djj = sqrt(s);
        L[j * n + j] = djj;
        for (size_t i = j + 1; i < n; i++) {
            double t = B[i * n + j];
            for (size_t k = 0; k < j; k++) {
                t -= L[i * n + k] * L[j * n + k];
            }
            L[i * n + j] = t / djj;
        }
    }
    return 0;
}

/* In-place solve L * X = B, L lower-triangular n x n (full storage,
 * strict upper triangle ignored), B = X row-major n x m. */
static void feast_trsm_left_lower(const double* L, double* B,
                                     size_t n, size_t m) {
    for (size_t i = 0; i < n; i++) {
        double inv_lii = 1.0 / L[i * n + i];
        for (size_t col = 0; col < m; col++) {
            double s = B[i * m + col];
            for (size_t k = 0; k < i; k++) {
                s -= L[i * n + k] * B[k * m + col];
            }
            B[i * m + col] = s * inv_lii;
        }
    }
}

/* In-place solve L^T * X = B, L lower-triangular n x n (so L^T is
 * upper-triangular), B = X row-major n x m.  Solve from bottom up. */
static void feast_trsm_left_lower_T(const double* L, double* B,
                                       size_t n, size_t m) {
    for (size_t ii = n; ii-- > 0; ) {
        double inv_lii = 1.0 / L[ii * n + ii];
        for (size_t col = 0; col < m; col++) {
            double s = B[ii * m + col];
            for (size_t k = ii + 1; k < n; k++) {
                s -= L[k * n + ii] * B[k * m + col];
            }
            B[ii * m + col] = s * inv_lii;
        }
    }
}

/* Build the eigenvalue/eigenvector result lists for the FEAST output,
 * where the # of returned eigenvalues (kept) may be smaller than the
 * matrix dimension n.  Mirrors direct_build_real_eigenvalue_list /
 * direct_build_real_eigenvector_list but separates "dim" from "count"
 * since FEAST returns m eigenpairs of n-dimensional vectors with m<=n.
 *
 * `vecs` is row-major dim x count (vecs[i*count + k] is the i-th
 * component of the k-th eigenvector); `perm` selects and orders the
 * columns. */
static Expr* feast_build_real_eigenvalue_list_subset(const double* vals,
                                                      size_t count,
                                                      const size_t* perm) {
    Expr** items = (Expr**)malloc(sizeof(Expr*) * count);
    for (size_t i = 0; i < count; i++) {
        items[i] = expr_new_real(vals[perm[i]]);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, count);
    free(items);
    return out;
}

static Expr* feast_build_real_eigenvector_list_subset(const double* vecs,
                                                       size_t dim,
                                                       size_t count,
                                                       const size_t* perm) {
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * count);
    for (size_t k = 0; k < count; k++) {
        size_t col = perm[k];
        Expr** comps = (Expr**)malloc(sizeof(Expr*) * dim);
        for (size_t i = 0; i < dim; i++) {
            comps[i] = expr_new_real(vecs[i * count + col]);
        }
        rows[k] = expr_new_function(expr_new_symbol(SYM_List), comps, dim);
        free(comps);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), rows, count);
    free(rows);
    return out;
}

/* Real-symmetric FEAST machine-precision kernel.
 *
 * Implements the spectral-projector contour integral
 *
 *     P_[a,b](A) = (1/2 pi i) integral_C (z I - A)^{-1} dz
 *
 * on the elliptic contour through (a, 0) and (b, 0).  Schwarz
 * symmetry (A is real) lets us evaluate the integral on the upper
 * half only and take the real part: the per-iteration projector
 * application is
 *
 *     Q = sum_e Re( w_e * X_e ),   w_e = (1/2) eta_e r exp(i t_e)
 *
 * where (xi_e, eta_e) are Gauss-Legendre nodes/weights on [-1, 1],
 * t_e = (pi/2)(1+xi_e), c = (a+b)/2, r = (b-a)/2, and X_e is the
 * solution to (z_e I - A) X_e = Y.
 *
 * Returns NULL on any failure (LU breakdown, Cholesky failure on B_q,
 * non-convergence within MaxIterations).  The caller then cascades to
 * direct_dispatch -- the user still gets all eigenvalues, just not
 * restricted to [a, b]. */
static Expr* feast_real_sym_machine(const MatD* A, MateigenWant want,
                                      Expr* k_spec,
                                      const FeastOpts* opts) {
    size_t n = A->n;
    if (n == 0) return NULL;
    if (opts->interval_high <= opts->interval_low) return NULL;

    /* Subspace size.  Auto: max(20, n/4), capped at n.  The user can
     * override with "SubspaceSize" -> integer. */
    size_t m0;
    if (opts->subspace_size > 0) {
        m0 = (size_t)opts->subspace_size;
    } else {
        size_t auto_m = n / 4;
        if (auto_m < 20) auto_m = 20;
        m0 = auto_m;
    }
    if (m0 > n) m0 = n;
    if (m0 == 0) return NULL;

    bool   want_Q   = (want & MATEIGEN_WANT_VECTORS) != 0;
    int64_t max_iter = opts->max_iterations > 0 ? opts->max_iterations : 20;
    double  tol      = (opts->tolerance_given && opts->tolerance > 0.0)
                         ? opts->tolerance : 1e-12;
    double  a_lo = opts->interval_low;
    double  b_hi = opts->interval_high;
    double  c    = 0.5 * (a_lo + b_hi);
    double  r    = 0.5 * (b_hi - a_lo);

    const double* gl_xi = NULL;
    const double* gl_w  = NULL;
    int64_t Ne = feast_gl_lookup(opts->contour_points, &gl_xi, &gl_w);

    /* Workspace allocations.  Y / Q / Xhat are n x m0 row-major; the
     * "small" buffers Aq / Bq / L / H / HQ are m0 x m0.  M_re / M_im
     * hold (z_e I - A) and are reused across nodes (LU is destroyed
     * in place but we rebuild M at each node anyway). */
    double* Y      = (double*)malloc(sizeof(double) * n * m0);
    double* Q      = (double*)malloc(sizeof(double) * n * m0);
    double* AQ     = (double*)malloc(sizeof(double) * n * m0);
    double* Xhat   = (double*)malloc(sizeof(double) * n * m0);
    double* Aq     = (double*)malloc(sizeof(double) * m0 * m0);
    double* Bq     = (double*)malloc(sizeof(double) * m0 * m0);
    double* Lchol  = (double*)malloc(sizeof(double) * m0 * m0);
    double* H      = (double*)malloc(sizeof(double) * m0 * m0);
    double* HQ     = (double*)malloc(sizeof(double) * m0 * m0);
    double* W      = (double*)malloc(sizeof(double) * m0 * m0);
    double* diag_h = (double*)malloc(sizeof(double) * m0);
    double* sub_h  = (double*)calloc(m0 ? m0 : 1, sizeof(double));
    double* u_h    = (double*)malloc(sizeof(double) * m0);
    double* p_h    = (double*)malloc(sizeof(double) * m0);
    double* q_h    = (double*)malloc(sizeof(double) * m0);
    double* M_re   = (double*)malloc(sizeof(double) * n * n);
    double* M_im   = (double*)malloc(sizeof(double) * n * n);
    int*    piv    = (int*)malloc(sizeof(int) * n);
    double* x_re   = (double*)malloc(sizeof(double) * n);
    double* x_im   = (double*)malloc(sizeof(double) * n);
    double* mu     = (double*)malloc(sizeof(double) * m0);
    size_t* keep   = (size_t*)malloc(sizeof(size_t) * m0);

    feast_seed_subspace(Y, n, m0, /*seed=*/ 1442695040888963407ull);

    Expr*  out      = NULL;
    bool   converged = false;
    size_t kept     = 0;
    const char* fail_reason = NULL;

    for (int64_t iter = 0; iter < max_iter && !converged; iter++) {
        /* ---- Step 1: Q = sum_e Re(w_e (z_e I - A)^{-1} Y) ---- */
        memset(Q, 0, sizeof(double) * n * m0);
        for (int64_t e = 0; e < Ne; e++) {
            double t  = (M_PI * 0.5) * (1.0 + gl_xi[e]);
            double ct = cos(t);
            double st = sin(t);
            double z_re = c + r * ct;
            double z_im =     r * st;
            /* w_e = (1/2) eta_e r (cos t + i sin t). */
            double w_re = 0.5 * gl_w[e] * r * ct;
            double w_im = 0.5 * gl_w[e] * r * st;

            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                    M_re[i * n + j] = -A->re[i * n + j];
                    M_im[i * n + j] = 0.0;
                }
                M_re[i * n + i] += z_re;
                M_im[i * n + i] += z_im;
            }
            if (feast_complex_lu_factor(M_re, M_im, n, piv) != 0) {
                fail_reason = "LU breakdown at a quadrature node";
                goto cleanup;
            }
            for (size_t col = 0; col < m0; col++) {
                for (size_t i = 0; i < n; i++) {
                    x_re[i] = Y[i * m0 + col];
                    x_im[i] = 0.0;
                }
                feast_complex_lu_solve(M_re, M_im, piv, n, x_re, x_im);
                /* Q[:,col] += Re(w * x) = w_re x_re - w_im x_im. */
                for (size_t i = 0; i < n; i++) {
                    Q[i * m0 + col] += w_re * x_re[i] - w_im * x_im[i];
                }
            }
        }

        /* ---- Step 2: A_q = Q^T A Q,  B_q = Q^T Q ---- */
        /* AQ = A * Q  (n x m0). */
        for (size_t i = 0; i < n; i++) {
            for (size_t col = 0; col < m0; col++) {
                double s = 0.0;
                for (size_t k2 = 0; k2 < n; k2++) {
                    s += A->re[i * n + k2] * Q[k2 * m0 + col];
                }
                AQ[i * m0 + col] = s;
            }
        }
        /* A_q = Q^T AQ. */
        for (size_t i = 0; i < m0; i++) {
            for (size_t j = 0; j < m0; j++) {
                double s = 0.0;
                for (size_t k2 = 0; k2 < n; k2++) {
                    s += Q[k2 * m0 + i] * AQ[k2 * m0 + j];
                }
                Aq[i * m0 + j] = s;
            }
        }
        /* B_q = Q^T Q. */
        for (size_t i = 0; i < m0; i++) {
            for (size_t j = 0; j < m0; j++) {
                double s = 0.0;
                for (size_t k2 = 0; k2 < n; k2++) {
                    s += Q[k2 * m0 + i] * Q[k2 * m0 + j];
                }
                Bq[i * m0 + j] = s;
            }
        }
        /* Symmetrise to suppress O(eps) drift before Cholesky. */
        for (size_t i = 0; i < m0; i++) {
            for (size_t j = i + 1; j < m0; j++) {
                double avg = 0.5 * (Aq[i * m0 + j] + Aq[j * m0 + i]);
                Aq[i * m0 + j] = Aq[j * m0 + i] = avg;
                avg = 0.5 * (Bq[i * m0 + j] + Bq[j * m0 + i]);
                Bq[i * m0 + j] = Bq[j * m0 + i] = avg;
            }
        }

        /* ---- Step 3: Cholesky B_q = L L^T,  H = L^{-1} A_q L^{-T} ---- */
        if (feast_cholesky(Bq, m0, Lchol) != 0) {
            fail_reason = "B_q not positive definite "
                          "(subspace likely undersized or rank-deficient)";
            goto cleanup;
        }
        memcpy(H, Aq, sizeof(double) * m0 * m0);
        feast_trsm_left_lower(Lchol, H, m0, m0);   /* H := L^{-1} A_q */
        /* Now apply L^{-1} from the right: H := H L^{-T}.  Equivalent
         * to (L^{-1} H^T)^T -- transpose, left-solve, transpose back. */
        for (size_t i = 0; i < m0; i++)
            for (size_t j = 0; j < m0; j++)
                HQ[i * m0 + j] = H[j * m0 + i];
        feast_trsm_left_lower(Lchol, HQ, m0, m0);
        for (size_t i = 0; i < m0; i++)
            for (size_t j = 0; j < m0; j++)
                H[i * m0 + j] = HQ[j * m0 + i];
        /* Re-symmetrise. */
        for (size_t i = 0; i < m0; i++) {
            for (size_t j = i + 1; j < m0; j++) {
                double avg = 0.5 * (H[i * m0 + j] + H[j * m0 + i]);
                H[i * m0 + j] = H[j * m0 + i] = avg;
            }
        }

        /* ---- Step 4: eigendecompose H (real symmetric, m0 x m0) ---- *
         * direct_tridiag_real_sym destroys H in place; W receives the   *
         * orthogonal eigenvector matrix of H.                           */
        direct_tridiag_real_sym(H, m0, diag_h, sub_h, W, true,
                                  u_h, p_h, q_h);
        if (direct_symtridiag_qr(diag_h, sub_h, m0, W, true) != 0) {
            fail_reason = "inner symmetric tridiagonal QR did not converge";
            goto cleanup;
        }
        /* Generalised eigenvectors of (A_q, B_q):  V = L^{-T} W. */
        feast_trsm_left_lower_T(Lchol, W, m0, m0);
        for (size_t k = 0; k < m0; k++) mu[k] = diag_h[k];

        /* Ritz vectors Xhat = Q V  (n x m0). */
        for (size_t i = 0; i < n; i++) {
            for (size_t col = 0; col < m0; col++) {
                double s = 0.0;
                for (size_t k2 = 0; k2 < m0; k2++) {
                    s += Q[i * m0 + k2] * W[k2 * m0 + col];
                }
                Xhat[i * m0 + col] = s;
            }
        }

        /* ---- Step 5: filter to [a, b] and compute residuals ---- */
        double pad = tol * fmax(fabs(a_lo), fabs(b_hi));
        if (pad < tol) pad = tol;
        size_t k_inside = 0;
        double max_rel_res = 0.0;
        for (size_t k = 0; k < m0; k++) {
            double muk = mu[k];
            if (muk < a_lo - pad) continue;
            if (muk > b_hi + pad) continue;
            /* residual: ||A xhat - mu xhat|| / max(|mu|, 1). */
            double res2 = 0.0;
            double xn2  = 0.0;
            for (size_t i = 0; i < n; i++) {
                double row = 0.0;
                for (size_t j = 0; j < n; j++) {
                    row += A->re[i * n + j] * Xhat[j * m0 + k];
                }
                double xv = Xhat[i * m0 + k];
                double diff = row - muk * xv;
                res2 += diff * diff;
                xn2  += xv * xv;
            }
            double res = sqrt(res2);
            double xn  = sqrt(xn2 > 0.0 ? xn2 : 1.0);
            double rel = res / (xn * fmax(fabs(muk), 1.0));
            if (rel > max_rel_res) max_rel_res = rel;
            keep[k_inside++] = k;
        }

        if (k_inside > 0 && max_rel_res < tol) {
            converged = true;
            kept = k_inside;
            break;
        }
        /* Not converged: continue.  Refresh Y from the current Ritz
         * basis so the next iteration polishes the same subspace.
         * Renormalise columns to unit length for numerical stability. */
        for (size_t col = 0; col < m0; col++) {
            double nrm2 = 0.0;
            for (size_t i = 0; i < n; i++) {
                double v = Xhat[i * m0 + col];
                nrm2 += v * v;
            }
            double inv = (nrm2 > 0.0) ? 1.0 / sqrt(nrm2) : 0.0;
            if (inv == 0.0) {
                /* Column collapsed; re-seed deterministically from the
                 * column index so different columns stay distinct. */
                for (size_t i = 0; i < n; i++) {
                    uint64_t s = (uint64_t)(col + 1) * 2654435761ull
                                + (uint64_t)i * 40503ull;
                    s ^= s >> 33;
                    s *= 0xff51afd7ed558ccdull;
                    s ^= s >> 33;
                    Y[i * m0 + col] = ((double)(uint32_t)(s >> 32)
                                         * (2.0 / 4294967295.0)) - 1.0;
                }
            } else {
                for (size_t i = 0; i < n; i++) {
                    Y[i * m0 + col] = Xhat[i * m0 + col] * inv;
                }
            }
        }
    }

    if (!converged) {
        if (!fail_reason) fail_reason = "MaxIterations exhausted";
        goto cleanup;
    }

    /* ---- Step 6: build the result ---- *
     * Pack the kept eigenpairs into compact arrays, sort descending by
     * |mu| to match the rest of the eigensolver, then emit. */
    {
        double* vals = (double*)malloc(sizeof(double) * kept);
        double* vecs = want_Q
            ? (double*)malloc(sizeof(double) * n * kept) : NULL;
        for (size_t i = 0; i < kept; i++) {
            size_t src = keep[i];
            vals[i] = mu[src];
            if (want_Q) {
                double nrm2 = 0.0;
                for (size_t r2 = 0; r2 < n; r2++) {
                    double v = Xhat[r2 * m0 + src];
                    nrm2 += v * v;
                }
                double inv = (nrm2 > 0.0) ? 1.0 / sqrt(nrm2) : 1.0;
                for (size_t r2 = 0; r2 < n; r2++) {
                    vecs[r2 * kept + i] = Xhat[r2 * m0 + src] * inv;
                }
            }
        }
        size_t* perm = (size_t*)malloc(sizeof(size_t) * kept);
        direct_sort_perm_desc_abs(vals, kept, perm);
        out = want_Q
            ? feast_build_real_eigenvector_list_subset(vecs, n, kept, perm)
            : feast_build_real_eigenvalue_list_subset(vals, kept, perm);
        free(perm);
        free(vals);
        if (vecs) free(vecs);
        out = direct_apply_k_spec_list(out, k_spec);
    }

cleanup:
    free(Y); free(Q); free(AQ); free(Xhat);
    free(Aq); free(Bq); free(Lchol); free(H); free(HQ); free(W);
    free(diag_h); free(sub_h); free(u_h); free(p_h); free(q_h);
    free(M_re); free(M_im); free(piv); free(x_re); free(x_im);
    free(mu); free(keep);

    if (!converged && fail_reason) feast_warn_fallback(fail_reason);
    return out;
}

/* --------- Complex Cholesky / triangular solves (Hermitian B) -------- *
 *                                                                       *
 *  feast_complex_hermitian_machine reduces the small generalised        *
 *  Hermitian-definite Rayleigh-Ritz problem A_q v = mu B_q v to a       *
 *  standard Hermitian eigenproblem via Cholesky B_q = L L^*.  The       *
 *  primitives below operate on paired re/im double arrays, matching    *
 *  the rest of mateigen's complex matrix convention.                    *
 * --------------------------------------------------------------------- */

/* Complex Hermitian Cholesky: B = L L^*, L lower-triangular with real-
 * positive diagonal.  B is row-major n*n; only the lower-triangle of
 * B is read (upper is implicit by Hermitian symmetry).  L is written
 * into fresh storage; strict-upper L is zeroed.  Returns 0 on success,
 * -1 on a non-positive pivot (rank-deficient or non-Hermitian B). */
static int feast_complex_cholesky(const double* B_re, const double* B_im,
                                    size_t n, double* L_re, double* L_im) {
    memset(L_re, 0, sizeof(double) * n * n);
    memset(L_im, 0, sizeof(double) * n * n);
    for (size_t j = 0; j < n; j++) {
        double s = B_re[j * n + j];   /* Hermitian: B_ii is real. */
        for (size_t k = 0; k < j; k++) {
            double lr = L_re[j * n + k];
            double li = L_im[j * n + k];
            s -= lr * lr + li * li;
        }
        if (s <= 0.0) return -1;
        double djj = sqrt(s);
        L_re[j * n + j] = djj;
        L_im[j * n + j] = 0.0;
        double inv_djj = 1.0 / djj;
        for (size_t i = j + 1; i < n; i++) {
            double tr = B_re[i * n + j];
            double ti = B_im[i * n + j];
            for (size_t k = 0; k < j; k++) {
                double lir = L_re[i * n + k];
                double lii = L_im[i * n + k];
                double ljr = L_re[j * n + k];
                double lji = L_im[j * n + k];
                /* L[i,k] * conj(L[j,k]) = (lir+i lii)(ljr - i lji). */
                tr -= lir * ljr + lii * lji;
                ti -= lii * ljr - lir * lji;
            }
            L_re[i * n + j] = tr * inv_djj;
            L_im[i * n + j] = ti * inv_djj;
        }
    }
    return 0;
}

/* In-place solve L X = B, L complex lower-triangular n*n (strict upper
 * triangle ignored), B = X row-major n x m complex. */
static void feast_complex_trsm_left_lower(const double* L_re,
                                            const double* L_im,
                                            double* B_re, double* B_im,
                                            size_t n, size_t m) {
    for (size_t i = 0; i < n; i++) {
        double lir = L_re[i * n + i];
        double lii = L_im[i * n + i];
        double denom = lir * lir + lii * lii;
        for (size_t col = 0; col < m; col++) {
            double sr = B_re[i * m + col];
            double si = B_im[i * m + col];
            for (size_t k = 0; k < i; k++) {
                double lr = L_re[i * n + k];
                double lim = L_im[i * n + k];
                double br = B_re[k * m + col];
                double bi = B_im[k * m + col];
                sr -= lr * br - lim * bi;
                si -= lr * bi + lim * br;
            }
            B_re[i * m + col] = (sr * lir + si * lii) / denom;
            B_im[i * m + col] = (si * lir - sr * lii) / denom;
        }
    }
}

/* In-place solve L^H X = B, L complex lower-triangular n*n (so L^H is
 * upper-triangular), B = X row-major n x m complex.  Solve bottom-up.
 * (L^H)[ii, k] = conj(L[k, ii]); (L^H)[ii, ii] = conj(L[ii, ii]). */
static void feast_complex_trsm_left_lower_H(const double* L_re,
                                              const double* L_im,
                                              double* B_re, double* B_im,
                                              size_t n, size_t m) {
    for (size_t ii = n; ii-- > 0; ) {
        double lir = L_re[ii * n + ii];
        double lii = L_im[ii * n + ii];
        double denom = lir * lir + lii * lii;
        for (size_t col = 0; col < m; col++) {
            double sr = B_re[ii * m + col];
            double si = B_im[ii * m + col];
            for (size_t k = ii + 1; k < n; k++) {
                double lr  = L_re[k * n + ii];
                double lim = L_im[k * n + ii];
                double br = B_re[k * m + col];
                double bi = B_im[k * m + col];
                /* (lr - i lim)(br + i bi) = (lr*br + lim*bi) + i(lr*bi - lim*br) */
                sr -= lr * br + lim * bi;
                si -= lr * bi - lim * br;
            }
            /* (sr + i si) / (lir - i lii). */
            B_re[ii * m + col] = (sr * lir - si * lii) / denom;
            B_im[ii * m + col] = (si * lir + sr * lii) / denom;
        }
    }
}

/* Build a List of List of complex eigenvectors from a subset of a dim
 * x count matrix.  Mirrors feast_build_real_eigenvector_list_subset
 * but emits Complex[re, im] heads for entries with non-zero imag
 * part, falling back to Real for purely-real components (matching
 * direct_build_complex_hermitian_eigvec_list's convention). */
static Expr* feast_build_complex_eigenvector_list_subset(const double* V_re,
                                                          const double* V_im,
                                                          size_t dim,
                                                          size_t count,
                                                          const size_t* perm) {
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * count);
    for (size_t k = 0; k < count; k++) {
        size_t col = perm[k];
        Expr** comps = (Expr**)malloc(sizeof(Expr*) * dim);
        for (size_t i = 0; i < dim; i++) {
            double rv = V_re[i * count + col];
            double iv = V_im[i * count + col];
            if (iv == 0.0) {
                comps[i] = expr_new_real(rv);
            } else {
                Expr** args = (Expr**)malloc(sizeof(Expr*) * 2);
                args[0] = expr_new_real(rv);
                args[1] = expr_new_real(iv);
                comps[i] = expr_new_function(expr_new_symbol(SYM_Complex), args, 2);
                free(args);
            }
        }
        rows[k] = expr_new_function(expr_new_symbol(SYM_List), comps, dim);
        free(comps);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), rows, count);
    free(rows);
    return out;
}

/* Complex-Hermitian FEAST machine-precision kernel.
 *
 * Mathematical structure parallels feast_real_sym_machine, with three
 * differences forced by A being complex Hermitian (not real
 * symmetric):
 *
 *   1.  No Schwarz symmetry.  For real A, complex conjugation of the
 *       resolvent gives (z_bar I - A)^{-1} = conj((zI - A)^{-1}), so
 *       the lower-half integral is the conjugate of the upper-half
 *       and only N_e/2 solves are needed; for Hermitian (not real) A
 *       that identity fails.  We therefore do GL quadrature on the
 *       *full* contour t in [0, 2pi] via t_e = pi(1 + xi_e), with
 *       w_e = (1/2) eta_e r exp(i t_e).  Same N_e solves per
 *       iteration as the real-symmetric kernel; slightly worse
 *       quadrature accuracy at fixed N_e (mitigated by the user
 *       bumping ContourPoints if needed).
 *
 *   2.  Q is complex (n x m0); A_q = Q^* A Q and B_q = Q^* Q are
 *       complex Hermitian m_0 x m_0.
 *
 *   3.  The inner Rayleigh-Ritz solve uses the existing complex-
 *       Hermitian Direct primitives (direct_tridiag_complex_hermitian
 *       + direct_phase_correct_tridiag + direct_symtridiag_qr +
 *       direct_compose_complex_Q_real_Z) instead of the real-only
 *       chain used by Phase 2.
 *
 * Returns NULL on any failure so the caller cascades to Direct. */
static Expr* feast_complex_hermitian_machine(const MatD* A,
                                               MateigenWant want,
                                               Expr* k_spec,
                                               const FeastOpts* opts) {
    size_t n = A->n;
    if (n == 0) return NULL;
    if (!A->is_complex) return NULL;
    if (opts->interval_high <= opts->interval_low) return NULL;

    size_t m0;
    if (opts->subspace_size > 0) {
        m0 = (size_t)opts->subspace_size;
    } else {
        size_t auto_m = n / 4;
        if (auto_m < 20) auto_m = 20;
        m0 = auto_m;
    }
    if (m0 > n) m0 = n;
    if (m0 == 0) return NULL;

    bool   want_Q   = (want & MATEIGEN_WANT_VECTORS) != 0;
    int64_t max_iter = opts->max_iterations > 0 ? opts->max_iterations : 20;
    double  tol      = (opts->tolerance_given && opts->tolerance > 0.0)
                         ? opts->tolerance : 1e-12;
    double  a_lo = opts->interval_low;
    double  b_hi = opts->interval_high;
    double  c    = 0.5 * (a_lo + b_hi);
    double  r    = 0.5 * (b_hi - a_lo);

    const double* gl_xi = NULL;
    const double* gl_w  = NULL;
    int64_t Ne = feast_gl_lookup(opts->contour_points, &gl_xi, &gl_w);

    /* Workspace.  Most buffers are doubled (re + im) compared to
     * Phase 2.  The inner-solve scratch (u, v, q, sub) follows the
     * shape used by direct_tridiag_complex_hermitian. */
    double* Y_re   = (double*)malloc(sizeof(double) * n * m0);
    double* Y_im   = (double*)malloc(sizeof(double) * n * m0);
    double* Q_re   = (double*)malloc(sizeof(double) * n * m0);
    double* Q_im   = (double*)malloc(sizeof(double) * n * m0);
    double* AQ_re  = (double*)malloc(sizeof(double) * n * m0);
    double* AQ_im  = (double*)malloc(sizeof(double) * n * m0);
    double* Xh_re  = (double*)malloc(sizeof(double) * n * m0);
    double* Xh_im  = (double*)malloc(sizeof(double) * n * m0);
    double* Aq_re  = (double*)malloc(sizeof(double) * m0 * m0);
    double* Aq_im  = (double*)malloc(sizeof(double) * m0 * m0);
    double* Bq_re  = (double*)malloc(sizeof(double) * m0 * m0);
    double* Bq_im  = (double*)malloc(sizeof(double) * m0 * m0);
    double* L_re   = (double*)malloc(sizeof(double) * m0 * m0);
    double* L_im   = (double*)malloc(sizeof(double) * m0 * m0);
    double* H_re   = (double*)malloc(sizeof(double) * m0 * m0);
    double* H_im   = (double*)malloc(sizeof(double) * m0 * m0);
    double* tmp_re = (double*)malloc(sizeof(double) * m0 * m0);
    double* tmp_im = (double*)malloc(sizeof(double) * m0 * m0);
    double* W_re   = (double*)malloc(sizeof(double) * m0 * m0);
    double* W_im   = (double*)malloc(sizeof(double) * m0 * m0);
    double* Z      = (double*)malloc(sizeof(double) * m0 * m0);   /* real */
    double* V_re   = (double*)malloc(sizeof(double) * m0 * m0);
    double* V_im   = (double*)malloc(sizeof(double) * m0 * m0);
    double* diag_h = (double*)malloc(sizeof(double) * m0);
    double* sub_re = (double*)calloc(m0 ? m0 : 1, sizeof(double));
    double* sub_im = (double*)calloc(m0 ? m0 : 1, sizeof(double));
    double* u_re   = (double*)malloc(sizeof(double) * m0);
    double* u_im   = (double*)malloc(sizeof(double) * m0);
    double* v_re   = (double*)malloc(sizeof(double) * m0);
    double* v_im   = (double*)malloc(sizeof(double) * m0);
    double* qs_re  = (double*)malloc(sizeof(double) * m0);
    double* qs_im  = (double*)malloc(sizeof(double) * m0);
    double* M_re   = (double*)malloc(sizeof(double) * n * n);
    double* M_im   = (double*)malloc(sizeof(double) * n * n);
    int*    piv    = (int*)malloc(sizeof(int) * n);
    double* x_re   = (double*)malloc(sizeof(double) * n);
    double* x_im   = (double*)malloc(sizeof(double) * n);
    double* mu     = (double*)malloc(sizeof(double) * m0);
    size_t* keep   = (size_t*)malloc(sizeof(size_t) * m0);

    /* Seed Y with independent re/im LCG streams.  Using the same seed
     * for both would make Y_re == Y_im, collapsing the effective
     * subspace dimension. */
    feast_seed_subspace(Y_re, n, m0, /*seed=*/ 1442695040888963407ull);
    feast_seed_subspace(Y_im, n, m0, /*seed=*/ 6364136223846793005ull);

    Expr*  out      = NULL;
    bool   converged = false;
    size_t kept     = 0;
    const char* fail_reason = NULL;

    for (int64_t iter = 0; iter < max_iter && !converged; iter++) {
        /* ---- Step 1: Q = sum_e w_e (z_e I - A)^{-1} Y  (full contour) ---- */
        memset(Q_re, 0, sizeof(double) * n * m0);
        memset(Q_im, 0, sizeof(double) * n * m0);
        for (int64_t e = 0; e < Ne; e++) {
            double t  = M_PI * (1.0 + gl_xi[e]);
            double ct = cos(t);
            double st = sin(t);
            double z_re = c + r * ct;
            double z_im =     r * st;
            double w_re = 0.5 * gl_w[e] * r * ct;
            double w_im = 0.5 * gl_w[e] * r * st;

            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                    M_re[i * n + j] = -A->re[i * n + j];
                    M_im[i * n + j] = -A->im[i * n + j];
                }
                M_re[i * n + i] += z_re;
                M_im[i * n + i] += z_im;
            }
            if (feast_complex_lu_factor(M_re, M_im, n, piv) != 0) {
                fail_reason = "LU breakdown at a quadrature node";
                goto cleanup;
            }
            for (size_t col = 0; col < m0; col++) {
                for (size_t i = 0; i < n; i++) {
                    x_re[i] = Y_re[i * m0 + col];
                    x_im[i] = Y_im[i * m0 + col];
                }
                feast_complex_lu_solve(M_re, M_im, piv, n, x_re, x_im);
                /* Q[:, col] += w_e * x_e (complex multiplication). */
                for (size_t i = 0; i < n; i++) {
                    double xr = x_re[i], xi = x_im[i];
                    Q_re[i * m0 + col] += w_re * xr - w_im * xi;
                    Q_im[i * m0 + col] += w_re * xi + w_im * xr;
                }
            }
        }

        /* ---- Step 2: AQ = A Q,  A_q = Q^* AQ,  B_q = Q^* Q ---- */
        for (size_t i = 0; i < n; i++) {
            for (size_t col = 0; col < m0; col++) {
                double sr = 0.0, si = 0.0;
                for (size_t k2 = 0; k2 < n; k2++) {
                    double ar = A->re[i * n + k2];
                    double ai = A->im[i * n + k2];
                    double qr = Q_re[k2 * m0 + col];
                    double qi = Q_im[k2 * m0 + col];
                    sr += ar * qr - ai * qi;
                    si += ar * qi + ai * qr;
                }
                AQ_re[i * m0 + col] = sr;
                AQ_im[i * m0 + col] = si;
            }
        }
        for (size_t i = 0; i < m0; i++) {
            for (size_t j = 0; j < m0; j++) {
                double sr = 0.0, si = 0.0;
                for (size_t k2 = 0; k2 < n; k2++) {
                    /* (Q^*)[i, k2] = conj(Q[k2, i]) = qr - i qi. */
                    double qr = Q_re[k2 * m0 + i];
                    double qi = Q_im[k2 * m0 + i];
                    double ar = AQ_re[k2 * m0 + j];
                    double ai = AQ_im[k2 * m0 + j];
                    /* (qr - i qi)(ar + i ai) = qr*ar + qi*ai + i(qr*ai - qi*ar) */
                    sr += qr * ar + qi * ai;
                    si += qr * ai - qi * ar;
                }
                Aq_re[i * m0 + j] = sr;
                Aq_im[i * m0 + j] = si;
            }
        }
        for (size_t i = 0; i < m0; i++) {
            for (size_t j = 0; j < m0; j++) {
                double sr = 0.0, si = 0.0;
                for (size_t k2 = 0; k2 < n; k2++) {
                    double qir = Q_re[k2 * m0 + i];
                    double qii = Q_im[k2 * m0 + i];
                    double qjr = Q_re[k2 * m0 + j];
                    double qji = Q_im[k2 * m0 + j];
                    /* conj(Q[k2,i]) * Q[k2,j] = (qir - i qii)(qjr + i qji) */
                    sr += qir * qjr + qii * qji;
                    si += qir * qji - qii * qjr;
                }
                Bq_re[i * m0 + j] = sr;
                Bq_im[i * m0 + j] = si;
            }
        }
        /* Hermitianise: A_ij <- (A_ij + conj(A_ji))/2.  Diagonal imag
         * parts are exactly zero for a Hermitian matrix; mop up
         * arithmetic drift. */
        for (size_t i = 0; i < m0; i++) {
            Aq_im[i * m0 + i] = 0.0;
            Bq_im[i * m0 + i] = 0.0;
            for (size_t j = i + 1; j < m0; j++) {
                double r2 = 0.5 * (Aq_re[i * m0 + j] + Aq_re[j * m0 + i]);
                Aq_re[i * m0 + j] = r2;
                Aq_re[j * m0 + i] = r2;
                double i2 = 0.5 * (Aq_im[i * m0 + j] - Aq_im[j * m0 + i]);
                Aq_im[i * m0 + j] = i2;
                Aq_im[j * m0 + i] = -i2;

                r2 = 0.5 * (Bq_re[i * m0 + j] + Bq_re[j * m0 + i]);
                Bq_re[i * m0 + j] = r2;
                Bq_re[j * m0 + i] = r2;
                i2 = 0.5 * (Bq_im[i * m0 + j] - Bq_im[j * m0 + i]);
                Bq_im[i * m0 + j] = i2;
                Bq_im[j * m0 + i] = -i2;
            }
        }

        /* ---- Step 3: Cholesky B_q = L L^*,  H := L^{-1} A_q L^{-*} ---- */
        if (feast_complex_cholesky(Bq_re, Bq_im, m0, L_re, L_im) != 0) {
            fail_reason = "B_q not positive definite "
                          "(subspace likely undersized or rank-deficient)";
            goto cleanup;
        }
        memcpy(H_re, Aq_re, sizeof(double) * m0 * m0);
        memcpy(H_im, Aq_im, sizeof(double) * m0 * m0);
        feast_complex_trsm_left_lower(L_re, L_im, H_re, H_im, m0, m0);
        /* Right-apply L^{-*}: H' = H L^{-*}, equivalent to
         *   (H L^{-*})^* = L^{-1} H^*  -> conj-transpose, left-solve,
         *   conj-transpose back. */
        for (size_t i = 0; i < m0; i++) {
            for (size_t j = 0; j < m0; j++) {
                tmp_re[i * m0 + j] =  H_re[j * m0 + i];
                tmp_im[i * m0 + j] = -H_im[j * m0 + i];
            }
        }
        feast_complex_trsm_left_lower(L_re, L_im, tmp_re, tmp_im, m0, m0);
        for (size_t i = 0; i < m0; i++) {
            for (size_t j = 0; j < m0; j++) {
                H_re[i * m0 + j] =  tmp_re[j * m0 + i];
                H_im[i * m0 + j] = -tmp_im[j * m0 + i];
            }
        }
        /* Hermitianise H. */
        for (size_t i = 0; i < m0; i++) {
            H_im[i * m0 + i] = 0.0;
            for (size_t j = i + 1; j < m0; j++) {
                double r2 = 0.5 * (H_re[i * m0 + j] + H_re[j * m0 + i]);
                H_re[i * m0 + j] = r2;
                H_re[j * m0 + i] = r2;
                double i2 = 0.5 * (H_im[i * m0 + j] - H_im[j * m0 + i]);
                H_im[i * m0 + j] = i2;
                H_im[j * m0 + i] = -i2;
            }
        }

        /* ---- Step 4: inner complex Hermitian eigendecomposition ----
         * direct_tridiag_complex_hermitian (Householder) -> Hermitian
         * tridiag; direct_phase_correct_tridiag rotates the sub-
         * diagonal to be real-positive; direct_symtridiag_qr finds
         * eigenvalues + a real orthogonal Z; direct_compose_complex_
         * Q_real_Z composes Q D and Z into V_h, the eigenvectors of H. */
        direct_tridiag_complex_hermitian(H_re, H_im, m0,
                                           diag_h, sub_re, sub_im,
                                           W_re, W_im, true,
                                           u_re, u_im, v_re, v_im,
                                           qs_re, qs_im);
        direct_phase_correct_tridiag(sub_re, sub_im, m0, W_re, W_im, true);
        for (size_t i = 0; i < m0; i++)
            for (size_t j = 0; j < m0; j++)
                Z[i * m0 + j] = (i == j) ? 1.0 : 0.0;
        if (direct_symtridiag_qr(diag_h, sub_re, m0, Z, true) != 0) {
            fail_reason = "inner symmetric tridiagonal QR did not converge";
            goto cleanup;
        }
        /* Eigenvectors of H in the original basis: V_h = (Q D) Z.   */
        direct_compose_complex_Q_real_Z(W_re, W_im, Z, m0, V_re, V_im);
        /* Generalised eigenvectors of (A_q, B_q): V = L^{-*} V_h. */
        feast_complex_trsm_left_lower_H(L_re, L_im, V_re, V_im, m0, m0);
        for (size_t k = 0; k < m0; k++) mu[k] = diag_h[k];

        /* Ritz vectors Xhat = Q V (complex, n x m0). */
        for (size_t i = 0; i < n; i++) {
            for (size_t col = 0; col < m0; col++) {
                double sr = 0.0, si = 0.0;
                for (size_t k2 = 0; k2 < m0; k2++) {
                    double qr = Q_re[i * m0 + k2];
                    double qi = Q_im[i * m0 + k2];
                    double vr = V_re[k2 * m0 + col];
                    double vi = V_im[k2 * m0 + col];
                    sr += qr * vr - qi * vi;
                    si += qr * vi + qi * vr;
                }
                Xh_re[i * m0 + col] = sr;
                Xh_im[i * m0 + col] = si;
            }
        }

        /* ---- Step 5: filter to [a, b] and compute residuals ---- */
        double pad = tol * fmax(fabs(a_lo), fabs(b_hi));
        if (pad < tol) pad = tol;
        size_t k_inside = 0;
        double max_rel_res = 0.0;
        for (size_t k = 0; k < m0; k++) {
            double muk = mu[k];
            if (muk < a_lo - pad) continue;
            if (muk > b_hi + pad) continue;
            double res2 = 0.0;
            double xn2  = 0.0;
            for (size_t i = 0; i < n; i++) {
                double rr = 0.0, ri = 0.0;
                for (size_t j = 0; j < n; j++) {
                    double ar = A->re[i * n + j];
                    double ai = A->im[i * n + j];
                    double xr = Xh_re[j * m0 + k];
                    double xi = Xh_im[j * m0 + k];
                    rr += ar * xr - ai * xi;
                    ri += ar * xi + ai * xr;
                }
                double xvr = Xh_re[i * m0 + k];
                double xvi = Xh_im[i * m0 + k];
                double dr = rr - muk * xvr;
                double di = ri - muk * xvi;
                res2 += dr * dr + di * di;
                xn2  += xvr * xvr + xvi * xvi;
            }
            double res = sqrt(res2);
            double xn  = sqrt(xn2 > 0.0 ? xn2 : 1.0);
            double rel = res / (xn * fmax(fabs(muk), 1.0));
            if (rel > max_rel_res) max_rel_res = rel;
            keep[k_inside++] = k;
        }

        if (k_inside > 0 && max_rel_res < tol) {
            converged = true;
            kept = k_inside;
            break;
        }
        /* Refresh Y from Xhat with per-column unit normalisation;
         * re-seed any column that collapsed to zero so distinct
         * columns stay distinct. */
        for (size_t col = 0; col < m0; col++) {
            double nrm2 = 0.0;
            for (size_t i = 0; i < n; i++) {
                double xr = Xh_re[i * m0 + col];
                double xi = Xh_im[i * m0 + col];
                nrm2 += xr * xr + xi * xi;
            }
            double inv = (nrm2 > 0.0) ? 1.0 / sqrt(nrm2) : 0.0;
            if (inv == 0.0) {
                for (size_t i = 0; i < n; i++) {
                    uint64_t s = (uint64_t)(col + 1) * 2654435761ull
                                + (uint64_t)i * 40503ull;
                    s ^= s >> 33;
                    s *= 0xff51afd7ed558ccdull;
                    s ^= s >> 33;
                    Y_re[i * m0 + col] = ((double)(uint32_t)(s >> 32)
                                            * (2.0 / 4294967295.0)) - 1.0;
                    s = (uint64_t)(col + 7) * 0x9E3779B97F4A7C15ull
                         + (uint64_t)i * 0xBF58476D1CE4E5B9ull;
                    s ^= s >> 27;
                    Y_im[i * m0 + col] = ((double)(uint32_t)(s >> 32)
                                            * (2.0 / 4294967295.0)) - 1.0;
                }
            } else {
                for (size_t i = 0; i < n; i++) {
                    Y_re[i * m0 + col] = Xh_re[i * m0 + col] * inv;
                    Y_im[i * m0 + col] = Xh_im[i * m0 + col] * inv;
                }
            }
        }
    }

    if (!converged) {
        if (!fail_reason) fail_reason = "MaxIterations exhausted";
        goto cleanup;
    }

    /* ---- Step 6: build the result ---- */
    {
        double* vals    = (double*)malloc(sizeof(double) * kept);
        double* vecs_re = want_Q
            ? (double*)malloc(sizeof(double) * n * kept) : NULL;
        double* vecs_im = want_Q
            ? (double*)malloc(sizeof(double) * n * kept) : NULL;
        for (size_t i = 0; i < kept; i++) {
            size_t src = keep[i];
            vals[i] = mu[src];
            if (want_Q) {
                double nrm2 = 0.0;
                for (size_t r2 = 0; r2 < n; r2++) {
                    double xr = Xh_re[r2 * m0 + src];
                    double xi = Xh_im[r2 * m0 + src];
                    nrm2 += xr * xr + xi * xi;
                }
                double inv = (nrm2 > 0.0) ? 1.0 / sqrt(nrm2) : 1.0;
                for (size_t r2 = 0; r2 < n; r2++) {
                    vecs_re[r2 * kept + i] = Xh_re[r2 * m0 + src] * inv;
                    vecs_im[r2 * kept + i] = Xh_im[r2 * m0 + src] * inv;
                }
            }
        }
        size_t* perm = (size_t*)malloc(sizeof(size_t) * kept);
        direct_sort_perm_desc_abs(vals, kept, perm);
        out = want_Q
            ? feast_build_complex_eigenvector_list_subset(vecs_re, vecs_im,
                                                           n, kept, perm)
            : feast_build_real_eigenvalue_list_subset(vals, kept, perm);
        free(perm);
        free(vals);
        if (vecs_re) free(vecs_re);
        if (vecs_im) free(vecs_im);
        out = direct_apply_k_spec_list(out, k_spec);
    }

cleanup:
    free(Y_re); free(Y_im);
    free(Q_re); free(Q_im);
    free(AQ_re); free(AQ_im);
    free(Xh_re); free(Xh_im);
    free(Aq_re); free(Aq_im);
    free(Bq_re); free(Bq_im);
    free(L_re);  free(L_im);
    free(H_re);  free(H_im);
    free(tmp_re); free(tmp_im);
    free(W_re); free(W_im);
    free(Z);
    free(V_re); free(V_im);
    free(diag_h); free(sub_re); free(sub_im);
    free(u_re); free(u_im);
    free(v_re); free(v_im);
    free(qs_re); free(qs_im);
    free(M_re); free(M_im);
    free(piv);
    free(x_re); free(x_im);
    free(mu); free(keep);

    if (!converged && fail_reason) feast_warn_fallback(fail_reason);
    return out;
}

/* Top-level machine-precision FEAST dispatcher.  Phases 2 and 3 cover
 * real-symmetric and complex-Hermitian inputs respectively; any other
 * shape (non-Hermitian or generalised) is rejected upstream by
 * feast_dispatch. */
static Expr* feast_dispatch_machine(Expr* m, Expr* a, int64_t n,
                                      MateigenWant want, Expr* k_spec,
                                      const FeastOpts* opts) {
    (void)a;
    MatD A;
    if (!matD_load(m, (size_t)n, &A)) return NULL;
    double norm = A.is_complex
        ? matD_norm_inf_complex(&A)
        : matD_norm_inf_real(A.re, A.n);
    double herm_tol = (norm + 1.0) * 1e-10;

    Expr* out = NULL;
    if (!A.is_complex) {
        if (matD_is_real_symmetric(&A, herm_tol)) {
            out = feast_real_sym_machine(&A, want, k_spec, opts);
        }
    } else if (matD_is_hermitian(&A, herm_tol)) {
        out = feast_complex_hermitian_machine(&A, want, k_spec, opts);
    }
    matD_free(&A);
    return out;
}
#ifdef USE_MPFR

/* =====================================================================
 *  Phase 4 -- MPFR-precision FEAST kernels
 *
 *  Mathematical structure mirrors the machine-precision Phase 2/3
 *  kernels exactly; each double operation is replaced by the
 *  corresponding mpfr_* call rounded MPFR_RNDN.  The two structural
 *  changes compared to the machine variants are:
 *
 *    1. Gauss-Legendre nodes/weights are computed at runtime via
 *       Newton iteration on the Legendre polynomial recurrence, at the
 *       working precision.  Shipping precomputed tables would couple
 *       FEAST output to a fixed bit-count.
 *
 *    2. Inner Hermitian eigendecomposition routes through the existing
 *       direct_tridiag_real_sym_M / direct_tridiag_complex_hermitian_M
 *       + direct_symtridiag_qr_M pipeline rather than the machine
 *       direct_tridiag_real_sym path.
 * ===================================================================== */

/* Complex LU with partial pivoting at MPFR precision.  Paired
 * (re, im) row-major n*n arrays; both must be pre-init2'd to `bits`.
 * Returns 0 on success, -1 on a zero pivot (singular at z_e). */
static int feast_complex_lu_factor_M(mpfr_t* A_re, mpfr_t* A_im, size_t n,
                                       mpfr_prec_t bits, int* piv) {
    if (n == 0) return 0;
    mpfr_t mag, cur, denom, t_re, t_im, fr, fi, xr, xi, prod;
    mpfr_init2(mag,   bits); mpfr_init2(cur,   bits);
    mpfr_init2(denom, bits);
    mpfr_init2(t_re,  bits); mpfr_init2(t_im,  bits);
    mpfr_init2(fr,    bits); mpfr_init2(fi,    bits);
    mpfr_init2(xr,    bits); mpfr_init2(xi,    bits);
    mpfr_init2(prod,  bits);

    int ret = 0;
    for (size_t k = 0; k < n; k++) {
        size_t pivot = k;
        mpfr_hypot(mag, A_re[k * n + k], A_im[k * n + k], MPFR_RNDN);
        for (size_t i = k + 1; i < n; i++) {
            mpfr_hypot(cur, A_re[i * n + k], A_im[i * n + k], MPFR_RNDN);
            if (mpfr_cmp(cur, mag) > 0) { mpfr_set(mag, cur, MPFR_RNDN); pivot = i; }
        }
        if (mpfr_zero_p(mag)) { ret = -1; goto done; }
        piv[k] = (int)pivot;
        if (pivot != k) {
            for (size_t j = 0; j < n; j++) {
                mpfr_swap(A_re[k * n + j], A_re[pivot * n + j]);
                mpfr_swap(A_im[k * n + j], A_im[pivot * n + j]);
            }
        }
        /* |a_kk|^2 = pkr^2 + pki^2. */
        mpfr_mul(prod,  A_re[k * n + k], A_re[k * n + k], MPFR_RNDN);
        mpfr_mul(denom, A_im[k * n + k], A_im[k * n + k], MPFR_RNDN);
        mpfr_add(denom, denom, prod, MPFR_RNDN);
        for (size_t i = k + 1; i < n; i++) {
            /* factor = a_ik / a_kk = (ar+i ai)*(pkr - i pki)/|pk|^2. */
            mpfr_mul(fr,   A_re[i * n + k], A_re[k * n + k], MPFR_RNDN);
            mpfr_mul(prod, A_im[i * n + k], A_im[k * n + k], MPFR_RNDN);
            mpfr_add(fr,   fr, prod, MPFR_RNDN);
            mpfr_mul(fi,   A_im[i * n + k], A_re[k * n + k], MPFR_RNDN);
            mpfr_mul(prod, A_re[i * n + k], A_im[k * n + k], MPFR_RNDN);
            mpfr_sub(fi,   fi, prod, MPFR_RNDN);
            mpfr_div(fr, fr, denom, MPFR_RNDN);
            mpfr_div(fi, fi, denom, MPFR_RNDN);
            mpfr_set(A_re[i * n + k], fr, MPFR_RNDN);
            mpfr_set(A_im[i * n + k], fi, MPFR_RNDN);
            for (size_t j = k + 1; j < n; j++) {
                /* A[i,j] -= f * A[k,j]; complex mult. */
                mpfr_set(xr, A_re[k * n + j], MPFR_RNDN);
                mpfr_set(xi, A_im[k * n + j], MPFR_RNDN);
                mpfr_mul(t_re, fr, xr, MPFR_RNDN);
                mpfr_mul(prod, fi, xi, MPFR_RNDN);
                mpfr_sub(t_re, t_re, prod, MPFR_RNDN);
                mpfr_mul(t_im, fr, xi, MPFR_RNDN);
                mpfr_mul(prod, fi, xr, MPFR_RNDN);
                mpfr_add(t_im, t_im, prod, MPFR_RNDN);
                mpfr_sub(A_re[i * n + j], A_re[i * n + j], t_re, MPFR_RNDN);
                mpfr_sub(A_im[i * n + j], A_im[i * n + j], t_im, MPFR_RNDN);
            }
        }
    }
done:
    mpfr_clear(mag); mpfr_clear(cur); mpfr_clear(denom);
    mpfr_clear(t_re); mpfr_clear(t_im);
    mpfr_clear(fr);  mpfr_clear(fi);
    mpfr_clear(xr);  mpfr_clear(xi);
    mpfr_clear(prod);
    return ret;
}

/* Solve LUx = Pb in place over a single complex RHS (length n) at MPFR
 * precision.  Unit-diagonal L occupies the strict lower triangle of LU;
 * U occupies the upper triangle including the diagonal. */
static void feast_complex_lu_solve_M(const mpfr_t* LU_re, const mpfr_t* LU_im,
                                       const int* piv, size_t n,
                                       mpfr_prec_t bits,
                                       mpfr_t* b_re, mpfr_t* b_im) {
    mpfr_t sr, si, prod, dr, di, denom;
    mpfr_init2(sr,    bits); mpfr_init2(si,    bits);
    mpfr_init2(prod,  bits);
    mpfr_init2(dr,    bits); mpfr_init2(di,    bits);
    mpfr_init2(denom, bits);

    for (size_t k = 0; k < n; k++) {
        int p = piv[k];
        if ((size_t)p != k) {
            mpfr_swap(b_re[k], b_re[p]);
            mpfr_swap(b_im[k], b_im[p]);
        }
    }
    /* Forward solve L y = Pb (unit-lower L). */
    for (size_t i = 1; i < n; i++) {
        mpfr_set(sr, b_re[i], MPFR_RNDN);
        mpfr_set(si, b_im[i], MPFR_RNDN);
        for (size_t j = 0; j < i; j++) {
            /* sr -= lr * yr - li * yi;  si -= lr * yi + li * yr. */
            mpfr_mul(prod, LU_re[i * n + j], b_re[j], MPFR_RNDN);
            mpfr_sub(sr, sr, prod, MPFR_RNDN);
            mpfr_mul(prod, LU_im[i * n + j], b_im[j], MPFR_RNDN);
            mpfr_add(sr, sr, prod, MPFR_RNDN);
            mpfr_mul(prod, LU_re[i * n + j], b_im[j], MPFR_RNDN);
            mpfr_sub(si, si, prod, MPFR_RNDN);
            mpfr_mul(prod, LU_im[i * n + j], b_re[j], MPFR_RNDN);
            mpfr_sub(si, si, prod, MPFR_RNDN);
        }
        mpfr_set(b_re[i], sr, MPFR_RNDN);
        mpfr_set(b_im[i], si, MPFR_RNDN);
    }
    /* Backward solve U x = y. */
    for (size_t ii = n; ii-- > 0; ) {
        mpfr_set(sr, b_re[ii], MPFR_RNDN);
        mpfr_set(si, b_im[ii], MPFR_RNDN);
        for (size_t j = ii + 1; j < n; j++) {
            mpfr_mul(prod, LU_re[ii * n + j], b_re[j], MPFR_RNDN);
            mpfr_sub(sr, sr, prod, MPFR_RNDN);
            mpfr_mul(prod, LU_im[ii * n + j], b_im[j], MPFR_RNDN);
            mpfr_add(sr, sr, prod, MPFR_RNDN);
            mpfr_mul(prod, LU_re[ii * n + j], b_im[j], MPFR_RNDN);
            mpfr_sub(si, si, prod, MPFR_RNDN);
            mpfr_mul(prod, LU_im[ii * n + j], b_re[j], MPFR_RNDN);
            mpfr_sub(si, si, prod, MPFR_RNDN);
        }
        mpfr_set(dr, LU_re[ii * n + ii], MPFR_RNDN);
        mpfr_set(di, LU_im[ii * n + ii], MPFR_RNDN);
        mpfr_mul(denom, dr, dr, MPFR_RNDN);
        mpfr_mul(prod,  di, di, MPFR_RNDN);
        mpfr_add(denom, denom, prod, MPFR_RNDN);
        /* (sr + i si) / (dr + i di) = (sr dr + si di + i(si dr - sr di))/|d|^2 */
        mpfr_mul(prod, sr, dr, MPFR_RNDN);
        mpfr_mul(b_re[ii], si, di, MPFR_RNDN);
        mpfr_add(b_re[ii], b_re[ii], prod, MPFR_RNDN);
        mpfr_div(b_re[ii], b_re[ii], denom, MPFR_RNDN);
        mpfr_mul(prod, si, dr, MPFR_RNDN);
        mpfr_mul(b_im[ii], sr, di, MPFR_RNDN);
        mpfr_sub(b_im[ii], prod, b_im[ii], MPFR_RNDN);
        mpfr_div(b_im[ii], b_im[ii], denom, MPFR_RNDN);
    }
    mpfr_clear(sr); mpfr_clear(si);
    mpfr_clear(prod);
    mpfr_clear(dr); mpfr_clear(di);
    mpfr_clear(denom);
}

/* Compute Gauss-Legendre nodes / weights of order Ne on [-1, 1] at the
 * working precision.  Newton iteration on the Legendre polynomial
 * recurrence:
 *
 *     P_0(x) = 1, P_1(x) = x,
 *     (k+1) P_{k+1}(x) = (2k+1) x P_k(x) - k P_{k-1}(x),
 *     P_n'(x) = n (x P_n(x) - P_{n-1}(x)) / (x^2 - 1),
 *     w_k = 2 / ((1 - x_k^2) (P_n'(x_k))^2).
 *
 * Initial guess for the k-th root (1-indexed): x_k = cos(pi (k - 1/4) /
 * (n + 1/2)), which converges quadratically.  `xi[]` / `w[]` must be
 * pre-init2'd to `bits` and have length Ne. */
static void feast_gl_compute_mpfr(int64_t Ne, mpfr_prec_t bits,
                                    mpfr_t* xi, mpfr_t* w) {
    if (Ne <= 0) return;
    mpfr_t x, x_new, P_km1, P_k, P_kp1, Pp, tmp, num, den;
    mpfr_t pi_val;
    mpfr_init2(x,     bits); mpfr_init2(x_new, bits);
    mpfr_init2(P_km1, bits); mpfr_init2(P_k,   bits);
    mpfr_init2(P_kp1, bits); mpfr_init2(Pp,    bits);
    mpfr_init2(tmp,   bits); mpfr_init2(num,   bits); mpfr_init2(den, bits);
    mpfr_init2(pi_val, bits);
    mpfr_const_pi(pi_val, MPFR_RNDN);

    /* Convergence threshold: ~ 2^(-bits/2 - 4).  Newton converges
     * quadratically, so a few iterations after this still buy precision. */
    mpfr_t conv_tol;
    mpfr_init2(conv_tol, bits);
    mpfr_set_ui(conv_tol, 1, MPFR_RNDN);
    mpfr_div_2si(conv_tol, conv_tol, (long)(bits / 2 + 4), MPFR_RNDN);

    const int max_iter = 100;

    for (int64_t k = 1; k <= Ne; k++) {
        /* Initial guess: x = cos(pi (k - 1/4) / (Ne + 1/2)). */
        mpfr_set_d(num, (double)k - 0.25, MPFR_RNDN);
        mpfr_set_d(den, (double)Ne + 0.5, MPFR_RNDN);
        mpfr_mul(tmp, pi_val, num, MPFR_RNDN);
        mpfr_div(tmp, tmp, den, MPFR_RNDN);
        mpfr_cos(x, tmp, MPFR_RNDN);

        for (int it = 0; it < max_iter; it++) {
            /* Evaluate P_Ne(x), P_{Ne-1}(x) via recurrence. */
            mpfr_set_ui(P_km1, 1, MPFR_RNDN);
            mpfr_set(P_k, x, MPFR_RNDN);
            for (int64_t i = 1; i < Ne; i++) {
                /* P_{i+1} = ((2i+1) x P_i - i P_{i-1}) / (i+1). */
                mpfr_mul(tmp, x, P_k, MPFR_RNDN);
                mpfr_mul_si(tmp, tmp, 2 * i + 1, MPFR_RNDN);
                mpfr_mul_si(num, P_km1, i, MPFR_RNDN);
                mpfr_sub(P_kp1, tmp, num, MPFR_RNDN);
                mpfr_div_si(P_kp1, P_kp1, i + 1, MPFR_RNDN);
                mpfr_set(P_km1, P_k, MPFR_RNDN);
                mpfr_set(P_k, P_kp1, MPFR_RNDN);
            }
            /* P_Ne = P_k; P_{Ne-1} = P_km1.
             * Pp = Ne (x P_Ne - P_{Ne-1}) / (x^2 - 1). */
            mpfr_mul(num, x, P_k, MPFR_RNDN);
            mpfr_sub(num, num, P_km1, MPFR_RNDN);
            mpfr_mul_si(num, num, Ne, MPFR_RNDN);
            mpfr_mul(den, x, x, MPFR_RNDN);
            mpfr_sub_ui(den, den, 1, MPFR_RNDN);
            mpfr_div(Pp, num, den, MPFR_RNDN);

            /* Newton step: x_new = x - P_Ne / Pp. */
            mpfr_div(tmp, P_k, Pp, MPFR_RNDN);
            mpfr_sub(x_new, x, tmp, MPFR_RNDN);

            mpfr_sub(tmp, x_new, x, MPFR_RNDN);
            mpfr_abs(tmp, tmp, MPFR_RNDN);
            mpfr_set(x, x_new, MPFR_RNDN);
            if (mpfr_cmp(tmp, conv_tol) < 0) break;
        }

        /* Re-evaluate P_Ne and its derivative at the converged root. */
        mpfr_set_ui(P_km1, 1, MPFR_RNDN);
        mpfr_set(P_k, x, MPFR_RNDN);
        for (int64_t i = 1; i < Ne; i++) {
            mpfr_mul(tmp, x, P_k, MPFR_RNDN);
            mpfr_mul_si(tmp, tmp, 2 * i + 1, MPFR_RNDN);
            mpfr_mul_si(num, P_km1, i, MPFR_RNDN);
            mpfr_sub(P_kp1, tmp, num, MPFR_RNDN);
            mpfr_div_si(P_kp1, P_kp1, i + 1, MPFR_RNDN);
            mpfr_set(P_km1, P_k, MPFR_RNDN);
            mpfr_set(P_k, P_kp1, MPFR_RNDN);
        }
        mpfr_mul(num, x, P_k, MPFR_RNDN);
        mpfr_sub(num, num, P_km1, MPFR_RNDN);
        mpfr_mul_si(num, num, Ne, MPFR_RNDN);
        mpfr_mul(den, x, x, MPFR_RNDN);
        mpfr_sub_ui(den, den, 1, MPFR_RNDN);
        mpfr_div(Pp, num, den, MPFR_RNDN);

        /* Roots from cos(.) are ordered with largest first (k=1 ->
         * x ~ +0.97...).  Store ascending: place at slot Ne - k. */
        mpfr_set(xi[Ne - k], x, MPFR_RNDN);
        /* w = 2 / ((1 - x^2) Pp^2).  Note (1 - x^2) = -den. */
        mpfr_mul(tmp, Pp, Pp, MPFR_RNDN);
        mpfr_neg(num, den, MPFR_RNDN);
        mpfr_mul(tmp, tmp, num, MPFR_RNDN);
        mpfr_ui_div(w[Ne - k], 2, tmp, MPFR_RNDN);
    }

    mpfr_clear(x);     mpfr_clear(x_new);
    mpfr_clear(P_km1); mpfr_clear(P_k);
    mpfr_clear(P_kp1); mpfr_clear(Pp);
    mpfr_clear(tmp);   mpfr_clear(num); mpfr_clear(den);
    mpfr_clear(pi_val);
    mpfr_clear(conv_tol);
}

/* Seed Y (n*m0 row-major, mpfr_t pre-init2'd) using the same Knuth
 * LCG stream as feast_seed_subspace, but writing directly into mpfr_t.
 * Same seed values keep MPFR FEAST output close to machine FEAST. */
static void feast_seed_subspace_M(mpfr_t* Y, size_t n, size_t m0,
                                    uint64_t seed) {
    uint64_t state = seed ? seed : 0xC0FFEE05u;
    for (size_t i = 0; i < n * m0; i++) {
        state = state * 6364136223846793005ull + 1442695040888963407ull;
        uint32_t v = (uint32_t)(state >> 32);
        double d = ((double)v * (2.0 / 4294967295.0)) - 1.0;
        mpfr_set_d(Y[i], d, MPFR_RNDN);
    }
}

/* Real SPD Cholesky factorisation at MPFR precision.  L lower-triangular
 * (full storage, strict upper triangle zeroed).  Returns 0 on success,
 * -1 on non-positive pivot. */
static int feast_cholesky_M(const mpfr_t* B, size_t n, mpfr_prec_t bits,
                              mpfr_t* L) {
    for (size_t i = 0; i < n * n; i++) mpfr_set_zero(L[i], 1);
    mpfr_t s, prod, djj;
    mpfr_init2(s,    bits);
    mpfr_init2(prod, bits);
    mpfr_init2(djj,  bits);
    int ret = 0;
    for (size_t j = 0; j < n; j++) {
        mpfr_set(s, B[j * n + j], MPFR_RNDN);
        for (size_t k = 0; k < j; k++) {
            mpfr_mul(prod, L[j * n + k], L[j * n + k], MPFR_RNDN);
            mpfr_sub(s, s, prod, MPFR_RNDN);
        }
        if (mpfr_sgn(s) <= 0) { ret = -1; goto done; }
        mpfr_sqrt(djj, s, MPFR_RNDN);
        mpfr_set(L[j * n + j], djj, MPFR_RNDN);
        for (size_t i = j + 1; i < n; i++) {
            mpfr_set(s, B[i * n + j], MPFR_RNDN);
            for (size_t k = 0; k < j; k++) {
                mpfr_mul(prod, L[i * n + k], L[j * n + k], MPFR_RNDN);
                mpfr_sub(s, s, prod, MPFR_RNDN);
            }
            mpfr_div(L[i * n + j], s, djj, MPFR_RNDN);
        }
    }
done:
    mpfr_clear(s); mpfr_clear(prod); mpfr_clear(djj);
    return ret;
}

/* In-place solve L X = B at MPFR precision; L lower-triangular n*n
 * (strict upper triangle ignored), B = X row-major n x m. */
static void feast_trsm_left_lower_M(const mpfr_t* L, mpfr_t* B,
                                      size_t n, size_t m, mpfr_prec_t bits) {
    mpfr_t s, prod;
    mpfr_init2(s, bits); mpfr_init2(prod, bits);
    for (size_t i = 0; i < n; i++) {
        for (size_t col = 0; col < m; col++) {
            mpfr_set(s, B[i * m + col], MPFR_RNDN);
            for (size_t k = 0; k < i; k++) {
                mpfr_mul(prod, L[i * n + k], B[k * m + col], MPFR_RNDN);
                mpfr_sub(s, s, prod, MPFR_RNDN);
            }
            mpfr_div(B[i * m + col], s, L[i * n + i], MPFR_RNDN);
        }
    }
    mpfr_clear(s); mpfr_clear(prod);
}

/* In-place solve L^T X = B at MPFR precision. */
static void feast_trsm_left_lower_T_M(const mpfr_t* L, mpfr_t* B,
                                        size_t n, size_t m, mpfr_prec_t bits) {
    mpfr_t s, prod;
    mpfr_init2(s, bits); mpfr_init2(prod, bits);
    for (size_t ii = n; ii-- > 0; ) {
        for (size_t col = 0; col < m; col++) {
            mpfr_set(s, B[ii * m + col], MPFR_RNDN);
            for (size_t k = ii + 1; k < n; k++) {
                mpfr_mul(prod, L[k * n + ii], B[k * m + col], MPFR_RNDN);
                mpfr_sub(s, s, prod, MPFR_RNDN);
            }
            mpfr_div(B[ii * m + col], s, L[ii * n + ii], MPFR_RNDN);
        }
    }
    mpfr_clear(s); mpfr_clear(prod);
}

/* Complex Hermitian Cholesky at MPFR precision: B = L L^*, with real-
 * positive diagonal of L.  Strict upper triangle of L is zeroed.
 * Returns 0 on success, -1 on non-positive pivot. */
static int feast_complex_cholesky_M(const mpfr_t* B_re, const mpfr_t* B_im,
                                      size_t n, mpfr_prec_t bits,
                                      mpfr_t* L_re, mpfr_t* L_im) {
    for (size_t i = 0; i < n * n; i++) {
        mpfr_set_zero(L_re[i], 1);
        mpfr_set_zero(L_im[i], 1);
    }
    mpfr_t s, prod, djj, tr, ti;
    mpfr_init2(s,    bits); mpfr_init2(prod, bits); mpfr_init2(djj, bits);
    mpfr_init2(tr,   bits); mpfr_init2(ti,   bits);
    int ret = 0;
    for (size_t j = 0; j < n; j++) {
        mpfr_set(s, B_re[j * n + j], MPFR_RNDN);
        for (size_t k = 0; k < j; k++) {
            /* s -= |L[j,k]|^2 */
            mpfr_mul(prod, L_re[j * n + k], L_re[j * n + k], MPFR_RNDN);
            mpfr_sub(s, s, prod, MPFR_RNDN);
            mpfr_mul(prod, L_im[j * n + k], L_im[j * n + k], MPFR_RNDN);
            mpfr_sub(s, s, prod, MPFR_RNDN);
        }
        if (mpfr_sgn(s) <= 0) { ret = -1; goto done; }
        mpfr_sqrt(djj, s, MPFR_RNDN);
        mpfr_set(L_re[j * n + j], djj, MPFR_RNDN);
        mpfr_set_zero(L_im[j * n + j], 1);
        for (size_t i = j + 1; i < n; i++) {
            mpfr_set(tr, B_re[i * n + j], MPFR_RNDN);
            mpfr_set(ti, B_im[i * n + j], MPFR_RNDN);
            for (size_t k = 0; k < j; k++) {
                /* L[i,k] * conj(L[j,k]) = (lir + i lii)(ljr - i lji). */
                mpfr_mul(prod, L_re[i * n + k], L_re[j * n + k], MPFR_RNDN);
                mpfr_sub(tr, tr, prod, MPFR_RNDN);
                mpfr_mul(prod, L_im[i * n + k], L_im[j * n + k], MPFR_RNDN);
                mpfr_sub(tr, tr, prod, MPFR_RNDN);
                mpfr_mul(prod, L_im[i * n + k], L_re[j * n + k], MPFR_RNDN);
                mpfr_sub(ti, ti, prod, MPFR_RNDN);
                mpfr_mul(prod, L_re[i * n + k], L_im[j * n + k], MPFR_RNDN);
                mpfr_add(ti, ti, prod, MPFR_RNDN);
            }
            mpfr_div(L_re[i * n + j], tr, djj, MPFR_RNDN);
            mpfr_div(L_im[i * n + j], ti, djj, MPFR_RNDN);
        }
    }
done:
    mpfr_clear(s); mpfr_clear(prod); mpfr_clear(djj);
    mpfr_clear(tr); mpfr_clear(ti);
    return ret;
}

/* Complex in-place solve L X = B at MPFR precision; L complex lower-
 * triangular n*n (strict upper triangle ignored), B = X row-major n x m
 * complex. */
static void feast_complex_trsm_left_lower_M(const mpfr_t* L_re,
                                              const mpfr_t* L_im,
                                              mpfr_t* B_re, mpfr_t* B_im,
                                              size_t n, size_t m,
                                              mpfr_prec_t bits) {
    mpfr_t sr, si, prod, denom;
    mpfr_init2(sr,    bits); mpfr_init2(si,    bits);
    mpfr_init2(prod,  bits); mpfr_init2(denom, bits);
    for (size_t i = 0; i < n; i++) {
        mpfr_mul(denom, L_re[i * n + i], L_re[i * n + i], MPFR_RNDN);
        mpfr_mul(prod,  L_im[i * n + i], L_im[i * n + i], MPFR_RNDN);
        mpfr_add(denom, denom, prod, MPFR_RNDN);
        for (size_t col = 0; col < m; col++) {
            mpfr_set(sr, B_re[i * m + col], MPFR_RNDN);
            mpfr_set(si, B_im[i * m + col], MPFR_RNDN);
            for (size_t k = 0; k < i; k++) {
                /* (lr + i li)(br + i bi). */
                mpfr_mul(prod, L_re[i * n + k], B_re[k * m + col], MPFR_RNDN);
                mpfr_sub(sr, sr, prod, MPFR_RNDN);
                mpfr_mul(prod, L_im[i * n + k], B_im[k * m + col], MPFR_RNDN);
                mpfr_add(sr, sr, prod, MPFR_RNDN);
                mpfr_mul(prod, L_re[i * n + k], B_im[k * m + col], MPFR_RNDN);
                mpfr_sub(si, si, prod, MPFR_RNDN);
                mpfr_mul(prod, L_im[i * n + k], B_re[k * m + col], MPFR_RNDN);
                mpfr_sub(si, si, prod, MPFR_RNDN);
            }
            /* (sr + i si) / (lir + i lii). */
            mpfr_mul(prod,            sr, L_re[i * n + i], MPFR_RNDN);
            mpfr_mul(B_re[i * m + col], si, L_im[i * n + i], MPFR_RNDN);
            mpfr_add(B_re[i * m + col], B_re[i * m + col], prod, MPFR_RNDN);
            mpfr_div(B_re[i * m + col], B_re[i * m + col], denom, MPFR_RNDN);
            mpfr_mul(prod,            si, L_re[i * n + i], MPFR_RNDN);
            mpfr_mul(B_im[i * m + col], sr, L_im[i * n + i], MPFR_RNDN);
            mpfr_sub(B_im[i * m + col], prod, B_im[i * m + col], MPFR_RNDN);
            mpfr_div(B_im[i * m + col], B_im[i * m + col], denom, MPFR_RNDN);
        }
    }
    mpfr_clear(sr); mpfr_clear(si);
    mpfr_clear(prod); mpfr_clear(denom);
}

/* Complex in-place solve L^H X = B at MPFR precision; L lower-triangular
 * (so L^H is upper-triangular).  Solve bottom-up. */
static void feast_complex_trsm_left_lower_H_M(const mpfr_t* L_re,
                                                const mpfr_t* L_im,
                                                mpfr_t* B_re, mpfr_t* B_im,
                                                size_t n, size_t m,
                                                mpfr_prec_t bits) {
    mpfr_t sr, si, prod, denom;
    mpfr_init2(sr,    bits); mpfr_init2(si,    bits);
    mpfr_init2(prod,  bits); mpfr_init2(denom, bits);
    for (size_t ii = n; ii-- > 0; ) {
        mpfr_mul(denom, L_re[ii * n + ii], L_re[ii * n + ii], MPFR_RNDN);
        mpfr_mul(prod,  L_im[ii * n + ii], L_im[ii * n + ii], MPFR_RNDN);
        mpfr_add(denom, denom, prod, MPFR_RNDN);
        for (size_t col = 0; col < m; col++) {
            mpfr_set(sr, B_re[ii * m + col], MPFR_RNDN);
            mpfr_set(si, B_im[ii * m + col], MPFR_RNDN);
            for (size_t k = ii + 1; k < n; k++) {
                /* (lr - i lim)(br + i bi). */
                mpfr_mul(prod, L_re[k * n + ii], B_re[k * m + col], MPFR_RNDN);
                mpfr_sub(sr, sr, prod, MPFR_RNDN);
                mpfr_mul(prod, L_im[k * n + ii], B_im[k * m + col], MPFR_RNDN);
                mpfr_sub(sr, sr, prod, MPFR_RNDN);
                mpfr_mul(prod, L_re[k * n + ii], B_im[k * m + col], MPFR_RNDN);
                mpfr_sub(si, si, prod, MPFR_RNDN);
                mpfr_mul(prod, L_im[k * n + ii], B_re[k * m + col], MPFR_RNDN);
                mpfr_add(si, si, prod, MPFR_RNDN);
            }
            /* (sr + i si) / (lir - i lii). */
            mpfr_mul(prod,             sr, L_re[ii * n + ii], MPFR_RNDN);
            mpfr_mul(B_re[ii * m + col], si, L_im[ii * n + ii], MPFR_RNDN);
            mpfr_sub(B_re[ii * m + col], prod, B_re[ii * m + col], MPFR_RNDN);
            mpfr_div(B_re[ii * m + col], B_re[ii * m + col], denom, MPFR_RNDN);
            mpfr_mul(prod,             si, L_re[ii * n + ii], MPFR_RNDN);
            mpfr_mul(B_im[ii * m + col], sr, L_im[ii * n + ii], MPFR_RNDN);
            mpfr_add(B_im[ii * m + col], B_im[ii * m + col], prod, MPFR_RNDN);
            mpfr_div(B_im[ii * m + col], B_im[ii * m + col], denom, MPFR_RNDN);
        }
    }
    mpfr_clear(sr); mpfr_clear(si);
    mpfr_clear(prod); mpfr_clear(denom);
}

/* Sort mpfr_t values descending by |.| -- subset variant.  perm[] is
 * filled in with the desired ordering. */
static void feast_sort_perm_desc_abs_M(const mpfr_t* vals, size_t n,
                                         size_t* perm) {
    for (size_t i = 0; i < n; i++) perm[i] = i;
    if (n == 0) return;
    mpfr_t ac, ap;
    mpfr_init2(ac, mpfr_get_prec(vals[0]));
    mpfr_init2(ap, mpfr_get_prec(vals[0]));
    for (size_t i = 1; i < n; i++) {
        size_t cur = perm[i];
        mpfr_abs(ac, vals[cur], MPFR_RNDN);
        size_t j = i;
        while (j > 0) {
            mpfr_abs(ap, vals[perm[j - 1]], MPFR_RNDN);
            int cmp = mpfr_cmp(ap, ac);
            if (cmp > 0 || (cmp == 0 && perm[j - 1] < cur)) break;
            perm[j] = perm[j - 1];
            j--;
        }
        perm[j] = cur;
    }
    mpfr_clear(ac); mpfr_clear(ap);
}

/* MPFR subset builders (m kept eigenpairs of dim n vectors). */
static Expr* feast_build_real_eigenvalue_list_subset_M(const mpfr_t* vals,
                                                        size_t count,
                                                        const size_t* perm) {
    Expr** items = (Expr**)malloc(sizeof(Expr*) * count);
    for (size_t i = 0; i < count; i++) {
        items[i] = expr_new_mpfr_copy(vals[perm[i]]);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), items, count);
    free(items);
    return out;
}

static Expr* feast_build_real_eigenvector_list_subset_M(const mpfr_t* vecs,
                                                         size_t dim,
                                                         size_t count,
                                                         const size_t* perm) {
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * count);
    for (size_t k = 0; k < count; k++) {
        size_t col = perm[k];
        Expr** comps = (Expr**)malloc(sizeof(Expr*) * dim);
        for (size_t i = 0; i < dim; i++) {
            comps[i] = expr_new_mpfr_copy(vecs[i * count + col]);
        }
        rows[k] = expr_new_function(expr_new_symbol(SYM_List), comps, dim);
        free(comps);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), rows, count);
    free(rows);
    return out;
}

static Expr* feast_build_complex_eigenvector_list_subset_M(const mpfr_t* V_re,
                                                            const mpfr_t* V_im,
                                                            size_t dim,
                                                            size_t count,
                                                            const size_t* perm) {
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * count);
    for (size_t k = 0; k < count; k++) {
        size_t col = perm[k];
        Expr** comps = (Expr**)malloc(sizeof(Expr*) * dim);
        for (size_t i = 0; i < dim; i++) {
            if (mpfr_zero_p(V_im[i * count + col])) {
                comps[i] = expr_new_mpfr_copy(V_re[i * count + col]);
            } else {
                Expr** args = (Expr**)malloc(sizeof(Expr*) * 2);
                args[0] = expr_new_mpfr_copy(V_re[i * count + col]);
                args[1] = expr_new_mpfr_copy(V_im[i * count + col]);
                comps[i] = expr_new_function(expr_new_symbol(SYM_Complex),
                                              args, 2);
                free(args);
            }
        }
        rows[k] = expr_new_function(expr_new_symbol(SYM_List), comps, dim);
        free(comps);
    }
    Expr* out = expr_new_function(expr_new_symbol(SYM_List), rows, count);
    free(rows);
    return out;
}

/* Default tolerance for MPFR FEAST: 2^(-0.7 * bits).  Matches the
 * planning doc; ~ 10^-28 at 128 bits, 10^-54 at 256 bits. */
static void feast_default_tolerance_M(mpfr_prec_t bits, mpfr_t tol) {
    mpfr_set_ui(tol, 1, MPFR_RNDN);
    long shift = (long)((double)bits * 0.7);
    if (shift < 1) shift = 1;
    mpfr_div_2si(tol, tol, shift, MPFR_RNDN);
}

/* Real-symmetric FEAST kernel at MPFR precision.
 *
 * Algorithmic structure is identical to feast_real_sym_machine; see the
 * commentary there.  The only differences are (a) every double op is an
 * mpfr_* call, and (b) the inner Rayleigh-Ritz eigendecomposition uses
 * the existing direct_tridiag_real_sym_M + direct_symtridiag_qr_M
 * pipeline. */
static Expr* feast_real_sym_mpfr(const MatM* A, MateigenWant want,
                                   Expr* k_spec, const FeastOpts* opts) {
    size_t n = A->n;
    mpfr_prec_t bits = A->bits;
    if (n == 0) return NULL;
    if (opts->interval_high <= opts->interval_low) return NULL;

    size_t m0;
    if (opts->subspace_size > 0) {
        m0 = (size_t)opts->subspace_size;
    } else {
        size_t auto_m = n / 4;
        if (auto_m < 20) auto_m = 20;
        m0 = auto_m;
    }
    if (m0 > n) m0 = n;
    if (m0 == 0) return NULL;

    bool   want_Q   = (want & MATEIGEN_WANT_VECTORS) != 0;
    int64_t max_iter = opts->max_iterations > 0 ? opts->max_iterations : 20;
    int64_t Ne = opts->contour_points > 0 ? opts->contour_points : 8;
    if (Ne != 2 && Ne != 4 && Ne != 8 && Ne != 16) Ne = 8;

    /* MPFR interval / contour centre / radius. */
    mpfr_t a_lo, b_hi, c_ctr, r_rad, tol;
    mpfr_init2(a_lo,  bits); mpfr_set_d(a_lo,  opts->interval_low,  MPFR_RNDN);
    mpfr_init2(b_hi,  bits); mpfr_set_d(b_hi,  opts->interval_high, MPFR_RNDN);
    mpfr_init2(c_ctr, bits); mpfr_init2(r_rad, bits);
    mpfr_add(c_ctr, a_lo, b_hi, MPFR_RNDN);
    mpfr_div_2si(c_ctr, c_ctr, 1, MPFR_RNDN);
    mpfr_sub(r_rad, b_hi, a_lo, MPFR_RNDN);
    mpfr_div_2si(r_rad, r_rad, 1, MPFR_RNDN);
    mpfr_init2(tol, bits);
    if (opts->tolerance_given && opts->tolerance > 0.0) {
        mpfr_set_d(tol, opts->tolerance, MPFR_RNDN);
    } else {
        feast_default_tolerance_M(bits, tol);
    }

    /* Gauss-Legendre nodes / weights at working precision. */
    mpfr_t* gl_xi = mpfr_array_alloc((size_t)Ne, bits);
    mpfr_t* gl_w  = mpfr_array_alloc((size_t)Ne, bits);
    feast_gl_compute_mpfr(Ne, bits, gl_xi, gl_w);

    /* Workspace. */
    mpfr_t* Y     = mpfr_array_alloc(n * m0, bits);
    mpfr_t* Q     = mpfr_array_alloc(n * m0, bits);
    mpfr_t* AQ    = mpfr_array_alloc(n * m0, bits);
    mpfr_t* Xhat  = mpfr_array_alloc(n * m0, bits);
    mpfr_t* Aq    = mpfr_array_alloc(m0 * m0, bits);
    mpfr_t* Bq    = mpfr_array_alloc(m0 * m0, bits);
    mpfr_t* Lchol = mpfr_array_alloc(m0 * m0, bits);
    mpfr_t* H     = mpfr_array_alloc(m0 * m0, bits);
    mpfr_t* HQ    = mpfr_array_alloc(m0 * m0, bits);
    mpfr_t* W     = mpfr_array_alloc(m0 * m0, bits);
    mpfr_t* diag_h= mpfr_array_alloc(m0, bits);
    mpfr_t* sub_h = mpfr_array_alloc(m0 ? m0 : 1, bits);
    mpfr_t* u_h   = mpfr_array_alloc(m0, bits);
    mpfr_t* p_h   = mpfr_array_alloc(m0, bits);
    mpfr_t* q_h   = mpfr_array_alloc(m0, bits);
    mpfr_t* qr_tmp= mpfr_array_alloc(12, bits);
    mpfr_t* M_re  = mpfr_array_alloc(n * n, bits);
    mpfr_t* M_im  = mpfr_array_alloc(n * n, bits);
    int*    piv   = (int*)malloc(sizeof(int) * n);
    mpfr_t* x_re  = mpfr_array_alloc(n, bits);
    mpfr_t* x_im  = mpfr_array_alloc(n, bits);
    mpfr_t* mu    = mpfr_array_alloc(m0, bits);
    size_t* keep  = (size_t*)malloc(sizeof(size_t) * m0);

    /* Inline scratch cells. */
    mpfr_t t_var, ct, st, z_re, z_im, w_re, w_im, prod, sum, sum2;
    mpfr_init2(t_var, bits); mpfr_init2(ct,   bits); mpfr_init2(st,   bits);
    mpfr_init2(z_re,  bits); mpfr_init2(z_im, bits);
    mpfr_init2(w_re,  bits); mpfr_init2(w_im, bits);
    mpfr_init2(prod,  bits); mpfr_init2(sum,  bits); mpfr_init2(sum2, bits);
    mpfr_t pi_val;
    mpfr_init2(pi_val, bits);
    mpfr_const_pi(pi_val, MPFR_RNDN);
    mpfr_t half_pi;
    mpfr_init2(half_pi, bits);
    mpfr_div_2si(half_pi, pi_val, 1, MPFR_RNDN);
    /* Padding for filter window. */
    mpfr_t pad, abs_a, abs_b;
    mpfr_init2(pad,   bits);
    mpfr_init2(abs_a, bits);
    mpfr_init2(abs_b, bits);
    mpfr_abs(abs_a, a_lo, MPFR_RNDN);
    mpfr_abs(abs_b, b_hi, MPFR_RNDN);
    if (mpfr_cmp(abs_a, abs_b) >= 0) mpfr_mul(pad, tol, abs_a, MPFR_RNDN);
    else                              mpfr_mul(pad, tol, abs_b, MPFR_RNDN);
    if (mpfr_cmp(pad, tol) < 0) mpfr_set(pad, tol, MPFR_RNDN);

    feast_seed_subspace_M(Y, n, m0, /*seed=*/ 1442695040888963407ull);

    Expr*  out      = NULL;
    bool   converged = false;
    size_t kept     = 0;
    const char* fail_reason = NULL;

    for (int64_t iter = 0; iter < max_iter && !converged; iter++) {
        /* Step 1: Q = sum_e Re(w_e (z_e I - A)^{-1} Y), upper-half. */
        for (size_t i = 0; i < n * m0; i++) mpfr_set_zero(Q[i], 1);
        for (int64_t e = 0; e < Ne; e++) {
            /* t = (pi/2)(1 + xi_e); z = c + r e^{i t}; w = (1/2) eta_e r e^{i t}. */
            mpfr_add_ui(t_var, gl_xi[e], 1, MPFR_RNDN);
            mpfr_mul(t_var, t_var, half_pi, MPFR_RNDN);
            mpfr_cos(ct, t_var, MPFR_RNDN);
            mpfr_sin(st, t_var, MPFR_RNDN);
            mpfr_mul(z_re, r_rad, ct, MPFR_RNDN);
            mpfr_add(z_re, z_re, c_ctr, MPFR_RNDN);
            mpfr_mul(z_im, r_rad, st, MPFR_RNDN);
            mpfr_mul(w_re, gl_w[e], r_rad, MPFR_RNDN);
            mpfr_div_2si(w_re, w_re, 1, MPFR_RNDN);
            mpfr_mul(w_im, w_re, st, MPFR_RNDN);
            mpfr_mul(w_re, w_re, ct, MPFR_RNDN);

            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                    mpfr_neg(M_re[i * n + j], A->re[i * n + j], MPFR_RNDN);
                    mpfr_set_zero(M_im[i * n + j], 1);
                }
                mpfr_add(M_re[i * n + i], M_re[i * n + i], z_re, MPFR_RNDN);
                mpfr_add(M_im[i * n + i], M_im[i * n + i], z_im, MPFR_RNDN);
            }
            if (feast_complex_lu_factor_M(M_re, M_im, n, bits, piv) != 0) {
                fail_reason = "LU breakdown at a quadrature node";
                goto cleanup;
            }
            for (size_t col = 0; col < m0; col++) {
                for (size_t i = 0; i < n; i++) {
                    mpfr_set(x_re[i], Y[i * m0 + col], MPFR_RNDN);
                    mpfr_set_zero(x_im[i], 1);
                }
                feast_complex_lu_solve_M(M_re, M_im, piv, n, bits, x_re, x_im);
                /* Q[:,col] += w_re x_re - w_im x_im. */
                for (size_t i = 0; i < n; i++) {
                    mpfr_mul(prod, w_re, x_re[i], MPFR_RNDN);
                    mpfr_add(Q[i * m0 + col], Q[i * m0 + col], prod, MPFR_RNDN);
                    mpfr_mul(prod, w_im, x_im[i], MPFR_RNDN);
                    mpfr_sub(Q[i * m0 + col], Q[i * m0 + col], prod, MPFR_RNDN);
                }
            }
        }

        /* Step 2: AQ = A Q, A_q = Q^T AQ, B_q = Q^T Q. */
        for (size_t i = 0; i < n; i++) {
            for (size_t col = 0; col < m0; col++) {
                mpfr_set_zero(sum, 1);
                for (size_t k2 = 0; k2 < n; k2++) {
                    mpfr_mul(prod, A->re[i * n + k2], Q[k2 * m0 + col], MPFR_RNDN);
                    mpfr_add(sum, sum, prod, MPFR_RNDN);
                }
                mpfr_set(AQ[i * m0 + col], sum, MPFR_RNDN);
            }
        }
        for (size_t i = 0; i < m0; i++) {
            for (size_t j = 0; j < m0; j++) {
                mpfr_set_zero(sum, 1);
                for (size_t k2 = 0; k2 < n; k2++) {
                    mpfr_mul(prod, Q[k2 * m0 + i], AQ[k2 * m0 + j], MPFR_RNDN);
                    mpfr_add(sum, sum, prod, MPFR_RNDN);
                }
                mpfr_set(Aq[i * m0 + j], sum, MPFR_RNDN);
            }
        }
        for (size_t i = 0; i < m0; i++) {
            for (size_t j = 0; j < m0; j++) {
                mpfr_set_zero(sum, 1);
                for (size_t k2 = 0; k2 < n; k2++) {
                    mpfr_mul(prod, Q[k2 * m0 + i], Q[k2 * m0 + j], MPFR_RNDN);
                    mpfr_add(sum, sum, prod, MPFR_RNDN);
                }
                mpfr_set(Bq[i * m0 + j], sum, MPFR_RNDN);
            }
        }
        /* Symmetrise. */
        for (size_t i = 0; i < m0; i++) {
            for (size_t j = i + 1; j < m0; j++) {
                mpfr_add(sum, Aq[i * m0 + j], Aq[j * m0 + i], MPFR_RNDN);
                mpfr_div_2si(sum, sum, 1, MPFR_RNDN);
                mpfr_set(Aq[i * m0 + j], sum, MPFR_RNDN);
                mpfr_set(Aq[j * m0 + i], sum, MPFR_RNDN);
                mpfr_add(sum, Bq[i * m0 + j], Bq[j * m0 + i], MPFR_RNDN);
                mpfr_div_2si(sum, sum, 1, MPFR_RNDN);
                mpfr_set(Bq[i * m0 + j], sum, MPFR_RNDN);
                mpfr_set(Bq[j * m0 + i], sum, MPFR_RNDN);
            }
        }

        /* Step 3: Cholesky B_q = L L^T,  H = L^{-1} A_q L^{-T}. */
        if (feast_cholesky_M(Bq, m0, bits, Lchol) != 0) {
            fail_reason = "B_q not positive definite "
                          "(subspace likely undersized or rank-deficient)";
            goto cleanup;
        }
        for (size_t i = 0; i < m0 * m0; i++) mpfr_set(H[i], Aq[i], MPFR_RNDN);
        feast_trsm_left_lower_M(Lchol, H, m0, m0, bits);
        /* Right-apply L^{-T}: transpose, left-solve, transpose back. */
        for (size_t i = 0; i < m0; i++)
            for (size_t j = 0; j < m0; j++)
                mpfr_set(HQ[i * m0 + j], H[j * m0 + i], MPFR_RNDN);
        feast_trsm_left_lower_M(Lchol, HQ, m0, m0, bits);
        for (size_t i = 0; i < m0; i++)
            for (size_t j = 0; j < m0; j++)
                mpfr_set(H[i * m0 + j], HQ[j * m0 + i], MPFR_RNDN);
        for (size_t i = 0; i < m0; i++) {
            for (size_t j = i + 1; j < m0; j++) {
                mpfr_add(sum, H[i * m0 + j], H[j * m0 + i], MPFR_RNDN);
                mpfr_div_2si(sum, sum, 1, MPFR_RNDN);
                mpfr_set(H[i * m0 + j], sum, MPFR_RNDN);
                mpfr_set(H[j * m0 + i], sum, MPFR_RNDN);
            }
        }

        /* Step 4: inner real-symmetric eigendecomposition at MPFR. */
        if (m0 == 1) {
            mpfr_set(diag_h[0], H[0], MPFR_RNDN);
            mpfr_set_ui(W[0], 1, MPFR_RNDN);
        } else {
            for (size_t i = 0; i < m0; i++) mpfr_set_zero(sub_h[i], 1);
            direct_tridiag_real_sym_M(H, m0, bits, diag_h, sub_h, W, true,
                                        u_h, p_h, q_h, qr_tmp);
            if (direct_symtridiag_qr_M(diag_h, sub_h, m0, bits,
                                          W, true, qr_tmp) != 0) {
                fail_reason = "inner symmetric tridiagonal QR did not converge";
                goto cleanup;
            }
        }
        /* V = L^{-T} W (in-place on W). */
        feast_trsm_left_lower_T_M(Lchol, W, m0, m0, bits);
        for (size_t k = 0; k < m0; k++) mpfr_set(mu[k], diag_h[k], MPFR_RNDN);

        /* Ritz vectors Xhat = Q V. */
        for (size_t i = 0; i < n; i++) {
            for (size_t col = 0; col < m0; col++) {
                mpfr_set_zero(sum, 1);
                for (size_t k2 = 0; k2 < m0; k2++) {
                    mpfr_mul(prod, Q[i * m0 + k2], W[k2 * m0 + col], MPFR_RNDN);
                    mpfr_add(sum, sum, prod, MPFR_RNDN);
                }
                mpfr_set(Xhat[i * m0 + col], sum, MPFR_RNDN);
            }
        }

        /* Step 5: filter to [a, b] and compute residuals. */
        size_t k_inside = 0;
        mpfr_t max_rel_res, res2, xn2, diff, rel, abs_mu, max_one;
        mpfr_init2(max_rel_res, bits); mpfr_set_zero(max_rel_res, 1);
        mpfr_init2(res2,        bits);
        mpfr_init2(xn2,         bits);
        mpfr_init2(diff,        bits);
        mpfr_init2(rel,         bits);
        mpfr_init2(abs_mu,      bits);
        mpfr_init2(max_one,     bits);
        mpfr_t lo_pad, hi_pad;
        mpfr_init2(lo_pad, bits); mpfr_init2(hi_pad, bits);
        mpfr_sub(lo_pad, a_lo, pad, MPFR_RNDN);
        mpfr_add(hi_pad, b_hi, pad, MPFR_RNDN);
        for (size_t k = 0; k < m0; k++) {
            if (mpfr_cmp(mu[k], lo_pad) < 0) continue;
            if (mpfr_cmp(mu[k], hi_pad) > 0) continue;
            mpfr_set_zero(res2, 1);
            mpfr_set_zero(xn2, 1);
            for (size_t i = 0; i < n; i++) {
                mpfr_set_zero(sum, 1);
                for (size_t j = 0; j < n; j++) {
                    mpfr_mul(prod, A->re[i * n + j], Xhat[j * m0 + k], MPFR_RNDN);
                    mpfr_add(sum, sum, prod, MPFR_RNDN);
                }
                mpfr_mul(prod, mu[k], Xhat[i * m0 + k], MPFR_RNDN);
                mpfr_sub(diff, sum, prod, MPFR_RNDN);
                mpfr_mul(prod, diff, diff, MPFR_RNDN);
                mpfr_add(res2, res2, prod, MPFR_RNDN);
                mpfr_mul(prod, Xhat[i * m0 + k], Xhat[i * m0 + k], MPFR_RNDN);
                mpfr_add(xn2, xn2, prod, MPFR_RNDN);
            }
            mpfr_sqrt(res2, res2, MPFR_RNDN);
            if (mpfr_zero_p(xn2)) mpfr_set_ui(sum, 1, MPFR_RNDN);
            else                   mpfr_sqrt(sum, xn2, MPFR_RNDN);
            mpfr_abs(abs_mu, mu[k], MPFR_RNDN);
            mpfr_set_ui(max_one, 1, MPFR_RNDN);
            if (mpfr_cmp(abs_mu, max_one) > 0) mpfr_set(max_one, abs_mu, MPFR_RNDN);
            mpfr_mul(sum2, sum, max_one, MPFR_RNDN);
            mpfr_div(rel, res2, sum2, MPFR_RNDN);
            if (mpfr_cmp(rel, max_rel_res) > 0)
                mpfr_set(max_rel_res, rel, MPFR_RNDN);
            keep[k_inside++] = k;
        }
        int converged_now = (k_inside > 0 && mpfr_cmp(max_rel_res, tol) < 0);
        mpfr_clear(max_rel_res); mpfr_clear(res2); mpfr_clear(xn2);
        mpfr_clear(diff); mpfr_clear(rel); mpfr_clear(abs_mu); mpfr_clear(max_one);
        mpfr_clear(lo_pad); mpfr_clear(hi_pad);

        if (converged_now) {
            converged = true;
            kept = k_inside;
            break;
        }
        /* Refresh Y from Xhat (column-normalised). */
        for (size_t col = 0; col < m0; col++) {
            mpfr_set_zero(sum, 1);
            for (size_t i = 0; i < n; i++) {
                mpfr_mul(prod, Xhat[i * m0 + col], Xhat[i * m0 + col], MPFR_RNDN);
                mpfr_add(sum, sum, prod, MPFR_RNDN);
            }
            if (mpfr_zero_p(sum)) {
                /* Reseed deterministically per column. */
                for (size_t i = 0; i < n; i++) {
                    uint64_t s = (uint64_t)(col + 1) * 2654435761ull
                                + (uint64_t)i * 40503ull;
                    s ^= s >> 33;
                    s *= 0xff51afd7ed558ccdull;
                    s ^= s >> 33;
                    double d = ((double)(uint32_t)(s >> 32)
                                  * (2.0 / 4294967295.0)) - 1.0;
                    mpfr_set_d(Y[i * m0 + col], d, MPFR_RNDN);
                }
            } else {
                mpfr_sqrt(sum, sum, MPFR_RNDN);
                for (size_t i = 0; i < n; i++) {
                    mpfr_div(Y[i * m0 + col], Xhat[i * m0 + col], sum, MPFR_RNDN);
                }
            }
        }
    }

    if (!converged) {
        if (!fail_reason) fail_reason = "MaxIterations exhausted";
        goto cleanup;
    }

    /* Step 6: emit result.  Pack and sort kept eigenpairs descending |mu|. */
    {
        mpfr_t* vals = mpfr_array_alloc(kept, bits);
        mpfr_t* vecs = want_Q ? mpfr_array_alloc(n * kept, bits) : NULL;
        for (size_t i = 0; i < kept; i++) {
            size_t src = keep[i];
            mpfr_set(vals[i], mu[src], MPFR_RNDN);
            if (want_Q) {
                mpfr_set_zero(sum, 1);
                for (size_t r2 = 0; r2 < n; r2++) {
                    mpfr_mul(prod, Xhat[r2 * m0 + src],
                                    Xhat[r2 * m0 + src], MPFR_RNDN);
                    mpfr_add(sum, sum, prod, MPFR_RNDN);
                }
                if (mpfr_zero_p(sum)) {
                    for (size_t r2 = 0; r2 < n; r2++)
                        mpfr_set(vecs[r2 * kept + i],
                                  Xhat[r2 * m0 + src], MPFR_RNDN);
                } else {
                    mpfr_sqrt(sum, sum, MPFR_RNDN);
                    for (size_t r2 = 0; r2 < n; r2++) {
                        mpfr_div(vecs[r2 * kept + i],
                                  Xhat[r2 * m0 + src], sum, MPFR_RNDN);
                    }
                }
            }
        }
        size_t* perm = (size_t*)malloc(sizeof(size_t) * kept);
        feast_sort_perm_desc_abs_M(vals, kept, perm);
        out = want_Q
            ? feast_build_real_eigenvector_list_subset_M(vecs, n, kept, perm)
            : feast_build_real_eigenvalue_list_subset_M(vals, kept, perm);
        free(perm);
        mpfr_array_free(vals, kept);
        if (vecs) mpfr_array_free(vecs, n * kept);
        out = direct_apply_k_spec_list(out, k_spec);
    }

cleanup:
    mpfr_array_free(Y, n * m0); mpfr_array_free(Q, n * m0);
    mpfr_array_free(AQ, n * m0); mpfr_array_free(Xhat, n * m0);
    mpfr_array_free(Aq, m0 * m0); mpfr_array_free(Bq, m0 * m0);
    mpfr_array_free(Lchol, m0 * m0); mpfr_array_free(H, m0 * m0);
    mpfr_array_free(HQ, m0 * m0); mpfr_array_free(W, m0 * m0);
    mpfr_array_free(diag_h, m0);
    mpfr_array_free(sub_h, m0 ? m0 : 1);
    mpfr_array_free(u_h, m0); mpfr_array_free(p_h, m0); mpfr_array_free(q_h, m0);
    mpfr_array_free(qr_tmp, 12);
    mpfr_array_free(M_re, n * n); mpfr_array_free(M_im, n * n);
    free(piv);
    mpfr_array_free(x_re, n); mpfr_array_free(x_im, n);
    mpfr_array_free(mu, m0);
    free(keep);
    mpfr_array_free(gl_xi, (size_t)Ne); mpfr_array_free(gl_w, (size_t)Ne);
    mpfr_clear(a_lo); mpfr_clear(b_hi); mpfr_clear(c_ctr); mpfr_clear(r_rad);
    mpfr_clear(tol);
    mpfr_clear(t_var); mpfr_clear(ct); mpfr_clear(st);
    mpfr_clear(z_re); mpfr_clear(z_im); mpfr_clear(w_re); mpfr_clear(w_im);
    mpfr_clear(prod); mpfr_clear(sum); mpfr_clear(sum2);
    mpfr_clear(pi_val); mpfr_clear(half_pi);
    mpfr_clear(pad); mpfr_clear(abs_a); mpfr_clear(abs_b);

    if (!converged && fail_reason) feast_warn_fallback(fail_reason);
    return out;
}

/* Complex-Hermitian FEAST kernel at MPFR precision.  Full-contour
 * quadrature (no Schwarz symmetry for complex A), complex subspace,
 * Hermitian Rayleigh-Ritz, inner direct_tridiag_complex_hermitian_M
 * pipeline. */
static Expr* feast_complex_hermitian_mpfr(const MatM* A, MateigenWant want,
                                            Expr* k_spec,
                                            const FeastOpts* opts) {
    size_t n = A->n;
    mpfr_prec_t bits = A->bits;
    if (n == 0) return NULL;
    if (!A->is_complex) return NULL;
    if (opts->interval_high <= opts->interval_low) return NULL;

    size_t m0;
    if (opts->subspace_size > 0) {
        m0 = (size_t)opts->subspace_size;
    } else {
        size_t auto_m = n / 4;
        if (auto_m < 20) auto_m = 20;
        m0 = auto_m;
    }
    if (m0 > n) m0 = n;
    if (m0 == 0) return NULL;

    bool   want_Q   = (want & MATEIGEN_WANT_VECTORS) != 0;
    int64_t max_iter = opts->max_iterations > 0 ? opts->max_iterations : 20;
    int64_t Ne = opts->contour_points > 0 ? opts->contour_points : 8;
    if (Ne != 2 && Ne != 4 && Ne != 8 && Ne != 16) Ne = 8;

    mpfr_t a_lo, b_hi, c_ctr, r_rad, tol;
    mpfr_init2(a_lo,  bits); mpfr_set_d(a_lo,  opts->interval_low,  MPFR_RNDN);
    mpfr_init2(b_hi,  bits); mpfr_set_d(b_hi,  opts->interval_high, MPFR_RNDN);
    mpfr_init2(c_ctr, bits); mpfr_init2(r_rad, bits);
    mpfr_add(c_ctr, a_lo, b_hi, MPFR_RNDN);
    mpfr_div_2si(c_ctr, c_ctr, 1, MPFR_RNDN);
    mpfr_sub(r_rad, b_hi, a_lo, MPFR_RNDN);
    mpfr_div_2si(r_rad, r_rad, 1, MPFR_RNDN);
    mpfr_init2(tol, bits);
    if (opts->tolerance_given && opts->tolerance > 0.0) {
        mpfr_set_d(tol, opts->tolerance, MPFR_RNDN);
    } else {
        feast_default_tolerance_M(bits, tol);
    }

    mpfr_t* gl_xi = mpfr_array_alloc((size_t)Ne, bits);
    mpfr_t* gl_w  = mpfr_array_alloc((size_t)Ne, bits);
    feast_gl_compute_mpfr(Ne, bits, gl_xi, gl_w);

    /* Workspace -- complex variants (re/im pairs). */
    mpfr_t* Y_re   = mpfr_array_alloc(n * m0, bits);
    mpfr_t* Y_im   = mpfr_array_alloc(n * m0, bits);
    mpfr_t* Q_re   = mpfr_array_alloc(n * m0, bits);
    mpfr_t* Q_im   = mpfr_array_alloc(n * m0, bits);
    mpfr_t* AQ_re  = mpfr_array_alloc(n * m0, bits);
    mpfr_t* AQ_im  = mpfr_array_alloc(n * m0, bits);
    mpfr_t* Xh_re  = mpfr_array_alloc(n * m0, bits);
    mpfr_t* Xh_im  = mpfr_array_alloc(n * m0, bits);
    mpfr_t* Aq_re  = mpfr_array_alloc(m0 * m0, bits);
    mpfr_t* Aq_im  = mpfr_array_alloc(m0 * m0, bits);
    mpfr_t* Bq_re  = mpfr_array_alloc(m0 * m0, bits);
    mpfr_t* Bq_im  = mpfr_array_alloc(m0 * m0, bits);
    mpfr_t* L_re   = mpfr_array_alloc(m0 * m0, bits);
    mpfr_t* L_im   = mpfr_array_alloc(m0 * m0, bits);
    mpfr_t* H_re   = mpfr_array_alloc(m0 * m0, bits);
    mpfr_t* H_im   = mpfr_array_alloc(m0 * m0, bits);
    mpfr_t* tmp_re = mpfr_array_alloc(m0 * m0, bits);
    mpfr_t* tmp_im = mpfr_array_alloc(m0 * m0, bits);
    mpfr_t* W_re   = mpfr_array_alloc(m0 * m0, bits);
    mpfr_t* W_im   = mpfr_array_alloc(m0 * m0, bits);
    mpfr_t* Z      = mpfr_array_alloc(m0 * m0, bits);
    mpfr_t* V_re   = mpfr_array_alloc(m0 * m0, bits);
    mpfr_t* V_im   = mpfr_array_alloc(m0 * m0, bits);
    mpfr_t* diag_h = mpfr_array_alloc(m0, bits);
    mpfr_t* sub_re = mpfr_array_alloc(m0 ? m0 : 1, bits);
    mpfr_t* sub_im = mpfr_array_alloc(m0 ? m0 : 1, bits);
    mpfr_t* u_re   = mpfr_array_alloc(m0, bits);
    mpfr_t* u_im   = mpfr_array_alloc(m0, bits);
    mpfr_t* v_re   = mpfr_array_alloc(m0, bits);
    mpfr_t* v_im   = mpfr_array_alloc(m0, bits);
    mpfr_t* qs_re  = mpfr_array_alloc(m0, bits);
    mpfr_t* qs_im  = mpfr_array_alloc(m0, bits);
    mpfr_t* qr_tmp = mpfr_array_alloc(12, bits);
    mpfr_t* M_re   = mpfr_array_alloc(n * n, bits);
    mpfr_t* M_im   = mpfr_array_alloc(n * n, bits);
    int*    piv    = (int*)malloc(sizeof(int) * n);
    mpfr_t* x_re   = mpfr_array_alloc(n, bits);
    mpfr_t* x_im   = mpfr_array_alloc(n, bits);
    mpfr_t* mu     = mpfr_array_alloc(m0, bits);
    size_t* keep   = (size_t*)malloc(sizeof(size_t) * m0);

    mpfr_t t_var, ct, st, z_re, z_im, w_re, w_im;
    mpfr_t prod, sum_r, sum_i, prod2;
    mpfr_init2(t_var, bits); mpfr_init2(ct,   bits); mpfr_init2(st,   bits);
    mpfr_init2(z_re,  bits); mpfr_init2(z_im, bits);
    mpfr_init2(w_re,  bits); mpfr_init2(w_im, bits);
    mpfr_init2(prod,  bits); mpfr_init2(sum_r, bits); mpfr_init2(sum_i, bits);
    mpfr_init2(prod2, bits);
    mpfr_t pi_val;
    mpfr_init2(pi_val, bits);
    mpfr_const_pi(pi_val, MPFR_RNDN);
    mpfr_t pad, abs_a, abs_b;
    mpfr_init2(pad,   bits);
    mpfr_init2(abs_a, bits);
    mpfr_init2(abs_b, bits);
    mpfr_abs(abs_a, a_lo, MPFR_RNDN);
    mpfr_abs(abs_b, b_hi, MPFR_RNDN);
    if (mpfr_cmp(abs_a, abs_b) >= 0) mpfr_mul(pad, tol, abs_a, MPFR_RNDN);
    else                              mpfr_mul(pad, tol, abs_b, MPFR_RNDN);
    if (mpfr_cmp(pad, tol) < 0) mpfr_set(pad, tol, MPFR_RNDN);

    feast_seed_subspace_M(Y_re, n, m0, /*seed=*/ 1442695040888963407ull);
    feast_seed_subspace_M(Y_im, n, m0, /*seed=*/ 6364136223846793005ull);

    Expr*  out      = NULL;
    bool   converged = false;
    size_t kept     = 0;
    const char* fail_reason = NULL;

    for (int64_t iter = 0; iter < max_iter && !converged; iter++) {
        /* Step 1: Q = sum_e w_e (z_e I - A)^{-1} Y (full contour). */
        for (size_t i = 0; i < n * m0; i++) {
            mpfr_set_zero(Q_re[i], 1);
            mpfr_set_zero(Q_im[i], 1);
        }
        for (int64_t e = 0; e < Ne; e++) {
            /* t = pi (1 + xi_e); z = c + r e^{i t}; w = (1/2) eta_e r e^{i t}. */
            mpfr_add_ui(t_var, gl_xi[e], 1, MPFR_RNDN);
            mpfr_mul(t_var, t_var, pi_val, MPFR_RNDN);
            mpfr_cos(ct, t_var, MPFR_RNDN);
            mpfr_sin(st, t_var, MPFR_RNDN);
            mpfr_mul(z_re, r_rad, ct, MPFR_RNDN);
            mpfr_add(z_re, z_re, c_ctr, MPFR_RNDN);
            mpfr_mul(z_im, r_rad, st, MPFR_RNDN);
            mpfr_mul(w_re, gl_w[e], r_rad, MPFR_RNDN);
            mpfr_div_2si(w_re, w_re, 1, MPFR_RNDN);
            mpfr_mul(w_im, w_re, st, MPFR_RNDN);
            mpfr_mul(w_re, w_re, ct, MPFR_RNDN);

            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                    mpfr_neg(M_re[i * n + j], A->re[i * n + j], MPFR_RNDN);
                    mpfr_neg(M_im[i * n + j], A->im[i * n + j], MPFR_RNDN);
                }
                mpfr_add(M_re[i * n + i], M_re[i * n + i], z_re, MPFR_RNDN);
                mpfr_add(M_im[i * n + i], M_im[i * n + i], z_im, MPFR_RNDN);
            }
            if (feast_complex_lu_factor_M(M_re, M_im, n, bits, piv) != 0) {
                fail_reason = "LU breakdown at a quadrature node";
                goto cleanup;
            }
            for (size_t col = 0; col < m0; col++) {
                for (size_t i = 0; i < n; i++) {
                    mpfr_set(x_re[i], Y_re[i * m0 + col], MPFR_RNDN);
                    mpfr_set(x_im[i], Y_im[i * m0 + col], MPFR_RNDN);
                }
                feast_complex_lu_solve_M(M_re, M_im, piv, n, bits, x_re, x_im);
                /* Q[:, col] += w_e * x_e (complex multiplication). */
                for (size_t i = 0; i < n; i++) {
                    /* xr = x_re[i], xi = x_im[i]. */
                    mpfr_mul(prod,  w_re, x_re[i], MPFR_RNDN);
                    mpfr_mul(prod2, w_im, x_im[i], MPFR_RNDN);
                    mpfr_sub(prod,  prod, prod2, MPFR_RNDN);
                    mpfr_add(Q_re[i * m0 + col], Q_re[i * m0 + col], prod, MPFR_RNDN);
                    mpfr_mul(prod,  w_re, x_im[i], MPFR_RNDN);
                    mpfr_mul(prod2, w_im, x_re[i], MPFR_RNDN);
                    mpfr_add(prod,  prod, prod2, MPFR_RNDN);
                    mpfr_add(Q_im[i * m0 + col], Q_im[i * m0 + col], prod, MPFR_RNDN);
                }
            }
        }

        /* Step 2: AQ = A Q (complex), A_q = Q^* AQ, B_q = Q^* Q. */
        for (size_t i = 0; i < n; i++) {
            for (size_t col = 0; col < m0; col++) {
                mpfr_set_zero(sum_r, 1); mpfr_set_zero(sum_i, 1);
                for (size_t k2 = 0; k2 < n; k2++) {
                    /* (ar + i ai)(qr + i qi). */
                    mpfr_mul(prod,  A->re[i * n + k2], Q_re[k2 * m0 + col], MPFR_RNDN);
                    mpfr_add(sum_r, sum_r, prod, MPFR_RNDN);
                    mpfr_mul(prod,  A->im[i * n + k2], Q_im[k2 * m0 + col], MPFR_RNDN);
                    mpfr_sub(sum_r, sum_r, prod, MPFR_RNDN);
                    mpfr_mul(prod,  A->re[i * n + k2], Q_im[k2 * m0 + col], MPFR_RNDN);
                    mpfr_add(sum_i, sum_i, prod, MPFR_RNDN);
                    mpfr_mul(prod,  A->im[i * n + k2], Q_re[k2 * m0 + col], MPFR_RNDN);
                    mpfr_add(sum_i, sum_i, prod, MPFR_RNDN);
                }
                mpfr_set(AQ_re[i * m0 + col], sum_r, MPFR_RNDN);
                mpfr_set(AQ_im[i * m0 + col], sum_i, MPFR_RNDN);
            }
        }
        for (size_t i = 0; i < m0; i++) {
            for (size_t j = 0; j < m0; j++) {
                mpfr_set_zero(sum_r, 1); mpfr_set_zero(sum_i, 1);
                for (size_t k2 = 0; k2 < n; k2++) {
                    /* conj(Q[k2,i]) * AQ[k2,j] = (qr - i qi)(ar + i ai). */
                    mpfr_mul(prod,  Q_re[k2 * m0 + i], AQ_re[k2 * m0 + j], MPFR_RNDN);
                    mpfr_add(sum_r, sum_r, prod, MPFR_RNDN);
                    mpfr_mul(prod,  Q_im[k2 * m0 + i], AQ_im[k2 * m0 + j], MPFR_RNDN);
                    mpfr_add(sum_r, sum_r, prod, MPFR_RNDN);
                    mpfr_mul(prod,  Q_re[k2 * m0 + i], AQ_im[k2 * m0 + j], MPFR_RNDN);
                    mpfr_add(sum_i, sum_i, prod, MPFR_RNDN);
                    mpfr_mul(prod,  Q_im[k2 * m0 + i], AQ_re[k2 * m0 + j], MPFR_RNDN);
                    mpfr_sub(sum_i, sum_i, prod, MPFR_RNDN);
                }
                mpfr_set(Aq_re[i * m0 + j], sum_r, MPFR_RNDN);
                mpfr_set(Aq_im[i * m0 + j], sum_i, MPFR_RNDN);
            }
        }
        for (size_t i = 0; i < m0; i++) {
            for (size_t j = 0; j < m0; j++) {
                mpfr_set_zero(sum_r, 1); mpfr_set_zero(sum_i, 1);
                for (size_t k2 = 0; k2 < n; k2++) {
                    /* (qir - i qii)(qjr + i qji). */
                    mpfr_mul(prod,  Q_re[k2 * m0 + i], Q_re[k2 * m0 + j], MPFR_RNDN);
                    mpfr_add(sum_r, sum_r, prod, MPFR_RNDN);
                    mpfr_mul(prod,  Q_im[k2 * m0 + i], Q_im[k2 * m0 + j], MPFR_RNDN);
                    mpfr_add(sum_r, sum_r, prod, MPFR_RNDN);
                    mpfr_mul(prod,  Q_re[k2 * m0 + i], Q_im[k2 * m0 + j], MPFR_RNDN);
                    mpfr_add(sum_i, sum_i, prod, MPFR_RNDN);
                    mpfr_mul(prod,  Q_im[k2 * m0 + i], Q_re[k2 * m0 + j], MPFR_RNDN);
                    mpfr_sub(sum_i, sum_i, prod, MPFR_RNDN);
                }
                mpfr_set(Bq_re[i * m0 + j], sum_r, MPFR_RNDN);
                mpfr_set(Bq_im[i * m0 + j], sum_i, MPFR_RNDN);
            }
        }
        /* Hermitianise A_q and B_q. */
        for (size_t i = 0; i < m0; i++) {
            mpfr_set_zero(Aq_im[i * m0 + i], 1);
            mpfr_set_zero(Bq_im[i * m0 + i], 1);
            for (size_t j = i + 1; j < m0; j++) {
                mpfr_add(sum_r, Aq_re[i * m0 + j], Aq_re[j * m0 + i], MPFR_RNDN);
                mpfr_div_2si(sum_r, sum_r, 1, MPFR_RNDN);
                mpfr_set(Aq_re[i * m0 + j], sum_r, MPFR_RNDN);
                mpfr_set(Aq_re[j * m0 + i], sum_r, MPFR_RNDN);
                mpfr_sub(sum_i, Aq_im[i * m0 + j], Aq_im[j * m0 + i], MPFR_RNDN);
                mpfr_div_2si(sum_i, sum_i, 1, MPFR_RNDN);
                mpfr_set(Aq_im[i * m0 + j], sum_i, MPFR_RNDN);
                mpfr_neg(Aq_im[j * m0 + i], sum_i, MPFR_RNDN);

                mpfr_add(sum_r, Bq_re[i * m0 + j], Bq_re[j * m0 + i], MPFR_RNDN);
                mpfr_div_2si(sum_r, sum_r, 1, MPFR_RNDN);
                mpfr_set(Bq_re[i * m0 + j], sum_r, MPFR_RNDN);
                mpfr_set(Bq_re[j * m0 + i], sum_r, MPFR_RNDN);
                mpfr_sub(sum_i, Bq_im[i * m0 + j], Bq_im[j * m0 + i], MPFR_RNDN);
                mpfr_div_2si(sum_i, sum_i, 1, MPFR_RNDN);
                mpfr_set(Bq_im[i * m0 + j], sum_i, MPFR_RNDN);
                mpfr_neg(Bq_im[j * m0 + i], sum_i, MPFR_RNDN);
            }
        }

        /* Step 3: Cholesky B_q = L L^*,  H = L^{-1} A_q L^{-*}. */
        if (feast_complex_cholesky_M(Bq_re, Bq_im, m0, bits, L_re, L_im) != 0) {
            fail_reason = "B_q not positive definite "
                          "(subspace likely undersized or rank-deficient)";
            goto cleanup;
        }
        for (size_t i = 0; i < m0 * m0; i++) {
            mpfr_set(H_re[i], Aq_re[i], MPFR_RNDN);
            mpfr_set(H_im[i], Aq_im[i], MPFR_RNDN);
        }
        feast_complex_trsm_left_lower_M(L_re, L_im, H_re, H_im, m0, m0, bits);
        /* Right-apply L^{-*}: conj-transpose, left-solve, conj-transpose back. */
        for (size_t i = 0; i < m0; i++) {
            for (size_t j = 0; j < m0; j++) {
                mpfr_set(tmp_re[i * m0 + j], H_re[j * m0 + i], MPFR_RNDN);
                mpfr_neg(tmp_im[i * m0 + j], H_im[j * m0 + i], MPFR_RNDN);
            }
        }
        feast_complex_trsm_left_lower_M(L_re, L_im, tmp_re, tmp_im, m0, m0, bits);
        for (size_t i = 0; i < m0; i++) {
            for (size_t j = 0; j < m0; j++) {
                mpfr_set(H_re[i * m0 + j], tmp_re[j * m0 + i], MPFR_RNDN);
                mpfr_neg(H_im[i * m0 + j], tmp_im[j * m0 + i], MPFR_RNDN);
            }
        }
        /* Hermitianise H. */
        for (size_t i = 0; i < m0; i++) {
            mpfr_set_zero(H_im[i * m0 + i], 1);
            for (size_t j = i + 1; j < m0; j++) {
                mpfr_add(sum_r, H_re[i * m0 + j], H_re[j * m0 + i], MPFR_RNDN);
                mpfr_div_2si(sum_r, sum_r, 1, MPFR_RNDN);
                mpfr_set(H_re[i * m0 + j], sum_r, MPFR_RNDN);
                mpfr_set(H_re[j * m0 + i], sum_r, MPFR_RNDN);
                mpfr_sub(sum_i, H_im[i * m0 + j], H_im[j * m0 + i], MPFR_RNDN);
                mpfr_div_2si(sum_i, sum_i, 1, MPFR_RNDN);
                mpfr_set(H_im[i * m0 + j], sum_i, MPFR_RNDN);
                mpfr_neg(H_im[j * m0 + i], sum_i, MPFR_RNDN);
            }
        }

        /* Step 4: inner complex Hermitian eigendecomposition at MPFR. */
        if (m0 == 1) {
            mpfr_set(diag_h[0], H_re[0], MPFR_RNDN);
            mpfr_set_ui(W_re[0], 1, MPFR_RNDN);
            mpfr_set_zero(W_im[0], 1);
            for (size_t i = 0; i < m0 * m0; i++) {
                mpfr_set_ui(Z[i], (i % (m0 + 1) == 0) ? 1 : 0, MPFR_RNDN);
            }
            mpfr_set(V_re[0], W_re[0], MPFR_RNDN);
            mpfr_set_zero(V_im[0], 1);
        } else {
            for (size_t i = 0; i < m0; i++) {
                mpfr_set_zero(sub_re[i], 1);
                mpfr_set_zero(sub_im[i], 1);
            }
            direct_tridiag_complex_hermitian_M(H_re, H_im, m0, bits,
                                                 diag_h, sub_re, sub_im,
                                                 W_re, W_im, true,
                                                 u_re, u_im, v_re, v_im,
                                                 qs_re, qs_im);
            direct_phase_correct_tridiag_M(sub_re, sub_im, m0, bits,
                                             W_re, W_im, true);
            for (size_t i = 0; i < m0; i++)
                for (size_t j = 0; j < m0; j++)
                    mpfr_set_si(Z[i * m0 + j], (i == j) ? 1 : 0, MPFR_RNDN);
            if (direct_symtridiag_qr_M(diag_h, sub_re, m0, bits,
                                          Z, true, qr_tmp) != 0) {
                fail_reason = "inner symmetric tridiagonal QR did not converge";
                goto cleanup;
            }
            direct_compose_complex_Q_real_Z_M(W_re, W_im, Z, m0, bits, V_re, V_im);
        }
        /* V = L^{-*} V. */
        feast_complex_trsm_left_lower_H_M(L_re, L_im, V_re, V_im, m0, m0, bits);
        for (size_t k = 0; k < m0; k++) mpfr_set(mu[k], diag_h[k], MPFR_RNDN);

        /* Ritz vectors Xhat = Q V (complex). */
        for (size_t i = 0; i < n; i++) {
            for (size_t col = 0; col < m0; col++) {
                mpfr_set_zero(sum_r, 1); mpfr_set_zero(sum_i, 1);
                for (size_t k2 = 0; k2 < m0; k2++) {
                    /* (qr + i qi)(vr + i vi). */
                    mpfr_mul(prod,  Q_re[i * m0 + k2], V_re[k2 * m0 + col], MPFR_RNDN);
                    mpfr_add(sum_r, sum_r, prod, MPFR_RNDN);
                    mpfr_mul(prod,  Q_im[i * m0 + k2], V_im[k2 * m0 + col], MPFR_RNDN);
                    mpfr_sub(sum_r, sum_r, prod, MPFR_RNDN);
                    mpfr_mul(prod,  Q_re[i * m0 + k2], V_im[k2 * m0 + col], MPFR_RNDN);
                    mpfr_add(sum_i, sum_i, prod, MPFR_RNDN);
                    mpfr_mul(prod,  Q_im[i * m0 + k2], V_re[k2 * m0 + col], MPFR_RNDN);
                    mpfr_add(sum_i, sum_i, prod, MPFR_RNDN);
                }
                mpfr_set(Xh_re[i * m0 + col], sum_r, MPFR_RNDN);
                mpfr_set(Xh_im[i * m0 + col], sum_i, MPFR_RNDN);
            }
        }

        /* Step 5: filter and residuals. */
        size_t k_inside = 0;
        mpfr_t max_rel_res, res2, xn2, dr, di, rel, abs_mu, max_one;
        mpfr_init2(max_rel_res, bits); mpfr_set_zero(max_rel_res, 1);
        mpfr_init2(res2,        bits);
        mpfr_init2(xn2,         bits);
        mpfr_init2(dr,          bits);
        mpfr_init2(di,          bits);
        mpfr_init2(rel,         bits);
        mpfr_init2(abs_mu,      bits);
        mpfr_init2(max_one,     bits);
        mpfr_t lo_pad, hi_pad;
        mpfr_init2(lo_pad, bits); mpfr_init2(hi_pad, bits);
        mpfr_sub(lo_pad, a_lo, pad, MPFR_RNDN);
        mpfr_add(hi_pad, b_hi, pad, MPFR_RNDN);
        for (size_t k = 0; k < m0; k++) {
            if (mpfr_cmp(mu[k], lo_pad) < 0) continue;
            if (mpfr_cmp(mu[k], hi_pad) > 0) continue;
            mpfr_set_zero(res2, 1);
            mpfr_set_zero(xn2,  1);
            for (size_t i = 0; i < n; i++) {
                mpfr_set_zero(sum_r, 1); mpfr_set_zero(sum_i, 1);
                for (size_t j = 0; j < n; j++) {
                    /* (ar + i ai)(xr + i xi). */
                    mpfr_mul(prod,  A->re[i * n + j], Xh_re[j * m0 + k], MPFR_RNDN);
                    mpfr_add(sum_r, sum_r, prod, MPFR_RNDN);
                    mpfr_mul(prod,  A->im[i * n + j], Xh_im[j * m0 + k], MPFR_RNDN);
                    mpfr_sub(sum_r, sum_r, prod, MPFR_RNDN);
                    mpfr_mul(prod,  A->re[i * n + j], Xh_im[j * m0 + k], MPFR_RNDN);
                    mpfr_add(sum_i, sum_i, prod, MPFR_RNDN);
                    mpfr_mul(prod,  A->im[i * n + j], Xh_re[j * m0 + k], MPFR_RNDN);
                    mpfr_add(sum_i, sum_i, prod, MPFR_RNDN);
                }
                mpfr_mul(prod, mu[k], Xh_re[i * m0 + k], MPFR_RNDN);
                mpfr_sub(dr,   sum_r, prod, MPFR_RNDN);
                mpfr_mul(prod, mu[k], Xh_im[i * m0 + k], MPFR_RNDN);
                mpfr_sub(di,   sum_i, prod, MPFR_RNDN);
                mpfr_mul(prod, dr, dr, MPFR_RNDN);
                mpfr_add(res2, res2, prod, MPFR_RNDN);
                mpfr_mul(prod, di, di, MPFR_RNDN);
                mpfr_add(res2, res2, prod, MPFR_RNDN);
                mpfr_mul(prod, Xh_re[i * m0 + k], Xh_re[i * m0 + k], MPFR_RNDN);
                mpfr_add(xn2, xn2, prod, MPFR_RNDN);
                mpfr_mul(prod, Xh_im[i * m0 + k], Xh_im[i * m0 + k], MPFR_RNDN);
                mpfr_add(xn2, xn2, prod, MPFR_RNDN);
            }
            mpfr_sqrt(res2, res2, MPFR_RNDN);
            if (mpfr_zero_p(xn2)) mpfr_set_ui(sum_r, 1, MPFR_RNDN);
            else                   mpfr_sqrt(sum_r, xn2, MPFR_RNDN);
            mpfr_abs(abs_mu, mu[k], MPFR_RNDN);
            mpfr_set_ui(max_one, 1, MPFR_RNDN);
            if (mpfr_cmp(abs_mu, max_one) > 0) mpfr_set(max_one, abs_mu, MPFR_RNDN);
            mpfr_mul(sum_i, sum_r, max_one, MPFR_RNDN);
            mpfr_div(rel, res2, sum_i, MPFR_RNDN);
            if (mpfr_cmp(rel, max_rel_res) > 0)
                mpfr_set(max_rel_res, rel, MPFR_RNDN);
            keep[k_inside++] = k;
        }
        int converged_now = (k_inside > 0 && mpfr_cmp(max_rel_res, tol) < 0);
        mpfr_clear(max_rel_res); mpfr_clear(res2); mpfr_clear(xn2);
        mpfr_clear(dr); mpfr_clear(di); mpfr_clear(rel);
        mpfr_clear(abs_mu); mpfr_clear(max_one);
        mpfr_clear(lo_pad); mpfr_clear(hi_pad);

        if (converged_now) {
            converged = true;
            kept = k_inside;
            break;
        }
        /* Refresh Y from Xhat. */
        for (size_t col = 0; col < m0; col++) {
            mpfr_set_zero(sum_r, 1);
            for (size_t i = 0; i < n; i++) {
                mpfr_mul(prod, Xh_re[i * m0 + col], Xh_re[i * m0 + col], MPFR_RNDN);
                mpfr_add(sum_r, sum_r, prod, MPFR_RNDN);
                mpfr_mul(prod, Xh_im[i * m0 + col], Xh_im[i * m0 + col], MPFR_RNDN);
                mpfr_add(sum_r, sum_r, prod, MPFR_RNDN);
            }
            if (mpfr_zero_p(sum_r)) {
                for (size_t i = 0; i < n; i++) {
                    uint64_t s = (uint64_t)(col + 1) * 2654435761ull
                                + (uint64_t)i * 40503ull;
                    s ^= s >> 33;
                    s *= 0xff51afd7ed558ccdull;
                    s ^= s >> 33;
                    double d = ((double)(uint32_t)(s >> 32)
                                  * (2.0 / 4294967295.0)) - 1.0;
                    mpfr_set_d(Y_re[i * m0 + col], d, MPFR_RNDN);
                    s = (uint64_t)(col + 7) * 0x9E3779B97F4A7C15ull
                         + (uint64_t)i * 0xBF58476D1CE4E5B9ull;
                    s ^= s >> 27;
                    d = ((double)(uint32_t)(s >> 32)
                            * (2.0 / 4294967295.0)) - 1.0;
                    mpfr_set_d(Y_im[i * m0 + col], d, MPFR_RNDN);
                }
            } else {
                mpfr_sqrt(sum_r, sum_r, MPFR_RNDN);
                for (size_t i = 0; i < n; i++) {
                    mpfr_div(Y_re[i * m0 + col],
                              Xh_re[i * m0 + col], sum_r, MPFR_RNDN);
                    mpfr_div(Y_im[i * m0 + col],
                              Xh_im[i * m0 + col], sum_r, MPFR_RNDN);
                }
            }
        }
    }

    if (!converged) {
        if (!fail_reason) fail_reason = "MaxIterations exhausted";
        goto cleanup;
    }

    /* Step 6: emit complex eigenpair list. */
    {
        mpfr_t* vals    = mpfr_array_alloc(kept, bits);
        mpfr_t* vecs_re = want_Q ? mpfr_array_alloc(n * kept, bits) : NULL;
        mpfr_t* vecs_im = want_Q ? mpfr_array_alloc(n * kept, bits) : NULL;
        for (size_t i = 0; i < kept; i++) {
            size_t src = keep[i];
            mpfr_set(vals[i], mu[src], MPFR_RNDN);
            if (want_Q) {
                mpfr_set_zero(sum_r, 1);
                for (size_t r2 = 0; r2 < n; r2++) {
                    mpfr_mul(prod, Xh_re[r2 * m0 + src],
                                    Xh_re[r2 * m0 + src], MPFR_RNDN);
                    mpfr_add(sum_r, sum_r, prod, MPFR_RNDN);
                    mpfr_mul(prod, Xh_im[r2 * m0 + src],
                                    Xh_im[r2 * m0 + src], MPFR_RNDN);
                    mpfr_add(sum_r, sum_r, prod, MPFR_RNDN);
                }
                if (mpfr_zero_p(sum_r)) {
                    for (size_t r2 = 0; r2 < n; r2++) {
                        mpfr_set(vecs_re[r2 * kept + i],
                                  Xh_re[r2 * m0 + src], MPFR_RNDN);
                        mpfr_set(vecs_im[r2 * kept + i],
                                  Xh_im[r2 * m0 + src], MPFR_RNDN);
                    }
                } else {
                    mpfr_sqrt(sum_r, sum_r, MPFR_RNDN);
                    for (size_t r2 = 0; r2 < n; r2++) {
                        mpfr_div(vecs_re[r2 * kept + i],
                                  Xh_re[r2 * m0 + src], sum_r, MPFR_RNDN);
                        mpfr_div(vecs_im[r2 * kept + i],
                                  Xh_im[r2 * m0 + src], sum_r, MPFR_RNDN);
                    }
                }
            }
        }
        size_t* perm = (size_t*)malloc(sizeof(size_t) * kept);
        feast_sort_perm_desc_abs_M(vals, kept, perm);
        out = want_Q
            ? feast_build_complex_eigenvector_list_subset_M(vecs_re, vecs_im,
                                                              n, kept, perm)
            : feast_build_real_eigenvalue_list_subset_M(vals, kept, perm);
        free(perm);
        mpfr_array_free(vals, kept);
        if (vecs_re) mpfr_array_free(vecs_re, n * kept);
        if (vecs_im) mpfr_array_free(vecs_im, n * kept);
        out = direct_apply_k_spec_list(out, k_spec);
    }

cleanup:
    mpfr_array_free(Y_re, n * m0); mpfr_array_free(Y_im, n * m0);
    mpfr_array_free(Q_re, n * m0); mpfr_array_free(Q_im, n * m0);
    mpfr_array_free(AQ_re, n * m0); mpfr_array_free(AQ_im, n * m0);
    mpfr_array_free(Xh_re, n * m0); mpfr_array_free(Xh_im, n * m0);
    mpfr_array_free(Aq_re, m0 * m0); mpfr_array_free(Aq_im, m0 * m0);
    mpfr_array_free(Bq_re, m0 * m0); mpfr_array_free(Bq_im, m0 * m0);
    mpfr_array_free(L_re, m0 * m0);  mpfr_array_free(L_im, m0 * m0);
    mpfr_array_free(H_re, m0 * m0);  mpfr_array_free(H_im, m0 * m0);
    mpfr_array_free(tmp_re, m0 * m0); mpfr_array_free(tmp_im, m0 * m0);
    mpfr_array_free(W_re, m0 * m0);  mpfr_array_free(W_im, m0 * m0);
    mpfr_array_free(Z, m0 * m0);
    mpfr_array_free(V_re, m0 * m0);  mpfr_array_free(V_im, m0 * m0);
    mpfr_array_free(diag_h, m0);
    mpfr_array_free(sub_re, m0 ? m0 : 1);
    mpfr_array_free(sub_im, m0 ? m0 : 1);
    mpfr_array_free(u_re, m0); mpfr_array_free(u_im, m0);
    mpfr_array_free(v_re, m0); mpfr_array_free(v_im, m0);
    mpfr_array_free(qs_re, m0); mpfr_array_free(qs_im, m0);
    mpfr_array_free(qr_tmp, 12);
    mpfr_array_free(M_re, n * n); mpfr_array_free(M_im, n * n);
    free(piv);
    mpfr_array_free(x_re, n); mpfr_array_free(x_im, n);
    mpfr_array_free(mu, m0);
    free(keep);
    mpfr_array_free(gl_xi, (size_t)Ne); mpfr_array_free(gl_w, (size_t)Ne);
    mpfr_clear(a_lo); mpfr_clear(b_hi); mpfr_clear(c_ctr); mpfr_clear(r_rad);
    mpfr_clear(tol);
    mpfr_clear(t_var); mpfr_clear(ct); mpfr_clear(st);
    mpfr_clear(z_re); mpfr_clear(z_im); mpfr_clear(w_re); mpfr_clear(w_im);
    mpfr_clear(prod); mpfr_clear(sum_r); mpfr_clear(sum_i); mpfr_clear(prod2);
    mpfr_clear(pi_val);
    mpfr_clear(pad); mpfr_clear(abs_a); mpfr_clear(abs_b);

    if (!converged && fail_reason) feast_warn_fallback(fail_reason);
    return out;
}

/* MPFR FEAST dispatcher.  Loads A into a MatM at working precision,
 * checks Hermitian shape, and routes to the right kernel.  Returns NULL
 * for any unsupported shape so the caller cascades to Direct (MPFR or
 * machine) -- see feast_dispatch. */
static Expr* feast_dispatch_mpfr(Expr* m, Expr* a, int64_t n,
                                   mpfr_prec_t bits,
                                   MateigenWant want, Expr* k_spec,
                                   const FeastOpts* opts) {
    (void)a;
    MatM A;
    if (!matM_load(m, (size_t)n, bits, &A)) return NULL;

    /* Hermitian tolerance: scale with matrix norm at working precision.
     * The factor mirrors feast_dispatch_machine's 1e-10 reference; at
     * MPFR we tighten to 2^(-bits + 16) so genuinely Hermitian inputs
     * with rounding noise are still recognised but a clearly non-
     * Hermitian matrix still cascades cleanly. */
    mpfr_t norm, herm_tol, scale;
    mpfr_init2(norm,     bits);
    mpfr_init2(herm_tol, bits);
    mpfr_init2(scale,    bits);
    if (A.is_complex) matM_norm_inf_complex(&A, norm);
    else              matM_norm_inf_real(A.re, A.n, bits, norm);
    mpfr_add_ui(scale, norm, 1, MPFR_RNDN);
    mpfr_set_ui(herm_tol, 1, MPFR_RNDN);
    long shift = (long)bits - 16;
    if (shift < 1) shift = 1;
    mpfr_div_2si(herm_tol, herm_tol, shift, MPFR_RNDN);
    mpfr_mul(herm_tol, herm_tol, scale, MPFR_RNDN);

    Expr* out = NULL;
    if (!A.is_complex) {
        if (matM_is_real_symmetric(&A, herm_tol)) {
            out = feast_real_sym_mpfr(&A, want, k_spec, opts);
        }
    } else if (matM_is_hermitian(&A, herm_tol)) {
        out = feast_complex_hermitian_mpfr(&A, want, k_spec, opts);
    }

    mpfr_clear(norm); mpfr_clear(herm_tol); mpfr_clear(scale);
    matM_free(&A);
    return out;
}
#endif /* USE_MPFR */

/* Top-level "FEAST" dispatcher.  Returns NULL when:
 *   - generalised eigenproblem (a != NULL): FEAST does not support
 *     generalised systems in Mathilda yet,
 *   - the user-supplied Interval is missing or malformed,
 *   - the requested precision kernel hasn't landed yet (Phase 1).
 * In every NULL case the caller falls through to direct_dispatch, so
 * the user still gets a result, just from the all-eigenvalues path. */
Expr* feast_dispatch(Expr* m, Expr* a, int64_t n,
                              MateigenWant want, Expr* k_spec,
                              Expr* method_value) {
    if (a != NULL) return NULL;
    if (n <= 0)    return NULL;

    FeastOpts opts;
    feast_parse_subopts(method_value, &opts);
    if (!opts.interval_given) return NULL;

    CommonInexactInfo info = common_scan_inexact(m);
    if (info.has_inexact && info.min_bits > 53) {
        Expr* out = feast_dispatch_mpfr(m, a, n,
                                          (mpfr_prec_t)info.min_bits,
                                          want, k_spec, &opts);
        if (out) return out;
    }
    return feast_dispatch_machine(m, a, n, want, k_spec, &opts);
}
