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


/* LAPACK-HOOK: replace with dsytrd when USE_LAPACK is set.
 *
 * Householder tridiagonalisation of a real symmetric n x n matrix `A`
 * (row-major, modified in place).  On return:
 *   diag[i]     = T_ii      for 0 <= i < n
 *   subdiag[i]  = T_{i+1,i} for 0 <= i < n-1
 *   Q (n*n)     = orthogonal accumulated reflectors, when `want_Q`.
 *
 * Algorithm: classic Householder reflectors applied symmetrically
 * via the rank-2 update A <- A - u q^T - q u^T where p = A u, K = u^T p / 2,
 * q = p - K u.  O(2 n^3 / 3) flops.  See Golub & Van Loan, Alg 8.3.1.
 *
 * Scratch buffers `u`, `p`, `q` (size n each) are caller-provided so
 * the inner loops never touch malloc/free.
 */
void direct_tridiag_real_sym(double* A, size_t n,
                                     double* diag, double* subdiag,
                                     double* Q, bool want_Q,
                                     double* u, double* p, double* q) {
    if (want_Q) {
        /* Q starts at identity; reflectors are applied from the right
         * so the columns of Q become the orthogonal eigenvectors of A
         * once symmetric QR finishes. */
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) Q[i * n + j] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (size_t k = 0; k + 2 < n; k++) {
        /* Compute Householder vector u for A[k+1:n, k]. */
        double sigma = 0.0;
        for (size_t i = k + 1; i < n; i++) {
            double v = A[i * n + k];
            sigma += v * v;
        }
        double xk1 = A[(k + 1) * n + k];
        if (sigma == 0.0) {
            subdiag[k] = xk1;
            continue;
        }
        /* sigma already contains xk1*xk1 (the loop runs i = k+1..n-1, and
         * xk1 = A[k+1, k] is the first such entry).  The cancellation-safe
         * form: alpha = -sign(xk1) * ||x||. */
        double norm_x = sqrt(sigma);
        double alpha = (xk1 >= 0.0) ? -norm_x : norm_x;

        /* u[k+1] = xk1 - alpha; u[i>k+1] = A[i, k];  then normalise. */
        u[k + 1] = xk1 - alpha;
        for (size_t i = k + 2; i < n; i++) u[i] = A[i * n + k];
        double unorm2 = u[k + 1] * u[k + 1];
        for (size_t i = k + 2; i < n; i++) unorm2 += u[i] * u[i];
        if (unorm2 == 0.0) {
            subdiag[k] = xk1;
            continue;
        }
        double unorm = sqrt(unorm2);
        for (size_t i = k + 1; i < n; i++) u[i] /= unorm;

        /* p = A_22 u   (only on the trailing sub-block) */
        for (size_t i = k + 1; i < n; i++) {
            double s = 0.0;
            for (size_t j = k + 1; j < n; j++) s += A[i * n + j] * u[j];
            p[i] = 2.0 * s;
        }
        /* K = u^T A u = (u^T p) / 2 (since p = 2 A u). */
        double K = 0.0;
        for (size_t i = k + 1; i < n; i++) K += u[i] * p[i];
        K *= 0.5;
        /* q = p - 2 K u.  The rank-2 update A <- A - u q^T - q u^T then
         * equals H A H for H = I - 2 u u^T (Golub & Van Loan Alg 8.3.1). */
        for (size_t i = k + 1; i < n; i++) q[i] = p[i] - 2.0 * K * u[i];

        /* A_22 -= u q^T + q u^T */
        for (size_t i = k + 1; i < n; i++) {
            for (size_t j = k + 1; j < n; j++) {
                A[i * n + j] -= u[i] * q[j] + q[i] * u[j];
            }
        }

        /* Set the new subdiagonal element and clear the eliminated
         * column / row entries explicitly (drift control). */
        subdiag[k] = alpha;
        A[(k + 1) * n + k] = alpha;
        A[k * n + (k + 1)] = alpha;
        for (size_t i = k + 2; i < n; i++) {
            A[i * n + k] = 0.0;
            A[k * n + i] = 0.0;
        }

        /* Q <- Q * H_k   (right multiplication; column i replaced by
         * Q_i - 2 (Q row . u) u_i).  Apply only to Q[*, k+1..n-1]. */
        if (want_Q) {
            for (size_t i = 0; i < n; i++) {
                double s = 0.0;
                for (size_t j = k + 1; j < n; j++) s += Q[i * n + j] * u[j];
                s *= 2.0;
                for (size_t j = k + 1; j < n; j++) Q[i * n + j] -= s * u[j];
            }
        }
    }

    /* Extract diagonal. */
    for (size_t i = 0; i < n; i++) diag[i] = A[i * n + i];
    /* The n-2 -> n-1 subdiagonal is already in A. */
    if (n >= 2) subdiag[n - 2] = A[(n - 1) * n + (n - 2)];
}

/* Implicit-shift symmetric tridiagonal QR with Wilkinson shift.
 *
 * LAPACK-HOOK: replace with dsteqr (or dstedc) when USE_LAPACK is set.
 *
 *   diag[0..n-1]    : in-place diagonal of the tridiagonal matrix
 *   sub [0..n-2]    : in-place sub/super-diagonal
 *   Q   (n*n)       : in-place orthogonal eigenvector accumulator,
 *                     when `want_Q`.  Caller initialises Q (typically
 *                     to the orthogonal matrix from the tridiag step).
 *
 * Iterates over the active sub-block, deflating when |sub[i]| falls
 * below |diag[i]|+|diag[i+1]| * relative_tol.  Returns 0 on success,
 * -1 if the maximum number of sweeps (30*n) is exceeded.  Stagnation
 * is exceptionally rare for symmetric tridiagonal inputs but we cap
 * to avoid theoretical hangs.
 */
int direct_symtridiag_qr(double* diag, double* sub, size_t n,
                                 double* Q, bool want_Q) {
    if (n == 0) return 0;
    const double rel_tol = 1e-14;   /* much tighter than chop threshold */
    const size_t max_sweeps = 30 * n;
    size_t sweeps = 0;

    size_t end = n;  /* active sub-block is [0..end-1]. */
    while (end > 1) {
        /* Find the largest m such that sub[m..end-2] are all "significant". */
        size_t m = end - 1;
        while (m > 0) {
            double tol = rel_tol * (fabs(diag[m - 1]) + fabs(diag[m]));
            if (fabs(sub[m - 1]) <= tol) { sub[m - 1] = 0.0; break; }
            m--;
        }
        if (m == end - 1) { end--; continue; }  /* deflated bottom */

        if (++sweeps > max_sweeps) return -1;

        /* Wilkinson shift on the trailing 2x2 block. */
        double d = (diag[end - 2] - diag[end - 1]) * 0.5;
        double e = sub[end - 2];
        double t = (d == 0.0) ? fabs(e)
                              : fabs(d) + sqrt(d * d + e * e);
        double sign_d = (d >= 0.0) ? 1.0 : -1.0;
        double mu = diag[end - 1] - sign_d * (e * e) / t;

        /* Implicit QR sweep on [m..end-1] using Givens rotations. */
        double x = diag[m] - mu;
        double z = sub[m];
        for (size_t k = m; k < end - 1; k++) {
            double c, s;
            double r = hypot(x, z);
            if (r == 0.0) { c = 1.0; s = 0.0; }
            else { c = x / r; s = z / r; }

            if (k > m) sub[k - 1] = r;

            double d_k    = diag[k];
            double d_k1   = diag[k + 1];
            double e_k    = sub[k];

            /* Two-sided Givens rotation Q^T A Q with Q = [c -s; s c]
             * (the rotation whose transpose annihilates z in [x; z]):
             *   d_k'   = c^2 d_k + 2 c s e_k + s^2 d_k1
             *   d_k1'  = s^2 d_k - 2 c s e_k + c^2 d_k1
             *   e_k'   = c s (d_k1 - d_k) + (c^2 - s^2) e_k
             */
            diag[k]     = c * c * d_k + 2.0 * c * s * e_k + s * s * d_k1;
            diag[k + 1] = s * s * d_k - 2.0 * c * s * e_k + c * c * d_k1;
            sub[k]      = c * s * (d_k1 - d_k) + (c * c - s * s) * e_k;

            /* Chase the bulge: rotation also affects sub[k+1] for k < end-2. */
            if (k + 1 < end - 1) {
                double t_next = sub[k + 1];
                x = sub[k];
                z = s * t_next;
                sub[k + 1] = c * t_next;
            }

            /* Update Q: post-multiply by the Givens rotation in columns
             * k and k+1.   Q_col_k  <- c Q_col_k + s Q_col_{k+1}
             *               Q_col_k1 <- -s Q_col_k + c Q_col_{k+1}    */
            if (want_Q) {
                for (size_t i = 0; i < n; i++) {
                    double qk  = Q[i * n + k];
                    double qk1 = Q[i * n + (k + 1)];
                    Q[i * n + k]       =  c * qk + s * qk1;
                    Q[i * n + (k + 1)] = -s * qk + c * qk1;
                }
            }
        }
    }
    return 0;
}

/* Sort eigenvalue indices into descending |lambda|, stable on ties.
 * Writes the permutation into `perm[0..n-1]`. */
void direct_sort_perm_desc_abs(const double* vals, size_t n,
                                       size_t* perm) {
    for (size_t i = 0; i < n; i++) perm[i] = i;
    /* Insertion sort: n is small for the matrix sizes we care about
     * and the comparison is cheap. */
    for (size_t i = 1; i < n; i++) {
        size_t cur = perm[i];
        double ac = fabs(vals[cur]);
        size_t j = i;
        while (j > 0) {
            double ap = fabs(vals[perm[j - 1]]);
            if (ap > ac || (ap == ac && perm[j - 1] < cur)) break;
            perm[j] = perm[j - 1];
            j--;
        }
        perm[j] = cur;
    }
}

/* Build the final List[Real,...] of eigenvalues from a sorted permutation. */
Expr* direct_build_real_eigenvalue_list(const double* vals, size_t n,
                                                const size_t* perm) {
    Expr** items = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        items[i] = expr_new_real(vals[perm[i]]);
    }
    Expr* out = expr_new_function(expr_new_symbol("List"), items, n);
    free(items);
    return out;
}

/* Build the final List[List[Real,...], ...] of eigenvectors from a sorted
 * permutation.  Q is the n x n column-major-of-eigenvectors matrix in
 * row-major storage: Q[i, p] is the i-th component of the p-th
 * eigenvector. */
Expr* direct_build_real_eigenvector_list(const double* Q, size_t n,
                                                 const size_t* perm) {
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t k = 0; k < n; k++) {
        size_t col = perm[k];
        Expr** comps = (Expr**)malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            comps[i] = expr_new_real(Q[i * n + col]);
        }
        rows[k] = expr_new_function(expr_new_symbol("List"), comps, n);
        free(comps);
    }
    Expr* out = expr_new_function(expr_new_symbol("List"), rows, n);
    free(rows);
    return out;
}

/* Apply k-spec selection to an already-sorted (descending |lambda|)
 * result.  Returns a fresh List with the trimmed entries.  Mirrors
 * eigen_apply_k_spec but for Expr trees we already own. */
Expr* direct_apply_k_spec_list(Expr* full_list, Expr* k_spec) {
    if (!k_spec) return full_list;
    size_t count = full_list->data.function.arg_count;
    size_t result_count = count;
    bool from_end = false;
    if (k_spec->type == EXPR_INTEGER) {
        int64_t k = k_spec->data.integer;
        if (k >= 0) {
            result_count = ((size_t)k < count) ? (size_t)k : count;
        } else {
            int64_t abs_k = -k;
            result_count = ((size_t)abs_k < count) ? (size_t)abs_k : count;
            from_end = true;
        }
    } else if (k_spec->type == EXPR_FUNCTION
        && k_spec->data.function.head->type == EXPR_SYMBOL
        && k_spec->data.function.head->data.symbol == SYM_UpTo
        && k_spec->data.function.arg_count == 1
        && k_spec->data.function.args[0]->type == EXPR_INTEGER) {
        int64_t k = k_spec->data.function.args[0]->data.integer;
        result_count = ((size_t)k < count) ? (size_t)k : count;
    }

    Expr** items = result_count
        ? (Expr**)malloc(sizeof(Expr*) * result_count) : NULL;
    if (from_end) {
        size_t start = count - result_count;
        for (size_t i = 0; i < start; i++)
            expr_free(full_list->data.function.args[i]);
        for (size_t i = 0; i < result_count; i++)
            items[i] = full_list->data.function.args[start + i];
    } else {
        for (size_t i = 0; i < result_count; i++)
            items[i] = full_list->data.function.args[i];
        for (size_t i = result_count; i < count; i++)
            expr_free(full_list->data.function.args[i]);
    }
    free(full_list->data.function.args);
    full_list->data.function.args = items;
    full_list->data.function.arg_count = result_count;
    return full_list;
}

/* ---------- Real non-symmetric Direct: Hessenberg + Francis QR ------ *
 *                                                                       *
 *  For a general real n x n matrix, eigenvalues are extracted by:       *
 *    1. Householder reduction to upper Hessenberg form.                 *
 *    2. Implicit double-shift QR (Francis) sweeps over the Hessenberg,  *
 *       deflating off-diagonal entries until the matrix is quasi-       *
 *       triangular (Schur form: 1x1 diagonal blocks for real            *
 *       eigenvalues, 2x2 blocks for complex conjugate pairs).           *
 *    3. Read eigenvalues off the Schur form.                            *
 *                                                                       *
 *  Eigenvectors require Q accumulation through both reduction and QR    *
 *  sweeps and arrive in step 2b (this commit covers eigenvalues only).  *
 *                                                                       *
 *  LAPACK-HOOK: dgehrd for step 1, dhseqr for step 2, dtrevc3 +         *
 *  dorghr for the eigenvector back-transform in step 2b.                *
 * --------------------------------------------------------------------- */

/* Householder reduction of an n x n matrix `A` (row-major, modified in
 * place) to upper Hessenberg form.  When `Q` is non-NULL the caller
 * passes a pre-initialised n x n identity matrix that this routine
 * post-multiplies by each Householder reflector, producing the
 * orthogonal back-transformation Q with Q^T A_in Q = A_out + (Schur).
 *
 * Scratch buffer `u` (size n) is caller-provided.
 *
 * LAPACK-HOOK: dgehrd (+ dorghr to materialise Q from the reflectors). */
void direct_hessenberg_real(double* A, size_t n, double* u, double* Q) {
    for (size_t k = 0; k + 2 < n; k++) {
        /* Householder vector for column k below the sub-diagonal. */
        double sigma = 0.0;
        for (size_t i = k + 1; i < n; i++) {
            double v = A[i * n + k];
            sigma += v * v;
        }
        if (sigma == 0.0) continue;
        double xk1 = A[(k + 1) * n + k];
        double norm_x = sqrt(sigma);
        double alpha = (xk1 >= 0.0) ? -norm_x : norm_x;
        u[k + 1] = xk1 - alpha;
        for (size_t i = k + 2; i < n; i++) u[i] = A[i * n + k];
        double unorm2 = u[k + 1] * u[k + 1];
        for (size_t i = k + 2; i < n; i++) unorm2 += u[i] * u[i];
        if (unorm2 == 0.0) continue;
        double unorm = sqrt(unorm2);
        for (size_t i = k + 1; i < n; i++) u[i] /= unorm;

        /* Left multiply: A[k+1..n-1, :] <- H * A[k+1..n-1, :] */
        for (size_t j = 0; j < n; j++) {
            double s = 0.0;
            for (size_t i = k + 1; i < n; i++) s += u[i] * A[i * n + j];
            s *= 2.0;
            for (size_t i = k + 1; i < n; i++) A[i * n + j] -= s * u[i];
        }
        /* Right multiply: A[:, k+1..n-1] <- A[:, k+1..n-1] * H */
        for (size_t i = 0; i < n; i++) {
            double s = 0.0;
            for (size_t j = k + 1; j < n; j++) s += A[i * n + j] * u[j];
            s *= 2.0;
            for (size_t j = k + 1; j < n; j++) A[i * n + j] -= s * u[j];
        }

        /* Q <- Q * H_k (right-multiply by reflector on cols k+1..n-1). */
        if (Q) {
            for (size_t i = 0; i < n; i++) {
                double s = 0.0;
                for (size_t j = k + 1; j < n; j++) s += Q[i * n + j] * u[j];
                s *= 2.0;
                for (size_t j = k + 1; j < n; j++) Q[i * n + j] -= s * u[j];
            }
        }

        /* Tidy up: explicit subdiag and zero below it. */
        A[(k + 1) * n + k] = alpha;
        for (size_t i = k + 2; i < n; i++) A[i * n + k] = 0.0;
    }
}

/* One Francis double-shift QR step applied to the active sub-block
 * H[q..p-1, q..p-1] of an upper Hessenberg matrix.  The shift pair is
 * derived from the eigenvalues of the trailing 2x2 block H[p-2:p-1,
 * p-2:p-1] (real or complex; real arithmetic throughout since the
 * shift pair is a conjugate pair when complex).
 *
 * Bulge chasing uses Householder reflectors of size 3 (size 2 at the
 * very end of the chase).  Sub-block size p - q must be >= 3.
 *
 * When `Q` is non-NULL each bulge-chase reflector also right-multiplies
 * Q so the caller can extract eigenvectors from the final Schur form.
 *
 * LAPACK-HOOK: this is the inner loop body of dhseqr / dlahqr. */
static void direct_francis_step(double* H, size_t n, size_t q, size_t p,
                                 double* Q) {
    double h11 = H[(p - 2) * n + (p - 2)];
    double h12 = H[(p - 2) * n + (p - 1)];
    double h21 = H[(p - 1) * n + (p - 2)];
    double h22 = H[(p - 1) * n + (p - 1)];
    double s = h11 + h22;                       /* trace of trailing 2x2 */
    double t = h11 * h22 - h12 * h21;           /* det of trailing 2x2   */

    /* First three entries of M's first column for the active block. */
    double g11 = H[q * n + q];
    double g12 = H[q * n + (q + 1)];
    double g21 = H[(q + 1) * n + q];
    double g22 = H[(q + 1) * n + (q + 1)];
    double g32 = H[(q + 2) * n + (q + 1)];
    double x = g11 * g11 + g12 * g21 - s * g11 + t;
    double y = g21 * (g11 + g22 - s);
    double z = g21 * g32;

    for (size_t k = q; k + 2 < p; k++) {
        /* 3-element Householder for (x, y, z). */
        double sig = x * x + y * y + z * z;
        if (sig == 0.0) {
            /* Skip degenerate step; advance and continue. */
            x = H[(k + 1) * n + k];
            y = H[(k + 2) * n + k];
            z = (k + 3 < p) ? H[(k + 3) * n + k] : 0.0;
            continue;
        }
        double norm_v = sqrt(sig);
        double a = (x >= 0.0) ? -norm_v : norm_v;
        double u1 = x - a;
        double u2 = y;
        double u3 = z;
        double un2 = u1 * u1 + u2 * u2 + u3 * u3;
        double un = sqrt(un2);
        u1 /= un; u2 /= un; u3 /= un;

        /* Apply P from the left to rows (k, k+1, k+2). */
        size_t col_start = (k > q) ? k - 1 : q;
        for (size_t j = col_start; j < n; j++) {
            double s0 = (u1 * H[k * n + j]
                       + u2 * H[(k + 1) * n + j]
                       + u3 * H[(k + 2) * n + j]) * 2.0;
            H[k * n + j]       -= s0 * u1;
            H[(k + 1) * n + j] -= s0 * u2;
            H[(k + 2) * n + j] -= s0 * u3;
        }
        /* Apply P from the right to cols (k, k+1, k+2). */
        size_t row_end = (k + 3 < p) ? (k + 3) : (p - 1);
        for (size_t i = 0; i <= row_end; i++) {
            double s0 = (H[i * n + k]       * u1
                       + H[i * n + (k + 1)] * u2
                       + H[i * n + (k + 2)] * u3) * 2.0;
            H[i * n + k]       -= s0 * u1;
            H[i * n + (k + 1)] -= s0 * u2;
            H[i * n + (k + 2)] -= s0 * u3;
        }
        /* Q <- Q * P (right-multiply by reflector on cols k, k+1, k+2). */
        if (Q) {
            for (size_t i = 0; i < n; i++) {
                double s0 = (Q[i * n + k]       * u1
                           + Q[i * n + (k + 1)] * u2
                           + Q[i * n + (k + 2)] * u3) * 2.0;
                Q[i * n + k]       -= s0 * u1;
                Q[i * n + (k + 1)] -= s0 * u2;
                Q[i * n + (k + 2)] -= s0 * u3;
            }
        }

        x = H[(k + 1) * n + k];
        y = H[(k + 2) * n + k];
        z = (k + 3 < p) ? H[(k + 3) * n + k] : 0.0;
    }

    /* Final 2-element Householder on (x, y) at rows (p-2, p-1). */
    {
        size_t k = p - 2;
        double norm_v = sqrt(x * x + y * y);
        if (norm_v == 0.0) return;
        double a = (x >= 0.0) ? -norm_v : norm_v;
        double u1 = x - a;
        double u2 = y;
        double un2 = u1 * u1 + u2 * u2;
        double un = sqrt(un2);
        u1 /= un; u2 /= un;
        size_t col_start = k - 1;
        for (size_t j = col_start; j < n; j++) {
            double s0 = (u1 * H[k * n + j]
                       + u2 * H[(k + 1) * n + j]) * 2.0;
            H[k * n + j]       -= s0 * u1;
            H[(k + 1) * n + j] -= s0 * u2;
        }
        for (size_t i = 0; i < p; i++) {
            double s0 = (H[i * n + k] * u1
                       + H[i * n + (k + 1)] * u2) * 2.0;
            H[i * n + k]       -= s0 * u1;
            H[i * n + (k + 1)] -= s0 * u2;
        }
        if (Q) {
            for (size_t i = 0; i < n; i++) {
                double s0 = (Q[i * n + k]       * u1
                           + Q[i * n + (k + 1)] * u2) * 2.0;
                Q[i * n + k]       -= s0 * u1;
                Q[i * n + (k + 1)] -= s0 * u2;
            }
        }
    }
}

/* When a trailing 2x2 block has REAL eigenvalues we need to triangularise
 * it (drive the sub-diagonal to zero) before deflating, otherwise the
 * Schur form is non-triangular and back-substitution returns the wrong
 * eigenvector.  Apply a Givens rotation whose first column is the
 * eigenvector for one of the two real eigenvalues; after the similarity
 * the (i+1, i) entry is zero and the two real eigenvalues sit on the
 * diagonal.  Caller deflates them as two 1x1 blocks on the next QR
 * iteration.
 *
 * The block is at rows / cols (p-2, p-1).  Q is right-multiplied when
 * non-NULL so the caller can recover eigenvectors via back-transform. */
static void direct_split_2x2_real(double* H, size_t n, size_t p, double* Q) {
    size_t i = p - 2;
    double a = H[i * n + i];
    double b = H[i * n + (i + 1)];
    double c = H[(i + 1) * n + i];
    double d = H[(i + 1) * n + (i + 1)];
    double tr = a + d;
    double det = a * d - b * c;
    double disc = tr * tr - 4.0 * det;
    if (disc < 0.0) return;             /* complex; caller handles it */
    double sq = sqrt(disc);
    double lam = (tr + sq) * 0.5;       /* larger eigenvalue */
    double v0 = lam - d;
    double v1 = c;
    double r = sqrt(v0 * v0 + v1 * v1);
    if (r == 0.0) return;
    double cs = v0 / r;
    double sn = v1 / r;
    /* Left-multiply rows (i, i+1) of H by G^T = [[cs, sn], [-sn, cs]]. */
    for (size_t j = 0; j < n; j++) {
        double r0 = H[i * n + j];
        double r1 = H[(i + 1) * n + j];
        H[i * n + j]       =  cs * r0 + sn * r1;
        H[(i + 1) * n + j] = -sn * r0 + cs * r1;
    }
    /* Right-multiply cols (i, i+1) of H by G = [[cs, -sn], [sn, cs]]. */
    for (size_t k = 0; k < n; k++) {
        double c0 = H[k * n + i];
        double c1 = H[k * n + (i + 1)];
        H[k * n + i]       =  cs * c0 + sn * c1;
        H[k * n + (i + 1)] = -sn * c0 + cs * c1;
    }
    if (Q) {
        for (size_t k = 0; k < n; k++) {
            double c0 = Q[k * n + i];
            double c1 = Q[k * n + (i + 1)];
            Q[k * n + i]       =  cs * c0 + sn * c1;
            Q[k * n + (i + 1)] = -sn * c0 + cs * c1;
        }
    }
    H[(i + 1) * n + i] = 0.0;
}

/* Drive Francis QR sweeps on an upper Hessenberg matrix `H` until it
 * reaches quasi-triangular (Schur) form.  Eigenvalues are written into
 * eval_re / eval_im in Schur position order (NOT in sort order -- the
 * caller is responsible for sorting via direct_sort_perm_desc_abs_complex).
 *
 * When `Q` is non-NULL each bulge-chase reflector right-multiplies Q so
 * the caller can extract eigenvectors from the final Schur form via
 * back-substitution + back-transformation.
 *
 * Returns 0 on success, -1 if the maximum number of total Francis
 * sweeps (30 * n) is exceeded without convergence.
 *
 * LAPACK-HOOK: dhseqr. */
int direct_qr_real_general(double* H, size_t n,
                                    double* eval_re, double* eval_im,
                                    double* Q) {
    const double eps = 1e-14;
    const size_t max_iter = 30 * n;
    size_t iter = 0;
    size_t p = n;

    while (p > 0) {
        if (p == 1) {
            eval_re[0] = H[0];
            eval_im[0] = 0.0;
            break;
        }
        if (iter++ > max_iter) return -1;

        /* Find largest q in [1..p-1] such that H[q, q-1] is negligible
         * (deflation point).  q = 0 means the whole [0..p-1] block is
         * unreduced. */
        size_t q = p - 1;
        while (q > 0) {
            double tol = eps * (fabs(H[(q - 1) * n + (q - 1)])
                              + fabs(H[q * n + q]));
            if (fabs(H[q * n + (q - 1)]) <= tol) {
                H[q * n + (q - 1)] = 0.0;
                break;
            }
            q -= 1;
        }

        if (q == p - 1) {
            /* Trailing 1x1 block deflated -- real eigenvalue. */
            eval_re[p - 1] = H[(p - 1) * n + (p - 1)];
            eval_im[p - 1] = 0.0;
            p -= 1;
            iter = 0;
            continue;
        }
        if (q == p - 2) {
            /* Trailing 2x2 block. */
            double a = H[(p - 2) * n + (p - 2)];
            double b = H[(p - 2) * n + (p - 1)];
            double c = H[(p - 1) * n + (p - 2)];
            double d = H[(p - 1) * n + (p - 1)];
            double tr = a + d;
            double det = a * d - b * c;
            double disc = tr * tr - 4.0 * det;
            if (disc < 0.0) {
                /* Complex conjugate pair: deflate the entire 2x2 block. */
                double sq = sqrt(-disc);
                eval_re[p - 2] = tr * 0.5;
                eval_im[p - 2] =  sq * 0.5;
                eval_re[p - 1] = tr * 0.5;
                eval_im[p - 1] = -sq * 0.5;
                p -= 2;
                iter = 0;
                continue;
            }
            /* Real eigenvalues.  Store the analytic values (more
             * accurate than re-reading after the Givens roundoff),
             * then triangularise the block via a Givens similarity
             * so the Schur form is properly upper triangular and
             * back-substitution can recover eigenvectors. */
            double sq2 = sqrt(disc);
            eval_re[p - 2] = (tr + sq2) * 0.5;
            eval_im[p - 2] = 0.0;
            eval_re[p - 1] = (tr - sq2) * 0.5;
            eval_im[p - 1] = 0.0;
            direct_split_2x2_real(H, n, p, Q);
            /* After the split H[p-2, p-2] = larger eigenvalue and
             * H[p-1, p-1] = smaller, with H[p-1, p-2] = 0.  Overwrite
             * the diagonal entries with the analytic eigenvalues so
             * back-substitution uses values matching eval_re. */
            H[(p - 2) * n + (p - 2)] = eval_re[p - 2];
            H[(p - 1) * n + (p - 1)] = eval_re[p - 1];
            p -= 2;
            iter = 0;
            continue;
        }

        direct_francis_step(H, n, q, p, Q);
    }
    return 0;
}

/* Sort permutation by descending |lambda| (stable) where each lambda
 * is given by its real and imaginary parts. */
void direct_sort_perm_desc_abs_complex(const double* re,
                                                const double* im,
                                                size_t n, size_t* perm) {
    for (size_t i = 0; i < n; i++) perm[i] = i;
    for (size_t i = 1; i < n; i++) {
        size_t cur = perm[i];
        double ac = hypot(re[cur], im[cur]);
        size_t j = i;
        while (j > 0) {
            size_t prev = perm[j - 1];
            double ap = hypot(re[prev], im[prev]);
            if (ap > ac || (ap == ac && prev < cur)) break;
            perm[j] = perm[j - 1];
            j--;
        }
        perm[j] = cur;
    }
}

/* Build a List of eigenvalues from real / imaginary parts in the order
 * given by `perm`.  Real eigenvalues become EXPR_REAL; complex pairs
 * become Complex[re, im]. */
Expr* direct_build_complex_eigenvalue_list(const double* re,
                                                    const double* im,
                                                    size_t n,
                                                    const size_t* perm) {
    Expr** items = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        size_t idx = perm[i];
        if (im[idx] == 0.0) {
            items[i] = expr_new_real(re[idx]);
        } else {
            Expr** comp_args = (Expr**)malloc(sizeof(Expr*) * 2);
            comp_args[0] = expr_new_real(re[idx]);
            comp_args[1] = expr_new_real(im[idx]);
            items[i] = expr_new_function(expr_new_symbol("Complex"),
                                          comp_args, 2);
            free(comp_args);
        }
    }
    Expr* out = expr_new_function(expr_new_symbol("List"), items, n);
    free(items);
    return out;
}

/* ---- Eigenvectors from Schur form ----------------------------------- *
 *                                                                       *
 *  After Hessenberg + Francis QR we have H = Q^T A Q in quasi-          *
 *  triangular Schur form (1x1 blocks for real eigenvalues, 2x2 blocks   *
 *  for complex conjugate pairs).  Eigenvectors of A are obtained by:    *
 *    1. Solving (H - lambda I) v_schur = 0 by back-substitution in the  *
 *       Schur basis (real and complex variants below).                  *
 *    2. Back-transforming: v_A = Q v_schur.                             *
 *    3. Normalising to unit 2-norm.                                     *
 *                                                                       *
 *  LAPACK-HOOK: dtrevc3 handles steps 1-2 in one call when USE_LAPACK   *
 *  is wired (dorghr materialises Q from the Hessenberg reflectors).     *
 * --------------------------------------------------------------------- */

/* Back-substitute for the eigenvector of a Schur quasi-triangular matrix
 * H corresponding to a real eigenvalue at Schur position k.  Writes the
 * Schur-basis eigenvector into v (length n; v[i] = 0 for i > k). */
static void schur_eigvec_real(const double* H, size_t n, size_t k,
                                double lambda, double* v) {
    for (size_t i = 0; i < n; i++) v[i] = 0.0;
    v[k] = 1.0;
    if (k == 0) return;

    size_t i = k;
    while (i > 0) {
        i--;
        bool is_2x2 = (i > 0) && (H[i * n + (i - 1)] != 0.0);
        if (is_2x2) {
            /* Solve 2x2 system coupling v[i-1] and v[i] (the rows of the
             * 2x2 block).  Right-hand side = -sum of already-known
             * components weighted by the off-diagonal entries. */
            double rhs1 = 0.0, rhs2 = 0.0;
            for (size_t j = i + 1; j <= k; j++) {
                rhs1 += H[(i - 1) * n + j] * v[j];
                rhs2 += H[i * n + j] * v[j];
            }
            rhs1 = -rhs1; rhs2 = -rhs2;
            double a = H[(i - 1) * n + (i - 1)] - lambda;
            double b = H[(i - 1) * n + i];
            double c = H[i * n + (i - 1)];
            double d = H[i * n + i] - lambda;
            double det = a * d - b * c;
            if (fabs(det) < 1e-300) {
                v[i - 1] = 0.0; v[i] = 0.0;
            } else {
                v[i - 1] = (d * rhs1 - b * rhs2) / det;
                v[i]     = (a * rhs2 - c * rhs1) / det;
            }
            i--;                            /* skip the row we just used */
        } else {
            double rhs = 0.0;
            for (size_t j = i + 1; j <= k; j++) rhs += H[i * n + j] * v[j];
            double diag = H[i * n + i] - lambda;
            if (fabs(diag) < 1e-300) {
                v[i] = 0.0;                 /* defective: leave zero */
            } else {
                v[i] = -rhs / diag;
            }
        }
    }
}

/* Back-substitute for the complex eigenvector of a Schur quasi-triangular
 * H corresponding to the complex eigenvalue lambda = a + i b (b > 0) at
 * Schur positions (k, k+1).  Writes real / imaginary parts of the Schur-
 * basis eigenvector into v_re / v_im (length n; both 0 above k+1). */
static void schur_eigvec_complex(const double* H, size_t n, size_t k,
                                   double a, double b,
                                   double* v_re, double* v_im) {
    for (size_t i = 0; i < n; i++) { v_re[i] = 0.0; v_im[i] = 0.0; }
    /* Initial values from the 2x2 block at (k, k+1):
     *   2x2 = [[alpha, beta], [gamma, delta]], eigenvalues lambda, conj(lambda)
     *   eigenvector for lambda: x_1 = 1, x_0 = (lambda - delta) / gamma.
     *   lambda - delta = (a - delta) + i b. */
    double delta = H[(k + 1) * n + (k + 1)];
    double gamma = H[(k + 1) * n + k];
    if (gamma == 0.0) {
        /* Shouldn't happen for an unreduced 2x2 block; bail. */
        v_re[k] = 1.0;
        v_re[k + 1] = 0.0;
        return;
    }
    v_re[k]     = (a - delta) / gamma;
    v_im[k]     = b / gamma;
    v_re[k + 1] = 1.0;
    v_im[k + 1] = 0.0;

    if (k == 0) return;
    size_t i = k;
    while (i > 0) {
        i--;
        bool is_2x2 = (i > 0) && (H[i * n + (i - 1)] != 0.0);
        if (is_2x2) {
            /* Two complex equations -> 2x2 complex system in
             *   (v[i-1], v[i]).  Solve via complex 2x2 inverse. */
            double rhs_u_top = 0.0, rhs_v_top = 0.0;
            double rhs_u_bot = 0.0, rhs_v_bot = 0.0;
            for (size_t j = i + 1; j <= k + 1; j++) {
                rhs_u_top += H[(i - 1) * n + j] * v_re[j];
                rhs_v_top += H[(i - 1) * n + j] * v_im[j];
                rhs_u_bot += H[i * n + j] * v_re[j];
                rhs_v_bot += H[i * n + j] * v_im[j];
            }
            /* Coeff matrix (complex):
             *   [(a11 - lambda)  a12 ]
             *   [   a21         (a22 - lambda)]
             * = [[A_re + i A_im, B_re + i B_im],
             *    [C_re + i C_im, D_re + i D_im]] */
            double A_re = H[(i - 1) * n + (i - 1)] - a, A_im = -b;
            double B_re = H[(i - 1) * n + i],            B_im = 0.0;
            double C_re = H[i * n + (i - 1)],            C_im = 0.0;
            double D_re = H[i * n + i] - a,              D_im = -b;
            /* det = A D - B C  (complex) */
            double det_re = (A_re * D_re - A_im * D_im)
                          - (B_re * C_re - B_im * C_im);
            double det_im = (A_re * D_im + A_im * D_re)
                          - (B_re * C_im + B_im * C_re);
            double det_mag2 = det_re * det_re + det_im * det_im;
            if (det_mag2 < 1e-300) {
                v_re[i - 1] = 0.0; v_im[i - 1] = 0.0;
                v_re[i] = 0.0; v_im[i] = 0.0;
            } else {
                /* RHS (negated). */
                double r1_re = -rhs_u_top, r1_im = -rhs_v_top;
                double r2_re = -rhs_u_bot, r2_im = -rhs_v_bot;
                /* v[i-1] = (D r1 - B r2) / det */
                double num1_re = (D_re * r1_re - D_im * r1_im)
                               - (B_re * r2_re - B_im * r2_im);
                double num1_im = (D_re * r1_im + D_im * r1_re)
                               - (B_re * r2_im + B_im * r2_re);
                /* v[i]   = (A r2 - C r1) / det */
                double num2_re = (A_re * r2_re - A_im * r2_im)
                               - (C_re * r1_re - C_im * r1_im);
                double num2_im = (A_re * r2_im + A_im * r2_re)
                               - (C_re * r1_im + C_im * r1_re);
                /* (p + iq) / (det_re + i det_im) = (p det_re + q det_im
                 *   + i (q det_re - p det_im)) / det_mag2 */
                v_re[i - 1] = (num1_re * det_re + num1_im * det_im) / det_mag2;
                v_im[i - 1] = (num1_im * det_re - num1_re * det_im) / det_mag2;
                v_re[i]     = (num2_re * det_re + num2_im * det_im) / det_mag2;
                v_im[i]     = (num2_im * det_re - num2_re * det_im) / det_mag2;
            }
            i--;
        } else {
            /* Single complex equation: (H[i,i] - lambda) v[i] = -rhs. */
            double rhs_u = 0.0, rhs_v = 0.0;
            for (size_t j = i + 1; j <= k + 1; j++) {
                rhs_u += H[i * n + j] * v_re[j];
                rhs_v += H[i * n + j] * v_im[j];
            }
            double diag_re = H[i * n + i] - a;
            double diag_im = -b;
            double denom = diag_re * diag_re + diag_im * diag_im;
            if (denom < 1e-300) {
                v_re[i] = 0.0; v_im[i] = 0.0;
            } else {
                double num_re = -rhs_u, num_im = -rhs_v;
                v_re[i] = (num_re * diag_re + num_im * diag_im) / denom;
                v_im[i] = (num_im * diag_re - num_re * diag_im) / denom;
            }
        }
    }
}

/* Compute eigenvectors of A from the Schur form H and back-transformation
 * Q, in sorted (descending |lambda|) order.  V_re / V_im are n x n
 * row-major arrays where row k = sorted-k-th eigenvector. */
void schur_compute_eigvecs(const double* H, const double* Q,
                                    size_t n,
                                    const double* eval_re, const double* eval_im,
                                    const size_t* perm,
                                    double* V_re, double* V_im) {
    double* v_schur_re = (double*)malloc(sizeof(double) * n);
    double* v_schur_im = (double*)malloc(sizeof(double) * n);
    double* w_re       = (double*)malloc(sizeof(double) * n);
    double* w_im       = (double*)malloc(sizeof(double) * n);

    size_t* inv_perm = (size_t*)malloc(sizeof(size_t) * n);
    for (size_t i = 0; i < n; i++) inv_perm[perm[i]] = i;

    size_t k = 0;
    while (k < n) {
        if (eval_im[k] == 0.0) {
            /* Real eigenvalue at Schur position k. */
            schur_eigvec_real(H, n, k, eval_re[k], v_schur_re);
            /* w = Q . v_schur. */
            for (size_t i = 0; i < n; i++) {
                double s = 0.0;
                for (size_t j = 0; j <= k; j++) s += Q[i * n + j] * v_schur_re[j];
                w_re[i] = s;
            }
            double norm2 = 0.0;
            for (size_t i = 0; i < n; i++) norm2 += w_re[i] * w_re[i];
            double inv = (norm2 > 0.0) ? 1.0 / sqrt(norm2) : 1.0;
            size_t sp = inv_perm[k];
            for (size_t i = 0; i < n; i++) {
                V_re[sp * n + i] = w_re[i] * inv;
                V_im[sp * n + i] = 0.0;
            }
            k++;
        } else {
            /* Complex pair at Schur positions k, k+1.  eval_im[k] is the
             * "+imag" root (the QR loop writes them in this order). */
            double a = eval_re[k];
            double b = fabs(eval_im[k]);
            schur_eigvec_complex(H, n, k, a, b, v_schur_re, v_schur_im);
            for (size_t i = 0; i < n; i++) {
                double s_re = 0.0, s_im = 0.0;
                for (size_t j = 0; j <= k + 1; j++) {
                    s_re += Q[i * n + j] * v_schur_re[j];
                    s_im += Q[i * n + j] * v_schur_im[j];
                }
                w_re[i] = s_re;
                w_im[i] = s_im;
            }
            double norm2 = 0.0;
            for (size_t i = 0; i < n; i++) {
                norm2 += w_re[i] * w_re[i] + w_im[i] * w_im[i];
            }
            double inv = (norm2 > 0.0) ? 1.0 / sqrt(norm2) : 1.0;
            size_t sp1 = inv_perm[k];       /* +imag eigenvalue */
            size_t sp2 = inv_perm[k + 1];   /* -imag eigenvalue */
            for (size_t i = 0; i < n; i++) {
                double r = w_re[i] * inv;
                double m = w_im[i] * inv;
                V_re[sp1 * n + i] = r;  V_im[sp1 * n + i] =  m;
                V_re[sp2 * n + i] = r;  V_im[sp2 * n + i] = -m;
            }
            k += 2;
        }
    }

    free(v_schur_re); free(v_schur_im);
    free(w_re); free(w_im);
    free(inv_perm);
}

/* Emit a List of List of (Real or Complex[re, im]) for the eigenvector
 * matrix V (rows = eigenvectors).  Entries with V_im[i,j] == 0 become
 * Real; others become Complex[re, im]. */
static Expr* direct_build_complex_eigenvector_list(const double* V_re,
                                                     const double* V_im,
                                                     size_t n) {
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t k = 0; k < n; k++) {
        Expr** comps = (Expr**)malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            double r = V_re[k * n + i];
            double m = V_im[k * n + i];
            if (m == 0.0) {
                comps[i] = expr_new_real(r);
            } else {
                Expr** args = (Expr**)malloc(sizeof(Expr*) * 2);
                args[0] = expr_new_real(r);
                args[1] = expr_new_real(m);
                comps[i] = expr_new_function(expr_new_symbol("Complex"), args, 2);
                free(args);
            }
        }
        rows[k] = expr_new_function(expr_new_symbol("List"), comps, n);
        free(comps);
    }
    Expr* out = expr_new_function(expr_new_symbol("List"), rows, n);
    free(rows);
    return out;
}

/* Top-level "Direct" kernel for real non-symmetric machine-precision
 * input.  When WANT_VALUES, returns a List of eigenvalues sorted by
 * descending |lambda| (Real or Complex[re, im]).  When WANT_VECTORS,
 * returns a List of List of eigenvector components in the same sorted
 * order.  Returns NULL on convergence failure -- caller falls back to
 * the symbolic path. */
static Expr* direct_real_general_machine(const MatD* A, MateigenWant want,
                                          Expr* k_spec) {
    size_t n = A->n;
    if (n == 0) return NULL;

    bool want_Q = (want & MATEIGEN_WANT_VECTORS) != 0;

    double* H = (double*)malloc(sizeof(double) * n * n);
    memcpy(H, A->re, sizeof(double) * n * n);
    double* Q = NULL;
    if (want_Q) {
        Q = (double*)malloc(sizeof(double) * n * n);
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++) Q[i * n + j] = (i == j) ? 1.0 : 0.0;
    }
    double* u_buf = (double*)malloc(sizeof(double) * n);
    direct_hessenberg_real(H, n, u_buf, Q);
    free(u_buf);

    double* eval_re = (double*)calloc(n, sizeof(double));
    double* eval_im = (double*)calloc(n, sizeof(double));
    int qr_status = direct_qr_real_general(H, n, eval_re, eval_im, Q);

    if (qr_status != 0) {
        free(H); free(eval_re); free(eval_im);
        if (Q) free(Q);
        return NULL;
    }

    size_t* perm = (size_t*)malloc(sizeof(size_t) * n);
    direct_sort_perm_desc_abs_complex(eval_re, eval_im, n, perm);

    Expr* out;
    if (want_Q) {
        double* V_re = (double*)malloc(sizeof(double) * n * n);
        double* V_im = (double*)malloc(sizeof(double) * n * n);
        schur_compute_eigvecs(H, Q, n, eval_re, eval_im, perm, V_re, V_im);
        out = direct_build_complex_eigenvector_list(V_re, V_im, n);
        free(V_re); free(V_im);
    } else {
        out = direct_build_complex_eigenvalue_list(eval_re, eval_im, n, perm);
    }

    free(H); free(eval_re); free(eval_im); free(perm);
    if (Q) free(Q);

    return direct_apply_k_spec_list(out, k_spec);
}

/* Top-level "Direct" kernel for real symmetric machine-precision input.
 * Returns a freshly-allocated List[Real,...] (eigenvalues) or
 * List[List[Real,...],...] (eigenvectors), or NULL when the input is
 * outside this kernel's supported domain so the dispatcher can fall
 * back to a different kernel / the symbolic path.  Caller owns the
 * result. */
static Expr* direct_real_sym_machine(const MatD* A, MateigenWant want,
                                       Expr* k_spec) {
    size_t n = A->n;
    if (n == 0) return NULL;

    /* Working copy of A (tridiag step modifies in place). */
    double* W = (double*)malloc(sizeof(double) * n * n);
    memcpy(W, A->re, sizeof(double) * n * n);

    double* diag    = (double*)malloc(sizeof(double) * n);
    double* sub     = (double*)calloc(n, sizeof(double));  /* length n-1, +1 slack */
    double* u       = (double*)malloc(sizeof(double) * n);
    double* p_buf   = (double*)malloc(sizeof(double) * n);
    double* q_buf   = (double*)malloc(sizeof(double) * n);
    bool want_Q     = (want & MATEIGEN_WANT_VECTORS) != 0;
    double* Q       = want_Q ? (double*)malloc(sizeof(double) * n * n) : NULL;

    direct_tridiag_real_sym(W, n, diag, sub, Q, want_Q, u, p_buf, q_buf);

    int qr_status = direct_symtridiag_qr(diag, sub, n, Q, want_Q);

    free(u); free(p_buf); free(q_buf); free(W);

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

/* ---------- Complex Hermitian Direct: Householder + phase + sym QR -- *
 *                                                                       *
 *  For a complex Hermitian n x n matrix A all eigenvalues are real and  *
 *  the eigenvectors form a complex unitary basis.  We reduce as:        *
 *                                                                       *
 *    1. Complex Householder reflectors zero the sub-column below the    *
 *       sub-diagonal at each step.  The result is a Hermitian tri-      *
 *       diagonal T with real diagonal but generally COMPLEX             *
 *       sub-diagonal.  Q (n x n complex) accumulates the reflectors.    *
 *    2. Diagonal phase correction D = diag(d_0, ..., d_{n-1}), |d_k|=1, *
 *       chosen so D^H T D has real-positive sub-diagonal.  T becomes a  *
 *       real symmetric tridiagonal.  Q is updated to Q D.               *
 *    3. The existing real symmetric tridiag QR (direct_symtridiag_qr)   *
 *       finds eigenvalues and the real orthogonal accumulator Z so that *
 *       Z^T T_real Z = Lambda.                                          *
 *    4. Final eigenvectors V = Q Z (complex Q times real Z).            *
 *                                                                       *
 *  This avoids needing a separate complex tridiagonal QR while costing  *
 *  only an O(n^2) diagonal-phase application step.  See Wilkinson 1965, *
 *  "The Algebraic Eigenvalue Problem", section 5.45.                    *
 *                                                                       *
 *  LAPACK-HOOK: this whole block maps to zhetrd (step 1) + the implicit *
 *  phase reduction (handled internally by LAPACK's zhetrd via the tau   *
 *  scalar) + dstedc/dsteqr (step 3), or simply zheevd as a one-call     *
 *  wrapper when USE_LAPACK is wired.                                    *
 * --------------------------------------------------------------------- */

/* Hermitian Householder tridiagonalisation.
 *
 * Input:  A (row-major n*n; A_re holds real parts, A_im holds imag parts).
 * Output: diag[i]     = real diagonal of T  (0 <= i < n)
 *         sub_re[k], sub_im[k] = complex sub-diagonal T_{k+1,k}, 0 <= k < n-1
 *         Q (complex, n*n via Q_re/Q_im) when want_Q: unitary accumulator
 *
 * Scratch buffers u, v, q (each length n complex via paired arrays) are
 * caller-provided so the inner loops never touch malloc/free. */
void direct_tridiag_complex_hermitian(double* A_re, double* A_im,
                                              size_t n,
                                              double* diag,
                                              double* sub_re, double* sub_im,
                                              double* Q_re, double* Q_im,
                                              bool want_Q,
                                              double* u_re, double* u_im,
                                              double* v_re, double* v_im,
                                              double* q_re, double* q_im) {
    if (want_Q) {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                Q_re[i * n + j] = (i == j) ? 1.0 : 0.0;
                Q_im[i * n + j] = 0.0;
            }
        }
    }

    for (size_t k = 0; k + 2 < n; k++) {
        /* Norm of column-k tail x = A[k+1:n, k] (complex). */
        double sigma = 0.0;
        for (size_t i = k + 1; i < n; i++) {
            double r = A_re[i * n + k];
            double m = A_im[i * n + k];
            sigma += r * r + m * m;
        }
        double xk1_re = A_re[(k + 1) * n + k];
        double xk1_im = A_im[(k + 1) * n + k];
        double xk1_abs = hypot(xk1_re, xk1_im);

        if (sigma == 0.0) {
            sub_re[k] = xk1_re;
            sub_im[k] = xk1_im;
            continue;
        }
        double norm_x = sqrt(sigma);

        /* alpha = -phase(x_{k+1}) * norm_x.  Phase = x/|x| (or 1 if x=0).
         * Choosing alpha "opposite" the leading entry maximises |u_0|^2 so
         * the reflector is numerically well-conditioned. */
        double alpha_re, alpha_im;
        if (xk1_abs == 0.0) {
            alpha_re = -norm_x;
            alpha_im = 0.0;
        } else {
            alpha_re = -xk1_re / xk1_abs * norm_x;
            alpha_im = -xk1_im / xk1_abs * norm_x;
        }

        /* u = x with u_{k+1} -= alpha.  Householder reflector
         * H = I - 2 u u^H / ||u||^2; we normalise u so that tau = 2. */
        u_re[k + 1] = xk1_re - alpha_re;
        u_im[k + 1] = xk1_im - alpha_im;
        for (size_t i = k + 2; i < n; i++) {
            u_re[i] = A_re[i * n + k];
            u_im[i] = A_im[i * n + k];
        }
        double unorm2 = 0.0;
        for (size_t i = k + 1; i < n; i++) {
            unorm2 += u_re[i] * u_re[i] + u_im[i] * u_im[i];
        }
        if (unorm2 == 0.0) {
            sub_re[k] = xk1_re;
            sub_im[k] = xk1_im;
            continue;
        }
        double inv_unorm = 1.0 / sqrt(unorm2);
        for (size_t i = k + 1; i < n; i++) {
            u_re[i] *= inv_unorm;
            u_im[i] *= inv_unorm;
        }

        /* v = A u   on the trailing sub-block (A is Hermitian; we read
         * the full sub-block to be agnostic to whether the upper triangle
         * has been kept in sync). */
        for (size_t i = k + 1; i < n; i++) {
            double s_re = 0.0, s_im = 0.0;
            for (size_t j = k + 1; j < n; j++) {
                double ar = A_re[i * n + j];
                double ai = A_im[i * n + j];
                double ur = u_re[j];
                double ui = u_im[j];
                /* s += A_ij * u_j */
                s_re += ar * ur - ai * ui;
                s_im += ar * ui + ai * ur;
            }
            v_re[i] = s_re;
            v_im[i] = s_im;
        }

        /* alpha_v = u^H v = sum conj(u_i) v_i.  Real for Hermitian A. */
        double alpha_v = 0.0;
        for (size_t i = k + 1; i < n; i++) {
            alpha_v += u_re[i] * v_re[i] + u_im[i] * v_im[i];
        }

        /* q = 2 v - 2 alpha_v u   so that  A <- A - u q^H - q u^H = H A H. */
        for (size_t i = k + 1; i < n; i++) {
            q_re[i] = 2.0 * v_re[i] - 2.0 * alpha_v * u_re[i];
            q_im[i] = 2.0 * v_im[i] - 2.0 * alpha_v * u_im[i];
        }

        /* Rank-2 update: A_ij -= u_i conj(q_j) + q_i conj(u_j). */
        for (size_t i = k + 1; i < n; i++) {
            for (size_t j = k + 1; j < n; j++) {
                double ur = u_re[i], ui = u_im[i];
                double qr = q_re[i], qi = q_im[i];
                double uj_r = u_re[j], uj_i = u_im[j];
                double qj_r = q_re[j], qj_i = q_im[j];
                /* u_i * conj(q_j) = (ur + i ui)(qj_r - i qj_i) */
                double t1_re = ur * qj_r + ui * qj_i;
                double t1_im = ui * qj_r - ur * qj_i;
                /* q_i * conj(u_j) = (qr + i qi)(uj_r - i uj_i) */
                double t2_re = qr * uj_r + qi * uj_i;
                double t2_im = qi * uj_r - qr * uj_i;
                A_re[i * n + j] -= t1_re + t2_re;
                A_im[i * n + j] -= t1_im + t2_im;
            }
        }

        /* Force the new sub-column / sub-row to their analytic values
         * to suppress drift.  Sub-column has only one non-zero entry
         * after H is applied: A[k+1, k] = alpha (complex).  All A[i, k]
         * for i > k+1 are exactly zero; conjugate row equally. */
        sub_re[k] = alpha_re;
        sub_im[k] = alpha_im;
        A_re[(k + 1) * n + k] = alpha_re;
        A_im[(k + 1) * n + k] = alpha_im;
        A_re[k * n + (k + 1)] = alpha_re;   /* T_{k,k+1} = conj(T_{k+1,k}) */
        A_im[k * n + (k + 1)] = -alpha_im;
        for (size_t i = k + 2; i < n; i++) {
            A_re[i * n + k] = 0.0;
            A_im[i * n + k] = 0.0;
            A_re[k * n + i] = 0.0;
            A_im[k * n + i] = 0.0;
        }

        /* Q <- Q H (right-multiply).  Column j -> col j - 2 (Q row * u) u_j.
         * In complex: Q_ij -= 2 (Q row . u)_i * conj(u_j) -- wait, this
         * isn't right.  H = I - 2 u u^H, so Q H = Q - 2 (Q u) u^H, which
         * in scalar form is Q_ij <- Q_ij - 2 (Q u)_i conj(u_j). */
        if (want_Q) {
            for (size_t i = 0; i < n; i++) {
                /* (Q u)_i = sum_j Q_ij u_j, restricted to j >= k+1. */
                double s_re = 0.0, s_im = 0.0;
                for (size_t j = k + 1; j < n; j++) {
                    double qr2 = Q_re[i * n + j];
                    double qi2 = Q_im[i * n + j];
                    double ur2 = u_re[j];
                    double ui2 = u_im[j];
                    s_re += qr2 * ur2 - qi2 * ui2;
                    s_im += qr2 * ui2 + qi2 * ur2;
                }
                s_re *= 2.0;
                s_im *= 2.0;
                /* Q_ij -= s * conj(u_j) = s * (u_re[j] - i u_im[j]) */
                for (size_t j = k + 1; j < n; j++) {
                    double ur2 = u_re[j];
                    double ui2 = u_im[j];
                    /* (s_re + i s_im)(ur2 - i ui2) */
                    double pr = s_re * ur2 + s_im * ui2;
                    double pi = s_im * ur2 - s_re * ui2;
                    Q_re[i * n + j] -= pr;
                    Q_im[i * n + j] -= pi;
                }
            }
        }
    }

    /* Extract diagonal (real for Hermitian, but defensively take the
     * real part). */
    for (size_t i = 0; i < n; i++) diag[i] = A_re[i * n + i];
    /* The (n-2 -> n-1) sub-diagonal entry hasn't been overwritten yet. */
    if (n >= 2) {
        sub_re[n - 2] = A_re[(n - 1) * n + (n - 2)];
        sub_im[n - 2] = A_im[(n - 1) * n + (n - 2)];
    }
}

/* Diagonal phase correction.  Multiplies the complex Hermitian tri-
 * diagonal T by a diagonal unitary D = diag(d_0, ..., d_{n-1}), |d_k|=1,
 * chosen so D^H T D has real-positive sub-diagonal.  Updates Q (the
 * Householder accumulator) by post-multiplication: Q <- Q D. */
void direct_phase_correct_tridiag(double* sub_re, double* sub_im,
                                           size_t n,
                                           double* Q_re, double* Q_im,
                                           bool want_Q) {
    /* d_0 = 1; d_{k+1} = d_k * sub[k] / |sub[k]|. */
    double d_re = 1.0, d_im = 0.0;
    /* k = 0 -> column 1 of Q gets multiplied by d_1, column 0 stays. */
    for (size_t k = 0; k + 1 < n; k++) {
        double sr = sub_re[k], si = sub_im[k];
        double mag = hypot(sr, si);
        double phase_re, phase_im;
        if (mag == 0.0) { phase_re = 1.0; phase_im = 0.0; }
        else            { phase_re = sr / mag; phase_im = si / mag; }
        /* d_{k+1} = d_k * phase */
        double new_d_re = d_re * phase_re - d_im * phase_im;
        double new_d_im = d_re * phase_im + d_im * phase_re;
        d_re = new_d_re;
        d_im = new_d_im;
        sub_re[k] = mag;
        sub_im[k] = 0.0;
        if (want_Q) {
            for (size_t i = 0; i < n; i++) {
                double qr = Q_re[i * n + (k + 1)];
                double qi = Q_im[i * n + (k + 1)];
                Q_re[i * n + (k + 1)] = qr * d_re - qi * d_im;
                Q_im[i * n + (k + 1)] = qr * d_im + qi * d_re;
            }
        }
    }
}

/* Compose complex Q (n*n) with real Z (n*n) into complex V (n*n) via
 * V = Q Z.  V_re/V_im are caller-allocated. */
void direct_compose_complex_Q_real_Z(const double* Q_re,
                                              const double* Q_im,
                                              const double* Z, size_t n,
                                              double* V_re, double* V_im) {
    for (size_t i = 0; i < n; i++) {
        for (size_t k = 0; k < n; k++) {
            double sr = 0.0, si = 0.0;
            for (size_t j = 0; j < n; j++) {
                double z = Z[j * n + k];
                sr += Q_re[i * n + j] * z;
                si += Q_im[i * n + j] * z;
            }
            V_re[i * n + k] = sr;
            V_im[i * n + k] = si;
        }
    }
}

/* Build a List of List of (Real or Complex[re, im]) from a complex
 * eigenvector matrix V (n*n; V[i,k] is the i-th component of the k-th
 * eigenvector) in the order given by `perm` (perm[r] -> column of V). */
Expr* direct_build_complex_hermitian_eigvec_list(const double* V_re,
                                                          const double* V_im,
                                                          size_t n,
                                                          const size_t* perm) {
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t r = 0; r < n; r++) {
        size_t col = perm[r];
        Expr** comps = (Expr**)malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            double rv = V_re[i * n + col];
            double iv = V_im[i * n + col];
            if (iv == 0.0) {
                comps[i] = expr_new_real(rv);
            } else {
                Expr** args = (Expr**)malloc(sizeof(Expr*) * 2);
                args[0] = expr_new_real(rv);
                args[1] = expr_new_real(iv);
                comps[i] = expr_new_function(expr_new_symbol("Complex"), args, 2);
                free(args);
            }
        }
        rows[r] = expr_new_function(expr_new_symbol("List"), comps, n);
        free(comps);
    }
    Expr* out = expr_new_function(expr_new_symbol("List"), rows, n);
    free(rows);
    return out;
}

/* Top-level "Direct" kernel for complex Hermitian machine-precision
 * input.  Eigenvalues are real (sorted by descending |lambda|);
 * eigenvectors are complex unitary.  Returns NULL on convergence
 * failure so the caller falls back to the symbolic path. */
static Expr* direct_complex_hermitian_machine(const MatD* A, MateigenWant want,
                                                Expr* k_spec) {
    size_t n = A->n;
    if (n == 0) return NULL;

    /* Working complex copy of A (tridiag step modifies in place). */
    double* W_re = (double*)malloc(sizeof(double) * n * n);
    double* W_im = (double*)malloc(sizeof(double) * n * n);
    memcpy(W_re, A->re, sizeof(double) * n * n);
    memcpy(W_im, A->im, sizeof(double) * n * n);

    double* diag   = (double*)malloc(sizeof(double) * n);
    double* sub_re = (double*)calloc(n, sizeof(double));   /* len n-1 + slack */
    double* sub_im = (double*)calloc(n, sizeof(double));

    double* u_re = (double*)malloc(sizeof(double) * n);
    double* u_im = (double*)malloc(sizeof(double) * n);
    double* v_re = (double*)malloc(sizeof(double) * n);
    double* v_im = (double*)malloc(sizeof(double) * n);
    double* q_re = (double*)malloc(sizeof(double) * n);
    double* q_im = (double*)malloc(sizeof(double) * n);

    bool want_Q   = (want & MATEIGEN_WANT_VECTORS) != 0;
    double* Q_re  = want_Q ? (double*)malloc(sizeof(double) * n * n) : NULL;
    double* Q_im  = want_Q ? (double*)malloc(sizeof(double) * n * n) : NULL;

    direct_tridiag_complex_hermitian(W_re, W_im, n,
                                      diag, sub_re, sub_im,
                                      Q_re, Q_im, want_Q,
                                      u_re, u_im, v_re, v_im, q_re, q_im);

    free(W_re); free(W_im);
    free(u_re); free(u_im); free(v_re); free(v_im); free(q_re); free(q_im);

    /* Phase-correct so sub-diagonal becomes real positive. */
    direct_phase_correct_tridiag(sub_re, sub_im, n, Q_re, Q_im, want_Q);

    /* Real symmetric tridiagonal QR.  We accumulate the real orthogonal
     * Z separately when eigenvectors are wanted -- it composes with the
     * complex Q from the Hermitian Householder step. */
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
        /* V = Q Z (complex). */
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

/* ---------- Complex general Direct: real block-embedding ------------- *
 *                                                                       *
 *  A complex n x n matrix A = R + i S is reduced to a real 2n x 2n      *
 *  block matrix                                                         *
 *                                                                       *
 *      M = [[ R, -S ],                                                  *
 *           [ S,  R ]]                                                  *
 *                                                                       *
 *  whose spectrum is spec(A) ∪ conj(spec(A)) as a multiset, with each   *
 *  pair (λ, conj λ) sharing a 2-dim real invariant subspace of M.       *
 *                                                                       *
 *  We dispatch M through the existing real Hessenberg + Francis QR      *
 *  pipeline (direct_hessenberg_real + direct_qr_real_general +          *
 *  schur_compute_eigvecs) and then unpair:                              *
 *                                                                       *
 *    - For each complex M-eigenvalue with positive imaginary part:      *
 *      find its conjugate partner among the M-eigenvalues, mark both    *
 *      as used, and emit the positive-imag one as A's eigenvalue.       *
 *    - For each real M-eigenvalue: pair it with the closest unused      *
 *      real M-eigenvalue, mark both as used, emit one copy.             *
 *                                                                       *
 *  The complex A-eigenvector is recovered from M's complex eigenvector  *
 *  w = w_re + i w_im (length 2n) via                                    *
 *                                                                       *
 *      x = (a - d) + i (b + c)                                          *
 *                                                                       *
 *  where a, b, c, d are the n-vector top/bottom splits of w_re / w_im   *
 *  respectively.  For real M-eigenvalues w_im = 0 and the formula       *
 *  collapses to x = top(w_re) + i * bot(w_re).  See test                *
 *  test_direct_general_complex_machine_2x2_diagonal_imag for a worked   *
 *  example.                                                             *
 *                                                                       *
 *  Cost: O((2n)^3) = 8 O(n^3) -- nominally 8x more flops than a native  *
 *  complex Hessenberg + complex QR at the same n.  In return the entire *
 *  numerical machinery is reused (real LAPACK-mappable kernels) and     *
 *  there are no separate complex tridiagonal / Givens implementations   *
 *  to maintain.  When USE_LAPACK is wired, the obvious replacement is   *
 *  a native zgehrd + zhseqr + ztrevc3 path that won't pay the 8x.       *
 *                                                                       *
 *  LAPACK-HOOK: zgehrd + zhseqr + ztrevc3 (or zgeev as a one-call       *
 *  wrapper) can replace this whole block when USE_LAPACK is set.        *
 * --------------------------------------------------------------------- */

static Expr* direct_complex_general_machine(const MatD* A, MateigenWant want,
                                              Expr* k_spec) {
    size_t n = A->n;
    if (n == 0) return NULL;
    size_t N = 2 * n;

    bool want_Q = (want & MATEIGEN_WANT_VECTORS) != 0;

    /* Build real 2n x 2n block embedding. */
    double* H = (double*)malloc(sizeof(double) * N * N);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            double r = A->re[i * n + j];
            double s = A->im[i * n + j];
            H[i * N + j]                   =  r;
            H[i * N + (j + n)]             = -s;
            H[(i + n) * N + j]             =  s;
            H[(i + n) * N + (j + n)]       =  r;
        }
    }

    /* Hessenberg reduction.  We ALWAYS accumulate Q -- M's eigenvectors
     * are needed for the values-only path too, to disambiguate which of
     * the conjugate pair (mu, conj mu) is in spec(A) vs spec(conj A). */
    double* Q = (double*)malloc(sizeof(double) * N * N);
    for (size_t i = 0; i < N; i++)
        for (size_t j = 0; j < N; j++) Q[i * N + j] = (i == j) ? 1.0 : 0.0;
    double* u_buf = (double*)malloc(sizeof(double) * N);
    direct_hessenberg_real(H, N, u_buf, Q);
    free(u_buf);

    /* Francis QR -> Schur form + complex-conjugate eigenvalue pairs. */
    double* M_eval_re = (double*)calloc(N, sizeof(double));
    double* M_eval_im = (double*)calloc(N, sizeof(double));
    int qr_status = direct_qr_real_general(H, N, M_eval_re, M_eval_im, Q);
    if (qr_status != 0) {
        free(H); free(Q); free(M_eval_re); free(M_eval_im);
        return NULL;
    }

    /* Eigenvectors of M (in the original basis), one row per eigenvalue. */
    double* M_evec_re = (double*)malloc(sizeof(double) * N * N);
    double* M_evec_im = (double*)malloc(sizeof(double) * N * N);
    /* calloc (not malloc) so GCC -O3 cross-procedural -Wmaybe-uninitialized
     * doesn't fire after schur_compute_eigvecs is inlined into the dispatcher. */
    size_t* identity_perm = (size_t*)calloc(N, sizeof(size_t));
    for (size_t i = 0; i < N; i++) identity_perm[i] = i;
    schur_compute_eigvecs(H, Q, N, M_eval_re, M_eval_im, identity_perm,
                            M_evec_re, M_evec_im);
    free(identity_perm);
    free(H); free(Q);

    /* Recover A's eigenvalues from M's spectrum:
     *
     * For each M-eigenvector w of eigenvalue mu, the vector
     *
     *     x = (a - d) + i (b + c)
     *
     * with w_re = [a; b] and w_im = [c; d] satisfies A x = mu x.  When mu
     * is in spec(A) the candidate x is non-zero; when mu is only in
     * spec(conj A) (i.e., conj mu in spec(A) but mu itself is not)
     * the candidate x is zero.  This disambiguates the conjugate-pair
     * doubling of spec(M) = spec(A) U conj(spec(A)).
     *
     * Grouped Gram-Schmidt handles algebraic multiplicity:
     *   1. Walk M-eigenvalues; group adjacent (after stable Schur order)
     *      values that fall within `group_tol`.
     *   2. Within each group, project each candidate x against the
     *      already-emitted vectors of the group and emit it iff the
     *      remaining norm stays above `extract_threshold`.
     *
     * This produces m_A(mu) ortho-normal eigenvectors per distinct mu --
     * the rank of the +J subspace in M's mu-eigenspace -- which sums to
     * exactly n. */
    double spec_norm = 0.0;
    for (size_t i = 0; i < N; i++) {
        double a = hypot(M_eval_re[i], M_eval_im[i]);
        if (a > spec_norm) spec_norm = a;
    }
    double group_tol = 1e-8 * (spec_norm == 0.0 ? 1.0 : spec_norm) * (double)N;
    /* Norm threshold: a "valid" candidate (in +J subspace) has |x| in
     * the range [eps, sqrt(2)]; a "wrong" candidate (in -J subspace) has
     * |x| ~ machine_eps.  1e-8 is safely between these. */
    double extract_threshold = sqrt((double)n) * 1e-9;

    int* used = (int*)calloc(N, sizeof(int));
    double* A_eval_re = (double*)malloc(sizeof(double) * n);
    double* A_eval_im = (double*)malloc(sizeof(double) * n);
    double* A_evec_re = (double*)calloc(n * n, sizeof(double));
    double* A_evec_im = (double*)calloc(n * n, sizeof(double));
    double* cand_re = (double*)malloc(sizeof(double) * n);
    double* cand_im = (double*)malloc(sizeof(double) * n);
    size_t out = 0;

    for (size_t i = 0; i < N && out < n; i++) {
        if (used[i]) continue;
        used[i] = 1;
        size_t group_start = out;

        /* Process group member j; iteration includes i itself. */
        for (size_t j = i; j < N && out < n; j++) {
            if (j != i) {
                if (used[j]) continue;
                double dr = M_eval_re[j] - M_eval_re[i];
                double di = M_eval_im[j] - M_eval_im[i];
                if (hypot(dr, di) > group_tol) continue;
                used[j] = 1;
            }
            /* Candidate x = (a - d) + i (b + c). */
            for (size_t l = 0; l < n; l++) {
                double a = M_evec_re[j * N + l];
                double b = M_evec_re[j * N + (l + n)];
                double c = M_evec_im[j * N + l];
                double d_im = M_evec_im[j * N + (l + n)];
                cand_re[l] = a - d_im;
                cand_im[l] = b + c;
            }
            /* Orthogonalise against already-emitted vectors in this
             * group (complex modified Gram-Schmidt).  Twice for numerical
             * stability per the "twice is enough" rule. */
            for (int pass = 0; pass < 2; pass++) {
                for (size_t f = group_start; f < out; f++) {
                    double pr = 0.0, pi = 0.0;
                    for (size_t l = 0; l < n; l++) {
                        double vr = A_evec_re[f * n + l];
                        double vi = A_evec_im[f * n + l];
                        /* conj(V_f) . cand = (vr - i vi)(cand_re + i cand_im) */
                        pr += vr * cand_re[l] + vi * cand_im[l];
                        pi += vr * cand_im[l] - vi * cand_re[l];
                    }
                    for (size_t l = 0; l < n; l++) {
                        double vr = A_evec_re[f * n + l];
                        double vi = A_evec_im[f * n + l];
                        /* (pr + i pi)(vr + i vi) */
                        double pvr = pr * vr - pi * vi;
                        double pvi = pr * vi + pi * vr;
                        cand_re[l] -= pvr;
                        cand_im[l] -= pvi;
                    }
                }
            }
            double norm2 = 0.0;
            for (size_t l = 0; l < n; l++) {
                norm2 += cand_re[l] * cand_re[l] + cand_im[l] * cand_im[l];
            }
            if (norm2 < extract_threshold * extract_threshold) continue;
            double inv = 1.0 / sqrt(norm2);
            for (size_t l = 0; l < n; l++) {
                A_evec_re[out * n + l] = cand_re[l] * inv;
                A_evec_im[out * n + l] = cand_im[l] * inv;
            }
            A_eval_re[out] = M_eval_re[i];
            A_eval_im[out] = M_eval_im[i];
            out++;
        }
    }
    free(used); free(M_eval_re); free(M_eval_im);
    free(M_evec_re); free(M_evec_im);
    free(cand_re); free(cand_im);

    if (out != n) {
        /* Extraction under-produced -- bail out so the symbolic path
         * can take over.  Most commonly hit when extract_threshold is
         * too tight or the matrix is wildly ill-conditioned. */
        free(A_eval_re); free(A_eval_im);
        free(A_evec_re); free(A_evec_im);
        return NULL;
    }

    /* Sort by descending |lambda| and emit. */
    size_t* perm = (size_t*)malloc(sizeof(size_t) * n);
    direct_sort_perm_desc_abs_complex(A_eval_re, A_eval_im, n, perm);

    Expr* result;
    if (want_Q) {
        /* Permute eigenvector rows into sort order. */
        double* V_re = (double*)malloc(sizeof(double) * n * n);
        double* V_im = (double*)malloc(sizeof(double) * n * n);
        for (size_t k = 0; k < n; k++) {
            size_t src = perm[k];
            for (size_t l = 0; l < n; l++) {
                V_re[k * n + l] = A_evec_re[src * n + l];
                V_im[k * n + l] = A_evec_im[src * n + l];
            }
        }
        result = direct_build_complex_eigenvector_list(V_re, V_im, n);
        free(V_re); free(V_im);
    } else {
        result = direct_build_complex_eigenvalue_list(A_eval_re, A_eval_im,
                                                        n, perm);
    }

    free(A_eval_re); free(A_eval_im); free(perm);
    free(A_evec_re); free(A_evec_im);

    return direct_apply_k_spec_list(result, k_spec);
}
#ifdef USE_MPFR

/* MPFR variant of `direct_tridiag_real_sym` -- see that routine for
 * the algorithm.  Workspace mpfr_t* arrays are caller-supplied and
 * pre-initialised; `tmp` holds 8 scratch cells used in the inner
 * loops. */
void direct_tridiag_real_sym_M(mpfr_t* A, size_t n, mpfr_prec_t bits,
                                       mpfr_t* diag, mpfr_t* sub,
                                       mpfr_t* Q, bool want_Q,
                                       mpfr_t* u, mpfr_t* p, mpfr_t* q,
                                       mpfr_t* tmp /* >= 8 cells */) {
    (void)bits;
    mpfr_t* sigma   = &tmp[0];
    mpfr_t* xk1     = &tmp[1];
    mpfr_t* norm_x  = &tmp[2];
    mpfr_t* alpha   = &tmp[3];
    mpfr_t* unorm2  = &tmp[4];
    mpfr_t* unorm   = &tmp[5];
    mpfr_t* K       = &tmp[6];
    mpfr_t* s_acc   = &tmp[7];

    if (want_Q) {
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                mpfr_set_si(Q[i * n + j], (i == j) ? 1 : 0, MPFR_RNDN);
    }

    for (size_t k = 0; k + 2 < n; k++) {
        /* sigma = sum_{i=k+1..n-1} A[i,k]^2 */
        mpfr_set_zero(*sigma, 1);
        for (size_t i = k + 1; i < n; i++) {
            mpfr_mul(*s_acc, A[i * n + k], A[i * n + k], MPFR_RNDN);
            mpfr_add(*sigma, *sigma, *s_acc, MPFR_RNDN);
        }
        mpfr_set(*xk1, A[(k + 1) * n + k], MPFR_RNDN);
        if (mpfr_zero_p(*sigma)) {
            mpfr_set(sub[k], *xk1, MPFR_RNDN);
            continue;
        }
        mpfr_sqrt(*norm_x, *sigma, MPFR_RNDN);
        /* alpha = -sign(xk1) * ||x||  (cancellation-safe) */
        if (mpfr_sgn(*xk1) >= 0) mpfr_neg(*alpha, *norm_x, MPFR_RNDN);
        else                     mpfr_set(*alpha, *norm_x, MPFR_RNDN);

        /* u[k+1] = xk1 - alpha;  u[i>k+1] = A[i,k] */
        mpfr_sub(u[k + 1], *xk1, *alpha, MPFR_RNDN);
        for (size_t i = k + 2; i < n; i++) mpfr_set(u[i], A[i * n + k], MPFR_RNDN);

        /* unorm2 = sum u^2 */
        mpfr_mul(*unorm2, u[k + 1], u[k + 1], MPFR_RNDN);
        for (size_t i = k + 2; i < n; i++) {
            mpfr_mul(*s_acc, u[i], u[i], MPFR_RNDN);
            mpfr_add(*unorm2, *unorm2, *s_acc, MPFR_RNDN);
        }
        if (mpfr_zero_p(*unorm2)) {
            mpfr_set(sub[k], *xk1, MPFR_RNDN);
            continue;
        }
        mpfr_sqrt(*unorm, *unorm2, MPFR_RNDN);
        for (size_t i = k + 1; i < n; i++)
            mpfr_div(u[i], u[i], *unorm, MPFR_RNDN);

        /* p[i] = 2 * sum_j A[i,j] u[j]   on the trailing block */
        for (size_t i = k + 1; i < n; i++) {
            mpfr_set_zero(*s_acc, 1);
            for (size_t j = k + 1; j < n; j++) {
                mpfr_t prod;
                mpfr_init2(prod, mpfr_get_prec(*s_acc));
                mpfr_mul(prod, A[i * n + j], u[j], MPFR_RNDN);
                mpfr_add(*s_acc, *s_acc, prod, MPFR_RNDN);
                mpfr_clear(prod);
            }
            mpfr_mul_2si(p[i], *s_acc, 1, MPFR_RNDN);   /* p[i] = 2 * s_acc */
        }
        /* K = (u^T p) / 2 */
        mpfr_set_zero(*K, 1);
        for (size_t i = k + 1; i < n; i++) {
            mpfr_mul(*s_acc, u[i], p[i], MPFR_RNDN);
            mpfr_add(*K, *K, *s_acc, MPFR_RNDN);
        }
        mpfr_div_2si(*K, *K, 1, MPFR_RNDN);
        /* q[i] = p[i] - 2 K u[i] */
        for (size_t i = k + 1; i < n; i++) {
            mpfr_mul(*s_acc, *K, u[i], MPFR_RNDN);
            mpfr_mul_2si(*s_acc, *s_acc, 1, MPFR_RNDN);
            mpfr_sub(q[i], p[i], *s_acc, MPFR_RNDN);
        }
        /* A_22 -= u q^T + q u^T */
        for (size_t i = k + 1; i < n; i++) {
            for (size_t j = k + 1; j < n; j++) {
                mpfr_mul(*s_acc, u[i], q[j], MPFR_RNDN);
                mpfr_sub(A[i * n + j], A[i * n + j], *s_acc, MPFR_RNDN);
                mpfr_mul(*s_acc, q[i], u[j], MPFR_RNDN);
                mpfr_sub(A[i * n + j], A[i * n + j], *s_acc, MPFR_RNDN);
            }
        }
        /* Set subdiag + clear eliminated entries (drift control). */
        mpfr_set(sub[k], *alpha, MPFR_RNDN);
        mpfr_set(A[(k + 1) * n + k], *alpha, MPFR_RNDN);
        mpfr_set(A[k * n + (k + 1)], *alpha, MPFR_RNDN);
        for (size_t i = k + 2; i < n; i++) {
            mpfr_set_zero(A[i * n + k], 1);
            mpfr_set_zero(A[k * n + i], 1);
        }

        if (want_Q) {
            for (size_t i = 0; i < n; i++) {
                mpfr_set_zero(*s_acc, 1);
                for (size_t j = k + 1; j < n; j++) {
                    mpfr_t prod;
                    mpfr_init2(prod, mpfr_get_prec(*s_acc));
                    mpfr_mul(prod, Q[i * n + j], u[j], MPFR_RNDN);
                    mpfr_add(*s_acc, *s_acc, prod, MPFR_RNDN);
                    mpfr_clear(prod);
                }
                mpfr_mul_2si(*s_acc, *s_acc, 1, MPFR_RNDN);  /* 2 * (Q_row . u) */
                for (size_t j = k + 1; j < n; j++) {
                    mpfr_t prod;
                    mpfr_init2(prod, mpfr_get_prec(*s_acc));
                    mpfr_mul(prod, *s_acc, u[j], MPFR_RNDN);
                    mpfr_sub(Q[i * n + j], Q[i * n + j], prod, MPFR_RNDN);
                    mpfr_clear(prod);
                }
            }
        }
    }

    for (size_t i = 0; i < n; i++) mpfr_set(diag[i], A[i * n + i], MPFR_RNDN);
    if (n >= 2) mpfr_set(sub[n - 2], A[(n - 1) * n + (n - 2)], MPFR_RNDN);
}

/* MPFR variant of `direct_symtridiag_qr`.  `tmp` provides 12+ scratch
 * cells.  rel_tol scales as 2^{-bits+3}: a few ULPs at the working
 * precision, mirroring the 1e-14 ~= 2^-46 (for bits=53) used in the
 * machine version. */
int direct_symtridiag_qr_M(mpfr_t* diag, mpfr_t* sub, size_t n,
                                    mpfr_prec_t bits,
                                    mpfr_t* Q, bool want_Q,
                                    mpfr_t* tmp /* >= 12 cells */) {
    if (n == 0) return 0;
    const size_t max_sweeps = 30 * n;
    size_t sweeps = 0;

    mpfr_t* tol     = &tmp[0];
    mpfr_t* d       = &tmp[1];
    mpfr_t* e       = &tmp[2];
    mpfr_t* t       = &tmp[3];
    mpfr_t* mu      = &tmp[4];
    mpfr_t* x       = &tmp[5];
    mpfr_t* z       = &tmp[6];
    mpfr_t* c       = &tmp[7];
    mpfr_t* s       = &tmp[8];
    mpfr_t* r       = &tmp[9];
    mpfr_t* scratch1= &tmp[10];
    mpfr_t* scratch2= &tmp[11];

    /* Relative tolerance: ~ 8 * 2^-bits == a few ULPs. */
    mpfr_t rel_tol;
    mpfr_init2(rel_tol, bits);
    mpfr_set_ui(rel_tol, 1, MPFR_RNDN);
    mpfr_div_2si(rel_tol, rel_tol, (long)bits - 3, MPFR_RNDN);

    size_t end = n;
    while (end > 1) {
        /* Find largest m s.t. sub[m..end-2] are all significant. */
        size_t m = end - 1;
        while (m > 0) {
            mpfr_abs(*scratch1, diag[m - 1], MPFR_RNDN);
            mpfr_abs(*scratch2, diag[m], MPFR_RNDN);
            mpfr_add(*scratch1, *scratch1, *scratch2, MPFR_RNDN);
            mpfr_mul(*tol, rel_tol, *scratch1, MPFR_RNDN);
            mpfr_abs(*scratch1, sub[m - 1], MPFR_RNDN);
            if (mpfr_cmp(*scratch1, *tol) <= 0) {
                mpfr_set_zero(sub[m - 1], 1);
                break;
            }
            m--;
        }
        if (m == end - 1) { end--; continue; }

        if (++sweeps > max_sweeps) {
            mpfr_clear(rel_tol);
            return -1;
        }

        /* Wilkinson shift on trailing 2x2 block. */
        mpfr_sub(*d, diag[end - 2], diag[end - 1], MPFR_RNDN);
        mpfr_div_2si(*d, *d, 1, MPFR_RNDN);                   /* d = (d_{e-2} - d_{e-1}) / 2 */
        mpfr_set(*e, sub[end - 2], MPFR_RNDN);
        if (mpfr_zero_p(*d)) {
            mpfr_abs(*t, *e, MPFR_RNDN);
        } else {
            mpfr_hypot(*scratch1, *d, *e, MPFR_RNDN);
            mpfr_abs(*scratch2, *d, MPFR_RNDN);
            mpfr_add(*t, *scratch2, *scratch1, MPFR_RNDN);
        }
        /* mu = d_{e-1} - sign(d) * e^2 / t */
        mpfr_mul(*scratch1, *e, *e, MPFR_RNDN);
        if (mpfr_zero_p(*t)) {
            mpfr_set_zero(*scratch1, 1);
        } else {
            mpfr_div(*scratch1, *scratch1, *t, MPFR_RNDN);
        }
        if (mpfr_sgn(*d) < 0) mpfr_neg(*scratch1, *scratch1, MPFR_RNDN);
        mpfr_sub(*mu, diag[end - 1], *scratch1, MPFR_RNDN);

        /* Implicit QR sweep using Givens rotations on [m..end-1]. */
        mpfr_sub(*x, diag[m], *mu, MPFR_RNDN);
        mpfr_set(*z, sub[m], MPFR_RNDN);
        for (size_t k = m; k < end - 1; k++) {
            mpfr_hypot(*r, *x, *z, MPFR_RNDN);
            if (mpfr_zero_p(*r)) {
                mpfr_set_ui(*c, 1, MPFR_RNDN);
                mpfr_set_zero(*s, 1);
            } else {
                mpfr_div(*c, *x, *r, MPFR_RNDN);
                mpfr_div(*s, *z, *r, MPFR_RNDN);
            }

            if (k > m) mpfr_set(sub[k - 1], *r, MPFR_RNDN);

            /* Snapshot d_k, d_k+1, e_k. */
            mpfr_t d_k, d_k1, e_k;
            mpfr_init2(d_k,  bits); mpfr_set(d_k,  diag[k],     MPFR_RNDN);
            mpfr_init2(d_k1, bits); mpfr_set(d_k1, diag[k + 1], MPFR_RNDN);
            mpfr_init2(e_k,  bits); mpfr_set(e_k,  sub[k],      MPFR_RNDN);

            /* diag[k]   = c^2 d_k + 2 c s e_k + s^2 d_k1 */
            mpfr_mul(*scratch1, *c, *c, MPFR_RNDN);
            mpfr_mul(*scratch1, *scratch1, d_k, MPFR_RNDN);
            mpfr_mul(*scratch2, *c, *s, MPFR_RNDN);
            mpfr_mul(*scratch2, *scratch2, e_k, MPFR_RNDN);
            mpfr_mul_2si(*scratch2, *scratch2, 1, MPFR_RNDN);
            mpfr_add(diag[k], *scratch1, *scratch2, MPFR_RNDN);
            mpfr_mul(*scratch1, *s, *s, MPFR_RNDN);
            mpfr_mul(*scratch1, *scratch1, d_k1, MPFR_RNDN);
            mpfr_add(diag[k], diag[k], *scratch1, MPFR_RNDN);

            /* diag[k+1] = s^2 d_k - 2 c s e_k + c^2 d_k1 */
            mpfr_mul(*scratch1, *s, *s, MPFR_RNDN);
            mpfr_mul(*scratch1, *scratch1, d_k, MPFR_RNDN);
            mpfr_mul(*scratch2, *c, *s, MPFR_RNDN);
            mpfr_mul(*scratch2, *scratch2, e_k, MPFR_RNDN);
            mpfr_mul_2si(*scratch2, *scratch2, 1, MPFR_RNDN);
            mpfr_sub(diag[k + 1], *scratch1, *scratch2, MPFR_RNDN);
            mpfr_mul(*scratch1, *c, *c, MPFR_RNDN);
            mpfr_mul(*scratch1, *scratch1, d_k1, MPFR_RNDN);
            mpfr_add(diag[k + 1], diag[k + 1], *scratch1, MPFR_RNDN);

            /* sub[k] = c s (d_k1 - d_k) + (c^2 - s^2) e_k */
            mpfr_sub(*scratch1, d_k1, d_k, MPFR_RNDN);
            mpfr_mul(*scratch1, *scratch1, *c, MPFR_RNDN);
            mpfr_mul(*scratch1, *scratch1, *s, MPFR_RNDN);
            mpfr_mul(*scratch2, *c, *c, MPFR_RNDN);
            mpfr_mul(d_k, *s, *s, MPFR_RNDN);                 /* reusing d_k as throwaway */
            mpfr_sub(*scratch2, *scratch2, d_k, MPFR_RNDN);
            mpfr_mul(*scratch2, *scratch2, e_k, MPFR_RNDN);
            mpfr_add(sub[k], *scratch1, *scratch2, MPFR_RNDN);

            mpfr_clear(d_k); mpfr_clear(d_k1); mpfr_clear(e_k);

            /* Bulge chase: next x, z. */
            if (k + 1 < end - 1) {
                mpfr_t t_next;
                mpfr_init2(t_next, bits);
                mpfr_set(t_next, sub[k + 1], MPFR_RNDN);
                mpfr_set(*x, sub[k], MPFR_RNDN);
                mpfr_mul(*z, *s, t_next, MPFR_RNDN);
                mpfr_mul(sub[k + 1], *c, t_next, MPFR_RNDN);
                mpfr_clear(t_next);
            }

            /* Q <- Q * Givens (post-multiply cols k, k+1). */
            if (want_Q) {
                for (size_t i = 0; i < n; i++) {
                    mpfr_t qk, qk1;
                    mpfr_init2(qk,  bits); mpfr_set(qk,  Q[i * n + k],     MPFR_RNDN);
                    mpfr_init2(qk1, bits); mpfr_set(qk1, Q[i * n + (k + 1)], MPFR_RNDN);
                    mpfr_mul(*scratch1, *c, qk,  MPFR_RNDN);
                    mpfr_mul(*scratch2, *s, qk1, MPFR_RNDN);
                    mpfr_add(Q[i * n + k], *scratch1, *scratch2, MPFR_RNDN);
                    mpfr_mul(*scratch1, *s, qk,  MPFR_RNDN);
                    mpfr_mul(*scratch2, *c, qk1, MPFR_RNDN);
                    mpfr_sub(Q[i * n + (k + 1)], *scratch2, *scratch1, MPFR_RNDN);
                    mpfr_clear(qk); mpfr_clear(qk1);
                }
            }
        }
    }
    mpfr_clear(rel_tol);
    return 0;
}

/* Sort permutation by descending |vals[i]|, stable. */
void direct_sort_perm_desc_abs_M(const mpfr_t* vals, size_t n,
                                         size_t* perm) {
    for (size_t i = 0; i < n; i++) perm[i] = i;
    for (size_t i = 1; i < n; i++) {
        size_t cur = perm[i];
        mpfr_t ac;
        mpfr_init2(ac, mpfr_get_prec(vals[0]));
        mpfr_abs(ac, vals[cur], MPFR_RNDN);
        size_t j = i;
        while (j > 0) {
            mpfr_t ap;
            mpfr_init2(ap, mpfr_get_prec(vals[0]));
            mpfr_abs(ap, vals[perm[j - 1]], MPFR_RNDN);
            int cmp = mpfr_cmp(ap, ac);
            mpfr_clear(ap);
            if (cmp > 0 || (cmp == 0 && perm[j - 1] < cur)) break;
            perm[j] = perm[j - 1];
            j--;
        }
        perm[j] = cur;
        mpfr_clear(ac);
    }
}

/* Build a List[MPFR, ...] of eigenvalues from a sorted permutation. */
Expr* direct_build_real_eigenvalue_list_M(const mpfr_t* vals, size_t n,
                                                  const size_t* perm) {
    Expr** items = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        items[i] = expr_new_mpfr_copy(vals[perm[i]]);
    }
    Expr* out = expr_new_function(expr_new_symbol("List"), items, n);
    free(items);
    return out;
}

/* Build the eigenvector list (rows of MPFR n-vectors). */
Expr* direct_build_real_eigenvector_list_M(const mpfr_t* Q, size_t n,
                                                   const size_t* perm) {
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t k = 0; k < n; k++) {
        size_t col = perm[k];
        Expr** comps = (Expr**)malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            comps[i] = expr_new_mpfr_copy(Q[i * n + col]);
        }
        rows[k] = expr_new_function(expr_new_symbol("List"), comps, n);
        free(comps);
    }
    Expr* out = expr_new_function(expr_new_symbol("List"), rows, n);
    free(rows);
    return out;
}

/* Real symmetric MPFR entry point. */
static Expr* direct_real_sym_mpfr(const MatM* A, MateigenWant want,
                                    Expr* k_spec) {
    size_t n = A->n;
    mpfr_prec_t bits = A->bits;
    bool want_Q = (want & MATEIGEN_WANT_VECTORS) != 0;

    if (n == 0) {
        Expr* empty = expr_new_function(expr_new_symbol("List"), NULL, 0);
        return direct_apply_k_spec_list(empty, k_spec);
    }

    /* Workspace: copy of A (gets destroyed), tridiag arrays, optional Q,
     * Householder vectors u/p/q, and scratch tmp[12]. */
    mpfr_t* Awork = mpfr_array_alloc(n * n, bits);
    for (size_t i = 0; i < n * n; i++) mpfr_set(Awork[i], A->re[i], MPFR_RNDN);

    mpfr_t* diag = mpfr_array_alloc(n, bits);
    mpfr_t* sub  = mpfr_array_alloc(n > 0 ? n - 1 + 1 : 0, bits);
    /* (allocate n cells so we can address sub[n-2] cleanly even for n==1) */
    mpfr_t* Q    = want_Q ? mpfr_array_alloc(n * n, bits) : NULL;
    mpfr_t* u    = mpfr_array_alloc(n, bits);
    mpfr_t* p    = mpfr_array_alloc(n, bits);
    mpfr_t* q    = mpfr_array_alloc(n, bits);
    mpfr_t* tmp  = mpfr_array_alloc(12, bits);

    /* Special case n == 1: trivial.  Skip tridiag. */
    if (n == 1) {
        mpfr_set(diag[0], Awork[0], MPFR_RNDN);
        if (want_Q) mpfr_set_ui(Q[0], 1, MPFR_RNDN);
    } else {
        direct_tridiag_real_sym_M(Awork, n, bits, diag, sub, Q, want_Q,
                                   u, p, q, tmp);
        direct_symtridiag_qr_M(diag, sub, n, bits, Q, want_Q, tmp);
    }

    size_t* perm = (size_t*)malloc(sizeof(size_t) * n);
    direct_sort_perm_desc_abs_M(diag, n, perm);

    Expr* result;
    if (want_Q) {
        result = direct_build_real_eigenvector_list_M(Q, n, perm);
    } else {
        result = direct_build_real_eigenvalue_list_M(diag, n, perm);
    }
    free(perm);

    mpfr_array_free(Awork, n * n);
    mpfr_array_free(diag, n);
    mpfr_array_free(sub, n > 0 ? n - 1 + 1 : 0);
    if (Q) mpfr_array_free(Q, n * n);
    mpfr_array_free(u, n);
    mpfr_array_free(p, n);
    mpfr_array_free(q, n);
    mpfr_array_free(tmp, 12);

    return direct_apply_k_spec_list(result, k_spec);
}

/* ===================================================================
 * 2d-B: Real general (non-symmetric) MPFR Direct kernel
 *
 * MPFR translation of direct_hessenberg_real / direct_francis_step /
 * direct_split_2x2_real / direct_qr_real_general / schur_eigvec_real /
 * schur_eigvec_complex / schur_compute_eigvecs.  See those routines for
 * algorithmic commentary -- this block is a one-to-one transcription,
 * with each double operation replaced by the corresponding mpfr_*
 * call rounded MPFR_RNDN, and inner-loop scratch lifted to caller-
 * supplied workspace cells where possible.
 * =================================================================== */

/* Householder reduction of A to upper Hessenberg form at MPFR precision.
 *
 * Workspace: `u` (length n), optional `Q` (n*n) for back-transformation,
 * `tmp` provides >= 8 scratch cells re-used inside the function. */
void direct_hessenberg_real_M(mpfr_t* A, size_t n, mpfr_prec_t bits,
                                       mpfr_t* u, mpfr_t* Q,
                                       mpfr_t* tmp /* >= 8 cells */) {
    (void)bits;
    mpfr_t* sigma  = &tmp[0];
    mpfr_t* xk1    = &tmp[1];
    mpfr_t* norm_x = &tmp[2];
    mpfr_t* alpha  = &tmp[3];
    mpfr_t* unorm2 = &tmp[4];
    mpfr_t* unorm  = &tmp[5];
    mpfr_t* s      = &tmp[6];
    mpfr_t* prod   = &tmp[7];

    for (size_t k = 0; k + 2 < n; k++) {
        /* sigma = sum_{i=k+1..n-1} A[i,k]^2 */
        mpfr_set_zero(*sigma, 1);
        for (size_t i = k + 1; i < n; i++) {
            mpfr_mul(*prod, A[i * n + k], A[i * n + k], MPFR_RNDN);
            mpfr_add(*sigma, *sigma, *prod, MPFR_RNDN);
        }
        if (mpfr_zero_p(*sigma)) continue;
        mpfr_set(*xk1, A[(k + 1) * n + k], MPFR_RNDN);
        mpfr_sqrt(*norm_x, *sigma, MPFR_RNDN);
        /* alpha = -sign(xk1) * ||x|| */
        if (mpfr_sgn(*xk1) >= 0) mpfr_neg(*alpha, *norm_x, MPFR_RNDN);
        else                     mpfr_set(*alpha, *norm_x, MPFR_RNDN);

        mpfr_sub(u[k + 1], *xk1, *alpha, MPFR_RNDN);
        for (size_t i = k + 2; i < n; i++) mpfr_set(u[i], A[i * n + k], MPFR_RNDN);
        mpfr_mul(*unorm2, u[k + 1], u[k + 1], MPFR_RNDN);
        for (size_t i = k + 2; i < n; i++) {
            mpfr_mul(*prod, u[i], u[i], MPFR_RNDN);
            mpfr_add(*unorm2, *unorm2, *prod, MPFR_RNDN);
        }
        if (mpfr_zero_p(*unorm2)) continue;
        mpfr_sqrt(*unorm, *unorm2, MPFR_RNDN);
        for (size_t i = k + 1; i < n; i++)
            mpfr_div(u[i], u[i], *unorm, MPFR_RNDN);

        /* Left: A[k+1..n-1, :] <- H * A[k+1..n-1, :] */
        for (size_t j = 0; j < n; j++) {
            mpfr_set_zero(*s, 1);
            for (size_t i = k + 1; i < n; i++) {
                mpfr_mul(*prod, u[i], A[i * n + j], MPFR_RNDN);
                mpfr_add(*s, *s, *prod, MPFR_RNDN);
            }
            mpfr_mul_2si(*s, *s, 1, MPFR_RNDN);
            for (size_t i = k + 1; i < n; i++) {
                mpfr_mul(*prod, *s, u[i], MPFR_RNDN);
                mpfr_sub(A[i * n + j], A[i * n + j], *prod, MPFR_RNDN);
            }
        }
        /* Right: A[:, k+1..n-1] <- A[:, k+1..n-1] * H */
        for (size_t i = 0; i < n; i++) {
            mpfr_set_zero(*s, 1);
            for (size_t j = k + 1; j < n; j++) {
                mpfr_mul(*prod, A[i * n + j], u[j], MPFR_RNDN);
                mpfr_add(*s, *s, *prod, MPFR_RNDN);
            }
            mpfr_mul_2si(*s, *s, 1, MPFR_RNDN);
            for (size_t j = k + 1; j < n; j++) {
                mpfr_mul(*prod, *s, u[j], MPFR_RNDN);
                mpfr_sub(A[i * n + j], A[i * n + j], *prod, MPFR_RNDN);
            }
        }
        /* Q <- Q * H_k */
        if (Q) {
            for (size_t i = 0; i < n; i++) {
                mpfr_set_zero(*s, 1);
                for (size_t j = k + 1; j < n; j++) {
                    mpfr_mul(*prod, Q[i * n + j], u[j], MPFR_RNDN);
                    mpfr_add(*s, *s, *prod, MPFR_RNDN);
                }
                mpfr_mul_2si(*s, *s, 1, MPFR_RNDN);
                for (size_t j = k + 1; j < n; j++) {
                    mpfr_mul(*prod, *s, u[j], MPFR_RNDN);
                    mpfr_sub(Q[i * n + j], Q[i * n + j], *prod, MPFR_RNDN);
                }
            }
        }

        mpfr_set(A[(k + 1) * n + k], *alpha, MPFR_RNDN);
        for (size_t i = k + 2; i < n; i++) mpfr_set_zero(A[i * n + k], 1);
    }
}

/* One Francis double-shift QR step on the active block H[q..p-1].
 * MPFR translation of direct_francis_step; `tmp` provides >= 14 cells.
 *
 * `exceptional`: when true, override the trailing-2x2 shift with an ad-
 * hoc shift derived from subdiagonal magnitudes.  This breaks QR
 * stalling cycles that arise at high precision when the standard shift
 * leaves the same active block unchanged after a sweep.  (See Stewart
 * "Matrix Algorithms II" 2.5.1; LAPACK dlahqr uses the same trick.) */
static void direct_francis_step_M(mpfr_t* H, size_t n, mpfr_prec_t bits,
                                    size_t q, size_t p, mpfr_t* Q,
                                    bool exceptional,
                                    mpfr_t* tmp /* >= 14 cells */) {
    mpfr_t* s_tr   = &tmp[0];   /* trace of trailing 2x2 */
    mpfr_t* t_det  = &tmp[1];   /* det of trailing 2x2   */
    mpfr_t* x      = &tmp[2];
    mpfr_t* y      = &tmp[3];
    mpfr_t* z      = &tmp[4];
    mpfr_t* sig    = &tmp[5];
    mpfr_t* nv     = &tmp[6];   /* norm_v */
    mpfr_t* u1     = &tmp[7];
    mpfr_t* u2     = &tmp[8];
    mpfr_t* u3     = &tmp[9];
    mpfr_t* un2    = &tmp[10];
    mpfr_t* un     = &tmp[11];
    mpfr_t* s0     = &tmp[12];
    mpfr_t* prod   = &tmp[13];

    /* Shift derived from trailing 2x2 block (or ad-hoc for exceptional). */
    if (!exceptional) {
        mpfr_t h11, h12, h21, h22;
        mpfr_init2(h11, bits); mpfr_set(h11, H[(p - 2) * n + (p - 2)], MPFR_RNDN);
        mpfr_init2(h12, bits); mpfr_set(h12, H[(p - 2) * n + (p - 1)], MPFR_RNDN);
        mpfr_init2(h21, bits); mpfr_set(h21, H[(p - 1) * n + (p - 2)], MPFR_RNDN);
        mpfr_init2(h22, bits); mpfr_set(h22, H[(p - 1) * n + (p - 1)], MPFR_RNDN);
        mpfr_add(*s_tr, h11, h22, MPFR_RNDN);
        mpfr_mul(*t_det, h11, h22, MPFR_RNDN);
        mpfr_mul(*prod,  h12, h21, MPFR_RNDN);
        mpfr_sub(*t_det, *t_det, *prod, MPFR_RNDN);
        mpfr_clear(h11); mpfr_clear(h12); mpfr_clear(h21); mpfr_clear(h22);
    } else {
        /* Exceptional shift: s_tr = 1.5 * (|H[p-1, p-2]| + |H[p-2, p-3]|),
         * t_det = s_tr^2 / 4 (so disc = 0, a real double shift centred
         * at s_tr/2).  Matches LAPACK dlahqr's "exceptional shift" path. */
        mpfr_t a1, a2;
        mpfr_init2(a1, bits); mpfr_init2(a2, bits);
        mpfr_abs(a1, H[(p - 1) * n + (p - 2)], MPFR_RNDN);
        if (p >= 3) mpfr_abs(a2, H[(p - 2) * n + (p - 3)], MPFR_RNDN);
        else        mpfr_set_zero(a2, 1);
        mpfr_add(*s_tr, a1, a2, MPFR_RNDN);
        mpfr_mul_d(*s_tr, *s_tr, 1.5, MPFR_RNDN);
        mpfr_mul(*t_det, *s_tr, *s_tr, MPFR_RNDN);
        mpfr_div_2si(*t_det, *t_det, 2, MPFR_RNDN);
        mpfr_clear(a1); mpfr_clear(a2);
    }

    /* First three entries of M's first column for the active block. */
    {
        mpfr_t g11, g12, g21, g22, g32;
        mpfr_init2(g11, bits); mpfr_set(g11, H[q * n + q],             MPFR_RNDN);
        mpfr_init2(g12, bits); mpfr_set(g12, H[q * n + (q + 1)],       MPFR_RNDN);
        mpfr_init2(g21, bits); mpfr_set(g21, H[(q + 1) * n + q],       MPFR_RNDN);
        mpfr_init2(g22, bits); mpfr_set(g22, H[(q + 1) * n + (q + 1)], MPFR_RNDN);
        mpfr_init2(g32, bits); mpfr_set(g32, H[(q + 2) * n + (q + 1)], MPFR_RNDN);

        /* x = g11^2 + g12 g21 - s_tr g11 + t_det */
        mpfr_mul(*x, g11, g11, MPFR_RNDN);
        mpfr_mul(*prod, g12, g21, MPFR_RNDN);
        mpfr_add(*x, *x, *prod, MPFR_RNDN);
        mpfr_mul(*prod, *s_tr, g11, MPFR_RNDN);
        mpfr_sub(*x, *x, *prod, MPFR_RNDN);
        mpfr_add(*x, *x, *t_det, MPFR_RNDN);
        /* y = g21 (g11 + g22 - s_tr) */
        mpfr_add(*prod, g11, g22, MPFR_RNDN);
        mpfr_sub(*prod, *prod, *s_tr, MPFR_RNDN);
        mpfr_mul(*y, g21, *prod, MPFR_RNDN);
        /* z = g21 g32 */
        mpfr_mul(*z, g21, g32, MPFR_RNDN);

        mpfr_clear(g11); mpfr_clear(g12); mpfr_clear(g21);
        mpfr_clear(g22); mpfr_clear(g32);
    }

    for (size_t k = q; k + 2 < p; k++) {
        /* sig = x^2 + y^2 + z^2 */
        mpfr_mul(*sig, *x, *x, MPFR_RNDN);
        mpfr_mul(*prod, *y, *y, MPFR_RNDN);
        mpfr_add(*sig, *sig, *prod, MPFR_RNDN);
        mpfr_mul(*prod, *z, *z, MPFR_RNDN);
        mpfr_add(*sig, *sig, *prod, MPFR_RNDN);
        if (mpfr_zero_p(*sig)) {
            mpfr_set(*x, H[(k + 1) * n + k], MPFR_RNDN);
            mpfr_set(*y, H[(k + 2) * n + k], MPFR_RNDN);
            if (k + 3 < p) mpfr_set(*z, H[(k + 3) * n + k], MPFR_RNDN);
            else mpfr_set_zero(*z, 1);
            continue;
        }
        mpfr_sqrt(*nv, *sig, MPFR_RNDN);
        /* a = -sign(x) * nv */
        mpfr_t a;
        mpfr_init2(a, bits);
        if (mpfr_sgn(*x) >= 0) mpfr_neg(a, *nv, MPFR_RNDN);
        else                   mpfr_set(a, *nv, MPFR_RNDN);
        mpfr_sub(*u1, *x, a, MPFR_RNDN);
        mpfr_set(*u2, *y, MPFR_RNDN);
        mpfr_set(*u3, *z, MPFR_RNDN);
        mpfr_clear(a);

        mpfr_mul(*un2, *u1, *u1, MPFR_RNDN);
        mpfr_mul(*prod, *u2, *u2, MPFR_RNDN);
        mpfr_add(*un2, *un2, *prod, MPFR_RNDN);
        mpfr_mul(*prod, *u3, *u3, MPFR_RNDN);
        mpfr_add(*un2, *un2, *prod, MPFR_RNDN);
        mpfr_sqrt(*un, *un2, MPFR_RNDN);
        mpfr_div(*u1, *u1, *un, MPFR_RNDN);
        mpfr_div(*u2, *u2, *un, MPFR_RNDN);
        mpfr_div(*u3, *u3, *un, MPFR_RNDN);

        /* Left: rows (k, k+1, k+2). */
        size_t col_start = (k > q) ? k - 1 : q;
        for (size_t j = col_start; j < n; j++) {
            mpfr_mul(*s0, *u1, H[k * n + j], MPFR_RNDN);
            mpfr_mul(*prod, *u2, H[(k + 1) * n + j], MPFR_RNDN);
            mpfr_add(*s0, *s0, *prod, MPFR_RNDN);
            mpfr_mul(*prod, *u3, H[(k + 2) * n + j], MPFR_RNDN);
            mpfr_add(*s0, *s0, *prod, MPFR_RNDN);
            mpfr_mul_2si(*s0, *s0, 1, MPFR_RNDN);
            mpfr_mul(*prod, *s0, *u1, MPFR_RNDN);
            mpfr_sub(H[k * n + j], H[k * n + j], *prod, MPFR_RNDN);
            mpfr_mul(*prod, *s0, *u2, MPFR_RNDN);
            mpfr_sub(H[(k + 1) * n + j], H[(k + 1) * n + j], *prod, MPFR_RNDN);
            mpfr_mul(*prod, *s0, *u3, MPFR_RNDN);
            mpfr_sub(H[(k + 2) * n + j], H[(k + 2) * n + j], *prod, MPFR_RNDN);
        }
        /* Right: cols (k, k+1, k+2). */
        size_t row_end = (k + 3 < p) ? (k + 3) : (p - 1);
        for (size_t i = 0; i <= row_end; i++) {
            mpfr_mul(*s0, H[i * n + k], *u1, MPFR_RNDN);
            mpfr_mul(*prod, H[i * n + (k + 1)], *u2, MPFR_RNDN);
            mpfr_add(*s0, *s0, *prod, MPFR_RNDN);
            mpfr_mul(*prod, H[i * n + (k + 2)], *u3, MPFR_RNDN);
            mpfr_add(*s0, *s0, *prod, MPFR_RNDN);
            mpfr_mul_2si(*s0, *s0, 1, MPFR_RNDN);
            mpfr_mul(*prod, *s0, *u1, MPFR_RNDN);
            mpfr_sub(H[i * n + k], H[i * n + k], *prod, MPFR_RNDN);
            mpfr_mul(*prod, *s0, *u2, MPFR_RNDN);
            mpfr_sub(H[i * n + (k + 1)], H[i * n + (k + 1)], *prod, MPFR_RNDN);
            mpfr_mul(*prod, *s0, *u3, MPFR_RNDN);
            mpfr_sub(H[i * n + (k + 2)], H[i * n + (k + 2)], *prod, MPFR_RNDN);
        }
        if (Q) {
            for (size_t i = 0; i < n; i++) {
                mpfr_mul(*s0, Q[i * n + k], *u1, MPFR_RNDN);
                mpfr_mul(*prod, Q[i * n + (k + 1)], *u2, MPFR_RNDN);
                mpfr_add(*s0, *s0, *prod, MPFR_RNDN);
                mpfr_mul(*prod, Q[i * n + (k + 2)], *u3, MPFR_RNDN);
                mpfr_add(*s0, *s0, *prod, MPFR_RNDN);
                mpfr_mul_2si(*s0, *s0, 1, MPFR_RNDN);
                mpfr_mul(*prod, *s0, *u1, MPFR_RNDN);
                mpfr_sub(Q[i * n + k], Q[i * n + k], *prod, MPFR_RNDN);
                mpfr_mul(*prod, *s0, *u2, MPFR_RNDN);
                mpfr_sub(Q[i * n + (k + 1)], Q[i * n + (k + 1)], *prod, MPFR_RNDN);
                mpfr_mul(*prod, *s0, *u3, MPFR_RNDN);
                mpfr_sub(Q[i * n + (k + 2)], Q[i * n + (k + 2)], *prod, MPFR_RNDN);
            }
        }

        mpfr_set(*x, H[(k + 1) * n + k], MPFR_RNDN);
        mpfr_set(*y, H[(k + 2) * n + k], MPFR_RNDN);
        if (k + 3 < p) mpfr_set(*z, H[(k + 3) * n + k], MPFR_RNDN);
        else mpfr_set_zero(*z, 1);
    }

    /* Final 2-element Householder on (x, y) at rows (p-2, p-1). */
    {
        size_t k = p - 2;
        mpfr_mul(*sig, *x, *x, MPFR_RNDN);
        mpfr_mul(*prod, *y, *y, MPFR_RNDN);
        mpfr_add(*sig, *sig, *prod, MPFR_RNDN);
        if (mpfr_zero_p(*sig)) return;
        mpfr_sqrt(*nv, *sig, MPFR_RNDN);
        mpfr_t a;
        mpfr_init2(a, bits);
        if (mpfr_sgn(*x) >= 0) mpfr_neg(a, *nv, MPFR_RNDN);
        else                   mpfr_set(a, *nv, MPFR_RNDN);
        mpfr_sub(*u1, *x, a, MPFR_RNDN);
        mpfr_set(*u2, *y, MPFR_RNDN);
        mpfr_clear(a);
        mpfr_mul(*un2, *u1, *u1, MPFR_RNDN);
        mpfr_mul(*prod, *u2, *u2, MPFR_RNDN);
        mpfr_add(*un2, *un2, *prod, MPFR_RNDN);
        mpfr_sqrt(*un, *un2, MPFR_RNDN);
        mpfr_div(*u1, *u1, *un, MPFR_RNDN);
        mpfr_div(*u2, *u2, *un, MPFR_RNDN);

        size_t col_start = k - 1;
        for (size_t j = col_start; j < n; j++) {
            mpfr_mul(*s0, *u1, H[k * n + j], MPFR_RNDN);
            mpfr_mul(*prod, *u2, H[(k + 1) * n + j], MPFR_RNDN);
            mpfr_add(*s0, *s0, *prod, MPFR_RNDN);
            mpfr_mul_2si(*s0, *s0, 1, MPFR_RNDN);
            mpfr_mul(*prod, *s0, *u1, MPFR_RNDN);
            mpfr_sub(H[k * n + j], H[k * n + j], *prod, MPFR_RNDN);
            mpfr_mul(*prod, *s0, *u2, MPFR_RNDN);
            mpfr_sub(H[(k + 1) * n + j], H[(k + 1) * n + j], *prod, MPFR_RNDN);
        }
        for (size_t i = 0; i < p; i++) {
            mpfr_mul(*s0, H[i * n + k], *u1, MPFR_RNDN);
            mpfr_mul(*prod, H[i * n + (k + 1)], *u2, MPFR_RNDN);
            mpfr_add(*s0, *s0, *prod, MPFR_RNDN);
            mpfr_mul_2si(*s0, *s0, 1, MPFR_RNDN);
            mpfr_mul(*prod, *s0, *u1, MPFR_RNDN);
            mpfr_sub(H[i * n + k], H[i * n + k], *prod, MPFR_RNDN);
            mpfr_mul(*prod, *s0, *u2, MPFR_RNDN);
            mpfr_sub(H[i * n + (k + 1)], H[i * n + (k + 1)], *prod, MPFR_RNDN);
        }
        if (Q) {
            for (size_t i = 0; i < n; i++) {
                mpfr_mul(*s0, Q[i * n + k], *u1, MPFR_RNDN);
                mpfr_mul(*prod, Q[i * n + (k + 1)], *u2, MPFR_RNDN);
                mpfr_add(*s0, *s0, *prod, MPFR_RNDN);
                mpfr_mul_2si(*s0, *s0, 1, MPFR_RNDN);
                mpfr_mul(*prod, *s0, *u1, MPFR_RNDN);
                mpfr_sub(Q[i * n + k], Q[i * n + k], *prod, MPFR_RNDN);
                mpfr_mul(*prod, *s0, *u2, MPFR_RNDN);
                mpfr_sub(Q[i * n + (k + 1)], Q[i * n + (k + 1)], *prod, MPFR_RNDN);
            }
        }
    }
}

/* MPFR variant of direct_split_2x2_real.  `tmp` provides >= 10 cells. */
static void direct_split_2x2_real_M(mpfr_t* H, size_t n, mpfr_prec_t bits,
                                      size_t p, mpfr_t* Q,
                                      mpfr_t* tmp /* >= 10 cells */) {
    (void)bits;
    mpfr_t* a    = &tmp[0];
    mpfr_t* b    = &tmp[1];
    mpfr_t* c    = &tmp[2];
    mpfr_t* d    = &tmp[3];
    mpfr_t* tr   = &tmp[4];
    mpfr_t* det  = &tmp[5];
    mpfr_t* disc = &tmp[6];
    mpfr_t* cs   = &tmp[7];
    mpfr_t* sn   = &tmp[8];
    mpfr_t* prod = &tmp[9];

    size_t i = p - 2;
    mpfr_set(*a, H[i * n + i],           MPFR_RNDN);
    mpfr_set(*b, H[i * n + (i + 1)],     MPFR_RNDN);
    mpfr_set(*c, H[(i + 1) * n + i],     MPFR_RNDN);
    mpfr_set(*d, H[(i + 1) * n + (i + 1)], MPFR_RNDN);
    mpfr_add(*tr,  *a, *d, MPFR_RNDN);
    mpfr_mul(*det, *a, *d, MPFR_RNDN);
    mpfr_mul(*prod, *b, *c, MPFR_RNDN);
    mpfr_sub(*det, *det, *prod, MPFR_RNDN);
    mpfr_mul(*disc, *tr, *tr, MPFR_RNDN);
    mpfr_mul_2si(*prod, *det, 2, MPFR_RNDN);    /* 4 * det */
    mpfr_sub(*disc, *disc, *prod, MPFR_RNDN);
    if (mpfr_sgn(*disc) < 0) return;

    mpfr_t sq, lam, v0, v1, r;
    mpfr_init2(sq,  bits);
    mpfr_init2(lam, bits);
    mpfr_init2(v0,  bits);
    mpfr_init2(v1,  bits);
    mpfr_init2(r,   bits);

    mpfr_sqrt(sq, *disc, MPFR_RNDN);
    mpfr_add(lam, *tr, sq, MPFR_RNDN);
    mpfr_div_2si(lam, lam, 1, MPFR_RNDN);
    mpfr_sub(v0, lam, *d, MPFR_RNDN);
    mpfr_set(v1, *c, MPFR_RNDN);
    mpfr_hypot(r, v0, v1, MPFR_RNDN);
    if (mpfr_zero_p(r)) {
        mpfr_clear(sq); mpfr_clear(lam); mpfr_clear(v0); mpfr_clear(v1); mpfr_clear(r);
        return;
    }
    mpfr_div(*cs, v0, r, MPFR_RNDN);
    mpfr_div(*sn, v1, r, MPFR_RNDN);
    mpfr_clear(sq); mpfr_clear(lam); mpfr_clear(v0); mpfr_clear(v1); mpfr_clear(r);

    /* Left: rows (i, i+1) <- G^T * rows. */
    for (size_t j = 0; j < n; j++) {
        mpfr_t r0, r1;
        mpfr_init2(r0, bits); mpfr_set(r0, H[i * n + j],       MPFR_RNDN);
        mpfr_init2(r1, bits); mpfr_set(r1, H[(i + 1) * n + j], MPFR_RNDN);
        mpfr_mul(*prod, *cs, r0, MPFR_RNDN);
        mpfr_mul(*a, *sn, r1, MPFR_RNDN);
        mpfr_add(H[i * n + j], *prod, *a, MPFR_RNDN);
        mpfr_mul(*prod, *sn, r0, MPFR_RNDN);
        mpfr_mul(*a, *cs, r1, MPFR_RNDN);
        mpfr_sub(H[(i + 1) * n + j], *a, *prod, MPFR_RNDN);
        mpfr_clear(r0); mpfr_clear(r1);
    }
    /* Right: cols (i, i+1) <- cols * G. */
    for (size_t k = 0; k < n; k++) {
        mpfr_t c0, c1;
        mpfr_init2(c0, bits); mpfr_set(c0, H[k * n + i],       MPFR_RNDN);
        mpfr_init2(c1, bits); mpfr_set(c1, H[k * n + (i + 1)], MPFR_RNDN);
        mpfr_mul(*prod, *cs, c0, MPFR_RNDN);
        mpfr_mul(*a, *sn, c1, MPFR_RNDN);
        mpfr_add(H[k * n + i], *prod, *a, MPFR_RNDN);
        mpfr_mul(*prod, *sn, c0, MPFR_RNDN);
        mpfr_mul(*a, *cs, c1, MPFR_RNDN);
        mpfr_sub(H[k * n + (i + 1)], *a, *prod, MPFR_RNDN);
        mpfr_clear(c0); mpfr_clear(c1);
    }
    if (Q) {
        for (size_t k = 0; k < n; k++) {
            mpfr_t c0, c1;
            mpfr_init2(c0, bits); mpfr_set(c0, Q[k * n + i],       MPFR_RNDN);
            mpfr_init2(c1, bits); mpfr_set(c1, Q[k * n + (i + 1)], MPFR_RNDN);
            mpfr_mul(*prod, *cs, c0, MPFR_RNDN);
            mpfr_mul(*a, *sn, c1, MPFR_RNDN);
            mpfr_add(Q[k * n + i], *prod, *a, MPFR_RNDN);
            mpfr_mul(*prod, *sn, c0, MPFR_RNDN);
            mpfr_mul(*a, *cs, c1, MPFR_RNDN);
            mpfr_sub(Q[k * n + (i + 1)], *a, *prod, MPFR_RNDN);
            mpfr_clear(c0); mpfr_clear(c1);
        }
    }
    mpfr_set_zero(H[(i + 1) * n + i], 1);
}

/* MPFR Francis QR sweep driver -- analogous to direct_qr_real_general.
 * Returns 0 on success, -1 if max sweeps exceeded.
 * Caller-supplied `tmp` workspace must have >= 14 cells (Francis step). */
int direct_qr_real_general_M(mpfr_t* H, size_t n, mpfr_prec_t bits,
                                      mpfr_t* eval_re, mpfr_t* eval_im,
                                      mpfr_t* Q, mpfr_t* tmp /* >= 14 cells */) {
    const size_t max_iter = 60 * n + 30;
    size_t iter = 0;
    size_t p = n;

    /* Relative tolerance: ~ 2^{-bits+10} == 1024 ULPs.  Non-symmetric
     * implicit-shift QR accumulates more per-sweep roundoff than the
     * symmetric tridiagonal QR (factor ~ n^2 from bulge chase + back-
     * transformation), so we deflate at a more conservative threshold
     * than the 2^{-bits+3} used in `direct_symtridiag_qr_M`. */
    mpfr_t rel_tol, tol_val, abs_a, abs_b;
    mpfr_init2(rel_tol, bits);
    mpfr_init2(tol_val, bits);
    mpfr_init2(abs_a,   bits);
    mpfr_init2(abs_b,   bits);
    mpfr_set_ui(rel_tol, 1, MPFR_RNDN);
    mpfr_div_2si(rel_tol, rel_tol, (long)bits - 10, MPFR_RNDN);

    while (p > 0) {
        if (p == 1) {
            mpfr_set(eval_re[0], H[0], MPFR_RNDN);
            mpfr_set_zero(eval_im[0], 1);
            break;
        }
        if (iter++ > max_iter) {
            mpfr_clear(rel_tol); mpfr_clear(tol_val);
            mpfr_clear(abs_a); mpfr_clear(abs_b);
            return -1;
        }

        size_t q = p - 1;
        while (q > 0) {
            mpfr_abs(abs_a, H[(q - 1) * n + (q - 1)], MPFR_RNDN);
            mpfr_abs(abs_b, H[q * n + q], MPFR_RNDN);
            mpfr_add(tol_val, abs_a, abs_b, MPFR_RNDN);
            mpfr_mul(tol_val, tol_val, rel_tol, MPFR_RNDN);
            mpfr_abs(abs_a, H[q * n + (q - 1)], MPFR_RNDN);
            if (mpfr_cmp(abs_a, tol_val) <= 0) {
                mpfr_set_zero(H[q * n + (q - 1)], 1);
                break;
            }
            q -= 1;
        }

        if (q == p - 1) {
            mpfr_set(eval_re[p - 1], H[(p - 1) * n + (p - 1)], MPFR_RNDN);
            mpfr_set_zero(eval_im[p - 1], 1);
            p -= 1;
            iter = 0;
            continue;
        }
        if (q == p - 2) {
            mpfr_t a, b, c, d, tr, det, disc, sq, half_tr;
            mpfr_init2(a, bits); mpfr_set(a, H[(p - 2) * n + (p - 2)], MPFR_RNDN);
            mpfr_init2(b, bits); mpfr_set(b, H[(p - 2) * n + (p - 1)], MPFR_RNDN);
            mpfr_init2(c, bits); mpfr_set(c, H[(p - 1) * n + (p - 2)], MPFR_RNDN);
            mpfr_init2(d, bits); mpfr_set(d, H[(p - 1) * n + (p - 1)], MPFR_RNDN);
            mpfr_init2(tr, bits);  mpfr_add(tr, a, d, MPFR_RNDN);
            mpfr_init2(det, bits); mpfr_mul(det, a, d, MPFR_RNDN);
            mpfr_t tmp1; mpfr_init2(tmp1, bits);
            mpfr_mul(tmp1, b, c, MPFR_RNDN);
            mpfr_sub(det, det, tmp1, MPFR_RNDN);
            mpfr_init2(disc, bits);
            mpfr_mul(disc, tr, tr, MPFR_RNDN);
            mpfr_mul_2si(tmp1, det, 2, MPFR_RNDN);
            mpfr_sub(disc, disc, tmp1, MPFR_RNDN);
            mpfr_init2(half_tr, bits);
            mpfr_div_2si(half_tr, tr, 1, MPFR_RNDN);
            mpfr_init2(sq, bits);
            if (mpfr_sgn(disc) < 0) {
                mpfr_neg(tmp1, disc, MPFR_RNDN);
                mpfr_sqrt(sq, tmp1, MPFR_RNDN);
                mpfr_div_2si(sq, sq, 1, MPFR_RNDN);
                mpfr_set(eval_re[p - 2], half_tr, MPFR_RNDN);
                mpfr_set(eval_im[p - 2], sq, MPFR_RNDN);
                mpfr_set(eval_re[p - 1], half_tr, MPFR_RNDN);
                mpfr_neg(eval_im[p - 1], sq, MPFR_RNDN);
                p -= 2;
                iter = 0;
            } else {
                mpfr_sqrt(sq, disc, MPFR_RNDN);
                mpfr_div_2si(sq, sq, 1, MPFR_RNDN);
                mpfr_add(eval_re[p - 2], half_tr, sq, MPFR_RNDN);
                mpfr_set_zero(eval_im[p - 2], 1);
                mpfr_sub(eval_re[p - 1], half_tr, sq, MPFR_RNDN);
                mpfr_set_zero(eval_im[p - 1], 1);
                direct_split_2x2_real_M(H, n, bits, p, Q, tmp);
                mpfr_set(H[(p - 2) * n + (p - 2)], eval_re[p - 2], MPFR_RNDN);
                mpfr_set(H[(p - 1) * n + (p - 1)], eval_re[p - 1], MPFR_RNDN);
                p -= 2;
                iter = 0;
            }
            mpfr_clear(a); mpfr_clear(b); mpfr_clear(c); mpfr_clear(d);
            mpfr_clear(tr); mpfr_clear(det); mpfr_clear(disc); mpfr_clear(sq);
            mpfr_clear(half_tr); mpfr_clear(tmp1);
            continue;
        }
        /* Inject an exceptional shift every 10th iter to break stalls. */
        bool exceptional = (iter > 0 && iter % 10 == 0);
        direct_francis_step_M(H, n, bits, q, p, Q, exceptional, tmp);
    }

    mpfr_clear(rel_tol); mpfr_clear(tol_val);
    mpfr_clear(abs_a); mpfr_clear(abs_b);
    return 0;
}

/* Sort permutation by descending |re + i*im| at MPFR precision. */
void direct_sort_perm_desc_abs_complex_M(const mpfr_t* re,
                                                  const mpfr_t* im,
                                                  size_t n, size_t* perm) {
    if (n == 0) return;
    mpfr_prec_t bits = mpfr_get_prec(re[0]);
    for (size_t i = 0; i < n; i++) perm[i] = i;
    for (size_t i = 1; i < n; i++) {
        size_t cur = perm[i];
        mpfr_t ac;
        mpfr_init2(ac, bits);
        mpfr_hypot(ac, re[cur], im[cur], MPFR_RNDN);
        size_t j = i;
        while (j > 0) {
            mpfr_t ap;
            mpfr_init2(ap, bits);
            mpfr_hypot(ap, re[perm[j - 1]], im[perm[j - 1]], MPFR_RNDN);
            int cmp = mpfr_cmp(ap, ac);
            mpfr_clear(ap);
            if (cmp > 0 || (cmp == 0 && perm[j - 1] < cur)) break;
            perm[j] = perm[j - 1];
            j--;
        }
        perm[j] = cur;
        mpfr_clear(ac);
    }
}

/* Build the eigenvalue list at MPFR precision, mixing MPFR reals and
 * Complex[mpfr, mpfr] entries. */
Expr* direct_build_complex_eigenvalue_list_M(const mpfr_t* re,
                                                      const mpfr_t* im,
                                                      size_t n,
                                                      const size_t* perm) {
    Expr** items = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t i = 0; i < n; i++) {
        size_t idx = perm[i];
        if (mpfr_zero_p(im[idx])) {
            items[i] = expr_new_mpfr_copy(re[idx]);
        } else {
            Expr** comp_args = (Expr**)malloc(sizeof(Expr*) * 2);
            comp_args[0] = expr_new_mpfr_copy(re[idx]);
            comp_args[1] = expr_new_mpfr_copy(im[idx]);
            items[i] = expr_new_function(expr_new_symbol("Complex"),
                                          comp_args, 2);
            free(comp_args);
        }
    }
    Expr* out = expr_new_function(expr_new_symbol("List"), items, n);
    free(items);
    return out;
}

/* MPFR back-substitution for a real eigenvector of the Schur form. */
static void schur_eigvec_real_M(const mpfr_t* H, size_t n, mpfr_prec_t bits,
                                  size_t k, const mpfr_t lambda, mpfr_t* v) {
    for (size_t i = 0; i < n; i++) mpfr_set_zero(v[i], 1);
    mpfr_set_ui(v[k], 1, MPFR_RNDN);
    if (k == 0) return;

    mpfr_t rhs, rhs1, rhs2, a, b, c, d, det, prod, abs_det, tiny;
    mpfr_init2(rhs,  bits); mpfr_init2(rhs1, bits); mpfr_init2(rhs2, bits);
    mpfr_init2(a, bits); mpfr_init2(b, bits); mpfr_init2(c, bits); mpfr_init2(d, bits);
    mpfr_init2(det, bits); mpfr_init2(prod, bits);
    mpfr_init2(abs_det, bits); mpfr_init2(tiny, bits);
    /* tiny = 2^{-bits-32} -- treat smaller dets as defective. */
    mpfr_set_ui(tiny, 1, MPFR_RNDN);
    mpfr_div_2si(tiny, tiny, (long)bits + 32, MPFR_RNDN);

    size_t i = k;
    while (i > 0) {
        i--;
        bool is_2x2 = (i > 0) && !mpfr_zero_p(H[i * n + (i - 1)]);
        if (is_2x2) {
            mpfr_set_zero(rhs1, 1); mpfr_set_zero(rhs2, 1);
            for (size_t j = i + 1; j <= k; j++) {
                mpfr_mul(prod, H[(i - 1) * n + j], v[j], MPFR_RNDN);
                mpfr_add(rhs1, rhs1, prod, MPFR_RNDN);
                mpfr_mul(prod, H[i * n + j], v[j], MPFR_RNDN);
                mpfr_add(rhs2, rhs2, prod, MPFR_RNDN);
            }
            mpfr_neg(rhs1, rhs1, MPFR_RNDN);
            mpfr_neg(rhs2, rhs2, MPFR_RNDN);
            mpfr_sub(a, H[(i - 1) * n + (i - 1)], lambda, MPFR_RNDN);
            mpfr_set(b, H[(i - 1) * n + i], MPFR_RNDN);
            mpfr_set(c, H[i * n + (i - 1)], MPFR_RNDN);
            mpfr_sub(d, H[i * n + i], lambda, MPFR_RNDN);
            mpfr_mul(det, a, d, MPFR_RNDN);
            mpfr_mul(prod, b, c, MPFR_RNDN);
            mpfr_sub(det, det, prod, MPFR_RNDN);
            mpfr_abs(abs_det, det, MPFR_RNDN);
            if (mpfr_cmp(abs_det, tiny) < 0) {
                mpfr_set_zero(v[i - 1], 1); mpfr_set_zero(v[i], 1);
            } else {
                mpfr_mul(prod, d, rhs1, MPFR_RNDN);
                mpfr_mul(abs_det, b, rhs2, MPFR_RNDN);
                mpfr_sub(prod, prod, abs_det, MPFR_RNDN);
                mpfr_div(v[i - 1], prod, det, MPFR_RNDN);
                mpfr_mul(prod, a, rhs2, MPFR_RNDN);
                mpfr_mul(abs_det, c, rhs1, MPFR_RNDN);
                mpfr_sub(prod, prod, abs_det, MPFR_RNDN);
                mpfr_div(v[i], prod, det, MPFR_RNDN);
            }
            i--;
        } else {
            mpfr_set_zero(rhs, 1);
            for (size_t j = i + 1; j <= k; j++) {
                mpfr_mul(prod, H[i * n + j], v[j], MPFR_RNDN);
                mpfr_add(rhs, rhs, prod, MPFR_RNDN);
            }
            mpfr_sub(a, H[i * n + i], lambda, MPFR_RNDN);
            mpfr_abs(abs_det, a, MPFR_RNDN);
            if (mpfr_cmp(abs_det, tiny) < 0) {
                mpfr_set_zero(v[i], 1);
            } else {
                mpfr_neg(rhs, rhs, MPFR_RNDN);
                mpfr_div(v[i], rhs, a, MPFR_RNDN);
            }
        }
    }

    mpfr_clear(rhs); mpfr_clear(rhs1); mpfr_clear(rhs2);
    mpfr_clear(a); mpfr_clear(b); mpfr_clear(c); mpfr_clear(d);
    mpfr_clear(det); mpfr_clear(prod);
    mpfr_clear(abs_det); mpfr_clear(tiny);
}

/* MPFR back-substitution for a complex eigenvector of the Schur form.
 * lambda = aval + i bval, bval > 0; writes v_re / v_im. */
static void schur_eigvec_complex_M(const mpfr_t* H, size_t n, mpfr_prec_t bits,
                                     size_t k,
                                     const mpfr_t aval, const mpfr_t bval,
                                     mpfr_t* v_re, mpfr_t* v_im) {
    for (size_t i = 0; i < n; i++) {
        mpfr_set_zero(v_re[i], 1);
        mpfr_set_zero(v_im[i], 1);
    }

    mpfr_t delta, gamma, tmp1, tmp2;
    mpfr_init2(delta, bits); mpfr_init2(gamma, bits);
    mpfr_init2(tmp1,  bits); mpfr_init2(tmp2,  bits);
    mpfr_set(delta, H[(k + 1) * n + (k + 1)], MPFR_RNDN);
    mpfr_set(gamma, H[(k + 1) * n + k],       MPFR_RNDN);
    if (mpfr_zero_p(gamma)) {
        mpfr_set_ui(v_re[k], 1, MPFR_RNDN);
        mpfr_clear(delta); mpfr_clear(gamma); mpfr_clear(tmp1); mpfr_clear(tmp2);
        return;
    }
    mpfr_sub(tmp1, aval, delta, MPFR_RNDN);
    mpfr_div(v_re[k], tmp1, gamma, MPFR_RNDN);
    mpfr_div(v_im[k], bval, gamma, MPFR_RNDN);
    mpfr_set_ui(v_re[k + 1], 1, MPFR_RNDN);
    mpfr_set_zero(v_im[k + 1], 1);
    if (k == 0) {
        mpfr_clear(delta); mpfr_clear(gamma); mpfr_clear(tmp1); mpfr_clear(tmp2);
        return;
    }

    mpfr_t rhs_u_top, rhs_v_top, rhs_u_bot, rhs_v_bot;
    mpfr_init2(rhs_u_top, bits); mpfr_init2(rhs_v_top, bits);
    mpfr_init2(rhs_u_bot, bits); mpfr_init2(rhs_v_bot, bits);
    mpfr_t A_re, A_im, B_re, B_im, C_re, C_im, D_re, D_im;
    mpfr_init2(A_re, bits); mpfr_init2(A_im, bits);
    mpfr_init2(B_re, bits); mpfr_init2(B_im, bits);
    mpfr_init2(C_re, bits); mpfr_init2(C_im, bits);
    mpfr_init2(D_re, bits); mpfr_init2(D_im, bits);
    mpfr_t det_re, det_im, det_mag2, num_re, num_im, r1_re, r1_im, r2_re, r2_im;
    mpfr_init2(det_re, bits); mpfr_init2(det_im, bits); mpfr_init2(det_mag2, bits);
    mpfr_init2(num_re, bits); mpfr_init2(num_im, bits);
    mpfr_init2(r1_re, bits); mpfr_init2(r1_im, bits);
    mpfr_init2(r2_re, bits); mpfr_init2(r2_im, bits);
    mpfr_t tiny;
    mpfr_init2(tiny, bits);
    mpfr_set_ui(tiny, 1, MPFR_RNDN);
    mpfr_div_2si(tiny, tiny, (long)bits + 32, MPFR_RNDN);

    size_t i = k;
    while (i > 0) {
        i--;
        bool is_2x2 = (i > 0) && !mpfr_zero_p(H[i * n + (i - 1)]);
        if (is_2x2) {
            mpfr_set_zero(rhs_u_top, 1); mpfr_set_zero(rhs_v_top, 1);
            mpfr_set_zero(rhs_u_bot, 1); mpfr_set_zero(rhs_v_bot, 1);
            for (size_t j = i + 1; j <= k + 1; j++) {
                mpfr_mul(tmp1, H[(i - 1) * n + j], v_re[j], MPFR_RNDN);
                mpfr_add(rhs_u_top, rhs_u_top, tmp1, MPFR_RNDN);
                mpfr_mul(tmp1, H[(i - 1) * n + j], v_im[j], MPFR_RNDN);
                mpfr_add(rhs_v_top, rhs_v_top, tmp1, MPFR_RNDN);
                mpfr_mul(tmp1, H[i * n + j], v_re[j], MPFR_RNDN);
                mpfr_add(rhs_u_bot, rhs_u_bot, tmp1, MPFR_RNDN);
                mpfr_mul(tmp1, H[i * n + j], v_im[j], MPFR_RNDN);
                mpfr_add(rhs_v_bot, rhs_v_bot, tmp1, MPFR_RNDN);
            }
            mpfr_sub(A_re, H[(i - 1) * n + (i - 1)], aval, MPFR_RNDN);
            mpfr_neg(A_im, bval, MPFR_RNDN);
            mpfr_set(B_re, H[(i - 1) * n + i], MPFR_RNDN);
            mpfr_set_zero(B_im, 1);
            mpfr_set(C_re, H[i * n + (i - 1)], MPFR_RNDN);
            mpfr_set_zero(C_im, 1);
            mpfr_sub(D_re, H[i * n + i], aval, MPFR_RNDN);
            mpfr_neg(D_im, bval, MPFR_RNDN);

            /* det = A D - B C  (complex) */
            mpfr_mul(tmp1, A_re, D_re, MPFR_RNDN);
            mpfr_mul(tmp2, A_im, D_im, MPFR_RNDN);
            mpfr_sub(det_re, tmp1, tmp2, MPFR_RNDN);
            mpfr_mul(tmp1, B_re, C_re, MPFR_RNDN);
            mpfr_mul(tmp2, B_im, C_im, MPFR_RNDN);
            mpfr_sub(tmp1, tmp1, tmp2, MPFR_RNDN);
            mpfr_sub(det_re, det_re, tmp1, MPFR_RNDN);
            mpfr_mul(tmp1, A_re, D_im, MPFR_RNDN);
            mpfr_mul(tmp2, A_im, D_re, MPFR_RNDN);
            mpfr_add(det_im, tmp1, tmp2, MPFR_RNDN);
            mpfr_mul(tmp1, B_re, C_im, MPFR_RNDN);
            mpfr_mul(tmp2, B_im, C_re, MPFR_RNDN);
            mpfr_add(tmp1, tmp1, tmp2, MPFR_RNDN);
            mpfr_sub(det_im, det_im, tmp1, MPFR_RNDN);
            mpfr_mul(det_mag2, det_re, det_re, MPFR_RNDN);
            mpfr_mul(tmp1, det_im, det_im, MPFR_RNDN);
            mpfr_add(det_mag2, det_mag2, tmp1, MPFR_RNDN);

            if (mpfr_cmp(det_mag2, tiny) < 0) {
                mpfr_set_zero(v_re[i - 1], 1); mpfr_set_zero(v_im[i - 1], 1);
                mpfr_set_zero(v_re[i], 1);     mpfr_set_zero(v_im[i], 1);
            } else {
                mpfr_neg(r1_re, rhs_u_top, MPFR_RNDN);
                mpfr_neg(r1_im, rhs_v_top, MPFR_RNDN);
                mpfr_neg(r2_re, rhs_u_bot, MPFR_RNDN);
                mpfr_neg(r2_im, rhs_v_bot, MPFR_RNDN);
                /* num1 = D r1 - B r2 (complex) */
                mpfr_mul(tmp1, D_re, r1_re, MPFR_RNDN);
                mpfr_mul(tmp2, D_im, r1_im, MPFR_RNDN);
                mpfr_sub(num_re, tmp1, tmp2, MPFR_RNDN);
                mpfr_mul(tmp1, B_re, r2_re, MPFR_RNDN);
                mpfr_mul(tmp2, B_im, r2_im, MPFR_RNDN);
                mpfr_sub(tmp1, tmp1, tmp2, MPFR_RNDN);
                mpfr_sub(num_re, num_re, tmp1, MPFR_RNDN);
                mpfr_mul(tmp1, D_re, r1_im, MPFR_RNDN);
                mpfr_mul(tmp2, D_im, r1_re, MPFR_RNDN);
                mpfr_add(num_im, tmp1, tmp2, MPFR_RNDN);
                mpfr_mul(tmp1, B_re, r2_im, MPFR_RNDN);
                mpfr_mul(tmp2, B_im, r2_re, MPFR_RNDN);
                mpfr_add(tmp1, tmp1, tmp2, MPFR_RNDN);
                mpfr_sub(num_im, num_im, tmp1, MPFR_RNDN);
                /* v[i-1] = num1 / det */
                mpfr_mul(tmp1, num_re, det_re, MPFR_RNDN);
                mpfr_mul(tmp2, num_im, det_im, MPFR_RNDN);
                mpfr_add(tmp1, tmp1, tmp2, MPFR_RNDN);
                mpfr_div(v_re[i - 1], tmp1, det_mag2, MPFR_RNDN);
                mpfr_mul(tmp1, num_im, det_re, MPFR_RNDN);
                mpfr_mul(tmp2, num_re, det_im, MPFR_RNDN);
                mpfr_sub(tmp1, tmp1, tmp2, MPFR_RNDN);
                mpfr_div(v_im[i - 1], tmp1, det_mag2, MPFR_RNDN);

                /* num2 = A r2 - C r1 (complex) */
                mpfr_mul(tmp1, A_re, r2_re, MPFR_RNDN);
                mpfr_mul(tmp2, A_im, r2_im, MPFR_RNDN);
                mpfr_sub(num_re, tmp1, tmp2, MPFR_RNDN);
                mpfr_mul(tmp1, C_re, r1_re, MPFR_RNDN);
                mpfr_mul(tmp2, C_im, r1_im, MPFR_RNDN);
                mpfr_sub(tmp1, tmp1, tmp2, MPFR_RNDN);
                mpfr_sub(num_re, num_re, tmp1, MPFR_RNDN);
                mpfr_mul(tmp1, A_re, r2_im, MPFR_RNDN);
                mpfr_mul(tmp2, A_im, r2_re, MPFR_RNDN);
                mpfr_add(num_im, tmp1, tmp2, MPFR_RNDN);
                mpfr_mul(tmp1, C_re, r1_im, MPFR_RNDN);
                mpfr_mul(tmp2, C_im, r1_re, MPFR_RNDN);
                mpfr_add(tmp1, tmp1, tmp2, MPFR_RNDN);
                mpfr_sub(num_im, num_im, tmp1, MPFR_RNDN);
                /* v[i] = num2 / det */
                mpfr_mul(tmp1, num_re, det_re, MPFR_RNDN);
                mpfr_mul(tmp2, num_im, det_im, MPFR_RNDN);
                mpfr_add(tmp1, tmp1, tmp2, MPFR_RNDN);
                mpfr_div(v_re[i], tmp1, det_mag2, MPFR_RNDN);
                mpfr_mul(tmp1, num_im, det_re, MPFR_RNDN);
                mpfr_mul(tmp2, num_re, det_im, MPFR_RNDN);
                mpfr_sub(tmp1, tmp1, tmp2, MPFR_RNDN);
                mpfr_div(v_im[i], tmp1, det_mag2, MPFR_RNDN);
            }
            i--;
        } else {
            mpfr_t rhs_u, rhs_v, diag_re, diag_im, denom;
            mpfr_init2(rhs_u, bits); mpfr_init2(rhs_v, bits);
            mpfr_init2(diag_re, bits); mpfr_init2(diag_im, bits);
            mpfr_init2(denom, bits);
            mpfr_set_zero(rhs_u, 1); mpfr_set_zero(rhs_v, 1);
            for (size_t j = i + 1; j <= k + 1; j++) {
                mpfr_mul(tmp1, H[i * n + j], v_re[j], MPFR_RNDN);
                mpfr_add(rhs_u, rhs_u, tmp1, MPFR_RNDN);
                mpfr_mul(tmp1, H[i * n + j], v_im[j], MPFR_RNDN);
                mpfr_add(rhs_v, rhs_v, tmp1, MPFR_RNDN);
            }
            mpfr_sub(diag_re, H[i * n + i], aval, MPFR_RNDN);
            mpfr_neg(diag_im, bval, MPFR_RNDN);
            mpfr_mul(denom, diag_re, diag_re, MPFR_RNDN);
            mpfr_mul(tmp1, diag_im, diag_im, MPFR_RNDN);
            mpfr_add(denom, denom, tmp1, MPFR_RNDN);
            if (mpfr_cmp(denom, tiny) < 0) {
                mpfr_set_zero(v_re[i], 1); mpfr_set_zero(v_im[i], 1);
            } else {
                mpfr_neg(num_re, rhs_u, MPFR_RNDN);
                mpfr_neg(num_im, rhs_v, MPFR_RNDN);
                mpfr_mul(tmp1, num_re, diag_re, MPFR_RNDN);
                mpfr_mul(tmp2, num_im, diag_im, MPFR_RNDN);
                mpfr_add(tmp1, tmp1, tmp2, MPFR_RNDN);
                mpfr_div(v_re[i], tmp1, denom, MPFR_RNDN);
                mpfr_mul(tmp1, num_im, diag_re, MPFR_RNDN);
                mpfr_mul(tmp2, num_re, diag_im, MPFR_RNDN);
                mpfr_sub(tmp1, tmp1, tmp2, MPFR_RNDN);
                mpfr_div(v_im[i], tmp1, denom, MPFR_RNDN);
            }
            mpfr_clear(rhs_u); mpfr_clear(rhs_v);
            mpfr_clear(diag_re); mpfr_clear(diag_im); mpfr_clear(denom);
        }
    }

    mpfr_clear(delta); mpfr_clear(gamma); mpfr_clear(tmp1); mpfr_clear(tmp2);
    mpfr_clear(rhs_u_top); mpfr_clear(rhs_v_top);
    mpfr_clear(rhs_u_bot); mpfr_clear(rhs_v_bot);
    mpfr_clear(A_re); mpfr_clear(A_im);
    mpfr_clear(B_re); mpfr_clear(B_im);
    mpfr_clear(C_re); mpfr_clear(C_im);
    mpfr_clear(D_re); mpfr_clear(D_im);
    mpfr_clear(det_re); mpfr_clear(det_im); mpfr_clear(det_mag2);
    mpfr_clear(num_re); mpfr_clear(num_im);
    mpfr_clear(r1_re); mpfr_clear(r1_im);
    mpfr_clear(r2_re); mpfr_clear(r2_im);
    mpfr_clear(tiny);
}

/* Build eigenvector matrix V from Schur form H and back-transformation Q
 * at MPFR precision.  Result is row-major n*n in V_re, V_im (rows indexed
 * by sorted permutation). */
void schur_compute_eigvecs_M(const mpfr_t* H, const mpfr_t* Q, size_t n,
                                      mpfr_prec_t bits,
                                      const mpfr_t* eval_re, const mpfr_t* eval_im,
                                      const size_t* perm,
                                      mpfr_t* V_re, mpfr_t* V_im) {
    mpfr_t* v_schur_re = mpfr_array_alloc(n, bits);
    mpfr_t* v_schur_im = mpfr_array_alloc(n, bits);
    mpfr_t* w_re       = mpfr_array_alloc(n, bits);
    mpfr_t* w_im       = mpfr_array_alloc(n, bits);

    size_t* inv_perm = (size_t*)malloc(sizeof(size_t) * n);
    for (size_t i = 0; i < n; i++) inv_perm[perm[i]] = i;

    mpfr_t norm2, inv, s_re, s_im, prod, abs_b;
    mpfr_init2(norm2, bits); mpfr_init2(inv,   bits);
    mpfr_init2(s_re,  bits); mpfr_init2(s_im,  bits);
    mpfr_init2(prod,  bits); mpfr_init2(abs_b, bits);

    size_t k = 0;
    while (k < n) {
        if (mpfr_zero_p(eval_im[k])) {
            schur_eigvec_real_M(H, n, bits, k, eval_re[k], v_schur_re);
            for (size_t i = 0; i < n; i++) {
                mpfr_set_zero(s_re, 1);
                for (size_t j = 0; j <= k; j++) {
                    mpfr_mul(prod, Q[i * n + j], v_schur_re[j], MPFR_RNDN);
                    mpfr_add(s_re, s_re, prod, MPFR_RNDN);
                }
                mpfr_set(w_re[i], s_re, MPFR_RNDN);
            }
            mpfr_set_zero(norm2, 1);
            for (size_t i = 0; i < n; i++) {
                mpfr_mul(prod, w_re[i], w_re[i], MPFR_RNDN);
                mpfr_add(norm2, norm2, prod, MPFR_RNDN);
            }
            if (mpfr_zero_p(norm2)) {
                mpfr_set_ui(inv, 1, MPFR_RNDN);
            } else {
                mpfr_sqrt(inv, norm2, MPFR_RNDN);
                mpfr_ui_div(inv, 1, inv, MPFR_RNDN);
            }
            size_t sp = inv_perm[k];
            for (size_t i = 0; i < n; i++) {
                mpfr_mul(V_re[sp * n + i], w_re[i], inv, MPFR_RNDN);
                mpfr_set_zero(V_im[sp * n + i], 1);
            }
            k++;
        } else {
            mpfr_abs(abs_b, eval_im[k], MPFR_RNDN);
            schur_eigvec_complex_M(H, n, bits, k, eval_re[k], abs_b,
                                     v_schur_re, v_schur_im);
            for (size_t i = 0; i < n; i++) {
                mpfr_set_zero(s_re, 1);
                mpfr_set_zero(s_im, 1);
                for (size_t j = 0; j <= k + 1; j++) {
                    mpfr_mul(prod, Q[i * n + j], v_schur_re[j], MPFR_RNDN);
                    mpfr_add(s_re, s_re, prod, MPFR_RNDN);
                    mpfr_mul(prod, Q[i * n + j], v_schur_im[j], MPFR_RNDN);
                    mpfr_add(s_im, s_im, prod, MPFR_RNDN);
                }
                mpfr_set(w_re[i], s_re, MPFR_RNDN);
                mpfr_set(w_im[i], s_im, MPFR_RNDN);
            }
            mpfr_set_zero(norm2, 1);
            for (size_t i = 0; i < n; i++) {
                mpfr_mul(prod, w_re[i], w_re[i], MPFR_RNDN);
                mpfr_add(norm2, norm2, prod, MPFR_RNDN);
                mpfr_mul(prod, w_im[i], w_im[i], MPFR_RNDN);
                mpfr_add(norm2, norm2, prod, MPFR_RNDN);
            }
            if (mpfr_zero_p(norm2)) {
                mpfr_set_ui(inv, 1, MPFR_RNDN);
            } else {
                mpfr_sqrt(inv, norm2, MPFR_RNDN);
                mpfr_ui_div(inv, 1, inv, MPFR_RNDN);
            }
            size_t sp1 = inv_perm[k];
            size_t sp2 = inv_perm[k + 1];
            for (size_t i = 0; i < n; i++) {
                mpfr_t rr, mm;
                mpfr_init2(rr, bits); mpfr_init2(mm, bits);
                mpfr_mul(rr, w_re[i], inv, MPFR_RNDN);
                mpfr_mul(mm, w_im[i], inv, MPFR_RNDN);
                mpfr_set(V_re[sp1 * n + i], rr, MPFR_RNDN);
                mpfr_set(V_im[sp1 * n + i], mm, MPFR_RNDN);
                mpfr_set(V_re[sp2 * n + i], rr, MPFR_RNDN);
                mpfr_neg(V_im[sp2 * n + i], mm, MPFR_RNDN);
                mpfr_clear(rr); mpfr_clear(mm);
            }
            k += 2;
        }
    }

    mpfr_array_free(v_schur_re, n); mpfr_array_free(v_schur_im, n);
    mpfr_array_free(w_re, n);       mpfr_array_free(w_im, n);
    free(inv_perm);
    mpfr_clear(norm2); mpfr_clear(inv);
    mpfr_clear(s_re); mpfr_clear(s_im);
    mpfr_clear(prod); mpfr_clear(abs_b);
}

/* Emit the eigenvector matrix V_re + i V_im as a List of List of MPFR
 * (or Complex[mpfr, mpfr]). */
static Expr* direct_build_complex_eigenvector_list_M(const mpfr_t* V_re,
                                                       const mpfr_t* V_im,
                                                       size_t n) {
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t k = 0; k < n; k++) {
        Expr** comps = (Expr**)malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            if (mpfr_zero_p(V_im[k * n + i])) {
                comps[i] = expr_new_mpfr_copy(V_re[k * n + i]);
            } else {
                Expr** args = (Expr**)malloc(sizeof(Expr*) * 2);
                args[0] = expr_new_mpfr_copy(V_re[k * n + i]);
                args[1] = expr_new_mpfr_copy(V_im[k * n + i]);
                comps[i] = expr_new_function(expr_new_symbol("Complex"),
                                              args, 2);
                free(args);
            }
        }
        rows[k] = expr_new_function(expr_new_symbol("List"), comps, n);
        free(comps);
    }
    Expr* out = expr_new_function(expr_new_symbol("List"), rows, n);
    free(rows);
    return out;
}

/* Top-level Direct kernel for real non-symmetric MPFR input. */
static Expr* direct_real_general_mpfr(const MatM* A, MateigenWant want,
                                        Expr* k_spec) {
    size_t n = A->n;
    mpfr_prec_t bits = A->bits;
    if (n == 0) return NULL;

    bool want_Q = (want & MATEIGEN_WANT_VECTORS) != 0;

    mpfr_t* H = mpfr_array_alloc(n * n, bits);
    for (size_t i = 0; i < n * n; i++) mpfr_set(H[i], A->re[i], MPFR_RNDN);
    mpfr_t* Q = NULL;
    if (want_Q) {
        Q = mpfr_array_alloc(n * n, bits);
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                mpfr_set_si(Q[i * n + j], (i == j) ? 1 : 0, MPFR_RNDN);
    }

    mpfr_t* u   = mpfr_array_alloc(n, bits);
    mpfr_t* tmp = mpfr_array_alloc(14, bits);

    if (n >= 3) direct_hessenberg_real_M(H, n, bits, u, Q, tmp);

    mpfr_t* eval_re = mpfr_array_alloc(n, bits);
    mpfr_t* eval_im = mpfr_array_alloc(n, bits);
    for (size_t i = 0; i < n; i++) {
        mpfr_set_zero(eval_re[i], 1);
        mpfr_set_zero(eval_im[i], 1);
    }
    int qr_status = direct_qr_real_general_M(H, n, bits, eval_re, eval_im, Q, tmp);

    if (qr_status != 0) {
        mpfr_array_free(H, n * n);
        if (Q) mpfr_array_free(Q, n * n);
        mpfr_array_free(u, n);
        mpfr_array_free(tmp, 14);
        mpfr_array_free(eval_re, n);
        mpfr_array_free(eval_im, n);
        return NULL;
    }

    size_t* perm = (size_t*)malloc(sizeof(size_t) * n);
    direct_sort_perm_desc_abs_complex_M(eval_re, eval_im, n, perm);

    Expr* out;
    if (want_Q) {
        mpfr_t* V_re = mpfr_array_alloc(n * n, bits);
        mpfr_t* V_im = mpfr_array_alloc(n * n, bits);
        schur_compute_eigvecs_M(H, Q, n, bits, eval_re, eval_im, perm, V_re, V_im);
        out = direct_build_complex_eigenvector_list_M(V_re, V_im, n);
        mpfr_array_free(V_re, n * n);
        mpfr_array_free(V_im, n * n);
    } else {
        out = direct_build_complex_eigenvalue_list_M(eval_re, eval_im, n, perm);
    }

    free(perm);
    mpfr_array_free(H, n * n);
    if (Q) mpfr_array_free(Q, n * n);
    mpfr_array_free(u, n);
    mpfr_array_free(tmp, 14);
    mpfr_array_free(eval_re, n);
    mpfr_array_free(eval_im, n);

    return direct_apply_k_spec_list(out, k_spec);
}

/* ===================================================================
 * 2d-C: Complex Hermitian MPFR Direct kernel
 *
 * MPFR translation of direct_tridiag_complex_hermitian /
 * direct_phase_correct_tridiag / direct_compose_complex_Q_real_Z /
 * direct_build_complex_hermitian_eigvec_list / direct_complex_hermitian_machine.
 * The post-tridiag eigenvalue extraction reuses direct_symtridiag_qr_M
 * from 2d-A.
 * =================================================================== */

/* Hermitian Householder tridiagonalisation at MPFR precision.  Workspace
 * arrays (u, v, q) are paired re/im, length n.  Optional Q is paired
 * re/im, length n*n. */
void direct_tridiag_complex_hermitian_M(mpfr_t* A_re, mpfr_t* A_im,
                                                 size_t n, mpfr_prec_t bits,
                                                 mpfr_t* diag,
                                                 mpfr_t* sub_re, mpfr_t* sub_im,
                                                 mpfr_t* Q_re, mpfr_t* Q_im,
                                                 bool want_Q,
                                                 mpfr_t* u_re, mpfr_t* u_im,
                                                 mpfr_t* v_re, mpfr_t* v_im,
                                                 mpfr_t* q_re, mpfr_t* q_im) {
    if (want_Q) {
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                mpfr_set_si(Q_re[i * n + j], (i == j) ? 1 : 0, MPFR_RNDN);
                mpfr_set_zero(Q_im[i * n + j], 1);
            }
        }
    }

    mpfr_t sigma, xk1_re, xk1_im, xk1_abs, norm_x, alpha_re, alpha_im;
    mpfr_t unorm2, inv_unorm, alpha_v;
    mpfr_t s_re, s_im, prod, t_re, t_im;
    mpfr_init2(sigma,    bits); mpfr_init2(xk1_re,   bits);
    mpfr_init2(xk1_im,   bits); mpfr_init2(xk1_abs,  bits);
    mpfr_init2(norm_x,   bits); mpfr_init2(alpha_re, bits);
    mpfr_init2(alpha_im, bits); mpfr_init2(unorm2,   bits);
    mpfr_init2(inv_unorm, bits); mpfr_init2(alpha_v, bits);
    mpfr_init2(s_re, bits); mpfr_init2(s_im, bits);
    mpfr_init2(prod, bits); mpfr_init2(t_re, bits); mpfr_init2(t_im, bits);

    for (size_t k = 0; k + 2 < n; k++) {
        /* sigma = sum_{i=k+1..n-1} |A[i, k]|^2 */
        mpfr_set_zero(sigma, 1);
        for (size_t i = k + 1; i < n; i++) {
            mpfr_mul(prod, A_re[i * n + k], A_re[i * n + k], MPFR_RNDN);
            mpfr_add(sigma, sigma, prod, MPFR_RNDN);
            mpfr_mul(prod, A_im[i * n + k], A_im[i * n + k], MPFR_RNDN);
            mpfr_add(sigma, sigma, prod, MPFR_RNDN);
        }
        mpfr_set(xk1_re, A_re[(k + 1) * n + k], MPFR_RNDN);
        mpfr_set(xk1_im, A_im[(k + 1) * n + k], MPFR_RNDN);
        mpfr_hypot(xk1_abs, xk1_re, xk1_im, MPFR_RNDN);
        if (mpfr_zero_p(sigma)) {
            mpfr_set(sub_re[k], xk1_re, MPFR_RNDN);
            mpfr_set(sub_im[k], xk1_im, MPFR_RNDN);
            continue;
        }
        mpfr_sqrt(norm_x, sigma, MPFR_RNDN);
        /* alpha = -phase(x_{k+1}) * norm_x. */
        if (mpfr_zero_p(xk1_abs)) {
            mpfr_neg(alpha_re, norm_x, MPFR_RNDN);
            mpfr_set_zero(alpha_im, 1);
        } else {
            mpfr_div(prod, xk1_re, xk1_abs, MPFR_RNDN);
            mpfr_mul(alpha_re, prod, norm_x, MPFR_RNDN);
            mpfr_neg(alpha_re, alpha_re, MPFR_RNDN);
            mpfr_div(prod, xk1_im, xk1_abs, MPFR_RNDN);
            mpfr_mul(alpha_im, prod, norm_x, MPFR_RNDN);
            mpfr_neg(alpha_im, alpha_im, MPFR_RNDN);
        }
        /* u = x with u_{k+1} -= alpha. */
        mpfr_sub(u_re[k + 1], xk1_re, alpha_re, MPFR_RNDN);
        mpfr_sub(u_im[k + 1], xk1_im, alpha_im, MPFR_RNDN);
        for (size_t i = k + 2; i < n; i++) {
            mpfr_set(u_re[i], A_re[i * n + k], MPFR_RNDN);
            mpfr_set(u_im[i], A_im[i * n + k], MPFR_RNDN);
        }
        mpfr_set_zero(unorm2, 1);
        for (size_t i = k + 1; i < n; i++) {
            mpfr_mul(prod, u_re[i], u_re[i], MPFR_RNDN);
            mpfr_add(unorm2, unorm2, prod, MPFR_RNDN);
            mpfr_mul(prod, u_im[i], u_im[i], MPFR_RNDN);
            mpfr_add(unorm2, unorm2, prod, MPFR_RNDN);
        }
        if (mpfr_zero_p(unorm2)) {
            mpfr_set(sub_re[k], xk1_re, MPFR_RNDN);
            mpfr_set(sub_im[k], xk1_im, MPFR_RNDN);
            continue;
        }
        mpfr_sqrt(inv_unorm, unorm2, MPFR_RNDN);
        mpfr_ui_div(inv_unorm, 1, inv_unorm, MPFR_RNDN);
        for (size_t i = k + 1; i < n; i++) {
            mpfr_mul(u_re[i], u_re[i], inv_unorm, MPFR_RNDN);
            mpfr_mul(u_im[i], u_im[i], inv_unorm, MPFR_RNDN);
        }
        /* v = A u  (complex; A Hermitian). */
        for (size_t i = k + 1; i < n; i++) {
            mpfr_set_zero(s_re, 1); mpfr_set_zero(s_im, 1);
            for (size_t j = k + 1; j < n; j++) {
                /* s += A_ij * u_j */
                mpfr_mul(prod, A_re[i * n + j], u_re[j], MPFR_RNDN);
                mpfr_add(s_re, s_re, prod, MPFR_RNDN);
                mpfr_mul(prod, A_im[i * n + j], u_im[j], MPFR_RNDN);
                mpfr_sub(s_re, s_re, prod, MPFR_RNDN);
                mpfr_mul(prod, A_re[i * n + j], u_im[j], MPFR_RNDN);
                mpfr_add(s_im, s_im, prod, MPFR_RNDN);
                mpfr_mul(prod, A_im[i * n + j], u_re[j], MPFR_RNDN);
                mpfr_add(s_im, s_im, prod, MPFR_RNDN);
            }
            mpfr_set(v_re[i], s_re, MPFR_RNDN);
            mpfr_set(v_im[i], s_im, MPFR_RNDN);
        }
        /* alpha_v = u^H v  (real). */
        mpfr_set_zero(alpha_v, 1);
        for (size_t i = k + 1; i < n; i++) {
            mpfr_mul(prod, u_re[i], v_re[i], MPFR_RNDN);
            mpfr_add(alpha_v, alpha_v, prod, MPFR_RNDN);
            mpfr_mul(prod, u_im[i], v_im[i], MPFR_RNDN);
            mpfr_add(alpha_v, alpha_v, prod, MPFR_RNDN);
        }
        /* q = 2 v - 2 alpha_v u. */
        for (size_t i = k + 1; i < n; i++) {
            mpfr_mul_2si(prod, v_re[i], 1, MPFR_RNDN);
            mpfr_mul(s_re, alpha_v, u_re[i], MPFR_RNDN);
            mpfr_mul_2si(s_re, s_re, 1, MPFR_RNDN);
            mpfr_sub(q_re[i], prod, s_re, MPFR_RNDN);
            mpfr_mul_2si(prod, v_im[i], 1, MPFR_RNDN);
            mpfr_mul(s_im, alpha_v, u_im[i], MPFR_RNDN);
            mpfr_mul_2si(s_im, s_im, 1, MPFR_RNDN);
            mpfr_sub(q_im[i], prod, s_im, MPFR_RNDN);
        }
        /* A_ij -= u_i conj(q_j) + q_i conj(u_j). */
        for (size_t i = k + 1; i < n; i++) {
            for (size_t j = k + 1; j < n; j++) {
                /* u_i * conj(q_j) = (u_re + i u_im)(q_re - i q_im) */
                mpfr_mul(t_re, u_re[i], q_re[j], MPFR_RNDN);
                mpfr_mul(prod, u_im[i], q_im[j], MPFR_RNDN);
                mpfr_add(t_re, t_re, prod, MPFR_RNDN);
                mpfr_mul(t_im, u_im[i], q_re[j], MPFR_RNDN);
                mpfr_mul(prod, u_re[i], q_im[j], MPFR_RNDN);
                mpfr_sub(t_im, t_im, prod, MPFR_RNDN);
                /* q_i * conj(u_j) = (q_re + i q_im)(u_re - i u_im) */
                mpfr_mul(s_re, q_re[i], u_re[j], MPFR_RNDN);
                mpfr_mul(prod, q_im[i], u_im[j], MPFR_RNDN);
                mpfr_add(s_re, s_re, prod, MPFR_RNDN);
                mpfr_mul(s_im, q_im[i], u_re[j], MPFR_RNDN);
                mpfr_mul(prod, q_re[i], u_im[j], MPFR_RNDN);
                mpfr_sub(s_im, s_im, prod, MPFR_RNDN);
                mpfr_add(prod, t_re, s_re, MPFR_RNDN);
                mpfr_sub(A_re[i * n + j], A_re[i * n + j], prod, MPFR_RNDN);
                mpfr_add(prod, t_im, s_im, MPFR_RNDN);
                mpfr_sub(A_im[i * n + j], A_im[i * n + j], prod, MPFR_RNDN);
            }
        }
        /* Force analytic sub-column / sub-row. */
        mpfr_set(sub_re[k], alpha_re, MPFR_RNDN);
        mpfr_set(sub_im[k], alpha_im, MPFR_RNDN);
        mpfr_set(A_re[(k + 1) * n + k], alpha_re, MPFR_RNDN);
        mpfr_set(A_im[(k + 1) * n + k], alpha_im, MPFR_RNDN);
        mpfr_set(A_re[k * n + (k + 1)], alpha_re, MPFR_RNDN);
        mpfr_neg(A_im[k * n + (k + 1)], alpha_im, MPFR_RNDN);
        for (size_t i = k + 2; i < n; i++) {
            mpfr_set_zero(A_re[i * n + k], 1);
            mpfr_set_zero(A_im[i * n + k], 1);
            mpfr_set_zero(A_re[k * n + i], 1);
            mpfr_set_zero(A_im[k * n + i], 1);
        }
        /* Q <- Q H. */
        if (want_Q) {
            for (size_t i = 0; i < n; i++) {
                /* (Q u)_i restricted to j >= k+1. */
                mpfr_set_zero(s_re, 1); mpfr_set_zero(s_im, 1);
                for (size_t j = k + 1; j < n; j++) {
                    mpfr_mul(prod, Q_re[i * n + j], u_re[j], MPFR_RNDN);
                    mpfr_add(s_re, s_re, prod, MPFR_RNDN);
                    mpfr_mul(prod, Q_im[i * n + j], u_im[j], MPFR_RNDN);
                    mpfr_sub(s_re, s_re, prod, MPFR_RNDN);
                    mpfr_mul(prod, Q_re[i * n + j], u_im[j], MPFR_RNDN);
                    mpfr_add(s_im, s_im, prod, MPFR_RNDN);
                    mpfr_mul(prod, Q_im[i * n + j], u_re[j], MPFR_RNDN);
                    mpfr_add(s_im, s_im, prod, MPFR_RNDN);
                }
                mpfr_mul_2si(s_re, s_re, 1, MPFR_RNDN);
                mpfr_mul_2si(s_im, s_im, 1, MPFR_RNDN);
                /* Q_ij -= s * conj(u_j) = (s_re + i s_im)(u_re - i u_im) */
                for (size_t j = k + 1; j < n; j++) {
                    mpfr_mul(prod, s_re, u_re[j], MPFR_RNDN);
                    mpfr_mul(t_re, s_im, u_im[j], MPFR_RNDN);
                    mpfr_add(prod, prod, t_re, MPFR_RNDN);
                    mpfr_sub(Q_re[i * n + j], Q_re[i * n + j], prod, MPFR_RNDN);
                    mpfr_mul(prod, s_im, u_re[j], MPFR_RNDN);
                    mpfr_mul(t_re, s_re, u_im[j], MPFR_RNDN);
                    mpfr_sub(prod, prod, t_re, MPFR_RNDN);
                    mpfr_sub(Q_im[i * n + j], Q_im[i * n + j], prod, MPFR_RNDN);
                }
            }
        }
    }
    /* Diagonal -- real part. */
    for (size_t i = 0; i < n; i++) mpfr_set(diag[i], A_re[i * n + i], MPFR_RNDN);
    if (n >= 2) {
        mpfr_set(sub_re[n - 2], A_re[(n - 1) * n + (n - 2)], MPFR_RNDN);
        mpfr_set(sub_im[n - 2], A_im[(n - 1) * n + (n - 2)], MPFR_RNDN);
    }

    mpfr_clear(sigma); mpfr_clear(xk1_re); mpfr_clear(xk1_im);
    mpfr_clear(xk1_abs); mpfr_clear(norm_x);
    mpfr_clear(alpha_re); mpfr_clear(alpha_im);
    mpfr_clear(unorm2); mpfr_clear(inv_unorm); mpfr_clear(alpha_v);
    mpfr_clear(s_re); mpfr_clear(s_im);
    mpfr_clear(prod); mpfr_clear(t_re); mpfr_clear(t_im);
}

/* Phase correction: rotate D = diag(d_k) onto the tridiagonal so the
 * subdiagonal becomes real-positive.  Updates Q by post-multiplying. */
void direct_phase_correct_tridiag_M(mpfr_t* sub_re, mpfr_t* sub_im,
                                             size_t n, mpfr_prec_t bits,
                                             mpfr_t* Q_re, mpfr_t* Q_im,
                                             bool want_Q) {
    mpfr_t d_re, d_im, sr, si, mag, phase_re, phase_im;
    mpfr_t new_d_re, new_d_im, qr, qi, prod;
    mpfr_init2(d_re, bits); mpfr_init2(d_im, bits);
    mpfr_init2(sr,   bits); mpfr_init2(si,   bits);
    mpfr_init2(mag,  bits);
    mpfr_init2(phase_re, bits); mpfr_init2(phase_im, bits);
    mpfr_init2(new_d_re, bits); mpfr_init2(new_d_im, bits);
    mpfr_init2(qr, bits); mpfr_init2(qi, bits); mpfr_init2(prod, bits);

    mpfr_set_ui(d_re, 1, MPFR_RNDN);
    mpfr_set_zero(d_im, 1);

    for (size_t k = 0; k + 1 < n; k++) {
        mpfr_set(sr, sub_re[k], MPFR_RNDN);
        mpfr_set(si, sub_im[k], MPFR_RNDN);
        mpfr_hypot(mag, sr, si, MPFR_RNDN);
        if (mpfr_zero_p(mag)) {
            mpfr_set_ui(phase_re, 1, MPFR_RNDN);
            mpfr_set_zero(phase_im, 1);
        } else {
            mpfr_div(phase_re, sr, mag, MPFR_RNDN);
            mpfr_div(phase_im, si, mag, MPFR_RNDN);
        }
        /* d_{k+1} = d_k * phase. */
        mpfr_mul(new_d_re, d_re, phase_re, MPFR_RNDN);
        mpfr_mul(prod, d_im, phase_im, MPFR_RNDN);
        mpfr_sub(new_d_re, new_d_re, prod, MPFR_RNDN);
        mpfr_mul(new_d_im, d_re, phase_im, MPFR_RNDN);
        mpfr_mul(prod, d_im, phase_re, MPFR_RNDN);
        mpfr_add(new_d_im, new_d_im, prod, MPFR_RNDN);
        mpfr_set(d_re, new_d_re, MPFR_RNDN);
        mpfr_set(d_im, new_d_im, MPFR_RNDN);
        mpfr_set(sub_re[k], mag, MPFR_RNDN);
        mpfr_set_zero(sub_im[k], 1);
        if (want_Q) {
            for (size_t i = 0; i < n; i++) {
                mpfr_set(qr, Q_re[i * n + (k + 1)], MPFR_RNDN);
                mpfr_set(qi, Q_im[i * n + (k + 1)], MPFR_RNDN);
                mpfr_mul(prod, qr, d_re, MPFR_RNDN);
                mpfr_mul(new_d_re, qi, d_im, MPFR_RNDN);
                mpfr_sub(Q_re[i * n + (k + 1)], prod, new_d_re, MPFR_RNDN);
                mpfr_mul(prod, qr, d_im, MPFR_RNDN);
                mpfr_mul(new_d_re, qi, d_re, MPFR_RNDN);
                mpfr_add(Q_im[i * n + (k + 1)], prod, new_d_re, MPFR_RNDN);
            }
        }
    }
    mpfr_clear(d_re); mpfr_clear(d_im);
    mpfr_clear(sr);   mpfr_clear(si);
    mpfr_clear(mag);  mpfr_clear(phase_re); mpfr_clear(phase_im);
    mpfr_clear(new_d_re); mpfr_clear(new_d_im);
    mpfr_clear(qr);   mpfr_clear(qi);  mpfr_clear(prod);
}

/* V = Q Z, where Q is complex n*n and Z is real n*n. */
void direct_compose_complex_Q_real_Z_M(const mpfr_t* Q_re,
                                                const mpfr_t* Q_im,
                                                const mpfr_t* Z, size_t n,
                                                mpfr_prec_t bits,
                                                mpfr_t* V_re, mpfr_t* V_im) {
    mpfr_t s_re, s_im, prod;
    mpfr_init2(s_re, bits); mpfr_init2(s_im, bits); mpfr_init2(prod, bits);
    for (size_t i = 0; i < n; i++) {
        for (size_t k = 0; k < n; k++) {
            mpfr_set_zero(s_re, 1); mpfr_set_zero(s_im, 1);
            for (size_t j = 0; j < n; j++) {
                mpfr_mul(prod, Q_re[i * n + j], Z[j * n + k], MPFR_RNDN);
                mpfr_add(s_re, s_re, prod, MPFR_RNDN);
                mpfr_mul(prod, Q_im[i * n + j], Z[j * n + k], MPFR_RNDN);
                mpfr_add(s_im, s_im, prod, MPFR_RNDN);
            }
            mpfr_set(V_re[i * n + k], s_re, MPFR_RNDN);
            mpfr_set(V_im[i * n + k], s_im, MPFR_RNDN);
        }
    }
    mpfr_clear(s_re); mpfr_clear(s_im); mpfr_clear(prod);
}

/* Build a List[List[...]] of Hermitian eigenvectors at MPFR precision. */
Expr* direct_build_complex_hermitian_eigvec_list_M(const mpfr_t* V_re,
                                                            const mpfr_t* V_im,
                                                            size_t n,
                                                            const size_t* perm) {
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * n);
    for (size_t r = 0; r < n; r++) {
        size_t col = perm[r];
        Expr** comps = (Expr**)malloc(sizeof(Expr*) * n);
        for (size_t i = 0; i < n; i++) {
            if (mpfr_zero_p(V_im[i * n + col])) {
                comps[i] = expr_new_mpfr_copy(V_re[i * n + col]);
            } else {
                Expr** args = (Expr**)malloc(sizeof(Expr*) * 2);
                args[0] = expr_new_mpfr_copy(V_re[i * n + col]);
                args[1] = expr_new_mpfr_copy(V_im[i * n + col]);
                comps[i] = expr_new_function(expr_new_symbol("Complex"),
                                              args, 2);
                free(args);
            }
        }
        rows[r] = expr_new_function(expr_new_symbol("List"), comps, n);
        free(comps);
    }
    Expr* out = expr_new_function(expr_new_symbol("List"), rows, n);
    free(rows);
    return out;
}

/* Top-level Direct kernel for complex Hermitian MPFR input. */
static Expr* direct_complex_hermitian_mpfr(const MatM* A, MateigenWant want,
                                             Expr* k_spec) {
    size_t n = A->n;
    mpfr_prec_t bits = A->bits;
    if (n == 0) return NULL;

    bool want_Q = (want & MATEIGEN_WANT_VECTORS) != 0;

    /* Working copy of A. */
    mpfr_t* W_re = mpfr_array_alloc(n * n, bits);
    mpfr_t* W_im = mpfr_array_alloc(n * n, bits);
    for (size_t i = 0; i < n * n; i++) {
        mpfr_set(W_re[i], A->re[i], MPFR_RNDN);
        mpfr_set(W_im[i], A->im[i], MPFR_RNDN);
    }

    mpfr_t* diag   = mpfr_array_alloc(n, bits);
    mpfr_t* sub_re = mpfr_array_alloc(n, bits);
    mpfr_t* sub_im = mpfr_array_alloc(n, bits);
    mpfr_t* u_re   = mpfr_array_alloc(n, bits);
    mpfr_t* u_im   = mpfr_array_alloc(n, bits);
    mpfr_t* v_re   = mpfr_array_alloc(n, bits);
    mpfr_t* v_im   = mpfr_array_alloc(n, bits);
    mpfr_t* q_re   = mpfr_array_alloc(n, bits);
    mpfr_t* q_im   = mpfr_array_alloc(n, bits);
    mpfr_t* Q_re   = want_Q ? mpfr_array_alloc(n * n, bits) : NULL;
    mpfr_t* Q_im   = want_Q ? mpfr_array_alloc(n * n, bits) : NULL;

    if (n == 1) {
        mpfr_set(diag[0], W_re[0], MPFR_RNDN);
        if (want_Q) {
            mpfr_set_ui(Q_re[0], 1, MPFR_RNDN);
            mpfr_set_zero(Q_im[0], 1);
        }
    } else {
        direct_tridiag_complex_hermitian_M(W_re, W_im, n, bits,
                                            diag, sub_re, sub_im,
                                            Q_re, Q_im, want_Q,
                                            u_re, u_im, v_re, v_im, q_re, q_im);
        direct_phase_correct_tridiag_M(sub_re, sub_im, n, bits,
                                          Q_re, Q_im, want_Q);
    }

    mpfr_array_free(W_re, n * n);
    mpfr_array_free(W_im, n * n);
    mpfr_array_free(u_re, n); mpfr_array_free(u_im, n);
    mpfr_array_free(v_re, n); mpfr_array_free(v_im, n);
    mpfr_array_free(q_re, n); mpfr_array_free(q_im, n);

    /* Real symmetric tridiagonal QR (sub is now real-positive). */
    mpfr_t* Z   = want_Q ? mpfr_array_alloc(n * n, bits) : NULL;
    if (want_Q) {
        for (size_t i = 0; i < n; i++)
            for (size_t j = 0; j < n; j++)
                mpfr_set_si(Z[i * n + j], (i == j) ? 1 : 0, MPFR_RNDN);
    }
    mpfr_t* tmp = mpfr_array_alloc(12, bits);
    int qr_status = (n >= 2)
        ? direct_symtridiag_qr_M(diag, sub_re, n, bits, Z, want_Q, tmp)
        : 0;
    mpfr_array_free(tmp, 12);
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

/* ===================================================================
 * 2d-D: Complex general (non-Hermitian) MPFR Direct kernel
 *
 * Reuses the 2d-B real-general MPFR pipeline via the same 2n x 2n
 * real block embedding M = [[R, -S], [S, R]] used by the machine
 * direct_complex_general_machine.  Grouped Gram-Schmidt at MPFR
 * precision extracts the n distinct A-eigenvectors from the doubled
 * M-eigenvectors.  See the comment above direct_complex_general_machine
 * for the algorithmic motivation.
 * =================================================================== */

static Expr* direct_complex_general_mpfr(const MatM* A, MateigenWant want,
                                           Expr* k_spec) {
    size_t n = A->n;
    mpfr_prec_t bits = A->bits;
    if (n == 0) return NULL;
    size_t N = 2 * n;

    bool want_Q = (want & MATEIGEN_WANT_VECTORS) != 0;

    /* Build the real 2n x 2n block matrix. */
    mpfr_t* H = mpfr_array_alloc(N * N, bits);
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            mpfr_set(H[i * N + j],                 A->re[i * n + j], MPFR_RNDN);
            mpfr_neg(H[i * N + (j + n)],           A->im[i * n + j], MPFR_RNDN);
            mpfr_set(H[(i + n) * N + j],           A->im[i * n + j], MPFR_RNDN);
            mpfr_set(H[(i + n) * N + (j + n)],     A->re[i * n + j], MPFR_RNDN);
        }
    }

    /* Always accumulate Q -- values-only path still needs eigenvectors
     * to disambiguate the doubled spectrum. */
    mpfr_t* Q = mpfr_array_alloc(N * N, bits);
    for (size_t i = 0; i < N; i++)
        for (size_t j = 0; j < N; j++)
            mpfr_set_si(Q[i * N + j], (i == j) ? 1 : 0, MPFR_RNDN);

    mpfr_t* u   = mpfr_array_alloc(N, bits);
    mpfr_t* tmp = mpfr_array_alloc(14, bits);
    if (N >= 3) direct_hessenberg_real_M(H, N, bits, u, Q, tmp);

    mpfr_t* M_eval_re = mpfr_array_alloc(N, bits);
    mpfr_t* M_eval_im = mpfr_array_alloc(N, bits);
    for (size_t i = 0; i < N; i++) {
        mpfr_set_zero(M_eval_re[i], 1);
        mpfr_set_zero(M_eval_im[i], 1);
    }
    int qr_status = direct_qr_real_general_M(H, N, bits, M_eval_re, M_eval_im,
                                               Q, tmp);
    if (qr_status != 0) {
        mpfr_array_free(H, N * N);
        mpfr_array_free(Q, N * N);
        mpfr_array_free(u, N);
        mpfr_array_free(tmp, 14);
        mpfr_array_free(M_eval_re, N);
        mpfr_array_free(M_eval_im, N);
        return NULL;
    }

    /* Eigenvectors of M (in the original basis). */
    mpfr_t* M_evec_re = mpfr_array_alloc(N * N, bits);
    mpfr_t* M_evec_im = mpfr_array_alloc(N * N, bits);
    size_t* identity_perm = (size_t*)malloc(sizeof(size_t) * N);
    for (size_t i = 0; i < N; i++) identity_perm[i] = i;
    schur_compute_eigvecs_M(H, Q, N, bits, M_eval_re, M_eval_im,
                              identity_perm, M_evec_re, M_evec_im);
    free(identity_perm);
    mpfr_array_free(H, N * N);
    mpfr_array_free(Q, N * N);
    mpfr_array_free(u, N);
    mpfr_array_free(tmp, 14);

    /* Tolerances scaled to ~|spec|. */
    mpfr_t spec_norm, group_tol, extract_threshold, mag, dr, di;
    mpfr_init2(spec_norm, bits); mpfr_set_zero(spec_norm, 1);
    mpfr_init2(group_tol, bits);
    mpfr_init2(extract_threshold, bits);
    mpfr_init2(mag, bits);
    mpfr_init2(dr, bits); mpfr_init2(di, bits);
    for (size_t i = 0; i < N; i++) {
        mpfr_hypot(mag, M_eval_re[i], M_eval_im[i], MPFR_RNDN);
        if (mpfr_cmp(mag, spec_norm) > 0) mpfr_set(spec_norm, mag, MPFR_RNDN);
    }
    if (mpfr_zero_p(spec_norm)) mpfr_set_ui(spec_norm, 1, MPFR_RNDN);
    /* group_tol = 2^{-bits/2 + 6} * spec_norm * N -- generous so that
     * any pair (lambda, conj lambda) groups together; tighter would be
     * brittle under non-symmetric QR roundoff. */
    mpfr_set_ui(group_tol, 1, MPFR_RNDN);
    mpfr_div_2si(group_tol, group_tol, (long)bits / 2 - 6, MPFR_RNDN);
    mpfr_mul(group_tol, group_tol, spec_norm, MPFR_RNDN);
    mpfr_mul_ui(group_tol, group_tol, (unsigned long)N, MPFR_RNDN);
    /* extract_threshold = sqrt(n) * 2^{-bits/2 + 6} -- midway between
     * the +J subspace norm ~ sqrt(2) and the -J subspace norm ~ 2^-bits. */
    mpfr_set_ui(extract_threshold, 1, MPFR_RNDN);
    mpfr_div_2si(extract_threshold, extract_threshold, (long)bits / 2 - 6,
                   MPFR_RNDN);
    mpfr_sqrt_ui(mag, (unsigned long)n, MPFR_RNDN);
    mpfr_mul(extract_threshold, extract_threshold, mag, MPFR_RNDN);
    mpfr_t threshold_sq;
    mpfr_init2(threshold_sq, bits);
    mpfr_mul(threshold_sq, extract_threshold, extract_threshold, MPFR_RNDN);

    int* used = (int*)calloc(N, sizeof(int));
    mpfr_t* A_eval_re = mpfr_array_alloc(n, bits);
    mpfr_t* A_eval_im = mpfr_array_alloc(n, bits);
    mpfr_t* A_evec_re = mpfr_array_alloc(n * n, bits);
    mpfr_t* A_evec_im = mpfr_array_alloc(n * n, bits);
    for (size_t i = 0; i < n; i++) {
        mpfr_set_zero(A_eval_re[i], 1); mpfr_set_zero(A_eval_im[i], 1);
    }
    for (size_t i = 0; i < n * n; i++) {
        mpfr_set_zero(A_evec_re[i], 1); mpfr_set_zero(A_evec_im[i], 1);
    }
    mpfr_t* cand_re = mpfr_array_alloc(n, bits);
    mpfr_t* cand_im = mpfr_array_alloc(n, bits);
    mpfr_t pr, pi, vr, vi, pvr, pvi, norm2, inv, prod;
    mpfr_init2(pr, bits); mpfr_init2(pi, bits);
    mpfr_init2(vr, bits); mpfr_init2(vi, bits);
    mpfr_init2(pvr, bits); mpfr_init2(pvi, bits);
    mpfr_init2(norm2, bits); mpfr_init2(inv, bits);
    mpfr_init2(prod, bits);

    size_t out = 0;
    for (size_t i = 0; i < N && out < n; i++) {
        if (used[i]) continue;
        used[i] = 1;
        size_t group_start = out;
        for (size_t j = i; j < N && out < n; j++) {
            if (j != i) {
                if (used[j]) continue;
                mpfr_sub(dr, M_eval_re[j], M_eval_re[i], MPFR_RNDN);
                mpfr_sub(di, M_eval_im[j], M_eval_im[i], MPFR_RNDN);
                mpfr_hypot(mag, dr, di, MPFR_RNDN);
                if (mpfr_cmp(mag, group_tol) > 0) continue;
                used[j] = 1;
            }
            /* Candidate x = (a - d) + i (b + c). */
            for (size_t l = 0; l < n; l++) {
                mpfr_sub(cand_re[l], M_evec_re[j * N + l],
                                       M_evec_im[j * N + (l + n)], MPFR_RNDN);
                mpfr_add(cand_im[l], M_evec_re[j * N + (l + n)],
                                       M_evec_im[j * N + l],       MPFR_RNDN);
            }
            /* Complex Gram-Schmidt (twice). */
            for (int pass = 0; pass < 2; pass++) {
                for (size_t f = group_start; f < out; f++) {
                    mpfr_set_zero(pr, 1); mpfr_set_zero(pi, 1);
                    for (size_t l = 0; l < n; l++) {
                        mpfr_set(vr, A_evec_re[f * n + l], MPFR_RNDN);
                        mpfr_set(vi, A_evec_im[f * n + l], MPFR_RNDN);
                        /* conj(V_f) . cand = (vr - i vi)(cand_re + i cand_im) */
                        mpfr_mul(prod, vr, cand_re[l], MPFR_RNDN);
                        mpfr_add(pr, pr, prod, MPFR_RNDN);
                        mpfr_mul(prod, vi, cand_im[l], MPFR_RNDN);
                        mpfr_add(pr, pr, prod, MPFR_RNDN);
                        mpfr_mul(prod, vr, cand_im[l], MPFR_RNDN);
                        mpfr_add(pi, pi, prod, MPFR_RNDN);
                        mpfr_mul(prod, vi, cand_re[l], MPFR_RNDN);
                        mpfr_sub(pi, pi, prod, MPFR_RNDN);
                    }
                    for (size_t l = 0; l < n; l++) {
                        mpfr_set(vr, A_evec_re[f * n + l], MPFR_RNDN);
                        mpfr_set(vi, A_evec_im[f * n + l], MPFR_RNDN);
                        /* (pr + i pi)(vr + i vi) */
                        mpfr_mul(pvr, pr, vr, MPFR_RNDN);
                        mpfr_mul(prod, pi, vi, MPFR_RNDN);
                        mpfr_sub(pvr, pvr, prod, MPFR_RNDN);
                        mpfr_mul(pvi, pr, vi, MPFR_RNDN);
                        mpfr_mul(prod, pi, vr, MPFR_RNDN);
                        mpfr_add(pvi, pvi, prod, MPFR_RNDN);
                        mpfr_sub(cand_re[l], cand_re[l], pvr, MPFR_RNDN);
                        mpfr_sub(cand_im[l], cand_im[l], pvi, MPFR_RNDN);
                    }
                }
            }
            mpfr_set_zero(norm2, 1);
            for (size_t l = 0; l < n; l++) {
                mpfr_mul(prod, cand_re[l], cand_re[l], MPFR_RNDN);
                mpfr_add(norm2, norm2, prod, MPFR_RNDN);
                mpfr_mul(prod, cand_im[l], cand_im[l], MPFR_RNDN);
                mpfr_add(norm2, norm2, prod, MPFR_RNDN);
            }
            if (mpfr_cmp(norm2, threshold_sq) < 0) continue;
            mpfr_sqrt(inv, norm2, MPFR_RNDN);
            mpfr_ui_div(inv, 1, inv, MPFR_RNDN);
            for (size_t l = 0; l < n; l++) {
                mpfr_mul(A_evec_re[out * n + l], cand_re[l], inv, MPFR_RNDN);
                mpfr_mul(A_evec_im[out * n + l], cand_im[l], inv, MPFR_RNDN);
            }
            mpfr_set(A_eval_re[out], M_eval_re[i], MPFR_RNDN);
            mpfr_set(A_eval_im[out], M_eval_im[i], MPFR_RNDN);
            out++;
        }
    }

    free(used);
    mpfr_array_free(M_eval_re, N); mpfr_array_free(M_eval_im, N);
    mpfr_array_free(M_evec_re, N * N); mpfr_array_free(M_evec_im, N * N);
    mpfr_array_free(cand_re, n); mpfr_array_free(cand_im, n);
    mpfr_clear(spec_norm); mpfr_clear(group_tol);
    mpfr_clear(extract_threshold); mpfr_clear(threshold_sq);
    mpfr_clear(mag); mpfr_clear(dr); mpfr_clear(di);
    mpfr_clear(pr); mpfr_clear(pi);
    mpfr_clear(vr); mpfr_clear(vi);
    mpfr_clear(pvr); mpfr_clear(pvi);
    mpfr_clear(norm2); mpfr_clear(inv); mpfr_clear(prod);

    if (out != n) {
        mpfr_array_free(A_eval_re, n); mpfr_array_free(A_eval_im, n);
        mpfr_array_free(A_evec_re, n * n); mpfr_array_free(A_evec_im, n * n);
        return NULL;
    }

    size_t* perm = (size_t*)malloc(sizeof(size_t) * n);
    direct_sort_perm_desc_abs_complex_M(A_eval_re, A_eval_im, n, perm);

    Expr* result;
    if (want_Q) {
        mpfr_t* V_re = mpfr_array_alloc(n * n, bits);
        mpfr_t* V_im = mpfr_array_alloc(n * n, bits);
        for (size_t k = 0; k < n; k++) {
            size_t src = perm[k];
            for (size_t l = 0; l < n; l++) {
                mpfr_set(V_re[k * n + l], A_evec_re[src * n + l], MPFR_RNDN);
                mpfr_set(V_im[k * n + l], A_evec_im[src * n + l], MPFR_RNDN);
            }
        }
        result = direct_build_complex_eigenvector_list_M(V_re, V_im, n);
        mpfr_array_free(V_re, n * n);
        mpfr_array_free(V_im, n * n);
    } else {
        result = direct_build_complex_eigenvalue_list_M(A_eval_re, A_eval_im,
                                                          n, perm);
    }

    mpfr_array_free(A_eval_re, n); mpfr_array_free(A_eval_im, n);
    mpfr_array_free(A_evec_re, n * n); mpfr_array_free(A_evec_im, n * n);
    free(perm);

    return direct_apply_k_spec_list(result, k_spec);
}

/* Dispatcher: route an MPFR-precision matrix through the appropriate
 * Direct kernel.  All four shapes (real symmetric, real general,
 * complex Hermitian, complex general) now have MPFR implementations. */
static Expr* direct_dispatch_mpfr(Expr* m, Expr* a, int64_t n,
                                    mpfr_prec_t bits,
                                    MateigenWant want, Expr* k_spec) {
    if (a != NULL) return NULL;
    if (n <= 0)    return NULL;

    MatM A;
    if (!matM_load(m, (size_t)n, bits, &A)) return NULL;

    Expr* out = NULL;
    /* Tolerance for Hermitian / symmetry detection: n * 2^{-bits+4} * ||A||_inf. */
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

    if (A.is_complex) {
        if (matM_is_hermitian(&A, tol))
            out = direct_complex_hermitian_mpfr(&A, want, k_spec);
        else
            out = direct_complex_general_mpfr(&A, want, k_spec);
    } else {
        if (matM_is_real_symmetric(&A, tol))
            out = direct_real_sym_mpfr(&A, want, k_spec);
        else
            out = direct_real_general_mpfr(&A, want, k_spec);
    }

    mpfr_clear(norm); mpfr_clear(tol); mpfr_clear(factor);
    matM_free(&A);
    return out;
}
#endif /* USE_MPFR */


/* Dispatcher entry point: route a numeric matrix through the
 * appropriate "Direct" kernel.  Returns NULL when the matrix shape
 * isn't yet supported by a numerical kernel so the caller can fall
 * back to the symbolic path.  This NULL return is also used for
 * Eigenvalues / Eigenvectors combined with a generalised pencil
 * ({m, a}) -- generalised numeric eigenvalues are not part of the
 * current numerical scope.
 *
 * Implemented kernels:
 *   - Real symmetric (machine precision):           values + vectors.
 *   - Real non-symmetric (machine precision):       values + vectors.
 *   - Complex Hermitian (machine precision):        values + vectors.
 *   - Complex non-Hermitian (machine precision):    values + vectors.
 *   - Real symmetric MPFR (step 2d-A):              values + vectors.
 *   - Real non-symmetric MPFR (step 2d-B):          values + vectors.
 *   - Complex Hermitian MPFR (step 2d-C):           values + vectors.
 *   - Complex non-Hermitian MPFR (step 2d-D):       values + vectors.
 *
 * Generalised pencils ({m, a}) still flow through the symbolic path.
 */
static Expr* direct_dispatch_machine(Expr* m, Expr* a, int64_t n,
                                       MateigenWant want, Expr* k_spec) {
    if (a != NULL) return NULL;          /* generalised: symbolic only */
    if (n <= 0)    return NULL;

    MatD A;
    if (!matD_load(m, (size_t)n, &A)) return NULL;

    Expr* out = NULL;
    if (A.is_complex) {
        double norm = matD_norm_inf_complex(&A);
        double herm_tol = 1e-12 * (norm == 0.0 ? 1.0 : norm) * (double)A.n;
        if (matD_is_hermitian(&A, herm_tol)) {
            out = direct_complex_hermitian_machine(&A, want, k_spec);
        } else {
            out = direct_complex_general_machine(&A, want, k_spec);
        }
    } else {
        double norm = matD_norm_inf_real(A.re, A.n);
        double sym_tol = 1e-12 * (norm == 0.0 ? 1.0 : norm) * (double)A.n;
        if (matD_is_real_symmetric(&A, sym_tol)) {
            out = direct_real_sym_machine(&A, want, k_spec);
        } else {
            out = direct_real_general_machine(&A, want, k_spec);
        }
    }
    matD_free(&A);
    return out;
}

/* Top-level "Direct" dispatcher.  Picks the MPFR kernel when any
 * input leaf carries MPFR precision (per common_scan_inexact); falls
 * back to the machine-precision kernel otherwise, or when the MPFR
 * kernel for the matrix shape isn't yet wired. */
Expr* direct_dispatch(Expr* m, Expr* a, int64_t n,
                               MateigenWant want, Expr* k_spec) {
    CommonInexactInfo info = common_scan_inexact(m);
    if (a) {
        CommonInexactInfo ia = common_scan_inexact(a);
        if (ia.has_inexact && (!info.has_inexact || ia.min_bits < info.min_bits))
            info = ia;
    }
    if (info.has_inexact && info.min_bits > 53) {
        Expr* out = direct_dispatch_mpfr(m, a, n, (mpfr_prec_t)info.min_bits,
                                          want, k_spec);
        if (out) return out;
        /* MPFR kernel not yet wired for this matrix shape -- fall
         * through to the machine kernel.  The machine kernel will
         * coerce the MPFR cells to doubles via eigen_leaf_to_double,
         * which is the closest behaviour-preserving fallback. */
    }
    return direct_dispatch_machine(m, a, n, want, k_spec);
}

#ifdef USE_MPFR
/* Public wrapper exposed in eigen.h.  Allocates workspace, runs
 * Hessenberg + Francis QR on a copy of the input, returns eigenvalues
 * as parallel (re, im) MPFR arrays.  Used by root_numeric.c as its
 * companion-matrix all-roots backend. */
int eigen_all_eigenvalues_real_mpfr(mpfr_t* A, size_t n, mpfr_prec_t bits,
                                     mpfr_t* eval_re, mpfr_t* eval_im) {
    if (n == 0) return 0;
    if (n == 1) {
        mpfr_set(eval_re[0], A[0], MPFR_RNDN);
        mpfr_set_zero(eval_im[0], 1);
        return 0;
    }

    mpfr_t* Q   = mpfr_array_alloc(n * n, bits);
    mpfr_t* u   = mpfr_array_alloc(n, bits);
    mpfr_t* tmp = mpfr_array_alloc(14, bits);

    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            mpfr_set_si(Q[i * n + j], (i == j) ? 1 : 0, MPFR_RNDN);

    if (n >= 3) direct_hessenberg_real_M(A, n, bits, u, Q, tmp);

    for (size_t i = 0; i < n; i++) {
        mpfr_set_zero(eval_re[i], 1);
        mpfr_set_zero(eval_im[i], 1);
    }
    int status = direct_qr_real_general_M(A, n, bits,
                                          eval_re, eval_im, Q, tmp);

    mpfr_array_free(Q,   n * n);
    mpfr_array_free(u,   n);
    mpfr_array_free(tmp, 14);
    return status;
}
#endif
