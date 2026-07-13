/* src/linalg/qrdecomp_machine.c
 *
 * Phase-3 machine-precision QR kernel.
 *
 * This file is the LAPACK fast path of QRDecomposition.  It is invoked
 * by qr_dispatch when common_scan_inexact reports an inexact input at
 * <= 53 bits of precision (i.e. the matrix is dominated by IEEE doubles).
 * On any failure -- USE_LAPACK off, a leaf the loader can't reduce to a
 * double, a non-zero LAPACK `info` -- it returns NULL and qr_dispatch
 * falls back to the symbolic Modified-Gram-Schmidt path.
 *
 * The kernel is deliberately a single file with no shared helpers other
 * than the lapack.c wrappers: this is small enough to keep all the
 * row/column-major bookkeeping in one place, and there is no other
 * builtin yet that wants the same load/store machinery.  When Inverse /
 * LinearSolve / Det / Eigenvalues acquire their own LAPACK fast paths
 * we can lift the load/store helpers into a shared `linalg/mach.c`.
 *
 * High-level flow:
 *
 *   1. mach_load_matrix     -- walk the n x p list-of-list expression
 *                              and copy every leaf into a column-major
 *                              double buffer A_cm of length n*p (real)
 *                              or 2*n*p (complex; interleaved re,im).
 *                              Returns false if any leaf isn't a
 *                              recognisable numeric value.
 *
 *   2. LAPACK QR            -- dgeqrf / dgeqp3 (real) or zgeqrf /
 *                              zgeqp3 (complex).  Overwrites A_cm
 *                              with the Householder reflectors in the
 *                              lower triangle and R in the upper.
 *
 *   3. Extract R            -- copy the (rank x p) row-major upper-
 *                              trapezoidal block out of A_cm before
 *                              dorgqr / zungqr overwrites it.  Rank is
 *                              determined from the R diagonal using the
 *                              standard LAPACK rank-revealing cutoff.
 *
 *   4. dorgqr / zungqr      -- form Q (n x k, k = min(n, p)) in place
 *                              over A_cm.
 *
 *   5. Build outputs        -- wrap Q (transposed to rank x n, with
 *                              optional Hermitian conjugation for the
 *                              complex case) and R into Mathilda Lists,
 *                              build the permutation matrix from jpvt
 *                              when pivoting was requested.
 *
 * Memory contract.  Standard builtin contract: this file does NOT call
 * expr_free on the input `m` -- the evaluator owns it (MEMORY.md /
 * SPEC.md §4.1).  Every malloc that is not handed to the caller is
 * matched by a free along every exit path.
 */

#include "qrdecomp_internal.h"
#include "linalg.h"
#include "lapack.h"
#include "expr.h"
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
 * of qr_machine_dispatch.  Wrapping them removes -Wunused-function noise
 * on builds compiled without LAPACK. */

/* ---------------------------------------------------------------------
 * One-shot warning helper.  Mirrors the matsol_warn_once / chop_warn_once
 * convention used elsewhere in linalg: a 64-bit counter passed by
 * pointer suppresses repeated reports of the same condition across a
 * REPL session. */
static void qr_machine_warn_once(uint64_t* counter, const char* msg)
{
    if (*counter) return;
    *counter = 1;
    fprintf(stderr, "%s", msg);
}

/* ---------------------------------------------------------------------
 * Numeric-leaf -> double conversion.
 *
 * Recognised forms: Integer, BigInt (via mpz_get_d), Real, MPFR, the
 * exact 2-arg Rational[p, q] head, and the 2-arg Complex[re, im] head.
 * Anything else triggers a fall-back to the symbolic kernel.
 *
 * `*out_im` is written even for purely-real inputs (set to 0.0); the
 * complex-detection scan reads it to decide whether to allocate an
 * imaginary column.
 *
 * Returns false on any unrecognised leaf, so the caller knows to bail
 * out and route through the symbolic dispatcher (which understands
 * symbolic content). */
static bool mach_leaf_to_double(Expr* e, double* out_re, double* out_im)
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
                    if (mach_leaf_to_double(e->data.function.args[0], &p, &dummy)
                        && mach_leaf_to_double(e->data.function.args[1], &q, &dummy)
                        && q != 0.0) {
                        *out_re = p / q;
                        return true;
                    }
                    return false;
                }
                if (h == SYM_Complex && e->data.function.arg_count == 2) {
                    double r, i, dummy;
                    if (mach_leaf_to_double(e->data.function.args[0], &r, &dummy)
                        && mach_leaf_to_double(e->data.function.args[1], &i, &dummy)) {
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

/* True iff `e` numerically represents a value with a non-zero imaginary
 * part.  Used by the first scan over the input matrix to decide whether
 * we need the z* (complex) LAPACK routines. */
static bool mach_leaf_is_complex(Expr* e)
{
    double r, i;
    return mach_leaf_to_double(e, &r, &i) && i != 0.0;
}

/* ---------------------------------------------------------------------
 * Load the n x p input matrix into a freshly-allocated column-major
 * double buffer.  For real inputs A_cm has length n*p; for complex
 * inputs A_cm has length 2*n*p with (re, im) interleaved per element
 * (matching the LAPACK `complex*16` layout).
 *
 * Returns true on success.  On failure (any leaf isn't a numeric form
 * we can convert), returns false with *out_A and *out_is_complex
 * undefined and nothing leaked. */
static bool mach_load_matrix(Expr* m, int n, int p,
                             double** out_A, bool* out_is_complex)
{
    if (m->type != EXPR_FUNCTION) return false;

    /* First pass: detect any complex entries. */
    bool is_complex = false;
    for (int i = 0; i < n && !is_complex; i++) {
        Expr* row = m->data.function.args[i];
        if (row->type != EXPR_FUNCTION) return false;
        for (int j = 0; j < p && !is_complex; j++) {
            if (mach_leaf_is_complex(row->data.function.args[j])) {
                is_complex = true;
            }
        }
    }

    size_t stride = is_complex ? 2 : 1;
    double* A = (double*)malloc(stride * (size_t)n * (size_t)p * sizeof(double));
    if (!A) return false;

    /* Second pass: column-major copy.  LAPACK addresses A_cm[i + j*lda]
     * for the (i, j) entry, lda = n.  For complex, each entry occupies
     * 2 consecutive doubles starting at 2*(i + j*n). */
    for (int i = 0; i < n; i++) {
        Expr* row = m->data.function.args[i];
        for (int j = 0; j < p; j++) {
            double re, im;
            if (!mach_leaf_to_double(row->data.function.args[j], &re, &im)) {
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

/* ---------------------------------------------------------------------
 * Determine the numerical rank of the (m x p) upper-trapezoidal R that
 * lives in the upper triangle of A_cm after dgeqrf / dgeqp3 / their z*
 * variants.  We use the standard LAPACK rank-revealing cutoff:
 *
 *     rank = #{ i : |R[i,i]| > max(m, p) * eps * |R[0,0]| }
 *
 * where the count starts at i = 0 and stops the moment one diagonal
 * falls under the threshold (R's diagonals are non-increasing in
 * magnitude when pivoting was used; for the no-pivot path the cutoff
 * still acts as a sanity guard against silent rank loss).
 *
 * For complex R the magnitude is sqrt(re^2 + im^2). */
static int mach_numerical_rank(const double* A_cm, int n, int p,
                                bool is_complex)
{
    int k = (n < p) ? n : p;
    if (k == 0) return 0;

    size_t stride = is_complex ? 2 : 1;

    double max_diag = 0.0;
    for (int i = 0; i < k; i++) {
        size_t off = stride * ((size_t)i + (size_t)i * (size_t)n);
        double mag = is_complex
            ? hypot(A_cm[off], A_cm[off + 1])
            : fabs(A_cm[off]);
        if (mag > max_diag) max_diag = mag;
    }
    if (max_diag == 0.0) return 0;

    /* LAPACK's documented "loose" rank-revealing tolerance.  Tight
     * enough to drop genuinely zero rows, loose enough not to discard
     * ill-conditioned-but-non-singular factors. */
    double tol = ((double)((n > p ? n : p))) * DBL_EPSILON * max_diag;

    int rank = 0;
    for (int i = 0; i < k; i++) {
        size_t off = stride * ((size_t)i + (size_t)i * (size_t)n);
        double mag = is_complex
            ? hypot(A_cm[off], A_cm[off + 1])
            : fabs(A_cm[off]);
        if (mag > tol) rank++;
        else break;
    }
    return rank;
}

/* ---------------------------------------------------------------------
 * Output construction helpers.  Each one builds the requested Mathilda
 * subtree from a freshly-allocated double buffer; the caller passes a
 * non-NULL imag buffer for complex outputs.  We never embed a
 * Complex[r, 0.0] head -- following the existing matD_load_M /
 * numericalize convention, a Complex with literal-zero imaginary part
 * is collapsed to a bare Real. */
static Expr* mach_make_scalar(double re, double im, bool is_complex)
{
    if (!is_complex || im == 0.0) return expr_new_real(re);
    Expr** args = (Expr**)malloc(sizeof(Expr*) * 2);
    args[0] = expr_new_real(re);
    args[1] = expr_new_real(im);
    Expr* z = expr_new_function(expr_new_symbol(SYM_Complex), args, 2);
    free(args);
    return z;
}

/* Build the `q` output (rank x n) from the LAPACK Q buffer (n x k
 * column-major, k = min(n, p)).  Only the first `rank` columns of Q
 * are used.
 *
 *   - Real input:   q[j, i] = Q[i, j]            (Transpose)
 *   - Complex input: q[j, i] = Conjugate(Q[i, j]) (ConjugateTranspose)
 *
 * The conjugation matches the existing symbolic-kernel convention in
 * build_q_from_Q (qrdecomp.c) -- it preserves the public identity
 * m == ConjugateTranspose[q] . r for both real and complex matrices,
 * and avoids emitting Conjugate[Conjugate[x]] for real inputs (which
 * the printer would render verbatim). */
static Expr* mach_build_q(const double* Q_cm, int n, int rank,
                          bool is_complex)
{
    size_t stride = is_complex ? 2 : 1;

    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)rank);
    for (int j = 0; j < rank; j++) {
        Expr** elems = NULL;
        if (n > 0) elems = (Expr**)malloc(sizeof(Expr*) * (size_t)n);
        for (int i = 0; i < n; i++) {
            size_t off = stride * ((size_t)i + (size_t)j * (size_t)n);
            double re = Q_cm[off];
            double im = is_complex ? Q_cm[off + 1] : 0.0;
            if (is_complex) im = -im;   /* Hermitian conjugate */
            elems[i] = mach_make_scalar(re, im, is_complex);
        }
        rows[j] = expr_new_function(expr_new_symbol(SYM_List),
                                     elems, (size_t)n);
        if (elems) free(elems);
    }
    Expr* q = expr_new_function(expr_new_symbol(SYM_List),
                                  rows, (size_t)rank);
    free(rows);
    return q;
}

/* Build the `r` output (rank x p) from the saved upper triangle of A_cm
 * captured before dorgqr ran.  `R_storage` is row-major, of size
 * rank * p * (is_complex ? 2 : 1).  Entries strictly below the leading
 * echelon must already be exact zero (the caller ensures this). */
static Expr* mach_build_r(const double* R_storage, int rank, int p,
                          bool is_complex)
{
    size_t stride = is_complex ? 2 : 1;

    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)rank);
    for (int j = 0; j < rank; j++) {
        Expr** elems = NULL;
        if (p > 0) elems = (Expr**)malloc(sizeof(Expr*) * (size_t)p);
        for (int k = 0; k < p; k++) {
            size_t off = stride * ((size_t)j * (size_t)p + (size_t)k);
            double re = R_storage[off];
            double im = is_complex ? R_storage[off + 1] : 0.0;
            elems[k] = mach_make_scalar(re, im, is_complex);
        }
        rows[j] = expr_new_function(expr_new_symbol(SYM_List),
                                     elems, (size_t)p);
        if (elems) free(elems);
    }
    Expr* r = expr_new_function(expr_new_symbol(SYM_List),
                                  rows, (size_t)rank);
    free(rows);
    return r;
}

/* Build a p x p Integer permutation matrix from LAPACK's 1-indexed
 * jpvt: the j-th column of m.P is column (jpvt[j] - 1) of m, so
 * P[jpvt[j] - 1, j] = 1.  Exactly matches the symbolic kernel's
 * build_perm_matrix output so the unit tests can compare against the
 * same Integer 0 / 1 entries. */
static Expr* mach_build_perm(const int* jpvt, int p)
{
    Expr** rows = (Expr**)malloc(sizeof(Expr*) * (size_t)p);
    for (int i = 0; i < p; i++) {
        Expr** elems = NULL;
        if (p > 0) elems = (Expr**)malloc(sizeof(Expr*) * (size_t)p);
        for (int j = 0; j < p; j++) {
            int src_col = jpvt[j] - 1;
            elems[j] = expr_new_integer(src_col == i ? 1 : 0);
        }
        rows[i] = expr_new_function(expr_new_symbol(SYM_List),
                                     elems, (size_t)p);
        if (elems) free(elems);
    }
    Expr* P = expr_new_function(expr_new_symbol(SYM_List), rows, (size_t)p);
    free(rows);
    return P;
}
#endif /* USE_LAPACK */

/* ---------------------------------------------------------------------
 * The kernel.
 *
 * Lives behind the qr_machine_dispatch entry point declared in
 * qrdecomp_internal.h.  Returns NULL on any failure path; the caller
 * (qr_dispatch) treats NULL as "fall back to symbolic" and never
 * frees the input.
 * ------------------------------------------------------------------ */
Expr* qr_machine_dispatch(Expr* m, int n, int p, const QrOpts* opts)
{
#ifndef USE_LAPACK
    (void)m; (void)n; (void)p; (void)opts;
    return NULL;
#else
    static uint64_t lapack_warn_counter = 0;
    /* Load the matrix into a column-major double buffer.  On a non-
     * numeric leaf we return NULL; qr_dispatch then re-runs the input
     * through the symbolic pipeline. */
    double* A_cm = NULL;
    bool is_complex = false;
    if (!mach_load_matrix(m, n, p, &A_cm, &is_complex)) return NULL;

    int k = (n < p) ? n : p;
    size_t stride = is_complex ? 2 : 1;
    int lda = n;

    /* tau holds k Householder reflector scalars (k complex doubles for
     * z* variants -> 2*k doubles). */
    double* tau = (double*)malloc(stride * (size_t)k * sizeof(double));
    if (!tau) { free(A_cm); return NULL; }

    /* jpvt: 1-indexed column-permutation array used when pivoting.  We
     * always allocate so we can read it back uniformly across the two
     * code paths; for the no-pivot case it stays untouched. */
    int* jpvt = NULL;
    if (opts->pivoting) {
        jpvt = (int*)malloc((size_t)p * sizeof(int));
        if (!jpvt) { free(A_cm); free(tau); return NULL; }
    }

    /* Factor via dgeqp3 / dgeqrf (real) or zgeqp3 / zgeqrf (complex). */
    int info = 0;
    if (is_complex) {
        info = opts->pivoting
            ? mat_lapack_zgeqp3(n, p, A_cm, lda, jpvt, tau)
            : mat_lapack_zgeqrf(n, p, A_cm, lda, tau);
    } else {
        info = opts->pivoting
            ? mat_lapack_dgeqp3(n, p, A_cm, lda, jpvt, tau)
            : mat_lapack_dgeqrf(n, p, A_cm, lda, tau);
    }
    if (info != 0) {
        if (info < 0) {
            /* info < 0 means either a stub (USE_LAPACK off, impossible
             * here) or a bad argument -- log once and fall back. */
            qr_machine_warn_once(&lapack_warn_counter,
                "QRDecomposition: LAPACK fast path unavailable; "
                "falling back to symbolic kernel.\n");
        } else {
            qr_machine_warn_once(&lapack_warn_counter,
                "QRDecomposition: LAPACK factorisation reported a "
                "numerical issue; falling back to symbolic kernel.\n");
        }
        free(A_cm); free(tau); if (jpvt) free(jpvt);
        return NULL;
    }

    /* Determine numerical rank from the R diagonal that now lives in
     * the upper triangle of A_cm. */
    int rank = mach_numerical_rank(A_cm, n, p, is_complex);

    /* Snapshot the (rank x p) upper trapezoidal R into row-major
     * storage *before* dorgqr overwrites A_cm with Q.  Entries strictly
     * below the leading echelon (i.e. R[j, k] with k < j) are zeroed --
     * the LAPACK lower-triangle of A actually holds the Householder
     * reflector data, which we must not leak into r. */
    double* R_row = NULL;
    if (rank > 0) {
        R_row = (double*)malloc(stride * (size_t)rank * (size_t)p
                                * sizeof(double));
        if (!R_row) {
            free(A_cm); free(tau); if (jpvt) free(jpvt);
            return NULL;
        }
        for (int j = 0; j < rank; j++) {
            for (int kcol = 0; kcol < p; kcol++) {
                size_t dst = stride * ((size_t)j * (size_t)p + (size_t)kcol);
                if (kcol < j) {
                    /* Strictly below the leading diagonal: literal zero. */
                    R_row[dst] = 0.0;
                    if (is_complex) R_row[dst + 1] = 0.0;
                } else {
                    size_t src = stride * ((size_t)j + (size_t)kcol * (size_t)n);
                    R_row[dst] = A_cm[src];
                    if (is_complex) R_row[dst + 1] = A_cm[src + 1];
                }
            }
        }
    }

    /* Form Q via dorgqr / zungqr.  We materialise the full leading
     * n x k block; rank-deficient columns (beyond `rank`) are produced
     * by LAPACK as extensions of the orthonormal basis and we just
     * ignore them when building `q`.  Pass k as the number of
     * reflectors -- always min(n, p), regardless of the numerical
     * rank. */
    if (k > 0) {
        if (is_complex) {
            info = mat_lapack_zungqr(n, k, k, A_cm, lda, tau);
        } else {
            info = mat_lapack_dorgqr(n, k, k, A_cm, lda, tau);
        }
        if (info != 0) {
            qr_machine_warn_once(&lapack_warn_counter,
                "QRDecomposition: LAPACK orthogonal-basis formation "
                "failed; falling back to symbolic kernel.\n");
            free(A_cm); free(tau); if (jpvt) free(jpvt);
            if (R_row) free(R_row);
            return NULL;
        }
    }

    free(tau);

    /* Build q and r.  Empty-rank case -> empty `{}` matrices, matching
     * the symbolic kernel's all-zero-input behaviour. */
    Expr* q;
    Expr* r;
    if (rank == 0) {
        q = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
        r = expr_new_function(expr_new_symbol(SYM_List), NULL, 0);
    } else {
        q = mach_build_q(A_cm, n, rank, is_complex);
        r = mach_build_r(R_row, rank, p, is_complex);
    }
    free(A_cm);
    if (R_row) free(R_row);

    /* Assemble the result list. */
    Expr* result;
    if (opts->pivoting) {
        Expr* P = mach_build_perm(jpvt, p);
        free(jpvt);
        Expr** items = (Expr**)malloc(sizeof(Expr*) * 3);
        items[0] = q; items[1] = r; items[2] = P;
        result = expr_new_function(expr_new_symbol(SYM_List), items, 3);
        free(items);
    } else {
        Expr** items = (Expr**)malloc(sizeof(Expr*) * 2);
        items[0] = q; items[1] = r;
        result = expr_new_function(expr_new_symbol(SYM_List), items, 2);
        free(items);
    }
    return result;
#endif /* USE_LAPACK */
}
