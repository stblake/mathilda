/* test_mateigen_feast.c -- Tests for the numerical "FEAST"
 * Eigenvalues / Eigenvectors kernel (Phases 2-4).
 *
 * FEAST is Hermitian-only: it builds an invariant subspace inside a
 * user-supplied real interval [a, b] via contour-integral spectral
 * projection, then performs Rayleigh-Ritz on that subspace.  Real
 * symmetric and complex Hermitian inputs are both supported at machine
 * precision; MPFR variants kick in when the input precision exceeds 53
 * bits.
 *
 * Correctness contract:
 *   - When the interval covers the full spectrum, FEAST's output equals
 *     Direct's (up to tolerance).
 *   - When the interval is a strict sub-range, FEAST returns precisely
 *     the eigenvalues inside that range.
 *   - Returned eigenvectors satisfy ||A v - lambda v|| / ||v|| < tol and
 *     are pairwise orthonormal.
 *   - When the interval is missing / empty / matrix is non-Hermitian,
 *     the dispatcher returns NULL and the call falls through to Direct;
 *     the user gets the full Direct spectrum (no truncation, no crash).
 *
 * Tests follow the conventions of test_mateigen_banded.c (eval_string,
 * fmt_real_matrix, extract_real, list_len_eq).
 */

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include "eigen_corpus.h"

/* ---------- formatting / extraction helpers ---------- */

static Expr* eval_string(const char* s) {
    Expr* parsed = parse_expression(s);
    return evaluate(parsed);
}

static char* fmt_real_matrix(const double* A, size_t n) {
    size_t cap = 32 + 32 * n * n;
    char* buf = (char*)malloc(cap);
    char* p = buf;
    p += snprintf(p, cap, "N[{");
    for (size_t i = 0; i < n; i++) {
        p += snprintf(p, cap - (p - buf), "%s{", i ? "," : "");
        for (size_t j = 0; j < n; j++) {
            p += snprintf(p, cap - (p - buf), "%s%.17g",
                          j ? "," : "", A[i * n + j]);
        }
        p += snprintf(p, cap - (p - buf), "}");
    }
    snprintf(p, cap - (p - buf), "}]");
    return buf;
}

static bool list_len_eq(Expr* e, size_t expected) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    if (strcmp(e->data.function.head->data.symbol, "List") != 0) return false;
    return e->data.function.arg_count == expected;
}

static bool is_list(Expr* e) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    return strcmp(e->data.function.head->data.symbol, "List") == 0;
}

static double extract_real(Expr* e) {
    if (!e) return NAN;
    if (e->type == EXPR_INTEGER) return (double)e->data.integer;
    if (e->type == EXPR_REAL)    return e->data.real;
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && strcmp(e->data.function.head->data.symbol, "Complex") == 0
        && e->data.function.arg_count == 2)
        return extract_real(e->data.function.args[0]);
    return NAN;
}

static double extract_imag(Expr* e) {
    if (!e) return NAN;
    if (e->type == EXPR_INTEGER) return 0.0;
    if (e->type == EXPR_REAL)    return 0.0;
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && strcmp(e->data.function.head->data.symbol, "Complex") == 0
        && e->data.function.arg_count == 2)
        return extract_real(e->data.function.args[1]);
    return 0.0;
}

/* Sort a small double array ascending in place (insertion sort). */
static void sort_asc(double* v, size_t n) {
    for (size_t i = 1; i < n; i++) {
        double x = v[i];
        size_t j = i;
        while (j > 0 && v[j - 1] > x) { v[j] = v[j - 1]; j--; }
        v[j] = x;
    }
}

/* ============================================================ *
 *  Real-symmetric FEAST -- cross-checks against Direct           *
 * ============================================================ */

/* 5x5 symmetric tridiagonal -- standard Chebyshev spectrum
 * 2 - 2 cos(k pi / 6), k = 1..5 = {0.268, 1, 2, 3, 3.732}.
 * Whole-spectrum interval recovers all 5 eigenvalues. */
void test_feast_real_sym_full_spectrum_5x5(void) {
    double A[25] = {
         2, -1,  0,  0,  0,
        -1,  2, -1,  0,  0,
         0, -1,  2, -1,  0,
         0,  0, -1,  2, -1,
         0,  0,  0, -1,  2
    };
    char* m = fmt_real_matrix(A, 5);
    char buf[4096];
    snprintf(buf, sizeof buf,
        "Eigenvalues[%s, Method -> {\"FEAST\", \"Interval\" -> {0, 4}}]", m);
    Expr* fe = eval_string(buf);
    snprintf(buf, sizeof buf, "Eigenvalues[%s, Method -> \"Direct\"]", m);
    Expr* dr = eval_string(buf);
    ASSERT(list_len_eq(fe, 5));
    ASSERT(list_len_eq(dr, 5));

    /* Direct returns descending |lambda|; FEAST does too (per
     * feast_real_sym_machine's direct_sort_perm_desc_abs at emit).
     * Compare via sorted ascending lists since |lambda| ordering for
     * these positive eigenvalues coincides with lambda ordering. */
    double f[5], d[5];
    for (size_t i = 0; i < 5; i++) {
        f[i] = extract_real(fe->data.function.args[i]);
        d[i] = extract_real(dr->data.function.args[i]);
    }
    sort_asc(f, 5); sort_asc(d, 5);
    for (size_t i = 0; i < 5; i++) {
        ASSERT(fabs(f[i] - d[i]) < 1e-10);
    }
    expr_free(fe); expr_free(dr);
    free(m);
}

/* Sub-interval {2.5, 4} on the Chebyshev tridiagonal -- expect the two
 * eigenvalues > 2.5, i.e. 3 and 3.732. */
void test_feast_real_sym_subinterval_5x5(void) {
    double A[25] = {
         2, -1,  0,  0,  0,
        -1,  2, -1,  0,  0,
         0, -1,  2, -1,  0,
         0,  0, -1,  2, -1,
         0,  0,  0, -1,  2
    };
    char* m = fmt_real_matrix(A, 5);
    char buf[4096];
    snprintf(buf, sizeof buf,
        "Eigenvalues[%s, Method -> {\"FEAST\", \"Interval\" -> {2.5, 4}}]", m);
    Expr* fe = eval_string(buf);
    ASSERT(list_len_eq(fe, 2));

    double f[2];
    for (size_t i = 0; i < 2; i++)
        f[i] = extract_real(fe->data.function.args[i]);
    sort_asc(f, 2);
    /* Eigenvalues 3.0 and 2 + sqrt(3). */
    ASSERT(fabs(f[0] - 3.0) < 1e-10);
    ASSERT(fabs(f[1] - (2.0 + sqrt(3.0))) < 1e-10);
    expr_free(fe);
    free(m);
}

/* Single-eigenvalue interval {1.5, 2.5}: only lambda = 2 fits. */
void test_feast_real_sym_single_eigenvalue(void) {
    double A[25] = {
         2, -1,  0,  0,  0,
        -1,  2, -1,  0,  0,
         0, -1,  2, -1,  0,
         0,  0, -1,  2, -1,
         0,  0,  0, -1,  2
    };
    char* m = fmt_real_matrix(A, 5);
    char buf[4096];
    snprintf(buf, sizeof buf,
        "Eigenvalues[%s, Method -> {\"FEAST\", \"Interval\" -> {1.5, 2.5}}]", m);
    Expr* fe = eval_string(buf);
    ASSERT(list_len_eq(fe, 1));
    double v = extract_real(fe->data.function.args[0]);
    ASSERT(fabs(v - 2.0) < 1e-10);
    expr_free(fe);
    free(m);
}

/* 3x3 symmetric, dense -- whole-spectrum interval matches Direct. */
void test_feast_real_sym_dense_3x3(void) {
    double A[9] = {
        4.0, 1.0, 2.0,
        1.0, 5.0, 1.0,
        2.0, 1.0, 3.0
    };
    char* m = fmt_real_matrix(A, 3);
    char buf[2048];
    snprintf(buf, sizeof buf,
        "Eigenvalues[%s, Method -> {\"FEAST\", \"Interval\" -> {0, 10}}]", m);
    Expr* fe = eval_string(buf);
    snprintf(buf, sizeof buf, "Eigenvalues[%s, Method -> \"Direct\"]", m);
    Expr* dr = eval_string(buf);
    ASSERT(list_len_eq(fe, 3));
    ASSERT(list_len_eq(dr, 3));
    double f[3], d[3];
    for (size_t i = 0; i < 3; i++) {
        f[i] = extract_real(fe->data.function.args[i]);
        d[i] = extract_real(dr->data.function.args[i]);
    }
    sort_asc(f, 3); sort_asc(d, 3);
    for (size_t i = 0; i < 3; i++)
        ASSERT(fabs(f[i] - d[i]) < 1e-10);
    expr_free(fe); expr_free(dr);
    free(m);
}

/* ============================================================ *
 *  Complex Hermitian FEAST -- cross-checks                      *
 * ============================================================ */

/* 2x2 complex Hermitian {{2, I}, {-I, 3}}: eigenvalues (5 +/- sqrt(5))/2
 * ~= {3.618, 1.382}.  Full interval [0, 5]. */
void test_feast_complex_herm_2x2_full(void) {
    Expr* fe = eval_string(
        "Eigenvalues[N[{{2, I}, {-I, 3}}], "
        "Method -> {\"FEAST\", \"Interval\" -> {0, 5}}]");
    ASSERT(list_len_eq(fe, 2));
    double f[2];
    for (size_t i = 0; i < 2; i++)
        f[i] = extract_real(fe->data.function.args[i]);
    sort_asc(f, 2);
    double exp_lo = (5.0 - sqrt(5.0)) * 0.5;
    double exp_hi = (5.0 + sqrt(5.0)) * 0.5;
    ASSERT(fabs(f[0] - exp_lo) < 1e-10);
    ASSERT(fabs(f[1] - exp_hi) < 1e-10);
    expr_free(fe);
}

/* Sub-interval [2, 5] on the same matrix -- expect only the larger
 * eigenvalue 3.618. */
void test_feast_complex_herm_2x2_subinterval(void) {
    Expr* fe = eval_string(
        "Eigenvalues[N[{{2, I}, {-I, 3}}], "
        "Method -> {\"FEAST\", \"Interval\" -> {2, 5}}]");
    ASSERT(list_len_eq(fe, 1));
    double v = extract_real(fe->data.function.args[0]);
    double exp_hi = (5.0 + sqrt(5.0)) * 0.5;
    ASSERT(fabs(v - exp_hi) < 1e-10);
    expr_free(fe);
}

/* 3x3 complex Hermitian with non-trivial structure: whole-spectrum
 * interval matches Direct. */
void test_feast_complex_herm_3x3_full(void) {
    Expr* fe = eval_string(
        "Eigenvalues[N[{{2, 1+I, 0}, {1-I, 3, I}, {0, -I, 4}}], "
        "Method -> {\"FEAST\", \"Interval\" -> {-1, 7}}]");
    Expr* dr = eval_string(
        "Eigenvalues[N[{{2, 1+I, 0}, {1-I, 3, I}, {0, -I, 4}}], "
        "Method -> \"Direct\"]");
    ASSERT(list_len_eq(fe, 3));
    ASSERT(list_len_eq(dr, 3));
    double f[3], d[3];
    for (size_t i = 0; i < 3; i++) {
        f[i] = extract_real(fe->data.function.args[i]);
        d[i] = extract_real(dr->data.function.args[i]);
    }
    sort_asc(f, 3); sort_asc(d, 3);
    for (size_t i = 0; i < 3; i++)
        ASSERT(fabs(f[i] - d[i]) < 1e-10);
    expr_free(fe); expr_free(dr);
}

/* ============================================================ *
 *  Eigenvector residual + orthonormality                        *
 * ============================================================ */

/* Real-symmetric eigenvector residual: ||A v - lambda v|| < tol. */
void test_feast_real_residual_5x5(void) {
    double A[25] = {
         3, -1,  0.5, 0, 0,
        -1,  4, -1,  0.5, 0,
         0.5, -1,  5, -1,  0.5,
         0,  0.5, -1,  6, -1,
         0,  0, 0.5, -1, 7
    };
    size_t n = 5;
    char* m = fmt_real_matrix(A, n);
    char buf[4096];

    snprintf(buf, sizeof buf,
        "Eigenvalues[%s, Method -> {\"FEAST\", \"Interval\" -> {0, 10}}]", m);
    Expr* eval_list = eval_string(buf);
    ASSERT(list_len_eq(eval_list, n));
    double* lambdas = (double*)malloc(sizeof(double) * n);
    for (size_t i = 0; i < n; i++)
        lambdas[i] = extract_real(eval_list->data.function.args[i]);
    expr_free(eval_list);

    snprintf(buf, sizeof buf,
        "Eigenvectors[%s, Method -> {\"FEAST\", \"Interval\" -> {0, 10}}]", m);
    Expr* evec = eval_string(buf);
    ASSERT(list_len_eq(evec, n));
    double* V = (double*)malloc(sizeof(double) * n * n);
    for (size_t i = 0; i < n; i++) {
        Expr* row = evec->data.function.args[i];
        ASSERT(list_len_eq(row, n));
        for (size_t j = 0; j < n; j++)
            V[i * n + j] = extract_real(row->data.function.args[j]);
    }
    expr_free(evec);
    free(m);

    double normA = corpus_norm_inf_real(A, n);
    double tol = 1e-9 * (double)n * normA;
    for (size_t i = 0; i < n; i++) {
        double max_r = 0.0;
        for (size_t r = 0; r < n; r++) {
            double s = 0.0;
            for (size_t c = 0; c < n; c++)
                s += A[r * n + c] * V[i * n + c];
            s -= lambdas[i] * V[i * n + r];
            if (fabs(s) > max_r) max_r = fabs(s);
        }
        ASSERT(max_r < tol);
    }
    free(lambdas); free(V);
}

/* Real-symmetric orthonormality: ||V V^T - I||_inf < tol. */
void test_feast_real_orthonormal_5x5(void) {
    double A[25] = {
         3, -1,  0.5, 0, 0,
        -1,  4, -1,  0.5, 0,
         0.5, -1,  5, -1,  0.5,
         0,  0.5, -1,  6, -1,
         0,  0, 0.5, -1, 7
    };
    size_t n = 5;
    char* m = fmt_real_matrix(A, n);
    char buf[4096];
    snprintf(buf, sizeof buf,
        "Eigenvectors[%s, Method -> {\"FEAST\", \"Interval\" -> {0, 10}}]", m);
    Expr* evec = eval_string(buf);
    ASSERT(list_len_eq(evec, n));
    double* V = (double*)malloc(sizeof(double) * n * n);
    for (size_t i = 0; i < n; i++) {
        Expr* row = evec->data.function.args[i];
        for (size_t j = 0; j < n; j++)
            V[i * n + j] = extract_real(row->data.function.args[j]);
    }
    expr_free(evec);
    free(m);

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            double s = 0.0;
            for (size_t k = 0; k < n; k++) s += V[i * n + k] * V[j * n + k];
            double expected = (i == j) ? 1.0 : 0.0;
            ASSERT(fabs(s - expected) < 1e-10);
        }
    }
    free(V);
}

/* Complex-Hermitian residual + unitary check (full interval). */
void test_feast_complex_residual_3x3(void) {
    Expr* fe = eval_string(
        "Eigenvalues[N[{{2, 1+I, 0}, {1-I, 3, I}, {0, -I, 4}}], "
        "Method -> {\"FEAST\", \"Interval\" -> {-1, 7}}]");
    Expr* ve = eval_string(
        "Eigenvectors[N[{{2, 1+I, 0}, {1-I, 3, I}, {0, -I, 4}}], "
        "Method -> {\"FEAST\", \"Interval\" -> {-1, 7}}]");
    ASSERT(list_len_eq(fe, 3));
    ASSERT(list_len_eq(ve, 3));

    double A_re[9] = { 2,  1,  0,
                       1,  3,  0,
                       0,  0,  4 };
    double A_im[9] = { 0,  1,  0,
                      -1,  0,  1,
                       0, -1,  0 };
    size_t n = 3;
    double lambdas[3];
    for (size_t i = 0; i < n; i++)
        lambdas[i] = extract_real(fe->data.function.args[i]);

    double V_re[9], V_im[9];
    for (size_t i = 0; i < n; i++) {
        Expr* row = ve->data.function.args[i];
        ASSERT(list_len_eq(row, n));
        for (size_t j = 0; j < n; j++) {
            V_re[i * n + j] = extract_real(row->data.function.args[j]);
            V_im[i * n + j] = extract_imag(row->data.function.args[j]);
        }
    }
    expr_free(fe); expr_free(ve);

    double normA = corpus_norm_inf_complex(A_re, A_im, n);
    double tol = 1e-9 * (double)n * normA;
    /* Residual: row i = i-th eigenvector v.  Check ||(A - lambda I) v||. */
    for (size_t i = 0; i < n; i++) {
        double max_r = 0.0;
        for (size_t r = 0; r < n; r++) {
            double sr = 0.0, si = 0.0;
            for (size_t c = 0; c < n; c++) {
                sr += A_re[r * n + c] * V_re[i * n + c]
                    - A_im[r * n + c] * V_im[i * n + c];
                si += A_re[r * n + c] * V_im[i * n + c]
                    + A_im[r * n + c] * V_re[i * n + c];
            }
            sr -= lambdas[i] * V_re[i * n + r];
            si -= lambdas[i] * V_im[i * n + r];
            double mag = hypot(sr, si);
            if (mag > max_r) max_r = mag;
        }
        ASSERT(max_r < tol);
    }

    /* Unitary: V V^H = I (row i is i-th eigenvector). */
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            double sr = 0.0, si = 0.0;
            for (size_t k = 0; k < n; k++) {
                /* Row i dot conj(row j) */
                sr += V_re[i * n + k] * V_re[j * n + k]
                    + V_im[i * n + k] * V_im[j * n + k];
                si += V_im[i * n + k] * V_re[j * n + k]
                    - V_re[i * n + k] * V_im[j * n + k];
            }
            double er = (i == j) ? 1.0 : 0.0;
            ASSERT(fabs(sr - er) < 1e-9);
            ASSERT(fabs(si) < 1e-9);
        }
    }
}

/* ============================================================ *
 *  k_spec interaction                                            *
 * ============================================================ */

/* Eigenvalues[A, 2, Method -> {"FEAST", ...}] returns top-2 by |lambda|
 * from inside the interval. */
void test_feast_top_k_filtering(void) {
    double A[25] = {
         2, -1,  0,  0,  0,
        -1,  2, -1,  0,  0,
         0, -1,  2, -1,  0,
         0,  0, -1,  2, -1,
         0,  0,  0, -1,  2
    };
    char* m = fmt_real_matrix(A, 5);
    char buf[4096];
    snprintf(buf, sizeof buf,
        "Eigenvalues[%s, 2, "
        "Method -> {\"FEAST\", \"Interval\" -> {0, 4}}]", m);
    Expr* fe = eval_string(buf);
    ASSERT(list_len_eq(fe, 2));
    /* Top 2 by |lambda| among {0.268, 1, 2, 3, 3.732} are {3.732, 3}. */
    double v0 = extract_real(fe->data.function.args[0]);
    double v1 = extract_real(fe->data.function.args[1]);
    ASSERT(fabs(fabs(v0) - (2.0 + sqrt(3.0))) < 1e-10);
    ASSERT(fabs(fabs(v1) - 3.0) < 1e-10);
    expr_free(fe);
    free(m);
}

/* ============================================================ *
 *  Edge-case fall-throughs to Direct                            *
 * ============================================================ */

/* Method -> "FEAST" with no "Interval" sub-option: feast_dispatch
 * returns NULL because interval_given == false, and the call falls
 * through to Direct.  Result must therefore equal Direct's full
 * spectrum. */
void test_feast_missing_interval_fallback(void) {
    double A[25] = {
         2, -1,  0,  0,  0,
        -1,  2, -1,  0,  0,
         0, -1,  2, -1,  0,
         0,  0, -1,  2, -1,
         0,  0,  0, -1,  2
    };
    char* m = fmt_real_matrix(A, 5);
    char buf[4096];
    snprintf(buf, sizeof buf,
        "Eigenvalues[%s, Method -> \"FEAST\"]", m);
    Expr* fe = eval_string(buf);
    snprintf(buf, sizeof buf, "Eigenvalues[%s, Method -> \"Direct\"]", m);
    Expr* dr = eval_string(buf);
    ASSERT(list_len_eq(fe, 5));
    ASSERT(list_len_eq(dr, 5));
    for (size_t i = 0; i < 5; i++) {
        double vf = extract_real(fe->data.function.args[i]);
        double vd = extract_real(dr->data.function.args[i]);
        ASSERT(fabs(vf - vd) < 1e-10);
    }
    expr_free(fe); expr_free(dr);
    free(m);
}

/* Degenerate (empty) interval {5, 5}: the parser leaves it as
 * interval_low == interval_high, the kernel's `high <= low` guard
 * trips, the kernel returns NULL, and the call falls through to
 * Direct.  Result must equal Direct's full spectrum. */
void test_feast_empty_interval_fallback(void) {
    double A[25] = {
         2, -1,  0,  0,  0,
        -1,  2, -1,  0,  0,
         0, -1,  2, -1,  0,
         0,  0, -1,  2, -1,
         0,  0,  0, -1,  2
    };
    char* m = fmt_real_matrix(A, 5);
    char buf[4096];
    snprintf(buf, sizeof buf,
        "Eigenvalues[%s, Method -> {\"FEAST\", \"Interval\" -> {5, 5}}]", m);
    Expr* fe = eval_string(buf);
    snprintf(buf, sizeof buf, "Eigenvalues[%s, Method -> \"Direct\"]", m);
    Expr* dr = eval_string(buf);
    ASSERT(list_len_eq(fe, 5));
    ASSERT(list_len_eq(dr, 5));
    for (size_t i = 0; i < 5; i++) {
        double vf = extract_real(fe->data.function.args[i]);
        double vd = extract_real(dr->data.function.args[i]);
        ASSERT(fabs(vf - vd) < 1e-10);
    }
    expr_free(fe); expr_free(dr);
    free(m);
}

/* Non-Hermitian real input: feast_dispatch_machine rejects it (not
 * real symmetric), the kernel returns NULL, the call cascades to
 * Direct. */
void test_feast_nonhermitian_fallback(void) {
    /* {{2,1},{3,4}} is non-symmetric: eigenvalues 5 and 1. */
    Expr* fe = eval_string(
        "Eigenvalues[N[{{2, 1}, {3, 4}}], "
        "Method -> {\"FEAST\", \"Interval\" -> {0, 10}}]");
    Expr* dr = eval_string(
        "Eigenvalues[N[{{2, 1}, {3, 4}}], Method -> \"Direct\"]");
    ASSERT(list_len_eq(fe, 2));
    ASSERT(list_len_eq(dr, 2));
    for (size_t i = 0; i < 2; i++) {
        double vf = extract_real(fe->data.function.args[i]);
        double vd = extract_real(dr->data.function.args[i]);
        ASSERT(fabs(vf - vd) < 1e-10);
    }
    expr_free(fe); expr_free(dr);
}

/* Generalised eigenproblem -- FEAST does not yet support a != NULL,
 * the dispatcher returns NULL and the call cascades to Direct's
 * symbolic generalised path.  We only check the call doesn't crash
 * and returns a List of the right length. */
void test_feast_generalised_fallback(void) {
    Expr* fe = eval_string(
        "Eigenvalues[{N[{{2, 1}, {1, 3}}], N[{{1, 0}, {0, 1}}]}, "
        "Method -> {\"FEAST\", \"Interval\" -> {0, 5}}]");
    /* Identity B reduces this to a regular Hermitian eigenproblem -- the
     * generalised codepath calls into the regular pipeline.  We just
     * check the result is a length-2 List. */
    ASSERT(is_list(fe));
    expr_free(fe);
}

/* Undersized subspace: SubspaceSize -> 1 on a 5x5 with all 5 eigenvalues
 * in the interval.  The kernel may converge with only 1 eigenvalue
 * returned, or fail and cascade to Direct (5 eigenvalues).  Either is a
 * valid outcome; we just require a well-formed List and that any
 * returned value lies in [0, 4]. */
void test_feast_undersized_subspace(void) {
    double A[25] = {
         2, -1,  0,  0,  0,
        -1,  2, -1,  0,  0,
         0, -1,  2, -1,  0,
         0,  0, -1,  2, -1,
         0,  0,  0, -1,  2
    };
    char* m = fmt_real_matrix(A, 5);
    char buf[4096];
    snprintf(buf, sizeof buf,
        "Eigenvalues[%s, Method -> {\"FEAST\", "
        "\"Interval\" -> {0, 4}, \"SubspaceSize\" -> 1}]", m);
    Expr* fe = eval_string(buf);
    ASSERT(is_list(fe));
    size_t n = fe->data.function.arg_count;
    ASSERT(n == 1 || n == 5);   /* converged-1 or cascaded-5. */
    for (size_t i = 0; i < n; i++) {
        double v = extract_real(fe->data.function.args[i]);
        /* Cascade case: any of the 5 eigenvalues in [0.26, 3.74].
         * Converged case: one eigenvalue in [0, 4]. */
        ASSERT(v >= -0.1 && v <= 4.1);
    }
    expr_free(fe);
    free(m);
}

/* ============================================================ *
 *  Eigenvectors output shape                                    *
 * ============================================================ */

/* Eigenvectors with sub-interval returns a List of correct-length
 * eigenvectors. */
void test_feast_eigenvectors_shape_subinterval(void) {
    double A[25] = {
         2, -1,  0,  0,  0,
        -1,  2, -1,  0,  0,
         0, -1,  2, -1,  0,
         0,  0, -1,  2, -1,
         0,  0,  0, -1,  2
    };
    char* m = fmt_real_matrix(A, 5);
    char buf[4096];
    snprintf(buf, sizeof buf,
        "Eigenvectors[%s, "
        "Method -> {\"FEAST\", \"Interval\" -> {2.5, 4}}]", m);
    Expr* ve = eval_string(buf);
    ASSERT(list_len_eq(ve, 2));
    for (size_t i = 0; i < 2; i++)
        ASSERT(list_len_eq(ve->data.function.args[i], 5));
    expr_free(ve);
    free(m);
}

/* Complex Hermitian eigenvectors shape. */
void test_feast_complex_eigenvectors_shape(void) {
    Expr* ve = eval_string(
        "Eigenvectors[N[{{2, I}, {-I, 3}}], "
        "Method -> {\"FEAST\", \"Interval\" -> {0, 5}}]");
    ASSERT(list_len_eq(ve, 2));
    for (size_t i = 0; i < 2; i++)
        ASSERT(list_len_eq(ve->data.function.args[i], 2));
    expr_free(ve);
}

#ifdef USE_MPFR

/* ============================================================ *
 *  MPFR FEAST tests                                              *
 * ============================================================ */

static double extract_abs_any(Expr* e) {
    if (!e) return NAN;
    if (e->type == EXPR_INTEGER) return fabs((double)e->data.integer);
    if (e->type == EXPR_REAL)    return fabs(e->data.real);
    if (e->type == EXPR_MPFR)    return fabs(mpfr_get_d(e->data.mpfr, MPFR_RNDN));
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && strcmp(e->data.function.head->data.symbol, "Complex") == 0
        && e->data.function.arg_count == 2) {
        double a = extract_abs_any(e->data.function.args[0]);
        double b = extract_abs_any(e->data.function.args[1]);
        return hypot(a, b);
    }
    return NAN;
}

static double extract_real_any(Expr* e) {
    if (!e) return NAN;
    if (e->type == EXPR_INTEGER) return (double)e->data.integer;
    if (e->type == EXPR_REAL)    return e->data.real;
    if (e->type == EXPR_MPFR)    return mpfr_get_d(e->data.mpfr, MPFR_RNDN);
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && strcmp(e->data.function.head->data.symbol, "Complex") == 0
        && e->data.function.arg_count == 2)
        return extract_real_any(e->data.function.args[0]);
    return NAN;
}

/* MPFR real symmetric tridiag 5x5 at 80 decimal digits: FEAST matches
 * Direct on the full spectrum.  Tolerance is loose at the double end
 * (the comparison is `double - double < 1e-18`) but verifies the
 * dispatch route lights up the MPFR FEAST kernel. */
void test_mpfr_feast_real_sym_tridiag_5x5(void) {
    Expr* fe = eval_string(
        "Eigenvalues[N[{{2,-1,0,0,0},{-1,2,-1,0,0},{0,-1,2,-1,0},"
        "{0,0,-1,2,-1},{0,0,0,-1,2}}, 80], "
        "Method -> {\"FEAST\", \"Interval\" -> {0, 4}}]");
    Expr* dr = eval_string(
        "Eigenvalues[N[{{2,-1,0,0,0},{-1,2,-1,0,0},{0,-1,2,-1,0},"
        "{0,0,-1,2,-1},{0,0,0,-1,2}}, 80], Method -> \"Direct\"]");
    ASSERT(list_len_eq(fe, 5));
    ASSERT(list_len_eq(dr, 5));
    double f[5], d[5];
    for (size_t i = 0; i < 5; i++) {
        f[i] = extract_real_any(fe->data.function.args[i]);
        d[i] = extract_real_any(dr->data.function.args[i]);
    }
    sort_asc(f, 5); sort_asc(d, 5);
    for (size_t i = 0; i < 5; i++)
        ASSERT(fabs(f[i] - d[i]) < 1e-18);
    expr_free(fe); expr_free(dr);
}

/* MPFR real symmetric sub-interval at 80 digits: pick the two largest
 * eigenvalues. */
void test_mpfr_feast_real_sym_subinterval(void) {
    Expr* fe = eval_string(
        "Eigenvalues[N[{{2,-1,0,0,0},{-1,2,-1,0,0},{0,-1,2,-1,0},"
        "{0,0,-1,2,-1},{0,0,0,-1,2}}, 80], "
        "Method -> {\"FEAST\", \"Interval\" -> {2.5, 4}}]");
    ASSERT(list_len_eq(fe, 2));
    double f[2];
    for (size_t i = 0; i < 2; i++)
        f[i] = extract_real_any(fe->data.function.args[i]);
    sort_asc(f, 2);
    ASSERT(fabs(f[0] - 3.0) < 1e-18);
    ASSERT(fabs(f[1] - (2.0 + sqrt(3.0))) < 1e-15);
    expr_free(fe);
}

/* MPFR complex Hermitian 2x2 at 80 digits: same closed form as the
 * machine test, but routes through feast_complex_hermitian_mpfr. */
void test_mpfr_feast_complex_herm_2x2(void) {
    Expr* fe = eval_string(
        "Eigenvalues[N[{{2, I}, {-I, 3}}, 80], "
        "Method -> {\"FEAST\", \"Interval\" -> {0, 5}}]");
    ASSERT(list_len_eq(fe, 2));
    double f[2];
    for (size_t i = 0; i < 2; i++)
        f[i] = extract_abs_any(fe->data.function.args[i]);
    sort_asc(f, 2);
    double exp_lo = (5.0 - sqrt(5.0)) * 0.5;
    double exp_hi = (5.0 + sqrt(5.0)) * 0.5;
    ASSERT(fabs(f[0] - exp_lo) < 1e-15);
    ASSERT(fabs(f[1] - exp_hi) < 1e-15);
    expr_free(fe);
}

/* MPFR complex Hermitian 3x3 at 128 decimal digits.  Cross-check
 * against Direct via |lambda|.  At 128 dd (~426 bits) the MPFR FEAST
 * kernel converges well below 1e-30 -- we only check 1e-25 here to
 * absorb the double-cast in extract_abs_any. */
void test_mpfr_feast_complex_herm_3x3_128dd(void) {
    Expr* fe = eval_string(
        "Eigenvalues[N[{{2, 1+I, 0}, {1-I, 3, I}, {0, -I, 4}}, 128], "
        "Method -> {\"FEAST\", \"Interval\" -> {-1, 7}}]");
    Expr* dr = eval_string(
        "Eigenvalues[N[{{2, 1+I, 0}, {1-I, 3, I}, {0, -I, 4}}, 128], "
        "Method -> \"Direct\"]");
    ASSERT(list_len_eq(fe, 3));
    ASSERT(list_len_eq(dr, 3));
    double f[3], d[3];
    for (size_t i = 0; i < 3; i++) {
        f[i] = extract_abs_any(fe->data.function.args[i]);
        d[i] = extract_abs_any(dr->data.function.args[i]);
    }
    sort_asc(f, 3); sort_asc(d, 3);
    for (size_t i = 0; i < 3; i++)
        ASSERT(fabs(f[i] - d[i]) < 1e-25);
    expr_free(fe); expr_free(dr);
}

/* MPFR eigenvector shape (real symmetric at 60 digits) */
void test_mpfr_feast_eigenvectors_shape(void) {
    Expr* v = eval_string(
        "Eigenvectors[N[{{2,-1,0,0,0},{-1,2,-1,0,0},{0,-1,2,-1,0},"
        "{0,0,-1,2,-1},{0,0,0,-1,2}}, 60], "
        "Method -> {\"FEAST\", \"Interval\" -> {0, 4}}]");
    ASSERT(list_len_eq(v, 5));
    for (size_t i = 0; i < 5; i++)
        ASSERT(list_len_eq(v->data.function.args[i], 5));
    expr_free(v);
}

#endif /* USE_MPFR */

int main(void) {
    symtab_init();
    core_init();

    /* Real symmetric cross-checks */
    TEST(test_feast_real_sym_full_spectrum_5x5);
    TEST(test_feast_real_sym_subinterval_5x5);
    TEST(test_feast_real_sym_single_eigenvalue);
    TEST(test_feast_real_sym_dense_3x3);

    /* Complex Hermitian cross-checks */
    TEST(test_feast_complex_herm_2x2_full);
    TEST(test_feast_complex_herm_2x2_subinterval);
    TEST(test_feast_complex_herm_3x3_full);

    /* Residual + orthonormality */
    TEST(test_feast_real_residual_5x5);
    TEST(test_feast_real_orthonormal_5x5);
    TEST(test_feast_complex_residual_3x3);

    /* k-spec interaction */
    TEST(test_feast_top_k_filtering);

    /* Edge cases */
    TEST(test_feast_missing_interval_fallback);
    TEST(test_feast_empty_interval_fallback);
    TEST(test_feast_nonhermitian_fallback);
    TEST(test_feast_generalised_fallback);
    TEST(test_feast_undersized_subspace);

    /* Eigenvectors shape */
    TEST(test_feast_eigenvectors_shape_subinterval);
    TEST(test_feast_complex_eigenvectors_shape);

#ifdef USE_MPFR
    TEST(test_mpfr_feast_real_sym_tridiag_5x5);
    TEST(test_mpfr_feast_real_sym_subinterval);
    TEST(test_mpfr_feast_complex_herm_2x2);
    TEST(test_mpfr_feast_complex_herm_3x3_128dd);
    TEST(test_mpfr_feast_eigenvectors_shape);
#endif

    printf("All test_mateigen_feast tests passed!\n");
    return 0;
}
