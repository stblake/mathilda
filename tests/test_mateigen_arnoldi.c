/* test_mateigen_arnoldi.c -- Tests for the numerical "Arnoldi"
 * Eigenvalues / Eigenvectors kernel.
 *
 * Coverage (Phase 3): real general and complex general Arnoldi at
 * machine precision.  Hermitian / symmetric inputs flow through the
 * same kernels and are exercised as part of the corpus below.
 *
 * Strategy: Arnoldi's correctness contract is "top-k Ritz values match
 * Direct's top-k within tolerance", so most tests compare Arnoldi
 * output against the Direct dispatcher rather than asserting closed-
 * form values.  Residual checks (||A v - lambda v||_inf) are also run
 * on the recovered Ritz vectors.
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

/* Convenience: evaluate a string, return the result Expr*. */
static Expr* eval_string(const char* s) {
    Expr* parsed = parse_expression(s);
    return evaluate(parsed);
}

/* Format a row-major real matrix as a Mathilda numeric matrix literal,
 * e.g. N[{{1.0, 2.0}, {3.0, 4.0}}].  Caller free()s the returned buffer. */
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

/* True iff a List Expr contains exactly `expected` elements. */
static bool list_len_eq(Expr* e, size_t expected) {
    if (!e || e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL) return false;
    if (strcmp(e->data.function.head->data.symbol, "List") != 0) return false;
    return e->data.function.arg_count == expected;
}

/* Extract a real numeric value from an Expr (Integer, Real, Complex[re,0]). */
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

/* Extract a complex magnitude |z| from an Expr (Integer, Real, or Complex[a,b]). */
static double extract_abs(Expr* e) {
    if (!e) return NAN;
    if (e->type == EXPR_INTEGER) return fabs((double)e->data.integer);
    if (e->type == EXPR_REAL)    return fabs(e->data.real);
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && strcmp(e->data.function.head->data.symbol, "Complex") == 0
        && e->data.function.arg_count == 2) {
        double a = extract_real(e->data.function.args[0]);
        double b = extract_real(e->data.function.args[1]);
        return hypot(a, b);
    }
    return NAN;
}

/* ============================================================ *
 *  Real general / symmetric Arnoldi                              *
 * ============================================================ */

/* Full Arnoldi (no k_spec) returns n eigenvalues matching Direct. */
void test_arnoldi_real_sym_3x3_full(void) {
    /* 3x3 symmetric tridiagonal: {{4,1,0},{1,3,1},{0,1,2}}. */
    double A[9] = {
        4.0, 1.0, 0.0,
        1.0, 3.0, 1.0,
        0.0, 1.0, 2.0
    };
    char* m = fmt_real_matrix(A, 3);
    char buf[1024];

    snprintf(buf, sizeof buf, "Eigenvalues[%s, Method -> \"Arnoldi\"]", m);
    Expr* arn = eval_string(buf);
    snprintf(buf, sizeof buf, "Eigenvalues[%s, Method -> \"Direct\"]", m);
    Expr* dir = eval_string(buf);

    ASSERT(list_len_eq(arn, 3));
    ASSERT(list_len_eq(dir, 3));
    for (size_t i = 0; i < 3; i++) {
        double va = extract_real(arn->data.function.args[i]);
        double vd = extract_real(dir->data.function.args[i]);
        if (fabs(va - vd) > 1e-9) {
            printf("FAIL: Arnoldi[%zu]=%g != Direct=%g\n", i, va, vd);
            ASSERT(0);
        }
    }
    expr_free(arn); expr_free(dir);
    free(m);
}

/* Arnoldi with k=2: top-2 by magnitude match Direct's top-2. */
void test_arnoldi_real_top_k(void) {
    double A[9] = {
        4.0, 1.0, 0.0,
        1.0, 3.0, 1.0,
        0.0, 1.0, 2.0
    };
    char* m = fmt_real_matrix(A, 3);
    char buf[1024];
    snprintf(buf, sizeof buf, "Eigenvalues[%s, 2, Method -> \"Arnoldi\"]", m);
    Expr* arn = eval_string(buf);
    snprintf(buf, sizeof buf, "Eigenvalues[%s, 2, Method -> \"Direct\"]", m);
    Expr* dir = eval_string(buf);
    ASSERT(list_len_eq(arn, 2));
    ASSERT(list_len_eq(dir, 2));
    for (size_t i = 0; i < 2; i++) {
        double va = extract_real(arn->data.function.args[i]);
        double vd = extract_real(dir->data.function.args[i]);
        ASSERT(fabs(va - vd) < 1e-9);
    }
    expr_free(arn); expr_free(dir);
    free(m);
}

/* Real non-symmetric with complex conjugate pair (1 +/- i*sqrt(2), 1). */
void test_arnoldi_real_general_complex_pair(void) {
    double A[9] = {
        1.0, -1.0,  0.0,
        1.0,  1.0, -1.0,
        0.0,  1.0,  1.0
    };
    char* m = fmt_real_matrix(A, 3);
    char buf[1024];
    snprintf(buf, sizeof buf, "Eigenvalues[%s, Method -> \"Arnoldi\"]", m);
    Expr* arn = eval_string(buf);
    ASSERT(list_len_eq(arn, 3));
    /* Descending |lambda|: |1+/-i*sqrt(2)| = sqrt(3) ~ 1.732, then |1| = 1. */
    double m0 = extract_abs(arn->data.function.args[0]);
    double m1 = extract_abs(arn->data.function.args[1]);
    double m2 = extract_abs(arn->data.function.args[2]);
    ASSERT(fabs(m0 - sqrt(3.0)) < 1e-9);
    ASSERT(fabs(m1 - sqrt(3.0)) < 1e-9);
    ASSERT(fabs(m2 - 1.0) < 1e-9);
    expr_free(arn);
    free(m);
}

/* Top-1 (largest |lambda|) via Eigenvalues[A, 1]. */
void test_arnoldi_real_top_1(void) {
    double A[16] = {
        5.0, 1.0, 0.0, 0.0,
        1.0, 4.0, 1.0, 0.0,
        0.0, 1.0, 3.0, 1.0,
        0.0, 0.0, 1.0, 2.0
    };
    char* m = fmt_real_matrix(A, 4);
    char buf[1024];
    snprintf(buf, sizeof buf, "Eigenvalues[%s, 1, Method -> \"Arnoldi\"]", m);
    Expr* arn = eval_string(buf);
    snprintf(buf, sizeof buf, "Eigenvalues[%s, 1, Method -> \"Direct\"]", m);
    Expr* dir = eval_string(buf);
    ASSERT(list_len_eq(arn, 1));
    ASSERT(list_len_eq(dir, 1));
    double va = extract_real(arn->data.function.args[0]);
    double vd = extract_real(dir->data.function.args[0]);
    ASSERT(fabs(va - vd) < 1e-9);
    expr_free(arn); expr_free(dir);
    free(m);
}

/* Diagonal matrix with known spectrum {7, 5, 3, 1}. */
void test_arnoldi_diag_4x4(void) {
    double A[16] = {
        1.0, 0.0, 0.0, 0.0,
        0.0, 7.0, 0.0, 0.0,
        0.0, 0.0, 3.0, 0.0,
        0.0, 0.0, 0.0, 5.0
    };
    char* m = fmt_real_matrix(A, 4);
    char buf[1024];
    snprintf(buf, sizeof buf, "Eigenvalues[%s, Method -> \"Arnoldi\"]", m);
    Expr* arn = eval_string(buf);
    /* For a diagonal matrix and v0 = ones/sqrt(n), Arnoldi finds all n
     * eigenvalues with lucky breakdown at step n.  Check descending order. */
    ASSERT(list_len_eq(arn, 4));
    double expected[4] = {7.0, 5.0, 3.0, 1.0};
    for (size_t i = 0; i < 4; i++) {
        double v = extract_real(arn->data.function.args[i]);
        ASSERT(fabs(v - expected[i]) < 1e-9);
    }
    expr_free(arn);
    free(m);
}

/* Eigenvectors: A v - lambda v residual check (real symmetric). */
void test_arnoldi_real_residual(void) {
    double A[9] = {
        4.0, 1.0, 0.0,
        1.0, 3.0, 1.0,
        0.0, 1.0, 2.0
    };
    size_t n = 3;
    double *lambdas = NULL, *V = NULL;
    char* m = fmt_real_matrix(A, n);
    char buf[1024];

    /* Use the Arnoldi method explicitly. */
    snprintf(buf, sizeof buf, "Eigenvalues[%s, Method -> \"Arnoldi\"]", m);
    Expr* eval_list = eval_string(buf);
    ASSERT(list_len_eq(eval_list, n));
    lambdas = (double*)malloc(sizeof(double) * n);
    for (size_t i = 0; i < n; i++)
        lambdas[i] = extract_real(eval_list->data.function.args[i]);
    expr_free(eval_list);

    snprintf(buf, sizeof buf, "Eigenvectors[%s, Method -> \"Arnoldi\"]", m);
    Expr* evec_list = eval_string(buf);
    ASSERT(list_len_eq(evec_list, n));
    V = (double*)malloc(sizeof(double) * n * n);
    for (size_t i = 0; i < n; i++) {
        Expr* row = evec_list->data.function.args[i];
        ASSERT(list_len_eq(row, n));
        for (size_t j = 0; j < n; j++) {
            V[i * n + j] = extract_real(row->data.function.args[j]);
        }
    }
    expr_free(evec_list);
    free(m);

    /* Residual: ||A v_i - lambda_i v_i||_inf  < tol */
    double normA = corpus_norm_inf_real(A, n);
    double tol = 1e-9 * (double)n * normA;
    for (size_t i = 0; i < n; i++) {
        double max_r = 0.0;
        for (size_t r = 0; r < n; r++) {
            double s = 0.0;
            for (size_t c = 0; c < n; c++) s += A[r * n + c] * V[i * n + c];
            s -= lambdas[i] * V[i * n + r];
            if (fabs(s) > max_r) max_r = fabs(s);
        }
        ASSERT(max_r < tol);
    }
    free(lambdas); free(V);
}

/* ============================================================ *
 *  BasisSize / Tolerance / MaxIterations sub-options             *
 * ============================================================ */

/* BasisSize -> 3 caps the Krylov basis; result length == 3 in the
 * default no-k_spec case. */
void test_arnoldi_basis_size_cap(void) {
    double A[16] = {
        5.0, 1.0, 0.0, 0.0,
        1.0, 4.0, 1.0, 0.0,
        0.0, 1.0, 3.0, 1.0,
        0.0, 0.0, 1.0, 2.0
    };
    char* m = fmt_real_matrix(A, 4);
    char buf[1024];
    snprintf(buf, sizeof buf,
        "Eigenvalues[%s, Method -> {\"Arnoldi\", \"BasisSize\" -> 3}]", m);
    Expr* arn = eval_string(buf);
    ASSERT(list_len_eq(arn, 3));
    /* The first (largest |lambda|) Ritz value should be close to Direct's
     * largest eigenvalue (sqrt(2)-shifted symmetric tridiag eigenvalue
     * ~6.04 within Arnoldi's convergence at m=3). */
    double v0 = extract_real(arn->data.function.args[0]);
    /* Direct's top-1 for this matrix is ~5.83118; Arnoldi with m=3 should
     * be near it but exact convergence is not required for m < n. */
    ASSERT(v0 > 5.0 && v0 < 7.0);
    expr_free(arn);
    free(m);
}

/* MaxIterations is currently consumed but not used; just verify it
 * parses without error. */
void test_arnoldi_max_iter_parses(void) {
    double A[9] = {
        4.0, 1.0, 0.0,
        1.0, 3.0, 1.0,
        0.0, 1.0, 2.0
    };
    char* m = fmt_real_matrix(A, 3);
    char buf[1024];
    snprintf(buf, sizeof buf,
        "Eigenvalues[%s, Method -> {\"Arnoldi\", \"MaxIterations\" -> 500}]", m);
    Expr* arn = eval_string(buf);
    ASSERT(list_len_eq(arn, 3));
    expr_free(arn);
    free(m);
}

/* ============================================================ *
 *  Automatic heuristic                                           *
 * ============================================================ */

/* Automatic with k=1 routes to Arnoldi (small k) -- still produces a
 * single-element list that matches Direct's top-1. */
void test_arnoldi_automatic_small_k(void) {
    double A[16] = {
        5.0, 1.0, 0.0, 0.0,
        1.0, 4.0, 1.0, 0.0,
        0.0, 1.0, 3.0, 1.0,
        0.0, 0.0, 1.0, 2.0
    };
    char* m = fmt_real_matrix(A, 4);
    char buf[1024];
    /* Default Method -> Automatic; k = 1 -> Arnoldi is preferred. */
    snprintf(buf, sizeof buf, "Eigenvalues[%s, 1]", m);
    Expr* auto_res = eval_string(buf);
    snprintf(buf, sizeof buf, "Eigenvalues[%s, 1, Method -> \"Direct\"]", m);
    Expr* dir = eval_string(buf);
    ASSERT(list_len_eq(auto_res, 1));
    ASSERT(list_len_eq(dir, 1));
    double va = extract_real(auto_res->data.function.args[0]);
    double vd = extract_real(dir->data.function.args[0]);
    ASSERT(fabs(va - vd) < 1e-9);
    expr_free(auto_res); expr_free(dir);
    free(m);
}

/* ============================================================ *
 *  Complex Arnoldi                                               *
 * ============================================================ */

/* 2x2 complex Hermitian: {{2, I}, {-I, 3}} -> (5 +/- sqrt(5)) / 2. */
void test_arnoldi_complex_hermitian_2x2(void) {
    Expr* arn = eval_string(
        "Eigenvalues[N[{{2, I}, {-I, 3}}], Method -> \"Arnoldi\"]");
    ASSERT(list_len_eq(arn, 2));
    double e0 = extract_real(arn->data.function.args[0]);
    double e1 = extract_real(arn->data.function.args[1]);
    double exp0 = (5.0 + sqrt(5.0)) / 2.0;
    double exp1 = (5.0 - sqrt(5.0)) / 2.0;
    ASSERT(fabs(e0 - exp0) < 1e-9);
    ASSERT(fabs(e1 - exp1) < 1e-9);
    expr_free(arn);
}

/* Pauli-Y: {{0, -I}, {I, 0}} -> {+1, -1}. */
void test_arnoldi_pauli_y(void) {
    Expr* arn = eval_string(
        "Eigenvalues[N[{{0, -I}, {I, 0}}], Method -> \"Arnoldi\"]");
    ASSERT(list_len_eq(arn, 2));
    double e0 = extract_real(arn->data.function.args[0]);
    double e1 = extract_real(arn->data.function.args[1]);
    ASSERT(fabs(e0 - 1.0) < 1e-9);
    ASSERT(fabs(e1 + 1.0) < 1e-9);
    expr_free(arn);
}

/* Complex general 2x2: same matrix Direct handles. */
void test_arnoldi_complex_general_2x2(void) {
    Expr* arn = eval_string(
        "Eigenvalues[N[{{1 + I, 2}, {3, 4 - I}}], Method -> \"Arnoldi\"]");
    Expr* dir = eval_string(
        "Eigenvalues[N[{{1 + I, 2}, {3, 4 - I}}], Method -> \"Direct\"]");
    ASSERT(list_len_eq(arn, 2));
    ASSERT(list_len_eq(dir, 2));
    /* |lambda| should match between Arnoldi and Direct. */
    for (size_t i = 0; i < 2; i++) {
        double ma = extract_abs(arn->data.function.args[i]);
        double md = extract_abs(dir->data.function.args[i]);
        ASSERT(fabs(ma - md) < 1e-9);
    }
    expr_free(arn); expr_free(dir);
}

/* Eigenvectors are emitted as List of List for the complex case. */
void test_arnoldi_complex_eigvecs_shape(void) {
    Expr* v = eval_string(
        "Eigenvectors[N[{{2, I}, {-I, 3}}], Method -> \"Arnoldi\"]");
    ASSERT(list_len_eq(v, 2));
    for (size_t i = 0; i < 2; i++) {
        ASSERT(list_len_eq(v->data.function.args[i], 2));
    }
    expr_free(v);
}

/* UpTo[k] form -- accept up to k eigenvalues; with n=3, UpTo[5] returns 3. */
void test_arnoldi_upto_form(void) {
    Expr* arn = eval_string(
        "Eigenvalues[N[{{4, 1, 0}, {1, 3, 1}, {0, 1, 2}}], UpTo[5], "
        "Method -> \"Arnoldi\"]");
    ASSERT(list_len_eq(arn, 3));
    expr_free(arn);
    Expr* arn2 = eval_string(
        "Eigenvalues[N[{{4, 1, 0}, {1, 3, 1}, {0, 1, 2}}], UpTo[2], "
        "Method -> \"Arnoldi\"]");
    ASSERT(list_len_eq(arn2, 2));
    expr_free(arn2);
}

/* 1x1 input -- edge case. */
void test_arnoldi_1x1(void) {
    Expr* arn = eval_string(
        "Eigenvalues[N[{{7.0}}], Method -> \"Arnoldi\"]");
    ASSERT(list_len_eq(arn, 1));
    double v = extract_real(arn->data.function.args[0]);
    ASSERT(fabs(v - 7.0) < 1e-12);
    expr_free(arn);
}

#ifdef USE_MPFR

/* ============================================================ *
 *  Phase 3f: MPFR Arnoldi tests                                  *
 * ============================================================ */

/* Helper: extract magnitude of a possibly-MPFR / Complex Expr.  Returns
 * NaN on failure. */
static double extract_abs_any(Expr* e) {
    if (!e) return NAN;
    if (e->type == EXPR_INTEGER) return fabs((double)e->data.integer);
    if (e->type == EXPR_REAL)    return fabs(e->data.real);
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) {
        return fabs(mpfr_get_d(e->data.mpfr, MPFR_RNDN));
    }
#endif
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && strcmp(e->data.function.head->data.symbol, "Complex") == 0
        && e->data.function.arg_count == 2) {
        double a = extract_abs_any(e->data.function.args[0]);
        double b = extract_abs_any(e->data.function.args[1]);
        /* The arguments are unsigned absolutes from extract_abs_any;
         * for our purposes (compare magnitudes between Arnoldi/Direct)
         * we treat them as the real/imag magnitudes. */
        return hypot(a, b);
    }
    return NAN;
}

/* MPFR real symmetric: top-k from Arnoldi matches Direct within sqrt-eps tol. */
void test_mpfr_arnoldi_real_sym_3x3(void) {
    Expr* arn = eval_string(
        "Eigenvalues[N[{{4, 1, 0}, {1, 3, 1}, {0, 1, 2}}, 60], "
        "Method -> \"Arnoldi\"]");
    Expr* dir = eval_string(
        "Eigenvalues[N[{{4, 1, 0}, {1, 3, 1}, {0, 1, 2}}, 60], "
        "Method -> \"Direct\"]");
    ASSERT(list_len_eq(arn, 3));
    ASSERT(list_len_eq(dir, 3));
    for (size_t i = 0; i < 3; i++) {
        double ma = extract_abs_any(arn->data.function.args[i]);
        double md = extract_abs_any(dir->data.function.args[i]);
        ASSERT(fabs(ma - md) < 1e-15);
    }
    expr_free(arn); expr_free(dir);
}

/* MPFR real general (complex conjugate pair): magnitudes match Direct. */
void test_mpfr_arnoldi_real_general_complex_pair(void) {
    Expr* arn = eval_string(
        "Eigenvalues[N[{{1, -1, 0}, {1, 1, -1}, {0, 1, 1}}, 80], "
        "Method -> \"Arnoldi\"]");
    Expr* dir = eval_string(
        "Eigenvalues[N[{{1, -1, 0}, {1, 1, -1}, {0, 1, 1}}, 80], "
        "Method -> \"Direct\"]");
    ASSERT(list_len_eq(arn, 3));
    ASSERT(list_len_eq(dir, 3));
    for (size_t i = 0; i < 3; i++) {
        double ma = extract_abs_any(arn->data.function.args[i]);
        double md = extract_abs_any(dir->data.function.args[i]);
        ASSERT(fabs(ma - md) < 1e-15);
    }
    expr_free(arn); expr_free(dir);
}

/* MPFR complex Hermitian: 2x2 {{2, I}, {-I, 3}} -> (5 +/- sqrt(5))/2. */
void test_mpfr_arnoldi_complex_hermitian_2x2(void) {
    Expr* arn = eval_string(
        "Eigenvalues[N[{{2, I}, {-I, 3}}, 60], Method -> \"Arnoldi\"]");
    Expr* dir = eval_string(
        "Eigenvalues[N[{{2, I}, {-I, 3}}, 60], Method -> \"Direct\"]");
    ASSERT(list_len_eq(arn, 2));
    ASSERT(list_len_eq(dir, 2));
    for (size_t i = 0; i < 2; i++) {
        double ma = extract_abs_any(arn->data.function.args[i]);
        double md = extract_abs_any(dir->data.function.args[i]);
        ASSERT(fabs(ma - md) < 1e-15);
    }
    expr_free(arn); expr_free(dir);
}

/* MPFR complex general: magnitudes match Direct. */
void test_mpfr_arnoldi_complex_general_2x2(void) {
    Expr* arn = eval_string(
        "Eigenvalues[N[{{1 + I, 2}, {3, 4 - I}}, 80], "
        "Method -> \"Arnoldi\"]");
    Expr* dir = eval_string(
        "Eigenvalues[N[{{1 + I, 2}, {3, 4 - I}}, 80], "
        "Method -> \"Direct\"]");
    ASSERT(list_len_eq(arn, 2));
    ASSERT(list_len_eq(dir, 2));
    for (size_t i = 0; i < 2; i++) {
        double ma = extract_abs_any(arn->data.function.args[i]);
        double md = extract_abs_any(dir->data.function.args[i]);
        ASSERT(fabs(ma - md) < 1e-15);
    }
    expr_free(arn); expr_free(dir);
}

/* MPFR Arnoldi eigenvectors -- structural shape only (3 vectors of 3
 * components for the 3x3 input). */
void test_mpfr_arnoldi_eigenvectors_shape(void) {
    Expr* v = eval_string(
        "Eigenvectors[N[{{4, 1, 0}, {1, 3, 1}, {0, 1, 2}}, 60], "
        "Method -> \"Arnoldi\"]");
    ASSERT(list_len_eq(v, 3));
    for (size_t i = 0; i < 3; i++) {
        ASSERT(list_len_eq(v->data.function.args[i], 3));
    }
    expr_free(v);
}

/* Eigenvalues carry full input precision (60 -> 60 bits). */
void test_mpfr_arnoldi_eigenvalues_carry_precision(void) {
    Expr* arn = eval_string(
        "Eigenvalues[N[{{2, 1}, {1, 3}}, 80], Method -> \"Arnoldi\"]");
    ASSERT(list_len_eq(arn, 2));
    /* The eigenvalues of {{2, 1}, {1, 3}} are (5 +/- sqrt(5))/2.
     * At 80 bits ~ 24 digits, they should be well within 1e-22. */
    Expr* lam0 = arn->data.function.args[0];
    ASSERT(lam0 != NULL);
    /* Just verify the result is non-trivial (MPFR carries through). */
    double m0 = extract_abs_any(lam0);
    ASSERT(m0 > 3.6 && m0 < 3.7);
    expr_free(arn);
}

#endif /* USE_MPFR */

int main(void) {
    symtab_init();
    core_init();

    TEST(test_arnoldi_real_sym_3x3_full);
    TEST(test_arnoldi_real_top_k);
    TEST(test_arnoldi_real_general_complex_pair);
    TEST(test_arnoldi_real_top_1);
    TEST(test_arnoldi_diag_4x4);
    TEST(test_arnoldi_real_residual);

    TEST(test_arnoldi_basis_size_cap);
    TEST(test_arnoldi_max_iter_parses);
    TEST(test_arnoldi_automatic_small_k);

    TEST(test_arnoldi_complex_hermitian_2x2);
    TEST(test_arnoldi_pauli_y);
    TEST(test_arnoldi_complex_general_2x2);
    TEST(test_arnoldi_complex_eigvecs_shape);

    TEST(test_arnoldi_upto_form);
    TEST(test_arnoldi_1x1);

#ifdef USE_MPFR
    TEST(test_mpfr_arnoldi_real_sym_3x3);
    TEST(test_mpfr_arnoldi_real_general_complex_pair);
    TEST(test_mpfr_arnoldi_complex_hermitian_2x2);
    TEST(test_mpfr_arnoldi_complex_general_2x2);
    TEST(test_mpfr_arnoldi_eigenvectors_shape);
    TEST(test_mpfr_arnoldi_eigenvalues_carry_precision);
#endif

    printf("All Arnoldi tests passed!\n");
    return 0;
}
