/* test_mateigen_banded.c -- Tests for the numerical "Banded"
 * Eigenvalues / Eigenvectors kernel (Phase 4).
 *
 * Banded is Hermitian-only: it reduces a real-symmetric / complex-
 * Hermitian banded matrix to symmetric tridiagonal via Givens with
 * bulge chasing, then runs the existing Phase 2 symmetric tridiag QR.
 *
 * Correctness contract: Banded must match Direct's spectrum (and
 * orthonormal eigenvector basis) within tolerance for every Hermitian
 * input, regardless of bandwidth.  Non-Hermitian / dense inputs make
 * the dispatcher return NULL and the call falls through to Direct --
 * we exercise that fallback transparently.
 *
 * Tests:
 *   - real symmetric tridiagonal / pentadiagonal / heptadiagonal
 *   - diagonal (b = 0)
 *   - dense fallback (b = n - 1)
 *   - complex Hermitian tridiagonal / banded
 *   - eigenvector residual + orthonormality
 *   - k-spec filtering
 *   - 1x1 edge case
 *   - MPFR variants
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

/* ---------- formatting / extraction helpers (mirrors test_mateigen_arnoldi.c) */

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

/* ============================================================ *
 *  Real symmetric Banded                                         *
 * ============================================================ */

/* 5x5 symmetric tridiagonal -- standard Chebyshev spectrum
 * 2 - 2 cos(k pi / 6), k = 1..5 = {3.732, 3, 2, 1, 0.268}. */
void test_banded_real_sym_tridiag_5x5(void) {
    double A[25] = {
         2, -1,  0,  0,  0,
        -1,  2, -1,  0,  0,
         0, -1,  2, -1,  0,
         0,  0, -1,  2, -1,
         0,  0,  0, -1,  2
    };
    char* m = fmt_real_matrix(A, 5);
    char buf[2048];
    snprintf(buf, sizeof buf, "Eigenvalues[%s, Method -> \"Banded\"]", m);
    Expr* bn = eval_string(buf);
    snprintf(buf, sizeof buf, "Eigenvalues[%s, Method -> \"Direct\"]", m);
    Expr* dr = eval_string(buf);
    ASSERT(list_len_eq(bn, 5));
    ASSERT(list_len_eq(dr, 5));
    for (size_t i = 0; i < 5; i++) {
        double vb = extract_real(bn->data.function.args[i]);
        double vd = extract_real(dr->data.function.args[i]);
        ASSERT(fabs(vb - vd) < 1e-11);
    }
    expr_free(bn); expr_free(dr);
    free(m);
}

/* 6x6 symmetric pentadiagonal -- b = 2. */
void test_banded_real_sym_pentadiag_6x6(void) {
    double A[36] = {
        4.0, 1.0, 0.5, 0, 0, 0,
        1.0, 4.0, 1.0, 0.5, 0, 0,
        0.5, 1.0, 4.0, 1.0, 0.5, 0,
        0, 0.5, 1.0, 4.0, 1.0, 0.5,
        0, 0, 0.5, 1.0, 4.0, 1.0,
        0, 0, 0, 0.5, 1.0, 4.0
    };
    char* m = fmt_real_matrix(A, 6);
    char buf[2048];
    snprintf(buf, sizeof buf, "Eigenvalues[%s, Method -> \"Banded\"]", m);
    Expr* bn = eval_string(buf);
    snprintf(buf, sizeof buf, "Eigenvalues[%s, Method -> \"Direct\"]", m);
    Expr* dr = eval_string(buf);
    ASSERT(list_len_eq(bn, 6));
    ASSERT(list_len_eq(dr, 6));
    for (size_t i = 0; i < 6; i++) {
        double vb = extract_real(bn->data.function.args[i]);
        double vd = extract_real(dr->data.function.args[i]);
        ASSERT(fabs(vb - vd) < 1e-11);
    }
    expr_free(bn); expr_free(dr);
    free(m);
}

/* 7x7 symmetric heptadiagonal -- b = 3.  Stress test for the chase
 * loop because q + b is large enough for multiple bulge propagations
 * per outer step. */
void test_banded_real_sym_heptadiag_7x7(void) {
    double A[49] = {
        5.0, 1.0, 0.5, 0.25, 0, 0, 0,
        1.0, 5.0, 1.0, 0.5, 0.25, 0, 0,
        0.5, 1.0, 5.0, 1.0, 0.5, 0.25, 0,
        0.25, 0.5, 1.0, 5.0, 1.0, 0.5, 0.25,
        0, 0.25, 0.5, 1.0, 5.0, 1.0, 0.5,
        0, 0, 0.25, 0.5, 1.0, 5.0, 1.0,
        0, 0, 0, 0.25, 0.5, 1.0, 5.0
    };
    char* m = fmt_real_matrix(A, 7);
    char buf[4096];
    snprintf(buf, sizeof buf, "Eigenvalues[%s, Method -> \"Banded\"]", m);
    Expr* bn = eval_string(buf);
    snprintf(buf, sizeof buf, "Eigenvalues[%s, Method -> \"Direct\"]", m);
    Expr* dr = eval_string(buf);
    ASSERT(list_len_eq(bn, 7));
    ASSERT(list_len_eq(dr, 7));
    for (size_t i = 0; i < 7; i++) {
        double vb = extract_real(bn->data.function.args[i]);
        double vd = extract_real(dr->data.function.args[i]);
        ASSERT(fabs(vb - vd) < 1e-10);
    }
    expr_free(bn); expr_free(dr);
    free(m);
}

/* Diagonal 4x4 (b = 0): no Givens rotations are applied; result must
 * still come back sorted by |lambda| descending. */
void test_banded_diag_4x4(void) {
    double A[16] = {
        5.0, 0, 0, 0,
        0, 2.0, 0, 0,
        0, 0, -3.0, 0,
        0, 0, 0, 1.0
    };
    char* m = fmt_real_matrix(A, 4);
    char buf[1024];
    snprintf(buf, sizeof buf, "Eigenvalues[%s, Method -> \"Banded\"]", m);
    Expr* bn = eval_string(buf);
    ASSERT(list_len_eq(bn, 4));
    /* Descending |lambda|: 5, -3, 2, 1 (sorted by abs value). */
    double exp[4] = {5.0, -3.0, 2.0, 1.0};
    for (size_t i = 0; i < 4; i++) {
        double v = extract_real(bn->data.function.args[i]);
        ASSERT(fabs(v - exp[i]) < 1e-12);
    }
    expr_free(bn);
    free(m);
}

/* Dense 3x3 (b = n - 1 = 2): banded_dispatch refuses, the wrapper
 * falls through to Direct.  Result must therefore equal Direct. */
void test_banded_dense_fallback_3x3(void) {
    double A[9] = {
        1.0, 2.0, 3.0,
        2.0, 4.0, 5.0,
        3.0, 5.0, 6.0
    };
    char* m = fmt_real_matrix(A, 3);
    char buf[1024];
    snprintf(buf, sizeof buf, "Eigenvalues[%s, Method -> \"Banded\"]", m);
    Expr* bn = eval_string(buf);
    snprintf(buf, sizeof buf, "Eigenvalues[%s, Method -> \"Direct\"]", m);
    Expr* dr = eval_string(buf);
    ASSERT(list_len_eq(bn, 3));
    ASSERT(list_len_eq(dr, 3));
    for (size_t i = 0; i < 3; i++) {
        double vb = extract_real(bn->data.function.args[i]);
        double vd = extract_real(dr->data.function.args[i]);
        ASSERT(fabs(vb - vd) < 1e-11);
    }
    expr_free(bn); expr_free(dr);
    free(m);
}

/* Non-symmetric real input: Banded refuses, falls back to Direct. */
void test_banded_nonsym_fallback(void) {
    double A[9] = {
        1.0, 2.0, 0,
        0,   3.0, 4.0,
        0,   0,   5.0
    };
    char* m = fmt_real_matrix(A, 3);
    char buf[1024];
    snprintf(buf, sizeof buf, "Eigenvalues[%s, Method -> \"Banded\"]", m);
    Expr* bn = eval_string(buf);
    snprintf(buf, sizeof buf, "Eigenvalues[%s, Method -> \"Direct\"]", m);
    Expr* dr = eval_string(buf);
    ASSERT(list_len_eq(bn, 3));
    ASSERT(list_len_eq(dr, 3));
    for (size_t i = 0; i < 3; i++) {
        double vb = extract_real(bn->data.function.args[i]);
        double vd = extract_real(dr->data.function.args[i]);
        ASSERT(fabs(vb - vd) < 1e-11);
    }
    expr_free(bn); expr_free(dr);
    free(m);
}

/* k-spec filtering: Eigenvalues[A, 2] returns the 2 largest. */
void test_banded_top_k(void) {
    double A[25] = {
         2, -1,  0,  0,  0,
        -1,  2, -1,  0,  0,
         0, -1,  2, -1,  0,
         0,  0, -1,  2, -1,
         0,  0,  0, -1,  2
    };
    char* m = fmt_real_matrix(A, 5);
    char buf[2048];
    snprintf(buf, sizeof buf, "Eigenvalues[%s, 2, Method -> \"Banded\"]", m);
    Expr* bn = eval_string(buf);
    snprintf(buf, sizeof buf, "Eigenvalues[%s, 2, Method -> \"Direct\"]", m);
    Expr* dr = eval_string(buf);
    ASSERT(list_len_eq(bn, 2));
    ASSERT(list_len_eq(dr, 2));
    for (size_t i = 0; i < 2; i++) {
        double vb = extract_real(bn->data.function.args[i]);
        double vd = extract_real(dr->data.function.args[i]);
        ASSERT(fabs(vb - vd) < 1e-11);
    }
    expr_free(bn); expr_free(dr);
    free(m);
}

/* Eigenvector residual: ||A v - lambda v||_inf < n * eps * ||A||_inf. */
void test_banded_real_residual(void) {
    double A[25] = {
         3, -1,  0.5, 0, 0,
        -1,  4, -1,  0.5, 0,
         0.5, -1,  5, -1,  0.5,
         0,  0.5, -1,  6, -1,
         0,  0, 0.5, -1, 7
    };
    size_t n = 5;
    char* m = fmt_real_matrix(A, n);
    char buf[2048];

    snprintf(buf, sizeof buf, "Eigenvalues[%s, Method -> \"Banded\"]", m);
    Expr* eval_list = eval_string(buf);
    ASSERT(list_len_eq(eval_list, n));
    double* lambdas = (double*)malloc(sizeof(double) * n);
    for (size_t i = 0; i < n; i++)
        lambdas[i] = extract_real(eval_list->data.function.args[i]);
    expr_free(eval_list);

    snprintf(buf, sizeof buf, "Eigenvectors[%s, Method -> \"Banded\"]", m);
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
    double tol = 1e-10 * (double)n * normA;
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

/* Orthonormality: ||V^T V - I||_inf < n * eps. */
void test_banded_real_orthonormal(void) {
    double A[25] = {
         3, -1,  0.5, 0, 0,
        -1,  4, -1,  0.5, 0,
         0.5, -1,  5, -1,  0.5,
         0,  0.5, -1,  6, -1,
         0,  0, 0.5, -1, 7
    };
    size_t n = 5;
    char* m = fmt_real_matrix(A, n);
    char buf[2048];
    snprintf(buf, sizeof buf, "Eigenvectors[%s, Method -> \"Banded\"]", m);
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

    /* V[i,*] = i-th eigenvector.  V V^T = I (rows form an orthonormal set). */
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            double s = 0.0;
            for (size_t k = 0; k < n; k++) s += V[i * n + k] * V[j * n + k];
            double expected = (i == j) ? 1.0 : 0.0;
            ASSERT(fabs(s - expected) < 1e-11);
        }
    }
    free(V);
}

/* ============================================================ *
 *  Complex Hermitian Banded                                      *
 * ============================================================ */

/* 3x3 complex Hermitian tridiagonal: real eigenvalues. */
void test_banded_complex_herm_tridiag_3x3(void) {
    Expr* bn = eval_string(
        "Eigenvalues[N[{{1, I, 0}, {-I, 2, I}, {0, -I, 3}}], "
        "Method -> \"Banded\"]");
    Expr* dr = eval_string(
        "Eigenvalues[N[{{1, I, 0}, {-I, 2, I}, {0, -I, 3}}], "
        "Method -> \"Direct\"]");
    ASSERT(list_len_eq(bn, 3));
    ASSERT(list_len_eq(dr, 3));
    for (size_t i = 0; i < 3; i++) {
        double vb = extract_real(bn->data.function.args[i]);
        double vd = extract_real(dr->data.function.args[i]);
        ASSERT(fabs(vb - vd) < 1e-11);
    }
    expr_free(bn); expr_free(dr);
}

/* 5x5 complex Hermitian pentadiagonal -- b = 2. */
void test_banded_complex_herm_pentadiag_5x5(void) {
    Expr* bn = eval_string(
        "Eigenvalues[N[{"
        "{4, I, 0.5, 0, 0},"
        "{-I, 4, I, 0.5, 0},"
        "{0.5, -I, 4, I, 0.5},"
        "{0, 0.5, -I, 4, I},"
        "{0, 0, 0.5, -I, 4}"
        "}], Method -> \"Banded\"]");
    Expr* dr = eval_string(
        "Eigenvalues[N[{"
        "{4, I, 0.5, 0, 0},"
        "{-I, 4, I, 0.5, 0},"
        "{0.5, -I, 4, I, 0.5},"
        "{0, 0.5, -I, 4, I},"
        "{0, 0, 0.5, -I, 4}"
        "}], Method -> \"Direct\"]");
    ASSERT(list_len_eq(bn, 5));
    ASSERT(list_len_eq(dr, 5));
    for (size_t i = 0; i < 5; i++) {
        double vb = extract_real(bn->data.function.args[i]);
        double vd = extract_real(dr->data.function.args[i]);
        ASSERT(fabs(vb - vd) < 1e-11);
    }
    expr_free(bn); expr_free(dr);
}

/* Eigenvectors shape for complex Hermitian banded -- 3 vectors of 3
 * components, each Real or Complex[re, im]. */
void test_banded_complex_eigvec_shape(void) {
    Expr* v = eval_string(
        "Eigenvectors[N[{{1, I, 0}, {-I, 2, I}, {0, -I, 3}}], "
        "Method -> \"Banded\"]");
    ASSERT(list_len_eq(v, 3));
    for (size_t i = 0; i < 3; i++) {
        ASSERT(list_len_eq(v->data.function.args[i], 3));
    }
    expr_free(v);
}

/* ============================================================ *
 *  Automatic dispatch + edge cases                              *
 * ============================================================ */

/* For a 12x12 symmetric tridiagonal (b = 1, n > 8), Automatic prefers
 * Banded.  Result must agree with explicit Banded. */
void test_banded_automatic_picks_for_tridiag(void) {
    /* 12x12 symmetric tridiagonal with random-ish entries. */
    size_t n = 12;
    double A[144] = { 0 };
    for (size_t i = 0; i < n; i++) {
        A[i * n + i] = 4.0 + 0.1 * (double)i;
        if (i + 1 < n) {
            A[i * n + (i + 1)] = 1.0 - 0.05 * (double)i;
            A[(i + 1) * n + i] = 1.0 - 0.05 * (double)i;
        }
    }
    char* m = fmt_real_matrix(A, n);
    char buf[8192];
    snprintf(buf, sizeof buf, "Eigenvalues[%s]", m);
    Expr* auto_res = eval_string(buf);
    snprintf(buf, sizeof buf, "Eigenvalues[%s, Method -> \"Banded\"]", m);
    Expr* bn = eval_string(buf);
    ASSERT(list_len_eq(auto_res, n));
    ASSERT(list_len_eq(bn, n));
    for (size_t i = 0; i < n; i++) {
        double va = extract_real(auto_res->data.function.args[i]);
        double vb = extract_real(bn->data.function.args[i]);
        ASSERT(fabs(va - vb) < 1e-10);
    }
    expr_free(auto_res); expr_free(bn);
    free(m);
}

/* 1x1 edge case. */
void test_banded_1x1(void) {
    Expr* bn = eval_string(
        "Eigenvalues[N[{{7.0}}], Method -> \"Banded\"]");
    ASSERT(list_len_eq(bn, 1));
    double v = extract_real(bn->data.function.args[0]);
    ASSERT(fabs(v - 7.0) < 1e-12);
    expr_free(bn);
}

/* UpTo[k] form. */
void test_banded_upto(void) {
    Expr* bn = eval_string(
        "Eigenvalues[N[{{4, 1, 0}, {1, 3, 1}, {0, 1, 2}}], UpTo[2], "
        "Method -> \"Banded\"]");
    ASSERT(list_len_eq(bn, 2));
    expr_free(bn);
}

#ifdef USE_MPFR

/* ============================================================ *
 *  MPFR Banded tests                                            *
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

/* MPFR real symmetric tridiag: Banded matches Direct at high precision. */
void test_mpfr_banded_real_tridiag_5x5(void) {
    Expr* bn = eval_string(
        "Eigenvalues[N[{{2,-1,0,0,0},{-1,2,-1,0,0},{0,-1,2,-1,0},"
        "{0,0,-1,2,-1},{0,0,0,-1,2}}, 80], Method -> \"Banded\"]");
    Expr* dr = eval_string(
        "Eigenvalues[N[{{2,-1,0,0,0},{-1,2,-1,0,0},{0,-1,2,-1,0},"
        "{0,0,-1,2,-1},{0,0,0,-1,2}}, 80], Method -> \"Direct\"]");
    ASSERT(list_len_eq(bn, 5));
    ASSERT(list_len_eq(dr, 5));
    for (size_t i = 0; i < 5; i++) {
        double mb = extract_abs_any(bn->data.function.args[i]);
        double md = extract_abs_any(dr->data.function.args[i]);
        ASSERT(fabs(mb - md) < 1e-20);
    }
    expr_free(bn); expr_free(dr);
}

/* MPFR real pentadiagonal -- b = 2 with chase. */
void test_mpfr_banded_real_pentadiag_6x6(void) {
    Expr* bn = eval_string(
        "Eigenvalues[N[{{4,1,1/2,0,0,0},{1,4,1,1/2,0,0},{1/2,1,4,1,1/2,0},"
        "{0,1/2,1,4,1,1/2},{0,0,1/2,1,4,1},{0,0,0,1/2,1,4}}, 80], "
        "Method -> \"Banded\"]");
    Expr* dr = eval_string(
        "Eigenvalues[N[{{4,1,1/2,0,0,0},{1,4,1,1/2,0,0},{1/2,1,4,1,1/2,0},"
        "{0,1/2,1,4,1,1/2},{0,0,1/2,1,4,1},{0,0,0,1/2,1,4}}, 80], "
        "Method -> \"Direct\"]");
    ASSERT(list_len_eq(bn, 6));
    ASSERT(list_len_eq(dr, 6));
    for (size_t i = 0; i < 6; i++) {
        double mb = extract_abs_any(bn->data.function.args[i]);
        double md = extract_abs_any(dr->data.function.args[i]);
        ASSERT(fabs(mb - md) < 1e-18);
    }
    expr_free(bn); expr_free(dr);
}

/* MPFR complex Hermitian tridiag. */
void test_mpfr_banded_complex_herm_tridiag_3x3(void) {
    Expr* bn = eval_string(
        "Eigenvalues[N[{{1, I, 0}, {-I, 2, I}, {0, -I, 3}}, 80], "
        "Method -> \"Banded\"]");
    Expr* dr = eval_string(
        "Eigenvalues[N[{{1, I, 0}, {-I, 2, I}, {0, -I, 3}}, 80], "
        "Method -> \"Direct\"]");
    ASSERT(list_len_eq(bn, 3));
    ASSERT(list_len_eq(dr, 3));
    for (size_t i = 0; i < 3; i++) {
        double mb = extract_abs_any(bn->data.function.args[i]);
        double md = extract_abs_any(dr->data.function.args[i]);
        ASSERT(fabs(mb - md) < 1e-20);
    }
    expr_free(bn); expr_free(dr);
}

/* MPFR eigenvector shape (real symmetric). */
void test_mpfr_banded_eigvec_shape(void) {
    Expr* v = eval_string(
        "Eigenvectors[N[{{2,-1,0,0,0},{-1,2,-1,0,0},{0,-1,2,-1,0},"
        "{0,0,-1,2,-1},{0,0,0,-1,2}}, 60], Method -> \"Banded\"]");
    ASSERT(list_len_eq(v, 5));
    for (size_t i = 0; i < 5; i++) {
        ASSERT(list_len_eq(v->data.function.args[i], 5));
    }
    expr_free(v);
}

#endif /* USE_MPFR */

int main(void) {
    symtab_init();
    core_init();

    TEST(test_banded_real_sym_tridiag_5x5);
    TEST(test_banded_real_sym_pentadiag_6x6);
    TEST(test_banded_real_sym_heptadiag_7x7);
    TEST(test_banded_diag_4x4);
    TEST(test_banded_dense_fallback_3x3);
    TEST(test_banded_nonsym_fallback);
    TEST(test_banded_top_k);
    TEST(test_banded_real_residual);
    TEST(test_banded_real_orthonormal);

    TEST(test_banded_complex_herm_tridiag_3x3);
    TEST(test_banded_complex_herm_pentadiag_5x5);
    TEST(test_banded_complex_eigvec_shape);

    TEST(test_banded_automatic_picks_for_tridiag);
    TEST(test_banded_1x1);
    TEST(test_banded_upto);

#ifdef USE_MPFR
    TEST(test_mpfr_banded_real_tridiag_5x5);
    TEST(test_mpfr_banded_real_pentadiag_6x6);
    TEST(test_mpfr_banded_complex_herm_tridiag_3x3);
    TEST(test_mpfr_banded_eigvec_shape);
#endif

    printf("All test_mateigen_banded tests passed!\n");
    return 0;
}
