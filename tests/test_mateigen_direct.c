/* test_mateigen_direct.c — Tests for the numerical "Direct" Eigenvalues
 * / Eigenvectors kernel.
 *
 * Step 2.1 coverage (this file): real symmetric matrices at machine
 * precision (double).  Subsequent commits extend this binary with
 * the non-symmetric, complex, and MPFR cases.
 *
 * Every test runs against the corpus assertions in eigen_corpus.h:
 *   - eigenvalue residual ||A v - lambda v||_inf
 *   - eigenvalue sum = trace(A)
 *   - orthonormality V V^T = I
 *   - eigenvalue sort order (descending |lambda|, stable on ties)
 *   - k-spec interaction (first k, last k, UpTo[k])
 */

#include <assert.h>
#include <math.h>
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

/* eps tolerance for residual / orthonormality checks.  c=64 absorbs
 * the Householder + symmetric-QR constant; tighter is fragile. */
static double machine_tol(const double* A, size_t n) {
    double normA = corpus_norm_inf_real(A, n);
    if (normA == 0.0) normA = 1.0;
    return 64.0 * (double)n * 2.220446049250313e-16 * normA;
}

/* Verify every (lambda_i, v_i) pair for a real symmetric matrix A. */
static void verify_full_spectrum(const double* A, size_t n) {
    double tol = machine_tol(A, n);
    double *lambdas = NULL, *V = NULL;
    size_t kL = corpus_eval_eigenvalues_real(A, n, &lambdas);
    size_t kV = corpus_eval_eigenvectors_real(A, n, &V);
    ASSERT(kL == n && kV == n);

    /* Sum of eigenvalues == trace(A). */
    double trace = 0.0;
    for (size_t i = 0; i < n; i++) trace += A[i * n + i];
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) sum += lambdas[i];
    if (fabs(sum - trace) > tol) {
        printf("FAIL: sum(lambda)=%g != trace=%g (tol %g)\n", sum, trace, tol);
        ASSERT(0);
    }

    /* Per-eigenpair residual. */
    for (size_t i = 0; i < n; i++) {
        corpus_assert_residual_real(A, n, lambdas[i], &V[i * n], tol);
    }

    /* Sort order: descending |lambda|, stable on ties. */
    for (size_t i = 1; i < n; i++) {
        ASSERT(fabs(lambdas[i - 1]) >= fabs(lambdas[i]) - tol);
    }

    /* Orthonormality (Hermitian / real-symmetric spectrum is
     * orthogonal; with our normalisation it's orthonormal). */
    corpus_assert_orthonormal_real(V, n, 64.0 * (double)n * 1e-15);

    free(lambdas); free(V);
}

/* --- 2x2 / 3x3 / 5x5 specific real symmetric matrices ---------------- */

void test_direct_diag_3x3(void) {
    /* Diagonal: eigenvalues = diagonal (sorted by |lambda| descending). */
    double A[9] = {
        3.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 5.0
    };
    double expected[3] = { 5.0, 3.0, 1.0 };
    char* in = corpus_matrix_to_eigenvalues_input(A, 3);
    corpus_check_real_eigenvalues("diag(3,1,5) machine", in, expected, 3, 1e-12);
    free(in);
    verify_full_spectrum(A, 3);
}

void test_direct_identity_4x4(void) {
    double A[16];
    for (size_t i = 0; i < 16; i++) A[i] = 0.0;
    for (size_t i = 0; i < 4; i++) A[i * 4 + i] = 1.0;
    double expected[4] = { 1.0, 1.0, 1.0, 1.0 };
    char* in = corpus_matrix_to_eigenvalues_input(A, 4);
    corpus_check_real_eigenvalues("I_4 machine", in, expected, 4, 1e-12);
    free(in);
    verify_full_spectrum(A, 4);
}

void test_direct_zero_3x3(void) {
    double A[9] = { 0.0 };
    double expected[3] = { 0.0, 0.0, 0.0 };
    char* in = corpus_matrix_to_eigenvalues_input(A, 3);
    corpus_check_real_eigenvalues("zero machine", in, expected, 3, 1e-12);
    free(in);
}

void test_direct_2x2_known(void) {
    /* [[4, 1], [1, 4]] -> eigenvalues 5 and 3. */
    double A[4] = { 4.0, 1.0, 1.0, 4.0 };
    double expected[2] = { 5.0, 3.0 };
    char* in = corpus_matrix_to_eigenvalues_input(A, 2);
    corpus_check_real_eigenvalues("[[4,1],[1,4]] machine", in, expected, 2, 1e-12);
    free(in);
    verify_full_spectrum(A, 2);
}

void test_direct_2x2_negative(void) {
    /* [[1, 2], [2, -1]] -> eigenvalues sqrt(5) and -sqrt(5). */
    double A[4] = { 1.0, 2.0, 2.0, -1.0 };
    double expected[2] = { sqrt(5.0), -sqrt(5.0) };  /* same |lambda|, stable. */
    char* in = corpus_matrix_to_eigenvalues_input(A, 2);
    corpus_check_real_eigenvalues("[[1,2],[2,-1]] machine", in, expected, 2, 1e-12);
    free(in);
    verify_full_spectrum(A, 2);
}

void test_direct_tridiag_3x3(void) {
    /* Toeplitz tridiag [[2,-1,0],[-1,2,-1],[0,-1,2]] -> 2 - sqrt(2), 2, 2 + sqrt(2). */
    double A[9] = {
         2.0, -1.0,  0.0,
        -1.0,  2.0, -1.0,
         0.0, -1.0,  2.0
    };
    double expected[3] = { 2.0 + sqrt(2.0), 2.0, 2.0 - sqrt(2.0) };
    char* in = corpus_matrix_to_eigenvalues_input(A, 3);
    corpus_check_real_eigenvalues("Toeplitz tridiag machine", in, expected, 3, 1e-12);
    free(in);
    verify_full_spectrum(A, 3);
}

void test_direct_tridiag_5x5(void) {
    double A[25] = {
         2.0, -1.0,  0.0,  0.0,  0.0,
        -1.0,  2.0, -1.0,  0.0,  0.0,
         0.0, -1.0,  2.0, -1.0,  0.0,
         0.0,  0.0, -1.0,  2.0, -1.0,
         0.0,  0.0,  0.0, -1.0,  2.0
    };
    /* Known closed form: lambda_k = 2 (1 - cos(k pi / (n+1))), k=1..n. */
    double expected[5];
    for (size_t k = 0; k < 5; k++) {
        double th = (double)(5 - k) * 3.14159265358979323846 / 6.0;
        expected[k] = 2.0 * (1.0 - cos(th));
    }
    char* in = corpus_matrix_to_eigenvalues_input(A, 5);
    corpus_check_real_eigenvalues("5x5 tridiag machine", in, expected, 5, 1e-11);
    free(in);
    verify_full_spectrum(A, 5);
}

void test_direct_hilbert_4(void) {
    /* Hilbert matrix: ill-conditioned but eigenvalues are well-defined. */
    double A[16];
    for (size_t i = 0; i < 4; i++)
        for (size_t j = 0; j < 4; j++)
            A[i * 4 + j] = 1.0 / (double)(i + j + 1);
    /* Reference eigenvalues (LAPACK / Mathematica) for n=4:
     *   1.5002142800592426, 0.1691412202214158,
     *   0.006738273605760491, 9.670230402035702e-5 */
    double expected[4] = {
        1.5002142800592426,
        0.1691412202214158,
        0.006738273605760491,
        9.670230402035702e-5
    };
    char* in = corpus_matrix_to_eigenvalues_input(A, 4);
    corpus_check_real_eigenvalues("Hilbert_4 machine", in, expected, 4, 1e-10);
    free(in);
    verify_full_spectrum(A, 4);
}

void test_direct_random_symmetric_5x5(void) {
    /* Deterministic "random" via a Lehmer LCG seeded at 12345.  Then
     * symmetrise.  No external randomness so failures reproduce. */
    unsigned long s = 12345UL;
    double A[25];
    for (size_t i = 0; i < 25; i++) {
        s = (s * 48271UL) & 0x7fffffffUL;
        A[i] = ((double)s / 2147483647.0 - 0.5) * 10.0;
    }
    /* Symmetrise: A <- (A + A^T) / 2. */
    double S[25];
    for (size_t i = 0; i < 5; i++)
        for (size_t j = 0; j < 5; j++)
            S[i * 5 + j] = 0.5 * (A[i * 5 + j] + A[j * 5 + i]);
    verify_full_spectrum(S, 5);
}

void test_direct_random_symmetric_10x10(void) {
    unsigned long s = 67890UL;
    double A[100];
    for (size_t i = 0; i < 100; i++) {
        s = (s * 48271UL) & 0x7fffffffUL;
        A[i] = ((double)s / 2147483647.0 - 0.5) * 10.0;
    }
    double S[100];
    for (size_t i = 0; i < 10; i++)
        for (size_t j = 0; j < 10; j++)
            S[i * 10 + j] = 0.5 * (A[i * 10 + j] + A[j * 10 + i]);
    verify_full_spectrum(S, 10);
}

void test_direct_repeated_eigenvalues(void) {
    /* Matrix with eigenvalues 2, 2, 5: 2*I + outer(v,v) with v = (1,1,1).
     * Construction:  M = 2 I + outer(v, v),  outer(v,v) has rank 1,
     * eigenvalues |v|^2 = 3 (once) and 0 (twice).  Adding 2 I yields
     * eigenvalues 5, 2, 2. */
    double A[9] = {
        3.0, 1.0, 1.0,
        1.0, 3.0, 1.0,
        1.0, 1.0, 3.0
    };
    double expected[3] = { 5.0, 2.0, 2.0 };
    char* in = corpus_matrix_to_eigenvalues_input(A, 3);
    corpus_check_real_eigenvalues("repeated eigenvalues machine",
                                   in, expected, 3, 1e-12);
    free(in);
    verify_full_spectrum(A, 3);
}

/* --- k-spec interaction --------------------------------------------- */

void test_direct_k_positive(void) {
    /* Toeplitz tridiag 5x5, k=2: top-2 |lambda|. */
    double A[25] = {
         2.0, -1.0,  0.0,  0.0,  0.0,
        -1.0,  2.0, -1.0,  0.0,  0.0,
         0.0, -1.0,  2.0, -1.0,  0.0,
         0.0,  0.0, -1.0,  2.0, -1.0,
         0.0,  0.0,  0.0, -1.0,  2.0
    };
    double expected[2];
    for (size_t k = 0; k < 2; k++) {
        double th = (double)(5 - k) * 3.14159265358979323846 / 6.0;
        expected[k] = 2.0 * (1.0 - cos(th));
    }
    char* in0 = corpus_matrix_to_eigenvalues_input(A, 5);
    /* Patch input: "Eigenvalues[<mat>]" -> "Eigenvalues[<mat>, 2]". */
    size_t blen = strlen(in0);
    char* in = (char*)malloc(blen + 8);
    snprintf(in, blen + 8, "%.*s, 2]", (int)(blen - 1), in0);
    corpus_check_real_eigenvalues("k=2 machine", in, expected, 2, 1e-11);
    free(in); free(in0);
}

void test_direct_k_negative(void) {
    double A[25] = {
         2.0, -1.0,  0.0,  0.0,  0.0,
        -1.0,  2.0, -1.0,  0.0,  0.0,
         0.0, -1.0,  2.0, -1.0,  0.0,
         0.0,  0.0, -1.0,  2.0, -1.0,
         0.0,  0.0,  0.0, -1.0,  2.0
    };
    /* Smallest 2 |lambda|: the last two of the n=5 spectrum. */
    double expected[2];
    for (size_t k = 0; k < 2; k++) {
        size_t idx = 3 + k;
        double th = (double)(5 - idx) * 3.14159265358979323846 / 6.0;
        expected[k] = 2.0 * (1.0 - cos(th));
    }
    char* in0 = corpus_matrix_to_eigenvalues_input(A, 5);
    size_t blen = strlen(in0);
    char* in = (char*)malloc(blen + 8);
    snprintf(in, blen + 8, "%.*s, -2]", (int)(blen - 1), in0);
    corpus_check_real_eigenvalues("k=-2 machine", in, expected, 2, 1e-11);
    free(in); free(in0);
}

void test_direct_upto_oversized(void) {
    /* UpTo[100] on a 3x3 must clamp to 3. */
    double A[9] = {
        3.0, 0.0, 0.0,
        0.0, 1.0, 0.0,
        0.0, 0.0, 5.0
    };
    double expected[3] = { 5.0, 3.0, 1.0 };
    corpus_check_real_eigenvalues(
        "UpTo[100] clamps machine",
        "Eigenvalues[{{3.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 5.0}}, UpTo[100]]",
        expected, 3, 1e-12);
    (void)A;
}

/* --- Method dispatch -------------------------------------------------- */

void test_direct_method_direct_equiv_automatic(void) {
    /* Method -> "Direct" must produce the same result as Method -> Automatic. */
    double A[9] = {
         2.0, -1.0,  0.0,
        -1.0,  2.0, -1.0,
         0.0, -1.0,  2.0
    };
    double *L1 = NULL, *L2 = NULL;
    size_t k1 = corpus_eval_eigenvalues_real(A, 3, &L1);
    ASSERT(k1 == 3);

    /* Method -> "Direct" path. */
    Expr* e = parse_expression(
        "Eigenvalues[{{2.0, -1.0, 0.0}, {-1.0, 2.0, -1.0}, "
        "{0.0, -1.0, 2.0}}, Method -> \"Direct\"]");
    ASSERT(e);
    Expr* r = evaluate(e);
    ASSERT(r);
    ASSERT(r->type == EXPR_FUNCTION);
    ASSERT(r->data.function.arg_count == 3);
    L2 = (double*)malloc(sizeof(double) * 3);
    for (size_t i = 0; i < 3; i++) {
        Expr* v = r->data.function.args[i];
        ASSERT(v->type == EXPR_REAL);
        L2[i] = v->data.real;
    }
    expr_free(r); expr_free(e);
    for (size_t i = 0; i < 3; i++) {
        if (fabs(L1[i] - L2[i]) > 1e-12) {
            printf("FAIL: Direct vs Automatic differ at i=%zu (%g vs %g)\n",
                   i, L1[i], L2[i]);
            ASSERT(0);
        }
    }
    printf("PASS: Direct == Automatic for symmetric 3x3 machine\n");
    free(L1); free(L2);
}

void test_direct_non_symmetric_falls_back(void) {
    /* Non-symmetric numeric matrix: step 2.1 dispatcher returns NULL
     * and the call falls back to the symbolic char-poly + Solve path.
     * Result should still be a List of three numeric values for this
     * upper-triangular case (eigenvalues are the diagonal). */
    Expr* e = parse_expression(
        "Eigenvalues[{{1.0, 2.0, 3.0}, {0.0, 4.0, 5.0}, {0.0, 0.0, 6.0}}]");
    ASSERT(e);
    Expr* r = evaluate(e);
    ASSERT(r);
    ASSERT(r->type == EXPR_FUNCTION);
    ASSERT(r->data.function.arg_count == 3);
    /* Sorted descending |lambda|: 6, 4, 1.  Allow integer or real leaves. */
    double got[3];
    for (size_t i = 0; i < 3; i++) {
        Expr* v = r->data.function.args[i];
        if (v->type == EXPR_REAL) got[i] = v->data.real;
        else if (v->type == EXPR_INTEGER) got[i] = (double)v->data.integer;
        else { ASSERT(0); }
    }
    ASSERT(fabs(got[0] - 6.0) < 1e-10);
    ASSERT(fabs(got[1] - 4.0) < 1e-10);
    ASSERT(fabs(got[2] - 1.0) < 1e-10);
    printf("PASS: non-symmetric matrix falls back to symbolic path\n");
    expr_free(r); expr_free(e);
}

void test_direct_symbolic_ignores_dispatch(void) {
    /* Pure-symbolic matrices must continue to use the symbolic
     * pipeline.  The dispatcher checks `inexact` so an Integer matrix
     * is never routed through the numerical kernel. */
    Expr* e = parse_expression(
        "Eigenvalues[{{2, -1, 0}, {-1, 2, -1}, {0, -1, 2}}]");
    ASSERT(e);
    Expr* r = evaluate(e);
    ASSERT(r);
    ASSERT(r->type == EXPR_FUNCTION);
    /* Symbolic result: at least one entry must not be EXPR_REAL --
     * either Integer 2, or Plus[2, Sqrt[2]], etc. */
    int saw_non_real = 0;
    for (size_t i = 0; i < r->data.function.arg_count; i++) {
        if (r->data.function.args[i]->type != EXPR_REAL) saw_non_real = 1;
    }
    ASSERT(saw_non_real);
    printf("PASS: integer matrix stays on the symbolic path\n");
    expr_free(r); expr_free(e);
}

/* --- Eigenvectors via the new dispatcher ---------------------------- */

void test_direct_eigenvectors_diag_4x4(void) {
    /* Diagonal: eigenvectors are the canonical basis (descending |lambda|
     * order: lambda=5 first => e2, then 3 => e0, then 2 => e3, then 1 => e1). */
    double A[16] = {
        3.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 5.0, 0.0,
        0.0, 0.0, 0.0, 2.0
    };
    double *V = NULL;
    size_t kV = corpus_eval_eigenvectors_real(A, 4, &V);
    ASSERT(kV == 4);
    /* First row corresponds to lambda=5 (max |.|): should be e2 = (0,0,1,0)
     * up to a sign.  Likewise:
     *   lambda=3 -> e0 = (1,0,0,0)
     *   lambda=2 -> e3 = (0,0,0,1)
     *   lambda=1 -> e1 = (0,1,0,0)
     */
    /* Sign-tolerant comparison via absolute values. */
    double expect[16] = {
        0,0,1,0,
        1,0,0,0,
        0,0,0,1,
        0,1,0,0
    };
    for (size_t i = 0; i < 16; i++) {
        ASSERT(fabs(fabs(V[i]) - expect[i]) < 1e-12);
    }
    /* Orthonormality. */
    corpus_assert_orthonormal_real(V, 4, 1e-13);
    printf("PASS: eigenvectors of diag(3,1,5,2) match canonical basis (sign-free)\n");
    free(V);
}

void test_direct_eigenvectors_tridiag_residual(void) {
    /* Larger tridiagonal: verify the residual closure A v_i = lambda_i v_i
     * for every (lambda, v) pair. */
    double A[36] = {
         2.0, -1.0,  0.0,  0.0,  0.0,  0.0,
        -1.0,  2.0, -1.0,  0.0,  0.0,  0.0,
         0.0, -1.0,  2.0, -1.0,  0.0,  0.0,
         0.0,  0.0, -1.0,  2.0, -1.0,  0.0,
         0.0,  0.0,  0.0, -1.0,  2.0, -1.0,
         0.0,  0.0,  0.0,  0.0, -1.0,  2.0
    };
    verify_full_spectrum(A, 6);
    printf("PASS: 6x6 tridiagonal residual + orthonormality\n");
}

int main(void) {
    symtab_init();
    core_init();
    printf("Running mateigen Direct (Phase 2.1: real symmetric machine) tests...\n");

    TEST(test_direct_diag_3x3);
    TEST(test_direct_identity_4x4);
    TEST(test_direct_zero_3x3);
    TEST(test_direct_2x2_known);
    TEST(test_direct_2x2_negative);
    TEST(test_direct_tridiag_3x3);
    TEST(test_direct_tridiag_5x5);
    TEST(test_direct_hilbert_4);
    TEST(test_direct_random_symmetric_5x5);
    TEST(test_direct_random_symmetric_10x10);
    TEST(test_direct_repeated_eigenvalues);
    TEST(test_direct_k_positive);
    TEST(test_direct_k_negative);
    TEST(test_direct_upto_oversized);
    TEST(test_direct_method_direct_equiv_automatic);
    TEST(test_direct_non_symmetric_falls_back);
    TEST(test_direct_symbolic_ignores_dispatch);
    TEST(test_direct_eigenvectors_diag_4x4);
    TEST(test_direct_eigenvectors_tridiag_residual);

    printf("All mateigen Direct (step 2.1) tests passed!\n");
    return 0;
}
