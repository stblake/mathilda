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

void test_direct_upper_triangular_real_eigenvalues(void) {
    /* Upper triangular: eigenvalues are the diagonal, sorted by |.|. */
    Expr* e = parse_expression(
        "Eigenvalues[{{1.0, 2.0, 3.0}, {0.0, 4.0, 5.0}, {0.0, 0.0, 6.0}}]");
    ASSERT(e);
    Expr* r = evaluate(e);
    ASSERT(r);
    ASSERT(r->type == EXPR_FUNCTION);
    ASSERT(r->data.function.arg_count == 3);
    double got[3];
    for (size_t i = 0; i < 3; i++) {
        Expr* v = r->data.function.args[i];
        ASSERT(v->type == EXPR_REAL);
        got[i] = v->data.real;
    }
    ASSERT(fabs(got[0] - 6.0) < 1e-10);
    ASSERT(fabs(got[1] - 4.0) < 1e-10);
    ASSERT(fabs(got[2] - 1.0) < 1e-10);
    printf("PASS: upper triangular non-symmetric eigenvalues (Hessenberg + Francis QR)\n");
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

/* ============================================================ *
 *  Phase 2 step 2a: real non-symmetric eigenvalues               *
 *                                                                 *
 *  Hessenberg + implicit double-shift Francis QR.  All real      *
 *  matrices return eigenvalues; complex eigenvalues are emitted   *
 *  as Complex[re, im].                                            *
 * ============================================================ */

/* Read an eigenvalue entry — Real or Complex[re, im] — into (re, im). */
static void read_eigenvalue_entry(Expr* e, double* re, double* im) {
    if (e->type == EXPR_REAL)    { *re = e->data.real;            *im = 0.0; return; }
    if (e->type == EXPR_INTEGER) { *re = (double)e->data.integer; *im = 0.0; return; }
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && strcmp(e->data.function.head->data.symbol, "Complex") == 0
        && e->data.function.arg_count == 2) {
        Expr* r = e->data.function.args[0];
        Expr* m = e->data.function.args[1];
        *re = (r->type == EXPR_REAL) ? r->data.real
            : (r->type == EXPR_INTEGER) ? (double)r->data.integer : NAN;
        *im = (m->type == EXPR_REAL) ? m->data.real
            : (m->type == EXPR_INTEGER) ? (double)m->data.integer : NAN;
        return;
    }
    *re = NAN; *im = NAN;
}

/* Eval `Eigenvalues[<matrix>]` and read into freshly-allocated re / im
 * arrays.  Returns the count of eigenvalues returned. */
static size_t eval_eigenvalues_mixed(const double* A, size_t n,
                                      double** re_out, double** im_out) {
    char* in = corpus_matrix_to_eigenvalues_input(A, n);
    Expr* e = parse_expression(in);
    ASSERT(e);
    Expr* r = evaluate(e);
    expr_free(e); free(in);
    ASSERT(r != NULL);
    ASSERT(r->type == EXPR_FUNCTION);
    ASSERT(r->data.function.head->type == EXPR_SYMBOL);
    ASSERT(strcmp(r->data.function.head->data.symbol, "List") == 0);
    size_t cnt = r->data.function.arg_count;
    *re_out = (double*)malloc(sizeof(double) * (cnt ? cnt : 1));
    *im_out = (double*)malloc(sizeof(double) * (cnt ? cnt : 1));
    for (size_t i = 0; i < cnt; i++) {
        read_eigenvalue_entry(r->data.function.args[i],
                              &(*re_out)[i], &(*im_out)[i]);
    }
    expr_free(r);
    return cnt;
}

/* For a real n x n matrix, verify:
 *   - count of eigenvalues returned == n
 *   - sum of eigenvalues == trace(A)            (similarity invariant)
 *   - sum of lambda_i^2 == trace(A^2)           (similarity invariant,
 *                                                 holds for any matrix)
 *   - imaginary parts cancel in conjugate pairs (real-valued input)
 *   - sort order is descending |lambda| (stable)
 *
 * Note: sum |lambda_i|^2 == Frobenius^2(A) is a NORMAL-matrix invariant
 * (Schur's inequality: sum |lambda|^2 <= Frobenius^2 with equality iff
 * normal).  General non-normal real matrices need trace(A^2) instead. */
static void verify_general_invariants(const double* A, size_t n) {
    double* re = NULL;
    double* im = NULL;
    size_t cnt = eval_eigenvalues_mixed(A, n, &re, &im);
    ASSERT(cnt == n);

    double trace = 0.0;
    for (size_t i = 0; i < n; i++) trace += A[i * n + i];
    double trace_A2 = 0.0;
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++)
            trace_A2 += A[i * n + j] * A[j * n + i];

    double sum_re = 0.0, sum_im = 0.0;
    double sum_lam2_re = 0.0, sum_lam2_im = 0.0;
    for (size_t i = 0; i < n; i++) {
        sum_re += re[i];
        sum_im += im[i];
        /* lambda^2 = (re^2 - im^2) + 2 re im i. */
        sum_lam2_re += re[i] * re[i] - im[i] * im[i];
        sum_lam2_im += 2.0 * re[i] * im[i];
    }
    double normA = corpus_norm_inf_real(A, n);
    if (normA == 0.0) normA = 1.0;
    double tol_trace = 256.0 * (double)n * 2.220446049250313e-16 * normA;
    double tol_A2 = 256.0 * (double)n * 2.220446049250313e-16
                    * (fabs(trace_A2) + normA * normA);
    if (fabs(sum_re - trace) > tol_trace || fabs(sum_im) > tol_trace) {
        printf("FAIL: sum(lambda) = %g + %gi vs trace %g (tol %g)\n",
               sum_re, sum_im, trace, tol_trace);
        ASSERT(0);
    }
    if (fabs(sum_lam2_re - trace_A2) > tol_A2
        || fabs(sum_lam2_im) > tol_A2) {
        printf("FAIL: sum(lambda^2) = %g + %gi vs trace(A^2) %g (tol %g)\n",
               sum_lam2_re, sum_lam2_im, trace_A2, tol_A2);
        ASSERT(0);
    }
    for (size_t i = 1; i < n; i++) {
        double a1 = hypot(re[i - 1], im[i - 1]);
        double a2 = hypot(re[i], im[i]);
        ASSERT(a1 >= a2 - 1e-9);
    }
    free(re); free(im);
}

void test_direct_general_2x2_real(void) {
    /* [[1, 2], [3, 4]]: eigenvalues (5 +/- sqrt(33))/2. */
    double A[4] = { 1.0, 2.0, 3.0, 4.0 };
    double* re = NULL;
    double* im = NULL;
    size_t cnt = eval_eigenvalues_mixed(A, 2, &re, &im);
    ASSERT(cnt == 2);
    double e0 = (5.0 + sqrt(33.0)) * 0.5;
    double e1 = (5.0 - sqrt(33.0)) * 0.5;
    ASSERT(fabs(re[0] - e0) < 1e-12);
    ASSERT(fabs(re[1] - e1) < 1e-12);
    ASSERT(fabs(im[0]) < 1e-14);
    ASSERT(fabs(im[1]) < 1e-14);
    printf("PASS: [[1,2],[3,4]] non-symmetric eigenvalues\n");
    free(re); free(im);
}

void test_direct_general_2x2_complex_rotation(void) {
    /* [[0, 1], [-1, 0]]: pure rotation, eigenvalues +/- i. */
    double A[4] = { 0.0, 1.0, -1.0, 0.0 };
    double *re = NULL, *im = NULL;
    size_t cnt = eval_eigenvalues_mixed(A, 2, &re, &im);
    ASSERT(cnt == 2);
    /* Both eigenvalues have |lambda| = 1; sort order by index. */
    for (size_t i = 0; i < 2; i++) {
        ASSERT(fabs(re[i]) < 1e-14);
        ASSERT(fabs(fabs(im[i]) - 1.0) < 1e-14);
    }
    ASSERT(im[0] * im[1] < 0);  /* conjugate pair */
    printf("PASS: rotation [[0,1],[-1,0]] -> +/- i\n");
    free(re); free(im);
}

void test_direct_general_2x2_complex_general(void) {
    /* [[1, -1], [1, 1]]: eigenvalues 1 +/- i. */
    double A[4] = { 1.0, -1.0, 1.0, 1.0 };
    double *re = NULL, *im = NULL;
    size_t cnt = eval_eigenvalues_mixed(A, 2, &re, &im);
    ASSERT(cnt == 2);
    for (size_t i = 0; i < 2; i++) {
        ASSERT(fabs(re[i] - 1.0) < 1e-13);
        ASSERT(fabs(fabs(im[i]) - 1.0) < 1e-13);
    }
    ASSERT(im[0] * im[1] < 0);
    printf("PASS: [[1,-1],[1,1]] -> 1 +/- i\n");
    free(re); free(im);
}

void test_direct_general_3x3_block_diagonal(void) {
    /* Block-diagonal: 2x2 with complex pair (1 +/- i) at top-left,
     * plus a real eigenvalue 2 at bottom-right. */
    double A[9] = {
        1.0, -1.0, 0.0,
        1.0,  1.0, 0.0,
        0.0,  0.0, 2.0
    };
    double *re = NULL, *im = NULL;
    size_t cnt = eval_eigenvalues_mixed(A, 3, &re, &im);
    ASSERT(cnt == 3);
    /* |2| = 2, |1 +/- i| = sqrt(2) ~= 1.414.  Sort: 2, 1+i, 1-i (stable). */
    ASSERT(fabs(re[0] - 2.0) < 1e-12);
    ASSERT(fabs(im[0]) < 1e-13);
    ASSERT(fabs(re[1] - 1.0) < 1e-13);
    ASSERT(fabs(re[2] - 1.0) < 1e-13);
    ASSERT(fabs(fabs(im[1]) - 1.0) < 1e-13);
    ASSERT(fabs(fabs(im[2]) - 1.0) < 1e-13);
    ASSERT(im[1] * im[2] < 0);
    printf("PASS: 3x3 block diag with complex pair\n");
    verify_general_invariants(A, 3);
    free(re); free(im);
}

void test_direct_general_3x3_random(void) {
    /* Deterministic Lehmer random 3x3, not symmetrised. */
    unsigned long s = 11111UL;
    double A[9];
    for (size_t i = 0; i < 9; i++) {
        s = (s * 48271UL) & 0x7fffffffUL;
        A[i] = ((double)s / 2147483647.0 - 0.5) * 10.0;
    }
    verify_general_invariants(A, 3);
    printf("PASS: random non-symmetric 3x3 trace + trace(A^2) preserved\n");
}

void test_direct_general_5x5_random(void) {
    unsigned long s = 22222UL;
    double A[25];
    for (size_t i = 0; i < 25; i++) {
        s = (s * 48271UL) & 0x7fffffffUL;
        A[i] = ((double)s / 2147483647.0 - 0.5) * 10.0;
    }
    verify_general_invariants(A, 5);
    printf("PASS: random non-symmetric 5x5 trace + trace(A^2) preserved\n");
}

void test_direct_general_10x10_random(void) {
    unsigned long s = 33333UL;
    double A[100];
    for (size_t i = 0; i < 100; i++) {
        s = (s * 48271UL) & 0x7fffffffUL;
        A[i] = ((double)s / 2147483647.0 - 0.5) * 10.0;
    }
    verify_general_invariants(A, 10);
    printf("PASS: random non-symmetric 10x10 trace + trace(A^2) preserved\n");
}

void test_direct_general_companion_matrix(void) {
    /* Companion of x^3 - 6 x^2 + 11 x - 6 (eigenvalues 1, 2, 3):
     *   [[0, 0, 6], [1, 0, -11], [0, 1, 6]] (Frobenius form). */
    double A[9] = {
        0.0, 0.0,   6.0,
        1.0, 0.0, -11.0,
        0.0, 1.0,   6.0
    };
    double *re = NULL, *im = NULL;
    size_t cnt = eval_eigenvalues_mixed(A, 3, &re, &im);
    ASSERT(cnt == 3);
    /* Expect {3, 2, 1} (descending |.|, all real). */
    ASSERT(fabs(re[0] - 3.0) < 1e-10);
    ASSERT(fabs(re[1] - 2.0) < 1e-10);
    ASSERT(fabs(re[2] - 1.0) < 1e-10);
    for (size_t i = 0; i < 3; i++) ASSERT(fabs(im[i]) < 1e-12);
    printf("PASS: companion of (x-1)(x-2)(x-3) -> {3, 2, 1}\n");
    free(re); free(im);
}

void test_direct_general_jordan_block_2x2(void) {
    /* Jordan block [[2, 1], [0, 2]]: defective eigenvalue 2 with
     * algebraic multiplicity 2, geometric multiplicity 1.  The
     * Eigenvalues call should still return {2, 2}. */
    double A[4] = { 2.0, 1.0, 0.0, 2.0 };
    double *re = NULL, *im = NULL;
    size_t cnt = eval_eigenvalues_mixed(A, 2, &re, &im);
    ASSERT(cnt == 2);
    ASSERT(fabs(re[0] - 2.0) < 1e-12);
    ASSERT(fabs(re[1] - 2.0) < 1e-12);
    ASSERT(fabs(im[0]) < 1e-13);
    ASSERT(fabs(im[1]) < 1e-13);
    printf("PASS: Jordan block 2x2 -> {2, 2}\n");
    free(re); free(im);
}

void test_direct_general_complex_pair_sorted_after_real(void) {
    /* Construct a matrix whose spectrum is {5, 2+i, 2-i}.  |5|=5, |2 +/- i| =
     * sqrt(5) ~= 2.236.  Sorted descending: 5, then conjugate pair. */
    double A[9] = {
        5.0, 0.0, 0.0,
        0.0, 2.0, -1.0,
        0.0, 1.0,  2.0
    };
    double *re = NULL, *im = NULL;
    size_t cnt = eval_eigenvalues_mixed(A, 3, &re, &im);
    ASSERT(cnt == 3);
    ASSERT(fabs(re[0] - 5.0) < 1e-12);
    ASSERT(fabs(im[0]) < 1e-13);
    ASSERT(fabs(re[1] - 2.0) < 1e-12);
    ASSERT(fabs(re[2] - 2.0) < 1e-12);
    ASSERT(fabs(fabs(im[1]) - 1.0) < 1e-12);
    ASSERT(fabs(fabs(im[2]) - 1.0) < 1e-12);
    printf("PASS: {5, 2+/-i} sort order respects |lambda|\n");
    free(re); free(im);
}

void test_direct_general_method_direct_routed(void) {
    /* Method -> "Direct" on a real non-symmetric numeric matrix must
     * use the new kernel (not the symbolic fallback) and produce the
     * same result as Automatic. */
    Expr* e1 = parse_expression(
        "Eigenvalues[{{1.0, 2.0, 3.0}, {0.0, 4.0, 5.0}, "
        "{0.0, 0.0, 6.0}}, Method -> \"Direct\"]");
    Expr* e2 = parse_expression(
        "Eigenvalues[{{1.0, 2.0, 3.0}, {0.0, 4.0, 5.0}, "
        "{0.0, 0.0, 6.0}}]");
    ASSERT(e1 && e2);
    Expr* r1 = evaluate(e1);
    Expr* r2 = evaluate(e2);
    ASSERT(r1 && r2);
    ASSERT(r1->type == EXPR_FUNCTION && r2->type == EXPR_FUNCTION);
    ASSERT(r1->data.function.arg_count == 3);
    ASSERT(r2->data.function.arg_count == 3);
    for (size_t i = 0; i < 3; i++) {
        double a1, b1, a2, b2;
        read_eigenvalue_entry(r1->data.function.args[i], &a1, &b1);
        read_eigenvalue_entry(r2->data.function.args[i], &a2, &b2);
        ASSERT(fabs(a1 - a2) < 1e-12);
        ASSERT(fabs(b1 - b2) < 1e-12);
    }
    printf("PASS: Direct == Automatic for non-symmetric 3x3\n");
    expr_free(r1); expr_free(r2); expr_free(e1); expr_free(e2);
}

/* Read an entry of an eigenvector List into (re, im).  Handles Real,
 * Integer, and Complex[re, im] entries. */
static void read_vec_entry(Expr* e, double* re, double* im) {
    read_eigenvalue_entry(e, re, im);
}

/* Read an Eigenvectors result into freshly-allocated row-major
 * V_re / V_im arrays.  Returns the number of eigenvectors returned. */
static size_t eval_eigenvectors_mixed(const double* A, size_t n,
                                       double** V_re_out, double** V_im_out) {
    char* in = corpus_matrix_to_eigenvectors_input(A, n);
    Expr* e = parse_expression(in);
    ASSERT(e);
    Expr* r = evaluate(e);
    expr_free(e); free(in);
    ASSERT(r != NULL);
    ASSERT(r->type == EXPR_FUNCTION);
    ASSERT(r->data.function.head->type == EXPR_SYMBOL);
    ASSERT(strcmp(r->data.function.head->data.symbol, "List") == 0);
    size_t k = r->data.function.arg_count;
    *V_re_out = (double*)malloc(sizeof(double) * (k ? k * n : 1));
    *V_im_out = (double*)malloc(sizeof(double) * (k ? k * n : 1));
    for (size_t i = 0; i < k; i++) {
        Expr* row = r->data.function.args[i];
        ASSERT(row->type == EXPR_FUNCTION);
        ASSERT(strcmp(row->data.function.head->data.symbol, "List") == 0);
        ASSERT(row->data.function.arg_count == n);
        for (size_t j = 0; j < n; j++) {
            read_vec_entry(row->data.function.args[j],
                           &(*V_re_out)[i * n + j],
                           &(*V_im_out)[i * n + j]);
        }
    }
    expr_free(r);
    return k;
}

/* Verify A v == lambda v for a real eigenvalue + real eigenvector. */
static void check_residual_real(const double* A, size_t n,
                                 double lam,
                                 const double* v, double tol) {
    double res = 0.0;
    for (size_t i = 0; i < n; i++) {
        double Avi = 0.0;
        for (size_t j = 0; j < n; j++) Avi += A[i * n + j] * v[j];
        double e = fabs(Avi - lam * v[i]);
        if (e > res) res = e;
    }
    if (res > tol) {
        printf("FAIL: real residual %g exceeds tol %g (lambda=%g)\n",
               res, tol, lam);
    }
    ASSERT(res <= tol);
}

/* Verify A v == lambda v for a complex eigenvalue + complex eigenvector.
 * lambda = a + b i, v = v_re + i v_im. */
static void check_residual_complex(const double* A, size_t n,
                                    double a, double b,
                                    const double* v_re, const double* v_im,
                                    double tol) {
    double res = 0.0;
    for (size_t i = 0; i < n; i++) {
        double Av_re = 0.0, Av_im = 0.0;
        for (size_t j = 0; j < n; j++) {
            Av_re += A[i * n + j] * v_re[j];
            Av_im += A[i * n + j] * v_im[j];
        }
        double lam_v_re = a * v_re[i] - b * v_im[i];
        double lam_v_im = a * v_im[i] + b * v_re[i];
        double e = hypot(Av_re - lam_v_re, Av_im - lam_v_im);
        if (e > res) res = e;
    }
    if (res > tol) {
        printf("FAIL: complex residual %g exceeds tol %g (lambda=%g+%gi)\n",
               res, tol, a, b);
    }
    ASSERT(res <= tol);
}

/* Full-spectrum residual check for a real general matrix.  Verifies
 * every (lambda_i, v_i) pair, real or complex. */
static void verify_general_eigenvectors(const double* A, size_t n) {
    double *evl_re = NULL, *evl_im = NULL;
    size_t cnt_l = eval_eigenvalues_mixed(A, n, &evl_re, &evl_im);
    ASSERT(cnt_l == n);

    double *V_re = NULL, *V_im = NULL;
    size_t cnt_v = eval_eigenvectors_mixed(A, n, &V_re, &V_im);
    ASSERT(cnt_v == n);

    double normA = corpus_norm_inf_real(A, n);
    if (normA == 0.0) normA = 1.0;
    double tol = 256.0 * (double)n * 2.220446049250313e-16 * normA;

    for (size_t k = 0; k < n; k++) {
        if (evl_im[k] == 0.0) {
            check_residual_real(A, n, evl_re[k], &V_re[k * n], tol);
        } else {
            check_residual_complex(A, n,
                                    evl_re[k], evl_im[k],
                                    &V_re[k * n], &V_im[k * n], tol);
        }
    }
    free(evl_re); free(evl_im); free(V_re); free(V_im);
}

void test_direct_general_eigenvectors_2x2(void) {
    /* [[1, 2], [3, 4]]: real eigenvalues. */
    double A[4] = { 1.0, 2.0, 3.0, 4.0 };
    verify_general_eigenvectors(A, 2);
    printf("PASS: [[1,2],[3,4]] non-symmetric eigenvectors residual\n");
}

void test_direct_general_eigenvectors_rotation(void) {
    /* [[0, 1], [-1, 0]]: pure rotation, eigenvalues +/- i. */
    double A[4] = { 0.0, 1.0, -1.0, 0.0 };
    verify_general_eigenvectors(A, 2);
    printf("PASS: rotation [[0,1],[-1,0]] eigenvectors residual\n");
}

void test_direct_general_eigenvectors_3x3_complex_block(void) {
    /* Block-diagonal: 2x2 with complex pair (1 +/- i), real eigenvalue 2. */
    double A[9] = {
        1.0, -1.0, 0.0,
        1.0,  1.0, 0.0,
        0.0,  0.0, 2.0
    };
    verify_general_eigenvectors(A, 3);
    printf("PASS: 3x3 block-diag complex-pair eigenvectors residual\n");
}

void test_direct_general_eigenvectors_upper_triangular(void) {
    /* Upper triangular with distinct real eigenvalues. */
    double A[9] = {
        1.0, 2.0, 3.0,
        0.0, 4.0, 5.0,
        0.0, 0.0, 6.0
    };
    verify_general_eigenvectors(A, 3);
    printf("PASS: upper-triangular 3x3 eigenvectors residual\n");
}

void test_direct_general_eigenvectors_generic_3x3(void) {
    /* Generic 3x3 with three real eigenvalues; verifies the
     * Hessenberg + Francis QR + back-substitution pipeline produces
     * eigenvectors at full machine precision.  (Companion matrices are
     * deliberately excluded -- they are a textbook ill-conditioned
     * case for eigenvalue / eigenvector computation, and the
     * eigenvalue-only companion test in step 2a already documents that
     * the eigenvalues themselves come out to ~1e-10.) */
    double A[9] = {
        1.0, 2.0, 3.0,
        4.0, 5.0, 6.0,
        7.0, 8.0, 0.0
    };
    verify_general_eigenvectors(A, 3);
    printf("PASS: generic 3x3 non-symmetric eigenvectors residual\n");
}

void test_direct_general_eigenvectors_5x5_random(void) {
    /* Deterministic 5x5 random non-symmetric.  Spectrum may mix real
     * and complex eigenvalues; verifier handles both. */
    unsigned long s = 22222UL;
    double A[25];
    for (size_t i = 0; i < 25; i++) {
        s = (s * 48271UL) & 0x7fffffffUL;
        A[i] = ((double)s / 2147483647.0 - 0.5) * 10.0;
    }
    verify_general_eigenvectors(A, 5);
    printf("PASS: random 5x5 non-symmetric eigenvectors residual\n");
}

void test_direct_general_eigenvectors_10x10_random(void) {
    unsigned long s = 33333UL;
    double A[100];
    for (size_t i = 0; i < 100; i++) {
        s = (s * 48271UL) & 0x7fffffffUL;
        A[i] = ((double)s / 2147483647.0 - 0.5) * 10.0;
    }
    verify_general_eigenvectors(A, 10);
    printf("PASS: random 10x10 non-symmetric eigenvectors residual\n");
}

void test_direct_general_eigenvectors_5_2plusi(void) {
    /* {5, 2+i, 2-i} spectrum -- block-diag mixing real and complex. */
    double A[9] = {
        5.0, 0.0, 0.0,
        0.0, 2.0, -1.0,
        0.0, 1.0,  2.0
    };
    verify_general_eigenvectors(A, 3);
    printf("PASS: {5, 2+/-i} eigenvectors residual\n");
}

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

/* --- Phase 2 step 2c-A: complex Hermitian -------------------------- */

/* Eigenvalue tolerance for a complex Hermitian matrix.  Same Householder
 * + Wilkinson-QR constants as the real path, plus a small fudge for the
 * extra phase-correction step. */
static double machine_tol_complex(const double* A_re, const double* A_im,
                                    size_t n) {
    double normA = corpus_norm_inf_complex(A_re, A_im, n);
    if (normA == 0.0) normA = 1.0;
    return 128.0 * (double)n * 2.220446049250313e-16 * normA;
}

/* Verify every (lambda_i, v_i) pair for a complex Hermitian matrix A.
 * Checks: eigenvalues real, sum == real(trace), per-pair residual,
 * descending |lambda| sort order, unitary orthonormality of V. */
static void verify_hermitian_spectrum(const double* A_re, const double* A_im,
                                        size_t n) {
    double tol = machine_tol_complex(A_re, A_im, n);
    double* lambdas = NULL;
    double *V_re = NULL, *V_im = NULL;
    size_t kL = corpus_eval_hermitian_eigenvalues(A_re, A_im, n, &lambdas);
    size_t kV = corpus_eval_eigenvectors_complex(A_re, A_im, n, &V_re, &V_im);
    ASSERT(kL == n && kV == n);

    /* Sum of eigenvalues == Re(trace A).  Imag(trace A) is exactly zero
     * for Hermitian A, but we recover it from the matrix data anyway. */
    double trace = 0.0;
    for (size_t i = 0; i < n; i++) trace += A_re[i * n + i];
    double sum = 0.0;
    for (size_t i = 0; i < n; i++) sum += lambdas[i];
    if (fabs(sum - trace) > tol) {
        printf("FAIL: sum(lambda)=%g != trace=%g (tol %g)\n", sum, trace, tol);
        ASSERT(0);
    }

    /* Residual.  Hermitian eigenvalues are real; vectors complex. */
    for (size_t i = 0; i < n; i++) {
        corpus_assert_residual_complex_real_lambda(A_re, A_im, n,
                                                     lambdas[i],
                                                     &V_re[i * n],
                                                     &V_im[i * n],
                                                     tol);
    }

    /* Descending |lambda| sort. */
    for (size_t i = 1; i < n; i++) {
        ASSERT(fabs(lambdas[i - 1]) >= fabs(lambdas[i]) - tol);
    }

    /* Unitary orthonormality. */
    corpus_assert_unitary(V_re, V_im, n,
                            128.0 * (double)n * 2.220446049250313e-16);

    free(lambdas); free(V_re); free(V_im);
}

void test_direct_hermitian_2x2_diag_imag(void) {
    /* [[2, i], [-i, 3]] -- Hermitian, eigenvalues (5 +/- sqrt(5))/2. */
    double Re[4] = { 2.0, 0.0,   0.0, 3.0 };
    double Im[4] = { 0.0, 1.0,  -1.0, 0.0 };
    verify_hermitian_spectrum(Re, Im, 2);
    printf("PASS: 2x2 Hermitian [[2,i],[-i,3]] residual + unitary\n");
}

void test_direct_hermitian_pauli_x(void) {
    /* Pauli-X: [[0,1],[1,0]] -- real, eigenvalues +/-1. */
    double Re[4] = { 0.0, 1.0,  1.0, 0.0 };
    double Im[4] = { 0.0, 0.0,  0.0, 0.0 };
    verify_hermitian_spectrum(Re, Im, 2);
    printf("PASS: Pauli-X 2x2 residual + unitary\n");
}

void test_direct_hermitian_pauli_y(void) {
    /* Pauli-Y: [[0, -i], [i, 0]] -- Hermitian, eigenvalues +/-1.
     * Forces the phase-correction path through a 2x2 Hermitian step. */
    double Re[4] = { 0.0,  0.0,  0.0, 0.0 };
    double Im[4] = { 0.0, -1.0,  1.0, 0.0 };
    verify_hermitian_spectrum(Re, Im, 2);
    printf("PASS: Pauli-Y 2x2 residual + unitary\n");
}

void test_direct_hermitian_pauli_z(void) {
    /* Pauli-Z: [[1,0],[0,-1]] -- real, eigenvalues +/-1. */
    double Re[4] = { 1.0, 0.0,  0.0, -1.0 };
    double Im[4] = { 0.0, 0.0,  0.0,  0.0 };
    verify_hermitian_spectrum(Re, Im, 2);
    printf("PASS: Pauli-Z 2x2 residual + unitary\n");
}

void test_direct_hermitian_3x3_known(void) {
    /* A canonical 3x3 Hermitian.
     *   A = [[1, 1+i, 0],
     *        [1-i,  3,   i],
     *        [ 0,  -i,   2]]
     * trace = 6.  Verified Hermitian by inspection. */
    double Re[9] = {
        1.0, 1.0,  0.0,
        1.0, 3.0,  0.0,
        0.0, 0.0,  2.0
    };
    double Im[9] = {
        0.0,  1.0, 0.0,
       -1.0,  0.0, 1.0,
        0.0, -1.0, 0.0
    };
    verify_hermitian_spectrum(Re, Im, 3);
    printf("PASS: 3x3 Hermitian (sum==trace==6) residual + unitary\n");
}

void test_direct_hermitian_diagonal_4x4(void) {
    /* Diagonal real, but loaded into the complex path via a zero im[] --
     * exercises the Hermitian dispatch + early-exit path (sigma == 0
     * at every step). */
    double Re[16] = {
        4.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 3.0, 0.0,
        0.0, 0.0, 0.0, 2.0
    };
    double Im[16] = {0};
    /* This routes through the real-symmetric path (no complex entries).
     * We assert directly that the eigenvalues come out in descending |lambda|
     * order: {4, 3, 2, 1}. */
    double* lambdas = NULL;
    size_t k = corpus_eval_hermitian_eigenvalues(Re, Im, 4, &lambdas);
    ASSERT(k == 4);
    const double expected[4] = { 4.0, 3.0, 2.0, 1.0 };
    for (size_t i = 0; i < 4; i++) {
        ASSERT(fabs(lambdas[i] - expected[i]) < 1e-12);
    }
    free(lambdas);
    printf("PASS: diag(4,1,3,2) Hermitian dispatch returns sorted real eigenvalues\n");
}

void test_direct_hermitian_5x5_random(void) {
    /* Build a 5x5 Hermitian matrix M = R + R^T  + i(S - S^T)/2.
     * Uses a small reproducible "random" seed via a linear congruential
     * generator so the test is independent of test order. */
    uint32_t s = 0xC0FFEE01u;
    #define LCG() (s = s * 1103515245u + 12345u, (double)(s >> 16) / 32768.0 - 1.0)
    double Re[25], Im[25];
    for (size_t i = 0; i < 5; i++) {
        for (size_t j = 0; j < 5; j++) {
            double a = LCG();
            double b = LCG();
            if (i == j) {
                Re[i * 5 + j] = 2.0 * a;
                Im[i * 5 + j] = 0.0;
            } else if (i < j) {
                Re[i * 5 + j] = a;
                Im[i * 5 + j] = b;
            } else {
                /* Conjugate of (i,j) is (j,i). */
                Re[i * 5 + j] = Re[j * 5 + i];
                Im[i * 5 + j] = -Im[j * 5 + i];
            }
        }
    }
    #undef LCG
    verify_hermitian_spectrum(Re, Im, 5);
    printf("PASS: random 5x5 Hermitian residual + unitary\n");
}

void test_direct_hermitian_repeated_eigvals(void) {
    /* I_4 + i * 0 = identity.  Eigenvalues all 1; eigenvectors should
     * form an orthonormal basis (any one will do since they live in a
     * 4-D degenerate eigenspace).  Exercises the kernel on a maximally
     * degenerate input. */
    double Re[16] = {0}; double Im[16] = {0};
    for (size_t i = 0; i < 4; i++) Re[i * 4 + i] = 1.0;
    double* lambdas = NULL;
    double *V_re = NULL, *V_im = NULL;
    size_t kL = corpus_eval_hermitian_eigenvalues(Re, Im, 4, &lambdas);
    size_t kV = corpus_eval_eigenvectors_complex(Re, Im, 4, &V_re, &V_im);
    ASSERT(kL == 4 && kV == 4);
    for (size_t i = 0; i < 4; i++) ASSERT(fabs(lambdas[i] - 1.0) < 1e-12);
    /* Per-eigenpair residual -- A v_i = v_i for any v_i. */
    for (size_t i = 0; i < 4; i++) {
        corpus_assert_residual_complex_real_lambda(Re, Im, 4, lambdas[i],
                                                     &V_re[i * 4], &V_im[i * 4],
                                                     1e-12);
    }
    corpus_assert_unitary(V_re, V_im, 4, 1e-12);
    free(lambdas); free(V_re); free(V_im);
    printf("PASS: 4x4 identity Hermitian residual + unitary (degenerate)\n");
}

void test_direct_hermitian_complex_offdiag_only(void) {
    /* [[0, 2i], [-2i, 0]] -- pure imag off-diagonal.  Eigenvalues +/-2,
     * eigenvectors must come out complex. */
    double Re[4] = { 0.0, 0.0,  0.0, 0.0 };
    double Im[4] = { 0.0, 2.0, -2.0, 0.0 };
    double* lambdas = NULL;
    double *V_re = NULL, *V_im = NULL;
    size_t kL = corpus_eval_hermitian_eigenvalues(Re, Im, 2, &lambdas);
    size_t kV = corpus_eval_eigenvectors_complex(Re, Im, 2, &V_re, &V_im);
    ASSERT(kL == 2 && kV == 2);
    /* Eigenvalues +/-2 sorted by abs desc: {2, -2} or {-2, 2}; the +2
     * is the first by stable tiebreak in our sort. */
    ASSERT(fabs(fabs(lambdas[0]) - 2.0) < 1e-12);
    ASSERT(fabs(fabs(lambdas[1]) - 2.0) < 1e-12);
    /* Verify each eigenvector via residual A v = lambda v. */
    for (size_t i = 0; i < 2; i++) {
        corpus_assert_residual_complex_real_lambda(Re, Im, 2, lambdas[i],
                                                     &V_re[i * 2], &V_im[i * 2],
                                                     1e-12);
    }
    /* Eigenvectors should be orthogonal in C^2. */
    corpus_assert_unitary(V_re, V_im, 2, 1e-12);
    free(lambdas); free(V_re); free(V_im);
    printf("PASS: 2x2 Hermitian [[0,2i],[-2i,0]] residual + unitary\n");
}

void test_direct_hermitian_eigenvalues_only(void) {
    /* Exercise the WANT_VALUES (no Q) path to make sure the kernel
     * doesn't crash when eigenvectors aren't requested.  Same matrix
     * as 3x3_known so we get the same trace = 6. */
    double Re[9] = {
        1.0, 1.0,  0.0,
        1.0, 3.0,  0.0,
        0.0, 0.0,  2.0
    };
    double Im[9] = {
        0.0,  1.0, 0.0,
       -1.0,  0.0, 1.0,
        0.0, -1.0, 0.0
    };
    double* lambdas = NULL;
    size_t k = corpus_eval_hermitian_eigenvalues(Re, Im, 3, &lambdas);
    ASSERT(k == 3);
    double sum = lambdas[0] + lambdas[1] + lambdas[2];
    ASSERT(fabs(sum - 6.0) < 1e-12);
    /* Eigenvalues should be sorted in descending |lambda|. */
    ASSERT(fabs(lambdas[0]) >= fabs(lambdas[1]) - 1e-12);
    ASSERT(fabs(lambdas[1]) >= fabs(lambdas[2]) - 1e-12);
    free(lambdas);
    printf("PASS: 3x3 Hermitian eigenvalues-only path\n");
}

/* --- Phase 2 step 2c-B: complex non-Hermitian ---------------------- */

/* Tolerance for the complex general path.  The block-embedding doubles
 * the matrix dimension and adds one Schur back-substitution pass, so we
 * loosen the constant beyond the Hermitian (128) path. */
static double machine_tol_complex_general(const double* A_re,
                                            const double* A_im, size_t n) {
    double normA = corpus_norm_inf_complex(A_re, A_im, n);
    if (normA == 0.0) normA = 1.0;
    return 256.0 * (double)n * 2.220446049250313e-16 * normA;
}

/* Check the eigenvalue residual for every (lambda_i, v_i) pair of a
 * complex (possibly non-Hermitian) matrix.  Unlike Hermitian, vectors
 * are NOT expected to be orthonormal -- only unit-norm.  Eigenvalues
 * are complex in general (no real check). */
static void verify_complex_general_spectrum(const double* A_re,
                                              const double* A_im, size_t n) {
    double tol = machine_tol_complex_general(A_re, A_im, n);

    double *eval_re = NULL, *eval_im = NULL;
    size_t kL = corpus_eval_eigenvalues_complex(A_re, A_im, n,
                                                  &eval_re, &eval_im);
    double *V_re = NULL, *V_im = NULL;
    size_t kV = corpus_eval_eigenvectors_complex(A_re, A_im, n, &V_re, &V_im);
    ASSERT(kL == n && kV == n);

    /* Eigenvalue sum == trace.  trace = sum_i A_ii. */
    double tr_re = 0.0, tr_im = 0.0;
    for (size_t i = 0; i < n; i++) {
        tr_re += A_re[i * n + i];
        tr_im += A_im[i * n + i];
    }
    double s_re = 0.0, s_im = 0.0;
    for (size_t i = 0; i < n; i++) { s_re += eval_re[i]; s_im += eval_im[i]; }
    if (hypot(s_re - tr_re, s_im - tr_im) > tol) {
        printf("FAIL: sum(lambda)=%g+%gi != trace=%g+%gi (tol %g)\n",
               s_re, s_im, tr_re, tr_im, tol);
        ASSERT(0);
    }

    /* Per-pair residual. */
    for (size_t i = 0; i < n; i++) {
        corpus_assert_residual_complex(A_re, A_im, n,
                                         eval_re[i], eval_im[i],
                                         &V_re[i * n], &V_im[i * n],
                                         tol);
    }

    /* Descending |lambda| sort. */
    for (size_t i = 1; i < n; i++) {
        double prev = hypot(eval_re[i - 1], eval_im[i - 1]);
        double cur  = hypot(eval_re[i],     eval_im[i]);
        ASSERT(prev >= cur - tol);
    }

    /* Eigenvectors unit-norm (NOT necessarily orthogonal). */
    for (size_t i = 0; i < n; i++) {
        double sq = 0.0;
        for (size_t j = 0; j < n; j++) {
            double r = V_re[i * n + j];
            double m = V_im[i * n + j];
            sq += r * r + m * m;
        }
        ASSERT(fabs(sq - 1.0) < 1e-10);
    }

    free(eval_re); free(eval_im); free(V_re); free(V_im);
}

void test_direct_complex_general_diagonal(void) {
    /* A = diag(1+i, 2+i) -- non-Hermitian (Im part is on the diagonal,
     * which violates Hermitian).  Eigenvalues 1+i, 2+i.  No mixing
     * (conj(spec) is disjoint from spec). */
    double Re[4] = { 1.0, 0.0,  0.0, 2.0 };
    double Im[4] = { 1.0, 0.0,  0.0, 1.0 };
    verify_complex_general_spectrum(Re, Im, 2);
    printf("PASS: diag(1+i, 2+i) complex general\n");
}

void test_direct_complex_general_mixing(void) {
    /* A = [[0, 1, 0], [-1, 0, 0], [0, 0, 2+i]] -- the canonical "mixing"
     * test: spec(A) = {i, -i, 2+i}, so conj(spec) and spec overlap on
     * {i, -i}.  Simple pair-up of M's eigenvalues would emit two copies
     * of i (or two of -i) here; grouped Gram-Schmidt is required. */
    double Re[9] = {
        0.0, 1.0, 0.0,
       -1.0, 0.0, 0.0,
        0.0, 0.0, 2.0
    };
    double Im[9] = {0};
    Im[2 * 3 + 2] = 1.0;
    verify_complex_general_spectrum(Re, Im, 3);
    printf("PASS: [[0,1,0],[-1,0,0],[0,0,2+i]] mixing case (grouped GS)\n");
}

void test_direct_complex_general_2x2_generic(void) {
    /* Generic 2x2 complex.  Char poly: lambda^2 - 5 lambda - 3 - i = 0.
     * Discriminant 37 + 4i, eigenvalues (5 +/- sqrt(37+4i))/2. */
    double Re[4] = { 1.0, 2.0,  3.0, 4.0 };
    double Im[4] = { 0.0, 1.0, -1.0, 0.0 };
    verify_complex_general_spectrum(Re, Im, 2);
    printf("PASS: generic complex 2x2 residual + sort + unit-norm\n");
}

void test_direct_complex_general_real_off_complex_diag(void) {
    /* A = [[1+i, 2], [3, 4]] -- only one complex entry; eigenvalues
     * generically complex. */
    double Re[4] = { 1.0, 2.0,  3.0, 4.0 };
    double Im[4] = { 1.0, 0.0,  0.0, 0.0 };
    verify_complex_general_spectrum(Re, Im, 2);
    printf("PASS: mixed-entries 2x2 complex general\n");
}

void test_direct_complex_general_3x3_random(void) {
    /* Reproducible "random" complex 3x3 via the LCG seed already used
     * elsewhere in this file. */
    uint32_t s = 0xC0FFEE07u;
    #define LCG() (s = s * 1103515245u + 12345u, (double)(s >> 16) / 32768.0 - 1.0)
    double Re[9], Im[9];
    for (size_t i = 0; i < 9; i++) { Re[i] = LCG(); Im[i] = LCG(); }
    #undef LCG
    verify_complex_general_spectrum(Re, Im, 3);
    printf("PASS: random complex 3x3 residual\n");
}

void test_direct_complex_general_5x5_random(void) {
    uint32_t s = 0xD15EA5Eu;
    #define LCG() (s = s * 1103515245u + 12345u, (double)(s >> 16) / 32768.0 - 1.0)
    double Re[25], Im[25];
    for (size_t i = 0; i < 25; i++) { Re[i] = LCG(); Im[i] = LCG(); }
    #undef LCG
    verify_complex_general_spectrum(Re, Im, 5);
    printf("PASS: random complex 5x5 residual\n");
}

void test_direct_complex_general_real_rotation_promoted(void) {
    /* [[0, 1], [-1, 0]] re-loaded into the complex path by adding a
     * zero-im entry that triggers is_complex=true.  Eigenvalues +/-i.
     * spec(A) and conj(spec(A)) coincide so we hit the mixing path. */
    double Re[4] = { 0.0, 1.0,  -1.0, 0.0 };
    double Im[4] = { 0.0, 0.0,  0.0, 0.0 };
    /* Stamp a complex 0 to force is_complex routing. */
    Im[0] = 1e-300;
    Im[0] = 0.0;
    /* The above doesn't actually mark complex; instead use a tiny but
     * non-zero imag entry on the diagonal so dispatch hits the complex
     * path.  We then assert the recovered eigenvalues are nearly +/-i. */
    Im[0] = 1e-15;
    double *eval_re = NULL, *eval_im = NULL;
    size_t k = corpus_eval_eigenvalues_complex(Re, Im, 2,
                                                 &eval_re, &eval_im);
    ASSERT(k == 2);
    /* One should be near +i, the other near -i. */
    double a0 = hypot(eval_re[0] - 0.0, eval_im[0] - 1.0);
    double b0 = hypot(eval_re[0] - 0.0, eval_im[0] + 1.0);
    double pick0 = (a0 < b0) ? 1.0 : -1.0;
    double pick1 = -pick0;
    ASSERT(fabs(eval_im[0] - pick0) < 1e-6);
    ASSERT(fabs(eval_im[1] - pick1) < 1e-6);
    ASSERT(fabs(eval_re[0]) < 1e-6);
    ASSERT(fabs(eval_re[1]) < 1e-6);
    free(eval_re); free(eval_im);
    printf("PASS: tiny-imag rotation matrix -> complex eigenvalues +/-i\n");
}

void test_direct_complex_general_eigenvalues_only(void) {
    /* Exercise the WANT_VALUES (no Q) path of the complex general
     * kernel.  Even values-only needs M's eigenvectors internally for
     * the +J / -J disambiguation, so this exercises the "compute Q
     * unconditionally then throw it away" branch. */
    double Re[4] = { 1.0, 2.0,  3.0, 4.0 };
    double Im[4] = { 0.0, 1.0, -1.0, 0.0 };
    double *eval_re = NULL, *eval_im = NULL;
    size_t k = corpus_eval_eigenvalues_complex(Re, Im, 2,
                                                 &eval_re, &eval_im);
    ASSERT(k == 2);
    /* Eigenvalues sum to trace = 5 (the imag part of trace is 0). */
    double sr = eval_re[0] + eval_re[1];
    double si = eval_im[0] + eval_im[1];
    ASSERT(fabs(sr - 5.0) < 1e-10);
    ASSERT(fabs(si - 0.0) < 1e-10);
    free(eval_re); free(eval_im);
    printf("PASS: complex 2x2 eigenvalues-only path\n");
}

#ifdef USE_MPFR
/* Helper: parse + evaluate `src`, assert the result prints as "True".
 *
 * Side note on Mathilda comparison semantics: `Less[mpfr_a, mpfr_b]`
 * will *not* always evaluate to True/False if either side is an exact
 * zero (machine or MPFR).  The tests below sidestep this by scaling
 * tiny MPFR residuals into the machine-precision range via `N[...]`
 * before comparing -- the scaled value is well clear of zero and Less
 * is happy. */
static void mpfr_assert_true(const char* label, const char* src) {
    Expr* e = parse_expression(src);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    char* s = expr_to_string(r);
    if (strcmp(s, "True") != 0) {
        printf("FAIL: %s\n  src: %s\n  got: %s\n", label, src, s);
        ASSERT(0);
    } else {
        printf("PASS: %s\n", label);
    }
    free(s);
    expr_free(r);
    expr_free(e);
}

/* Helper: assert that the head of the evaluated `src` is `head_name`. */
static void mpfr_assert_head(const char* label, const char* src,
                              const char* head_name) {
    Expr* e = parse_expression(src);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    bool ok = (r->type == EXPR_FUNCTION
               && r->data.function.head->type == EXPR_SYMBOL
               && strcmp(r->data.function.head->data.symbol, head_name) == 0);
    if (!ok) {
        char* s = expr_to_string(r);
        printf("FAIL: %s (head)\n  src: %s\n  got: %s\n", label, src, s);
        free(s);
        ASSERT(0);
    }
    printf("PASS: %s (head=%s)\n", label, head_name);
    expr_free(r); expr_free(e);
}

/* Verify residual ||A v - lambda v|| is accurate to at least `digits`
 * decimal places.  Strategy: scale the residual by 10^digits and then
 * compare its machine-precision norm to 1.0 -- a true MPFR residual of
 * 1e-{digits+1} becomes ~0.1 after scaling, exposing any drift.
 *
 * NOTE: identifiers passed to Mathilda's parser MUST NOT contain
 * underscores; in Mathematica syntax `max_norm` parses as
 * `Pattern[max, Blank[norm]]`, which silently breaks Module
 * declarations.  Use mixedCase / camelCase identifiers throughout. */
static void mpfr_assert_eigenpair_residuals(const char* label,
                                              const char* matrix_src,
                                              int digits) {
    char buf[2048];
    snprintf(buf, sizeof(buf),
        "Module[{m, vals, vecs, n, k, scaled, mxN, normK},\n"
        "  m = %s;\n"
        "  vals = Eigenvalues[m];\n"
        "  vecs = Eigenvectors[m];\n"
        "  n = Length[vals];\n"
        "  mxN = 0.0;\n"
        "  k = 1;\n"
        "  While[k <= n,\n"
        "    scaled = (Dot[m, vecs[[k]]] - vals[[k]] vecs[[k]]) * 10^%d;\n"
        "    normK = N[Norm[scaled]];\n"
        "    If[normK > mxN, mxN = normK];\n"
        "    k = k + 1];\n"
        "  TrueQ[mxN < 1.0]]",
        matrix_src, digits);
    mpfr_assert_true(label, buf);
}

/* Verify orthonormality: scale ||V V^T - I||_inf by 10^digits and
 * confirm it's below 1.0 at machine precision. */
static void mpfr_assert_orthonormal(const char* label,
                                      const char* matrix_src,
                                      int digits) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
        "Module[{m, vecs, n, prod, mxErr, scaled, k, jj, entry},\n"
        "  m = %s;\n"
        "  vecs = Eigenvectors[m];\n"
        "  n = Length[vecs];\n"
        "  prod = Dot[vecs, Transpose[vecs]];\n"
        "  mxErr = 0.0;\n"
        "  k = 1;\n"
        "  While[k <= n,\n"
        "    jj = 1;\n"
        "    While[jj <= n,\n"
        "      scaled = (prod[[k]][[jj]] - If[k == jj, 1, 0]) * 10^%d;\n"
        "      entry = N[Abs[scaled]];\n"
        "      If[entry > mxErr, mxErr = entry];\n"
        "      jj = jj + 1];\n"
        "    k = k + 1];\n"
        "  TrueQ[mxErr < 1.0]]",
        matrix_src, digits);
    mpfr_assert_true(label, buf);
}

/* Verify Precision[lambda_i] > min_decimal_prec for every eigenvalue,
 * i.e. MPFR routing happened (a slip to machine would yield
 * Precision == MachinePrecision ~ 15.95 digits).  Mathilda's NumberQ /
 * Head don't fire on MPFR values, so we lean on Precision alone -- it
 * returns a Real and Less is well-behaved there. */
static void mpfr_assert_eigenvalues_have_precision(const char* label,
                                                     const char* matrix_src,
                                                     int n,
                                                     int min_decimal_prec) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
        "Module[{m, vs, ok, k},\n"
        "  m = %s;\n"
        "  vs = Eigenvalues[m];\n"
        "  ok = TrueQ[Length[vs] == %d];\n"
        "  k = 1;\n"
        "  While[k <= Length[vs] && ok,\n"
        "    ok = ok && TrueQ[Precision[vs[[k]]] > %d];\n"
        "    k = k + 1];\n"
        "  ok]",
        matrix_src, n, min_decimal_prec);
    mpfr_assert_true(label, buf);
}

/* --- MPFR real symmetric tests (step 2d-A) -----------------------------
 *
 * All "small residual" tests use the scale-by-10^k + N[]-to-machine
 * pattern: multiply the (MPFR) residual by 10^k so its magnitude moves
 * into the comfortable machine-precision range, convert via N[] to a
 * plain double, then compare to 1.0.  This sidesteps Mathilda's
 * Less[machine_zero, mpfr_value] non-evaluation quirk and gives us a
 * sharp digit count: passing scale=40 with threshold 1.0 means the
 * residual is bounded by 10^-40. */

void test_mpfr_sym_2x2_golden_ratio(void) {
    /* {{1,1},{1,2}} -- eigenvalues are phi and 1/phi == (3 +/- Sqrt[5])/2.
     * At 50 decimal digits we expect ~10^-49 residual. */
    const char* src =
        "Module[{m, vs, e1, e2, d1, d2},\n"
        "  m = SetPrecision[{{1, 1}, {1, 2}}, 50];\n"
        "  vs = Eigenvalues[m];\n"
        "  e1 = SetPrecision[(3 + Sqrt[5])/2, 50];\n"
        "  e2 = SetPrecision[(3 - Sqrt[5])/2, 50];\n"
        "  d1 = N[(vs[[1]] - e1) * 10^45];\n"
        "  d2 = N[(vs[[2]] - e2) * 10^45];\n"
        "  TrueQ[N[Abs[d1]] < 1.0] && TrueQ[N[Abs[d2]] < 1.0]]";
    mpfr_assert_true("MPFR sym 2x2 golden ratio (50 digits, scale=45)", src);
}

void test_mpfr_sym_3x3_tridiag_known(void) {
    /* Tridiag(2,3,2) with off-diag 1 -- analytic eigenvalues 4, 2, 1. */
    const char* src =
        "Module[{m, vs, d1, d2, d3},\n"
        "  m = SetPrecision[{{2, 1, 0}, {1, 3, 1}, {0, 1, 2}}, 80];\n"
        "  vs = Eigenvalues[m];\n"
        "  d1 = N[(vs[[1]] - 4) * 10^70];\n"
        "  d2 = N[(vs[[2]] - 2) * 10^70];\n"
        "  d3 = N[(vs[[3]] - 1) * 10^70];\n"
        "  TrueQ[N[Abs[d1]] < 1.0]\n"
        " && TrueQ[N[Abs[d2]] < 1.0]\n"
        " && TrueQ[N[Abs[d3]] < 1.0]]";
    mpfr_assert_true("MPFR sym 3x3 tridiag(2,3,2) -> {4,2,1} at 80 digits", src);
}

void test_mpfr_sym_diagonal_identity(void) {
    /* Identity at high precision -- degenerate spectrum; eigenvalues all 1. */
    const char* src =
        "Module[{m, vs, d},\n"
        "  m = SetPrecision[{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, 50];\n"
        "  vs = Eigenvalues[m];\n"
        "  Length[vs] == 3\n"
        " && TrueQ[N[Abs[(vs[[1]] - 1) * 10^45]] < 1.0]\n"
        " && TrueQ[N[Abs[(vs[[2]] - 1) * 10^45]] < 1.0]\n"
        " && TrueQ[N[Abs[(vs[[3]] - 1) * 10^45]] < 1.0]]";
    mpfr_assert_true("MPFR sym 3x3 identity -> {1,1,1}", src);
}

void test_mpfr_sym_5x5_residual(void) {
    /* Symmetric 5x5 at 60 digits via (B + B^T)/2; verify per-eigenpair
     * residual and orthonormality at the 45-digit level. */
    const char* matrix_src =
        "Module[{b, sym},\n"
        "  b = SetPrecision[{{2, 1, 3, 0, 1},\n"
        "                     {1, 5, 2, 1, 0},\n"
        "                     {3, 2, 7, 1, 2},\n"
        "                     {0, 1, 1, 4, 1},\n"
        "                     {1, 0, 2, 1, 6}}, 60];\n"
        "  sym = (b + Transpose[b])/2;\n"
        "  sym]";
    mpfr_assert_eigenpair_residuals(
        "MPFR sym 5x5 eigenpair residual at 45 digits (input 60)",
        matrix_src, 45);
    mpfr_assert_orthonormal(
        "MPFR sym 5x5 orthonormality at 45 digits (input 60)",
        matrix_src, 45);
}

void test_mpfr_sym_eigenvalues_carry_precision(void) {
    /* Confirms MPFR routing fired -- a slip to the machine kernel would
     * yield Precision[lambda] == MachinePrecision (~15.95). */
    mpfr_assert_eigenvalues_have_precision(
        "MPFR sym 2x2 eigenvalues carry > 49 digits",
        "SetPrecision[{{1, 1}, {1, 2}}, 50]", 2, 49);
}

void test_mpfr_sym_eigenvectors_are_lists(void) {
    mpfr_assert_head("MPFR sym Eigenvectors returns a List",
        "Eigenvectors[SetPrecision[{{1, 1}, {1, 2}}, 40]]", "List");
}

void test_mpfr_sym_high_precision_200(void) {
    /* Stress: 200 decimal digits.  Verify residual at 180 digits. */
    const char* matrix_src = "SetPrecision[{{4, 1, 0}, {1, 4, 1}, {0, 1, 4}}, 200]";
    mpfr_assert_eigenpair_residuals(
        "MPFR sym 3x3 at 200 digits, residual < 10^-180",
        matrix_src, 180);
}

void test_mpfr_sym_k_spec(void) {
    /* k-spec interacts with the MPFR kernel: top 2 of diag(5,3,1). */
    const char* src =
        "Module[{m, vs},\n"
        "  m = SetPrecision[{{5, 0, 0}, {0, 3, 0}, {0, 0, 1}}, 50];\n"
        "  vs = Eigenvalues[m, 2];\n"
        "  Length[vs] == 2\n"
        " && TrueQ[N[Abs[(vs[[1]] - 5) * 10^45]] < 1.0]\n"
        " && TrueQ[N[Abs[(vs[[2]] - 3) * 10^45]] < 1.0]]";
    mpfr_assert_true("MPFR sym k-spec: top 2 of diag(5,3,1)", src);
}

void test_mpfr_sym_repeated_eigvals(void) {
    /* Near-degenerate 2I + small symmetric perturbation: three real
     * eigenvalues all close to 2. */
    const char* src =
        "Module[{m, vs, ok, i, ds},\n"
        "  m = SetPrecision[{{2, 1/100, 0}, {1/100, 2, 1/100}, {0, 1/100, 2}}, 60];\n"
        "  vs = Eigenvalues[m];\n"
        "  ok = TrueQ[Length[vs] == 3];\n"
        "  For[i = 1, i <= 3, i++,\n"
        "    ds = N[Abs[vs[[i]] - 2]];\n"
        "    ok = ok && TrueQ[ds < 0.05]];\n"
        "  ok]";
    mpfr_assert_true("MPFR sym near-degenerate {2,2,2}+epsilon", src);
}

void test_mpfr_sym_1x1(void) {
    /* Edge case: 1x1 -- bypasses tridiag entirely. */
    const char* src =
        "Module[{m, vs, vc},\n"
        "  m = SetPrecision[{{3}}, 50];\n"
        "  vs = Eigenvalues[m];\n"
        "  vc = Eigenvectors[m];\n"
        "  Length[vs] == 1\n"
        " && TrueQ[N[Abs[(vs[[1]] - 3) * 10^45]] < 1.0]\n"
        " && Length[vc] == 1 && Length[vc[[1]]] == 1\n"
        " && TrueQ[N[Abs[(Abs[vc[[1]][[1]]] - 1) * 10^45]] < 1.0]]";
    mpfr_assert_true("MPFR sym 1x1 edge case", src);
}
#endif  /* USE_MPFR */

int main(void) {
    symtab_init();
    core_init();
    printf("Running mateigen Direct (Phase 2: machine-precision) tests...\n");

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
    TEST(test_direct_upper_triangular_real_eigenvalues);
    TEST(test_direct_symbolic_ignores_dispatch);

    /* Phase 2 step 2a: real non-symmetric eigenvalues. */
    TEST(test_direct_general_2x2_real);
    TEST(test_direct_general_2x2_complex_rotation);
    TEST(test_direct_general_2x2_complex_general);
    TEST(test_direct_general_3x3_block_diagonal);
    TEST(test_direct_general_3x3_random);
    TEST(test_direct_general_5x5_random);
    TEST(test_direct_general_10x10_random);
    TEST(test_direct_general_companion_matrix);
    TEST(test_direct_general_jordan_block_2x2);
    TEST(test_direct_general_complex_pair_sorted_after_real);
    TEST(test_direct_general_method_direct_routed);

    /* Phase 2 step 2b: real non-symmetric eigenvectors (real and complex
     * eigenvalues, residual-checked through the actual matrix A). */
    TEST(test_direct_general_eigenvectors_2x2);
    TEST(test_direct_general_eigenvectors_rotation);
    TEST(test_direct_general_eigenvectors_3x3_complex_block);
    TEST(test_direct_general_eigenvectors_upper_triangular);
    TEST(test_direct_general_eigenvectors_generic_3x3);
    TEST(test_direct_general_eigenvectors_5x5_random);
    TEST(test_direct_general_eigenvectors_10x10_random);
    TEST(test_direct_general_eigenvectors_5_2plusi);
    TEST(test_direct_eigenvectors_diag_4x4);
    TEST(test_direct_eigenvectors_tridiag_residual);

    /* Phase 2 step 2c-A: complex Hermitian. */
    TEST(test_direct_hermitian_2x2_diag_imag);
    TEST(test_direct_hermitian_pauli_x);
    TEST(test_direct_hermitian_pauli_y);
    TEST(test_direct_hermitian_pauli_z);
    TEST(test_direct_hermitian_3x3_known);
    TEST(test_direct_hermitian_diagonal_4x4);
    TEST(test_direct_hermitian_5x5_random);
    TEST(test_direct_hermitian_repeated_eigvals);
    TEST(test_direct_hermitian_complex_offdiag_only);
    TEST(test_direct_hermitian_eigenvalues_only);

    /* Phase 2 step 2c-B: complex non-Hermitian. */
    TEST(test_direct_complex_general_diagonal);
    TEST(test_direct_complex_general_mixing);
    TEST(test_direct_complex_general_2x2_generic);
    TEST(test_direct_complex_general_real_off_complex_diag);
    TEST(test_direct_complex_general_3x3_random);
    TEST(test_direct_complex_general_5x5_random);
    TEST(test_direct_complex_general_real_rotation_promoted);
    TEST(test_direct_complex_general_eigenvalues_only);

#ifdef USE_MPFR
    /* Phase 2 step 2d-A: real symmetric MPFR. */
    TEST(test_mpfr_sym_2x2_golden_ratio);
    TEST(test_mpfr_sym_3x3_tridiag_known);
    TEST(test_mpfr_sym_diagonal_identity);
    TEST(test_mpfr_sym_5x5_residual);
    TEST(test_mpfr_sym_eigenvalues_carry_precision);
    TEST(test_mpfr_sym_eigenvectors_are_lists);
    TEST(test_mpfr_sym_high_precision_200);
    TEST(test_mpfr_sym_k_spec);
    TEST(test_mpfr_sym_repeated_eigvals);
    TEST(test_mpfr_sym_1x1);
#endif

    printf("All mateigen Direct (machine-precision) tests passed!\n");
    return 0;
}
