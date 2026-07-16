/*
 * test_lapack_builtin.c -- unit tests for the LAPACK` context builtins.
 *
 * Distinct from test_lapack.c (which tests the low-level linkage probe).
 * Results are cross-checked in-language against Mathilda's own high-level
 * operators (LinearSolve, Inverse, Eigenvalues, Dot, ...) via a `closeQ`
 * tolerance predicate, so the tests are independent of float printing.
 *
 * When built without LAPACK the LAPACK` builtins are not registered and the
 * suite skips (exit 0).
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

static char* eval_str(const char* expr) {
    Expr* p = parse_expression(expr);
    ASSERT(p != NULL);
    Expr* e = evaluate(p);
    char* s = expr_to_string(e);
    expr_free(p);
    expr_free(e);
    return s;
}

static void run(const char* expr) { char* s = eval_str(expr); free(s); }

static void check_true(const char* expr) {
    char* s = eval_str(expr);
    if (strcmp(s, "True") != 0) {
        fprintf(stderr, "FAIL: %s\n  expected True, got: %s\n", expr, s);
        free(s);
        exit(1);
    }
    free(s);
}

static int lapack_available(void) {
    char* s = eval_str("Head[LAPACK`dgesv[{{2,0},{0,2}},{2,2}]]");
    int ok = (strcmp(s, "NDArray") == 0 || strcmp(s, "List") == 0);
    free(s);
    return ok;
}

static void test_solve(void) {
    check_true("closeQ[LAPACK`dgesv[{{3,2},{1,2}},{7,5}], "
               "LinearSolve[{{3,2},{1,2}},{7,5}]]");
    /* Matrix RHS. */
    check_true("closeQ[LAPACK`dgesv[{{3,2},{1,2}},{{7},{5}}], "
               "LinearSolve[{{3,2},{1,2}},{{7},{5}}]]");
    /* Triangular solve: A.X == B. */
    check_true("With[{A={{2,1},{0,3}}},"
               "closeQ[A.Normal[LAPACK`dtrtrs[A,{1,2}]], {1,2}]]");
    /* Least squares (overdetermined) matches the normal equations. */
    check_true("With[{A={{1,0},{0,1},{1,1}},b={1,2,4}},"
               "closeQ[LAPACK`dgels[A,b], "
               "LinearSolve[Transpose[A].A, Transpose[A].b]]]");
}

static void test_factor(void) {
    /* LU of a diagonal matrix needs no pivoting. */
    check_true("closeQ[Normal[First[LAPACK`dgetrf[{{2,0},{0,3}}]]], {{2,0},{0,3}}]");
    check_true("Last[LAPACK`dgetrf[{{2,0},{0,3}}]] == {1,2}");
    /* Inverse matches Inverse[]. */
    check_true("closeQ[LAPACK`dgetri[{{4,3},{6,3}}], Inverse[{{4,3},{6,3}}]]");
    /* QR: A == Q.R and Q has orthonormal columns. */
    check_true("Module[{q,r},{q,r}=Normal/@LAPACK`dgeqrf[{{1,2},{3,4}}];"
               "closeQ[q.r,{{1,2},{3,4}}] && "
               "closeQ[Transpose[q].q, IdentityMatrix[2]]]");
    /* Cholesky: L.L^T == A for an SPD matrix. */
    check_true("Module[{l},l=Normal[LAPACK`dpotrf[{{4,2},{2,3}}]];"
               "closeQ[l.Transpose[l], {{4,2},{2,3}}]]");
}

static void test_svd(void) {
    check_true("Module[{u,d,vt},{u,d,vt}=Normal/@LAPACK`dgesdd[{{1,2},{3,4}}];"
               "closeQ[u.DiagonalMatrix[d].vt, {{1,2},{3,4}}]]");
    check_true("Module[{u,d,vt},{u,d,vt}=Normal/@LAPACK`dgesvd[{{2,0},{0,3}}];"
               "closeQ[u.DiagonalMatrix[d].vt, {{2,0},{0,3}}]]");
    /* Singular values are nonnegative and descending. */
    check_true("Module[{d},d=Normal[LAPACK`dgesdd[{{1,2},{3,4}}][[2]]];"
               "d[[1]] >= d[[2]] && d[[2]] >= 0]");
}

static void test_eigen(void) {
    /* General eigenvalues (sorted) match Eigenvalues. */
    check_true("closeQ[Sort[First[LAPACK`dgeev[{{2,1},{1,2}}]]], "
               "Sort[N[Eigenvalues[{{2,1},{1,2}}]]]]");
    /* A rotation has purely imaginary eigenvalues +/- I. */
    check_true("closeQ[Sort[First[LAPACK`dgeev[{{0,-1},{1,0}}]]], Sort[{I,-I}]]");
    /* Symmetric eigenvalues come out ascending and match Eigenvalues. */
    check_true("closeQ[First[LAPACK`dsyev[{{2,1},{1,2}}]], "
               "Sort[Eigenvalues[{{2,1},{1,2}}]]]");
    /* dgeev eigenvectors satisfy A.v == lambda v. */
    check_true("Module[{vals,vecs},{vals,vecs}=LAPACK`dgeev[{{2,1},{1,2}}];"
               "closeQ[{{2,1},{1,2}}.vecs[[1]], vals[[1]] vecs[[1]]]]");
}

static void test_norm(void) {
    check_true("LAPACK`dlange[{{3,4}}] == 5");
    check_true("closeQ[LAPACK`dlange[{{1,2},{3,4}}], Sqrt[30]]");
}

static void test_complex(void) {
    check_true("closeQ[LAPACK`zgesv[{{Complex[1,1],0},{0,2}},{Complex[1,1],4}], "
               "LinearSolve[{{Complex[1,1],0},{0,2}},{Complex[1,1],4}]]");
    /* Hermitian eigenvalues are real and match Eigenvalues. */
    check_true("closeQ[First[LAPACK`zheev[{{2,Complex[0,1]},{Complex[0,-1],2}}]], "
               "Sort[N[Eigenvalues[{{2,Complex[0,1]},{Complex[0,-1],2}}]]]]");
    /* Complex QR reconstructs A. */
    check_true("Module[{q,r},{q,r}=LAPACK`zgeqrf[{{Complex[1,1],2},{3,Complex[0,4]}}];"
               "closeQ[q.r, {{Complex[1,1],2},{3,Complex[0,4]}}]]");
    check_true("closeQ[LAPACK`zgetri[{{Complex[2,0],0},{0,Complex[0,1]}}], "
               "Inverse[{{Complex[2,0],0},{0,Complex[0,1]}}]]");
}

static void test_unevaluated(void) {
    assert_eval_eq("LAPACK`dgesv[{{a,b},{c,d}},{e,f}]",
                   "LAPACK`dgesv[List[List[a, b], List[c, d]], List[e, f]]", 1);
    /* Singular system stays unevaluated. */
    assert_eval_eq("LAPACK`dgetri[{{1,2},{2,4}}]",
                   "LAPACK`dgetri[List[List[1, 2], List[2, 4]]]", 1);
}

int main(void) {
    symtab_init();
    core_init();

    if (!lapack_available()) {
        printf("SKIP: LAPACK not available in this build (USE_LAPACK off)\n");
        return 0;
    }

    run("closeQ[a_,b_]:=Max[Abs[Flatten[Normal[a]-Normal[b]]]]<1/100000000");

    printf("Running LAPACK builtin tests...\n");
    TEST(test_solve);
    TEST(test_factor);
    TEST(test_svd);
    TEST(test_eigen);
    TEST(test_norm);
    TEST(test_complex);
    TEST(test_unevaluated);
    printf("All LAPACK builtin tests passed.\n");
    return 0;
}
