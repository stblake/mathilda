/*
 * test_linearsolve.c -- unit tests for LinearSolve.
 *
 * Strategy
 * --------
 *  Numerical / exact-integer cases compare against a fixed canonical
 *  FullForm string, since the result is unique up to printing.
 *
 *  Symbolic and floating-point cases instead verify the defining
 *  relation
 *       m . LinearSolve[m, b]  ==  b
 *  (under Simplify when necessary).  That avoids brittle dependence
 *  on whatever ordering Mathilda's Plus / Times canonicaliser
 *  happens to produce on a given platform.
 *
 *  Inconsistent / shape-error cases verify that LinearSolve returns
 *  unevaluated (i.e. prints back as itself, in FullForm).  The
 *  LinearSolve::nosol / ::lvec1 / ::matrix messages they print to
 *  stderr are not inspected here -- they are explicitly informational
 *  output that the REPL relies on.
 */

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

/* Suppress LinearSolve::nosol / ::lvec1 / ::matrix messages while
 * running the negative-path tests.  We still want the *behaviour*
 * (return unevaluated) to be tested, just not the noise on stderr. */
static void silence_stderr(void) { fflush(stderr); freopen("/dev/null", "w", stderr); }

static void run_test(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    char* s = expr_to_string_fullform(r);
    if (strcmp(s, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n",
               input, expected, s);
        ASSERT(0);
    } else {
        printf("PASS: %s -> %s\n", input, s);
    }
    free(s);
    expr_free(r);
    expr_free(e);
}

/* Evaluate `input` and assert the printed result is the string
 * "True".  Used for verification-by-substitution tests. */
static void run_check_true(const char* input) {
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    char* s = expr_to_string(r);
    if (strcmp(s, "True") != 0) {
        printf("FAIL (expected True): %s\n  got: %s\n", input, s);
        ASSERT(0);
    } else {
        printf("PASS: %s -> True\n", input);
    }
    free(s);
    expr_free(r);
    expr_free(e);
}

/* -----------------------------------------------------------------
 * Diagonal / trivial systems
 * ----------------------------------------------------------------- */
static void test_linearsolve_identity_and_diagonal(void) {
    /* Identity: x = b. */
    run_test("LinearSolve[{{1, 0}, {0, 1}}, {3, 5}]",
             "List[3, 5]");
    run_test("LinearSolve[IdentityMatrix[3], {1, 2, 3}]",
             "List[1, 2, 3]");

    /* Diagonal scaling. */
    run_test("LinearSolve[{{2, 0}, {0, 3}}, {6, 12}]",
             "List[3, 4]");
    run_test("LinearSolve[{{2, 0, 0}, {0, 3, 0}, {0, 0, 5}}, {4, 9, 25}]",
             "List[2, 3, 5]");

    /* 1x1 system. */
    run_test("LinearSolve[{{5}}, {15}]", "List[3]");
    run_test("LinearSolve[{{a}}, {b}]",
             "List[Times[Power[a, -1], b]]");
}

/* -----------------------------------------------------------------
 * Exact-integer / exact-rational square systems
 * ----------------------------------------------------------------- */
static void test_linearsolve_exact_square(void) {
    /* 2x2 with rational solution.  x1 + 2 x2 = 5, 3 x1 + 4 x2 = 6  =>
     * x1 = -4, x2 = 9/2. */
    run_test("LinearSolve[{{1, 2}, {3, 4}}, {5, 6}]",
             "List[-4, Rational[9, 2]]");

    /* Matrix RHS: LinearSolve[{{1,2},{3,4}}, {{5,6},{7,8}}] = {{-3,-4},{4,5}} */
    run_test("LinearSolve[{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}]",
             "List[List[-3, -4], List[4, 5]]");

    /* 3x3 exact integer matrix. */
    run_check_true(
        "CompoundExpression[Set[mm, {{2, 3, 2}, {4, 9, 2}, {7, 2, 4}}], "
        "Set[bb, {1, 2, 3}], "
        "Equal[Dot[mm, LinearSolve[mm, bb]], bb]]");

    /* 3x3 exact integer matrix, second exact-integer instance. */
    run_check_true(
        "CompoundExpression[Set[mm, {{1, 2, 3}, {4, 2, 2}, {5, 1, 7}}], "
        "Set[bb, {6, 7, 8}], "
        "Equal[Dot[mm, LinearSolve[mm, bb]], bb]]");

    /* 4x4 invertible system. */
    run_check_true(
        "CompoundExpression[Set[mm, {{1, 0, 2, 1}, {3, 1, 0, 1}, "
        "{1, 2, 1, 0}, {0, 1, 3, 2}}], "
        "Set[bb, {2, 5, 3, 7}], "
        "Equal[Dot[mm, LinearSolve[mm, bb]], bb]]");
}

/* -----------------------------------------------------------------
 * Symbolic systems
 * ----------------------------------------------------------------- */
static void test_linearsolve_symbolic(void) {
    /* 2x2 fully symbolic.  We don't pin the canonical form -- just
     * insist the answer back-substitutes. */
    run_check_true(
        "CompoundExpression[Set[mm, {{r, s}, {t, u}}], "
        "Set[bb, {y, z}], "
        "Equal[Together[Dot[mm, LinearSolve[mm, bb]] - bb], {0, 0}]]");

    /* 3x3 symbolic with constant RHS. */
    run_check_true(
        "CompoundExpression[Set[mm, {{a, b, c}, {d, e, f}, {g, h, i}}], "
        "Set[bb, {3, 2, 1}], "
        "Equal[Together[Dot[mm, LinearSolve[mm, bb]] - bb], {0, 0, 0}]]");
}

/* -----------------------------------------------------------------
 * Rectangular systems
 * ----------------------------------------------------------------- */
static void test_linearsolve_rectangular(void) {
    /* Overdetermined but consistent: spec example.
     *   m is 4x2, b is length 4.  Answer is {-1, 2}. */
    run_test("LinearSolve[{{1, 5}, {2, 6}, {3, 7}, {4, 8}}, {9, 10, 11, 12}]",
             "List[-1, 2]");

    /* Underdetermined (rank 2 in 3 unknowns).  Mathematica
     * convention: pivots take the b value, free vars are 0. */
    run_test("LinearSolve[{{1, 2, 3}, {4, 5, 6}}, {6, 15}]",
             "List[0, 3, 0]");

    /* 3x4 underdetermined exact-rational case.  Spec:
     * LinearSolve[{{1,2,3,4},{5,6,7,8},{9,10,11,12}}, {13,14,15}]
     *   == {-25/2, 51/4, 0, 0}. */
    run_test("LinearSolve[{{1, 2, 3, 4}, {5, 6, 7, 8}, "
             "{9, 10, 11, 12}}, {13, 14, 15}]",
             "List[Rational[-25, 2], Rational[51, 4], 0, 0]");

    /* Verify the rectangular answer back-substitutes. */
    run_check_true(
        "CompoundExpression[Set[mm, {{1, 5}, {2, 6}, {3, 7}, {4, 8}}], "
        "Set[bb, {9, 10, 11, 12}], "
        "Equal[Dot[mm, LinearSolve[mm, bb]], bb]]");

    /* Matrix-RHS rectangular case: 3x3 m, 3x2 b. */
    run_test("LinearSolve[{{1, 1, 1}, {1, 2, 3}, {1, 4, 9}}, "
             "{{1, 2}, {3, 4}, {5, 6}}]",
             "List[List[-2, -1], List[4, 4], List[-1, -1]]");
}

/* -----------------------------------------------------------------
 * Complex / non-integer numeric systems
 * ----------------------------------------------------------------- */
static void test_linearsolve_complex(void) {
    /* Spec example: LinearSolve over Gaussian integers. */
    run_test("LinearSolve[{{1 + I, 2, 3 - 2 I}, {0, 4, 5 I}, "
             "{1 + I, 6, 3 + 3 I}}, {1, 2, 3}]",
             "List[0, Rational[1, 2], 0]");

    /* Verify a 2x2 complex result by back-substitution. */
    run_check_true(
        "CompoundExpression[Set[mm, {{1 + I, 2}, {3, 4 - I}}], "
        "Set[bb, {5, 6 + I}], "
        "Equal[Dot[mm, LinearSolve[mm, bb]], bb]]");
}

/* -----------------------------------------------------------------
 * Floating-point input
 * ----------------------------------------------------------------- */
static void test_linearsolve_float(void) {
    /* Symmetric verification: m . LinearSolve[m, b] should equal b
     * componentwise after a Norm.  We use a strict equality test
     * with carefully chosen exact-representable floats so the
     * round-trip is exact. */
    run_check_true(
        "CompoundExpression[Set[mm, {{1.0, 0.0}, {0.0, 4.0}}], "
        "Set[bb, {2.0, 8.0}], "
        "Equal[Dot[mm, LinearSolve[mm, bb]], bb]]");
}

/* -----------------------------------------------------------------
 * Singular / inconsistent / shape-error cases
 * ----------------------------------------------------------------- */
static void test_linearsolve_no_solution(void) {
    silence_stderr();

    /* Singular 3x3 with inconsistent b: LinearSolve::nosol. */
    run_test("LinearSolve[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, {1, 1, 2}]",
             "LinearSolve[List[List[1, 2, 3], List[4, 5, 6], List[7, 8, 9]], "
             "List[1, 1, 2]]");

    /* Same singular matrix but consistent b -- should give a
     * particular solution {-1, 1, 0}. */
    run_test("LinearSolve[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, {1, 1, 1}]",
             "List[-1, 1, 0]");

    /* Identity matrix on the RHS with a singular matrix -> no
     * solution (spec example). */
    run_test("LinearSolve[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, "
             "{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}]",
             "LinearSolve[List[List[1, 2, 3], List[4, 5, 6], List[7, 8, 9]], "
             "List[List[1, 0, 0], List[0, 1, 0], List[0, 0, 1]]]");

    /* Shape mismatch: vector b is the wrong length. */
    run_test("LinearSolve[{{1, 2}, {3, 4}}, {5, 6, 7}]",
             "LinearSolve[List[List[1, 2], List[3, 4]], List[5, 6, 7]]");

    /* Argument m is not a matrix at all (scalar). */
    run_test("LinearSolve[5, {1}]", "LinearSolve[5, List[1]]");
}

/* -----------------------------------------------------------------
 * Higher-rank tensor inputs.
 *
 * Mathematica's LinearSolve interprets a rank-N matrix `m` with
 * dimensions {d_1, ..., d_{N-1}, n} as a (d_1*...*d_{N-1}) x n
 * linear system, and a rank-(N-1+p) `b` with dimensions
 * {d_1, ..., d_{N-1}, e_1, ..., e_p} produces a result of shape
 * {n, e_1, ..., e_p}.  The two natural cross-checks are:
 *   - Dot[m, x] == b should hold when the system is consistent;
 *   - Dimensions[x] should equal {n, e_1, ..., e_p}.
 * ----------------------------------------------------------------- */
static void test_linearsolve_array_inputs(void) {
    /* Rank-2 m, rank-3 b.  m = {{1,0},{0,2}}, b is 2x2x3.
     * Expected solution shape: {2, 2, 3}. */
    run_check_true(
        "CompoundExpression["
        "Set[mm, {{1, 0}, {0, 2}}], "
        "Set[bb, {{{1, 2, 3}, {4, 5, 6}}, {{2, 4, 6}, {8, 10, 12}}}], "
        "Set[xx, LinearSolve[mm, bb]], "
        "And[Equal[Dimensions[xx], {2, 2, 3}], "
        "Equal[Dot[mm, xx], bb]]]");

    /* Rank-3 m (square: leading 2x2 -> 4 equations, last dim 2 unknowns
     * but here we use {{{1,0},{0,2}},{{2,1},{1,2}}}), rank-2 b matching
     * leading dims exactly (so trail_rank == 0 -> flat vector result). */
    run_test(
        "LinearSolve["
        "{{{1, 0}, {0, 2}}, {{2, 1}, {1, 2}}}, "
        "{{3, 10}, {11, 13}}]",
        "List[3, 5]");

    /* Rank-3 m with rank-3 b: trail_rank == 1, result has shape
     * {n, k}. */
    run_test(
        "LinearSolve["
        "{{{1, 0}, {0, 2}}, {{2, 1}, {1, 2}}}, "
        "{{{3, 6, 9}, {10, 20, 30}}, {{11, 22, 33}, {13, 26, 39}}}]",
        "List[List[3, 6, 9], List[5, 10, 15]]");

    /* The bug-report case: rank-3 square m and rank-4 b, verified by
     * back-substitution.  m has dims {2, 3, 6}, b has dims
     * {2, 3, 4, 5}, result must have dims {6, 4, 5}. */
    run_check_true(
        "CompoundExpression["
        "Set[mm, {{{2, -1, 0, -3, 0, 2}, {1, 0, -2, 2, -1, 2}, "
        "         {-3, 3, 0, -2, 2, -2}}, "
        "        {{3, 1, 1, 1, 3, -3}, {-1, 0, 2, 3, 0, -1}, "
        "         {1, 0, -1, 1, 3, 3}}}], "
        "Set[bb, Table[i + j + k + l, "
        "  {i, 2}, {j, 3}, {k, 4}, {l, 5}]], "
        "Set[xx, LinearSolve[mm, bb]], "
        "And[Equal[Dimensions[xx], {6, 4, 5}], "
        "Equal[Dot[mm, xx], bb]]]");

    /* Leading-dim mismatch must error: m's leading dims {2, 2} don't
     * match b's first dim {2, 3}. */
    silence_stderr();
    run_test(
        "LinearSolve["
        "{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}}, "
        "{{1, 2, 3}, {4, 5, 6}}]",
        "LinearSolve[List[List[List[1, 2], List[3, 4]], "
        "List[List[5, 6], List[7, 8]]], "
        "List[List[1, 2, 3], List[4, 5, 6]]]");
}

/* -----------------------------------------------------------------
 * Inverse via LinearSolve[m, IdentityMatrix[n]]
 * ----------------------------------------------------------------- */
static void test_linearsolve_inverse_check(void) {
    /* For an invertible m, LinearSolve[m, IdentityMatrix[n]] should
     * yield the same result as Inverse[m]. */
    run_check_true(
        "CompoundExpression[Set[mm, {{1, 2}, {3, 4}}], "
        "Equal[LinearSolve[mm, IdentityMatrix[2]], Inverse[mm]]]");

    run_check_true(
        "CompoundExpression[Set[mm, {{1, 4, 2, -9}, {4, 12, 2, 5}, "
        "{6, 7, -11, 9}, {5, 15, 10, 12}}], "
        "Equal[LinearSolve[mm, IdentityMatrix[4]], Inverse[mm]]]");
}

int main(void) {
    /* test_utils.h sets alarm(60); under valgrind that's tight, so
     * extend it.  Native runs need no headroom -- the suite runs in
     * well under a second. */
    alarm(600);

    symtab_init();
    core_init();

    printf("Running LinearSolve tests...\n");
    TEST(test_linearsolve_identity_and_diagonal);
    TEST(test_linearsolve_exact_square);
    TEST(test_linearsolve_symbolic);
    TEST(test_linearsolve_rectangular);
    TEST(test_linearsolve_complex);
    TEST(test_linearsolve_float);
    TEST(test_linearsolve_no_solution);
    TEST(test_linearsolve_array_inputs);
    TEST(test_linearsolve_inverse_check);

    printf("All LinearSolve tests passed!\n");
    symtab_clear();
    return 0;
}
