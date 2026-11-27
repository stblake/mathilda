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
static void silence_stderr(void) { fflush(stderr); test_mute_messages(); }

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

/* A small machine-real matrix is packed to reach the buffer path. When that
 * path declined (non-square Det, singular Inverse or LinearSolve) it delisted
 * and re-evaluated the same builtin, which packed it again -- to
 * $RecursionLimit, with rows of Hold[List][...] in the result. Each now emits
 * its one message and stays unevaluated, or answers. */
static void test_linalg_machine_decline_no_recursion(void) {
    silence_stderr();
    run_check_true("Quiet[Det[{{1., 2., 3.}, {4., 5., 6.}}]] === Det[{{1., 2., 3.}, {4., 5., 6.}}]");
    run_check_true("FreeQ[Quiet[Det[{{1., 2., 3.}, {4., 5., 6.}}]], Hold]");
    run_check_true("Check[Det[{{1., 2., 3.}, {4., 5., 6.}}], failed] === failed");
    run_check_true("FreeQ[Quiet[Det[NDArray[{{1., 2., 3.}, {4., 5., 6.}}]]], Hold]");
    run_check_true("Quiet[Inverse[{{1., 2.}, {2., 4.}}]] === Inverse[{{1., 2.}, {2., 4.}}]");
    run_check_true("Check[Inverse[{{1., 2.}, {2., 4.}}], failed, Inverse::sing] === failed");
    run_check_true("FreeQ[Quiet[Inverse[NDArray[{{1., 2., 3.}, {4., 5., 6.}}]]], Hold]");
    run_check_true("Check[LinearSolve[{{1., 2.}, {2., 4.}}, {1., 3.}], failed, LinearSolve::nosol] === failed");
    /* The nonsingular machine cases keep the fast path. */
    run_check_true("Det[{{1., 2.}, {3., 4.}}] == -2.");
    run_check_true("LinearSolve[{{1., 2.}, {3., 4.}}, {1., 1.}] == {-1., 1.}");
}

/* Values the solver invents -- a free variable's 0, a pivot normalised to 1 --
 * take the input's exactness, so a machine-real input answers in Reals only. */
static void test_linalg_singular_machine_heads(void) {
    /* singular machine system: the minimum-norm solution, as in the reference
     * ({0.6, 1.2}; v0.6196) */
    run_check_true("Max[Abs[LinearSolve[{{1., 2.}, {2., 4.}}, {3., 6.}] - {0.6, 1.2}]] < 10^-14");
    run_check_true("Union[Head /@ LinearSolve[{{1., 2., 3.}, {2., 4., 6.}, {1., 1., 1.}}, {1., 2., 3.}]] === {Real}");
    run_check_true("Union[Head /@ Flatten[RowReduce[{{1., 2., 1.}, {2., 4., 2.}}]]] === {Real}");
    run_check_true("RowReduce[{{1., 2., 1.}, {2., 4., 2.}}] == {{1., 2., 1.}, {0., 0., 0.}}");
    /* Exact input stays exact. */
    run_check_true("LinearSolve[{{1, 2}, {2, 4}}, {3, 6}] === {3, 0}");
    run_check_true("RowReduce[{{1, 2, 1}, {2, 4, 2}}] === {{1, 2, 1}, {0, 0, 0}}");
    run_check_true("RowReduce[{{x, 2 x}, {1, 3}}] === {{1, 0}, {0, 1}}");
}

/* Algebraic-number entries (I Sqrt[3], Sqrt[2], Root objects): elimination
 * over the number field their atoms generate (flint_qqbar_alg_mat_*).  The
 * symbolic fraction-free path treated each atom as an opaque indeterminate,
 * and this 12x12 system over Q(I Sqrt[3]) did not finish in 40 s (Inverse
 * likewise; RowReduce took 2 s); each is now milliseconds.  Every answer is
 * checked exactly against its defining identity. */
#define ALG12_A \
    "{{0, (30091839012864*I)*Sqrt[3], 0, 1, 0, 0, 0, 0, 0, 0, 0, 0}, " \
    "{0, 0, (60183678025728*I)*Sqrt[3], 0, 1, 0, 0, 0, 0, 0, 0, 0}, " \
    "{(-90275517038592*I)*Sqrt[3], 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0}, " \
    "{0, 70214291030016 + (6687075336192*I)*Sqrt[3], 0, 1, 0, 0, 1, 0, 0, 0, 0, 0}, " \
    "{0, 0, 140428582060032 + (103649667710976*I)*Sqrt[3], 0, 1, 0, 0, 1, 0, 0, 0, 0}, " \
    "{-210642873090048 - (200612260085760*I)*Sqrt[3], 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0}, " \
    "{0, -70214291030016 - (113680280715264*I)*Sqrt[3], 0, 0, 0, 0, 1, 0, 0, 1, 0, 0}, " \
    "{0, 0, 70214291030016 - (26748301344768*I)*Sqrt[3], 0, 0, 0, 0, 1, 0, 0, 1, 0}, " \
    "{-210642873090048 - (60183678025728*I)*Sqrt[3], 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1}, " \
    "{0, -140428582060032 - (40122452017152*I)*Sqrt[3], 0, 0, 0, 0, 0, 0, 0, 1, 0, 0}, " \
    "{0, 0, -70214291030016 - (20061226008576*I)*Sqrt[3], 0, 0, 0, 0, 0, 0, 0, 1, 0}, " \
    "{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1}}"
#define ALG12_B \
    "{0, 90275517038592 - (30091839012864*I)*Sqrt[3], 0, 0, " \
    "85353086582784 - (122132000931840*I)*Sqrt[3], 0, 0, " \
    "-75043845439488 - (21918746935296*I)*Sqrt[3], 0, 0, " \
    "5108182548480 + (45044882472960*I)*Sqrt[3], 0}"

static void test_linearsolve_algebraic_field(void) {
    /* timing gates: 5 s each, against > 40 s (LinearSolve, Inverse) before */
    run_check_true("TimeConstrained[Module[{a = " ALG12_A ", b = " ALG12_B ", x},"
                   " x = LinearSolve[a, b]; Expand[a . x - b] === ConstantArray[0, 12]], 5, False]");
    run_check_true("TimeConstrained[Module[{a = " ALG12_A ", ia},"
                   " ia = Inverse[a]; Expand[a . ia] === IdentityMatrix[12]], 5, False]");
    run_check_true("TimeConstrained[RowReduce[" ALG12_A "] === IdentityMatrix[12], 5, False]");
    run_check_true("TimeConstrained[Module[{a = " ALG12_A ", d}, d = Det[a];"
                   " Abs[d - Det[N[a, 60]]] < 10^-30 Abs[d]], 5, False]");
    /* every Method answers the same */
    run_check_true("TimeConstrained[Module[{a = " ALG12_A ", b = " ALG12_B "},"
                   " LinearSolve[a, b, Method -> \"OneStepRowReduction\"] === LinearSolve[a, b]], 5, False]");

    /* rectangular consistent / inconsistent, singular, mixed atoms, Root */
    run_check_true("LinearSolve[{{1, Sqrt[2]}, {Sqrt[2], 2}, {I, I Sqrt[2]}}, {1, Sqrt[2], I}] === {1, 0}");
    silence_stderr();
    run_check_true("Check[LinearSolve[{{1, Sqrt[2]}, {Sqrt[2], 2}, {I, I Sqrt[2]}}, {1, 2, 3}], "
                   "failed, LinearSolve::nosol] === failed");
    run_check_true("Check[Inverse[{{1, Sqrt[2]}, {Sqrt[2], 2}}], failed, Inverse::sing] === failed");
    run_check_true("RowReduce[{{1, Sqrt[2]}, {Sqrt[2], 2}, {I, I Sqrt[2]}}] === {{1, Sqrt[2]}, {0, 0}, {0, 0}}");
    run_check_true("Det[{{1, Sqrt[2]}, {Sqrt[3], 1}}] === 1 - Sqrt[6]");
    run_check_true("Module[{a = {{1, Sqrt[2], Sqrt[3]}, {Sqrt[6], 1, 0}, {0, 1, 1}}, x},"
                   " x = LinearSolve[a, {1, 0, 0}]; Expand[a . x] === {1, 0, 0}]");
    run_check_true("Module[{r = Root[#^3 - # - 1 &, 1], a, x}, a = {{r, 1}, {1, 2}};"
                   " x = LinearSolve[a, {1, 0}]; RootReduce[a . x] === {1, 0}]");
}

/* A square system over Q(k, x): the coordinates of the GeneralisedGoursat
 * cube-root descent over a parameter (3 x 3, entries of degree 12).  The
 * division-free loop multiplied and exactly divided unreduced rational
 * functions through the evaluator: 26 s, and an answer whose check by Together
 * did not finish in a minute.  The reduced row echelon form of [m | b] on
 * fmpz_mpoly (RowReduce's path) takes milliseconds. */
#define PARAM3_A \
    "{{(2 + 14*k - 30*k^2 + 16*k^3 + 16*k^4 - 30*k^5 + 14*k^6 + 2*k^7 - 2*x + 4*k*x - 2*k^2*x +" \
    " 2*k^3*x - 2*k^4*x + 4*k^5*x - 2*k^6*x)/(4 - 12*k - 3*k^2 + 26*k^3 - 3*k^4 - 12*k^5 + 4*k^" \
    "6 + 4*k*x - 6*k^2*x - 6*k^3*x + 4*k^4*x + k^2*x^2), (-14 + 52*k - 212*k^2 + 384*k^3 - 214*" \
    "k^4 - 214*k^5 + 384*k^6 - 212*k^7 + 52*k^8 - 14*k^9 + 28*k*x - 68*k^2*x + 52*k^3*x - 28*k^" \
    "4*x + 52*k^5*x - 68*k^6*x + 28*k^7*x - 2*x^2 + 4*k*x^2 - 2*k^2*x^2 - 2*k^5*x^2 + 4*k^6*x^2" \
    " - 2*k^7*x^2)/(8 - 36*k + 18*k^2 + 105*k^3 - 99*k^4 - 99*k^5 + 105*k^6 + 18*k^7 - 36*k^8 +" \
    " 8*k^9 + 12*k*x - 36*k^2*x - 9*k^3*x + 78*k^4*x - 9*k^5*x - 36*k^6*x + 12*k^7*x + 6*k^2*x^" \
    "2 - 9*k^3*x^2 - 9*k^4*x^2 + 6*k^5*x^2 + k^3*x^3), (-15 - 168*k + 939*k^2 - 2034*k^3 + 2358" \
    "*k^4 - 1080*k^5 - 1080*k^6 + 2358*k^7 - 2034*k^8 + 939*k^9 - 168*k^10 - 15*k^11 + 23*x - 9" \
    "9*k*x + 206*k^2*x - 350*k^3*x + 500*k^4*x - 560*k^5*x + 500*k^6*x - 350*k^7*x + 206*k^8*x " \
    "- 99*k^9*x + 23*k^10*x - x^2 - 10*k*x^2 + 42*k^2*x^2 - 64*k^3*x^2 + 32*k^4*x^2 + 32*k^5*x^" \
    "2 - 64*k^6*x^2 + 42*k^7*x^2 - 10*k^8*x^2 - k^9*x^2 + x^3 - 3*k*x^3 + 3*k^2*x^3 - k^3*x^3 -" \
    " k^4*x^3 - k^5*x^3 + 3*k^6*x^3 - 3*k^7*x^3 + k^8*x^3)/(8 - 36*k + 18*k^2 + 105*k^3 - 99*k^" \
    "4 - 99*k^5 + 105*k^6 + 18*k^7 - 36*k^8 + 8*k^9 + 12*k*x - 36*k^2*x - 9*k^3*x + 78*k^4*x - " \
    "9*k^5*x - 36*k^6*x + 12*k^7*x + 6*k^2*x^2 - 9*k^3*x^2 - 9*k^4*x^2 + 6*k^5*x^2 + k^3*x^3)}," \
    " {(-10 + 2*k + 8*k^2 - 8*k^3 + 8*k^4 + 2*k^5 - 10*k^6 + 2*x - 4*k*x - 4*k^4*x + 2*k^5*x)/(" \
    "4 - 12*k - 3*k^2 + 26*k^3 - 3*k^4 - 12*k^5 + 4*k^6 + 4*k*x - 6*k^2*x - 6*k^3*x + 4*k^4*x +" \
    " k^2*x^2), (14 + 36*k - 64*k^2 - 32*k^3 + 108*k^4 - 32*k^5 - 64*k^6 + 36*k^7 + 14*k^8 - 8*" \
    "x - 4*k*x + 40*k^2*x - 24*k^3*x - 24*k^4*x + 40*k^5*x - 4*k^6*x - 8*k^7*x + 2*x^2 - 4*k*x^" \
    "2 + 2*k^2*x^2 + 2*k^4*x^2 - 4*k^5*x^2 + 2*k^6*x^2)/(8 - 36*k + 18*k^2 + 105*k^3 - 99*k^4 -" \
    " 99*k^5 + 105*k^6 + 18*k^7 - 36*k^8 + 8*k^9 + 12*k*x - 36*k^2*x - 9*k^3*x + 78*k^4*x - 9*k" \
    "^5*x - 36*k^6*x + 12*k^7*x + 6*k^2*x^2 - 9*k^3*x^2 - 9*k^4*x^2 + 6*k^5*x^2 + k^3*x^3), (87" \
    " - 144*k - 75*k^2 + 219*k^3 - 12*k^4 - 150*k^5 - 12*k^6 + 219*k^7 - 75*k^8 - 144*k^9 + 87*" \
    "k^10 - 27*x + 87*k*x - 54*k^2*x - 39*k^3*x + 33*k^4*x + 33*k^5*x - 39*k^6*x - 54*k^7*x + 8" \
    "7*k^8*x - 27*k^9*x + 5*x^2 - 6*k*x^2 - 10*k^2*x^2 + 43*k^3*x^2 - 60*k^4*x^2 + 43*k^5*x^2 -" \
    " 10*k^6*x^2 - 6*k^7*x^2 + 5*k^8*x^2 - x^3 + 3*k*x^3 - 3*k^2*x^3 + 2*k^3*x^3 + 2*k^4*x^3 - " \
    "3*k^5*x^3 + 3*k^6*x^3 - k^7*x^3)/(8 - 36*k + 18*k^2 + 105*k^3 - 99*k^4 - 99*k^5 + 105*k^6 " \
    "+ 18*k^7 - 36*k^8 + 8*k^9 + 12*k*x - 36*k^2*x - 9*k^3*x + 78*k^4*x - 9*k^5*x - 36*k^6*x + " \
    "12*k^7*x + 6*k^2*x^2 - 9*k^3*x^2 - 9*k^4*x^2 + 6*k^5*x^2 + k^3*x^3)}, {(4 - 2*k^2 - 2*k^3 " \
    "+ 4*k^5 + 2*k^2*x)/(4 - 12*k - 3*k^2 + 26*k^3 - 3*k^4 - 12*k^5 + 4*k^6 + 4*k*x - 6*k^2*x -" \
    " 6*k^3*x + 4*k^4*x + k^2*x^2), (-4 - 28*k + 60*k^2 - 32*k^3 - 32*k^4 + 60*k^5 - 28*k^6 - 4" \
    "*k^7 + 4*x - 8*k*x + 4*k^2*x - 4*k^3*x + 4*k^4*x - 8*k^5*x + 4*k^6*x)/(8 - 36*k + 18*k^2 +" \
    " 105*k^3 - 99*k^4 - 99*k^5 + 105*k^6 + 18*k^7 - 36*k^8 + 8*k^9 + 12*k*x - 36*k^2*x - 9*k^3" \
    "*x + 78*k^4*x - 9*k^5*x - 36*k^6*x + 12*k^7*x + 6*k^2*x^2 - 9*k^3*x^2 - 9*k^4*x^2 + 6*k^5*" \
    "x^2 + k^3*x^3), (-42 + 114*k - 174*k^2 + 243*k^3 - 141*k^4 - 141*k^5 + 243*k^6 - 174*k^7 +" \
    " 114*k^8 - 42*k^9 + 4*x - 8*k*x - 32*k^2*x + 89*k^3*x - 106*k^4*x + 89*k^5*x - 32*k^6*x - " \
    "8*k^7*x + 4*k^8*x - 2*x^2 + 6*k*x^2 - 12*k^2*x^2 + 7*k^3*x^2 + 7*k^4*x^2 - 12*k^5*x^2 + 6*" \
    "k^6*x^2 - 2*k^7*x^2 - k^3*x^3)/(8 - 36*k + 18*k^2 + 105*k^3 - 99*k^4 - 99*k^5 + 105*k^6 + " \
    "18*k^7 - 36*k^8 + 8*k^9 + 12*k*x - 36*k^2*x - 9*k^3*x + 78*k^4*x - 9*k^5*x - 36*k^6*x + 12" \
    "*k^7*x + 6*k^2*x^2 - 9*k^3*x^2 - 9*k^4*x^2 + 6*k^5*x^2 + k^3*x^3)}}"
#define PARAM3_B \
    "{(66 - 48*k - 252*k^2 + 204*k^3 + 300*k^4 - 270*k^5 - 270*k^8 + 300*k^9 + 204*k^10 - 252*k" \
    "^11 - 48*k^12 + 66*k^13 - 38*x + 146*k*x - 358*k^2*x + 1450*k^3*x - 4170*k^4*x + 7152*k^5*" \
    "x - 8364*k^6*x + 7152*k^7*x - 4170*k^8*x + 1450*k^9*x - 358*k^10*x + 146*k^11*x - 38*k^12*" \
    "x + 6*x^2 + 4*k*x^2 - 130*k^2*x^2 + 380*k^3*x^2 - 450*k^4*x^2 + 190*k^5*x^2 + 190*k^6*x^2 " \
    "- 450*k^7*x^2 + 380*k^8*x^2 - 130*k^9*x^2 + 4*k^10*x^2 + 6*k^11*x^2 - 2*x^3 + 10*k*x^3 - 2" \
    "0*k^2*x^3 + 20*k^3*x^3 - 10*k^4*x^3 + 8*k^5*x^3 - 10*k^6*x^3 + 20*k^7*x^3 - 20*k^8*x^3 + 1" \
    "0*k^9*x^3 - 2*k^10*x^3)/(16 - 96*k + 120*k^2 + 280*k^3 - 639*k^4 - 180*k^5 + 1014*k^6 - 18" \
    "0*k^7 - 639*k^8 + 280*k^9 + 120*k^10 - 96*k^11 + 16*k^12 + 32*k*x - 144*k^2*x + 72*k^3*x +" \
    " 420*k^4*x - 396*k^5*x - 396*k^6*x + 420*k^7*x + 72*k^8*x - 144*k^9*x + 32*k^10*x + 24*k^2" \
    "*x^2 - 72*k^3*x^2 - 18*k^4*x^2 + 156*k^5*x^2 - 18*k^6*x^2 - 72*k^7*x^2 + 24*k^8*x^2 + 8*k^" \
    "3*x^3 - 12*k^4*x^3 - 12*k^5*x^3 + 8*k^6*x^3 + k^4*x^4), (-162 + 216*k + 432*k^2 - 648*k^3 " \
    "- 270*k^4 + 432*k^5 + 432*k^7 - 270*k^8 - 648*k^9 + 432*k^10 + 216*k^11 - 162*k^12 + 62*x " \
    "- 178*k*x - 154*k^2*x + 726*k^3*x - 480*k^4*x + 24*k^5*x + 24*k^6*x - 480*k^7*x + 726*k^8*" \
    "x - 154*k^9*x - 178*k^10*x + 62*k^11*x - 14*x^2 + 44*k*x^2 + 6*k^2*x^2 - 248*k^3*x^2 + 440" \
    "*k^4*x^2 - 456*k^5*x^2 + 440*k^6*x^2 - 248*k^7*x^2 + 6*k^8*x^2 + 44*k^9*x^2 - 14*k^10*x^2 " \
    "+ 2*x^3 - 10*k*x^3 + 20*k^2*x^3 - 20*k^3*x^3 + 4*k^4*x^3 + 4*k^5*x^3 - 20*k^6*x^3 + 20*k^7" \
    "*x^3 - 10*k^8*x^3 + 2*k^9*x^3)/(16 - 96*k + 120*k^2 + 280*k^3 - 639*k^4 - 180*k^5 + 1014*k" \
    "^6 - 180*k^7 - 639*k^8 + 280*k^9 + 120*k^10 - 96*k^11 + 16*k^12 + 32*k*x - 144*k^2*x + 72*" \
    "k^3*x + 420*k^4*x - 396*k^5*x - 396*k^6*x + 420*k^7*x + 72*k^8*x - 144*k^9*x + 32*k^10*x +" \
    " 24*k^2*x^2 - 72*k^3*x^2 - 18*k^4*x^2 + 156*k^5*x^2 - 18*k^6*x^2 - 72*k^7*x^2 + 24*k^8*x^2" \
    " + 8*k^3*x^3 - 12*k^4*x^3 - 12*k^5*x^3 + 8*k^6*x^3 + k^4*x^4), (60 - 36*k - 276*k^2 + 204*" \
    "k^3 + 504*k^4 - 456*k^5 - 456*k^6 + 504*k^7 + 204*k^8 - 276*k^9 - 36*k^10 + 60*k^11 - 16*x" \
    " + 32*k*x + 208*k^2*x - 672*k^3*x + 780*k^4*x - 664*k^5*x + 780*k^6*x - 672*k^7*x + 208*k^" \
    "8*x + 32*k^9*x - 16*k^10*x + 4*x^2 - 20*k*x^2 + 40*k^2*x^2 - 8*k^3*x^2 - 16*k^4*x^2 - 16*k" \
    "^5*x^2 - 8*k^6*x^2 + 40*k^7*x^2 - 20*k^8*x^2 + 4*k^9*x^2 + 4*k^4*x^3)/(16 - 96*k + 120*k^2" \
    " + 280*k^3 - 639*k^4 - 180*k^5 + 1014*k^6 - 180*k^7 - 639*k^8 + 280*k^9 + 120*k^10 - 96*k^" \
    "11 + 16*k^12 + 32*k*x - 144*k^2*x + 72*k^3*x + 420*k^4*x - 396*k^5*x - 396*k^6*x + 420*k^7" \
    "*x + 72*k^8*x - 144*k^9*x + 32*k^10*x + 24*k^2*x^2 - 72*k^3*x^2 - 18*k^4*x^2 + 156*k^5*x^2" \
    " - 18*k^6*x^2 - 72*k^7*x^2 + 24*k^8*x^2 + 8*k^3*x^3 - 12*k^4*x^3 - 12*k^5*x^3 + 8*k^6*x^3 " \
    "+ k^4*x^4)}"

static void test_linearsolve_rational_function_field(void) {
    run_check_true("TimeConstrained[Module[{a = " PARAM3_A ", b = " PARAM3_B ", s},"
                   " s = LinearSolve[a, b]; Together[a . s - b] === {0, 0, 0}], 5, False]");
    /* rectangular consistent, free variables 0, and an inconsistent system */
    run_check_true("Module[{a = {{1/(k + 1), k}, {2/(k + 1), 2 k}, {0, 0}}, s},"
                   " s = LinearSolve[a, {1/k, 2/k, 0}]; Together[a . s - {1/k, 2/k, 0}] === {0, 0, 0} && s[[2]] === 0]");
    run_check_true("Module[{a = {{k, 1/(1 - k)}, {1, k^2}}, s}, s = LinearSolve[a, {1, x}];"
                   " Together[a . s - {1, x}] === {0, 0}]");
    silence_stderr();
    run_check_true("Check[LinearSolve[{{k, 1}, {2 k, 2}}, {1, 3}], failed, LinearSolve::nosol] === failed");
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
    TEST(test_linalg_machine_decline_no_recursion);
    TEST(test_linalg_singular_machine_heads);
    TEST(test_linearsolve_algebraic_field);
    TEST(test_linearsolve_rational_function_field);

    printf("All LinearSolve tests passed!\n");
    symtab_clear();
    return 0;
}
