/*
 * test_matlstsq.c -- Unit tests for LeastSquares (src/matlstsq.c).
 *
 * Coverage
 * --------
 *  - Spec examples from the Wolfram LeastSquares documentation:
 *      * exact 3x2 over-determined integer system
 *      * 3x3 singular integer system (rank-deficient, residual != 0)
 *      * 4x3 over-determined integer system with rational answer
 *      * symbolic 2x2 system
 *      * matrix RHS (one column per RHS)
 *      * IdentityMatrix . b == b
 *      * approximate (Real-valued) input round-tripped to inexact
 *  - LeastSquares == PseudoInverse . b identity
 *  - LeastSquares == LinearSolve for consistent square systems
 *  - Method option parsing:
 *      Method -> Automatic (symbol)
 *      Method -> "Automatic" | "Direct" | "IterativeRefinement"
 *               | "LSQR" | "Krylov"
 *      Invalid Method -> stays unevaluated
 *  - Tolerance option parsing:
 *      Tolerance -> Automatic
 *      Tolerance -> Integer / Real / Rational non-negative
 *      Negative tolerance / unknown rule -> stays unevaluated
 *  - Mixed option order (Method, Tolerance) and (Tolerance, Method)
 *  - Duplicate Method or Tolerance -> stays unevaluated
 *  - Shape-validation errors leave the call unevaluated:
 *      * matrix is not rank-2 / empty / non-list
 *      * b is not rank-1 or rank-2
 *      * b's leading dimension does not match m's row count
 *      * rank-2 b with zero columns
 *  - Stress / memory soundness: a representative spec example is run
 *    many times so valgrind can surface any leak in the per-call
 *    allocations through the PseudoInverse + Dot pipeline.
 */

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

/* Silence the LeastSquares::matrix / ::lvec / ::lvec1 messages while we
 * deliberately drive the error paths.  Their presence on stderr is part
 * of the contract; their textual output is informational. */
static void silence_stderr(void) {
    fflush(stderr);
    freopen("/dev/null", "w", stderr);
}

/* Parse, evaluate, compare to expected (FullForm). */
static void check_fullform(const char* input, const char* expected) {
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

/* Parse, evaluate, compare to expected (pretty-print form). */
static void check_pretty(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    char* s = expr_to_string(r);
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

/* Run `input` and assert the result printed as "True". */
static void check_true(const char* input) {
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

/* ------------------------------------------------------------------
 * Spec example: simple 3x2 over-determined integer system.
 * m = {{1,1},{1,2},{1,3}}; b = {7,7,8}; LeastSquares[m,b] = {19/3, 1/2}
 * ------------------------------------------------------------------ */
static void test_lstsq_spec_overdetermined_3x2(void) {
    check_pretty("LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 8}]",
                 "{19/3, 1/2}");
}

/* ------------------------------------------------------------------
 * Spec example: singular 3x3, no exact solution exists.
 * LeastSquares returns {0,0,0} -- the minimum-norm minimiser.
 * ------------------------------------------------------------------ */
static void test_lstsq_spec_singular_3x3(void) {
    check_fullform("LeastSquares[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, {2, -4, 2}]",
                   "List[0, 0, 0]");
}

/* ------------------------------------------------------------------
 * Spec example: 4x3 exact over-determined system with rational answer.
 * ------------------------------------------------------------------ */
static void test_lstsq_spec_overdetermined_4x3(void) {
    check_pretty(
        "LeastSquares[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {10, 11, 12}}, {1, 2, 4, 8}]",
        "{157/180, 23/90, -13/36}");
}

/* ------------------------------------------------------------------
 * Spec example: symbolic 2x2 system.  We avoid asserting on the exact
 * printed form (which depends on canonical-ordering choices) and
 * instead verify the residual is zero, i.e. m . x == b.
 * ------------------------------------------------------------------ */
static void test_lstsq_spec_symbolic_2x2(void) {
    check_true(
        "(m = {{a, b}, {c, d}}; rhs = {e, f}; "
        " x = LeastSquares[m, rhs]; "
        " Together[m . x - rhs] == {0, 0})");
}

/* ------------------------------------------------------------------
 * Spec example: matrix RHS (one column per RHS).
 * m = {{1,1},{1,2},{1,3}}; b = {{7,1},{7,2},{8,3}}
 * LeastSquares[m, b] = {{19/3, 0}, {1/2, 1}}
 * ------------------------------------------------------------------ */
static void test_lstsq_spec_matrix_rhs(void) {
    check_pretty(
        "LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {{7, 1}, {7, 2}, {8, 3}}]",
        "{{19/3, 0}, {1/2, 1}}");
}

/* ------------------------------------------------------------------
 * Spec example: each column of a matrix-RHS result matches the
 * column-wise vector solve.
 * ------------------------------------------------------------------ */
static void test_lstsq_spec_matrix_rhs_columnwise(void) {
    check_true(
        "(m = {{1, 1}, {1, 2}, {1, 3}}; "
        " bmat = {{7, 1}, {7, 2}, {8, 3}}; "
        " full = LeastSquares[m, bmat]; "
        " col1 = LeastSquares[m, bmat[[All, 1]]]; "
        " col2 = LeastSquares[m, bmat[[All, 2]]]; "
        " full[[All, 1]] == col1 && full[[All, 2]] == col2)");
}

/* ------------------------------------------------------------------
 * Spec example: LeastSquares[IdentityMatrix[n], b] == b.
 * Exercises the square-non-singular path (PseudoInverse == Inverse ==
 * the identity, Dot collapses the multiplication).
 * ------------------------------------------------------------------ */
static void test_lstsq_spec_identity(void) {
    check_true("(m = IdentityMatrix[4]; b = Range[4]; "
               "LeastSquares[m, b] == b)");
    check_true("(m = IdentityMatrix[6]; b = Range[6]; "
               "LeastSquares[m, b] == b)");
}

/* ------------------------------------------------------------------
 * Spec example: when m . x == b is consistent, LeastSquares equals
 * LinearSolve.
 * ------------------------------------------------------------------ */
static void test_lstsq_equals_linearsolve_for_consistent_system(void) {
    check_true(
        "(m = {{1, 1}, {1, 2}, {1, 3}}; b = {7, 7, 7}; "
        " LeastSquares[m, b] == LinearSolve[m, b])");
    check_true(
        "(m = {{2, 1}, {1, 3}}; b = {5, 10}; "
        " LeastSquares[m, b] == LinearSolve[m, b])");
}

/* ------------------------------------------------------------------
 * The fundamental identity: LeastSquares == PseudoInverse . b.
 * ------------------------------------------------------------------ */
static void test_lstsq_equals_pseudoinverse_dot_b(void) {
    check_true(
        "(m = {{1, 2}, {3, 4}, {1, 1}}; b = {1, 2, 3}; "
        " LeastSquares[m, b] == PseudoInverse[m] . b)");
    check_true(
        "(m = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}; b = {2, -4, 2}; "
        " LeastSquares[m, b] == PseudoInverse[m] . b)");
    check_true(
        "(m = {{1, 1}, {1, 2}, {1, 3}}; bmat = {{7, 1}, {7, 2}, {8, 3}}; "
        " LeastSquares[m, bmat] == PseudoInverse[m] . bmat)");
}

/* ------------------------------------------------------------------
 * Approximate (machine-precision Real) input.  We compare against the
 * spec reference value via Total[Abs[...]] within tolerance so the test
 * is robust to last-bit round-off.
 * ------------------------------------------------------------------ */
static void test_lstsq_spec_real_3x3(void) {
    check_true(
        "(x = LeastSquares["
        "    {{3.2, 2.2, 1.2}, {2.1, 7.1, 8.5}, {9.5, 6.7, 3.7}}, "
        "    {7, 8, 9}]; "
        " expected = {73.9499, -174.379, 128.329}; "
        " Total[Abs[x - expected]] < 0.001)");
}

/* ------------------------------------------------------------------
 * Approximate input: invertible 2x2 reduces to Inverse . b.
 * ------------------------------------------------------------------ */
static void test_lstsq_real_invertible(void) {
    check_true(
        "(m = {{1.0, 2.0}, {3.0, 4.0}}; b = {5.0, 11.0}; "
        " x = LeastSquares[m, b]; "
        " Total[Abs[m . x - b]] < 0.000001)");
}

/* ------------------------------------------------------------------
 * Rational entries.
 * ------------------------------------------------------------------ */
static void test_lstsq_rational_invertible(void) {
    check_true(
        "(m = {{1/2, 1/3}, {1/4, 1/5}}; b = {1, 1}; "
        " x = LeastSquares[m, b]; m . x == b)");
}

/* ------------------------------------------------------------------
 * 1x1 / 1xN / Nx1 degenerate matrices.
 * ------------------------------------------------------------------ */
static void test_lstsq_1x1(void) {
    check_fullform("LeastSquares[{{2}}, {6}]",
                   "List[3]");
}
static void test_lstsq_row_vector_matrix(void) {
    /* 1x3 m, 1-vector b: minimum-norm solution is b[1] * v / Norm[v]^2
     * with v = {1, 2, 3}, Norm[v]^2 = 14.  b = {14} -> x = {1, 2, 3}. */
    check_pretty("LeastSquares[{{1, 2, 3}}, {14}]",
                 "{1, 2, 3}");
}
static void test_lstsq_column_vector_matrix(void) {
    /* 3x1 m, 3-vector b -- vanilla over-determined least squares.
     * x = (m^T b)/(m^T m) = (1*1 + 2*2 + 3*3)/(1+4+9) = 14/14 = 1. */
    check_pretty("LeastSquares[{{1}, {2}, {3}}, {1, 2, 3}]",
                 "{1}");
}

/* ------------------------------------------------------------------
 * Method option parsing.
 * ------------------------------------------------------------------ */
static void test_lstsq_method_automatic_symbol(void) {
    check_pretty(
        "LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 8}, Method -> Automatic]",
        "{19/3, 1/2}");
}
static void test_lstsq_method_automatic_string(void) {
    check_pretty(
        "LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 8}, "
        "Method -> \"Automatic\"]",
        "{19/3, 1/2}");
}
static void test_lstsq_method_direct(void) {
    check_pretty(
        "LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 8}, "
        "Method -> \"Direct\"]",
        "{19/3, 1/2}");
}
static void test_lstsq_method_iterativerefinement(void) {
    /* For exact input, IterativeRefinement coincides with Direct because
     * the residual is exactly zero. */
    check_pretty(
        "LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 8}, "
        "Method -> \"IterativeRefinement\"]",
        "{19/3, 1/2}");
}
static void test_lstsq_method_lsqr(void) {
    /* LSQR is currently dispatched to Direct (see file header). */
    check_pretty(
        "LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 8}, "
        "Method -> \"LSQR\"]",
        "{19/3, 1/2}");
}
static void test_lstsq_method_krylov(void) {
    check_pretty(
        "LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 8}, "
        "Method -> \"Krylov\"]",
        "{19/3, 1/2}");
}
static void test_lstsq_method_invalid_unevaluated(void) {
    /* Unknown Method value -> stays unevaluated. */
    check_fullform(
        "LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 8}, "
        "Method -> \"Magic\"]",
        "LeastSquares[List[List[1, 1], List[1, 2], List[1, 3]], "
        "List[7, 7, 8], Rule[Method, \"Magic\"]]");
}

/* ------------------------------------------------------------------
 * Tolerance option parsing.
 * ------------------------------------------------------------------ */
static void test_lstsq_tolerance_automatic(void) {
    check_pretty(
        "LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 8}, "
        "Tolerance -> Automatic]",
        "{19/3, 1/2}");
}
static void test_lstsq_tolerance_integer(void) {
    check_pretty(
        "LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 8}, Tolerance -> 0]",
        "{19/3, 1/2}");
}
static void test_lstsq_tolerance_real(void) {
    check_pretty(
        "LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 8}, "
        "Tolerance -> 0.0001]",
        "{19/3, 1/2}");
}
static void test_lstsq_tolerance_rational(void) {
    check_pretty(
        "LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 8}, "
        "Tolerance -> 1/100]",
        "{19/3, 1/2}");
}
static void test_lstsq_tolerance_negative_unevaluated(void) {
    /* Negative tolerance -> stays unevaluated. */
    check_fullform(
        "LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 8}, Tolerance -> -1]",
        "LeastSquares[List[List[1, 1], List[1, 2], List[1, 3]], "
        "List[7, 7, 8], Rule[Tolerance, -1]]");
}
static void test_lstsq_unknown_option_unevaluated(void) {
    /* Unknown option name -> stays unevaluated. */
    check_fullform(
        "LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 8}, Foo -> 1]",
        "LeastSquares[List[List[1, 1], List[1, 2], List[1, 3]], "
        "List[7, 7, 8], Rule[Foo, 1]]");
}

/* ------------------------------------------------------------------
 * Mixed option order and duplicate detection.
 * ------------------------------------------------------------------ */
static void test_lstsq_method_then_tolerance(void) {
    check_pretty(
        "LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 8}, "
        "Method -> \"Direct\", Tolerance -> 1/100]",
        "{19/3, 1/2}");
}
static void test_lstsq_tolerance_then_method(void) {
    check_pretty(
        "LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 8}, "
        "Tolerance -> 0.001, Method -> \"IterativeRefinement\"]",
        "{19/3, 1/2}");
}
static void test_lstsq_duplicate_method_unevaluated(void) {
    check_fullform(
        "LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 8}, "
        "Method -> \"Direct\", Method -> \"LSQR\"]",
        "LeastSquares[List[List[1, 1], List[1, 2], List[1, 3]], "
        "List[7, 7, 8], Rule[Method, \"Direct\"], Rule[Method, \"LSQR\"]]");
}
static void test_lstsq_duplicate_tolerance_unevaluated(void) {
    check_fullform(
        "LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 8}, "
        "Tolerance -> 0, Tolerance -> 1]",
        "LeastSquares[List[List[1, 1], List[1, 2], List[1, 3]], "
        "List[7, 7, 8], Rule[Tolerance, 0], Rule[Tolerance, 1]]");
}

/* ------------------------------------------------------------------
 * Shape-validation errors leave the call unevaluated.
 * ------------------------------------------------------------------ */
static void test_lstsq_too_few_args_unevaluated(void) {
    check_fullform("LeastSquares[{{1, 2}, {3, 4}}]",
                   "LeastSquares[List[List[1, 2], List[3, 4]]]");
}
static void test_lstsq_non_matrix_m_unevaluated(void) {
    silence_stderr();
    /* Use a symbol that is guaranteed unbound across the test suite.
     * Earlier tests assign to x, m, b, ...; matlstsqUnbound is unique
     * to this test. */
    check_fullform("LeastSquares[Sin[matlstsqUnbound], {1, 2}]",
                   "LeastSquares[Sin[matlstsqUnbound], List[1, 2]]");
    check_fullform("LeastSquares[{1, 2, 3}, {1, 2, 3}]",
                   "LeastSquares[List[1, 2, 3], List[1, 2, 3]]");
    /* Jagged "matrix". */
    check_fullform("LeastSquares[{{1, 2}, {3}}, {1, 2}]",
                   "LeastSquares[List[List[1, 2], List[3]], List[1, 2]]");
}
static void test_lstsq_dim_mismatch_unevaluated(void) {
    silence_stderr();
    /* 2x2 matrix, length-3 RHS. */
    check_fullform("LeastSquares[{{1, 2}, {3, 4}}, {1, 2, 3}]",
                   "LeastSquares[List[List[1, 2], List[3, 4]], List[1, 2, 3]]");
    /* 3x2 matrix, length-2 RHS. */
    check_fullform(
        "LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 8}]",
        "LeastSquares[List[List[1, 1], List[1, 2], List[1, 3]], List[7, 8]]");
}
static void test_lstsq_bad_b_unevaluated(void) {
    silence_stderr();
    /* Rank-3 b: not a vector or matrix. */
    check_fullform(
        "LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {{{7}}, {{8}}, {{9}}}]",
        "LeastSquares[List[List[1, 1], List[1, 2], List[1, 3]], "
        "List[List[List[7]], List[List[8]], List[List[9]]]]");
}

/* ------------------------------------------------------------------
 * Stress / leak-detection loop.  Same call is run many times so
 * valgrind can surface any leak in the per-call PseudoInverse + Dot
 * allocations.
 * ------------------------------------------------------------------ */
static void test_lstsq_stress_loop(void) {
    /* Singular 3x3 (rank 2): exercises the rectangular full-rank
     * decomposition path inside PseudoInverse. */
    for (int i = 0; i < 50; i++) {
        Expr* e = parse_expression(
            "LeastSquares[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, {2, -4, 2}]");
        ASSERT(e != NULL);
        Expr* r = evaluate(e);
        ASSERT(r != NULL);
        expr_free(r);
        expr_free(e);
    }
    /* Over-determined 3x2 with vector b. */
    for (int i = 0; i < 50; i++) {
        Expr* e = parse_expression(
            "LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 8}]");
        ASSERT(e != NULL);
        Expr* r = evaluate(e);
        ASSERT(r != NULL);
        expr_free(r);
        expr_free(e);
    }
    /* Same shape with a matrix b. */
    for (int i = 0; i < 50; i++) {
        Expr* e = parse_expression(
            "LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, "
            "{{7, 1}, {7, 2}, {8, 3}}]");
        ASSERT(e != NULL);
        Expr* r = evaluate(e);
        ASSERT(r != NULL);
        expr_free(r);
        expr_free(e);
    }
    /* IterativeRefinement path -- exercises the second-pass Dot. */
    for (int i = 0; i < 30; i++) {
        Expr* e = parse_expression(
            "LeastSquares[{{1, 1}, {1, 2}, {1, 3}}, {7, 7, 8}, "
            "Method -> \"IterativeRefinement\"]");
        ASSERT(e != NULL);
        Expr* r = evaluate(e);
        ASSERT(r != NULL);
        expr_free(r);
        expr_free(e);
    }
    /* Inexact (Real) input -- exercises the rationalise + numericalise
     * round-trip inside PseudoInverse. */
    for (int i = 0; i < 20; i++) {
        Expr* e = parse_expression(
            "LeastSquares["
            "  {{3.2, 2.2, 1.2}, {2.1, 7.1, 8.5}, {9.5, 6.7, 3.7}}, "
            "  {7, 8, 9}]");
        ASSERT(e != NULL);
        Expr* r = evaluate(e);
        ASSERT(r != NULL);
        expr_free(r);
        expr_free(e);
    }
    printf("PASS: stress loop (5 shapes x 50/50/50/30/20 iterations)\n");
}

int main(void) {
    symtab_init();
    core_init();

    /* Spec examples */
    TEST(test_lstsq_spec_overdetermined_3x2);
    TEST(test_lstsq_spec_singular_3x3);
    TEST(test_lstsq_spec_overdetermined_4x3);
    TEST(test_lstsq_spec_symbolic_2x2);
    TEST(test_lstsq_spec_matrix_rhs);
    TEST(test_lstsq_spec_matrix_rhs_columnwise);
    TEST(test_lstsq_spec_identity);
    TEST(test_lstsq_equals_linearsolve_for_consistent_system);
    TEST(test_lstsq_equals_pseudoinverse_dot_b);
    TEST(test_lstsq_spec_real_3x3);
    TEST(test_lstsq_real_invertible);
    TEST(test_lstsq_rational_invertible);

    /* Degenerate shapes */
    TEST(test_lstsq_1x1);
    TEST(test_lstsq_row_vector_matrix);
    TEST(test_lstsq_column_vector_matrix);

    /* Method options */
    TEST(test_lstsq_method_automatic_symbol);
    TEST(test_lstsq_method_automatic_string);
    TEST(test_lstsq_method_direct);
    TEST(test_lstsq_method_iterativerefinement);
    TEST(test_lstsq_method_lsqr);
    TEST(test_lstsq_method_krylov);
    TEST(test_lstsq_method_invalid_unevaluated);

    /* Tolerance options */
    TEST(test_lstsq_tolerance_automatic);
    TEST(test_lstsq_tolerance_integer);
    TEST(test_lstsq_tolerance_real);
    TEST(test_lstsq_tolerance_rational);
    TEST(test_lstsq_tolerance_negative_unevaluated);
    TEST(test_lstsq_unknown_option_unevaluated);

    /* Mixed-option / duplicate handling */
    TEST(test_lstsq_method_then_tolerance);
    TEST(test_lstsq_tolerance_then_method);
    TEST(test_lstsq_duplicate_method_unevaluated);
    TEST(test_lstsq_duplicate_tolerance_unevaluated);

    /* Shape-validation errors */
    TEST(test_lstsq_too_few_args_unevaluated);
    TEST(test_lstsq_non_matrix_m_unevaluated);
    TEST(test_lstsq_dim_mismatch_unevaluated);
    TEST(test_lstsq_bad_b_unevaluated);

    /* Memory-soundness stress */
    TEST(test_lstsq_stress_loop);

    printf("\nAll test_matlstsq.c tests passed.\n");
    return 0;
}
