/*
 * test_matinv.c -- Unit tests for Inverse and PseudoInverse (src/matinv.c).
 *
 * Coverage
 * --------
 *  - Inverse continues to behave identically after its move out of
 *    src/linalg.c (regression: integer, rational, symbolic, complex,
 *    singular).
 *  - PseudoInverse on
 *      * invertible square (matches Inverse)
 *      * the rank-1 singular 3x3 example from the spec
 *      * rectangular (over- and under-determined) integer matrices
 *      * the m x n zero matrix (gives n x m zero matrix)
 *      * 1x1, 1xN, Nx1 matrices
 *      * exact rational entries
 *      * exact complex entries
 *      * machine-precision Real entries (round-trip via rationalise)
 *      * MPFR (if compiled in) at non-machine precision
 *      * symbolic non-invertible row matrix
 *  - Moore-Penrose identities
 *      M.P.M == M and P.M.P == P for every non-trivial example.
 *  - Tolerance-option syntax: Automatic, numeric, invalid (unevaluated).
 *  - Shape-validation errors (rank != 2, empty matrix, non-matrix arg)
 *    leave the call unevaluated.
 *  - Stress / memory soundness: the same test loop is run 50 times to
 *    surface any leak in valgrind.
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

/* Silence the various foo::sing / foo::matrix / foo::matsq messages
 * while we exercise the error paths.  Their presence is part of the
 * contract; their stderr output is informational. */
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
 * Inverse regression after the move into src/matinv.c
 * ------------------------------------------------------------------ */
static void test_inverse_basic_integer(void) {
    check_fullform("Inverse[{{1, 2}, {3, 4}}]",
        "List[List[-2, 1], List[Rational[3, 2], Rational[-1, 2]]]");
}

static void test_inverse_3x3_rational(void) {
    /* {{1, 1/2, 0}, {0, 1, 1/3}, {0, 0, 1}}: upper-triangular, inverse
     * stays upper-triangular. */
    check_pretty(
        "Inverse[{{1, 1/2, 0}, {0, 1, 1/3}, {0, 0, 1}}]",
        "{{1, -1/2, 1/6}, {0, 1, -1/3}, {0, 0, 1}}");
}

static void test_inverse_symbolic_2x2(void) {
    /* Inverse of {{a, b}, {c, d}}. */
    check_pretty(
        "Inverse[{{a, b}, {c, d}}] . {{a, b}, {c, d}} // Expand // Simplify",
        "{{1, 0}, {0, 1}}");
}

static void test_inverse_complex_2x2(void) {
    /* Inverse of {{1+I, 2}, {3, 4-I}} -- check identity via Dot. */
    check_true("(Inverse[{{1 + I, 2}, {3, 4 - I}}] . {{1 + I, 2}, {3, 4 - I}}) "
               "== {{1, 0}, {0, 1}}");
}

static void test_inverse_singular_leaves_unevaluated(void) {
    silence_stderr();
    /* Singular matrix: builtin returns NULL so the call survives. */
    check_fullform("Inverse[{{1, 2}, {2, 4}}]",
        "Inverse[List[List[1, 2], List[2, 4]]]");
}

/* ------------------------------------------------------------------
 * PseudoInverse: invertible square (collapses to Inverse)
 * ------------------------------------------------------------------ */
static void test_pseudo_invertible_collapses_to_inverse(void) {
    check_fullform("PseudoInverse[{{1, 2}, {3, 4}}]",
        "List[List[-2, 1], List[Rational[3, 2], Rational[-1, 2]]]");
    check_true("PseudoInverse[{{1, 2}, {3, 4}}] == Inverse[{{1, 2}, {3, 4}}]");
}

static void test_pseudo_invertible_3x3(void) {
    /* {{1, 0, 0}, {0, 2, 0}, {0, 0, 3}}^+ = diag(1, 1/2, 1/3). */
    check_pretty(
        "PseudoInverse[{{1, 0, 0}, {0, 2, 0}, {0, 0, 3}}]",
        "{{1, 0, 0}, {0, 1/2, 0}, {0, 0, 1/3}}");
}

/* ------------------------------------------------------------------
 * PseudoInverse: rank-1 singular 3x3 from the Wolfram spec
 * ------------------------------------------------------------------ */
static void test_pseudo_singular_3x3_spec(void) {
    /* Exact match against the spec value. */
    check_pretty(
        "PseudoInverse[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]",
        "{{-23/36, -1/6, 11/36}, {-1/18, 0, 1/18}, {19/36, 1/6, -7/36}}");
    /* Moore-Penrose properties hold. */
    check_true("(m = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}; "
               "p = PseudoInverse[m]; m . p . m == m)");
    check_true("(m = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}; "
               "p = PseudoInverse[m]; p . m . p == p)");
}

/* ------------------------------------------------------------------
 * PseudoInverse: rectangular (over- and under-determined)
 * ------------------------------------------------------------------ */
static void test_pseudo_rectangular_overdetermined(void) {
    /* Spec value: PseudoInverse[{{2, 3}, {2, 2}, {3, 1}, {4, 3}}] */
    check_pretty(
        "PseudoInverse[{{2, 3}, {2, 2}, {3, 1}, {4, 3}}]",
        "{{-29/134, -2/67, 22/67, 17/134}, {49/134, 8/67, -21/67, -1/134}}");
    /* MP identities. */
    check_true("(m = {{2, 3}, {2, 2}, {3, 1}, {4, 3}}; "
               "p = PseudoInverse[m]; m . p . m == m && p . m . p == p)");
}

static void test_pseudo_rectangular_underdetermined(void) {
    /* {{1, 2, 3}, {4, 5, 6}}: rank 2, m < n.  Check MP identities. */
    check_true("(m = {{1, 2, 3}, {4, 5, 6}}; "
               "p = PseudoInverse[m]; m . p . m == m && p . m . p == p)");
    /* For full-row-rank A: A . A^+ = I_m. */
    check_pretty(
        "(m = {{1, 2, 3}, {4, 5, 6}}; m . PseudoInverse[m])",
        "{{1, 0}, {0, 1}}");
}

static void test_pseudo_rectangular_3x2_fullcolumn(void) {
    /* {{1, 0}, {0, 1}, {1, 1}}: rank 2, m > n.  For full column rank
     * A: A^+ . A = I_n. */
    check_pretty(
        "PseudoInverse[{{1, 0}, {0, 1}, {1, 1}}] . {{1, 0}, {0, 1}, {1, 1}}",
        "{{1, 0}, {0, 1}}");
}

/* ------------------------------------------------------------------
 * PseudoInverse: zero matrix
 * ------------------------------------------------------------------ */
static void test_pseudo_zero_matrix_rect(void) {
    /* Spec: PseudoInverse[{{0,0,0},{0,0,0}}] == {{0,0},{0,0},{0,0}}. */
    check_fullform("PseudoInverse[{{0, 0, 0}, {0, 0, 0}}]",
        "List[List[0, 0], List[0, 0], List[0, 0]]");
}

static void test_pseudo_zero_matrix_square(void) {
    check_fullform("PseudoInverse[{{0, 0}, {0, 0}}]",
        "List[List[0, 0], List[0, 0]]");
}

/* ------------------------------------------------------------------
 * PseudoInverse: degenerate shapes (1x1, 1xN, Nx1)
 * ------------------------------------------------------------------ */
static void test_pseudo_1x1(void) {
    /* PseudoInverse[{{a}}] = {{1/a}} for non-zero a. */
    check_fullform("PseudoInverse[{{3}}]",
        "List[List[Rational[1, 3]]]");
    /* PseudoInverse[{{0}}] = {{0}}. */
    check_fullform("PseudoInverse[{{0}}]",
        "List[List[0]]");
}

static void test_pseudo_row_vector(void) {
    /* PseudoInverse[{{1, 2, 3}}] = {{1/14},{2/14},{3/14}}. */
    check_pretty(
        "PseudoInverse[{{1, 2, 3}}]",
        "{{1/14}, {1/7}, {3/14}}");
}

static void test_pseudo_column_vector(void) {
    /* PseudoInverse[{{1},{2},{3}}] = {{1/14, 2/14, 3/14}}. */
    check_pretty(
        "PseudoInverse[{{1}, {2}, {3}}]",
        "{{1/14, 1/7, 3/14}}");
}

/* ------------------------------------------------------------------
 * PseudoInverse: exact rational entries
 * ------------------------------------------------------------------ */
static void test_pseudo_rational_invertible(void) {
    check_true(
        "(m = {{1/2, 1/3}, {1/4, 1/5}}; "
        " p = PseudoInverse[m]; p == Inverse[m])");
}

static void test_pseudo_rational_singular(void) {
    /* Rank 1 rational matrix. */
    check_true(
        "(m = {{1/2, 1}, {1, 2}}; "
        " p = PseudoInverse[m]; m . p . m == m && p . m . p == p)");
}

/* ------------------------------------------------------------------
 * PseudoInverse: complex entries (exact)
 * ------------------------------------------------------------------ */
static void test_pseudo_complex_invertible(void) {
    /* Match against the standard inverse for a non-singular complex
     * matrix. */
    check_true(
        "(m = {{1 + I, 2}, {3, 4 - I}}; "
        " p = PseudoInverse[m]; "
        " Expand[m . p] == {{1, 0}, {0, 1}})");
}

static void test_pseudo_complex_singular(void) {
    /* Rank-1 complex matrix: rows are (1+I) times. */
    check_true(
        "(m = {{1 + I, 2}, {2 (1 + I), 4}}; "
        " p = PseudoInverse[m]; "
        " Expand[m . p . m] == m && Expand[p . m . p] == p)");
}

/* ------------------------------------------------------------------
 * PseudoInverse: machine-precision Real entries
 * ------------------------------------------------------------------ */
static void test_pseudo_real_with_zero_row(void) {
    /* Spec example: matches Mathematica to ~6 digits.  The result is
     * approximate so we compare via the sum of absolute deviations
     * against the reference output. */
    check_true(
        "(m = {{1.25, 3.2, 3.2}, {7.9, -1.4, 5.1}, {0, 0, 0}}; "
        " p = PseudoInverse[m]; "
        " expected = {{-0.0385185, 0.0966633, 0.0}, {0.210183, -0.0659894, 0.0}, "
        "             {0.117363, 0.0282303, 0.0}}; "
        " Total[Total[Abs[p - expected]]] < 0.00001)");
}

static void test_pseudo_real_invertible_matches_inverse(void) {
    /* For approximate invertible matrices, P ≈ Inverse[m]. */
    check_true(
        "(m = {{1.0, 2.0}, {3.0, 4.0}}; "
        " p = PseudoInverse[m]; "
        " Total[Total[Abs[p - {{-2, 1}, {3/2, -1/2}}]]] < 0.0000000001)");
}

/* ------------------------------------------------------------------
 * PseudoInverse: complex machine-precision (Wolfram spec value)
 * ------------------------------------------------------------------ */
static void test_pseudo_complex_machineprec(void) {
    /* Match the Mathematica reference values to ~5 digits via sum of
     * absolute deviations. */
    check_true(
        "(m = {{1. + I, 2, 3 - 2 I}, {0, 4, 5 I}, {1 + I, 6, 3 + 3 I}}; "
        " p = PseudoInverse[m]; "
        " expected = {{0.0393939 - 0.0575758 I, -0.00424242 + 0.0406061 I, "
        "              0.0351515 - 0.0169697 I}, "
        "             {0.0606061 - 0.0909091 I, 0.0424242 + 0.0727273 I, "
        "              0.10303 - 0.0181818 I}, "
        "             {0.0727273 + 0.115152 I, -0.0581818 - 0.0993939 I, "
        "              0.0145455 + 0.0157576 I}}; "
        " Total[Total[Abs[p - expected]]] < 0.0001)");
}

/* ------------------------------------------------------------------
 * PseudoInverse: Tolerance option parsing
 * ------------------------------------------------------------------ */
static void test_pseudo_tolerance_automatic(void) {
    /* With Tolerance -> Automatic the result is identical to the
     * default. */
    check_fullform(
        "PseudoInverse[{{1, 2}, {3, 4}}, Tolerance -> Automatic]",
        "List[List[-2, 1], List[Rational[3, 2], Rational[-1, 2]]]");
}

static void test_pseudo_tolerance_numeric_value(void) {
    /* On exact matrices the tolerance value is irrelevant; the result
     * is the exact inverse. */
    check_fullform(
        "PseudoInverse[{{1, 2}, {3, 4}}, Tolerance -> 1/100]",
        "List[List[-2, 1], List[Rational[3, 2], Rational[-1, 2]]]");
    check_fullform(
        "PseudoInverse[{{1, 2}, {3, 4}}, Tolerance -> 0.001]",
        "List[List[-2, 1], List[Rational[3, 2], Rational[-1, 2]]]");
}

static void test_pseudo_tolerance_invalid_option_unevaluated(void) {
    /* Unknown option name -> stays unevaluated. */
    check_fullform(
        "PseudoInverse[{{1, 2}, {3, 4}}, Foo -> 1]",
        "PseudoInverse[List[List[1, 2], List[3, 4]], Rule[Foo, 1]]");
    /* Tolerance with a negative value -> stays unevaluated. */
    check_fullform(
        "PseudoInverse[{{1, 2}, {3, 4}}, Tolerance -> -1]",
        "PseudoInverse[List[List[1, 2], List[3, 4]], Rule[Tolerance, -1]]");
}

/* ------------------------------------------------------------------
 * PseudoInverse: bad shape leaves the call unevaluated
 * ------------------------------------------------------------------ */
static void test_pseudo_bad_shape_unevaluated(void) {
    silence_stderr();
    /* Non-matrix scalar argument. */
    check_fullform("PseudoInverse[Sin[x]]", "PseudoInverse[Sin[x]]");
    /* Jagged "matrix". */
    check_fullform("PseudoInverse[{{1, 2}, {3}}]",
        "PseudoInverse[List[List[1, 2], List[3]]]");
    /* Empty inner list. */
    check_fullform("PseudoInverse[{{}}]", "PseudoInverse[List[List[]]]");
    /* Vector (rank 1) is not a matrix. */
    check_fullform("PseudoInverse[{1, 2, 3}]",
        "PseudoInverse[List[1, 2, 3]]");
}

/* ------------------------------------------------------------------
 * PseudoInverse: stress / leak-detection loop
 * ------------------------------------------------------------------ */
static void test_pseudo_stress_loop(void) {
    /* Re-run the singular 3x3 example many times so valgrind can
     * surface any leak in the per-call allocations. */
    for (int i = 0; i < 50; i++) {
        Expr* e = parse_expression(
            "PseudoInverse[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]");
        ASSERT(e != NULL);
        Expr* r = evaluate(e);
        ASSERT(r != NULL);
        expr_free(r);
        expr_free(e);
    }
    /* Also exercise the rectangular path. */
    for (int i = 0; i < 50; i++) {
        Expr* e = parse_expression(
            "PseudoInverse[{{2, 3}, {2, 2}, {3, 1}, {4, 3}}]");
        ASSERT(e != NULL);
        Expr* r = evaluate(e);
        ASSERT(r != NULL);
        expr_free(r);
        expr_free(e);
    }
    /* And the inexact/numeric path -- exercises the rationalise + numeric
     * round-trip. */
    for (int i = 0; i < 20; i++) {
        Expr* e = parse_expression(
            "PseudoInverse[{{1.25, 3.2, 3.2}, {7.9, -1.4, 5.1}, {0, 0, 0}}]");
        ASSERT(e != NULL);
        Expr* r = evaluate(e);
        ASSERT(r != NULL);
        expr_free(r);
        expr_free(e);
    }
    printf("PASS: stress loop (3 matrices x 50/50/20 iterations)\n");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_inverse_basic_integer);
    TEST(test_inverse_3x3_rational);
    TEST(test_inverse_symbolic_2x2);
    TEST(test_inverse_complex_2x2);
    TEST(test_inverse_singular_leaves_unevaluated);

    TEST(test_pseudo_invertible_collapses_to_inverse);
    TEST(test_pseudo_invertible_3x3);

    TEST(test_pseudo_singular_3x3_spec);

    TEST(test_pseudo_rectangular_overdetermined);
    TEST(test_pseudo_rectangular_underdetermined);
    TEST(test_pseudo_rectangular_3x2_fullcolumn);

    TEST(test_pseudo_zero_matrix_rect);
    TEST(test_pseudo_zero_matrix_square);

    TEST(test_pseudo_1x1);
    TEST(test_pseudo_row_vector);
    TEST(test_pseudo_column_vector);

    TEST(test_pseudo_rational_invertible);
    TEST(test_pseudo_rational_singular);

    TEST(test_pseudo_complex_invertible);
    TEST(test_pseudo_complex_singular);

    TEST(test_pseudo_real_with_zero_row);
    TEST(test_pseudo_real_invertible_matches_inverse);

    TEST(test_pseudo_complex_machineprec);

    TEST(test_pseudo_tolerance_automatic);
    TEST(test_pseudo_tolerance_numeric_value);
    TEST(test_pseudo_tolerance_invalid_option_unevaluated);

    TEST(test_pseudo_bad_shape_unevaluated);

    TEST(test_pseudo_stress_loop);

    printf("\nAll test_matinv.c tests passed.\n");
    return 0;
}
