/*
 * test_eigen.c -- unit tests for Eigenvalues and Eigenvectors.
 *
 * The Mathematica reference output for several of the spec examples uses
 * a different canonical form than Mathilda's printer (e.g.
 * `3/2 (5 + Sqrt[33])` vs `1/2 (15 + 3 Sqrt[33])`), so most tests check
 * either the exact Mathilda-canonical string or verify the
 * mathematical relationship m.v == lambda*v / Det[m - lambda I] == 0.
 *
 * Numerical eigenvector tests for n >= 3 are deliberately omitted: the
 * exact-symbolic Cardano roots that come out of Solve are too costly
 * to substitute into m - lambda I and RowReduce, and we have no purely
 * numerical eigenvector routine yet.  Limitations are documented in the
 * Eigenvalues changelog.
 */
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

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
        printf("PASS: %s\n", input);
    }
    free(s);
    expr_free(r);
    expr_free(e);
}

/* Check that evaluating `input` yields True. */
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

/* Approximate-equality check for a numerical eigenvalue list result. */
static void run_check_numeric(const char* input,
                              const double* expected,
                              size_t n) {
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    ASSERT(r->type == EXPR_FUNCTION);
    ASSERT(r->data.function.arg_count == n);
    for (size_t i = 0; i < n; i++) {
        Expr* v = r->data.function.args[i];
        double d = 0;
        if (v->type == EXPR_REAL) d = v->data.real;
        else if (v->type == EXPR_INTEGER) d = (double)v->data.integer;
        else {
            char* s = expr_to_string(v);
            printf("FAIL: %s -- entry %zu is not numerical: %s\n",
                   input, i, s);
            free(s);
            ASSERT(0);
        }
        double err = fabs(d - expected[i]);
        double tol = 1e-4 * (1 + fabs(expected[i]));
        if (err > tol) {
            printf("FAIL: %s -- entry %zu: expected %g, got %g (err %g)\n",
                   input, i, expected[i], d, err);
            ASSERT(0);
        }
    }
    printf("PASS: %s -> numeric match\n", input);
    expr_free(r);
    expr_free(e);
}

/* ============================ Tests ============================ */

/* Spec: Eigenvalues of an exact matrix. */
void test_eigenvalues_exact_3x3(void) {
    run_test("Eigenvalues[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]",
             "List[Times[Rational[1, 2], Plus[15, Times[3, Power[33, Rational[1, 2]]]]], "
             "Times[Rational[1, 2], Plus[15, Times[-3, Power[33, Rational[1, 2]]]]], 0]");
}

/* Spec: symbolic 2x2 eigenvalues -- verify via char polynomial.
 * Uses Simplify because the Sqrt argument doesn't auto-expand. */
void test_eigenvalues_symbolic_2x2(void) {
    run_check_true(
        "(mt = {{a, b}, {c, d}}; "
        " evs = Eigenvalues[mt]; "
        " r = And @@ (Expand[Det[mt - #*IdentityMatrix[2]]] == 0 & /@ evs); "
        " Clear[mt, evs]; "
        " Simplify[r])");
}

/* Spec: machine-precision 3x3 eigenvalues with chop. */
void test_eigenvalues_numeric_3x3(void) {
    double expected[3] = { 6.60674, 4.52536, 0.667901 };
    run_check_numeric(
        "Eigenvalues[{{1.1, 2.2, 3.25}, {0.76, 4.6, 5}, {0.1, 0.1, 6.1}}]",
        expected, 3);
}

/* Spec: repeated eigenvalues. */
void test_eigenvalues_repeated(void) {
    run_test("Eigenvalues[{{1, 0, 1}, {0, 1, 0}, {0, 0, 1}}]",
             "List[1, 1, 1]");
}

/* Spec: IdentityMatrix gives all-one eigenvalues. */
void test_eigenvalues_identity(void) {
    run_test("Eigenvalues[IdentityMatrix[3]]", "List[1, 1, 1]");
    run_test("Eigenvalues[IdentityMatrix[5]]",
             "List[1, 1, 1, 1, 1]");
    run_test("Eigenvalues[IdentityMatrix[12]]",
             "List[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]");
}

/* Spec: 4x4 block matrix, eigenvalues sorted by descending |lambda|. */
void test_eigenvalues_4x4_block(void) {
    run_test("Eigenvalues[{{7/2, 0, 1/2, 0}, {0, 3, 0, 1}, "
             "{1/2, 0, 7/2, 0}, {0, 1, 0, 3}}]",
             "List[4, 4, 3, 2]");
}

/* Spec: Eigenvalues[m, k] -- first k by absolute value. */
void test_eigenvalues_k_positive(void) {
    run_test("Eigenvalues[{{7/2, 0, 1/2, 0}, {0, 3, 0, 1}, "
             "{1/2, 0, 7/2, 0}, {0, 1, 0, 3}}, 3]",
             "List[4, 4, 3]");
}

/* Spec: Eigenvalues[m, -k] -- smallest k by absolute value. */
void test_eigenvalues_k_negative(void) {
    run_test("Eigenvalues[{{7/2, 0, 1/2, 0}, {0, 3, 0, 1}, "
             "{1/2, 0, 7/2, 0}, {0, 1, 0, 3}}, -3]",
             "List[4, 3, 2]");
}

/* Spec: Eigenvalues[m, UpTo[k]]. */
void test_eigenvalues_upto(void) {
    run_test("Eigenvalues[{{7/2, 0, 1/2, 0}, {0, 3, 0, 1}, "
             "{1/2, 0, 7/2, 0}, {0, 1, 0, 3}}, UpTo[2]]",
             "List[4, 4]");
    run_test("Eigenvalues[{{7/2, 0, 1/2, 0}, {0, 3, 0, 1}, "
             "{1/2, 0, 7/2, 0}, {0, 1, 0, 3}}, UpTo[10]]",
             "List[4, 4, 3, 2]");
}

/* Spec: generalised eigenvalues with finite + Infinity components. */
void test_eigenvalues_generalized_exact(void) {
    run_test("Eigenvalues[{{{1, 1, 1}, {1, 0, 1}, {0, 0, 1}}, "
             "{{0, 1, 1}, {0, 1, 1}, {1, 0, 0}}}]",
             "List[Infinity, "
             "Times[Rational[1, 2], Plus[1, Power[5, Rational[1, 2]]]], "
             "Times[Rational[1, 2], Plus[1, Times[-1, Power[5, Rational[1, 2]]]]]]");
}

/* Spec: generalised eigenvalues symbolic 2x2 -- verify via char poly. */
void test_eigenvalues_generalized_symbolic(void) {
    run_check_true(
        "(am = {{x, 1 + x}, {1 - x, x}}; "
        " bm = {{1, 1}, {1, 2 x}}; "
        " evs = Eigenvalues[{am, bm}]; "
        " r = And @@ (Simplify[Together[Det[am - #*bm]]] == 0 & /@ evs); "
        " Clear[am, bm, evs]; "
        " r)");
}

/* Spec: generalised eigenvalues take smallest 2.
 * This case has a degenerate char poly: one finite root, one zero root,
 * and one Infinity (since det of b's lead is 0).
 * Per the spec the answer is {1, 0}.  We verify a weaker invariant: the
 * computed list has length 2 with at least a 0 entry and a 1 entry. */
void test_eigenvalues_generalized_k_negative(void) {
    run_test("Eigenvalues[{{{1,2,3},{4,5,6},{7,8,9}}, "
             "{{11,12,13},{1,15,16},{17,18,19}}}, -2]",
             "List[1, 0]");
}

/* Eigenvectors: defective matrix gets zero padding in-line. */
void test_eigenvectors_defective(void) {
    run_test("Eigenvectors[{{2, 1, 0}, {0, 2, 0}, {0, 0, 1}}]",
             "List[List[1, 0, 0], List[0, 0, 0], List[0, 0, 1]]");
}

/* Eigenvectors: identity-defective case (eigenvalue 1 with mult 3, rank 1
 * defect). */
void test_eigenvectors_identity_defective(void) {
    run_test("Eigenvectors[{{1, 0, 1}, {0, 1, 0}, {0, 0, 1}}]",
             "List[List[1, 0, 0], List[0, 1, 0], List[0, 0, 0]]");
}

/* Eigenvectors: m.v == lambda*v for the 4x4 block case. */
void test_eigenvectors_4x4_verify(void) {
    run_check_true(
        "(mt = {{7/2, 0, 1/2, 0}, {0, 3, 0, 1}, {1/2, 0, 7/2, 0}, {0, 1, 0, 3}}; "
        " evs = Eigenvalues[mt]; "
        " evecs = Eigenvectors[mt]; "
        " r = And @@ Table[mt . evecs[[i]] == evs[[i]] * evecs[[i]], {i, 1, 4}]; "
        " Clear[mt, evs, evecs]; "
        " r)");
}

/* Eigenvectors: m.v == lambda*v for the singular 3x3 matrix. */
void test_eigenvectors_3x3_singular_verify(void) {
    run_check_true(
        "(mt = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}; "
        " evs = Eigenvalues[mt]; "
        " evecs = Eigenvectors[mt]; "
        " r = And @@ Table["
        "    Together[Expand[mt . evecs[[i]] - evs[[i]]*evecs[[i]]]] == "
        "      {0, 0, 0}, {i, 1, 3}]; "
        " Clear[mt, evs, evecs]; "
        " r)");
}

/* Spec: norm of unit-vector eigenvectors for inexact 2x2. */
void test_eigenvectors_numeric_2x2_normalized(void) {
    run_test("Norm /@ Eigenvectors[{{1., 2.}, {2., 1.}}]",
             "List[1.0, 1.0]");
}

/* Eigenvectors: exact 2x2 not normalized. */
void test_eigenvectors_2x2_exact_norm(void) {
    run_test("Norm /@ Eigenvectors[{{1, 2}, {2, 1}}]",
             "List[Power[2, Rational[1, 2]], Power[2, Rational[1, 2]]]");
}

/* Spec: generalised eigenvectors -- verify m.v == lambda*a.v. */
void test_eigenvectors_generalized_verify(void) {
    run_check_true(
        "(am = {{1, 1, 1}, {1, 0, 1}, {0, 0, 1}}; "
        " bm = {{0, 1, 1}, {0, 1, 1}, {1, 0, 0}}; "
        " evs = Eigenvalues[{am, bm}]; "
        " evecs = Eigenvectors[{am, bm}]; "
        " r = And @@ Table["
        "    If[evs[[i]] === Infinity, True, "
        "       Expand[am . evecs[[i]] - evs[[i]]*(bm . evecs[[i]])] == "
        "         {0, 0, 0}], "
        "    {i, 1, 3}]; "
        " Clear[am, bm, evs, evecs]; "
        " r)");
}

/* Generalised: take 2 smallest from 3x3. */
void test_eigenvectors_generalized_k(void) {
    run_test("Eigenvectors[{{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, "
             "{{11, 12, 13}, {1, 15, 16}, {17, 18, 19}}}, -2]",
             "List[List[0, 0, 0], List[1, -2, 1]]");
}

/* Eigenvalues of a triangular matrix are its diagonal entries. */
void test_eigenvalues_triangular(void) {
    run_test("Eigenvalues[{{5, 1, 2, 3}, {0, 3, 7, 4}, "
             "{0, 0, 2, 9}, {0, 0, 0, 1}}]",
             "List[5, 3, 2, 1]");
}

/* Diagonal matrix eigenvalues. */
void test_eigenvalues_diagonal(void) {
    run_test("Eigenvalues[DiagonalMatrix[{7, 3, 1, 5}]]",
             "List[7, 5, 3, 1]");
}

/* Eigenvalues on a non-square matrix should be unevaluated. */
void test_eigenvalues_unevaluated(void) {
    run_test("Eigenvalues[{{1, 2, 3}, {4, 5, 6}}]",
             "Eigenvalues[List[List[1, 2, 3], List[4, 5, 6]]]");
    run_test("Eigenvalues[x]", "Eigenvalues[x]");
}

/* Docstrings registered? */
void test_docstring(void) {
    run_test("Head[Information[\"Eigenvalues\"]]", "String");
    run_test("Head[Information[\"Eigenvectors\"]]", "String");
}

/* ==================== Additional thorough tests ==================== */

/* 1x1 numeric: trivially the single matrix element / unit vector. */
void test_eigenvalues_1x1_numeric(void) {
    run_test("Eigenvalues[{{5}}]", "List[5]");
    run_test("Eigenvectors[{{5}}]", "List[List[1]]");
}

/* 1x1 symbolic: the lone entry is the eigenvalue. */
void test_eigenvalues_1x1_symbolic(void) {
    run_test("Eigenvalues[{{a}}]", "List[a]");
    run_test("Eigenvectors[{{a}}]", "List[List[1]]");
}

/* Complex eigenvalues from a 2D rotation: roots are +-I. */
void test_eigenvalues_complex_rotation(void) {
    run_test("Eigenvalues[{{0, -1}, {1, 0}}]",
             "List[Complex[0, -1], Complex[0, 1]]");
}

/* Complex eigenvectors verify m.v == lambda*v. */
void test_eigenvectors_complex_rotation_verify(void) {
    run_check_true(
        "(mt = {{0, -1}, {1, 0}}; "
        " evs = Eigenvalues[mt]; "
        " evecs = Eigenvectors[mt]; "
        " r = And @@ Table[Expand[mt . evecs[[i]] - evs[[i]] * evecs[[i]]] "
        "                    == {0, 0}, {i, 1, 2}]; "
        " Clear[mt, evs, evecs]; r)");
}

/* Zero matrix: zero eigenvalues with full identity eigenvector basis. */
void test_eigenvalues_zero_matrix(void) {
    run_test("Eigenvalues[{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}]",
             "List[0, 0, 0]");
    run_test("Eigenvectors[{{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}]",
             "List[List[1, 0, 0], List[0, 1, 0], List[0, 0, 1]]");
}

/* Nilpotent shift matrix: all eigenvalues zero, single Jordan block. */
void test_eigenvalues_nilpotent_3x3(void) {
    run_test("Eigenvalues[{{0, 1, 0}, {0, 0, 1}, {0, 0, 0}}]",
             "List[0, 0, 0]");
    run_test("Eigenvectors[{{0, 1, 0}, {0, 0, 1}, {0, 0, 0}}]",
             "List[List[1, 0, 0], List[0, 0, 0], List[0, 0, 0]]");
}

/* Defective 2x2 Jordan block: lambda=2 with multiplicity 2 but only one
 * independent eigenvector; expect zero-padded second slot. */
void test_eigenvalues_jordan_2x2(void) {
    run_test("Eigenvalues[{{2, 1}, {0, 2}}]", "List[2, 2]");
    run_test("Eigenvectors[{{2, 1}, {0, 2}}]",
             "List[List[1, 0], List[0, 0]]");
}

/* Trace identity: Sum[eigenvalues] == Tr[m] (after simplification). */
void test_eigenvalues_trace_identity(void) {
    run_check_true(
        "Simplify[Plus @@ Eigenvalues[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]] "
        "  == Tr[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]");
    run_check_true(
        "Simplify[Plus @@ Eigenvalues[{{7/2, 0, 1/2, 0}, {0, 3, 0, 1}, "
        "                              {1/2, 0, 7/2, 0}, {0, 1, 0, 3}}]] "
        "  == Tr[{{7/2, 0, 1/2, 0}, {0, 3, 0, 1}, "
        "         {1/2, 0, 7/2, 0}, {0, 1, 0, 3}}]");
}

/* Determinant identity: Product[eigenvalues] == Det[m]. */
void test_eigenvalues_det_identity(void) {
    run_check_true(
        "Together[Times @@ Eigenvalues[{{1, 2}, {3, 4}}]] == "
        "  Det[{{1, 2}, {3, 4}}]");
    run_check_true(
        "Times @@ Eigenvalues[{{7/2, 0, 1/2, 0}, {0, 3, 0, 1}, "
        "                      {1/2, 0, 7/2, 0}, {0, 1, 0, 3}}] == "
        "  Det[{{7/2, 0, 1/2, 0}, {0, 3, 0, 1}, "
        "       {1/2, 0, 7/2, 0}, {0, 1, 0, 3}}]");
}

/* Lower-triangular matrix: eigenvalues are diagonal, descending |lambda|. */
void test_eigenvalues_lower_triangular(void) {
    run_test("Eigenvalues[{{5, 0, 0, 0}, {1, 3, 0, 0}, "
             "{2, 7, 2, 0}, {3, 4, 9, 1}}]",
             "List[5, 3, 2, 1]");
}

/* k larger than n saturates to the full list. */
void test_eigenvalues_k_saturating(void) {
    run_test("Eigenvalues[{{7/2, 0, 1/2, 0}, {0, 3, 0, 1}, "
             "{1/2, 0, 7/2, 0}, {0, 1, 0, 3}}, 10]",
             "List[4, 4, 3, 2]");
    run_test("Eigenvalues[{{7/2, 0, 1/2, 0}, {0, 3, 0, 1}, "
             "{1/2, 0, 7/2, 0}, {0, 1, 0, 3}}, -10]",
             "List[4, 4, 3, 2]");
}

/* k = 0 / UpTo[0] yield empty lists. */
void test_eigenvalues_k_zero(void) {
    run_test("Eigenvalues[{{1, 2}, {3, 4}}, 0]", "List[]");
    run_test("Eigenvalues[{{1, 2}, {3, 4}}, UpTo[0]]", "List[]");
}

/* N[Eigenvalues[m]] returns a fully-numerical list. */
void test_eigenvalues_N_numericalize(void) {
    double expected[3] = { 16.1168, -1.11684, 0.0 };
    run_check_numeric(
        "N[Eigenvalues[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]]",
        expected, 3);
}

/* Eigenvectors of IdentityMatrix[n] are the standard basis. */
void test_eigenvectors_identity_basis(void) {
    run_test("Eigenvectors[IdentityMatrix[3]]",
             "List[List[1, 0, 0], List[0, 1, 0], List[0, 0, 1]]");
}

/* Eigenvectors of a DiagonalMatrix: basis vectors in descending-|lambda|
 * order matching the eigenvalue sort. */
void test_eigenvectors_diagonal_basis(void) {
    run_test("Eigenvectors[DiagonalMatrix[{3, 5, 7}]]",
             "List[List[0, 0, 1], List[0, 1, 0], List[1, 0, 0]]");
}

/* Symmetric tridiagonal 3x3: closed-form Sqrt[2] spectrum. */
void test_eigenvalues_symmetric_tridiagonal(void) {
    run_test("Eigenvalues[{{2, -1, 0}, {-1, 2, -1}, {0, -1, 2}}]",
             "List[Times[Rational[1, 2], Plus[4, Times[2, Power[2, Rational[1, 2]]]]], "
             "2, "
             "Times[Rational[1, 2], Plus[4, Times[-2, Power[2, Rational[1, 2]]]]]]");
}

/* Symmetric tridiagonal 3x3: corresponding eigenvectors verify m.v == lambda*v. */
void test_eigenvectors_symmetric_tridiagonal_verify(void) {
    run_check_true(
        "(mt = {{2, -1, 0}, {-1, 2, -1}, {0, -1, 2}}; "
        " evs = Eigenvalues[mt]; "
        " evecs = Eigenvectors[mt]; "
        " r = And @@ Table[Expand[mt . evecs[[i]] - evs[[i]] * evecs[[i]]] "
        "                    == {0, 0, 0}, {i, 1, 3}]; "
        " Clear[mt, evs, evecs]; r)");
}

/* Cubics->False keeps closed-form answer when one root is 0 (degree drops
 * to a quadratic). */
void test_eigenvalues_cubics_option_singular(void) {
    run_test("Eigenvalues[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, Cubics -> False]",
             "List[Times[Rational[1, 2], Plus[15, Times[3, Power[33, Rational[1, 2]]]]], "
             "Times[Rational[1, 2], Plus[15, Times[-3, Power[33, Rational[1, 2]]]]], 0]");
}

/* Quartics->False: 4x4 singular matrix degenerates to a quadratic, so the
 * closed form survives.  Exercises the option-parsing pipeline. */
void test_eigenvalues_quartics_option_singular(void) {
    run_test("Eigenvalues[{{1,2,3,4},{5,6,7,8},{9,10,11,12},{13,14,15,16}}, "
             "Quartics -> False]",
             "List[Times[Rational[1, 2], Plus[34, Times[6, Power[41, Rational[1, 2]]]]], "
             "Times[Rational[1, 2], Plus[34, Times[-6, Power[41, Rational[1, 2]]]]], 0, 0]");
}

/* Non-square arg with k stays unevaluated. */
void test_eigenvalues_nonsquare_with_k_unevaluated(void) {
    run_test("Eigenvalues[{{1, 2, 3}, {4, 5, 6}}, 2]",
             "Eigenvalues[List[List[1, 2, 3], List[4, 5, 6]], 2]");
}

/* Empty list / non-matrix scalars stay unevaluated. */
void test_eigenvalues_empty_and_scalar_unevaluated(void) {
    run_test("Eigenvalues[{}]", "Eigenvalues[List[]]");
    run_test("Eigenvalues[5]", "Eigenvalues[5]");
}

/* Bigint matrix entries: the closed-form radical does not auto-simplify back
 * to the diagonal entries, but the spectrum still obeys Tr / Det identities,
 * which is what we verify. */
void test_eigenvalues_bigint_diagonal(void) {
    run_check_true(
        "Simplify[Plus @@ Eigenvalues[DiagonalMatrix[{12345678901, 98765432109}]]] "
        "  == 12345678901 + 98765432109");
    run_check_true(
        "Simplify[Times @@ Eigenvalues[DiagonalMatrix[{12345678901, 98765432109}]]] "
        "  == 12345678901 * 98765432109");
}

/* Rational-entry matrix produces a closed-form spectrum. */
void test_eigenvalues_rational_2x2(void) {
    run_test("Eigenvalues[{{1/2, 1/3}, {1/4, 1/5}}]",
             "List[Times[Rational[1, 120], Plus[42, Times[2, Power[381, Rational[1, 2]]]]], "
             "Times[Rational[1, 120], Plus[42, Times[-2, Power[381, Rational[1, 2]]]]]]");
}

/* 5x5 upper triangular: spectrum is the diagonal (descending |lambda|). */
void test_eigenvalues_upper_triangular_5x5(void) {
    run_test("Eigenvalues[{{1,2,3,4,5},{0,2,3,4,5},{0,0,3,4,5},"
             "{0,0,0,4,5},{0,0,0,0,5}}]",
             "List[5, 4, 3, 2, 1]");
}

/* -------- Phase 1: Method option parsing (no kernel dispatch yet) --------
 *
 * Method -> Automatic and the new Direct/Arnoldi/Banded/FEAST values must
 * all parse cleanly.  Until Phases 2-5 land, every numeric-matrix call
 * with an explicit Method falls back to the symbolic char-poly path; the
 * result must equal the no-Method result so Phase 1 introduces no
 * behavioural regression.
 *
 * Unknown Method values are tolerated (one warning, same fallback).
 * Symbolic matrices ignore Method entirely. */
static void run_pair_same(const char* a, const char* b) {
    Expr* ea = parse_expression(a);
    Expr* eb = parse_expression(b);
    ASSERT(ea && eb);
    Expr* ra = evaluate(ea);
    Expr* rb = evaluate(eb);
    char* sa = expr_to_string_fullform(ra);
    char* sb = expr_to_string_fullform(rb);
    if (strcmp(sa, sb) != 0) {
        printf("FAIL: %s\n   not equal to: %s\n   got:    %s\n   vs:     %s\n",
               a, b, sa, sb);
        ASSERT(0);
    } else {
        printf("PASS: %s == %s\n", a, b);
    }
    free(sa); free(sb);
    expr_free(ra); expr_free(rb);
    expr_free(ea); expr_free(eb);
}

void test_eigenvalues_method_automatic_numeric(void) {
    run_pair_same(
        "Eigenvalues[{{1.0, 2.0}, {3.0, 4.0}}, Method -> Automatic]",
        "Eigenvalues[{{1.0, 2.0}, {3.0, 4.0}}]");
}

void test_eigenvalues_method_direct_numeric_fallback(void) {
    run_pair_same(
        "Eigenvalues[{{1.0, 2.0}, {3.0, 4.0}}, Method -> \"Direct\"]",
        "Eigenvalues[{{1.0, 2.0}, {3.0, 4.0}}]");
}

void test_eigenvalues_method_arnoldi_numeric_fallback(void) {
    run_pair_same(
        "Eigenvalues[{{1.0, 2.0}, {3.0, 4.0}}, Method -> \"Arnoldi\"]",
        "Eigenvalues[{{1.0, 2.0}, {3.0, 4.0}}]");
}

void test_eigenvalues_method_banded_numeric_fallback(void) {
    run_pair_same(
        "Eigenvalues[{{2.0, -1.0, 0.0}, {-1.0, 2.0, -1.0}, "
        "{0.0, -1.0, 2.0}}, Method -> \"Banded\"]",
        "Eigenvalues[{{2.0, -1.0, 0.0}, {-1.0, 2.0, -1.0}, "
        "{0.0, -1.0, 2.0}}]");
}

void test_eigenvalues_method_feast_numeric_fallback(void) {
    run_pair_same(
        "Eigenvalues[{{2.0, -1.0}, {-1.0, 2.0}}, "
        "Method -> {\"FEAST\", Interval -> {0.0, 5.0}}]",
        "Eigenvalues[{{2.0, -1.0}, {-1.0, 2.0}}]");
}

void test_eigenvalues_method_unknown_string_numeric_fallback(void) {
    /* An unknown Method string falls back to the symbolic path with a
     * warning -- it must not leave the call unevaluated. */
    run_pair_same(
        "Eigenvalues[{{1.0, 0.0}, {0.0, 2.0}}, Method -> \"NoSuchMethod\"]",
        "Eigenvalues[{{1.0, 0.0}, {0.0, 2.0}}]");
}

void test_eigenvalues_method_ignored_for_symbolic(void) {
    /* Symbolic matrices ignore Method and emit no warning. */
    run_pair_same(
        "Eigenvalues[{{1, 2}, {3, 4}}, Method -> \"Direct\"]",
        "Eigenvalues[{{1, 2}, {3, 4}}]");
    run_pair_same(
        "Eigenvalues[{{a, b}, {c, d}}, Method -> \"Arnoldi\"]",
        "Eigenvalues[{{a, b}, {c, d}}]");
}

void test_eigenvalues_method_value_list_with_suboptions(void) {
    /* {Method -> {"Arnoldi", MaxIterations -> 20}} must parse: the head
     * element classifies the method, sub-options are accepted (and
     * ignored in Phase 1).  Result equals the no-Method baseline. */
    run_pair_same(
        "Eigenvalues[{{1.0, 2.0}, {3.0, 4.0}}, "
        "Method -> {\"Arnoldi\", MaxIterations -> 20, Tolerance -> 1.0*^-8}]",
        "Eigenvalues[{{1.0, 2.0}, {3.0, 4.0}}]");
}

void test_eigenvectors_method_fallback(void) {
    /* Eigenvectors mirrors Eigenvalues' Method behaviour.
     *
     * Note: "Direct" and the no-Method (Automatic) path return string-
     * identical results on this diagonal input.  "Arnoldi" now uses a
     * real numerical kernel (Phase 3) whose output is mathematically
     * equivalent but not string-equal (CGS roundoff noise + sign
     * convention on the second eigenvector), so we limit the strict
     * string comparison to Direct; Arnoldi's eigenvector residuals are
     * verified by tests/test_mateigen_arnoldi.c. */
    run_pair_same(
        "Eigenvectors[{{2.0, 0.0}, {0.0, 3.0}}, Method -> \"Direct\"]",
        "Eigenvectors[{{2.0, 0.0}, {0.0, 3.0}}]");
}

void test_eigenvalues_method_combines_with_k_and_cubics(void) {
    /* Method must coexist with the other options and the k-spec. */
    run_pair_same(
        "Eigenvalues[{{1.0, 2.0, 3.0}, {0.0, 4.0, 5.0}, {0.0, 0.0, 6.0}}, "
        "2, Method -> \"Direct\", Cubics -> True]",
        "Eigenvalues[{{1.0, 2.0, 3.0}, {0.0, 4.0, 5.0}, {0.0, 0.0, 6.0}}, "
        "2, Cubics -> True]");
}

int main(void) {
    symtab_init();
    core_init();

    printf("Running Eigenvalues/Eigenvectors tests...\n");
    TEST(test_eigenvalues_exact_3x3);
    TEST(test_eigenvalues_symbolic_2x2);
    TEST(test_eigenvalues_numeric_3x3);
    TEST(test_eigenvalues_repeated);
    TEST(test_eigenvalues_identity);
    TEST(test_eigenvalues_4x4_block);
    TEST(test_eigenvalues_k_positive);
    TEST(test_eigenvalues_k_negative);
    TEST(test_eigenvalues_upto);
    TEST(test_eigenvalues_generalized_exact);
    TEST(test_eigenvalues_generalized_symbolic);
    TEST(test_eigenvalues_generalized_k_negative);
    TEST(test_eigenvectors_defective);
    TEST(test_eigenvectors_identity_defective);
    TEST(test_eigenvectors_4x4_verify);
    TEST(test_eigenvectors_3x3_singular_verify);
    TEST(test_eigenvectors_numeric_2x2_normalized);
    TEST(test_eigenvectors_2x2_exact_norm);
    TEST(test_eigenvectors_generalized_verify);
    TEST(test_eigenvectors_generalized_k);
    TEST(test_eigenvalues_triangular);
    TEST(test_eigenvalues_diagonal);
    TEST(test_eigenvalues_unevaluated);
    TEST(test_docstring);

    /* Additional thorough tests */
    TEST(test_eigenvalues_1x1_numeric);
    TEST(test_eigenvalues_1x1_symbolic);
    TEST(test_eigenvalues_complex_rotation);
    TEST(test_eigenvectors_complex_rotation_verify);
    TEST(test_eigenvalues_zero_matrix);
    TEST(test_eigenvalues_nilpotent_3x3);
    TEST(test_eigenvalues_jordan_2x2);
    TEST(test_eigenvalues_trace_identity);
    TEST(test_eigenvalues_det_identity);
    TEST(test_eigenvalues_lower_triangular);
    TEST(test_eigenvalues_k_saturating);
    TEST(test_eigenvalues_k_zero);
    TEST(test_eigenvalues_N_numericalize);
    TEST(test_eigenvectors_identity_basis);
    TEST(test_eigenvectors_diagonal_basis);
    TEST(test_eigenvalues_symmetric_tridiagonal);
    TEST(test_eigenvectors_symmetric_tridiagonal_verify);
    TEST(test_eigenvalues_cubics_option_singular);
    TEST(test_eigenvalues_quartics_option_singular);
    TEST(test_eigenvalues_nonsquare_with_k_unevaluated);
    TEST(test_eigenvalues_empty_and_scalar_unevaluated);
    TEST(test_eigenvalues_bigint_diagonal);
    TEST(test_eigenvalues_rational_2x2);
    TEST(test_eigenvalues_upper_triangular_5x5);

    /* Phase 1: Method option parsing + warning fallback. */
    TEST(test_eigenvalues_method_automatic_numeric);
    TEST(test_eigenvalues_method_direct_numeric_fallback);
    TEST(test_eigenvalues_method_arnoldi_numeric_fallback);
    TEST(test_eigenvalues_method_banded_numeric_fallback);
    TEST(test_eigenvalues_method_feast_numeric_fallback);
    TEST(test_eigenvalues_method_unknown_string_numeric_fallback);
    TEST(test_eigenvalues_method_ignored_for_symbolic);
    TEST(test_eigenvalues_method_value_list_with_suboptions);
    TEST(test_eigenvectors_method_fallback);
    TEST(test_eigenvalues_method_combines_with_k_and_cubics);

    printf("All eigen tests passed!\n");
    return 0;
}
