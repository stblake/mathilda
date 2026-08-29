/*
 * test_schurdecomp.c -- unit tests for SchurDecomposition.
 *
 * SchurDecomposition[m] returns {q, t} with m == q . t . ConjugateTranspose[q],
 * q orthonormal and t (quasi-)upper-triangular; SchurDecomposition[{m, a}]
 * returns {q, s, p, t} with m == q.s.p^H and a == q.t.p^H.  The Schur form is
 * not unique (q and t are determined only up to the eigenvalue ordering and
 * unitary freedom within degenerate blocks), so the tests verify the defining
 * RELATIONS -- reconstruction residual, orthonormality, and triangularity --
 * rather than fixed factor values.
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

/* Evaluate `input` and require the result to be True. */
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

/* ---- Standard form: real ---- */

static void test_real_reconstruction(void) {
    check_true("Module[{m,q,t},"
        "m={{2.7,4.8,8.1},{-0.6,0,0},{0.1,0,0.3}};"
        "{q,t}=SchurDecomposition[m];"
        "Max[Abs[Flatten[m-q.t.ConjugateTranspose[q]]]]<10^-9]");
}
static void test_real_orthonormal(void) {
    check_true("Module[{m,q,t},"
        "m={{2.7,4.8,8.1},{-0.6,0,0},{0.1,0,0.3}};"
        "{q,t}=SchurDecomposition[m];"
        "Max[Abs[Flatten[q.ConjugateTranspose[q]-IdentityMatrix[3]]]]<10^-9]");
}
static void test_real_triangular(void) {
    /* This matrix has an all-real spectrum, so t is fully upper-triangular. */
    check_true("Module[{m,q,t},"
        "m={{2.7,4.8,8.1},{-0.6,0,0},{0.1,0,0.3}};"
        "{q,t}=SchurDecomposition[m];"
        "UpperTriangularMatrixQ[Chop[t,10^-9]]]");
}
static void test_machine_reconstruction(void) {
    check_true("Module[{m,q,t},"
        "m={{-1.2,2.7,3.8},{4.2,4.4,5.3},{3.5,7.6,6.8}};"
        "{q,t}=SchurDecomposition[m];"
        "Max[Abs[Flatten[m-q.t.ConjugateTranspose[q]]]]<10^-9]");
}

/* ---- Standard form: real matrix with complex eigenvalues (block form) ---- */

static void test_block_reconstruction(void) {
    check_true("Module[{m,q,t},"
        "m={{1.81066,0.31066,1.5},{-0.53033,2.03033,0.43934},"
        "{-0.96967,-0.53033,2.56066}};"
        "{q,t}=SchurDecomposition[m];"
        "Max[Abs[Flatten[m-q.t.ConjugateTranspose[q]]]]<10^-7]");
}
static void test_block_real_and_structure(void) {
    /* Default RealBlockDiagonalForm -> True: q, t are real, t is block upper-
     * triangular (nonzero first subdiagonal) but not upper-triangular. */
    check_true("Module[{m,q,t},"
        "m={{1.81066,0.31066,1.5},{-0.53033,2.03033,0.43934},"
        "{-0.96967,-0.53033,2.56066}};"
        "{q,t}=SchurDecomposition[m];"
        "FreeQ[Chop[q,10^-9],Complex] && UpperTriangularMatrixQ[t,-1] &&"
        " !UpperTriangularMatrixQ[t]]");
}
static void test_rbdf_false(void) {
    /* RealBlockDiagonalForm -> False makes t complex upper-triangular. */
    check_true("Module[{m,q,t},"
        "m={{1.81066,0.31066,1.5},{-0.53033,2.03033,0.43934},"
        "{-0.96967,-0.53033,2.56066}};"
        "{q,t}=SchurDecomposition[m,RealBlockDiagonalForm->False];"
        "Max[Abs[Flatten[m-q.t.ConjugateTranspose[q]]]]<10^-7 &&"
        " UpperTriangularMatrixQ[Chop[t,10^-9]]]");
}

/* ---- Standard form: complex (with a symbolic constant to numericalise) ---- */

static void test_complex_reconstruction(void) {
    check_true("Module[{m,q,t},"
        "m={{Pi,0.3},{I,1+1.5 I}};"
        "{q,t}=SchurDecomposition[m];"
        "Max[Abs[Flatten[m-q.t.ConjugateTranspose[q]]]]<10^-9]");
}
static void test_complex_triangular(void) {
    check_true("Module[{m,q,t},"
        "m={{Pi,0.3},{I,1+1.5 I}};"
        "{q,t}=SchurDecomposition[m];"
        "UpperTriangularMatrixQ[Chop[t,10^-9]]]");
}

/* ---- Pivoting ---- */

static void test_pivoting(void) {
    /* m . d == d . q . t . ConjugateTranspose[q] with the scaling/permutation d. */
    check_true("Module[{m,q,t,d},"
        "m=N[{{0,0,10^6},{0,10^-6,0},{-1,0,0}}];"
        "{q,t,d}=SchurDecomposition[m,Pivoting->True];"
        "Max[Abs[Flatten[m.d-d.q.t.ConjugateTranspose[q]]]]<10^-4]");
}

/* ---- Arbitrary precision (MPFR real-Schur path) ---- */

static void test_mpfr_reconstruction(void) {
    /* A 25-digit real matrix must reconstruct far tighter than machine could. */
    check_true("Module[{m,q,t},"
        "m=N[{{27/10,48/10,81/10},{-6/10,0,0},{1/10,0,3/10}},25];"
        "{q,t}=SchurDecomposition[m];"
        "Max[Abs[Flatten[m-q.t.ConjugateTranspose[q]]]]<10^-20]");
}
static void test_mpfr_complex_eigenvalue_block(void) {
    /* A high-precision real matrix with a complex-conjugate eigenvalue pair
     * exercises the MPFR 2x2-block path. */
    check_true("Module[{m,q,t},"
        "m=N[{{181066,31066,150000},{-53033,203033,43934},"
        "{-96967,-53033,256066}}/100000,25];"
        "{q,t}=SchurDecomposition[m];"
        "Max[Abs[Flatten[m-q.t.ConjugateTranspose[q]]]]<10^-18 &&"
        " UpperTriangularMatrixQ[Chop[t,10^-18],-1]]");
}

/* ---- Generalized (QZ) ---- */

static void test_generalized_real(void) {
    check_true("Module[{m,a,q,s,p,t},"
        "m={{0.5,1},{1.5,2}};a={{2.5,3},{3.5,4}};"
        "{q,s,p,t}=SchurDecomposition[{m,a}];"
        "Max[Abs[Flatten[m-q.s.ConjugateTranspose[p]]]]<10^-8 &&"
        " Max[Abs[Flatten[a-q.t.ConjugateTranspose[p]]]]<10^-8]");
}
static void test_generalized_complex(void) {
    check_true("Module[{m,a,q,s,p,t},"
        "m={{1.+2. I,0.3},{0.5,2.-1. I}};a={{2.,1. I},{1.,3.}};"
        "{q,s,p,t}=SchurDecomposition[{m,a}];"
        "Max[Abs[Flatten[m-q.s.ConjugateTranspose[p]]]]<10^-8 &&"
        " Max[Abs[Flatten[a-q.t.ConjugateTranspose[p]]]]<10^-8]");
}
static void test_generalized_triangular(void) {
    /* s and t are upper-triangular for a complex pencil. */
    check_true("Module[{m,a,q,s,p,t},"
        "m={{1.+2. I,0.3},{0.5,2.-1. I}};a={{2.,1. I},{1.,3.}};"
        "{q,s,p,t}=SchurDecomposition[{m,a}];"
        "UpperTriangularMatrixQ[Chop[s,10^-9]] &&"
        " UpperTriangularMatrixQ[Chop[t,10^-9]]]");
}

/* ---- Surfaces & edge cases ---- */

static void test_ndarray_input(void) {
    check_true("Module[{m,q,t},"
        "m=NDArray[{{-1.2,2.7,3.8},{4.2,4.4,5.3},{3.5,7.6,6.8}}];"
        "{q,t}=SchurDecomposition[m];"
        "Max[Abs[Flatten[Normal[m]-q.t.ConjugateTranspose[q]]]]<10^-9]");
}
static void test_exact_integer_numeric(void) {
    /* An exact numeric matrix has no closed-form Schur decomposition; it is
     * computed at machine precision and returns a numeric factor pair. */
    check_true("MatrixQ[First[SchurDecomposition[{{1,2},{3,4}}]], NumericQ]");
}
static void test_one_by_one(void) {
    check_true("Module[{q,t},{q,t}=SchurDecomposition[{{5.0}}];"
        "q=={{1.}} && t=={{5.}}]");
}
static void test_targetstructure_dense(void) {
    /* TargetStructure -> \"Structured\" returns dense matrices (no structured
     * type); the reconstruction still holds. */
    check_true("Module[{m,q,t},"
        "m={{-1.2,2.7,3.8},{4.2,4.4,5.3},{3.5,7.6,6.8}};"
        "{q,t}=SchurDecomposition[m,TargetStructure->\"Structured\"];"
        "Max[Abs[Flatten[m-q.t.ConjugateTranspose[q]]]]<10^-9]");
}

/* ---- Unevaluated / error cases ---- */

static void test_symbolic_unevaluated(void) {
    check_true("Head[SchurDecomposition[{{a,b},{c,d}}]]===SchurDecomposition");
}
static void test_arg_error(void) {
    check_true("Head[SchurDecomposition[]]===SchurDecomposition");
}
static void test_nonsquare(void) {
    check_true("Head[SchurDecomposition[{{1.,2.,3.},{4.,5.,6.}}]]===SchurDecomposition");
}
static void test_unknown_option(void) {
    check_true("Head[SchurDecomposition[{{1.,2.},{3.,4.}},Foo->True]]"
        "===SchurDecomposition");
}
static void test_docstring(void) {
    const char* doc = symtab_get_docstring("SchurDecomposition");
    ASSERT(doc != NULL);
    ASSERT(strstr(doc, "SchurDecomposition[m]") != NULL);
    printf("PASS: docstring present\n");
}

int main(void) {
    symtab_init();
    core_init();

    printf("Running SchurDecomposition tests...\n");
    TEST(test_real_reconstruction);
    TEST(test_real_orthonormal);
    TEST(test_real_triangular);
    TEST(test_machine_reconstruction);
    TEST(test_block_reconstruction);
    TEST(test_block_real_and_structure);
    TEST(test_rbdf_false);
    TEST(test_complex_reconstruction);
    TEST(test_complex_triangular);
    TEST(test_pivoting);
    TEST(test_mpfr_reconstruction);
    TEST(test_mpfr_complex_eigenvalue_block);
    TEST(test_generalized_real);
    TEST(test_generalized_complex);
    TEST(test_generalized_triangular);
    TEST(test_ndarray_input);
    TEST(test_exact_integer_numeric);
    TEST(test_one_by_one);
    TEST(test_targetstructure_dense);
    TEST(test_symbolic_unevaluated);
    TEST(test_arg_error);
    TEST(test_nonsquare);
    TEST(test_unknown_option);
    TEST(test_docstring);

    printf("\nAll SchurDecomposition tests passed.\n");
    return 0;
}
