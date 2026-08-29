/*
 * test_jordandecomp.c -- unit tests for JordanDecomposition[m].
 *
 * JordanDecomposition[m] returns {s, j} with m == s . j . Inverse[s].  The
 * generalized eigenvectors in s are not canonical (any valid Jordan basis is
 * accepted), so the exact tests verify the *relation* m.s == s.j together with
 * the deterministic Jordan form j, rather than a fixed s.  Numeric tests check
 * the residual m.s - s.j and the block structure; exact printed values are not
 * asserted (see the JordanDecomposition changelog).
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

/* Evaluate `input` and require its FullForm to equal `expected`. */
static void check_fullform(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    char* s = expr_to_string_fullform(r);
    if (strcmp(s, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, s);
        ASSERT(0);
    } else {
        printf("PASS: %s\n", input);
    }
    free(s);
    expr_free(r);
    expr_free(e);
}

/* ---- Exact / symbolic ---- */

static void test_exact_3x3_relation(void) {
    check_true("Module[{s,j,m},m={{27,48,81},{-6,0,0},{1,0,3}};"
               "{s,j}=JordanDecomposition[m];"
               "Simplify[m.s-s.j]==ConstantArray[0,{3,3}]]");
}
static void test_exact_3x3_jordan_form(void) {
    check_true("Module[{s,j},{s,j}=JordanDecomposition[{{27,48,81},{-6,0,0},{1,0,3}}];"
               "j=={{6,0,0},{0,12,1},{0,0,12}}]");
}
static void test_deficient_relation(void) {
    check_true("Module[{s,j,m},m={{-103,-191,-255},{110,190,222},{9,9,33}};"
               "{s,j}=JordanDecomposition[m];"
               "Simplify[m.s-s.j]==ConstantArray[0,{3,3}]]");
}
static void test_deficient_jordan_form(void) {
    check_true("Module[{s,j},{s,j}=JordanDecomposition[{{-103,-191,-255},{110,190,222},{9,9,33}}];"
               "j=={{24,0,0},{0,48,1},{0,0,48}}]");
}
static void test_deficient_superdiagonal(void) {
    /* The lone 1 marks the deficient eigenspace at 48. */
    check_true("Module[{s,j},{s,j}=JordanDecomposition[{{-103,-191,-255},{110,190,222},{9,9,33}}];"
               "Diagonal[j,1]=={0,1}]");
}
static void test_exact_4x4_single_chain(void) {
    /* A 4x4 with one size-3 block and one size-1 block for eigenvalue 1. */
    check_true("Module[{s,j,m},m={{1,0,0,0},{0,1,0,0},{1,-1,1,0},{1,-1,1,1}};"
               "{s,j}=JordanDecomposition[m];"
               "j=={{1,0,0,0},{0,1,1,0},{0,0,1,1},{0,0,0,1}} && "
               "Simplify[m.s-s.j]==ConstantArray[0,{4,4}]]");
}
static void test_symbolic_2x2(void) {
    check_true("Module[{s,j,m},m={{a,b},{c,d}};"
               "{s,j}=JordanDecomposition[m];"
               "Simplify[m.s-s.j]==ConstantArray[0,{2,2}]]");
}
static void test_diagonal_distinct(void) {
    check_fullform("JordanDecomposition[{{2,0},{0,3}}]",
                   "List[List[List[1, 0], List[0, 1]], List[List[2, 0], List[0, 3]]]");
}
static void test_1x1(void) {
    check_fullform("JordanDecomposition[{{5}}]",
                   "List[List[List[1]], List[List[5]]]");
}
static void test_nilpotent_2x2(void) {
    /* {{0,1},{0,0}} is a single 2x2 nilpotent Jordan block for eigenvalue 0. */
    check_true("Module[{s,j},{s,j}=JordanDecomposition[{{0,1},{0,0}}];"
               "j=={{0,1},{0,0}} && Diagonal[j,1]=={1}]");
}

/* ---- Numeric ---- */

static void test_numeric_real_residual(void) {
    check_true("Module[{s,j,m},m={{-1.2,2.7,3.8},{4.2,4.4,5.3},{3.5,7.6,6.8}};"
               "{s,j}=JordanDecomposition[m];"
               "Max[Abs[Flatten[N[m.s-s.j]]]]<10^-8]");
}
static void test_numeric_diagonalizable_form(void) {
    /* Distinct real spectrum -> j is diagonal (exact zeros off-diagonal). */
    check_true("Module[{s,j},{s,j}=JordanDecomposition[{{-1.2,2.7,3.8},{4.2,4.4,5.3},{3.5,7.6,6.8}}];"
               "Max[Abs[Flatten[j-DiagonalMatrix[Diagonal[j]]]]]<10^-12]");
}
static void test_numeric_complex_residual(void) {
    check_true("Module[{s,j,m},m={{3.14,0.3},{I,1+1.5 I}};"
               "{s,j}=JordanDecomposition[m];"
               "Max[Abs[Flatten[N[m.s-s.j]]]]<10^-8]");
}
static void test_mpfr_residual(void) {
    check_true("Module[{s,j,m},m=N[{{1,2,0},{0,3,1},{2,0,1}},20];"
               "{s,j}=JordanDecomposition[m];"
               "Max[Abs[Flatten[N[m.s-s.j]]]]<10^-15]");
}
static void test_numeric_defective_block(void) {
    /* N of an exactly defective matrix rationalizes back to a genuine block. */
    check_true("Module[{s,j,m},m=N[{{1,0,0,0},{0,1,0,0},{1,-1,1,0},{1,-1,1,1}}];"
               "{s,j}=JordanDecomposition[m];"
               "Max[Abs[Diagonal[j,1]-{0,1,1}]]<10^-6 && Max[Abs[Flatten[N[m.s-s.j]]]]<10^-6]");
}

/* ---- Surfaces ---- */

static void test_ndarray_input(void) {
    check_true("Module[{r,m},m={{2.,0.},{0.,3.}};"
               "r=JordanDecomposition[NDArray[m]];"
               "Max[Abs[Flatten[N[m.r[[1]]-r[[1]].r[[2]]]]]]<10^-8]");
}

/* ---- Errors / unevaluated ---- */

static void test_arg_error(void) {
    check_fullform("JordanDecomposition[]", "JordanDecomposition[]");
}
static void test_two_args(void) {
    check_fullform("JordanDecomposition[{{1,2},{3,4}}, x]",
                   "JordanDecomposition[List[List[1, 2], List[3, 4]], x]");
}
static void test_nonsquare(void) {
    check_fullform("JordanDecomposition[{{1,2,3},{4,5,6}}]",
                   "JordanDecomposition[List[List[1, 2, 3], List[4, 5, 6]]]");
}

static void test_docstring(void) {
    const char* doc = symtab_get_docstring("JordanDecomposition");
    ASSERT(doc != NULL);
    ASSERT(strstr(doc, "JordanDecomposition[m]") != NULL);
    printf("PASS: docstring present\n");
}

int main(void) {
    symtab_init();
    core_init();

    printf("Running JordanDecomposition tests...\n");
    TEST(test_exact_3x3_relation);
    TEST(test_exact_3x3_jordan_form);
    TEST(test_deficient_relation);
    TEST(test_deficient_jordan_form);
    TEST(test_deficient_superdiagonal);
    TEST(test_exact_4x4_single_chain);
    TEST(test_symbolic_2x2);
    TEST(test_diagonal_distinct);
    TEST(test_1x1);
    TEST(test_nilpotent_2x2);
    TEST(test_numeric_real_residual);
    TEST(test_numeric_diagonalizable_form);
    TEST(test_numeric_complex_residual);
    TEST(test_mpfr_residual);
    TEST(test_numeric_defective_block);
    TEST(test_ndarray_input);
    TEST(test_arg_error);
    TEST(test_two_args);
    TEST(test_nonsquare);
    TEST(test_docstring);

    printf("\nAll JordanDecomposition tests passed.\n");
    return 0;
}
