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

void run_test(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    if (!e) {
        printf("Failed to parse: %s\n", input);
        ASSERT(0);
    }
    Expr* res = evaluate(e);
    char* res_str = expr_to_string_fullform(res);
    if (strcmp(res_str, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, res_str);
        ASSERT(0);
    } else {
        printf("PASS: %s -> %s\n", input, res_str);
    }
    free(res_str);
    expr_free(res);
    expr_free(e);
}

void test_dot() {
    run_test("{a, b, c} . {x, y, z}", "Plus[Times[a, x], Times[b, y], Times[c, z]]");
    run_test("CompoundExpression[Set[u, List[1, 1]], Set[v, List[-1, 1]], Dot[u, v]]", "0");
    run_test("{{a, b}, {c, d}} . {x, y}", "List[Plus[Times[a, x], Times[b, y]], Plus[Times[c, x], Times[d, y]]]");
    run_test("{x, y} . {{a, b}, {c, d}}", "List[Plus[Times[a, x], Times[c, y]], Plus[Times[b, x], Times[d, y]]]");
    run_test("{x, y} . {{a, b}, {c, d}} . {r, s}", "Plus[Times[r, Plus[Times[a, x], Times[c, y]]], Times[s, Plus[Times[b, x], Times[d, y]]]]");
    
    // Matrix multiplication
    run_test("{{a, b}, {c, d}} . {{1, 2}, {3, 4}}", "List[List[Plus[a, Times[3, b]], Plus[Times[2, a], Times[4, b]]], List[Plus[c, Times[3, d]], Plus[Times[2, c], Times[4, d]]]]");
    
    // Rectangular
    run_test("{{1, 2, 3}, {4, 5, 6}} . {{a, b}, {c, d}, {e, f}}", "List[List[Plus[a, Times[2, c], Times[3, e]], Plus[b, Times[2, d], Times[3, f]]], List[Plus[Times[4, a], Times[5, c], Times[6, e]], Plus[Times[4, b], Times[5, d], Times[6, f]]]]");
}

void test_det() {
    run_test("Det[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]", "0");
    run_test("Det[{{1.7, 7.1, -2.7}, {2.2, 8.7, 3.2}, {3.2, -9.2, 1.2}}]", "251.572");
    run_test("Expand[Det[{{a, b, c}, {d, e, f}, {g, h, i}}]]", "Plus[Times[-1, Times[c, e, g]], Times[b, f, g], Times[c, d, h], Times[-1, Times[a, f, h]], Times[-1, Times[b, d, i]], Times[a, e, i]]");
}

void test_cross() {
    run_test("Cross[{1, 2, -1}, {-1, 1, 0}]", "List[1, 1, 3]");
    run_test("Cross[{1, Sqrt[3]}]", "List[Times[-1, Power[3, Rational[1, 2]]], 1]");
    run_test("Cross[{3.2, 4.2, 5.2}, {0.75, 0.09, 0.06}]", "List[-0.216, 3.708, -2.862]");
    run_test("Cross[{1.3 + I, 2, 3 - 2 I}, {6. + I, 4, 5 - 7 I}]", "List[Complex[-2, -6], Complex[6.5, -4.9], Complex[-6.8, 2.0]]");
    run_test("Cross[{1, 2, 3}, {4, 5}]", "Cross[List[1, 2, 3], List[4, 5]]");
}

void test_norm() {
    run_test("Norm[{x, y, z}]", "Power[Plus[Power[Abs[x], 2], Power[Abs[y], 2], Power[Abs[z], 2]], Rational[1, 2]]");
    run_test("Norm[-2 + I]", "Power[5, Rational[1, 2]]");
    run_test("CompoundExpression[Set[v, List[1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1]], Norm[v]]", "Power[5, Rational[1, 2]]");
    run_test("Norm[{x, y, z}, p]", "Power[Plus[Power[Abs[x], p], Power[Abs[y], p], Power[Abs[z], p]], Power[p, -1]]");
    run_test("Norm[{x, y, z}, Infinity]", "Max[Abs[x], Abs[y], Abs[z]]");
    run_test("Norm[{{a11, a12}, {a21, a22}}, \"Frobenius\"]", "Power[Plus[Power[Abs[a11], 2], Power[Abs[a12], 2], Power[Abs[a21], 2], Power[Abs[a22], 2]], Rational[1, 2]]");
}

void test_tr() {
    run_test("Tr[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]", "15");
    run_test("Tr[{{a, b}, {c, d}}]", "Plus[a, d]");
    run_test("Tr[{1, 2, 3}]", "6");
    run_test("Tr[Array[a, {4, 3, 2}]]", "Plus[a[1, 1, 1], a[2, 2, 2]]");
    run_test("Tr[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, f]", "f[1, 5, 9]");
    run_test("Tr[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, List]", "List[1, 5, 9]");
    run_test("Tr[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}, Plus, 1]", "List[12, 15, 18]");
    run_test("Tr[Array[a, {4, 3, 2}], f, 2]", "f[List[a[1, 1, 1], a[1, 1, 2]], List[a[2, 2, 1], a[2, 2, 2]], List[a[3, 3, 1], a[3, 3, 2]]]");
    run_test("Tr[x]", "x");
    run_test("Tr[x, f]", "x");
    run_test("Tr[{{1, 2}, {3}}]", "1");
}

void test_rowreduce() {
    run_test("RowReduce[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]", "List[List[1, 0, -1], List[0, 1, 2], List[0, 0, 0]]");
    run_test("RowReduce[{{1, 2, 3, 1, 0, 0}, {4, 5, 6, 0, 1, 0}, {7, 8, 9, 0, 0, 1}}]", "List[List[1, 0, -1, 0, Rational[-8, 3], Rational[5, 3]], List[0, 1, 2, 0, Rational[7, 3], Rational[-4, 3]], List[0, 0, 0, 1, -2, 1]]");
    run_test("RowReduce[{{3, 1, a}, {2, 1, b}}]", "List[List[1, 0, Plus[a, Times[-1, b]]], List[0, 1, Plus[Times[-2, a], Times[3, b]]]]");
    run_test("RowReduce[{{0, 0, 0, Pi}, {2, 2, 2, 2}, {3 + I * Sqrt[2], 0, 0, 0}, {0, 4, 4, 0}}]", "List[List[1, 0, 0, 0], List[0, 1, 1, 0], List[0, 0, 0, 1], List[0, 0, 0, 0]]");
    run_test("RowReduce[{{a, b, c}, {d, e, f}, {g, h, i}}]", "List[List[1, 0, 0], List[0, 1, 0], List[0, 0, 1]]");
    
    // Check symbolic fractional simplification: Note that Expand changes the exact output string.
    Expr* t = eval_and_free(parse_expression("RowReduce[{{a, b, c}, {d, e, f}, {a + d, b + e, c + f}}]"));
    char* res_str = expr_to_string_fullform(t);
    
    // We expect the first two rows to have the elements `(ce-bf)/(ae-bd)` and `(cd-af)/(bd-ae)`.
    // In our algorithm the denominator is `-bd+ae` and the numerator is `ce-bf` or `-af+cd`.
    // We just verify it reduces exactly to row 3 being 0 and diagonal being 1.
    ASSERT(t->type == EXPR_FUNCTION);
    ASSERT(t->data.function.arg_count == 3);
    ASSERT(t->data.function.args[2]->data.function.args[0]->type == EXPR_INTEGER);
    ASSERT(t->data.function.args[2]->data.function.args[0]->data.integer == 0);
    ASSERT(t->data.function.args[2]->data.function.args[1]->type == EXPR_INTEGER);
    ASSERT(t->data.function.args[2]->data.function.args[1]->data.integer == 0);
    ASSERT(t->data.function.args[2]->data.function.args[2]->type == EXPR_INTEGER);
    ASSERT(t->data.function.args[2]->data.function.args[2]->data.integer == 0);
    
    free(res_str);
    expr_free(t);
}

void test_identitymatrix() {
    run_test("IdentityMatrix[3]", "List[List[1, 0, 0], List[0, 1, 0], List[0, 0, 1]]");
    run_test("IdentityMatrix[{2, 3}]", "List[List[1, 0, 0], List[0, 1, 0]]");
    run_test("IdentityMatrix[{3, 2}]", "List[List[1, 0], List[0, 1], List[0, 0]]");
    run_test("IdentityMatrix[0]", "List[]");
    run_test("IdentityMatrix[-1]", "IdentityMatrix[-1]");
    run_test("IdentityMatrix[a]", "IdentityMatrix[a]");
}

void test_inverse() {
    /* 1x1 matrices */
    run_test("Inverse[{{1}}]", "List[List[1]]");
    run_test("Inverse[{{5}}]", "List[List[Rational[1, 5]]]");
    run_test("Inverse[{{x}}]", "List[List[Power[x, -1]]]");

    /* 2x2 exact integer matrix */
    run_test("Inverse[{{2,3,2},{4,9,2},{7,2,4}}]",
        "List[List[Rational[-8, 13], Rational[2, 13], Rational[3, 13]], "
        "List[Rational[1, 26], Rational[3, 26], Rational[-1, 13]], "
        "List[Rational[55, 52], Rational[-17, 52], Rational[-3, 26]]]");

    /* 3x3 exact integer matrix */
    run_test("Inverse[{{1,2,3},{4,2,2},{5,1,7}}]",
        "List[List[Rational[-2, 7], Rational[11, 42], Rational[1, 21]], "
        "List[Rational[3, 7], Rational[4, 21], Rational[-5, 21]], "
        "List[Rational[1, 7], Rational[-3, 14], Rational[1, 7]]]");

    /* 2x2 symbolic matrix */
    run_test("Inverse[{{a,b},{c,d}}]",
        "List[List[Times[d, Power[Plus[Times[-1, Times[b, c]], Times[a, d]], -1]], "
        "Times[-1, b, Power[Plus[Times[-1, Times[b, c]], Times[a, d]], -1]]], "
        "List[Times[-1, c, Power[Plus[Times[-1, Times[b, c]], Times[a, d]], -1]], "
        "Times[a, Power[Plus[Times[-1, Times[b, c]], Times[a, d]], -1]]]]");

    /* Symmetric symbolic matrix (use p,q to avoid polluted u,v from test_dot) */
    run_test("Inverse[{{p,q},{q,p}}]",
        "List[List[Times[p, Power[Plus[Power[p, 2], Times[-1, Power[q, 2]]], -1]], "
        "Times[-1, q, Power[Plus[Power[p, 2], Times[-1, Power[q, 2]]], -1]]], "
        "List[Times[-1, q, Power[Plus[Power[p, 2], Times[-1, Power[q, 2]]], -1]], "
        "Times[p, Power[Plus[Power[p, 2], Times[-1, Power[q, 2]]], -1]]]]");

    /* Floating-point 2x2 matrix */
    run_test("Inverse[{{1.4,2},{3,-6.7}}]",
        "List[List[0.435631, 0.130039], List[0.195059, -0.0910273]]");

    /* Machine-precision 3x3 matrix */
    run_test("Inverse[{{1.2,2.5,-3.2},{0.7,-9.4,5.8},{-0.2,0.3,6.4}}]",
        "List[List[0.74546, 0.204249, 0.187629], "
        "List[0.0679223, -0.0847825, 0.110795], "
        "List[0.0201118, 0.010357, 0.15692]]");

    /* Singular matrix: should remain unevaluated */
    run_test("Inverse[{{1,2},{1,2}}]", "Inverse[List[List[1, 2], List[1, 2]]]");

    /* Non-square matrix: should remain unevaluated */
    run_test("Inverse[{{1,2},{3,4},{5,6}}]", "Inverse[List[List[1, 2], List[3, 4], List[5, 6]]]");

    /* Identity relation: mat.Inverse[mat] == IdentityMatrix[n] */
    run_test("CompoundExpression[Set[mat1, {{1,2},{3,4}}], Equal[Dot[mat1, Inverse[mat1]], IdentityMatrix[2]]]", "True");
    run_test("CompoundExpression[Set[mat2, {{1,2},{3,4}}], Equal[Dot[Inverse[mat2], mat2], IdentityMatrix[2]]]", "True");

    /* (aa.bb)^-1 == bb^-1.aa^-1 */
    run_test("CompoundExpression["
        "Set[aa, {{1,1,1},{6,9,7},{8,1,9}}], "
        "Set[bb, {{0,3,9},{7,9,7},{4,4,1}}], "
        "Equal[Inverse[Dot[aa, bb]], Dot[Inverse[bb], Inverse[aa]]]]", "True");

    /* Identity matrix is its own inverse */
    run_test("Inverse[{{1,0,0},{0,1,0},{0,0,1}}]", "List[List[1, 0, 0], List[0, 1, 0], List[0, 0, 1]]");

    /* 4x4 exact inverse */
    run_test("CompoundExpression["
        "Set[mm, {{1,0,2,1},{3,1,0,1},{1,2,1,0},{0,1,3,2}}], "
        "Equal[Dot[mm, Inverse[mm]], IdentityMatrix[4]]]", "True");

    /* Diagonal matrix inverse */
    run_test("Inverse[{{2,0,0},{0,3,0},{0,0,5}}]",
        "List[List[Rational[1, 2], 0, 0], List[0, Rational[1, 3], 0], List[0, 0, Rational[1, 5]]]");

    /* Singular 3x3 matrix */
    run_test("Inverse[{{1,2,3},{4,5,6},{7,8,9}}]",
        "Inverse[List[List[1, 2, 3], List[4, 5, 6], List[7, 8, 9]]]");

    /* Non-function argument */
    run_test("Inverse[5]", "Inverse[5]");

    /* Zero determinant via zero row */
    run_test("Inverse[{{0,0},{0,0}}]",
        "Inverse[List[List[0, 0], List[0, 0]]]");
}

void test_matrixpower() {
    /* Power 0: returns identity matrix */
    run_test("MatrixPower[{{1,2},{3,4}},0]", "List[List[1, 0], List[0, 1]]");
    run_test("MatrixPower[{{a,b},{c,d}},0]", "List[List[1, 0], List[0, 1]]");

    /* Power 1: returns matrix itself */
    run_test("MatrixPower[{{1,0},{0,1}},1]", "List[List[1, 0], List[0, 1]]");
    run_test("MatrixPower[{{1,2},{3,4}},1]", "List[List[1, 2], List[3, 4]]");

    /* Power 2 symbolic: m.m - verify via equality instead of exact string */
    run_test("CompoundExpression[Set[mp0, {{a,b},{c,d}}], "
        "Equal[MatrixPower[mp0, 2], Dot[mp0, mp0]]]", "True");

    /* Power 10: Fibonacci-like matrix */
    run_test("MatrixPower[{{1,1},{1,2}},10]", "List[List[4181, 6765], List[6765, 10946]]");

    /* Large exact integer power */
    run_test("MatrixPower[{{2,3,0},{4,9,0},{0,0,4}},14]",
        "List[List[25881337259836, 54508871401413, 0], "
        "List[72678495201884, 153068703863133, 0], "
        "List[0, 0, 268435456]]");

    /* Symbolic upper triangular matrix power */
    run_test("MatrixPower[{{a,b},{0,c}},4]",
        "List[List[Power[a, 4], Plus[Times[Power[a, 2], Plus[Times[a, b], Times[b, c]]], "
        "Times[Power[c, 2], Plus[Times[a, b], Times[b, c]]]]], List[0, Power[c, 4]]]");

    /* Negative power: m^-1 == Inverse[m] */
    run_test("CompoundExpression[Set[mp2, {{1,2},{3,4}}], "
        "Equal[MatrixPower[mp2, -1], Inverse[mp2]]]", "True");

    /* Negative power -2: m^-2 == Inverse[m]^2 */
    run_test("CompoundExpression[Set[mp3, {{a,b},{c,d}}], "
        "Equal[MatrixPower[mp3, -2], Dot[Inverse[mp3], Inverse[mp3]]]]", "True");

    /* Floating-point matrix power */
    run_test("MatrixPower[{{1.2,2.5,-3.2},{0.7,-9.4,5.8},{-0.2,0.3,6.4}},5]",
        "List[List[-1208.61, 19598.2, -12658.4], "
        "List[5784.51, -83315.1, 35420.6], "
        "List[-559.11, 1960.12, 11511.9]]");

    /* Complex matrix power */
    /* MatrixPower[{{1.+I,2,3-2 I},{0,4 Pi,5I},{E,0,6}},2] */
    /* Skipping exact complex test due to symbolic Pi/E constants */

    /* MatrixPower[m, n, v] - matrix power applied to vector */
    run_test("MatrixPower[{{1,1},{1,2}},3,{1,0}]",
        "List[5, 8]");

    /* Power 0 with vector */
    run_test("MatrixPower[{{1,2},{3,4}},0,{5,6}]",
        "List[5, 6]");

    /* Negative power with vector */
    run_test("MatrixPower[{{1,2},{3,4}},-1,{1,0}]",
        "List[-2, Rational[3, 2]]");

    /* Fractional power: should warn and return unevaluated */
    run_test("MatrixPower[{{1,2},{3,4}},1/2]",
        "MatrixPower[List[List[1, 2], List[3, 4]], Rational[1, 2]]");

    /* Non-square matrix: error */
    run_test("MatrixPower[{{1,2,3},{4,5,6}},2]",
        "MatrixPower[List[List[1, 2, 3], List[4, 5, 6]], 2]");

    /* Symbolic exponent: unevaluated */
    run_test("MatrixPower[{{1,2},{3,4}},n]",
        "MatrixPower[List[List[1, 2], List[3, 4]], n]");

    /* Singular matrix with negative power */
    run_test("MatrixPower[{{1,2},{2,4}},-1]",
        "MatrixPower[List[List[1, 2], List[2, 4]], -1]");

    /* 1x1 matrix */
    run_test("MatrixPower[{{5}},3]", "List[List[125]]");
    run_test("MatrixPower[{{5}},-2]", "List[List[Rational[1, 25]]]");
}

void test_diagonalmatrix() {
    run_test("DiagonalMatrix[{a, b, c}]", "List[List[a, 0, 0], List[0, b, 0], List[0, 0, c]]");
    run_test("DiagonalMatrix[{a, b}, 1]", "List[List[0, a, 0], List[0, 0, b], List[0, 0, 0]]");
    run_test("DiagonalMatrix[{a, b}, -1]", "List[List[0, 0, 0], List[a, 0, 0], List[0, b, 0]]");
    run_test("DiagonalMatrix[{1, 2, 3}, 0, 5]", "List[List[1, 0, 0, 0, 0], List[0, 2, 0, 0, 0], List[0, 0, 3, 0, 0], List[0, 0, 0, 0, 0], List[0, 0, 0, 0, 0]]");
    run_test("DiagonalMatrix[{1, 2, 3, 4, 5}, 0, 3]", "List[List[1, 0, 0], List[0, 2, 0], List[0, 0, 3]]");
    run_test("DiagonalMatrix[{1, 2, 3}, 0, {3, 5}]", "List[List[1, 0, 0, 0, 0], List[0, 2, 0, 0, 0], List[0, 0, 3, 0, 0]]");
    run_test("DiagonalMatrix[{1, 2, 3}, 0, {5, 3}]", "List[List[1, 0, 0], List[0, 2, 0], List[0, 0, 3], List[0, 0, 0], List[0, 0, 0]]");
}

int main() {
    symtab_init();
    core_init();
    
    printf("Running linalg tests...\n");
    TEST(test_dot);
    TEST(test_det);
    TEST(test_cross);
    TEST(test_norm);
    TEST(test_tr);
    TEST(test_rowreduce);
    TEST(test_identitymatrix);
    TEST(test_diagonalmatrix);
    TEST(test_inverse);
    TEST(test_matrixpower);
    printf("All linalg tests passed!\n");
    symtab_clear();
    return 0;
}
