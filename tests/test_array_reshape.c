#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "ndarray.h"
#include "pack.h"
#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Compare an evaluated result (delisted so a packed buffer reads as a List) in
 * OutputForm against `expected`. */
static void run(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* r = test_delist(evaluate(e));
    char* s = expr_to_string(r);
    ASSERT_MSG(strcmp(s, expected) == 0,
               "ArrayReshape %s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

/* ---------- Core reshaping ---------- */

static void test_ar_matrix(void) {
    run("ArrayReshape[{a,b,c,d,e,f},{2,3}]", "{{a, b, c}, {d, e, f}}");
}
static void test_ar_depth3(void) {
    run("ArrayReshape[Range[24],{2,3,4}]",
        "{{{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}}, "
        "{{13, 14, 15, 16}, {17, 18, 19, 20}, {21, 22, 23, 24}}}");
}
static void test_ar_single_int_dims(void) {
    run("ArrayReshape[{1,2,3,4},3]", "{1, 2, 3}");
}
static void test_ar_nested_input(void) {
    /* The input is fully flattened first. */
    run("ArrayReshape[{{1,2,3},{4,5,6},{7,8,9},{10,11,12}},{3,4}]",
        "{{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}}");
    run("ArrayReshape[{{1,2,3},{4,5,6},{7,8,9},{10,11,12}},{2,2,3}]",
        "{{{1, 2, 3}, {4, 5, 6}}, {{7, 8, 9}, {10, 11, 12}}}");
}

/* ---------- Dropping extra elements ---------- */

static void test_ar_drop(void) {
    run("ArrayReshape[Range[100],{2,3,4}]",
        "{{{1, 2, 3, 4}, {5, 6, 7, 8}, {9, 10, 11, 12}}, "
        "{{13, 14, 15, 16}, {17, 18, 19, 20}, {21, 22, 23, 24}}}");
}

/* ---------- Padding a short input ---------- */

static void test_ar_pad_default_zero(void) {
    run("ArrayReshape[{1,2,3},{5}]", "{1, 2, 3, 0, 0}");
}
static void test_ar_pad_constant(void) {
    run("ArrayReshape[{1,2,3,4,5,6,7},{5,3},x]",
        "{{1, 2, 3}, {4, 5, 6}, {7, x, x}, {x, x, x}, {x, x, x}}");
}
static void test_ar_pad_cyclic(void) {
    run("ArrayReshape[{1,2,3},{6},{a,b}]", "{1, 2, 3, a, b, a}");
}
static void test_ar_pad_periodic(void) {
    run("ArrayReshape[{1,2,3},{7},\"Periodic\"]", "{1, 2, 3, 1, 2, 3, 1}");
}

/* ---------- Zero dimensions ---------- */

static void test_ar_zero_dim(void) {
    run("ArrayReshape[{1,2,3,4,5,6},{0,3}]", "{}");
}

/* ---------- Packed input ---------- */

static void test_ar_packed(void) {
    run("ArrayReshape[Range[6],{2,3}]", "{{1, 2, 3}, {4, 5, 6}}");
    run("ArrayReshape[Range[7],{3,2}]", "{{1, 2}, {3, 4}, {5, 6}}");
    /* packed with an exact-int pad on the tail */
    run("ArrayReshape[Range[3],{5}]", "{1, 2, 3, 0, 0}");
}
static void test_ar_packed_real(void) {
    run("ArrayReshape[Range[1.,6.],{2,3}]",
        "{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}}");
}

/* ---------- Attributes / documentation ---------- */

static void test_ar_protected(void) {
    run("MemberQ[Attributes[ArrayReshape], Protected]", "True");
}
static void test_ar_docstring(void) {
    SymbolDef* def = symtab_get_def("ArrayReshape");
    ASSERT_MSG(def && def->docstring && def->docstring[0],
               "ArrayReshape should have a non-empty docstring");
}

/* ---------- Unevaluated / error cases ---------- */

static void test_ar_arity(void) {
    run("ArrayReshape[]", "ArrayReshape[]");
    run("ArrayReshape[x]", "ArrayReshape[x]");
}
static void test_ar_nonlist(void) {
    run("ArrayReshape[5,{2,2}]", "ArrayReshape[5, {2, 2}]");
}
static void test_ar_bad_dims(void) {
    run("ArrayReshape[{1,2,3},{2,x}]", "ArrayReshape[{1, 2, 3}, {2, x}]");
    run("ArrayReshape[{1,2,3},-1]", "ArrayReshape[{1, 2, 3}, -1]");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_ar_matrix);
    TEST(test_ar_depth3);
    TEST(test_ar_single_int_dims);
    TEST(test_ar_nested_input);
    TEST(test_ar_drop);
    TEST(test_ar_pad_default_zero);
    TEST(test_ar_pad_constant);
    TEST(test_ar_pad_cyclic);
    TEST(test_ar_pad_periodic);
    TEST(test_ar_zero_dim);
    TEST(test_ar_packed);
    TEST(test_ar_packed_real);
    TEST(test_ar_protected);
    TEST(test_ar_docstring);
    TEST(test_ar_arity);
    TEST(test_ar_nonlist);
    TEST(test_ar_bad_dims);

    printf("All ArrayReshape tests passed!\n");
    return 0;
}
