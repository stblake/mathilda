#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void run_full(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* r = evaluate(e);
    char* s = expr_to_string_fullform(r);
    ASSERT_MSG(strcmp(s, expected) == 0,
               "ConstantArray %s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

/* ---------- Flat vectors ---------- */

static void test_ca_vector_symbol(void) {
    run_full("ConstantArray[c, 10]",
             "List[c, c, c, c, c, c, c, c, c, c]");
}

static void test_ca_vector_exact_zero(void) {
    run_full("ConstantArray[0, 10]",
             "List[0, 0, 0, 0, 0, 0, 0, 0, 0, 0]");
}

static void test_ca_vector_machine_zero(void) {
    run_full("ConstantArray[0., 10]",
             "List[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]");
}

static void test_ca_vector_exact_one(void) {
    run_full("ConstantArray[1, 10]",
             "List[1, 1, 1, 1, 1, 1, 1, 1, 1, 1]");
}

static void test_ca_vector_machine_one(void) {
    run_full("ConstantArray[1., 10]",
             "List[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]");
}

static void test_ca_vector_short(void) {
    run_full("ConstantArray[x, 1]", "List[x]");
    run_full("ConstantArray[x, 3]", "List[x, x, x]");
}

/* ---------- Matrices / nested arrays ---------- */

static void test_ca_matrix_symbol(void) {
    run_full("ConstantArray[c, {3, 4}]",
             "List[List[c, c, c, c], List[c, c, c, c], List[c, c, c, c]]");
}

static void test_ca_matrix_exact_zero(void) {
    run_full("ConstantArray[0, {3, 3}]",
             "List[List[0, 0, 0], List[0, 0, 0], List[0, 0, 0]]");
}

static void test_ca_matrix_machine_zero(void) {
    run_full("ConstantArray[0., {3, 3}]",
             "List[List[0.0, 0.0, 0.0], List[0.0, 0.0, 0.0], List[0.0, 0.0, 0.0]]");
}

static void test_ca_deeply_nested(void) {
    /* ConstantArray[x, {2,1,2,1,2}]
       -> {{{{{x,x}},{{x,x}}}},{{{{x,x}},{{x,x}}}}} */
    run_full("ConstantArray[x, {2, 1, 2, 1, 2}]",
             "List[List[List[List[List[x, x]], List[List[x, x]]]], "
             "List[List[List[List[x, x]], List[List[x, x]]]]]");
}

static void test_ca_compound_element(void) {
    /* The constant element is copied verbatim, including compound structure. */
    run_full("ConstantArray[{{1, 2}, {3, 4}}, {2, 2}]",
             "List["
             "List[List[List[1, 2], List[3, 4]], List[List[1, 2], List[3, 4]]], "
             "List[List[List[1, 2], List[3, 4]], List[List[1, 2], List[3, 4]]]]");
}

/* ---------- Zero dimensions ---------- */

static void test_ca_zero_dim_flat(void) {
    run_full("ConstantArray[c, 0]", "List[]");
}

static void test_ca_zero_dim_nested(void) {
    /* A zero at any level truncates that level to an empty list. */
    run_full("ConstantArray[c, {0, 3}]", "List[]");
    run_full("ConstantArray[c, {2, 0}]", "List[List[], List[]]");
}

/* ---------- Attributes & documentation ---------- */

static void test_ca_attributes_protected(void) {
    run_full("MemberQ[Attributes[ConstantArray], Protected]", "True");
}

static void test_ca_docstring_present(void) {
    SymbolDef* def = symtab_get_def("ConstantArray");
    ASSERT_MSG(def != NULL && def->docstring != NULL && def->docstring[0] != '\0',
               "ConstantArray should have a non-empty docstring");
}

/* ---------- Unevaluated / error cases ---------- */

static void test_ca_wrong_arity(void) {
    /* 0 or 1 args -> arity error, stays unevaluated. */
    run_full("ConstantArray[]", "ConstantArray[]");
    run_full("ConstantArray[c]", "ConstantArray[c]");
}

static void test_ca_bad_dims(void) {
    /* Non-integer, symbolic, negative, or empty dims -> unevaluated. */
    run_full("ConstantArray[c, x]", "ConstantArray[c, x]");
    run_full("ConstantArray[c, -1]", "ConstantArray[c, -1]");
    run_full("ConstantArray[c, 2.5]", "ConstantArray[c, 2.5]");
    run_full("ConstantArray[c, {2, x}]", "ConstantArray[c, List[2, x]]");
    run_full("ConstantArray[c, {}]", "ConstantArray[c, List[]]");
}

int main(void) {
    symtab_init();
    core_init();

    /* Vectors */
    TEST(test_ca_vector_symbol);
    TEST(test_ca_vector_exact_zero);
    TEST(test_ca_vector_machine_zero);
    TEST(test_ca_vector_exact_one);
    TEST(test_ca_vector_machine_one);
    TEST(test_ca_vector_short);

    /* Matrices / nested */
    TEST(test_ca_matrix_symbol);
    TEST(test_ca_matrix_exact_zero);
    TEST(test_ca_matrix_machine_zero);
    TEST(test_ca_deeply_nested);
    TEST(test_ca_compound_element);

    /* Zero dimensions */
    TEST(test_ca_zero_dim_flat);
    TEST(test_ca_zero_dim_nested);

    /* Attributes & docs */
    TEST(test_ca_attributes_protected);
    TEST(test_ca_docstring_present);

    /* Unevaluated cases */
    TEST(test_ca_wrong_arity);
    TEST(test_ca_bad_dims);

    printf("All ConstantArray tests passed!\n");
    return 0;
}
