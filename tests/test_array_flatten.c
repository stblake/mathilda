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
               "ArrayFlatten %s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(e);
    expr_free(r);
}

/* ---------- Block matrices ---------- */

static void test_af_block_repeat(void) {
    /* A 2x3 grid of the 2x2 block m -> a 4x6 matrix. */
    run_full("m = {{1, 2}, {3, 4}}; ArrayFlatten[{{m, m, m}, {m, m, m}}]",
             "List[List[1, 2, 1, 2, 1, 2], List[3, 4, 3, 4, 3, 4], "
             "List[1, 2, 1, 2, 1, 2], List[3, 4, 3, 4, 3, 4]]");
}

static void test_af_zero_blocks(void) {
    /* 0 is a scalar, replicated to a 2x2 zero block. */
    run_full("m = {{1, 2}, {3, 4}}; ArrayFlatten[{{0, 0, m}, {m, m, 0}}]",
             "List[List[0, 0, 0, 0, 1, 2], List[0, 0, 0, 0, 3, 4], "
             "List[1, 2, 1, 2, 0, 0], List[3, 4, 3, 4, 0, 0]]");
}

/* ---------- Scalar grids ---------- */

static void test_af_all_scalar_grid(void) {
    /* Every element is a scalar (depth 0 < 2): each fills a 1x1 block, so the
       grid passes through unchanged. */
    run_full("ArrayFlatten[{{1, 2}, {3, 4}}]",
             "List[List[1, 2], List[3, 4]]");
}

static void test_af_symbolic_scalar_replication(void) {
    /* A symbol s (depth 0) is replicated to fill the 2x2 block demanded by the
       m blocks sharing its row/column. */
    run_full("m = {{1, 2}, {3, 4}}; ArrayFlatten[{{s, m}, {m, s}}]",
             "List[List[s, s, 1, 2], List[s, s, 3, 4], "
             "List[1, 2, s, s], List[3, 4, s, s]]");
}

/* ---------- Rank / Dimensions (general r) ---------- */

static void test_af_dims_rank4(void) {
    /* ArrayFlatten[Array[a, {1,2,3,4}]] // Dimensions -> {3, 8}. */
    run_full("Dimensions[ArrayFlatten[Array[a, {1, 2, 3, 4}]]]",
             "List[3, 8]");
}

static void test_af_dims_rank6_default(void) {
    /* Default r=2 flattens the first two level-pairs; trailing axes carried
       along: {3, 8, 5, 6}. */
    run_full("Dimensions[ArrayFlatten[Array[a, {1, 2, 3, 4, 5, 6}]]]",
             "List[3, 8, 5, 6]");
}

static void test_af_dims_rank6_r3(void) {
    /* r=3 flattens three level-pairs of a rank-6 array -> rank 3: {4, 10, 18}. */
    run_full("Dimensions[ArrayFlatten[Array[a, {1, 2, 3, 4, 5, 6}], 3]]",
             "List[4, 10, 18]");
}

/* ---------- Attributes & documentation ---------- */

static void test_af_attributes_protected(void) {
    run_full("MemberQ[Attributes[ArrayFlatten], Protected]", "True");
}

static void test_af_docstring_present(void) {
    SymbolDef* def = symtab_get_def("ArrayFlatten");
    ASSERT_MSG(def != NULL && def->docstring != NULL && def->docstring[0] != '\0',
               "ArrayFlatten should have a non-empty docstring");
}

/* ---------- Unevaluated / error cases ---------- */

static void test_af_wrong_arity(void) {
    run_full("ArrayFlatten[]", "ArrayFlatten[]");
    run_full("ArrayFlatten[a, 2, 3]", "ArrayFlatten[a, 2, 3]");
}

static void test_af_non_array(void) {
    /* A non-grid first argument leaves the call unevaluated. */
    run_full("ArrayFlatten[5]", "ArrayFlatten[5]");
    run_full("ArrayFlatten[x]", "ArrayFlatten[x]");
}

static void test_af_bad_level(void) {
    /* Non-integer or < 1 level count -> unevaluated. */
    run_full("ArrayFlatten[{{1, 2}, {3, 4}}, 0]",
             "ArrayFlatten[List[List[1, 2], List[3, 4]], 0]");
    run_full("ArrayFlatten[{{1, 2}, {3, 4}}, x]",
             "ArrayFlatten[List[List[1, 2], List[3, 4]], x]");
}

static void test_af_mismatched_blocks(void) {
    /* Two blocks in the same grid row disagree on their first dimension
       (2 vs 3), so they do not fit and the call stays unevaluated. */
    run_full("ArrayFlatten[{{{{1, 2}, {3, 4}}, {{1, 2}, {3, 4}, {5, 6}}}}]",
             "ArrayFlatten[List[List[List[List[1, 2], List[3, 4]], "
             "List[List[1, 2], List[3, 4], List[5, 6]]]]]");
}

int main(void) {
    symtab_init();
    core_init();

    /* Block matrices */
    TEST(test_af_block_repeat);
    TEST(test_af_zero_blocks);

    /* Scalar grids */
    TEST(test_af_all_scalar_grid);
    TEST(test_af_symbolic_scalar_replication);

    /* Rank / Dimensions */
    TEST(test_af_dims_rank4);
    TEST(test_af_dims_rank6_default);
    TEST(test_af_dims_rank6_r3);

    /* Attributes & docs */
    TEST(test_af_attributes_protected);
    TEST(test_af_docstring_present);

    /* Unevaluated cases */
    TEST(test_af_wrong_arity);
    TEST(test_af_non_array);
    TEST(test_af_bad_level);
    TEST(test_af_mismatched_blocks);

    printf("All ArrayFlatten tests passed!\n");
    return 0;
}
