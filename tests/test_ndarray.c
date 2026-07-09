/* End-to-end tests for the NDArray[...] dense N-dimensional array type.
 *
 * Each assert_eval_eq drives the full pipeline: parse -> evaluate -> print,
 * so these exercise the parser, the evaluator/builtins, and the printer
 * together. */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include <stdio.h>

/* ---------- Construction, round-trip, printing ---------- */

void test_ndarray_construct_rank2() {
    assert_eval_eq("NDArray[{{1, 2}, {3, 4}}]",
                   "NDArray[{{1.0, 2.0}, {3.0, 4.0}}]", 0);
}

void test_ndarray_construct_rank1() {
    assert_eval_eq("NDArray[{1, 2, 3}]", "NDArray[{1.0, 2.0, 3.0}]", 0);
}

void test_ndarray_construct_rank3() {
    assert_eval_eq("NDArray[{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}}]",
                   "NDArray[{{{1.0, 2.0}, {3.0, 4.0}}, {{5.0, 6.0}, {7.0, 8.0}}}]", 0);
}

void test_ndarray_construct_rank4() {
    assert_eval_eq("Dimensions[NDArray[{{{{1,2,3}}}, {{{4,5,6}}}}]]", "{2, 1, 1, 3}", 0);
}

void test_ndarray_construct_reals() {
    assert_eval_eq("NDArray[{{1.5, 2.5}, {3.5, 4.5}}]",
                   "NDArray[{{1.5, 2.5}, {3.5, 4.5}}]", 0);
}

void test_ndarray_reject_jagged() {
    /* Ragged nested lists can't form a rectangular array: a NDArray::ragged
     * warning is printed (to stderr, not checked here) and the expression is
     * left unevaluated. Covers unequal sublist lengths and mixed scalar/list
     * siblings. */
    assert_eval_eq("NDArray[{{1, 2}, {3}}]", "NDArray[{{1, 2}, {3}}]", 0);
    assert_eval_eq("NDArray[{1, {2}}]", "NDArray[{1, {2}}]", 0);
    assert_eval_eq("NDArray[{1, {2, 3}}]", "NDArray[{1, {2, 3}}]", 0);
}

void test_ndarray_reject_symbolic() {
    /* A symbolic entry can't be packed to machine precision -- unevaluated. */
    assert_eval_eq("NDArray[{{1, x}, {3, 4}}]", "NDArray[{{1, x}, {3, 4}}]", 0);
}

void test_ndarray_reject_non_list() {
    assert_eval_eq("NDArray[5]", "NDArray[5]", 0);
}

void test_ndarray_normal_roundtrip() {
    /* Normal[NDArray[list]] gives back the original list (as reals). */
    assert_eval_eq("Normal[NDArray[{{1, 2}, {3, 4}}]]",
                   "{{1.0, 2.0}, {3.0, 4.0}}", 0);
    assert_eval_eq("Normal[NDArray[{1, 2, 3}]]", "{1.0, 2.0, 3.0}", 0);
    /* Rank-3 round-trip. */
    assert_eval_eq("Normal[NDArray[{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}}]]",
                   "{{{1.0, 2.0}, {3.0, 4.0}}, {{5.0, 6.0}, {7.0, 8.0}}}", 0);
}

/* ---------- Head / Dimensions / Depth / Length / predicates ---------- */

void test_ndarray_head() {
    assert_eval_eq("Head[NDArray[{{1, 2}, {3, 4}}]]", "NDArray", 0);
}

void test_ndarray_dimensions() {
    assert_eval_eq("Dimensions[NDArray[{{1, 2}, {3, 4}}]]", "{2, 2}", 0);
    assert_eval_eq("Dimensions[NDArray[{1, 2, 3}]]", "{3}", 0);
    assert_eval_eq("Dimensions[NDArray[{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}}]]",
                   "{2, 2, 2}", 0);
}

void test_ndarray_depth() {
    /* Depth of a rank-r NDArray matches the equivalent nested list: r + 1. */
    assert_eval_eq("Depth[NDArray[{1, 2, 3}]]", "2", 0);
    assert_eval_eq("Depth[NDArray[{{1, 2}, {3, 4}}]]", "3", 0);
    assert_eval_eq("Depth[NDArray[{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}}]]", "4", 0);
}

void test_ndarray_length() {
    /* Length is the leading-axis length (numpy len). */
    assert_eval_eq("Length[NDArray[{1, 2, 3}]]", "3", 0);
    assert_eval_eq("Length[NDArray[{{1, 2}, {3, 4}}]]", "2", 0);
    assert_eval_eq("Length[NDArray[{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}}]]", "2", 0);
}

void test_ndarray_ndarrayq() {
    assert_eval_eq("NDArrayQ[NDArray[{{1, 2}, {3, 4}}]]", "True", 0);
    assert_eval_eq("NDArrayQ[{{1, 2}, {3, 4}}]", "False", 0);
    assert_eval_eq("NDArrayQ[5]", "False", 0);
}

void test_ndarray_matrixq() {
    assert_eval_eq("MatrixQ[NDArray[{{1, 2}, {3, 4}}]]", "True", 0);
    assert_eval_eq("MatrixQ[NDArray[{1, 2, 3}]]", "False", 0);
}

void test_ndarray_vectorq() {
    assert_eval_eq("VectorQ[NDArray[{1, 2, 3}]]", "True", 0);
    assert_eval_eq("VectorQ[NDArray[{{1, 2}, {3, 4}}]]", "False", 0);
}

void test_ndarray_listq_always_false() {
    /* Visibly distinct from List: ListQ must never call an NDArray a List,
     * even though a plain nested-list equivalent answers True. */
    assert_eval_eq("ListQ[NDArray[{{1, 2}, {3, 4}}]]", "False", 0);
    assert_eval_eq("ListQ[{{1, 2}, {3, 4}}]", "True", 0);
}

/* ---------- Equality / structural comparison ---------- */

void test_ndarray_equality() {
    /* Structurally-equal NDArray values are decided True via expr_eq, same
     * as identical Lists. Unequal ones stay symbolically unevaluated --
     * matching e.g. {1, 2} == {1, 3}, which is also not a decidable
     * numeric comparison and is left as Equal[...]. */
    assert_eval_eq("NDArray[{{1, 2}, {3, 4}}] == NDArray[{{1, 2}, {3, 4}}]", "True", 0);
    assert_eval_eq("NDArray[{{1, 2}, {3, 4}}] == NDArray[{{1, 2}, {3, 5}}]",
                   "NDArray[{{1.0, 2.0}, {3.0, 4.0}}] == NDArray[{{1.0, 2.0}, {3.0, 5.0}}]", 0);
}

void test_ndarray_sameq() {
    assert_eval_eq("NDArray[{{1, 2}, {3, 4}}] === NDArray[{{1, 2}, {3, 4}}]", "True", 0);
    /* Not structurally identical to the equivalent nested List -- distinct types. */
    assert_eval_eq("NDArray[{{1, 2}, {3, 4}}] === {{1.0, 2.0}, {3.0, 4.0}}", "False", 0);
}

/* ---------- Dot fast path ---------- */

void test_ndarray_dot_matrix_matrix() {
    assert_eval_eq("Dot[NDArray[{{1, 2}, {3, 4}}], NDArray[{{5, 6}, {7, 8}}]]",
                   "NDArray[{{19.0, 22.0}, {43.0, 50.0}}]", 0);
}

void test_ndarray_dot_vector_vector() {
    assert_eval_eq("Dot[NDArray[{1, 2, 3}], NDArray[{4, 5, 6}]]", "32.0", 0);
}

void test_ndarray_dot_matrix_vector() {
    assert_eval_eq("Dot[NDArray[{{1, 2}, {3, 4}}], NDArray[{5, 6}]]",
                   "NDArray[{17.0, 39.0}]", 0);
}

void test_ndarray_dot_matches_list_path() {
    /* The NDArray fast path and the generic nested-list path must agree
     * (modulo Integer vs. machine-Real representation). */
    assert_eval_eq("Normal[Dot[NDArray[{{1, 2}, {3, 4}}], NDArray[{{5, 6}, {7, 8}}]]]",
                   "{{19.0, 22.0}, {43.0, 50.0}}", 0);
    assert_eval_eq("Dot[{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}]", "{{19, 22}, {43, 50}}", 0);
}

void test_ndarray_dot_shape_mismatch() {
    /* Incompatible inner dimensions: Dot stays unevaluated (error printed
     * to stderr, not checked here). */
    assert_eval_eq("Dot[NDArray[{{1, 2}, {3, 4}}], NDArray[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]]",
                   "Dot[NDArray[{{1.0, 2.0}, {3.0, 4.0}}], "
                   "NDArray[{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}}]]", 0);
}

/* ---------- Elementwise Plus / Times fast paths ---------- */

void test_ndarray_plus_elementwise() {
    assert_eval_eq("NDArray[{{1, 2}, {3, 4}}] + NDArray[{{5, 6}, {7, 8}}]",
                   "NDArray[{{6.0, 8.0}, {10.0, 12.0}}]", 0);
}

void test_ndarray_times_elementwise() {
    assert_eval_eq("NDArray[{{1, 2}, {3, 4}}] * NDArray[{{5, 6}, {7, 8}}]",
                   "NDArray[{{5.0, 12.0}, {21.0, 32.0}}]", 0);
}

void test_ndarray_plus_matches_list_path() {
    assert_eval_eq("{{1, 2}, {3, 4}} + {{5, 6}, {7, 8}}", "{{6, 8}, {10, 12}}", 0);
    assert_eval_eq("Normal[NDArray[{{1, 2}, {3, 4}}] + NDArray[{{5, 6}, {7, 8}}]]",
                   "{{6.0, 8.0}, {10.0, 12.0}}", 0);
}

void test_ndarray_times_matches_list_path() {
    assert_eval_eq("{{1, 2}, {3, 4}} * {{5, 6}, {7, 8}}", "{{5, 12}, {21, 32}}", 0);
    assert_eval_eq("Normal[NDArray[{{1, 2}, {3, 4}}] * NDArray[{{5, 6}, {7, 8}}]]",
                   "{{5.0, 12.0}, {21.0, 32.0}}", 0);
}

void test_ndarray_plus_shape_mismatch_degrades() {
    /* Mismatched shapes: falls through and stays symbolic/unevaluated. */
    assert_eval_eq("NDArray[{1, 2}] + NDArray[{1, 2, 3}]",
                   "NDArray[{1.0, 2.0}] + NDArray[{1.0, 2.0, 3.0}]", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_ndarray_construct_rank2);
    TEST(test_ndarray_construct_rank1);
    TEST(test_ndarray_construct_rank3);
    TEST(test_ndarray_construct_rank4);
    TEST(test_ndarray_construct_reals);
    TEST(test_ndarray_reject_jagged);
    TEST(test_ndarray_reject_symbolic);
    TEST(test_ndarray_reject_non_list);
    TEST(test_ndarray_normal_roundtrip);

    TEST(test_ndarray_head);
    TEST(test_ndarray_dimensions);
    TEST(test_ndarray_depth);
    TEST(test_ndarray_length);
    TEST(test_ndarray_ndarrayq);
    TEST(test_ndarray_matrixq);
    TEST(test_ndarray_vectorq);
    TEST(test_ndarray_listq_always_false);

    TEST(test_ndarray_equality);
    TEST(test_ndarray_sameq);

    TEST(test_ndarray_dot_matrix_matrix);
    TEST(test_ndarray_dot_vector_vector);
    TEST(test_ndarray_dot_matrix_vector);
    TEST(test_ndarray_dot_matches_list_path);
    TEST(test_ndarray_dot_shape_mismatch);

    TEST(test_ndarray_plus_elementwise);
    TEST(test_ndarray_times_elementwise);
    TEST(test_ndarray_plus_matches_list_path);
    TEST(test_ndarray_times_matches_list_path);
    TEST(test_ndarray_plus_shape_mismatch_degrades);

    printf("All NDArray tests passed.\n");
    return 0;
}
