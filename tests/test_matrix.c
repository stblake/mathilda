/* End-to-end tests for the Matrix[...] dense ndarray type.
 *
 * Each assert_eval_eq drives the full pipeline: parse -> evaluate -> print,
 * so these exercise the parser, the evaluator/builtins, and the printer
 * together. */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"
#include <stdio.h>

/* ---------- Construction, round-trip, printing ---------- */

void test_matrix_construct_rank2() {
    assert_eval_eq("Matrix[{{1, 2}, {3, 4}}]",
                   "Matrix[{{1.0, 2.0}, {3.0, 4.0}}]", 0);
}

void test_matrix_construct_rank1() {
    assert_eval_eq("Matrix[{1, 2, 3}]", "Matrix[{1.0, 2.0, 3.0}]", 0);
}

void test_matrix_construct_rank3() {
    assert_eval_eq("Matrix[{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}}]",
                   "Matrix[{{{1.0, 2.0}, {3.0, 4.0}}, {{5.0, 6.0}, {7.0, 8.0}}}]", 0);
}

void test_matrix_construct_reals() {
    assert_eval_eq("Matrix[{{1.5, 2.5}, {3.5, 4.5}}]",
                   "Matrix[{{1.5, 2.5}, {3.5, 4.5}}]", 0);
}

void test_matrix_reject_jagged() {
    /* Jagged nested list can't be packed -- stays unevaluated as Matrix[...]. */
    assert_eval_eq("Matrix[{{1, 2}, {3}}]", "Matrix[{{1, 2}, {3}}]", 0);
}

void test_matrix_reject_symbolic() {
    /* A symbolic entry can't be packed to machine precision -- unevaluated. */
    assert_eval_eq("Matrix[{{1, x}, {3, 4}}]", "Matrix[{{1, x}, {3, 4}}]", 0);
}

void test_matrix_reject_non_list() {
    assert_eval_eq("Matrix[5]", "Matrix[5]", 0);
}

void test_matrix_normal_roundtrip() {
    /* Normal[Matrix[list]] gives back the original list (as reals). */
    assert_eval_eq("Normal[Matrix[{{1, 2}, {3, 4}}]]",
                   "{{1.0, 2.0}, {3.0, 4.0}}", 0);
    assert_eval_eq("Normal[Matrix[{1, 2, 3}]]", "{1.0, 2.0, 3.0}", 0);
}

/* ---------- Head / Dimensions / predicates ---------- */

void test_matrix_head() {
    assert_eval_eq("Head[Matrix[{{1, 2}, {3, 4}}]]", "Matrix", 0);
}

void test_matrix_dimensions() {
    assert_eval_eq("Dimensions[Matrix[{{1, 2}, {3, 4}}]]", "{2, 2}", 0);
    assert_eval_eq("Dimensions[Matrix[{1, 2, 3}]]", "{3}", 0);
    assert_eval_eq("Dimensions[Matrix[{{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}}]]",
                   "{2, 2, 2}", 0);
}

void test_matrix_matrixq() {
    assert_eval_eq("MatrixQ[Matrix[{{1, 2}, {3, 4}}]]", "True", 0);
    assert_eval_eq("MatrixQ[Matrix[{1, 2, 3}]]", "False", 0);
}

void test_matrix_vectorq() {
    assert_eval_eq("VectorQ[Matrix[{1, 2, 3}]]", "True", 0);
    assert_eval_eq("VectorQ[Matrix[{{1, 2}, {3, 4}}]]", "False", 0);
}

void test_matrix_listq_always_false() {
    /* Visibly distinct from List: ListQ must never call a Matrix a List,
     * even though a plain nested-list equivalent answers True. */
    assert_eval_eq("ListQ[Matrix[{{1, 2}, {3, 4}}]]", "False", 0);
    assert_eval_eq("ListQ[{{1, 2}, {3, 4}}]", "True", 0);
}

/* ---------- Equality / structural comparison ---------- */

void test_matrix_equality() {
    /* Structurally-equal Matrix values are decided True via expr_eq, same
     * as identical Lists. Unequal ones stay symbolically unevaluated --
     * matching e.g. {1, 2} == {1, 3}, which is also not a decidable
     * numeric comparison and is left as Equal[...]. */
    assert_eval_eq("Matrix[{{1, 2}, {3, 4}}] == Matrix[{{1, 2}, {3, 4}}]", "True", 0);
    assert_eval_eq("Matrix[{{1, 2}, {3, 4}}] == Matrix[{{1, 2}, {3, 5}}]",
                   "Matrix[{{1.0, 2.0}, {3.0, 4.0}}] == Matrix[{{1.0, 2.0}, {3.0, 5.0}}]", 0);
}

void test_matrix_sameq() {
    assert_eval_eq("Matrix[{{1, 2}, {3, 4}}] === Matrix[{{1, 2}, {3, 4}}]", "True", 0);
    /* Not structurally identical to the equivalent nested List -- distinct types. */
    assert_eval_eq("Matrix[{{1, 2}, {3, 4}}] === {{1.0, 2.0}, {3.0, 4.0}}", "False", 0);
}

/* ---------- Dot fast path ---------- */

void test_matrix_dot_matrix_matrix() {
    assert_eval_eq("Dot[Matrix[{{1, 2}, {3, 4}}], Matrix[{{5, 6}, {7, 8}}]]",
                   "Matrix[{{19.0, 22.0}, {43.0, 50.0}}]", 0);
}

void test_matrix_dot_vector_vector() {
    assert_eval_eq("Dot[Matrix[{1, 2, 3}], Matrix[{4, 5, 6}]]", "32.0", 0);
}

void test_matrix_dot_matrix_vector() {
    assert_eval_eq("Dot[Matrix[{{1, 2}, {3, 4}}], Matrix[{5, 6}]]",
                   "Matrix[{17.0, 39.0}]", 0);
}

void test_matrix_dot_matches_list_path() {
    /* The Matrix fast path and the generic nested-list path must agree
     * (modulo Integer vs. machine-Real representation). */
    assert_eval_eq("Normal[Dot[Matrix[{{1, 2}, {3, 4}}], Matrix[{{5, 6}, {7, 8}}]]]",
                   "{{19.0, 22.0}, {43.0, 50.0}}", 0);
    assert_eval_eq("Dot[{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}]", "{{19, 22}, {43, 50}}", 0);
}

void test_matrix_dot_shape_mismatch() {
    /* Incompatible inner dimensions: Dot stays unevaluated (error printed
     * to stderr, not checked here). */
    assert_eval_eq("Dot[Matrix[{{1, 2}, {3, 4}}], Matrix[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]]",
                   "Dot[Matrix[{{1.0, 2.0}, {3.0, 4.0}}], "
                   "Matrix[{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}}]]", 0);
}

/* ---------- Elementwise Plus / Times fast paths ---------- */

void test_matrix_plus_elementwise() {
    assert_eval_eq("Matrix[{{1, 2}, {3, 4}}] + Matrix[{{5, 6}, {7, 8}}]",
                   "Matrix[{{6.0, 8.0}, {10.0, 12.0}}]", 0);
}

void test_matrix_times_elementwise() {
    assert_eval_eq("Matrix[{{1, 2}, {3, 4}}] * Matrix[{{5, 6}, {7, 8}}]",
                   "Matrix[{{5.0, 12.0}, {21.0, 32.0}}]", 0);
}

void test_matrix_plus_matches_list_path() {
    assert_eval_eq("{{1, 2}, {3, 4}} + {{5, 6}, {7, 8}}", "{{6, 8}, {10, 12}}", 0);
    assert_eval_eq("Normal[Matrix[{{1, 2}, {3, 4}}] + Matrix[{{5, 6}, {7, 8}}]]",
                   "{{6.0, 8.0}, {10.0, 12.0}}", 0);
}

void test_matrix_times_matches_list_path() {
    assert_eval_eq("{{1, 2}, {3, 4}} * {{5, 6}, {7, 8}}", "{{5, 12}, {21, 32}}", 0);
    assert_eval_eq("Normal[Matrix[{{1, 2}, {3, 4}}] * Matrix[{{5, 6}, {7, 8}}]]",
                   "{{5.0, 12.0}, {21.0, 32.0}}", 0);
}

void test_matrix_plus_shape_mismatch_degrades() {
    /* Mismatched shapes: falls through and stays symbolic/unevaluated. */
    assert_eval_eq("Matrix[{1, 2}] + Matrix[{1, 2, 3}]",
                   "Matrix[{1.0, 2.0}] + Matrix[{1.0, 2.0, 3.0}]", 0);
}

int main() {
    symtab_init();
    core_init();

    TEST(test_matrix_construct_rank2);
    TEST(test_matrix_construct_rank1);
    TEST(test_matrix_construct_rank3);
    TEST(test_matrix_construct_reals);
    TEST(test_matrix_reject_jagged);
    TEST(test_matrix_reject_symbolic);
    TEST(test_matrix_reject_non_list);
    TEST(test_matrix_normal_roundtrip);

    TEST(test_matrix_head);
    TEST(test_matrix_dimensions);
    TEST(test_matrix_matrixq);
    TEST(test_matrix_vectorq);
    TEST(test_matrix_listq_always_false);

    TEST(test_matrix_equality);
    TEST(test_matrix_sameq);

    TEST(test_matrix_dot_matrix_matrix);
    TEST(test_matrix_dot_vector_vector);
    TEST(test_matrix_dot_matrix_vector);
    TEST(test_matrix_dot_matches_list_path);
    TEST(test_matrix_dot_shape_mismatch);

    TEST(test_matrix_plus_elementwise);
    TEST(test_matrix_times_elementwise);
    TEST(test_matrix_plus_matches_list_path);
    TEST(test_matrix_times_matches_list_path);
    TEST(test_matrix_plus_shape_mismatch_degrades);

    printf("All Matrix tests passed.\n");
    return 0;
}
