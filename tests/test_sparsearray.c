/* Normal[SparseArray[...]]: the dense nested List a sparse specification
 * denotes (src/sparsearray.c). Normal used to hand the SparseArray back
 * unchanged, so every caller that densified one silently kept the spec. */
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"
#include "core.h"
#include "test_utils.h"
#include <stdio.h>

static void test_explicit_positions(void) {
    assert_eval_eq("Normal[SparseArray[{1 -> 1, 3 -> 2}]]", "{1, 0, 2}", 0);
    assert_eval_eq("Normal[SparseArray[{1 -> 1, 3 -> 2}, 5]]", "{1, 0, 2, 0, 0}", 0);
    assert_eval_eq("Normal[SparseArray[{{1, 1} -> 1, {2, 3} -> 2}]]",
                   "{{1, 0, 0}, {0, 0, 2}}", 0);
    assert_eval_eq("Normal[SparseArray[{{1, 1} -> 1, {2, 3} -> 2}, {3, 3}, x]]",
                   "{{1, x, x}, {x, x, 2}, {x, x, x}}", 0);
    assert_eval_eq("Normal[SparseArray[{{1, 1}, {2, 2}} -> {a, b}, {2, 2}]]",
                   "{{a, 0}, {0, b}}", 0);
    assert_eval_eq("Normal[SparseArray[{{1, 2, 1} -> 7}]]", "{{{0}, {7}}}", 0);
    /* first rule wins; negative indices count from the end */
    assert_eval_eq("Normal[SparseArray[{1 -> a, 1 -> b}]]", "{a}", 0);
    assert_eval_eq("Normal[SparseArray[{-1 -> 9}, 4]]", "{0, 0, 0, 9}", 0);
    assert_eval_eq("Normal[SparseArray[{}, {2, 2}]]", "{{0, 0}, {0, 0}}", 0);
    assert_eval_eq("Normal[SparseArray[{2 -> 1.5}, 3, 0.]]", "{0.0, 1.5, 0.0}", 0);
}

static void test_patterns_and_bands(void) {
    assert_eval_eq("Normal[SparseArray[{{i_, i_} -> 1}, {3, 3}]]",
                   "{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}", 0);
    assert_eval_eq("Normal[SparseArray[{i_, j_} :> i + 10 j, {2, 3}]]",
                   "{{11, 21, 31}, {12, 22, 32}}", 0);
    assert_eval_eq("Normal[SparseArray[{{i_, j_} /; i > j -> 1}, {3, 3}]]",
                   "{{0, 0, 0}, {1, 0, 0}, {1, 1, 0}}", 0);
    assert_eval_eq("Normal[SparseArray[{i_} :> i^2, 4]]", "{1, 4, 9, 16}", 0);
    assert_eval_eq("Normal[SparseArray[Band[{1, 1}] -> 1, {3, 3}]]",
                   "{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}", 0);
    assert_eval_eq("Normal[SparseArray[{Band[{1, 2}] -> {a, b}, Band[{2, 1}] -> -1}, {3, 3}]]",
                   "{{0, a, 0}, {-1, 0, b}, {0, -1, 0}}", 0);
    /* a pattern rule without dims determines no size: left alone */
    assert_eval_eq("Normal[SparseArray[{{i_, i_} -> 1}]]",
                   "Normal[SparseArray[List[Rule[List[Pattern[i, Blank[]], Pattern[i, Blank[]]], 1]]]]", 1);
}

static void test_dense_and_internal_forms(void) {
    assert_eval_eq("Normal[SparseArray[{{1, 2}, {3, 4}}]]", "{{1, 2}, {3, 4}}", 0);
    assert_eval_eq("Normal[SparseArray[{1, 2}, 4]]", "{1, 2, 0, 0}", 0);
    assert_eval_eq("Normal[SparseArray[{1, 2, 3}, 2]]", "{1, 2}", 0);
    /* InputForm of Mathematica's own SparseArray */
    assert_eval_eq("Normal[SparseArray[Automatic, {3, 3}, 0, "
                   "{1, {{0, 1, 2, 3}, {{1}, {2}, {3}}}, {1, 1, 1}}]]",
                   "{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}", 0);
    assert_eval_eq("Normal[SparseArray[Automatic, {5}, 0, {1, {{0, 2}, {{1}, {3}}}, {a, b}}]]",
                   "{a, 0, b, 0, 0}", 0);
    /* nested inside another expression */
    assert_eval_eq("Normal[{SparseArray[{2 -> 1}], 3}]", "{{0, 1}, 3}", 0);
    /* the SparseArray itself stays inert */
    assert_eval_eq("SparseArray[{1 -> 1}]", "SparseArray[List[Rule[1, 1]]]", 1);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_explicit_positions);
    TEST(test_patterns_and_bands);
    TEST(test_dense_and_internal_forms);
    printf("All SparseArray tests passed!\n");
    return 0;
}
