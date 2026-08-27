/* test_quantile_family.c -- Quantile, InterquartileRange, MeanDeviation,
 * MedianDeviation (ticket STATS-1). One assert_eval_eq per acceptance-criteria
 * row in thoughts/shared/tickets/STATS-1/plan.md; AC ids in the comments. */

#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"

/* --- Quantile: default parameters (Wolfram {{0,0},{1,0}}) --- */

void test_quantile_default_half(void) {
    /* AC-1: left-continuous type-1 result, NOT the median */
    assert_eval_eq("Quantile[{1, 2, 3, 4}, 1/2]", "2", 0);
}

void test_quantile_default_quarter(void) {
    /* AC-2 */
    assert_eval_eq("Quantile[{1, 2, 3, 4}, 1/4]", "1", 0);
}

void test_quantile_q_list(void) {
    /* AC-3 */
    assert_eval_eq("Quantile[{1, 2, 3, 4}, {1/4, 3/4}]", "{1, 3}", 0);
}

void test_quantile_edges(void) {
    /* AC-4a / AC-4b: q = 1 and q = 0 clamp to max/min */
    assert_eval_eq("Quantile[{1, 2, 3, 4}, 1]", "4", 0);
    assert_eval_eq("Quantile[{1, 2, 3, 4}, 0]", "1", 0);
}

void test_quantile_explicit_params(void) {
    /* AC-5: Quartiles' parameters reproduce the interpolated median */
    assert_eval_eq("Quantile[{1, 2, 3, 4}, 1/2, {{1/2, 0}, {0, 1}}]", "5/2", 0);
}

void test_quantile_reals(void) {
    /* AC-6: real input, real result */
    assert_eval_eq("Quantile[{1., 2., 3., 4.}, 1/2]", "2.0", 0);
}

void test_quantile_unsorted(void) {
    /* AC-7: sorts first */
    assert_eval_eq("Quantile[{3, 1, 4, 2}, 1/2]", "2", 0);
}

void test_quantile_symbolic_declines(void) {
    /* AC-11a / AC-11b: symbolic data or q stays unevaluated */
    assert_eval_eq("Quantile[x, 1/2]", "Quantile[x, 1/2]", 0);
    assert_eval_eq("Quantile[{1, 2}, q]", "Quantile[{1, 2}, q]", 0);
}

void test_quantile_matrix_columnwise(void) {
    /* AC-13 */
    assert_eval_eq("Quantile[{{1, 2}, {3, 4}}, 1/2]", "{1, 2}", 0);
}

void test_quantile_empty_declines(void) {
    /* AC-15 */
    assert_eval_eq("Quantile[{}, 1/2]", "Quantile[{}, 1/2]", 0);
}

void test_quantile_out_of_range_q(void) {
    /* AC-16: message (stdout) + unevaluated */
    assert_eval_eq("Quantile[{1, 2, 3}, 2]", "Quantile[{1, 2, 3}, 2]", 0);
}

void test_quantile_ndarray_exact(void) {
    /* AC-18: visible NDArray materialised to the exact path */
    assert_eval_eq("Quantile[NDArray[{1, 2, 3, 4}, DataType->\"int64\"], 1/2]", "2", 0);
}

void test_quantile_singleton(void) {
    /* edge: single element at every q */
    assert_eval_eq("Quantile[{7}, 0]", "7", 0);
    assert_eval_eq("Quantile[{7}, 1/2]", "7", 0);
    assert_eval_eq("Quantile[{7}, 1]", "7", 0);
}

/* --- InterquartileRange, MeanDeviation, MedianDeviation --- */

void test_iqr_basic(void) {
    /* AC-8 */
    assert_eval_eq("InterquartileRange[{1, 2, 3, 4, 5, 6, 7, 8}]", "4", 0);
}

void test_iqr_matrix_columnwise(void) {
    /* AC-17: per-column IQR, including the 3-column collision case the plan
     * review flagged (a 3-list of column triples must not be mistaken for a
     * vector's own {q1,q2,q3}). */
    assert_eval_eq("InterquartileRange[{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}, {10, 11, 12}}]", "{6, 6, 6}", 0);
}

void test_mean_deviation(void) {
    /* AC-9, AC-14 */
    assert_eval_eq("MeanDeviation[{1, 2, 3, 4}]", "1", 0);
    assert_eval_eq("MeanDeviation[{1/2, 3/2}]", "1/2", 0);
}

void test_median_deviation(void) {
    /* AC-10 */
    assert_eval_eq("MedianDeviation[{1, 2, 3, 4}]", "1", 0);
    assert_eval_eq("MedianDeviation[{1, 2, 3, 10}]", "1", 0);
}

void test_meandeviation_ndarray(void) {
    /* AC-19: visible NDArray materialised, real result */
    assert_eval_eq("MeanDeviation[NDArray[{1., 2., 3., 4.}, DataType->\"float64\"]]", "1.0", 0);
}

void test_deviation_declines(void) {
    /* symbolic + empty declines for the phase-2 heads */
    assert_eval_eq("MeanDeviation[{a, b}]", "MeanDeviation[{a, b}]", 0);
    assert_eval_eq("InterquartileRange[{}]", "InterquartileRange[{}]", 0);
    assert_eval_eq("MedianDeviation[x]", "MedianDeviation[x]", 0);
}

/* --- Quartiles regression + the integer-h correction --- */

void test_quartiles_regression(void) {
    /* AC-12: default Quartiles unchanged after the engine extraction */
    assert_eval_eq("Quartiles[{1, 2, 3, 4}]", "{3/2, 5/2, 7/2}", 0);
    assert_eval_eq("Quartiles[{1, 2, 3, 4, 5, 6, 7, 8}]", "{5/2, 9/2, 13/2}", 0);
}

void test_quartiles_param_integer_h(void) {
    /* Parameterized Quartiles at integer h now takes the Wolfram Floor/Ceiling
     * collapse (was s[[h+1]] before STATS-1's shared engine): with {{0,0},{1,0}}
     * h = n q = 1, 2, 3 -> elements 1, 2, 3. Behavior change is deliberate and
     * recorded in the plan's deviation notes. */
    assert_eval_eq("Quartiles[{1, 2, 3, 4}, {{0, 0}, {1, 0}}]", "{1, 2, 3}", 0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_quantile_default_half);
    TEST(test_quantile_default_quarter);
    TEST(test_quantile_q_list);
    TEST(test_quantile_edges);
    TEST(test_quantile_explicit_params);
    TEST(test_quantile_reals);
    TEST(test_quantile_unsorted);
    TEST(test_quantile_symbolic_declines);
    TEST(test_quantile_matrix_columnwise);
    TEST(test_quantile_empty_declines);
    TEST(test_quantile_out_of_range_q);
    TEST(test_quantile_ndarray_exact);
    TEST(test_quantile_singleton);
    TEST(test_iqr_basic);
    TEST(test_iqr_matrix_columnwise);
    TEST(test_mean_deviation);
    TEST(test_median_deviation);
    TEST(test_meandeviation_ndarray);
    TEST(test_deviation_declines);
    TEST(test_quartiles_regression);
    TEST(test_quartiles_param_integer_h);
    printf("All quantile family tests passed (phase 1 set).\n");
    return 0;
}
