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

/* --- Regressions pinned by the STATS-1 adversarial review --- */

void test_quantile_selects_not_recomputes(void) {
    /* The default parameters (c=1, d=0) must SELECT x_(Ceiling[nq]), not evaluate
     * A[j] + 1*(A[j+1]-A[j]): that identity does not hold in floating point.
     * With a huge negative first element the subtraction rounds and the sum
     * returned 0. instead of the element. */
    assert_eval_eq("Quantile[{-1.0*10^308, 2.0, 3.0, 5.0}, 3/10]", "2.0", 0);
    /* ... and it must not drag exact data to Real when q is inexact. */
    assert_eval_eq("Quantile[{1, 2, 3, 4}, 0.3]", "2", 0);
    assert_eval_eq("Quantile[{11, 21, 31}, 0.5]", "21", 0);
}

void test_quantile_exact_irrational_q(void) {
    /* An exact irrational q is NumericQ: it must answer as its N[] form does,
     * not silently decline. n=4, h = 4 Pi/4 = Pi -> Ceiling = 4. */
    assert_eval_eq("Quantile[{1, 2, 3, 4}, Pi/4]", "4", 0);
}

void test_quantile_bad_params_decline(void) {
    /* Symbolic parameters must decline, never yield Indeterminate or a
     * half-evaluated expression; a non-List head is not a parameter matrix. */
    assert_eval_eq("Quantile[{1, 2, 3}, 1/2, {{a, 0}, {1, 0}}]",
                   "Quantile[{1, 2, 3}, 1/2, {{a, 0}, {1, 0}}]", 0);
    assert_eval_eq("Quantile[{1, 2, 3}, 1/2, f[{0, 0}, {1, 0}]]",
                   "Quantile[{1, 2, 3}, 1/2, f[{0, 0}, {1, 0}]]", 0);
}

void test_deviation_nonfinite_declines(void) {
    /* Infinity is NumericQ and free of I, so it passes the element gate; the
     * composed tree cannot reduce, so these must decline rather than hand back
     * a half-evaluated expression. */
    assert_eval_eq("MeanDeviation[{1, 2, Infinity}]", "MeanDeviation[{1, 2, Infinity}]", 0);
    assert_eval_eq("MedianDeviation[{1, 2, ComplexInfinity}]",
                   "MedianDeviation[{1, 2, ComplexInfinity}]", 0);
}

/* --- Quartiles regression + the integer-h correction --- */

void test_quartiles_regression(void) {
    /* AC-12: default Quartiles unchanged after the engine extraction */
    assert_eval_eq("Quartiles[{1, 2, 3, 4}]", "{3/2, 5/2, 7/2}", 0);
    assert_eval_eq("Quartiles[{1, 2, 3, 4, 5, 6, 7, 8}]", "{5/2, 9/2, 13/2}", 0);
}

void test_quartiles_mixed_exactness_at_integer_h(void) {
    /* DELIBERATE, RECORDED deviation from the pre-STATS-1 behavior: at integer h
     * the Floor/Ceiling neighbours coincide, so the upper element is never
     * consulted. Before, A[j] + 0*(A[j+1]-A[j]) let a Real neighbour turn an
     * exact element into a Real ({2.0, 3.5, 5} here). The documented definition
     * gives the element itself. */
    assert_eval_eq("Quartiles[{1, 2, 3., 4, 5, 6}]", "{2, 3.5, 5}", 0);
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
    TEST(test_quantile_selects_not_recomputes);
    TEST(test_quantile_exact_irrational_q);
    TEST(test_quantile_bad_params_decline);
    TEST(test_deviation_nonfinite_declines);
    TEST(test_quartiles_mixed_exactness_at_integer_h);
    TEST(test_quartiles_regression);
    TEST(test_quartiles_param_integer_h);
    printf("All quantile family tests passed (phase 1 set).\n");
    return 0;
}
