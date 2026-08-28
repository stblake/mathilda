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

/* --- Follow-up: siblings of the "identity that is not an identity in floating
       point" bug the first adversarial pass found --- */

void test_quantile_interpolation_no_overflow(void) {
    /* The interpolating branch must not form A[j+1] - A[j]. Two neighbours of
     * opposite sign near the top of the double range make that difference
     * round to Infinity, and Infinity survives the rest of the expression:
     * this returned inf.0 where the answer is exactly 0. The convex form
     * (1-w) A[j] + w A[j+1] never leaves the interval it interpolates. */
    assert_eval_eq("Quantile[{-1.0*10^308, 1.0*10^308}, 1/2, {{1/2, 0}, {0, 1}}]", "0.0", 0);
    /* The same list through Quartiles, which reaches the engine's fractional
     * weight on its default parameters. The overflow condition is a property of
     * the NEIGHBOUR PAIR the engine lands on -- opposite-sign and huge, ADJACENT
     * after sorting -- not of the list, so the list has to be exactly this: a
     * four-element list with small values in the middle never overflows and
     * would pin nothing. Middle entry was inf.0 before. */
    assert_eval_eq("Quartiles[{-1.0*10^308, 1.0*10^308}]", "{-1e+308, 0.0, 1e+308}", 0);
    /* The NDArray kernel (ndreduce.c) carries its own copy of the same formula
     * and must not disagree with the boxed path on the same data. */
    assert_eval_eq("Quartiles[NDArray[{-1.0*10^308, 1.0*10^308}, DataType->\"float64\"]]",
                   "NDArray[{-1e+308, 0.0, 1e+308}]", 0);
    assert_eval_eq("Quartiles[NDArray[{1., 2., 3., 4.}, DataType->\"float64\"]]",
                   "NDArray[{1.5, 2.5, 3.5}]", 0);
}

void test_quantile_extrapolating_weight_no_nan(void) {
    /* w = c + d g is user-controlled through {{a,b},{c,d}} and nothing bounds it
     * to [0,1]. Outside that interval the convex form is the WRONG one: (1-w)
     * and w have opposite signs, each product can overflow on its own, and
     * inf + -inf is NaN -- on input the difference form gets exactly right,
     * because equal neighbours make the difference zero. Both of these returned
     * nan.0 while the convex form was used unconditionally. */
    assert_eval_eq("Quantile[{1.0*10^308, 1.0*10^308}, 1/2, {{1/2, 0}, {0, 100}}]", "1e+308", 0);
    assert_eval_eq("Quantile[{1.0*10^308, 1.0*10^308}, 1/2, {{1/2, 0}, {-100, 0}}]", "1e+308", 0);
    /* An ordinary extrapolating weight on ordinary data still extrapolates. */
    assert_eval_eq("Quantile[{1, 2, 3, 4}, 1/2, {{1/2, 0}, {0, 4}}]", "4", 0);
}

void test_stats_gate_accepts_complex_with_zero_imaginary(void) {
    /* Complex[x, 0] is a real number wearing a Complex head, and it really does
     * reach the gate: builtin_complex normalises an int64 or double zero away
     * but not an MPFR one, so z + Conjugate[z] at 30 digits is Complex[4.0, 0.0].
     * Rejecting on the head alone declined a genuinely real value. */
    assert_eval_eq("Median[{1, N[2 + I, 30] + N[Conjugate[2 + I], 30], 3}]", "3", 0);
}

void test_quantile_interpolation_exact_unchanged(void) {
    /* ... and the rewrite must not disturb the exact path it shares: a
     * fractional weight over exact rationals still returns an exact rational,
     * with the same value as before. */
    assert_eval_eq("Quantile[{1/3, 2/3, 1, 4/3}, 1/2, {{1/2, 0}, {0, 1}}]", "5/6", 0);
    assert_eval_eq("Quantile[{1, 2, 3, 4}, 1/2, {{1/2, 0}, {0, 1}}]", "5/2", 0);
}

void test_stats_gate_rejects_evaluated_complex(void) {
    /* The element gate promises "a rectangular array of real numbers", and used
     * to test that by searching for the literal symbol I. An evaluated complex
     * number does not contain it -- 2 + I is Complex[2, 1] -- so complex data
     * walked straight through and every head in the subsystem answered anyway:
     * Median returned 3, Quartiles returned a complex third quartile. All of
     * these must decline instead. */
    assert_eval_eq("Quantile[{1, 2 + I, 3}, 1/2]", "Quantile[{1, 2 + I, 3}, 1/2]", 0);
    assert_eval_eq("Median[{1, 2 + I, 3}]", "Median[{1, 2 + I, 3}]", 0);
    assert_eval_eq("Quartiles[{1, 2 + I, 3, 4}]", "Quartiles[{1, 2 + I, 3, 4}]", 0);
    assert_eval_eq("MeanDeviation[{1, 2 + I, 3}]", "MeanDeviation[{1, 2 + I, 3}]", 0);
    assert_eval_eq("InterquartileRange[{1, 2 + I, 3, 4}]", "InterquartileRange[{1, 2 + I, 3, 4}]", 0);
    /* The literal-I form kept working all along; pin it so a future rewrite of
     * the gate cannot trade one hole for the other. */
    assert_eval_eq("Quantile[{1, I, 3}, 1/2]", "Quantile[{1, I, 3}, 1/2]", 0);
}

void test_stats_gate_still_accepts_exact_reals(void) {
    /* The leaf fast path added to that gate must accept exactly what the
     * evaluator round-trip accepted: machine integers, machine reals, exact
     * rationals, bignums, and (the long way round) exact irrationals. */
    assert_eval_eq("Quantile[{1/2, 3/2, 5/2, 7/2}, 1/2]", "3/2", 0);
    assert_eval_eq("Quantile[{10^30, 2*10^30, 3*10^30, 4*10^30}, 1/2, {{1/2, 0}, {0, 1}}]",
                   "2500000000000000000000000000000", 0);
    assert_eval_eq("Quantile[{1., 2., 3., 4.}, 1/4]", "1.0", 0);
    /* An exact irrational is not a numeric leaf, so it still takes the
     * evaluator round-trip; it must still be accepted as real. (Asserted on a
     * constant list so this does not also pin Sort's ordering of symbolic
     * reals, which is a separate question.) */
    assert_eval_eq("Median[{Pi, Pi, Pi}]", "Pi", 0);
}

/* --- Coverage the STATS-1 phase-1 set did not reach --- */

void test_deviation_matrix_columnwise(void) {
    /* Only InterquartileRange had a matrix test; MeanDeviation and
     * MedianDeviation take the same columnwise branch and had none. */
    assert_eval_eq("MeanDeviation[{{1, 2}, {3, 4}, {5, 6}}]", "{4/3, 4/3}", 0);
    assert_eval_eq("MedianDeviation[{{1, 2}, {3, 10}, {5, 6}}]", "{2, 4}", 0);
}

void test_quantile_qlist_out_of_range_declines(void) {
    /* A q-LIST with one bad entry must decline as a whole, not return a
     * partially built list. Only the scalar-q form was pinned. */
    assert_eval_eq("Quantile[{1, 2, 3}, {1/2, 2}]", "Quantile[{1, 2, 3}, {1/2, 2}]", 0);
    assert_eval_eq("Quantile[{1, 2, 3}, {-1/2, 1/2}]", "Quantile[{1, 2, 3}, {-1/2, 1/2}]", 0);
}

void test_iqr_reals(void) {
    /* IQR had exact-integer and matrix coverage only. */
    assert_eval_eq("InterquartileRange[{1., 2., 3., 4., 5., 6., 7., 8.}]", "4.0", 0);
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
    TEST(test_quantile_interpolation_no_overflow);
    TEST(test_quantile_extrapolating_weight_no_nan);
    TEST(test_stats_gate_accepts_complex_with_zero_imaginary);
    TEST(test_quantile_interpolation_exact_unchanged);
    TEST(test_stats_gate_rejects_evaluated_complex);
    TEST(test_stats_gate_still_accepts_exact_reals);
    TEST(test_deviation_matrix_columnwise);
    TEST(test_quantile_qlist_out_of_range_declines);
    TEST(test_iqr_reals);
    printf("All quantile family tests passed (phase 1 set).\n");
    return 0;
}
