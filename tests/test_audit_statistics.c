/* Builtin-audit regressions: statistics (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* Median#1 (high): "Median sorts in canonical order, not numeric order: wrong
 * median when data contains exact numeric non-rationals (Sqrt, -Pi)".  The
 * order statistics now sort by certified numeric value (src/numeric_order.c). */
static void test_median_numeric_order(void) {
    audit_true("Median#1", "Median[{Sqrt[2], 1, 3}] === Sqrt[2]");
    audit_true("Median#1", "Median[{-Pi, 0, -1}] === -1");
    audit_true("Median#1", "Median[{Log[3], 10, Sqrt[5]}] === Sqrt[5]");
    /* Even length: the mean of the two middle values by value. */
    audit_true("Median#1", "PossibleZeroQ[Median[{Pi, 1, E, 4}] - (E + Pi)/2]");
    /* -Infinity ranks below every real. */
    audit_true("Median#1", "Median[{-Infinity, 5, 1}] === 1");
    /* Ordinary data is unchanged. */
    audit_true("Median#1", "Median[{3, 1, 2}] === 2");
    audit_true("Median#1", "Median[{1., 2, 1/2, 4}] === 1.5");
    audit_true("Median#1", "Median[{{1, 2}, {3, 4}}] === {2, 3}");
}

/* Median#2 (high): "Median picks the wrong element when an exact constant is
 * mixed with reals". */
static void test_median_constant_with_reals(void) {
    audit_true("Median#2", "Median[{2.72, E, 2.71}] === E");
    audit_true("Median#2", "Median[{3.2, Pi, 3.1}] === Pi");
    /* packed == unpacked for machine data */
    audit_true("Median#2",
               "Module[{v = RandomReal[1, 1001]}, Median[v] === Median[Developer`FromPackedArray[v]]]");
}

/* MedianDeviation#1 (high): "MedianDeviation gives the wrong value for exact
 * non-rational numeric data (orders by canonical form, not by value)". */
static void test_mediandeviation_numeric_order(void) {
    audit_true("MedianDeviation#1", "PossibleZeroQ[MedianDeviation[{Pi, E, Sqrt[2]}] - (Pi - E)]");
    audit_true("MedianDeviation#1", "PossibleZeroQ[MedianDeviation[{Pi, E, 1}] - (Pi - E)]");
    audit_true("MedianDeviation#1",
               "Abs[N[MedianDeviation[{Pi, E, Sqrt[2]}]] - MedianDeviation[N[{Pi, E, Sqrt[2]}]]] < 10^-12");
}

/* InterquartileRange#1 (high): "Wrong IQR for exact irrational data: values
 * ordered canonically instead of numerically". */
static void test_iqr_numeric_order(void) {
    audit_true("InterquartileRange#1", "InterquartileRange[{Sqrt[2], 1, 3}] === 3/2");
    audit_true("InterquartileRange#1",
               "And @@ PossibleZeroQ[Quartiles[{Sqrt[2], 1, 3}] - {3/4 + Sqrt[2]/4, Sqrt[2], 9/4 + Sqrt[2]/4}]");
    audit_true("InterquartileRange#1",
               "Abs[InterquartileRange[N[{Sqrt[2], 1, 3}]] - 3/2] < 10^-12");
}

/* InterquartileRange#2 (high): "Ragged or mixed-depth list hits
 * $RecursionLimit with ~500 Transpose::nmtx messages". */
static void test_iqr_ragged(void) {
    audit_chk("InterquartileRange#2",
              "Check[InterquartileRange[{{1, 2}, 3}], fired, InterquartileRange::rectn]", "fired");
    audit_chk("InterquartileRange#2",
              "Check[InterquartileRange[{{1, 2}, {3, 4}, {5}}], fired, InterquartileRange::rectn]", "fired");
    /* unevaluated, and no Transpose / recursion messages on the way */
    audit_true("InterquartileRange#2",
               "MatchQ[Check[Quiet[InterquartileRange[{{1, 2}, 3}], InterquartileRange::rectn], fired, "
               "{Transpose::nmtx, $RecursionLimit::reclim}], HoldPattern[InterquartileRange[{{1, 2}, 3}]]]");
    /* The same shape guard serves Median and MedianDeviation. */
    audit_chk("InterquartileRange#2", "Check[Median[{{1, 2}, 3}], fired, Median::rectn]", "fired");
    audit_chk("InterquartileRange#2",
              "Check[MedianDeviation[{{1, 2}, {3}}], fired, MedianDeviation::rectn]", "fired");
    /* A rectangular matrix still works columnwise. */
    audit_true("InterquartileRange#2", "InterquartileRange[{{1, 2}, {3, 4}, {5, 6}, {7, 8}}] === {4, 4}");
}

/* Quantile#1 (high): "Quantile orders exact numeric data canonically, not by
 * value, so it returns wrong elements (even the minimum)". */
static void test_quantile_numeric_order(void) {
    audit_true("Quantile#1", "Quantile[{2, Sqrt[3]}, 0] === Sqrt[3]");
    audit_true("Quantile#1", "Quantile[{10, Log[3], Sqrt[5]}, 0] === Log[3]");
    audit_true("Quantile#1",
               "Quantile[{Pi, 3.5, E, 3, Sqrt[2]}, {0, 1/5, 2/5, 3/5, 4/5, 1}] === {Sqrt[2], Sqrt[2], E, 3, Pi, 3.5}");
    audit_true("Quantile#1", "Quantile[{5, 1, 3}, 1/2] === 3");
    audit_true("Quantile#1",
               "Module[{v = RandomReal[1, 999]}, Quantile[v, {1/4, 1/2}] === Quantile[Developer`FromPackedArray[v], {1/4, 1/2}]]");
}

/* Quantile#2 (high): "Quantile on ragged nested list recurses until
 * $RecursionLimit, printing ~1000 Transpose::nmtx messages". */
static void test_quantile_ragged(void) {
    audit_chk("Quantile#2", "Check[Quantile[{{1, 2}, {3}}, 1/2], fired, Quantile::rectn]", "fired");
    audit_true("Quantile#2",
               "MatchQ[Check[Quiet[Quantile[{{1, 2}, {3}}, 1/2], Quantile::rectn], fired, "
               "{Transpose::nmtx, $RecursionLimit::reclim}], HoldPattern[Quantile[{{1, 2}, {3}}, 1/2]]]");
    audit_true("Quantile#2", "Quantile[{{1, 2}, {3, 4}, {5, 7}}, 1/2] === {3, 4}");
}

/* Quartiles#1 (high): "Quartiles on a ragged or mixed nested list recurses
 * until $RecursionLimit, emitting ~1000 Transpose messages". */
static void test_quartiles_ragged(void) {
    audit_chk("Quartiles#1", "Check[Quartiles[{{1, 2}, {3}}], fired, Quartiles::rectn]", "fired");
    audit_chk("Quartiles#1", "Check[Quartiles[{{1, 2}, 3}], fired, Quartiles::rectn]", "fired");
    audit_true("Quartiles#1",
               "MatchQ[Check[Quiet[Quartiles[{{1, 2}, {3}}], Quartiles::rectn], fired, "
               "{Transpose::nmtx, $RecursionLimit::reclim}], HoldPattern[Quartiles[{{1, 2}, {3}}]]]");
    audit_true("Quartiles#1",
               "Quartiles[{{1, 2}, {3, 4}, {5, 6}}] === {{3/2, 3, 9/2}, {5/2, 4, 11/2}}");
}

/* Covariance#1 (high): "Covariance of two matrices with different numbers of
 * rows returns a zero matrix instead of staying unevaluated". */
static void test_covariance_row_mismatch(void) {
    audit_chk("Covariance#1",
              "Check[Covariance[{{1, 2}, {3, 4}, {5, 6}}, {{1, 2}, {3, 4}}], fired, Covariance::vctmat]", "fired");
    audit_true("Covariance#1",
               "MatchQ[Quiet[Covariance[{{1, 2}, {3, 4}}, {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]], "
               "HoldPattern[Covariance[{{1, 2}, {3, 4}}, {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}]]]");
    audit_chk("Covariance#1",
              "Check[Correlation[{{1, 2}, {3, 4}}, {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}}], fired, Correlation::vctmat]", "fired");
    /* the packed path re-evaluates through the same guard */
    audit_chk("Covariance#1",
              "Check[Covariance[{{1., 2.}, {3., 4.}, {5., 6.}}, {{1., 2.}, {3., 4.}}], fired, Covariance::vctmat]", "fired");
    audit_chk("Covariance#1", "Check[Covariance[{1, 2, 3}, {1, 2}], fired, Covariance::vctmat]", "fired");
    /* Equal row counts are unchanged. */
    audit_true("Covariance#1",
               "Covariance[{{1, 2}, {3, 5}, {5, 6}}, {{1, 2}, {3, 4}, {4, 4}}] === {{3, 2}, {19/6, 7/3}}");
    audit_true("Covariance#1", "Covariance[{1, 2, 4}, {1, 3, 2}] === 1/2");
}

/* Wave 3 S12/S13: declared argument counts and option checks
 * (src/arity.h, src/arity_table.c).  The evaluator checks a call
 * against its head's declared arity before the C code runs and
 * prints argx/argr/argrx/argt(u)/argb(u)/argm(u), nonopt or optx;
 * the call stays unevaluated. */
static void test_wave3_arity_options(void) {
    audit_chk("CentralMoment#3", "Quiet[Check[CentralMoment[{1, 2, 3}], fired, CentralMoment::argr]]", "fired");
    audit_chk("CentralMoment#3", "Quiet[CentralMoment[{1, 2, 3}]]", "CentralMoment[{1, 2, 3}]");
    audit_chk("Covariance#2", "Quiet[Check[Covariance[{1, 2}, {3, 4}, {5, 6}], fired, Covariance::argt]]", "fired");
    audit_chk("Covariance#2", "Quiet[Covariance[{1, 2}, {3, 4}, {5, 6}]]", "Covariance[{1, 2}, {3, 4}, {5, 6}]");
}

/* ================================================================== */
/* Wave 3, shared fix S11: Indeterminate / infinity / pole values     */
/* (src/infinity_rules.c and the per-head rules it does not cover).   */
/* ================================================================== */

static void test_s11_special_values(void) {
    /* CentralMoment#1: the 0th central moment is 1 */
    audit_true("CentralMoment#1", "CentralMoment[{1, 2, 3}, 0] === 1");
    audit_true("CentralMoment#1", "CentralMoment[{1., 2., 3.}, 0] === 1.");
    audit_true("CentralMoment#1", "CentralMoment[{1, 2, 3}, 2] === 2/3");
}

int main(void) {
    audit_init();
    TEST(test_median_numeric_order);
    TEST(test_median_constant_with_reals);
    TEST(test_mediandeviation_numeric_order);
    TEST(test_iqr_numeric_order);
    TEST(test_iqr_ragged);
    TEST(test_quantile_numeric_order);
    TEST(test_quantile_ragged);
    TEST(test_quartiles_ragged);
    TEST(test_covariance_row_mismatch);
    TEST(test_wave3_arity_options);
    TEST(test_s11_special_values);
    return audit_finish("audit_statistics");
}
