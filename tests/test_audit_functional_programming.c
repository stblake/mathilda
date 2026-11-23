/* Builtin-audit regressions: functional programming (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* FixedPointList#1 (high): "FixedPointList segfaults on long symbolic
 * iterations (e.g. FixedPointList[f, x, 100000]); run time grows
 * quadratically".  SameQ now sits behind the memoised-hash pre-check, so each
 * step is O(1) for a growing iterate. */
static void test_fixedpointlist_long_symbolic(void) {
    audit_chk("FixedPointList#1", "Length[FixedPointList[f, x, 100000]]", "100001");
    audit_true("FixedPointList#1",
               "Last[FixedPointList[f, x, 100000]] === Nest[f, x, 100000]");
    audit_true("FixedPointList#1",
               "FixedPoint[g, x, 100000] === Nest[g, x, 100000]");
    /* Convergence is still detected by SameQ. */
    audit_chk("FixedPointList#1", "FixedPointList[Floor[#/2] &, 20]", "{20, 10, 5, 2, 1, 0, 0}");
    audit_chk("FixedPointList#1", "FixedPoint[Floor[#/2] &, 2^100]", "0");
    audit_chk("FixedPointList#1", "FixedPointList[# /. a -> b &, {a, a}]", "{{a, a}, {b, b}, {b, b}}");
    /* 0. and -0. are SameQ, so they hash alike and still converge. */
    audit_chk("FixedPointList#1", "Length[FixedPointList[-# &, 0., 10]]", "2");
    audit_true("FixedPointList#1", "Hash[0.] === Hash[-0.]");
    audit_chk("FixedPointList#1", "FixedPointList[f, x, SameTest -> (True &)]", "{x, f[x]}");
}

/* TakeLargest#1 (high): "TakeLargest ranks symbolic numeric constants (Pi, E,
 * -Pi, Sqrt[2]) above every explicit number".  Keys rank in certified numeric
 * order (src/numeric_order.c), not canonical order. */
static void test_takelargest_numeric_constants(void) {
    audit_true("TakeLargest#1", "TakeLargest[{1, Pi, 4}, 1] === {4}");
    audit_true("TakeLargest#1", "TakeLargest[{-Pi, 1}, 1] === {1}");
    audit_true("TakeLargest#1", "TakeLargest[{E, 3, Pi}, 3] === {Pi, 3, E}");
    audit_true("TakeLargest#1", "TakeSmallest[{Pi, 3, Sqrt[2], 1/2}, 2] === {1/2, Sqrt[2]}");
    /* close values are certified, not read through doubles */
    audit_true("TakeLargest#1", "TakeLargest[{314159265358979/10^14, Pi}, 1] === {Pi}");
    audit_true("TakeLargest#1", "RankedMin[{Pi, 3, E}, 2] === 3");
    audit_true("TakeLargest#1", "RankedMax[{Pi, 314159265358979/10^14, E}, 1] === Pi");
    /* ties settle by position; machine data unchanged */
    audit_true("TakeLargest#1", "TakeLargest[{3, 3., 1}, 2] === {3, 3.}");
    audit_true("TakeLargest#1", "TakeLargest[{5, 1, 4, 2}, 2] === {5, 4}");
    audit_true("TakeLargest#1",
               "Module[{v = RandomReal[1, 500]}, TakeLargest[v, 5] === TakeLargest[Developer`FromPackedArray[v], 5]]");
}

/* TakeLargest#2 (high): "TakeLargest treats -Infinity as larger than finite
 * numbers". */
static void test_takelargest_minus_infinity(void) {
    audit_true("TakeLargest#2", "TakeLargest[{-Infinity, 3}, 1] === {3}");
    audit_true("TakeLargest#2", "TakeLargest[{-Infinity, Infinity, 3}, 2] === {Infinity, 3}");
    audit_true("TakeLargest#2", "TakeSmallest[{-Infinity, Infinity, 3}, 1] === {-Infinity}");
}

/* TakeLargestBy#1 (high): "Exact numeric keys (Pi, E, Sqrt[2]) are ranked by
 * canonical order, not numeric value". */
static void test_takelargestby_numeric_keys(void) {
    audit_true("TakeLargestBy#1", "TakeLargestBy[{1, 2, 3}, {Pi, 3, E}[[#]] &, 2] === {1, 2}");
    audit_true("TakeLargestBy#1", "TakeLargestBy[{1, 2, 3}, {Sqrt[2], 1.5, 1}[[#]] &, 1] === {2}");
    audit_true("TakeLargestBy#1", "TakeLargestBy[{Pi, 3, E}, N, 2] === {Pi, 3}");
}

/* TakeSmallestBy#1 (high): "TakeSmallestBy ranks exact numeric keys like
 * Sqrt[2] and E by canonical order, not numeric value". */
static void test_takesmallestby_numeric_keys(void) {
    audit_true("TakeSmallestBy#1", "TakeSmallestBy[{Sqrt[2], 1.5, 1}, Identity, 2] === {1, Sqrt[2]}");
    audit_true("TakeSmallestBy#1", "TakeSmallestBy[{E, 2.8, 3}, Identity, 1] === {E}");
    audit_true("TakeSmallestBy#1", "TakeSmallestBy[{Pi, 3, 4}, Identity, 1] === {3}");
}

int main(void) {
    audit_init();
    TEST(test_fixedpointlist_long_symbolic);
    TEST(test_takelargest_numeric_constants);
    TEST(test_takelargest_minus_infinity);
    TEST(test_takelargestby_numeric_keys);
    TEST(test_takesmallestby_numeric_keys);
    return audit_finish("audit_functional_programming");
}
