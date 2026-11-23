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

int main(void) {
    audit_init();
    TEST(test_fixedpointlist_long_symbolic);
    return audit_finish("audit_functional_programming");
}
