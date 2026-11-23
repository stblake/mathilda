/* Builtin-audit regressions: numerical calculus (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* NSeries#1 (critical): an order of 2^31 truncated to a negative int and the
 * sampler overran its buffer.  Orders past the documented maximum are
 * refused with NSeries::ivar and the call stays unevaluated. */
static void test_nseries_huge_order(void) {
    audit_chk("NSeries#1", "Check[NSeries[Exp[x], {x, 0, 2^31}], fired, NSeries::ivar]", "fired");
    audit_chk("NSeries#1", "Quiet[NSeries[Exp[x], {x, 0, 2^31}]]", "NSeries[E^x, {x, 0, 2147483648}]");
    audit_chk("NSeries#1", "Quiet[NSeries[Exp[x], {x, 0, 2^64}]]",
              "NSeries[E^x, {x, 0, 18446744073709551616}]");
    audit_true("NSeries#1",
               "Abs[SeriesCoefficient[NSeries[Exp[x], {x, 0, 4}], 3] - 1/6] < 10^-10");
}

int main(void) {
    audit_init();
    TEST(test_nseries_huge_order);
    return audit_finish("audit_numerical_calculus");
}
