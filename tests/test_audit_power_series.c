/* Builtin-audit regressions: power series (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* Series#1 (critical): Series[x, {x, 0, 2^62}] sized its dense coefficient
 * buffer unchecked and segfaulted.  An order no machine can hold aborts with
 * General::nomem; ordinary orders are unchanged. */
static void test_series_huge_order(void) {
    audit_chk("Series#1", "Quiet[CheckAbort[Series[x, {x, 0, 2^62}], aborted]]", "aborted");
    audit_chk("Series#1",
              "Quiet[Check[CheckAbort[Series[Sin[x], {x, 0, 2^62}], 0], fired, General::nomem]]",
              "fired");
    audit_chk("Series#1", "Quiet[CheckAbort[Series[x, {x, 0, 9223372036854775807}], aborted]]",
              "aborted");
    audit_true("Series#1", "Normal[Series[Sin[x], {x, 0, 5}]] === x - x^3/6 + x^5/120");
    /* A constant is returned verbatim whatever the order (no buffer). */
    audit_chk("Series#1", "Series[y, {x, 0, 2^62}]", "y");
}

int main(void) {
    audit_init();
    TEST(test_series_huge_order);
    return audit_finish("audit_power_series");
}
