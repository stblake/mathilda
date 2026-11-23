/* Builtin-audit regressions: lists and iteration (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* UnitVector#1 (critical): UnitVector[opt] read args[1] past the end. */
static void test_unit_vector_option_only(void) {
    audit_chk("UnitVector#1",
              "Check[UnitVector[WorkingPrecision -> MachinePrecision], fired, UnitVector::argt]",
              "fired");
    audit_chk("UnitVector#1", "Quiet[UnitVector[WorkingPrecision -> MachinePrecision]]",
              "UnitVector[WorkingPrecision -> MachinePrecision]");
    audit_chk("UnitVector#1", "UnitVector[3, 2]", "{0, 1, 0}");
    audit_true("UnitVector#1", "UnitVector[2, WorkingPrecision -> MachinePrecision] === {0., 1.}");
    /* A length past memory aborts with General::nomem, never mallocs it. */
    audit_chk("UnitVector#1", "Quiet[CheckAbort[UnitVector[10^15, 1], aborted]]", "aborted");
    audit_chk("UnitVector#1",
              "Quiet[Check[CheckAbort[UnitVector[10^15, 1], 0], fired, General::nomem]]",
              "fired");
}

int main(void) {
    audit_init();
    TEST(test_unit_vector_option_only);
    return audit_finish("audit_lists_and_iteration");
}
