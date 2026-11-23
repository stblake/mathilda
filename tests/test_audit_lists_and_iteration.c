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

/* Tuples#1 (high): "Tuples[{a}, 10^6] crashes with SIGSEGV (unbounded
 * recursion in tuples_rec)".  The index tuples come from an odometer loop. */
static void test_tuples_long_tuple(void) {
    audit_chk("Tuples#1", "Length[Tuples[{a}, 10^6]]", "1");
    audit_chk("Tuples#1", "Length[First[Tuples[{a}, 10^6]]]", "1000000");
    audit_true("Tuples#1", "First[Tuples[{a}, 10^6]] === ConstantArray[a, 10^6]");
    /* Lexicographic order and the other forms are unchanged. */
    audit_chk("Tuples#1", "Tuples[{0, 1}, 3]",
              "{{0, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 1}, {1, 0, 0}, {1, 0, 1}, {1, 1, 0}, {1, 1, 1}}");
    audit_chk("Tuples#1", "Tuples[{{a, b}, {1, 2, 3}}]",
              "{{a, 1}, {a, 2}, {a, 3}, {b, 1}, {b, 2}, {b, 3}}");
    audit_chk("Tuples#1", "Tuples[{{a, b}, {}}]", "{}");
    audit_chk("Tuples#1", "Tuples[{}]", "{{}}");
    audit_chk("Tuples#1", "Tuples[{a, b}, {2, 1}]",
              "{{{a}, {a}}, {{a}, {b}}, {{b}, {a}}, {{b}, {b}}}");
}

int main(void) {
    audit_init();
    TEST(test_unit_vector_option_only);
    TEST(test_tuples_long_tuple);
    return audit_finish("audit_lists_and_iteration");
}
