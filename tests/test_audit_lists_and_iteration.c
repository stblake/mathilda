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

/* Array#1 (high): Array[f, 10^15] segfaulted (unchecked malloc of the result
 * slots, which also wrapped at 2^61).  Every level's slot count is checked. */
static void test_array_huge_dimension(void) {
    audit_chk("Array#1", "Quiet[CheckAbort[Array[f, 10^15], aborted]]", "aborted");
    audit_chk("Array#1", "Quiet[Check[CheckAbort[Array[f, 10^15], 0], fired, General::nomem]]",
              "fired");
    audit_chk("Array#1", "Quiet[CheckAbort[Array[f, 2^61], aborted]]", "aborted");
    audit_chk("Array#1", "Quiet[CheckAbort[Array[f, {10^15, 0}], aborted]]", "aborted");
    audit_chk("Array#1", "Quiet[CheckAbort[Array[f, {2^32, 2^32}], aborted]]", "aborted");
    audit_chk("Array#1", "Array[f, 3]", "{f[1], f[2], f[3]}");
    audit_chk("Array#1", "Array[f, {2, 2}, 0]", "{{f[0, 0], f[0, 1]}, {f[1, 0], f[1, 1]}}");
    audit_chk("Array#1", "Array[f, {0, 10^15}]", "{}");
    /* The same class in ArrayReshape (List and packed paths). */
    audit_chk("Array#1", "Quiet[CheckAbort[ArrayReshape[{1}, {10^15, 0}], aborted]]", "aborted");
    audit_chk("Array#1", "Quiet[CheckAbort[ArrayReshape[Range[100], 2^61], aborted]]", "aborted");
    audit_chk("Array#1", "ArrayReshape[Range[6], {2, 3}]", "{{1, 2, 3}, {4, 5, 6}}");
    audit_chk("Array#1", "ArrayReshape[{1, 2}, {2, 2}]", "{{1, 2}, {0, 0}}");
}

int main(void) {
    audit_init();
    TEST(test_unit_vector_option_only);
    TEST(test_tuples_long_tuple);
    TEST(test_array_huge_dimension);
    return audit_finish("audit_lists_and_iteration");
}
