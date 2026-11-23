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

/* Range#1 (high): "Range[0., 100., 0.01] drops the final element (10000
 * instead of 10001 elements)".  The count is Floor[(imax - imin)/di] + 1 with
 * a relative slack, and element k is imin + k di (never past imax). */
static void test_range_real_count(void) {
    audit_chk("Range#1", "Length[Range[0., 100., 0.01]]", "10001");
    audit_true("Range#1", "Last[Range[0., 100., 0.01]] === 100.");
    audit_true("Range#1", "Range[0., 0.3, 0.1] === {0., 0.1, 0.2, 0.3}");
    audit_true("Range#1", "Last[Range[0., 1., 0.1]] === 1.");
    audit_true("Range#1", "Range[0.1, 1., 0.3] === {0.1, 0.4, 0.7, 0.1 + 3*0.3}");
    audit_true("Range#1", "Length[Range[0., 0.2999999999999, 0.1]] === 3");
    audit_true("Range#1", "Range[1., 0., -0.5] === {1., 0.5, 0.}");
    audit_true("Range#1", "Range[1, 2, 0.25] === {1., 1.25, 1.5, 1.75, 2.}");
    audit_true("Range#1", "Range[0, 1, 1/3] === {0, 1/3, 2/3, 1}");
    audit_true("Range#1", "Range[5., 1., 1.] === {}");
    audit_true("Range#1", "Length[Range[-1.*^6, 0., 0.1]] === 10000001");
}

/* ArrayPad#1 (high): "ArrayPad silently erases the array's data when the
 * per-level spec has more levels than the array". */
static void test_arraypad_spec_too_deep(void) {
    audit_chk("ArrayPad#1", "Check[ArrayPad[{1, 2}, {{1, 1}, {1, 1}}], fired, ArrayPad::depth]", "fired");
    audit_true("ArrayPad#1",
               "MatchQ[Quiet[ArrayPad[{1, 2}, {{1, 1}, {1, 1}}]], HoldPattern[ArrayPad[{1, 2}, {{1, 1}, {1, 1}}]]]");
    audit_chk("ArrayPad#1",
              "Check[ArrayPad[{{1, 2}, {3, 4}}, {{1, 1}, {1, 1}, {1, 1}}], fired, ArrayPad::depth]", "fired");
    /* a ragged level is not an array level */
    audit_chk("ArrayPad#1", "Check[ArrayPad[{1, {2, 3}}, {{1, 1}, {1, 1}}], fired, ArrayPad::depth]", "fired");
    audit_true("ArrayPad#1", "ArrayPad[{1, 2}, 1] === {0, 1, 2, 0}");
    audit_true("ArrayPad#1", "ArrayPad[{1, 2}, {{1, 2}}] === {0, 1, 2, 0, 0}");
    audit_true("ArrayPad#1",
               "ArrayPad[{{1, 2}, {3, 4}}, {{1, 0}, {0, 1}}] === {{0, 0, 0}, {1, 2, 0}, {3, 4, 0}}");
    audit_true("ArrayPad#1",
               "ArrayPad[{{1, 2}, {3, 4}}, 1] === {{0, 0, 0, 0}, {0, 1, 2, 0}, {0, 3, 4, 0}, {0, 0, 0, 0}}");
}

/* Permutations#1 (high): "Permutations[list, n] returns permutations of
 * exactly n elements instead of at most n". */
static void test_permutations_at_most_n(void) {
    audit_chk("Permutations#1", "Permutations[{a, b, c}, 2]",
              "{{}, {a}, {b}, {c}, {a, b}, {a, c}, {b, a}, {b, c}, {c, a}, {c, b}}");
    audit_chk("Permutations#1", "Permutations[{a, b}, Infinity]", "{{}, {a}, {b}, {a, b}, {b, a}}");
    audit_chk("Permutations#1", "Permutations[{a, a, b}, 2]",
              "{{}, {a}, {b}, {a, a}, {a, b}, {b, a}}");
    /* {n} is the exactly-n form */
    audit_chk("Permutations#1", "Permutations[{a, b, c}, {2}]",
              "{{a, b}, {a, c}, {b, a}, {b, c}, {c, a}, {c, b}}");
    audit_chk("Permutations#1", "Permutations[{a, b}, 0]", "{{}}");
    audit_chk("Permutations#1", "Length[Permutations[Range[4], 9]]", "65");
    audit_chk("Permutations#1", "Permutations[{a, b}]", "{{a, b}, {b, a}}");
}

int main(void) {
    audit_init();
    TEST(test_unit_vector_option_only);
    TEST(test_tuples_long_tuple);
    TEST(test_array_huge_dimension);
    TEST(test_range_real_count);
    TEST(test_arraypad_spec_too_deep);
    TEST(test_permutations_at_most_n);
    return audit_finish("audit_lists_and_iteration");
}
