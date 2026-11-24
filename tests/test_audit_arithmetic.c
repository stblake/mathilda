/* Builtin-audit regressions: arithmetic (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* ComplexExpand#1 (critical): the TargetFunctions -> Conjugate path freed its
 * real and imaginary parts twice after moving them into the result; the
 * recycled nodes then turned up inside the caller's expression (and the
 * kernel crashed).  The results are also checked for the mathematics. */
static void test_complex_expand_conjugate(void) {
    audit_chk("ComplexExpand#1", "ComplexExpand[Exp[x], TargetFunctions -> Conjugate]", "E^x");
    audit_chk("ComplexExpand#1", "Head[ComplexExpand[Exp[x], TargetFunctions -> Conjugate]]", "Power");
    audit_chk("ComplexExpand#1", "f[ComplexExpand[2^x, TargetFunctions -> Conjugate]]", "f[2^x]");
    audit_chk("ComplexExpand#1", "ComplexExpand[x + 1, TargetFunctions -> Conjugate]", "1 + x");
    audit_chk("ComplexExpand#1", "ComplexExpand[Sin[x], TargetFunctions -> Conjugate]", "Sin[x]");
    /* z complex: re + I im must reassemble z itself. */
    audit_true("ComplexExpand#1",
               "PossibleZeroQ[ComplexExpand[z^2, {z}, TargetFunctions -> Conjugate] - z^2]");
    audit_true("ComplexExpand#1",
               "PossibleZeroQ[ComplexExpand[Sqrt[z], {z}, TargetFunctions -> Conjugate] - Sqrt[z]]");
    /* The Re / Im / Abs wrappers keep their (previously correct) answers. */
    audit_true("ComplexExpand#1",
               "PossibleZeroQ[ComplexExpand[Re[z^2], {z}, TargetFunctions -> Conjugate] - (z^2 + Conjugate[z]^2)/2]");
}

/* Differences#1 (high): Differences[list, 1, -2^63] crashed -- -s overflowed
 * at INT64_MIN, so the |s| >= length guard failed.  The result length is
 * l - n|s|, empty here. */
static void test_differences_int64_min_step(void) {
    audit_chk("Differences#1", "Differences[{1, 2, 3}, 1, -2^63]", "{}");
    audit_chk("Differences#1", "Differences[{1, 2, 3}, 1, -2^62]", "{}");
    audit_chk("Differences#1", "Differences[LinkedList[{1, 2, 3}], 1, -2^63] // Normal", "{}");
    audit_chk("Differences#1", "Differences[{1, 2, 4, 8}, 1, -2]", "{-3, -6}");
    audit_chk("Differences#1", "Differences[{1, 2, 3}, 1, 2]", "{2}");
}

/* Plus#1 (high): "Infinity + I Infinity returns Infinity instead of
 * Indeterminate".  Infinities in different directions have no sum; the
 * complex-directed term was absorbed as if it were finite. */
static void test_plus_infinity_directions(void) {
    audit_chk("Plus#1", "Quiet[Infinity + I Infinity]", "Indeterminate");
    audit_chk("Plus#1", "Quiet[Infinity - I Infinity]", "Indeterminate");
    audit_chk("Plus#1", "Quiet[-Infinity + (1 + I) Infinity]", "Indeterminate");
    audit_chk("Plus#1", "Quiet[ComplexInfinity + I Infinity]", "Indeterminate");
    audit_chk("Plus#1", "Check[Infinity + I Infinity, fired, Infinity::indet]", "fired");
    /* Unchanged: a finite term beside a complex-directed infinity, and the
     * real-direction rules. */
    audit_chk("Plus#1", "3 + I Infinity", "3 + I Infinity");
    audit_chk("Plus#1", "Infinity + 5", "Infinity");
    audit_chk("Plus#1", "Quiet[Infinity - Infinity]", "Indeterminate");
}

/* Power#1 (high): "Finite numeric base raised to Infinity/-Infinity stays
 * unevaluated".  b^Infinity is decided by the certified sign of |b| - 1:
 * 0 inside the unit circle, Infinity (positive real b) or ComplexInfinity
 * outside it, Indeterminate on it; b^-Infinity = (1/b)^Infinity. */
static void test_power_finite_base_infinity(void) {
    audit_chk("Power#1", "2^Infinity", "Infinity");
    audit_chk("Power#1", "(1/2)^Infinity", "0");
    audit_chk("Power#1", "2^(-Infinity)", "0");
    audit_chk("Power#1", "(1/2)^(-Infinity)", "Infinity");
    audit_chk("Power#1", "E^Infinity", "Infinity");
    audit_chk("Power#1", "E^(-Infinity)", "0");
    audit_chk("Power#1", "3^Infinity", "Infinity");
    audit_chk("Power#1", "2.^Infinity", "Infinity");
    audit_chk("Power#1", "(-2)^Infinity", "ComplexInfinity");
    audit_chk("Power#1", "(-1/2)^Infinity", "0");
    audit_chk("Power#1", "(-1/2)^(-Infinity)", "ComplexInfinity");
    audit_chk("Power#1", "(1 + I)^Infinity", "ComplexInfinity");
    audit_chk("Power#1", "(I/2)^Infinity", "0");
    audit_chk("Power#1", "Quiet[1^Infinity]", "Indeterminate");
    audit_chk("Power#1", "Quiet[I^Infinity]", "Indeterminate");
    audit_chk("Power#1", "Quiet[(-1)^Infinity]", "Indeterminate");
    audit_chk("Power#1", "Check[(-1)^Infinity, fired, Infinity::indet]", "fired");
    audit_chk("Power#1", "(Sqrt[2] - 1)^Infinity", "0");
    /* An exact base on the unit circle is not taken from rounding noise. */
    audit_chk("Power#1", "Quiet[(Sqrt[5 + 2 Sqrt[6]] - Sqrt[2] - Sqrt[3] + 1)^Infinity]",
              "Indeterminate");
    /* Symbolic bases and the existing table are unchanged. */
    audit_chk("Power#1", "x^Infinity", "x^Infinity");
    audit_chk("Power#1", "0^Infinity", "0");
    audit_chk("Power#1", "Infinity^2", "Infinity");
}

/* IntervalIntersection#1 (high): "Intervals with symbolic endpoints are
 * reported as disjoint (Interval[]), even X intersected with itself".  An
 * overlap that cannot be decided leaves the call unevaluated. */
static void test_interval_intersection_symbolic(void) {
    audit_chk("IntervalIntersection#1",
              "IntervalIntersection[Interval[{a, b}], Interval[{a, b}]]",
              "IntervalIntersection[Interval[{a, b}], Interval[{a, b}]]");
    audit_chk("IntervalIntersection#1",
              "IntervalIntersection[Interval[{0, a}], Interval[{1, 2}]]",
              "IntervalIntersection[Interval[{0, a}], Interval[{1, 2}]]");
    audit_chk("IntervalIntersection#1",
              "IntervalIntersection[Interval[{0, a}], Interval[{1, 2}]] /. a -> 3/2",
              "Interval[{1, 3/2}]");
    /* Numeric endpoints: overlaps and proven disjointness as before. */
    audit_chk("IntervalIntersection#1",
              "IntervalIntersection[Interval[{0, 3}], Interval[{1, 5}]]", "Interval[{1, 3}]");
    audit_chk("IntervalIntersection#1",
              "IntervalIntersection[Interval[{0, 1}], Interval[{2, 3}]]", "Interval[]");
    audit_chk("IntervalIntersection#1",
              "IntervalIntersection[Interval[{0, 1}, {2, 5}], Interval[{Pi, 7}]]",
              "Interval[{Pi, 5}]");
}

/* Ratios#1 (high): "Level spec deeper than the list recurses into
 * Rational/Times internals".  Below the top level only a List is a level of
 * the array; a deeper spec draws Ratios::lvldep and stays unevaluated.  The
 * same walk in Differences had the same defect. */
static void test_ratios_level_depth(void) {
    audit_chk("Ratios#1", "Quiet[Ratios[{1, 2, 3}, {1, 1}]]", "Ratios[{1, 2, 3}, {1, 1}]");
    audit_chk("Ratios#1", "Quiet[Ratios[{a, b, c}, {1, 1}]]", "Ratios[{a, b, c}, {1, 1}]");
    audit_chk("Ratios#1", "Check[Ratios[{1, 2, 3}, {1, 1}], fired, Ratios::lvldep]", "fired");
    audit_chk("Ratios#1", "Quiet[Differences[{a, b^2, c}, {1, 1}]]",
              "Differences[{a, b^2, c}, {1, 1}]");
    /* Arrays deep enough, and zero-order trailing levels, are unchanged. */
    audit_chk("Ratios#1", "Ratios[{{1, 2, 4}, {3, 6, 9}}, {1, 1}]", "{{1, 3/4}}");
    audit_chk("Ratios#1", "Ratios[{{1, 2, 4}, {3, 6, 9}}, {0, 1}]", "{{2, 2}, {2, 3/2}}");
    audit_chk("Ratios#1", "Ratios[{1, 2, 4}, {1, 0}]", "{2, 2}");
    audit_chk("Ratios#1", "Ratios[{1, 2, 4}, {1}]", "{2, 2}");
    audit_chk("Ratios#1", "Differences[{{1, 2, 4}, {3, 6, 9}}, {1, 1}]", "{{2, 1}}");
}

/* Wave 3 S12/S13: declared argument counts and option checks
 * (src/arity.h, src/arity_table.c).  The evaluator checks a call
 * against its head's declared arity before the C code runs and
 * prints argx/argr/argrx/argt(u)/argb(u)/argm(u), nonopt or optx;
 * the call stays unevaluated. */
static void test_wave3_arity_options(void) {
    audit_chk("Total#3", "Quiet[Check[Total[], fired, Total::argt]]", "fired");
    audit_chk("Total#3", "Quiet[Check[Total[{1, 2}, {1}, 3], fired, Total::nonopt]]", "fired");
    audit_chk("Total#3", "Quiet[Check[Total[{1, 2}, x], fired, Total::level]]", "fired");
    audit_chk("Abs#7", "Quiet[Check[Abs[], fired, Abs::argx]]", "fired");
    audit_chk("Abs#7", "Quiet[Abs[]]", "Abs[]");
    audit_chk("Minus#2", "Quiet[Check[Minus[Sequence[1, 2]], fired, Minus::argx]]", "fired");
    audit_chk("Minus#2", "Quiet[Minus[Sequence[1, 2]]]", "Minus[1, 2]");
    audit_chk("Arg#4", "Quiet[Check[Arg[1, 2], fired, Arg::argx]]", "fired");
    audit_chk("Arg#4", "Quiet[Arg[1, 2]]", "Arg[1, 2]");
    audit_chk("Accumulate#2", "Quiet[Check[Accumulate[], fired, Accumulate::argx]]", "fired");
    audit_chk("Accumulate#2", "Quiet[Accumulate[{1., 2., 3.}, Method -> \"CompensatedSummation\"]]", "{1.0, 3.0, 6.0}");
    audit_chk("Factorial#4", "Quiet[Check[Factorial[1, 2], fired, Factorial::argx]]", "fired");
    audit_chk("Factorial#4", "Quiet[Factorial[1, 2]]", "Factorial[1, 2]");
    audit_chk("FactorialPower#3", "Quiet[Check[FactorialPower[5], fired, FactorialPower::argtu]]", "fired");
    audit_chk("FactorialPower#3", "Quiet[FactorialPower[5]]", "FactorialPower[5]");
    audit_chk("IntervalMemberQ#3", "Quiet[Check[IntervalMemberQ[Interval[{1, 5}]], fired, IntervalMemberQ::argr]]", "fired");
    audit_chk("IntervalMemberQ#3", "Quiet[IntervalMemberQ[Interval[{1, 5}]]]", "IntervalMemberQ[Interval[{1, 5}]]");
    audit_chk("Precision#2", "Quiet[Check[Precision[1., 2.], fired, Precision::argx]]", "fired");
    audit_chk("ReIm#4", "Quiet[Check[ReIm[1, 2], fired, ReIm::argx]]", "fired");
    audit_chk("ReIm#4", "Quiet[ReIm[1, 2]]", "ReIm[1, 2]");
    audit_chk("SetAccuracy#6", "Quiet[Check[SetAccuracy[1.5], fired, SetAccuracy::argr]]", "fired");
    audit_chk("Sign#5", "Quiet[Check[Sign[1, 2], fired, Sign::argx]]", "fired");
    audit_chk("Sign#5", "Quiet[Sign[1, 2]]", "Sign[1, 2]");
}

int main(void) {
    audit_init();
    TEST(test_complex_expand_conjugate);
    TEST(test_differences_int64_min_step);
    TEST(test_plus_infinity_directions);
    TEST(test_power_finite_base_infinity);
    TEST(test_interval_intersection_symbolic);
    TEST(test_ratios_level_depth);
    TEST(test_wave3_arity_options);
    return audit_finish("audit_arithmetic");
}
