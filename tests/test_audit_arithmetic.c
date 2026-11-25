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
    audit_chk("Total#2", "Total[{1, 2}, Method -> \"CompensatedSummation\"]", "3");
    audit_chk("Total#2", "Total[{{1, 2}, {3, 4}}, 2, Method -> Automatic]", "10");
    audit_chk("Total#2", "Total[Join[{1.}, Table[10.^-16, {10^4}]], Method -> \"CompensatedSummation\"] == 1.000000000001", "True");
    audit_chk("Total#2", "Quiet[Check[Total[{1., 2.}, Method -> \"CompensatedSummation\"], fired]]", "3.0");
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

/* ================================================================== */
/* Wave 3, shared fix S11: Indeterminate / infinity / pole values     */
/* (src/infinity_rules.c and the per-head rules it does not cover).   */
/* ================================================================== */

static void test_s11_special_values(void) {
    /* Abs#2: Abs of the infinities reads the direction; Indeterminate propagates */
    audit_true("Abs#2", "Abs[Infinity] === Infinity");
    audit_true("Abs#2", "Check[Abs[-Infinity], \"m\"] === Infinity");
    audit_true("Abs#2", "Abs[ComplexInfinity] === Infinity");
    audit_true("Abs#2", "Abs[Indeterminate] === Indeterminate");
    audit_true("Abs#2", "Abs[DirectedInfinity[1 + I]] === Infinity");
    /* Minus#1: Minus of a complex-directed infinity flips the direction */
    audit_true("Minus#1", "Minus[DirectedInfinity[I]] === DirectedInfinity[-I]");
    audit_true("Minus#1", "Minus[DirectedInfinity[1 + I]] === DirectedInfinity[-1 - I]");
    audit_true("Minus#1", "Minus[DirectedInfinity[1 + I]] === (-1 - I)/Sqrt[2] Infinity");
    /* Times#3: ComplexInfinity absorbs every finite factor */
    audit_true("Times#3", "x ComplexInfinity === ComplexInfinity");
    audit_true("Times#3", "Quiet[x/0] === ComplexInfinity");
    audit_true("Times#3", "Check[x/0, \"fired\", Power::infy] === \"fired\"");
    audit_true("Times#3", "Quiet[0 x ComplexInfinity] === Indeterminate");
    /* Arg#1: Arg of a directed infinity is the argument of its direction */
    audit_true("Arg#1", "Arg[DirectedInfinity[1 + I]] === Pi/4");
    audit_true("Arg#1", "Arg[Infinity] === 0");
    audit_true("Arg#1", "Check[Arg[-Infinity], \"m\"] === Pi");
    audit_true("Arg#1", "Arg[I Infinity] === Pi/2");
    /* Conjugate#1: Conjugate of the infinities, no spurious message */
    audit_true("Conjugate#1", "Conjugate[Infinity] === Infinity");
    audit_true("Conjugate#1", "Check[Conjugate[-Infinity], \"m\"] === -Infinity");
    audit_true("Conjugate#1", "Check[Conjugate[I Infinity], \"m\"] === DirectedInfinity[-I]");
    audit_true("Conjugate#1", "Conjugate[ComplexInfinity] === ComplexInfinity");
    /* Im#1: Im of a real infinity is 0 */
    audit_true("Im#1", "Im[Infinity] === 0");
    audit_true("Im#1", "Check[Im[-Infinity], \"m\"] === 0");
    audit_true("Im#1", "Im[I Infinity] === Infinity");
    /* Re#3: Re of the infinities and of Indeterminate */
    audit_true("Re#3", "Re[Infinity] === Infinity");
    audit_true("Re#3", "Check[Re[-Infinity], \"m\"] === -Infinity");
    audit_true("Re#3", "Re[Indeterminate] === Indeterminate");
    audit_true("Re#3", "Re[I Infinity] === 0");
    /* Factorial#1: an inexact pole of Factorial is ComplexInfinity, no Overflow[] */
    audit_true("Factorial#1", "Check[Factorial[-1.], \"m\"] === ComplexInfinity");
    audit_true("Factorial#1", "Check[Factorial[-2.], \"m\"] === ComplexInfinity");
    audit_true("Factorial#1", "Check[Factorial[N[-1, 30]], \"m\"] === ComplexInfinity");
    audit_true("Factorial#1", "Abs[Factorial[-0.5] - Sqrt[Pi]] < 10.^-14");
    /* ReIm#2: ReIm of the infinities and of Indeterminate */
    audit_true("ReIm#2", "ReIm[Infinity] === {Infinity, 0}");
    audit_true("ReIm#2", "Check[ReIm[-Infinity], \"m\"] === {-Infinity, 0}");
    audit_true("ReIm#2", "ReIm[Indeterminate] === {Indeterminate, Indeterminate}");
    /* Sign#1: Sign of a directed infinity is its unit direction */
    audit_true("Sign#1", "Sign[Infinity] === 1");
    audit_true("Sign#1", "Sign[-Infinity] === -1");
    audit_true("Sign#1", "Sign[I Infinity] === I");
    audit_true("Sign#1", "Sign[DirectedInfinity[1 + I]] === (1 + I)/Sqrt[2]");
    /* Sign#3: Sign[Indeterminate], Sign[ComplexInfinity] */
    audit_true("Sign#3", "Sign[Indeterminate] === Indeterminate");
    audit_true("Sign#3", "Sign[ComplexInfinity] === Indeterminate");
    /* Power#4: 0^(I b) has no value */
    audit_true("Power#4", "Check[0^I, \"fired\", Power::indet] === \"fired\"");
    audit_true("Power#4", "Quiet[0^I] === Indeterminate");
    audit_true("Power#4", "Quiet[0^(-2 I)] === Indeterminate");
    audit_true("Power#4", "0^2 === 0");
    /* Sqrt#4: Sqrt[-Infinity] is I Infinity */
    audit_true("Sqrt#4", "Sqrt[-Infinity] === DirectedInfinity[I]");
    audit_true("Sqrt#4", "Sqrt[Infinity] === Infinity");
    audit_true("Sqrt#4", "(-Infinity)^(-1/2) === 0");
    /* Arg#3: Arg[Indeterminate] */
    audit_true("Arg#3", "Arg[Indeterminate] === Indeterminate");
    /* Conjugate#2: Conjugate[Indeterminate] */
    audit_true("Conjugate#2", "Conjugate[Indeterminate] === Indeterminate");
    /* Im#3: Im[Indeterminate], Im[ComplexInfinity] */
    audit_true("Im#3", "Im[Indeterminate] === Indeterminate");
    audit_true("Im#3", "Im[ComplexInfinity] === Indeterminate");
    /* Factorial#3: Factorial[Infinity], Factorial[Indeterminate] */
    audit_true("Factorial#3", "Factorial[Infinity] === Infinity");
    audit_true("Factorial#3", "Factorial[Indeterminate] === Indeterminate");
    /* Multinomial#1: an Indeterminate argument */
    audit_true("Multinomial#1", "Multinomial[Indeterminate, 2] === Indeterminate");
    audit_true("Multinomial#1", "Multinomial[Indeterminate, 2, 3] === Indeterminate");
    audit_true("Multinomial#1", "Multinomial[2, 3] === 10");
}

static void test_wave4_arith_misc(void) {
    audit_chk("Sqrt#2", "Sqrt[2*10^40]", "100000000000000000000 Sqrt[2]");
    audit_chk("Sqrt#2", "(2^70*3)^(1/3)", "8388608 6^(1/3)");
    audit_chk("Total#1", "Total[f[1,2,3]]", "6");
    audit_chk("Total#1", "Total[1/2]", "1/2");
    audit_chk("Arg#2", "Arg[E^(I Pi/3)]", "1/3 Pi");
    audit_chk("ComplexExpand#2", "ComplexExpand[Re[f[x]]]", "Re[f[x]]");
    audit_chk("ComplexExpand#2", "ComplexExpand[Im[f[x]]]", "Im[f[x]]");
    audit_chk("ComplexExpand#2", "ComplexExpand[Re[Gamma[x]]]", "Gamma[x]");
    audit_chk("Im#2", "Im[(-1)^(1/3)]", "1/2 Sqrt[3]");
    audit_chk("Im#2", "Im[(-8)^(1/3)]", "Sqrt[3]");
    audit_chk("Im#2", "Im[(-1)^(1/4)]", "1/Sqrt[2]");
    audit_chk("Re#1", "Re[(-1)^(1/3)]", "1/2");
    audit_chk("Re#1", "Re[E^(I Pi/4)]", "1/Sqrt[2]");
    audit_chk("Re#1", "Re[(-8)^(1/3)]", "1");
    audit_chk("Re#2", "Re[1/(1 + I Sqrt[3])]", "1/4");
    audit_chk("Re#2", "Re[Exp[1 + I]]", "E Cos[1]");
    audit_chk("ReIm#1", "ReIm[x]", "{Re[x], Im[x]}");
    audit_chk("ReIm#1", "ReIm[{1+I, 2, x}]", "{{1, 1}, {2, 0}, {Re[x], Im[x]}}");
    audit_chk("ReIm#3", "ReIm[(-1)^(1/3)]", "{1/2, 1/2 Sqrt[3]}");
    audit_chk("ReIm#3", "ReIm[Exp[I Pi/3]]", "{1/2, 1/2 Sqrt[3]}");
    audit_chk("Quotient#2", "Quotient[Pi, 1]", "3");
    audit_chk("Quotient#2", "Quotient[Sqrt[2], 1/10]", "14");
    audit_chk("QuotientRemainder#2", "QuotientRemainder[10, Pi]", "{3, 10 - 3 Pi}");
    audit_chk("Sign#2", "Sign[-x]", "-Sign[x]");
    audit_chk("Sign#2", "Sign[2 x]", "Sign[x]");
    audit_chk("Sign#2", "Sign[-3 x]", "-Sign[x]");
    audit_chk("Sign#2", "Sign[Pi x]", "Sign[x]");
    audit_chk("Factorial#2", "D[x!, x]", "Gamma[1 + x] PolyGamma[0, 1 + x]");
    audit_chk("FactorialPower#1", "FactorialPower[1/2, 33] == FactorialPower[1/2, 32] (1/2 - 32)", "True");
    audit_chk("FactorialPower#1", "Head[FactorialPower[2.5, 33]]", "Real");
    audit_chk("Accumulate#1", "Accumulate[{10^20, 1, -10^20}, Method->\"CompensatedSummation\"]", "{100000000000000000000, 100000000000000000001, 1}");
    audit_chk("Accumulate#1", "Accumulate[{1,2,3}, Method->\"CompensatedSummation\"]", "{1, 3, 6}");
    audit_chk("ListGradient#1", "Check[ListGradient[{{1,2},{3,4,5}}], \"fired\", ListGradient::rect]", "\"fired\"");
    audit_chk("ListGradient#1", "Quiet[ListGradient[{{1,2},{3,4,5}}]]", "ListGradient[{{1, 2}, {3, 4, 5}}]");
    audit_chk("MantissaExponent#1", "MantissaExponent[Pi]", "{1/10 Pi, 1}");
    audit_chk("DigitCount#1", "DigitCount[1.*^-20]", "{1, 0, 0, 0, 0, 0, 0, 0, 0, 15}");
    audit_chk("DigitCount#1", "RealDigits[1.*^-20]", "{{1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, -19}");
    audit_chk("DigitCount#1", "RealDigits[0.1 + 0.2]", "{{3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 0}");
    audit_chk("IntervalIntersection#2", "IntervalIntersection[Interval[{1,2}], Interval[]]", "Interval[]");
    audit_chk("IntervalIntersection#2", "IntervalIntersection[Interval[], Interval[]]", "Interval[]");
    audit_chk("IntervalMemberQ#1", "IntervalMemberQ[Interval[], 1]", "False");
    audit_chk("IntervalMemberQ#2", "IntervalMemberQ[Interval[{1, 5}], 2 + I]", "False");
    audit_chk("IntervalUnion#1", "IntervalUnion[Interval[], Interval[{1,2}]]", "Interval[{1, 2}]");
    audit_chk("IntervalUnion#1", "IntervalUnion[Interval[]]", "Interval[]");
    audit_chk("IntervalUnion#1", "IntervalUnion[IntervalIntersection[Interval[{1,2}], Interval[{3,4}]], Interval[{5,6}]]", "Interval[{5, 6}]");
    audit_chk("Interval#1", "Interval[{1, 2}]^Infinity", "Interval[{1, 2}]^Infinity");
    audit_chk("NValues#1", "N[nvw4f[x_], p_] := p; N[nvw4f[1], 20]", "20.000000000000000000");
    audit_chk("NValues#2", "N[nvw4a] = 3; N[nvw4a] =.; NValues[nvw4a]", "{}");
    audit_true("SetAccuracy#2", "Min[List @@ Map[Accuracy, SetAccuracy[f[1000.5, 0.001], 20]]] >= 20");
    audit_chk("SetAccuracy#2", "Map[Accuracy, SetAccuracy[f[1000.5, 0.001], 20]][[2]]", "20.0");
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
    TEST(test_s11_special_values);
    TEST(test_wave4_arith_misc);
    return audit_finish("audit_arithmetic");
}
