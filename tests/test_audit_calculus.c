/* Builtin-audit regressions: calculus (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* D#1: a number or numeric constant is not a differentiation variable. */
static void test_d_numeric_variable(void) {
    audit_chk("D#1", "Check[D[x^2, 2], fired, D::ivar]", "fired");
    audit_chk("D#1", "Quiet[D[x^2, 2]]", "D[x^2, 2]");
    audit_chk("D#1", "Quiet[D[Pi^x, Pi]]", "D[Pi^x, Pi]");
    audit_chk("D#1", "Quiet[D[x^2, I]]", "D[x^2, I]");
    audit_chk("D#1", "Quiet[D[x^2, x, 0]]", "D[x^2, x, 0]");
    audit_chk("D#1", "Quiet[D[x, 1.5]]", "D[x, 1.5]");
    audit_chk("D#1", "Quiet[D[x^2, {2, 1}]]", "D[x^2, {2, 1}]");
    /* ordinary specs unchanged */
    audit_chk("D#1", "D[x^2, {x, 0}]", "x^2");
    audit_chk("D#1", "D[x^2 y, x, y]", "2 x");
    audit_chk("D#1", "D[x^3, {x, 2}]", "6 x");
}

/* Dt#1: Constants -> {a} was read as a differentiation variable. */
static void test_dt_constants_option(void) {
    audit_chk("Dt#1", "Dt[x^2, x, Constants -> {a}]", "2 x");
    audit_chk("Dt#1", "Dt[a x^2, x, Constants -> {a}]", "2 a x");
    audit_chk("Dt#1", "Dt[a x, Constants -> a]", "a Dt[x]");
    audit_chk("Dt#1", "Dt[a x^2 + y, x, Constants -> {a}]", "2 a x + Dt[y, x]");
    audit_chk("Dt#1", "Dt[x^2 y, x]", "2 x y + x^2 Dt[y, x]");
}

/* Integrate#2: an indefinite variable among several specs. */
static void test_integrate_mixed_specs(void) {
    audit_true("Integrate#2", "PossibleZeroQ[Integrate[x, x, y] - x^2 y/2]");
    audit_true("Integrate#2", "PossibleZeroQ[D[Integrate[x, x, y], x, y] - x]");
    audit_true("Integrate#2", "PossibleZeroQ[Integrate[x y, {x, 0, 1}, y] - y^2/4]");
    audit_true("Integrate#2", "PossibleZeroQ[Integrate[x y, x, {y, 0, 2}] - x^2]");
    audit_true("Integrate#2", "PossibleZeroQ[D[Integrate[Sin[x] Cos[y], x, y], x, y] - Sin[x] Cos[y]]");
    audit_chk("Integrate#2", "Integrate[x y, {x, 0, 1}, {y, 0, 2}]", "1");
}

/* Residue#1: a logarithmic branch point has no residue. */
static void test_residue_branch_point(void) {
    audit_chk("Residue#1", "Residue[Log[z]/z, {z, 0}]", "Residue[Log[z]/z, {z, 0}]");
    audit_chk("Residue#1", "Residue[Log[z], {z, 0}]", "Residue[Log[z], {z, 0}]");
    audit_chk("Residue#1", "Residue[Log[z]/z^2, {z, 0}]", "Residue[Log[z]/z^2, {z, 0}]");
    audit_chk("Residue#1", "Residue[Log[1 + z]/z^2, {z, 0}]", "1");
    audit_chk("Residue#1", "Residue[Exp[z]/z^2, {z, 0}]", "1");
    audit_chk("Residue#1", "Residue[1/(z^2 + 1), {z, I}]", "-1/2*I");
}

/* LaplaceTransform#1: the transform integrates from 0-, so DiracDelta[t] counts in full. */
static void test_laplace_delta_origin(void) {
    audit_chk("LaplaceTransform#1", "LaplaceTransform[DiracDelta[t], t, s]", "1");
    audit_chk("LaplaceTransform#1", "LaplaceTransform[DiracDelta[2 t], t, s]", "1/2");
    audit_chk("LaplaceTransform#1", "LaplaceTransform[Cos[t] DiracDelta[t], t, s]", "1");
    audit_chk("LaplaceTransform#1", "LaplaceTransform[E^(3 t) DiracDelta[t], t, s]", "1");
    audit_true("LaplaceTransform#1",
               "PossibleZeroQ[LaplaceTransform[3 DiracDelta[t] + t, t, s] - (3 + 1/s^2)]");
    audit_chk("LaplaceTransform#1", "LaplaceTransform[DiracDelta[t - 2], t, s]", "E^(-2 s)");
}

/* Maximize#1: x == I has no real solution. */
static void test_maximize_complex_equality(void) {
    audit_chk("Maximize#1", "Check[Maximize[{-x^2, x == I}, x], fired, Maximize::infeas]", "fired");
    audit_chk("Maximize#1", "Quiet[Maximize[{-x^2, x == I}, x]]", "{-Infinity, {x -> Indeterminate}}");
    audit_chk("Maximize#1", "Quiet[Minimize[{x^2, x == I}, x]]", "{Infinity, {x -> Indeterminate}}");
    audit_chk("Maximize#1", "Maximize[{-x^2, x == 1}, x]", "{-1, {x -> 1}}");
    audit_chk("Maximize#1", "Minimize[{x^2 + y^2, x == 1}, {x, y}]", "{1, {x -> 1, y -> 0}}");
}

/* PolynomialQuotientRemainder#1: p == Q q + R must hold. */
static void test_pqr_nonpolynomial_terms(void) {
    audit_chk("PolynomialQuotientRemainder#1",
              "PolynomialQuotientRemainder[1/x + x^2, x + 1, x]", "{-1 + x, 1 + 1/x}");
    audit_chk("PolynomialQuotientRemainder#1",
              "PolynomialQuotientRemainder[Sqrt[x] + x^2, x + 1, x]", "{-1 + x, 1 + Sqrt[x]}");
    audit_true("PolynomialQuotientRemainder#1",
               "Module[{qr = PolynomialQuotientRemainder[1/x + x^2, x + 1, x]}, "
               "PossibleZeroQ[qr[[1]] (x + 1) + qr[[2]] - (1/x + x^2)]]");
    audit_chk("PolynomialQuotientRemainder#1",
              "PolynomialQuotientRemainder[x^2 + 1, x + 1, x]", "{-1 + x, 2}");
    audit_chk("PolynomialQuotientRemainder#1",
              "PolynomialQuotientRemainder[a x^2 + 1, x + 1, x]", "{-a + a x, 1 + a}");
}

/* Wave 3 S12/S13: declared argument counts and option checks
 * (src/arity.h, src/arity_table.c).  The evaluator checks a call
 * against its head's declared arity before the C code runs and
 * prints argx/argr/argrx/argt(u)/argb(u)/argm(u), nonopt or optx;
 * the call stays unevaluated. */
static void test_wave3_arity_options(void) {
    audit_chk("SubresultantPolynomialRemainders#3", "Quiet[Check[SubresultantPolynomialRemainders[x^2 - 1, x + 1], fired, SubresultantPolynomialRemainders::argrx]]", "fired");
    audit_chk("SubresultantPolynomialRemainders#3", "Quiet[SubresultantPolynomialRemainders[x^2 - 1, x + 1]]", "SubresultantPolynomialRemainders[-1 + x^2, 1 + x]");
    audit_chk("Limit#3", "Quiet[Check[Limit[x^2, x -> 1, Foo -> 1], fired, Limit::optx]]", "fired");
    audit_chk("Limit#3", "Quiet[Limit[x^2, x -> 1, Foo -> 1]]", "Limit[x^2, x -> 1, Foo -> 1]");
    audit_chk("DSolve#2", "Quiet[Check[DSolve[y'[x] == y[x], y[x], x, Foo -> 1], fired, DSolve::optx]]", "fired");
    audit_chk("Maximize#4", "Quiet[Check[Maximize[-x^2 + 1, x, Foo -> 1], fired, Maximize::optx]]", "fired");
    audit_chk("Maximize#4", "Quiet[Maximize[-x^2 + 1, x, Foo -> 1]]", "Maximize[1 - x^2, x, Foo -> 1]");
    audit_chk("Limit#3", "Limit[x^2, x -> 1, Direction -> \"FromAbove\"]", "1");
}

/* ================================================================== */
/* Wave 3, shared fix S11: Indeterminate / infinity / pole values     */
/* (src/infinity_rules.c and the per-head rules it does not cover).   */
/* ================================================================== */

static void test_s11_special_values(void) {
    /* Cot#1: Cot at an inexact pole */
    audit_true("Cot#1", "Check[Cot[0.], \"m\"] === ComplexInfinity");
    audit_true("Cot#1", "Check[Cot[-0.], \"m\"] === ComplexInfinity");
    audit_true("Cot#1", "Check[Cot[Complex[0., 0.]], \"m\"] === ComplexInfinity");
    audit_true("Cot#1", "Abs[Cot[1.] - 0.6420926159343306] < 10.^-15");
    /* Exp#1: Exp along the imaginary axis / at ComplexInfinity has no limit */
    audit_true("Exp#1", "Exp[ComplexInfinity] === Indeterminate");
    audit_true("Exp#1", "Exp[I Infinity] === Indeterminate");
    audit_true("Exp#1", "N[Exp[Infinity I]] === Indeterminate");
    audit_true("Exp#1", "2^ComplexInfinity === Indeterminate");
    audit_true("Exp#1", "Exp[-Infinity] === 0");
    /* Log#1: Log of every infinity is Infinity */
    audit_true("Log#1", "Log[-Infinity] === Infinity");
    audit_true("Log#1", "Log[I Infinity] === Infinity");
    audit_true("Log#1", "Log[DirectedInfinity[1 + I]] === Infinity");
    /* Log#2: Log[0, 0] and Log[1, 1] are 0/0 forms */
    audit_true("Log#2", "Quiet[Log[0, 0]] === Indeterminate");
    audit_true("Log#2", "Quiet[Log[1, 1]] === Indeterminate");
    audit_true("Log#2", "Check[Log[1, 1], \"fired\"] === \"fired\"");
    audit_true("Log#2", "Log[2, 2] === 1");
    /* Csc#1: Csc at an inexact pole */
    audit_true("Csc#1", "Check[Csc[0.], \"m\"] === ComplexInfinity");
    audit_true("Csc#1", "Check[Csc[0. + 0. I], \"m\"] === ComplexInfinity");
    audit_true("Csc#1", "Check[Csc[0.`20], \"m\"] === ComplexInfinity");
    /* Csch#1: Csch at an inexact pole */
    audit_true("Csch#1", "Check[Csch[0.], \"m\"] === ComplexInfinity");
    audit_true("Csch#1", "Check[Csch[-0.], \"m\"] === ComplexInfinity");
    audit_true("Csch#1", "Check[Csch[0. + 0. I], \"m\"] === ComplexInfinity");
    /* Csch#2: Csch at an imaginary-axis pole */
    audit_true("Csch#2", "Csch[I Pi] === ComplexInfinity");
    audit_true("Csch#2", "Csch[2 I Pi] === ComplexInfinity");
    /* Csch#3: Csch[ComplexInfinity], Csch[Indeterminate] */
    audit_true("Csch#3", "Csch[ComplexInfinity] === Indeterminate");
    audit_true("Csch#3", "Csch[Indeterminate] === Indeterminate");
    /* FractionalPart#1: FractionalPart at the infinities ranges over the unit interval */
    audit_true("FractionalPart#1", "FractionalPart[Infinity] === Interval[{0, 1}]");
    audit_true("FractionalPart#1", "FractionalPart[-Infinity] === Interval[{-1, 0}]");
    audit_true("FractionalPart#1", "FractionalPart[ComplexInfinity] === Indeterminate");
    /* Sech#1: Sech[ComplexInfinity], Sech[Indeterminate] */
    audit_true("Sech#1", "Sech[ComplexInfinity] === Indeterminate");
    audit_true("Sech#1", "Sech[Indeterminate] === Indeterminate");
    audit_true("Sech#1", "Sech[Infinity] === 0");
    /* Cos#4: Cos at the infinities */
    audit_true("Cos#4", "Cos[Infinity] === Interval[{-1, 1}]");
    audit_true("Cos#4", "Cos[ComplexInfinity] === Indeterminate");
    audit_true("Cos#4", "Cos[Indeterminate] === Indeterminate");
    audit_true("Cos#4", "Cos[-Infinity] === Interval[{-1, 1}]");
    /* Cot#4: Cot[Indeterminate] */
    audit_true("Cot#4", "Cot[Indeterminate] === Indeterminate");
    /* Exp#2: Exp of DirectedInfinity[+-1] */
    audit_true("Exp#2", "Exp[DirectedInfinity[1]] === Infinity");
    audit_true("Exp#2", "Exp[DirectedInfinity[-1]] === 0");
    /* Log#5: Log[1, z] is ComplexInfinity */
    audit_true("Log#5", "Quiet[Log[1, 5]] === ComplexInfinity");
    audit_true("Log#5", "Check[Log[1, 5], \"fired\", Power::infy] === \"fired\"");
    /* Log#6: Log of an arbitrary-precision zero */
    audit_true("Log#6", "Check[Log[0``10], \"m\"] === Indeterminate");
    audit_true("Log#6", "Log[0.] === Indeterminate");
    audit_true("Log#6", "Log[0] === -Infinity");
    /* Sin#4: Sin at the infinities */
    audit_true("Sin#4", "Sin[Indeterminate] === Indeterminate");
    audit_true("Sin#4", "Sin[ComplexInfinity] === Indeterminate");
    audit_true("Sin#4", "Sin[Infinity] === Interval[{-1, 1}]");
    audit_true("Sin#4", "Sin[I Infinity] === DirectedInfinity[I]");
    /* Tan#4: Tan[Indeterminate] */
    audit_true("Tan#4", "Tan[Indeterminate] === Indeterminate");
    audit_true("Tan#4", "Tan[Infinity] === Indeterminate");
    /* ArcTan#3: ArcTan[0, 0] */
    audit_true("ArcTan#3", "Check[ArcTan[0, 0], \"fired\", ArcTan::indet] === \"fired\"");
    audit_true("ArcTan#3", "Quiet[ArcTan[0, 0]] === Indeterminate");
    audit_true("ArcTan#3", "ArcTan[0, 1] === Pi/2");
    /* ArcTan#4: two-argument ArcTan with an infinite coordinate */
    audit_true("ArcTan#4", "ArcTan[1, Infinity] === Pi/2");
    audit_true("ArcTan#4", "ArcTan[Infinity, 1] === 0");
    audit_true("ArcTan#4", "ArcTan[-3, -Infinity] === -Pi/2");
    /* Csc#3: Csc[Indeterminate], Csc[ComplexInfinity] */
    audit_true("Csc#3", "Csc[Indeterminate] === Indeterminate");
    audit_true("Csc#3", "Csc[ComplexInfinity] === Indeterminate");
    /* Sec#2: Sec[Indeterminate], Sec[ComplexInfinity] */
    audit_true("Sec#2", "Sec[Indeterminate] === Indeterminate");
    audit_true("Sec#2", "Sec[ComplexInfinity] === Indeterminate");
}

int main(void) {
    audit_init();
    TEST(test_d_numeric_variable);
    TEST(test_dt_constants_option);
    TEST(test_integrate_mixed_specs);
    TEST(test_residue_branch_point);
    TEST(test_laplace_delta_origin);
    TEST(test_maximize_complex_equality);
    TEST(test_pqr_nonpolynomial_terms);
    TEST(test_wave3_arity_options);
    TEST(test_s11_special_values);
    return audit_finish("audit_calculus");
}
