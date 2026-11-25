/* Builtin-audit regressions: elementary functions (review/issues/<Head>.md; tests/test_audit.h). */
#include "test_audit.h"

/* UnitStep#1 (high): "UnitStep gives 0 for an exact zero argument, because
 * numerical noise is taken as a certified negative sign".  The sign comes
 * from the certified sign Positive uses; an exact zero is >= 0. */
static void test_unitstep_exact_zero(void) {
    audit_chk("UnitStep#1", "UnitStep[Sqrt[5 + 2 Sqrt[6]] - Sqrt[2] - Sqrt[3]]", "1");
    audit_chk("UnitStep#1", "UnitStep[Sin[Pi/7]^2 + Cos[Pi/7]^2 - 1]", "1");
    audit_chk("UnitStep#1", "UnitStep[Log[2] + Log[3] - Log[6]]", "1");
    audit_chk("UnitStep#1", "UnitStep[Pi - 3, Sqrt[5 + 2 Sqrt[6]] - Sqrt[2] - Sqrt[3]]", "1");
    /* Genuine signs, and the machine path. */
    audit_chk("UnitStep#1", "UnitStep[Pi - 3]", "1");
    audit_chk("UnitStep#1", "UnitStep[3 - Pi]", "0");
    audit_chk("UnitStep#1", "UnitStep[Pi - 80143857/25510582]", "1");
    audit_chk("UnitStep#1", "UnitStep[1 - Zeta[3]]", "0");
    audit_chk("UnitStep#1", "UnitStep[{-1., 0., 2.}]", "{0, 1, 1}");
    audit_chk("UnitStep#1", "UnitStep[x]", "UnitStep[x]");
    audit_chk("UnitStep#1", "UnitStep[I]", "UnitStep[I]");
}

/* Piecewise#1 (high): "Values and default of a symbolic Piecewise are never
 * evaluated".  The values and the default kept in the returned form are
 * evaluated; a side effect in one runs once. */
static void test_piecewise_values_evaluated(void) {
    audit_chk("Piecewise#1", "FullForm[Piecewise[{{1 + 1, x > 0}}, 2 + 2]]",
              "Piecewise[List[List[2, Greater[x, 0]]], 4]");
    audit_chk("Piecewise#1", "Piecewise[{{1 + 1, x > 0}, {Sin[0], x < -1}}]",
              "Piecewise[{{2, x > 0}, {0, x < -1}}, 0]");
    audit_run("auditPwN = 0");
    audit_run("auditPw = Piecewise[{{auditPwN++; 7, x > 0}}, 1]");
    audit_chk("Piecewise#1", "auditPw", "Piecewise[{{7, x > 0}}, 1]");
    audit_chk("Piecewise#1", "auditPwN", "1");
    /* A clause that is decided still yields its value alone. */
    audit_chk("Piecewise#1", "Piecewise[{{1 + 1, 1 > 0}}, 2 + 2]", "2");
    audit_chk("Piecewise#1", "Piecewise[{{1 + 1, 1 < 0}}, 2 + 2]", "4");
    audit_chk("Piecewise#1", "Piecewise[{{x^2, x > 0}}, -x] /. x -> 3", "9");
}

/* Wave 3 S12/S13: declared argument counts and option checks
 * (src/arity.h, src/arity_table.c).  The evaluator checks a call
 * against its head's declared arity before the C code runs and
 * prints argx/argr/argrx/argt(u)/argb(u)/argm(u), nonopt or optx;
 * the call stays unevaluated. */
static void test_wave3_arity_options(void) {
    audit_chk("TrigReduce#3", "Quiet[Check[TrigReduce[], fired, TrigReduce::argx]]", "fired");
    audit_chk("TrigReduce#3", "Quiet[TrigReduce[]]", "TrigReduce[]");
    audit_chk("Clip#3", "Quiet[Check[Clip[1, {1, 2}, {3, 4}, 5], fired, Clip::argb]]", "fired");
    audit_chk("Clip#3", "Quiet[Clip[1, {1, 2}, {3, 4}, 5]]", "Clip[1, {1, 2}, {3, 4}, 5]");
    audit_chk("Ramp#4", "Quiet[Check[Ramp[1, 2], fired, Ramp::argx]]", "fired");
    audit_chk("Ramp#4", "Quiet[Ramp[1, 2]]", "Ramp[1, 2]");
    audit_chk("Rationalize#1", "Quiet[Check[Rationalize[1, 2, 3], fired, Rationalize::argt]]", "fired");
    audit_chk("Rationalize#1", "Quiet[Rationalize[1, 2, 3]]", "Rationalize[1, 2, 3]");
}

/* ================================================================== */
/* Wave 3, shared fix S11: Indeterminate / infinity / pole values     */
/* (src/infinity_rules.c and the per-head rules it does not cover).   */
/* ================================================================== */

static void test_s11_special_values(void) {
    /* ArcCot#2: ArcCot at an inexact pole agrees with the exact one */
    audit_true("ArcCot#2", "Check[ArcCot[1.0 I], \"m\"] === ArcCot[I]");
    audit_true("ArcCot#2", "Check[ArcCot[N[1, 30] I], \"m\"] === -I Infinity");
    /* ArcCsc#1: ArcCsc[0] */
    audit_true("ArcCsc#1", "ArcCsc[0] === ComplexInfinity");
    audit_true("ArcCsc#1", "ArcCsc[0.] === ComplexInfinity");
    /* ArcSec#1: ArcSec[0] */
    audit_true("ArcSec#1", "ArcSec[0] === ComplexInfinity");
    /* ArcSec#2: ArcSec of an arbitrary-precision zero */
    audit_true("ArcSec#2", "Check[ArcSec[0.0`20], \"m\"] === ComplexInfinity");
    audit_true("ArcSec#2", "ArcSec[0.] === ComplexInfinity");
    /* ArcTanh#1: ArcTanh at an inexact +-1 is the signed infinity of the exact call */
    audit_true("ArcTanh#1", "Check[ArcTanh[1.], \"m\"] === Infinity");
    audit_true("ArcTanh#1", "Check[ArcTanh[-1.], \"m\"] === -Infinity");
    audit_true("ArcTanh#1", "Check[ArcTanh[N[1, 30]], \"m\"] === Infinity");
    audit_true("ArcTanh#1", "Abs[ArcTanh[0.5] - 0.5493061443340549] < 10.^-15");
    /* Cosh#1: Cosh[Indeterminate], Cosh[ComplexInfinity] */
    audit_true("Cosh#1", "Cosh[Indeterminate] === Indeterminate");
    audit_true("Cosh#1", "Cosh[ComplexInfinity] === Indeterminate");
    audit_true("Cosh#1", "Cosh[-Infinity] === Infinity");
    /* Log10#1: Log10[-Infinity] */
    audit_true("Log10#1", "Log10[-Infinity] === Infinity");
    audit_true("Log10#1", "Log10[1000] === 3");
    /* Log2#1: Log2 of the infinities */
    audit_true("Log2#1", "Log2[-Infinity] === Infinity");
    audit_true("Log2#1", "Log2[I Infinity] === Infinity");
    /* Ceiling#2: Ceiling[Indeterminate] */
    audit_true("Ceiling#2", "Ceiling[Indeterminate] === Indeterminate");
    /* Floor#2: Floor[Indeterminate] */
    audit_true("Floor#2", "Floor[Indeterminate] === Indeterminate");
    /* ArcTanh#2: ArcTanh[Indeterminate] */
    audit_true("ArcTanh#2", "ArcTanh[Indeterminate] === Indeterminate");
    /* Log10#3: Log10 of an arbitrary-precision zero */
    audit_true("Log10#3", "Check[Log10[0``10], \"m\"] === Indeterminate");
    audit_true("Log10#3", "Check[Log2[0``10], \"m\"] === Indeterminate");
    /* Ramp#3: Ramp[Indeterminate] */
    audit_true("Ramp#3", "Ramp[Indeterminate] === Indeterminate");
    /* Sinh#1: Sinh[Indeterminate] */
    audit_true("Sinh#1", "Sinh[Indeterminate] === Indeterminate");
    /* Sinh#2: Sinh[ComplexInfinity] */
    audit_true("Sinh#2", "Sinh[ComplexInfinity] === Indeterminate");
    /* IntegerPart#2: IntegerPart[Indeterminate] (category bitwise) */
    audit_true("IntegerPart#2", "IntegerPart[Indeterminate] === Indeterminate");
}

/* ================================================================== */
/* Wave 4 G1-specfun: elementary-functions medium findings.            */
/* ================================================================== */

static void test_wave4_elementary(void) {
    /* Floor/Ceiling of exact complex numbers act on both parts */
    audit_true("Ceiling#1", "Ceiling[I Sqrt[2]] === 2 I");
    audit_true("Ceiling#1", "Ceiling[Pi + I E] === 4 + 3 I");
    audit_true("Ceiling#1", "Ceiling[Sqrt[2] + I] === 2 + I");
    audit_true("Floor#1", "Floor[Sqrt[-2]] === I");
    audit_true("Floor#1", "Floor[Pi + E I] === 3 + 2 I");
    audit_true("Floor#1", "Floor[-1/2 + Pi I] === -1 + 3 I");
    /* a trailing {v, True} clause is the default */
    audit_chk("Piecewise#2", "FullForm[Piecewise[{{a, x > 0}, {b, True}}]]",
              "Piecewise[List[List[a, Greater[x, 0]]], b]");
    /* TrigExpand expands reciprocal-function denominators */
    audit_true("TrigExpand#1", "TrigExpand[Sec[2 x]] === 1/(Cos[x]^2 - Sin[x]^2)");
    audit_true("TrigExpand#1", "TrigExpand[Csc[x + y]] === 1/(Cos[y] Sin[x] + Cos[x] Sin[y])");
    audit_true("TrigExpand#1", "TrigExpand[Sech[2 x]] === 1/(Cosh[x]^2 + Sinh[x]^2)");
    audit_true("TrigExpand#2", "TrigExpand[Tan[2 x]] === (2 Cos[x] Sin[x])/(Cos[x]^2 - Sin[x]^2)");
    audit_true("TrigExpand#2", "FreeQ[TrigExpand[{Cot[2 x], Tanh[x + y]}], Csc | Sec | Sech]");
    audit_true("TrigExpand#2", "PossibleZeroQ[TrigExpand[Tanh[x + y]] - Tanh[x + y]]");
    /* TrigReduce: negated collapse orientation; held code untouched */
    audit_true("TrigReduce#1", "TrigReduce[{Tan[x] + Cot[y], Tanh[x] - Coth[y]}] === "
                               "{Cos[x - y] Csc[y] Sec[x], -Cosh[x - y] Csch[y] Sech[x]}");
    audit_true("TrigReduce#2", "TrigReduce[Hold[Sin[x]^2]] === Hold[Sin[x]^2]");
    /* ArcCsc/ArcSec invert their own table values */
    audit_true("ArcCsc#2", "{ArcCsc[Csc[Pi/10]], ArcCsc[Csc[3 Pi/10]], ArcCsc[Csc[Pi/12]], "
                           "ArcCsc[Csc[5 Pi/12]], ArcCsc[Csc[Pi/5]]} === {Pi/10, 3 Pi/10, Pi/12, 5 Pi/12, Pi/5}");
    /* Clip: exact comparison inside double resolution */
    audit_true("Clip#1", "Clip[10^20 + 1, {0, 10^20}] === 10^20");
    audit_true("Clip#2", "Clip[1/3 + 10^-30, {0, 1/3}] === 1/3");
    audit_true("Clip#2", "Clip[Sqrt[2] + 10^-20, {0, Sqrt[2]}] === Sqrt[2]");
    audit_true("Clip#2", "Clip[Pi, {0, 314159265358979323846/10^20}] === 157079632679489661923/50000000000000000000");
    /* ExpToTrig */
    audit_true("ExpToTrig#1", "ExpToTrig[Exp[I Pi/3]] == 1/2 + I Sqrt[3]/2 && FreeQ[ExpToTrig[Exp[I Pi/3]], Power[-1, _]]");
    audit_true("ExpToTrig#2", "ExpToTrig[(x^2 - 1)/(x - 1)] === (x^2 - 1)/(x - 1)");
    audit_true("ExpToTrig#2", "ExpToTrig[1/(x - 1) + 1/(x + 1)] === 1/(x - 1) + 1/(x + 1)");
    audit_true("ExpToTrig#2", "ExpToTrig[(x + 1)^2] === (x + 1)^2");
    audit_true("ExpToTrig#3", "Chop[N[ExpToTrig[Log[a + I b] - Log[a - I b]] - (Log[a + I b] - Log[a - I b]) /. {a -> -1, b -> 1}]] == 0");
    audit_true("ExpToTrig#3", "Chop[N[ExpToTrig[Log[a + b] - Log[a - b]] - (Log[a + b] - Log[a - b]) /. {a -> -1, b -> 3}]] == 0");
    audit_true("ExpToTrig#3", "ExpToTrig[Log[1 + I b] - Log[1 - I b]] === 2 I ArcTan[b]");
    /* Ramp: derivative and definite integral */
    audit_true("Ramp#1", "(D[Ramp[x], x] /. x -> 1) === 1");
    audit_true("Ramp#1", "(D[Ramp[x], x] /. x -> -1) === 0");
    audit_true("Ramp#2", "Integrate[Ramp[x], {x, -1, 2}] === 2");
    /* TrigFactor: sum-to-product, 1 - Sin^2, bivariate degree 10 */
    audit_true("TrigFactor#1", "TrigFactor[Sin[x] + Sin[y]] === 2 Cos[x/2 - y/2] Sin[x/2 + y/2]");
    audit_true("TrigFactor#1", "TrigFactor[Cos[x] + Cos[y]] === 2 Cos[x/2 - y/2] Cos[x/2 + y/2]");
    audit_true("TrigFactor#2", "TrigFactor[Sec[x]^2 - Tan[x]^2] === 1");
    audit_true("TrigFactor#2", "TrigFactor[Csc[x]^2 - Cot[x]^2] === 1");
    audit_true("TrigFactor#3", "Head[TrigFactor[10*Cos[x]^9*Sin[x] - 120*Cos[x]^7*Sin[x]^3 + "
                               "252*Cos[x]^5*Sin[x]^5 - 120*Cos[x]^3*Sin[x]^7 + 10*Cos[x]*Sin[x]^9]] === Times");
}

int main(void) {
    audit_init();
    TEST(test_unitstep_exact_zero);
    TEST(test_piecewise_values_evaluated);
    TEST(test_wave3_arity_options);
    TEST(test_s11_special_values);
    TEST(test_wave4_elementary);
    return audit_finish("audit_elementary_functions");
}
