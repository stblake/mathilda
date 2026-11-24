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

int main(void) {
    audit_init();
    TEST(test_unitstep_exact_zero);
    TEST(test_piecewise_values_evaluated);
    TEST(test_wave3_arity_options);
    return audit_finish("audit_elementary_functions");
}
