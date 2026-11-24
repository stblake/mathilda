/* Builtin-audit regressions: simplification (review/issues/<Head>.md; tests/test_audit.h). */
#include "test_audit.h"

/* Factorial2#1 (high): "Factorial2 of a negative odd integer below -1
 * returns ComplexInfinity instead of the finite value".  The recurrence
 * n!! = n (n-2)!! run down from (-1)!! = 1 gives (-2k-1)!! = (-1)^k/(2k-1)!!;
 * only negative even integers are poles. */
static void test_factorial2_negative_odd(void) {
    audit_chk("Factorial2#1", "Factorial2[-3]", "-1");
    audit_chk("Factorial2#1", "Factorial2[-5]", "1/3");
    audit_chk("Factorial2#1", "Factorial2[-7]", "-1/15");
    audit_chk("Factorial2#1", "Factorial2[-9]", "1/105");
    audit_true("Factorial2#1",
               "And @@ Table[Factorial2[n] == n Factorial2[n - 2], {n, -41, 15, 2}]");
    audit_true("Factorial2#1", "Abs[Factorial2[-41] - Factorial2[-41.]] < 10^-12 Abs[Factorial2[-41]]");
    audit_true("Factorial2#1", "Abs[N[Factorial2[-5]] - Factorial2[-5.]] < 10^-14");
    /* Poles and the non-negative values are unchanged. */
    audit_chk("Factorial2#1", "Factorial2[-2]", "ComplexInfinity");
    audit_chk("Factorial2#1", "Factorial2[-10]", "ComplexInfinity");
    audit_chk("Factorial2#1", "Factorial2[-1]", "1");
    audit_chk("Factorial2#1", "Factorial2[7]", "105");
    audit_chk("Factorial2#1", "Factorial2[-10^30]", "ComplexInfinity");
    audit_chk("Factorial2#1", "Factorial2[-10^30 - 1]", "Factorial2[-1000000000000000000000000000001]");
}

/* Wave 3 S12/S13: declared argument counts and option checks
 * (src/arity.h, src/arity_table.c).  The evaluator checks a call
 * against its head's declared arity before the C code runs and
 * prints argx/argr/argrx/argt(u)/argb(u)/argm(u), nonopt or optx;
 * the call stays unevaluated. */
static void test_wave3_arity_options(void) {
    audit_chk("Refine#4", "Quiet[Check[Refine[], fired, Refine::argt]]", "fired");
    audit_chk("Refine#4", "Quiet[Refine[]]", "Refine[]");
    audit_chk("Assuming#1", "Quiet[Check[Assuming[x > 0], fired, Assuming::argr]]", "fired");
    audit_chk("Assuming#1", "Quiet[Assuming[x > 0]]", "Assuming[x > 0]");
}

/* ================================================================== */
/* Wave 3, shared fix S11: Indeterminate / infinity / pole values     */
/* (src/infinity_rules.c and the per-head rules it does not cover).   */
/* ================================================================== */

static void test_s11_special_values(void) {
    /* ArcCos#1: ArcCos at +-Infinity (branch of ArcCos[10.^8] = 0 + 19.1 I) */
    audit_true("ArcCos#1", "ArcCos[Infinity] === DirectedInfinity[I]");
    audit_true("ArcCos#1", "ArcCos[-Infinity] === DirectedInfinity[-I]");
    audit_true("ArcCos#1", "Im[ArcCos[10.^8]] > 19");
    audit_true("ArcCos#1", "Im[ArcCos[-10.^8]] < -19");
    /* ArcCoth#1: ArcCoth[0] is the principal I Pi/2 */
    audit_true("ArcCoth#1", "ArcCoth[0] === I Pi/2");
    audit_true("ArcCoth#1", "Abs[ArcCoth[0.] - N[I Pi/2]] < 10.^-15");
    /* ArcCoth#2: ArcCoth at an inexact pole */
    audit_true("ArcCoth#2", "Check[ArcCoth[1.], \"m\"] === Infinity");
    audit_true("ArcCoth#2", "ArcCoth[1] === Infinity");
    audit_true("ArcCoth#2", "Check[ArcCoth[-1.], \"m\"] === -Infinity");
    /* Coth#1: Coth at an inexact pole */
    audit_true("Coth#1", "Check[Coth[0.], \"m\"] === ComplexInfinity");
    audit_true("Coth#1", "Check[Coth[0. + 0. I], \"m\"] === ComplexInfinity");
    /* Coth#2: Coth at an imaginary-axis pole */
    audit_true("Coth#2", "Coth[I Pi] === ComplexInfinity");
    /* Tanh#1: Tanh at an imaginary-axis pole */
    audit_true("Tanh#1", "Tanh[I Pi/2] === ComplexInfinity");
    audit_true("Tanh#1", "Tanh[3 I Pi/2] === ComplexInfinity");
    /* ArcCos#2: ArcCos along the imaginary axis */
    audit_true("ArcCos#2", "ArcCos[I Infinity] === DirectedInfinity[-I]");
    audit_true("ArcCos#2", "ArcCos[-I Infinity] === DirectedInfinity[I]");
    audit_true("ArcCos#2", "Im[ArcCos[10.^8 I]] < -19");
    /* ArcCos#3: ArcCos[Indeterminate] */
    audit_true("ArcCos#3", "ArcCos[Indeterminate] === Indeterminate");
    /* Coth#3: Coth[ComplexInfinity], Coth[Indeterminate] */
    audit_true("Coth#3", "Coth[ComplexInfinity] === Indeterminate");
    audit_true("Coth#3", "Coth[Indeterminate] === Indeterminate");
    audit_true("Coth#3", "Coth[-Infinity] === -1");
    /* Element#3: an infinity is in no numeric domain */
    audit_true("Element#3", "Element[Infinity, Reals] === False");
    audit_true("Element#3", "Element[ComplexInfinity, Complexes] === False");
    audit_true("Element#3", "Element[I Infinity, Complexes] === False");
    audit_true("Element#3", "Element[2, Reals] === True");
    /* Factorial2#2: Factorial2 at an inexact pole */
    audit_true("Factorial2#2", "Check[Factorial2[-2.], \"m\"] === ComplexInfinity");
    audit_true("Factorial2#2", "Check[Factorial2[-4.], \"m\"] === ComplexInfinity");
    audit_true("Factorial2#2", "Factorial2[-2] === ComplexInfinity");
    /* Tanh#2: Tanh[Indeterminate], Tanh[ComplexInfinity] */
    audit_true("Tanh#2", "Tanh[Indeterminate] === Indeterminate");
    audit_true("Tanh#2", "Tanh[ComplexInfinity] === Indeterminate");
    audit_true("Tanh#2", "Tanh[Infinity] === 1");
}

int main(void) {
    audit_init();
    TEST(test_factorial2_negative_odd);
    TEST(test_wave3_arity_options);
    TEST(test_s11_special_values);
    return audit_finish("audit_simplification");
}
