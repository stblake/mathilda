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

int main(void) {
    audit_init();
    TEST(test_factorial2_negative_odd);
    TEST(test_wave3_arity_options);
    return audit_finish("audit_simplification");
}
