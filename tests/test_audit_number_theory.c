/* Builtin-audit regressions: number theory (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* PrimitiveRootList#1 (high): the root buffer, sizeof(mpz_t) * (phi + 1),
 * wrapped for phi near 2^60 (80 bytes for n = 3*2^60 + 5, a prime) and the
 * enumeration overran it.  The buffer is now sized by the exact count
 * phi(phi(n)), checked: an impossible list aborts with General::nomem. */
static void test_primitive_root_list_huge_phi(void) {
    audit_chk("PrimitiveRootList#1",
              "Quiet[CheckAbort[TimeConstrained[PrimitiveRootList[3458764513820540933], 10], aborted]]",
              "aborted");
    audit_chk("PrimitiveRootList#1",
              "Quiet[Check[CheckAbort[PrimitiveRootList[3458764513820540933], 0], fired, General::nomem]]",
              "fired");
    audit_chk("PrimitiveRootList#1", "PrimitiveRootList[13]", "{2, 6, 7, 11}");
    audit_chk("PrimitiveRootList#1", "PrimitiveRootList[54]", "{5, 11, 23, 29, 41, 47}");
    audit_true("PrimitiveRootList#1",
               "Length[PrimitiveRootList[1000003]] == EulerPhi[EulerPhi[1000003]]");
    audit_true("PrimitiveRootList#1",
               "Length[PrimitiveRootList[2 * 3^7]] == EulerPhi[EulerPhi[2 * 3^7]]");
}

/* Builtin-audit regressions: number theory (review/issues/<Head>.md; tests/test_audit.h). */
#include "test_audit.h"

/* PowerMod#1, #2, #3 (high): PowerMod[a, 1/r, m] falsely reported "no
 * integer solutions" for square roots mod large powers of 2 (#1), for a
 * divisible by p (#2) and for cube roots mod a large p = 1 (mod 3) (#3).
 * Every returned root is checked by raising it back; small moduli are
 * checked against brute force in both directions. */
static void test_powermod_roots(void) {
    audit_chk("PowerMod#1", "Mod[PowerMod[9, 1/2, 2^21]^2 - 9, 2^21]", "0");
    audit_chk("PowerMod#1", "Mod[PowerMod[4, 1/2, 2^30]^2 - 4, 2^30]", "0");
    audit_chk("PowerMod#1", "Mod[PowerMod[17, 1/2, 2^30]^2 - 17, 2^30]", "0");
    audit_chk("PowerMod#1", "Mod[PowerMod[3^5, 1/5, 2^40]^5 - 3^5, 2^40]", "0");
    audit_chk("PowerMod#1", "Check[PowerMod[3, 1/2, 2^30], fired, PowerMod::root]", "fired");
    audit_chk("PowerMod#2", "Mod[PowerMod[1009^2, 1/2, 1009^3]^2 - 1009^2, 1009^3]", "0");
    audit_chk("PowerMod#2", "Mod[PowerMod[4*1009^2, 1/2, 1009^3]^2 - 4*1009^2, 1009^3]", "0");
    audit_chk("PowerMod#2", "Check[PowerMod[1009, 1/2, 1009^3], fired, PowerMod::root]", "fired");
    audit_chk("PowerMod#3", "PowerMod[1, 1/3, 1000003]", "1");
    audit_chk("PowerMod#3", "Mod[PowerMod[125, 1/3, 1000003]^3 - 125, 1000003]", "0");
    audit_chk("PowerMod#3", "Mod[PowerMod[8, 1/3, 1000003*7]^3 - 8, 1000003*7]", "0");
    audit_chk("PowerMod#3", "Check[PowerMod[2, 1/3, 1000003], fired, PowerMod::root]", "fired");
    /* Random roots of random powers over several modulus shapes. */
    audit_run("SeedRandom[11]; auditPmBad = 0; "
              "Do[With[{x = RandomInteger[{0, m - 1}], r = RandomChoice[{2, 3, 4, 6, 8, 9, 16}]}, "
              "With[{a = PowerMod[x, r, m]}, With[{pm = Quiet[PowerMod[a, 1/r, m]]}, "
              "If[!IntegerQ[pm] || PowerMod[pm, r, m] != a, auditPmBad++]]]], "
              "{m, {2^61, 3^20, 1009^3, 1000003^2, 2^10 3^7 1000003, 998244353, 7^9 2^20}}, {20}]");
    audit_chk("PowerMod#3", "auditPmBad", "0");
    /* Brute force: a root is returned exactly when one exists. */
    audit_run("auditPmBrute = 0; "
              "Do[With[{sol = Select[Range[0, m - 1], Mod[#^r - a, m] == 0 &, 1], "
              "pm = Quiet[PowerMod[a, 1/r, m]]}, "
              "If[sol === {}, If[IntegerQ[pm], auditPmBrute++], "
              "If[!IntegerQ[pm] || Mod[pm^r - a, m] != 0, auditPmBrute++]]], "
              "{m, 2, 40}, {r, 2, 5}, {a, 0, m - 1}]");
    audit_chk("PowerMod#1", "auditPmBrute", "0");
    /* Integer exponents unchanged. */
    audit_chk("PowerMod#1", "PowerMod[3, 100, 7]", "4");
    audit_chk("PowerMod#1", "PowerMod[3, -1, 7]", "5");
}

/* ContinuedFraction#1 (high): "ContinuedFraction of an inexact real drops
 * the last term when the expansion terminates (0.1 gives {0})".  A value
 * integral to the input's precision ends the expansion with that integer. */
static void test_continued_fraction_inexact_last_term(void) {
    audit_chk("ContinuedFraction#1", "ContinuedFraction[0.1]", "{0, 10}");
    audit_chk("ContinuedFraction#1", "ContinuedFraction[1/3.]", "{0, 3}");
    audit_chk("ContinuedFraction#1", "ContinuedFraction[0.7]", "{0, 1, 2, 3}");
    audit_chk("ContinuedFraction#1", "ContinuedFraction[N[1/10, 30]]", "{0, 10}");
    audit_chk("ContinuedFraction#1", "ContinuedFraction[0.1, 3]", "{0, 10}");
    audit_chk("ContinuedFraction#1", "ContinuedFraction[1.1]", "{1, 10}");
    audit_chk("ContinuedFraction#1", "ContinuedFraction[10.1]", "{10, 10}");
    audit_chk("ContinuedFraction#1", "ContinuedFraction[0.01]", "{0, 100}");
    audit_true("ContinuedFraction#1", "FromContinuedFraction[ContinuedFraction[0.7]] == 7/10");
    /* Expansions that exhaust the precision are unchanged. */
    audit_chk("ContinuedFraction#1", "ContinuedFraction[N[Pi]]",
              "{3, 7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14}");
    audit_chk("ContinuedFraction#1", "ContinuedFraction[0.75]", "{0, 1, 3}");
    audit_chk("ContinuedFraction#1", "ContinuedFraction[2.]", "{2}");
}

/* Wave 3 S12/S13: declared argument counts and option checks
 * (src/arity.h, src/arity_table.c).  The evaluator checks a call
 * against its head's declared arity before the C code runs and
 * prints argx/argr/argrx/argt(u)/argb(u)/argm(u), nonopt or optx;
 * the call stays unevaluated. */
static void test_wave3_arity_options(void) {
    audit_chk("FromContinuedFraction#3", "Quiet[Check[FromContinuedFraction[{1}, {2}], fired, FromContinuedFraction::argx]]", "fired");
    audit_chk("FromContinuedFraction#3", "Quiet[FromContinuedFraction[{1}, {2}]]", "FromContinuedFraction[{1}, {2}]");
}

/* ================================================================== */
/* Wave 3, shared fix S11: Indeterminate / infinity / pole values     */
/* (src/infinity_rules.c and the per-head rules it does not cover).   */
/* ================================================================== */

static void test_s11_special_values(void) {
    /* Divisible#1: an infinite divisor is not decided */
    audit_true("Divisible#1", "Divisible[2, Infinity] === Divisible[2, Infinity]");
    audit_true("Divisible#1", "Divisible[2, -Infinity] =!= True");
    audit_true("Divisible#1", "Divisible[2, ComplexInfinity] =!= True");
    audit_true("Divisible#1", "Divisible[6, 3] === True");
    /* DivisorSigma#1: an infinite index is no number of divisor powers */
    audit_true("DivisorSigma#1", "Quiet[DivisorSigma[Infinity, 6]] =!= 1");
    audit_true("DivisorSigma#1", "Quiet[DivisorSigma[ComplexInfinity, 6]] === Indeterminate");
    audit_true("DivisorSigma#1", "DivisorSigma[2, 6] === 50");
    /* FromContinuedFraction#2: an Infinity term ends the fraction */
    audit_true("FromContinuedFraction#2", "Check[FromContinuedFraction[{1, 2, Infinity}], \"m\"] === 3/2");
    audit_true("FromContinuedFraction#2", "FromContinuedFraction[{1, Infinity, 5}] === 1");
    audit_true("FromContinuedFraction#2", "FromContinuedFraction[{1, 2}] === 3/2");
}

static void test_wave4_numbertheory(void) {
    audit_chk("FromContinuedFraction#1", "FromContinuedFraction[{{-1}}]", "1/2 (-1 - Sqrt[5])");
    audit_chk("FromContinuedFraction#1", "FromContinuedFraction[{-2,{-3}}]", "1/2 (-1 - Sqrt[13])");
    audit_chk("FromContinuedFraction#1", "FromContinuedFraction[{{1}}]", "1/2 (1 + Sqrt[5])");
    audit_chk("GCD#2", "GCD[-1/2]", "1/2");
    audit_chk("GCD#2", "GCD[-2^70/3]", "1180591620717411303424/3");
    audit_chk("PrimitiveRoot#1", "PrimitiveRoot[2000000000005408000000000197359]", "11");
}

int main(void) {
    audit_init();
    TEST(test_primitive_root_list_huge_phi);
    TEST(test_powermod_roots);
    TEST(test_continued_fraction_inexact_last_term);
    TEST(test_wave3_arity_options);
    TEST(test_s11_special_values);
    TEST(test_wave4_numbertheory);
    return audit_finish("audit_number_theory");
}
