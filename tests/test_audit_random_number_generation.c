/* Builtin-audit regressions: random number generation (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

static void test_wave4_random(void) {
    audit_chk("RandomChoice#1", "SeedRandom[1]; RandomChoice[{10^400, 1} -> {a, b}, 10]", "{a, a, a, a, a, a, a, a, a, a}");
    audit_chk("RandomChoice#1", "SeedRandom[1]; Union[RandomChoice[{1.5*^308, 1.5*^308} -> {a, b}, 40]]", "{a, b}");
    audit_chk("RandomSample#1", "RandomSample[{}]", "{}");
    audit_chk("RandomSample#1", "RandomSample[{}, 0]", "{}");
    audit_chk("RandomSample#1", "RandomSample[{}, UpTo[3]]", "{}");
    audit_chk("SeedRandom#1", "SeedRandom[2^70]; e1w4 = RandomReal[]; SeedRandom[-2^70]; e1w4 == RandomReal[]", "False");
}

/* Wave 5 G7 (arith-misc): low-severity findings. Message counts are read
 * off $MessageList inside a Block (each diagnostic once per evaluation). */
static void test_wave5_g7_random_number_generation(void) {
    audit_chk("RandomChoice#2", "RandomChoice[{10^-400, 0} -> {a, b}, 5]", "{a, a, a, a, a}");
    audit_chk("RandomChoice#3", "Quiet[Check[RandomChoice[], fired, RandomChoice::argt]]", "fired");
    audit_chk("RandomChoice#3", "Quiet[Check[RandomChoice[{a, b}, 3, 4], fired, RandomChoice::argt]]", "fired");
    audit_chk("RandomChoice#3", "Quiet[Check[RandomChoice[{}], fired, RandomChoice::lrwl]]", "fired");
    audit_chk("RandomChoice#3", "Quiet[Check[RandomChoice[x], fired, RandomChoice::lrwl]]", "fired");
    audit_chk("RandomComplex#1", "Precision[Re[RandomComplex[1 + I, WorkingPrecision -> 10]]]", "10.0");
    audit_chk("RandomComplex#2", "Quiet[Check[RandomComplex[1 + I, -1], fired, RandomComplex::array]]", "fired");
    audit_chk("RandomComplex#2", "Quiet[Check[RandomComplex[1 + I, 2.5], fired, RandomComplex::array]]", "fired");
    audit_chk("RandomInteger#1", "Quiet[Check[RandomInteger[1.5], fired, RandomInteger::unifr]]", "fired");
    audit_chk("RandomInteger#1", "Quiet[Check[RandomInteger[10, -1], fired, RandomInteger::array]]", "fired");
    audit_chk("RandomInteger#1", "Quiet[Check[RandomInteger[1, 2, 3], fired, RandomInteger::argb]]", "fired");
    audit_chk("RandomInteger#1", "Length[RandomInteger[10, Developer`ToPackedArray[{2}]]]", "2");
    audit_chk("RandomPrime#1", "ArrayDepth[RandomPrime[10, Table[1, {64}]]]", "64");
    audit_chk("RandomPrime#1", "ArrayDepth[RandomPrime[10, Table[1, {63}]]]", "63");
    audit_chk("RandomReal#1", "Quiet[Check[RandomReal[x], fired, RandomReal::unifr]]", "fired");
    audit_chk("RandomReal#1", "Quiet[Check[RandomReal[{1, 2}, -1], fired, RandomReal::array]]", "fired");
    audit_chk("RandomReal#1", "Quiet[Check[RandomReal[{1, 2}, 2.5], fired, RandomReal::array]]", "fired");
    audit_chk("RandomReal#1", "Quiet[Check[RandomReal[1, WorkingPrecision -> -1], fired, RandomReal::precw]]", "fired");
    audit_chk("RandomSample#2", "Sort[RandomSample[{a, b, c}, UpTo[10^30]]]", "{a, b, c}");
    audit_chk("RandomSample#3", "Quiet[Check[RandomSample[{a, b, c}, 5], fired, RandomSample::smplen]]", "fired");
    audit_chk("RandomSample#3", "Quiet[RandomSample[{a, b, c}, 5]]", "RandomSample[{a, b, c}, 5]");
    audit_chk("SeedRandom#2", "SeedRandom[2^70 + 1]; w5a = RandomReal[]; SeedRandom[2^70 + 2^65 + 1]; w5a == RandomReal[]", "False");
    audit_chk("SeedRandom#3", "SeedRandom[-5]; w5c = RandomInteger[10^30]; SeedRandom[2^64 - 5]; w5c == RandomInteger[10^30]", "False");
    audit_chk("SeedRandom#3", "SeedRandom[-5]; w5c = RandomReal[]; SeedRandom[2^64 - 5]; w5c == RandomReal[]", "False");
    audit_chk("SeedRandom#4", "SeedRandom[\"abc\"]; w5s = RandomReal[]; SeedRandom[\"abc\"]; w5s == RandomReal[]", "True");
    audit_chk("SeedRandom#4", "SeedRandom[\"abc\"]; w5s = RandomReal[]; SeedRandom[\"abd\"]; w5s == RandomReal[]", "False");
    audit_chk("SeedRandom#4", "SeedRandom[\"abc\"]", "Null");
}

int main(void) {
    audit_init();
    TEST(test_wave4_random);
    TEST(test_wave5_g7_random_number_generation);
    return audit_finish("audit_random_number_generation");
}
