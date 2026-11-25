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

int main(void) {
    audit_init();
    TEST(test_wave4_random);
    return audit_finish("audit_random_number_generation");
}
