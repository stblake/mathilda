/* Builtin-audit regressions: bitwise (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

static void test_wave4_integerpart_complex(void) {
    audit_chk("IntegerPart#1", "IntegerPart[Sqrt[2] + I Pi]", "1 + 3*I");
    audit_chk("IntegerPart#1", "IntegerPart[Sqrt[-2]]", "I");
    audit_chk("IntegerPart#1", "IntegerPart[1/3 + I Sqrt[2]]", "I");
}

int main(void) {
    audit_init();
    TEST(test_wave4_integerpart_complex);
    return audit_finish("audit_bitwise");
}
