/* Builtin-audit regressions: bitwise (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

static void test_wave4_integerpart_complex(void) {
    audit_chk("IntegerPart#1", "IntegerPart[Sqrt[2] + I Pi]", "1 + 3*I");
    audit_chk("IntegerPart#1", "IntegerPart[Sqrt[-2]]", "I");
    audit_chk("IntegerPart#1", "IntegerPart[1/3 + I Sqrt[2]]", "I");
}

/* Wave 5 G7 (arith-misc): low-severity findings. Message counts are read
 * off $MessageList inside a Block (each diagnostic once per evaluation). */
static void test_wave5_g7_bitwise(void) {
    audit_chk("BitShiftLeft#1", "Block[{$MessageList = {}}, BitShiftLeft[{1, 2, 3}, 2.]; Length[$MessageList]]", "3");
    audit_chk("BitXor#1", "Block[{$MessageList = {}}, BitXor[x, 2.5]; Length[$MessageList]]", "1");
    audit_chk("BitXor#1", "Block[{$MessageList = {}}, BitXor[2.5, 3]; Length[$MessageList]]", "1");
}

int main(void) {
    audit_init();
    TEST(test_wave4_integerpart_complex);
    TEST(test_wave5_g7_bitwise);
    return audit_finish("audit_bitwise");
}
