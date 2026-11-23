/* Builtin-audit regressions: arithmetic (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* ComplexExpand#1 (critical): the TargetFunctions -> Conjugate path freed its
 * real and imaginary parts twice after moving them into the result; the
 * recycled nodes then turned up inside the caller's expression (and the
 * kernel crashed).  The results are also checked for the mathematics. */
static void test_complex_expand_conjugate(void) {
    audit_chk("ComplexExpand#1", "ComplexExpand[Exp[x], TargetFunctions -> Conjugate]", "E^x");
    audit_chk("ComplexExpand#1", "Head[ComplexExpand[Exp[x], TargetFunctions -> Conjugate]]", "Power");
    audit_chk("ComplexExpand#1", "f[ComplexExpand[2^x, TargetFunctions -> Conjugate]]", "f[2^x]");
    audit_chk("ComplexExpand#1", "ComplexExpand[x + 1, TargetFunctions -> Conjugate]", "1 + x");
    audit_chk("ComplexExpand#1", "ComplexExpand[Sin[x], TargetFunctions -> Conjugate]", "Sin[x]");
    /* z complex: re + I im must reassemble z itself. */
    audit_true("ComplexExpand#1",
               "PossibleZeroQ[ComplexExpand[z^2, {z}, TargetFunctions -> Conjugate] - z^2]");
    audit_true("ComplexExpand#1",
               "PossibleZeroQ[ComplexExpand[Sqrt[z], {z}, TargetFunctions -> Conjugate] - Sqrt[z]]");
    /* The Re / Im / Abs wrappers keep their (previously correct) answers. */
    audit_true("ComplexExpand#1",
               "PossibleZeroQ[ComplexExpand[Re[z^2], {z}, TargetFunctions -> Conjugate] - (z^2 + Conjugate[z]^2)/2]");
}

int main(void) {
    audit_init();
    TEST(test_complex_expand_conjugate);
    return audit_finish("audit_arithmetic");
}
