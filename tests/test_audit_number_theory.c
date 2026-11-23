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

int main(void) {
    audit_init();
    TEST(test_primitive_root_list_huge_phi);
    return audit_finish("audit_number_theory");
}
