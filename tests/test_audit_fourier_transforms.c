/* Builtin-audit regressions: fourier transforms (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

static void test_wave4_exact_constants(void) {
    audit_chk("FourierDCT#1", "FourierDCT[{1, Pi}] == FourierDCT[N[{1, Pi}]]", "True");
    audit_chk("FourierDCT#1", "Head /@ FourierDCT[{1, Sqrt[2]}]", "{Real, Real}");
    audit_chk("FourierDST#1", "FourierDST[{Pi, E}] == FourierDST[N[{Pi, E}]]", "True");
    audit_chk("FourierDST#1", "FourierDST[{1.5, Pi}] == FourierDST[{1.5, N[Pi]}]", "True");
    audit_chk("InverseFourier#1", "InverseFourier[{Pi,E}] == InverseFourier[N[{Pi,E}]]", "True");
    audit_chk("InverseFourier#1", "Head /@ InverseFourier[{Sqrt[2],1}]", "{Real, Real}");
    audit_chk("InverseFourier#1", "Fourier[{x, 1}]", "{(1 + x)/Sqrt[2], (-1 + x)/Sqrt[2]}");
}

int main(void) {
    audit_init();
    TEST(test_wave4_exact_constants);
    return audit_finish("audit_fourier_transforms");
}
