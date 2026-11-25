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

/* Wave 5 G7 (arith-misc): low-severity findings. Message counts are read
 * off $MessageList inside a Block (each diagnostic once per evaluation). */
static void test_wave5_g7_fourier_transforms(void) {
    audit_chk("InverseFourier#2", "InverseFourier[{1, 2, 3, 4}, FourierParameters -> {0, 1}, FourierParameters -> {-1, 1}]", "{5.0, -1.0 + 1.0*I, -1.0, -1.0 - 1.0*I}");
    audit_chk("InverseFourier#3", "Quiet[Check[InverseFourier[{1, 2, 3, 4}, FourierParameters -> {0, 0}], fired, InverseFourier::fpopt]]", "fired");
    audit_chk("InverseFourier#4", "Quiet[Check[InverseFourier[{}], fired, InverseFourier::fftl]]", "fired");
    audit_chk("InverseFourier#4", "Quiet[Check[InverseFourier[{{1, 2}, {3}}], fired, InverseFourier::fftl]]", "fired");
    audit_chk("InverseFourier#4", "Quiet[Check[InverseFourier[5], fired, InverseFourier::fftl]]", "fired");
    audit_chk("ListConvolve#2", "Quiet[Check[ListConvolve[{1, 2}, {1, 2, 3}, 3], fired, ListConvolve::klist]]", "fired");
    audit_chk("ListConvolve#2", "ListConvolve[{1, 2}, {1, 2, 3}, 2]", "{4, 7, 7}");
}

int main(void) {
    audit_init();
    TEST(test_wave4_exact_constants);
    TEST(test_wave5_g7_fourier_transforms);
    return audit_finish("audit_fourier_transforms");
}
