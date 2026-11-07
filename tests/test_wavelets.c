/* Tests for HaarWavelet / DaubechiesWavelet, WaveletPhi, WaveletPsi and
 * WaveletFilterCoefficients.  Expected values from the reference
 * (wolframscript): the Daubechies values are the cascade's dyadic values at
 * spacing 2^-8, linearly interpolated -- WaveletPhi[DaubechiesWavelet[2], 1.3]
 * = 0.3406737739116784, and 0.34067377391167836263936642769158893255 at 40
 * digits (also an independent mpmath cascade). */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void check(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    char* s = expr_to_string(r);
    if (strcmp(s, expected) != 0) {
        printf("Wavelet test failed: %s\n  expected: %s\n  got:      %s\n", input, expected, s);
        ASSERT(0);
    }
    free(s);
    expr_free(e);
    expr_free(r);
}

/* |value - ref| < tol, printed True */
#define CLOSE(expr, ref, tol) check("Max[Abs[(" expr ") - (" ref ")]] < " tol, "True")

static void test_haar(void) {
    check("WaveletPhi[HaarWavelet[], 0.3]", "1.0");
    check("WaveletPhi[HaarWavelet[], 1/2]", "1.0");
    check("WaveletPhi[HaarWavelet[], 1]", "0.0");
    check("WaveletPhi[HaarWavelet[], 0.]", "1.0");
    check("WaveletPsi[HaarWavelet[], 0.3]", "1.0");
    check("WaveletPsi[HaarWavelet[], 0.7]", "-1.0");
    check("WaveletPsi[HaarWavelet[], 1/2]", "0.0");
    check("WaveletPsi[HaarWavelet[], 0.]", "0.0");
    check("WaveletPhi[HaarWavelet[], x] /. x -> 0.5", "1");
    check("Head[WaveletPsi[HaarWavelet[], x]]", "Piecewise");
}

static void test_daubechies_values(void) {
    CLOSE("WaveletPhi[DaubechiesWavelet[2], 1.]", "1.3660254037844386", "1*^-15");
    CLOSE("WaveletPhi[DaubechiesWavelet[2], 0.5]", "0.9330127018922193", "1*^-15");
    CLOSE("WaveletPhi[DaubechiesWavelet[2], 1.3]", "0.3406737739116784", "1*^-15");
    CLOSE("WaveletPhi[DaubechiesWavelet[2], 1/3]", "0.704616143501375", "1*^-15");
    CLOSE("WaveletPsi[DaubechiesWavelet[2], 1.3]", "-0.012500835765420517", "1*^-15");
    CLOSE("WaveletPsi[DaubechiesWavelet[2], 0.5]", "1.732050807568877", "1*^-14");
    check("WaveletPhi[DaubechiesWavelet[2], 4.]", "0.0");
    check("WaveletPhi[DaubechiesWavelet[2], -0.1]", "0.0");
    CLOSE("WaveletPhi[DaubechiesWavelet[3], 1.7]", "0.20033345258187266", "1*^-14");
    CLOSE("WaveletPsi[DaubechiesWavelet[3], 0.3]", "0.6058574872971068", "1*^-14");
    CLOSE("WaveletPsi[DaubechiesWavelet[3], -1.5]", "0.064081996982659", "1*^-14");
    CLOSE("WaveletPhi[DaubechiesWavelet[4], 2.2]", "-0.29182565197467675", "1*^-14");
    CLOSE("WaveletPhi[DaubechiesWavelet[1], 0.999]", "0.256", "1*^-13");
    check("WaveletPsi[DaubechiesWavelet[1], 0.5]", "-1.0");
    CLOSE("WaveletPhi[DaubechiesWavelet[], 0.25]", "0.6372595264191645", "1*^-15");
    CLOSE("WaveletPhi[DaubechiesWavelet[2], 1.3, MaxRecursion -> 4]", "0.3323557158514988", "1*^-15");
    CLOSE("WaveletPhi[DaubechiesWavelet[2], 1.3, MaxRecursion -> 12]", "0.34080374356886867", "1*^-15");
}

static void test_precision(void) {
    check("Abs[WaveletPhi[DaubechiesWavelet[2], 1.3`40, WorkingPrecision -> 40] - "
          "0.34067377391167836263936642769158893255`40] < 10^-37", "True");
    check("Precision[WaveletPhi[DaubechiesWavelet[2], 13/10, WorkingPrecision -> 40]] >= 39", "True");
    check("MachineNumberQ[WaveletPhi[DaubechiesWavelet[2], 1.3, WorkingPrecision -> 30]]", "True");
    check("MachineNumberQ[WaveletPhi[DaubechiesWavelet[2], 1/3]]", "True");
}

static void test_symbolic(void) {
    check("w = WaveletPhi[DaubechiesWavelet[2], x]; {Head[w], w[[1, 1, 2]]}", "{Piecewise, 0 <= x <= 3}");
    check("Abs[(WaveletPhi[DaubechiesWavelet[2], x] /. x -> 1.3) - 0.3406737739116784] < 1*^-15", "True");
    check("WaveletPsi[DaubechiesWavelet[2], x][[1, 1, 2]]", "-1 <= x <= 2");
}

static void test_filters(void) {
    CLOSE("WaveletFilterCoefficients[DaubechiesWavelet[2], \"PrimalLowpass\"][[All, 2]]",
          "{(1 + Sqrt[3])/8, (3 + Sqrt[3])/8, (3 - Sqrt[3])/8, (1 - Sqrt[3])/8}", "1*^-15");
    check("WaveletFilterCoefficients[DaubechiesWavelet[2], \"PrimalHighpass\"][[All, 1]]", "{-2, -1, 0, 1}");
    check("Max[Abs[WaveletFilterCoefficients[DaubechiesWavelet[3], \"PrimalLowpass\", WorkingPrecision -> 25][[All, 2]] - "
          "{0.23523360389208184037666455378558522915`25, 0.5705584579157218128800628148297107791`25, "
          "0.32518250026311626425346741451708064159`25, -0.0954672077841636807533291075711704583`25, "
          "-0.06041610415519810463013196830266587074`25, 0.02490874986844186787326629274145967921`25}]] < "
          "10^-23", "True");
    check("Total[WaveletFilterCoefficients[DaubechiesWavelet[6]][[All, 2]]]", "1.0");
    check("DaubechiesWavelet[]", "DaubechiesWavelet[2]");
}

static void test_attributes(void) {
    check("Attributes[WaveletPhi]", "{Protected}");
    check("Attributes[WaveletPsi]", "{Protected}");
    check("WaveletPhi[foo, 1.]", "WaveletPhi[foo, 1.0]");
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_haar);
    TEST(test_daubechies_values);
    TEST(test_precision);
    TEST(test_symbolic);
    TEST(test_filters);
    TEST(test_attributes);
    printf("All wavelet tests passed!\n");
    return 0;
}
