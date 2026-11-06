/* Tests for the zero and characteristic-value heads:
 *
 *   ZetaZero[k]                        nontrivial zeros of Zeta
 *   AiryAiZero / AiryAiPrimeZero /
 *   AiryBiZero / AiryBiPrimeZero [k]   zeros of the Airy functions
 *   MathieuCharacteristicA / B [r, q]  Mathieu characteristic values
 *
 * Reference values, each from a published table or an independent derivation:
 *   - zeta zeros: Odlyzko's tables (gamma_1, gamma_2, gamma_10, gamma_1000 to
 *     30 digits, gamma_{10^6});
 *   - Airy zeros: DLMF Table 9.9.1 (a_1, a_2, a_10, a'_1, b_1, b'_1) and a_500
 *     from the asymptotic series T(3 pi (4k - 1)/8) to t^-10 (DLMF 9.9.6);
 *   - Mathieu: Abramowitz & Stegun Table 20.1 (q = 1, 5, 25), the small-q
 *     series a_nu = nu^2 + q^2/(2 (nu^2 - 1)) + (5 nu^2 + 7) q^4 /
 *     (32 (nu^2 - 1)^3 (nu^2 - 4)) (DLMF 28.6.21 for non-integer nu), the
 *     characteristic value as the eigenvalue of the Fourier matrix computed by
 *     Eigenvalues (a different algorithm), and the defining property: at
 *     a = a_r(q), MathieuC is pi- or 2 pi-periodic in z.
 */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "attr.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- helpers --------------------------------------------------------- */

static double eval_real(const char* input) {
    Expr* e = parse_expression(input);
    ASSERT(e != NULL);
    Expr* r = evaluate(e);
    expr_free(e);
    ASSERT_MSG(r->type == EXPR_REAL, "%s: expected a Real result", input);
    double v = r->data.real;
    expr_free(r);
    return v;
}

static void assert_close(const char* input, double expected, double tol) {
    double v = eval_real(input);
    ASSERT_MSG(fabs(v - expected) <= tol * fmax(1.0, fabs(expected)),
               "%s: expected %.17g, got %.17g", input, expected, v);
}

static void assert_true(const char* input) {
    assert_eval_eq(input, "True", 0);
}

/* ---- ZetaZero ---------------------------------------------------------- */

void test_zetazero_symbolic() {
    assert_eval_eq("ZetaZero[1]", "ZetaZero[1]", 0);
    assert_eval_eq("ZetaZero[k]", "ZetaZero[k]", 0);
    assert_eval_eq("ZetaZero[{1, 2}]", "{ZetaZero[1], ZetaZero[2]}", 0);
    /* a non-integer (or zero) index: ZetaZero::intnz, through the funnel */
    assert_eval_eq("Check[ZetaZero[2.5], fail]", "fail", 0);
    assert_eval_eq("Check[ZetaZero[0], fail]", "fail", 0);
    assert_eval_eq("Quiet[ZetaZero[1.]] === Unevaluated[ZetaZero[1.]]", "True", 0);
    assert_eval_eq("Check[ZetaZero[3], fail]", "ZetaZero[3]", 0);
}

void test_zetazero_machine() {
    assert_close("Im[N[ZetaZero[1]]]", 14.134725141734693790, 1e-15);
    assert_close("Re[N[ZetaZero[1]]]", 0.5, 0);
    assert_close("Im[N[ZetaZero[2]]]", 21.022039638771554993, 1e-15);
    assert_close("Im[N[ZetaZero[10]]]", 49.773832477672302182, 1e-15);
    /* the conjugate below the axis */
    assert_close("Im[N[ZetaZero[-2]]]", -21.022039638771554993, 1e-15);
    /* Odlyzko: zero number 10^6 */
    assert_close("Im[N[ZetaZero[10^6]]]", 600269.677012444955, 1e-14);
    /* it is a zero of Zeta, and of RiemannSiegelZ */
    assert_true("Abs[Zeta[N[ZetaZero[7]]]] < 10^-12");
    assert_true("Abs[RiemannSiegelZ[Im[N[ZetaZero[7]]]]] < 10^-12");
}

void test_zetazero_arbitrary() {
    /* Odlyzko: gamma_1000 = 1419.42248094599568646598903807991 */
    assert_true("Abs[N[Im[ZetaZero[1000]], 30] - "
                "1419.42248094599568646598903807991`35] < 10^-26");
    assert_true("Precision[N[Im[ZetaZero[1000]], 30]] >= 30");
    assert_true("Abs[N[Im[ZetaZero[1]], 50] - "
                "14.134725141734693790457251983562470270784257115699`55] < 10^-48");
    assert_true("Abs[Zeta[N[ZetaZero[3], 40]]] < 10^-36");
}

/* ---- Airy zeros ----------------------------------------------------------- */

void test_airyzero_symbolic() {
    assert_eval_eq("AiryAiZero[1]", "AiryAiZero[1]", 0);
    assert_eval_eq("AiryBiPrimeZero[k]", "AiryBiPrimeZero[k]", 0);
    assert_eval_eq("AiryAiZero[2.] === Unevaluated[AiryAiZero[2.]]", "True", 0);
    assert_eval_eq("N[AiryAiZero[0]]", "AiryAiZero[0]", 0);
    assert_eval_eq("AiryAiZero[{1, 2}]", "{AiryAiZero[1], AiryAiZero[2]}", 0);
}

void test_airyzero_values() {
    /* DLMF Table 9.9.1 */
    assert_close("N[AiryAiZero[1]]", -2.338107410459767039, 1e-15);
    assert_close("N[AiryAiZero[2]]", -4.087949444130970617, 1e-15);
    assert_close("N[AiryAiZero[10]]", -12.828776752865757200, 1e-15);
    assert_close("N[AiryAiPrimeZero[1]]", -1.018792971647471089, 1e-15);
    assert_close("N[AiryBiZero[1]]", -1.173713222709127925, 1e-15);
    assert_close("N[AiryBiPrimeZero[1]]", -2.294439682614123247, 1e-15);
    /* each is a zero of its function */
    assert_true("Abs[AiryAi[N[AiryAiZero[37]]]] < 10^-12");
    assert_true("Abs[AiryBiPrime[N[AiryBiPrimeZero[5]]]] < 10^-12");
    /* arbitrary precision: a_500 from the asymptotic series (26 digits agree) */
    assert_true("Abs[N[AiryAiZero[500], 30] - "
                "(-177.00925104647920682546615344202`35)] < 10^-25");
    assert_true("Precision[N[AiryAiZero[500], 30]] >= 30");
    assert_true("Abs[AiryAi[N[AiryAiZero[3], 40]]] < 10^-37");
}

/* ---- Mathieu characteristic values --------------------------------------- */

void test_mathieu_char_exact() {
    assert_eval_eq("MathieuCharacteristicA[0, 25]", "MathieuCharacteristicA[0, 25]", 0);
    assert_eval_eq("MathieuCharacteristicA[r, q]", "MathieuCharacteristicA[r, q]", 0);
    /* q = 0: a_r = b_r = r^2 */
    assert_eval_eq("MathieuCharacteristicA[3, 0]", "9", 0);
    assert_eval_eq("MathieuCharacteristicB[5/2, 0]", "25/4", 0);
    /* no b_0 */
    assert_eval_eq("Check[MathieuCharacteristicB[0, 2.], fail]", "fail", 0);
    assert_eval_eq("Head[Quiet[MathieuCharacteristicB[0, 2.]]]", "MathieuCharacteristicB", 0);
}

void test_mathieu_char_machine() {
    /* Abramowitz & Stegun Table 20.1 */
    assert_close("MathieuCharacteristicA[0, 1.]", -0.45513860410741354823, 1e-13);
    assert_close("MathieuCharacteristicA[0, 5.]", -5.80004602085150841910, 1e-13);
    assert_close("MathieuCharacteristicA[1, 5.]", 1.85818754154775088515, 1e-13);
    assert_close("MathieuCharacteristicA[2, 5.]", 7.44910973952917781248, 1e-13);
    assert_close("MathieuCharacteristicA[3, 1.]", 9.07836884720310199254, 1e-13);
    assert_close("MathieuCharacteristicB[1, 5.]", -5.79008060, 1e-8);
    assert_close("MathieuCharacteristicB[2, 5.]", 2.09946045, 1e-8);
    assert_close("MathieuCharacteristicB[2, 1.]", 3.91702477, 1e-8);
    assert_close("MathieuCharacteristicA[0, 25.]", -40.256779546566787276, 1e-13);
    /* symmetries: a_{-r} = a_r, b_{-r} = b_r, a_{2m+1}(-q) = b_{2m+1}(q),
     * a_{2m}(-q) = a_{2m}(q) */
    assert_true("MathieuCharacteristicA[-3, 2.] == MathieuCharacteristicA[3, 2.]");
    assert_true("MathieuCharacteristicB[-4, 2.] == MathieuCharacteristicB[4, 2.]");
    assert_true("Abs[MathieuCharacteristicA[3, -2.] - MathieuCharacteristicB[3, 2.]] < 10^-13");
    assert_true("Abs[MathieuCharacteristicA[4, -2.] - MathieuCharacteristicA[4, 2.]] < 10^-13");
    /* exact r with a machine q, and a machine r */
    assert_true("MathieuCharacteristicA[2., 1.] == MathieuCharacteristicA[2, 1.]");
    /* non-integer order: A = B, even in r, and the small-q series */
    assert_true("MathieuCharacteristicA[1.5, 2.] == MathieuCharacteristicB[1.5, 2.]");
    assert_true("MathieuCharacteristicA[-1.5, 2.] == MathieuCharacteristicA[1.5, 2.]");
    assert_close("MathieuCharacteristicA[0.5, 0.01]",
                 0.25 - 0.0001 / 1.5 + 8.25e-8 / 50.625, 1e-12);
    assert_close("MathieuCharacteristicA[2.3, 0.02]",
                 2.3 * 2.3 + 0.0004 / (2 * (2.3 * 2.3 - 1)), 1e-9);
    /* continuity across an integer order from both sides */
    assert_true("Abs[MathieuCharacteristicA[2 + 10^-9, 3.] - MathieuCharacteristicA[2, 3.]] < 10^-6");
    assert_true("Abs[MathieuCharacteristicA[3 - 10^-9, 3.] - MathieuCharacteristicB[3, 3.]] < 10^-6");
    /* the eigenvalue of the Fourier matrix, by Eigenvalues */
    assert_true("Abs[MathieuCharacteristicA[0, 25.] - Min[Eigenvalues[N[Table[Which[i == j, "
                "(2 (i - 1))^2, Abs[i - j] == 1, If[Min[i, j] == 1, Sqrt[2] 25, 25], True, 0], "
                "{i, 40}, {j, 40}]]]]] < 10^-11");
    /* the defining property: MathieuC at a_r(q) has period pi (even r) or
     * 2 pi (odd r), MathieuS at b_r(q) likewise */
    assert_true("With[{a = MathieuCharacteristicA[2, 1.5]}, "
                "Abs[MathieuC[a, 1.5, 0.3 + Pi] - MathieuC[a, 1.5, 0.3]] < 10^-10]");
    assert_true("With[{a = MathieuCharacteristicA[3, 1.5]}, "
                "Abs[MathieuC[a, 1.5, 0.3 + 2 Pi] - MathieuC[a, 1.5, 0.3]] < 10^-10]");
    assert_true("With[{b = MathieuCharacteristicB[4, 2.5]}, "
                "Abs[MathieuS[b, 2.5, 0.4 + Pi] - MathieuS[b, 2.5, 0.4]] < 10^-10]");
    /* large q: the asymptotic expansion (DLMF 28.8.1, s = 2r + 1 = 1)
     * -2q + 2 s Sqrt[q] - (s^2 + 1)/8 - (s^3 + 3 s)/(2^7 Sqrt[q])
     * - (5 s^4 + 34 s^2 + 9)/(2^12 q), next term ~ 10^-10 at q = 10^4 */
    assert_true("Abs[MathieuCharacteristicA[0, 10000.] - (-20000 + 200 - 1/4 - 4/12800 "
                "- 48/40960000)] < 10^-8");
    /* large r: the series nu^2 + q^2/(2 (nu^2 - 1)) + (5 nu^2 + 7) q^4 /
     * (32 (nu^2 - 1)^3 (nu^2 - 4)), next term ~ 10^-11 at r = 50, q = 10 */
    assert_true("Abs[MathieuCharacteristicA[50, 10.] - (2500 + 100/4998 + "
                "12507 10^4/(32 2499^3 2496))] < 10^-10");
    /* complex q stays unevaluated */
    assert_eval_eq("Head[MathieuCharacteristicA[2, 1. + I]]", "MathieuCharacteristicA", 0);
}

void test_mathieu_char_arbitrary() {
    /* a_0(25): two methods in the benchmark notes agree to 30 digits */
    assert_true("Abs[N[MathieuCharacteristicA[0, 25], 30] - "
                "(-40.2567795465667872758017193249`35)] < 10^-27");
    assert_true("Precision[N[MathieuCharacteristicA[0, 25], 30]] >= 30");
    /* against Eigenvalues of the 60 x 60 Fourier matrix at 60 digits */
    assert_true("Abs[N[MathieuCharacteristicA[2, 25], 50] - Sort[Eigenvalues[N[Table[Which[i == j, "
                "(2 (i - 1))^2, Abs[i - j] == 1, If[Min[i, j] == 1, Sqrt[2] 25, 25], True, 0], "
                "{i, 60}, {j, 60}], 60]]][[2]]] < 10^-47");
    /* non-integer order at 40 digits against the machine value */
    assert_true("Abs[N[MathieuCharacteristicA[7/3, 5], 40] - MathieuCharacteristicA[7/3, 5.]] < 10^-13");
    /* b_7(3): the 4th eigenvalue of the sin((2k+1) z) matrix */
    assert_true("Abs[N[MathieuCharacteristicB[7, 3], 30] - Sort[Eigenvalues[N[Table[Which[i == j, "
                "(2 i - 1)^2 - If[i == 1, 3, 0], Abs[i - j] == 1, 3, True, 0], "
                "{i, 50}, {j, 50}], 50]]][[4]]] < 10^-27");
    assert_true("Precision[MathieuCharacteristicA[1, 2.`40]] >= 39");
}

void test_mathieu_char_surfaces() {
    /* scalar == unpacked List == packed List == NDArray == Compile, exactly */
    assert_true("With[{v = Range[0., 20., 0.5]}, "
                "MathieuCharacteristicA[3, v] === Table[MathieuCharacteristicA[3, x], {x, v}]]");
    assert_true("PackedArrayQ[MathieuCharacteristicA[3, Range[0., 20., 0.5]]]");
    assert_true("With[{v = Range[0., 20., 0.5]}, "
                "MathieuCharacteristicB[2, FromPackedArray[v]] === MathieuCharacteristicB[2, v]]");
    assert_true("With[{v = Range[0., 3., 0.25]}, Normal[MathieuCharacteristicA[1, NDArray[v]]] "
                "=== MathieuCharacteristicA[1, FromPackedArray[v]]]");
    assert_eval_eq("Head[MathieuCharacteristicA[1, NDArray[{1.5, 2.5}]]]", "NDArray", 0);
    assert_eval_eq("CompileDiagnostics[{{x, _Real}}, MathieuCharacteristicA[2, x]][[1]]",
                   "\"Compiled\" -> True", 0);
    assert_eval_eq("CompileDiagnostics[{{v, _Real, 1}}, MathieuCharacteristicB[2, v]][[1]]",
                   "\"Compiled\" -> True", 0);
    assert_true("Compile[{{x, _Real}}, MathieuCharacteristicA[2, x]][3.25] === "
                "MathieuCharacteristicA[2, 3.25]");
    assert_true("With[{w = Range[0., 5., 0.5]}, Compile[{{v, _Real, 1}}, "
                "MathieuCharacteristicB[2, v]][w] === MathieuCharacteristicB[2, w]]");
    /* an int64 buffer of exact q is exact: unevaluated element-wise */
    assert_eval_eq("MathieuCharacteristicA[2, NDArray[{1, 2}, DataType -> \"int64\"]]",
                   "{MathieuCharacteristicA[2, 1], MathieuCharacteristicA[2, 2]}", 0);
}

/* ---- BesselYZero ------------------------------------------------------------ */

void test_besselyzero() {
    /* exact stays symbolic (as BesselJZero) */
    assert_eval_eq("BesselYZero[0, 1]", "BesselYZero[0, 1]", 0);
    /* DLMF Table 10.21.1 (and mpmath besselyzero to 30 digits) */
    assert_close("N[BesselYZero[0, 1]]", 0.89357696627916752158, 1e-15);
    assert_close("N[BesselYZero[0, 2]]", 3.95767841931485786838, 1e-15);
    assert_close("N[BesselYZero[1, 1]]", 2.19714132603101703515, 1e-15);
    assert_close("N[BesselYZero[1, 2]]", 5.42968104079413513277, 1e-15);
    /* each is a zero of Y_n: small order, large order, large k (McMahon) */
    assert_true("Abs[BesselY[0, N[BesselYZero[0, 7]]]] < 10^-14");
    assert_true("Abs[BesselY[25/2, N[BesselYZero[25/2, 3]]]] < 10^-14");
    assert_true("Abs[BesselY[2, N[BesselYZero[2, 1000]]]] < 10^-14");
    /* consecutive zeros interlace with those of J (DLMF 10.21.3):
     * j_{n,k} < y_{n,k+1} < j_{n,k+1} */
    assert_true("N[BesselJZero[3, 4]] < N[BesselYZero[3, 5]] < N[BesselJZero[3, 5]]");
    /* arbitrary precision */
    assert_true("Abs[BesselY[0, N[BesselYZero[0, 1], 40]]] < 10^-38");
    assert_true("Abs[N[BesselYZero[0, 1], 30] - 0.893576966279167521584887102058`30] < 10^-28");
    /* surfaces: scalar == packed == Compile */
    assert_true("With[{v = N[Range[1, 30]]}, BesselYZero[1, v] === Table[BesselYZero[1, k], {k, v}]]");
    assert_eval_eq("CompileDiagnostics[{{k, _Real}}, BesselYZero[1, k]][[1]]",
                   "\"Compiled\" -> True", 0);
    assert_true("Compile[{{k, _Real}}, BesselYZero[2, k]][9.] === BesselYZero[2, 9.]");
}

void test_attributes_docs() {
    const char* heads[] = { "ZetaZero", "AiryAiZero", "AiryAiPrimeZero", "AiryBiZero",
                            "AiryBiPrimeZero", "MathieuCharacteristicA", "MathieuCharacteristicB" };
    for (size_t i = 0; i < sizeof heads / sizeof heads[0]; i++) {
        SymbolDef* d = symtab_get_def(heads[i]);
        ASSERT_MSG(d && d->docstring && strstr(d->docstring, heads[i]), "%s docstring", heads[i]);
        ASSERT_MSG((d->attributes & ATTR_LISTABLE) && (d->attributes & ATTR_NHOLDFIRST) &&
                   (d->attributes & ATTR_NUMERICFUNCTION) && (d->attributes & ATTR_PROTECTED),
                   "%s attributes", heads[i]);
    }
}

int main() {
    symtab_init();
    core_init();

    TEST(test_zetazero_symbolic);
    TEST(test_zetazero_machine);
    TEST(test_zetazero_arbitrary);
    TEST(test_airyzero_symbolic);
    TEST(test_airyzero_values);
    TEST(test_mathieu_char_exact);
    TEST(test_mathieu_char_machine);
    TEST(test_mathieu_char_arbitrary);
    TEST(test_mathieu_char_surfaces);
    TEST(test_besselyzero);
    TEST(test_attributes_docs);

    printf("All special-zero and Mathieu characteristic-value tests passed.\n");
    return 0;
}
