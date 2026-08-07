#include "test_utils.h"
#include "symtab.h"
#include "core.h"

/*
 * PossibleZeroQ assumptions STRESS corpus (auto-derived from the hardening
 * sweep). Each identity holds only on its assumed region (branch cuts of
 * Sqrt/Log/inverse-trig, special functions under integer order, poles,
 * boundary/Equal-pinned ranges, multi-symbol mixed domains) and is paired
 * with the matching non-identity that must stay False. Assumptions use the
 * Element[x,dom] function form (parser has no \[Element] operator).
 *
 * assert_eval_eq uses libc assert() (a no-op under -DNDEBUG); run the binary
 * and grep for FAIL:. Known limitations are listed at the bottom, NOT asserted
 * (asserting current behaviour would lock in a wrong verdict).
 */

static void test_pzq_stress_00(void) {
    assert_eval_eq("PossibleZeroQ[Sqrt[z^2] - z, Assumptions -> Re[z] > 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Sqrt[z^2] + z, Assumptions -> Re[z] < 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Sqrt[z^2] - z, Assumptions -> Re[z] < 0]", "False", 0);
    assert_eval_eq("PossibleZeroQ[Sqrt[x^2] - x, Assumptions -> x >= 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Sqrt[x^2] + x, Assumptions -> x <= 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Sqrt[x^2] - x, Assumptions -> x < 0]", "False", 0);
    assert_eval_eq("PossibleZeroQ[Sqrt[x^2] - x, Assumptions -> 0 <= x]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Sqrt[z^2] - z, Assumptions -> 0 < Re[z]]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Sqrt[a] Sqrt[b] - Sqrt[a b], Assumptions -> a > 0 && b > 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Sqrt[a] Sqrt[b] - Sqrt[a b], Assumptions -> a < 0 && b < 0]", "False", 0);
    assert_eval_eq("PossibleZeroQ[Sqrt[a b] - Sqrt[a] Sqrt[b], Assumptions -> a >= 0 && b >= 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Sqrt[1/z] - 1/Sqrt[z], Assumptions -> Re[z] > 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Conjugate[Sqrt[z]] - Sqrt[Conjugate[z]], Assumptions -> Re[z] > 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Sqrt[x - 3] Sqrt[x - 3] - (x - 3), Assumptions -> x > 3]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Sqrt[1 - Cos[t]^2] - Sin[t], Assumptions -> 0 <= t <= Pi]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Sqrt[1 - Cos[t]^2] - Sin[t], Assumptions -> -Pi <= t <= 0]", "False", 0);
}

static void test_pzq_stress_01(void) {
    assert_eval_eq("PossibleZeroQ[Log[z^2] - 2 Log[z], Assumptions -> Re[z] > 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Log[a b] - Log[a] - Log[b], Assumptions -> a > 0 && b > 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Log[a b] - Log[a] - Log[b], Assumptions -> a < 0 && b < 0]", "False", 0);
    assert_eval_eq("PossibleZeroQ[Log[1/z] + Log[z], Assumptions -> Re[z] > 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Log[Exp[z]] - z, Assumptions -> -Pi < Im[z] <= Pi]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Log[Exp[z]] - z, Assumptions -> Im[z] > Pi]", "False", 0);
    assert_eval_eq("PossibleZeroQ[Exp[Log[z]] - z, Assumptions -> Re[z] > 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Log[x] - Log[Abs[x]], Assumptions -> x > 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Im[Log[x]] - Pi, Assumptions -> x < 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Im[Log[x]], Assumptions -> x > 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Log[z] - Log[Abs[z]], Assumptions -> Element[z, Complexes]]", "False", 0);
    assert_eval_eq("PossibleZeroQ[ArcSin[Sin[x]] - x, Assumptions -> -Pi/2 <= x <= Pi/2]", "True", 0);
    assert_eval_eq("PossibleZeroQ[ArcSin[Sin[x]] - x, Assumptions -> Pi/2 < x < 3 Pi/2]", "False", 0);
    assert_eval_eq("PossibleZeroQ[ArcCos[Cos[x]] - x, Assumptions -> 0 <= x <= Pi]", "True", 0);
    assert_eval_eq("PossibleZeroQ[ArcTan[Tan[x]] - x, Assumptions -> -Pi/2 < x < Pi/2]", "True", 0);
    assert_eval_eq("PossibleZeroQ[ArcTan[z] - I/2 Log[(I + z)/(I - z)], Assumptions -> Im[z] == 0 && -1 < Re[z] < 1]", "True", 0);
}

static void test_pzq_stress_02(void) {
    assert_eval_eq("PossibleZeroQ[ArcTan[1/x] + ArcTan[x] - Pi/2, Assumptions -> x > 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[ArcTan[1/x] + ArcTan[x] - Pi/2, Assumptions -> x < 0]", "False", 0);
    assert_eval_eq("PossibleZeroQ[Abs[x] - x, Assumptions -> x >= 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Abs[x] + x, Assumptions -> x <= 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Abs[Exp[I t]] - 1, Assumptions -> Im[t] == 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Abs[Exp[I t]] - 1, Assumptions -> Element[t, Complexes]]", "False", 0);
    assert_eval_eq("PossibleZeroQ[Abs[z]^2 - z Conjugate[z], Assumptions -> Element[z, Complexes]]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Conjugate[z] - Re[z] + I Im[z], Assumptions -> Element[z, Complexes]]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Re[z] - (z + Conjugate[z])/2, Assumptions -> Element[z, Complexes]]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Sign[x] - x/Abs[x], Assumptions -> x > 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Sign[x] + 1, Assumptions -> x < 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Arg[x] - Pi, Assumptions -> x < 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Floor[x] - x, Assumptions -> Element[x, Integers]]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Ceiling[x] - x, Assumptions -> Element[x, Integers]]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Round[x] - x, Assumptions -> Element[x, Integers]]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Mod[x, 1], Assumptions -> Element[x, Integers]]", "True", 0);
}

static void test_pzq_stress_03(void) {
    assert_eval_eq("PossibleZeroQ[x - Floor[x], Assumptions -> Element[x, Integers]]", "True", 0);
    assert_eval_eq("PossibleZeroQ[x - Floor[x], Assumptions -> x > 0]", "False", 0);
    assert_eval_eq("PossibleZeroQ[Mod[x^2 - x, 2], Assumptions -> Element[x, Integers]]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Gamma[n] - Factorial[n - 1], Assumptions -> Element[n, Integers] && n > 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Gamma[n + 1] - Factorial[n], Assumptions -> Element[n, Integers] && n >= 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Zeta[-2 n], Assumptions -> Element[n, Integers] && n > 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Zeta[-2 n], Assumptions -> Element[n, Integers]]", "False", 0);
    assert_eval_eq("PossibleZeroQ[Zeta[-2 n], Assumptions -> Element[n, PositiveIntegers]]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Sin[n Pi], Assumptions -> Element[n, PositiveIntegers]]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Zeta[1 - 2 n] + BernoulliB[2 n]/(2 n), Assumptions -> Element[n, Integers] && n > 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Fibonacci[n + 2] - Fibonacci[n + 1] - Fibonacci[n], Assumptions -> Element[n, Integers]]", "True", 0);
    assert_eval_eq("PossibleZeroQ[LucasL[n] - Fibonacci[n - 1] - Fibonacci[n + 1], Assumptions -> Element[n, Integers]]", "True", 0);
    assert_eval_eq("PossibleZeroQ[BesselJ[-n, x] - (-1)^n BesselJ[n, x], Assumptions -> Element[n, Integers]]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Binomial[n, 0] - 1, Assumptions -> Element[n, Integers] && n >= 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Binomial[n, n] - 1, Assumptions -> Element[n, Integers] && n >= 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Sin[n Pi/2], Assumptions -> Element[n, Integers]]", "False", 0);
}

static void test_pzq_stress_04(void) {
    assert_eval_eq("PossibleZeroQ[1/Gamma[x], Assumptions -> Element[x, Integers] && x <= 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Sin[x + I y] - (Sin[x] Cosh[y] + I Cos[x] Sinh[y]), Assumptions -> Im[x] == 0 && Im[y] == 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Cos[x + I y] - (Cos[x] Cosh[y] - I Sin[x] Sinh[y]), Assumptions -> Im[x] == 0 && Im[y] == 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Re[Sin[x + I y]] - Sin[x] Cosh[y], Assumptions -> Im[x] == 0 && Im[y] == 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Re[Sqrt[x]] - Sqrt[x], Assumptions -> x >= 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Im[Sqrt[x]], Assumptions -> x >= 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Re[Cos[I y]] - Cosh[y], Assumptions -> Im[y] == 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[(x^a)^b - x^(a b), Assumptions -> Im[x] == 0 && Re[x] > 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Sqrt[x]^2 - x, Assumptions -> Element[x, Complexes]]", "True", 0);
    assert_eval_eq("PossibleZeroQ[(x^2)^(1/2) - Abs[x], Assumptions -> Im[x] == 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[(x^a)^b - x^(a b), Assumptions -> Element[x, Complexes]]", "False", 0);
    assert_eval_eq("PossibleZeroQ[Sqrt[(-1 + I)^2] + (-1 + I)]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Sqrt[(-1 + I)^2] - (-1 + I)]", "False", 0);
    assert_eval_eq("PossibleZeroQ[Zeta[-2 n - 1], Assumptions -> Element[n, Integers] && n > 0]", "False", 0);
    assert_eval_eq("PossibleZeroQ[1/Gamma[x], Assumptions -> Element[x, Integers] && x > 0]", "False", 0);
    assert_eval_eq("PossibleZeroQ[Gamma[x], Assumptions -> Element[x, Integers] && x <= 0]", "False", 0);
}

static void test_pzq_stress_05(void) {
    assert_eval_eq("PossibleZeroQ[BesselJ[n, x] - BesselJ[n + 1, x], Assumptions -> Element[n, Integers]]", "False", 0);
    assert_eval_eq("PossibleZeroQ[Fibonacci[n] - n, Assumptions -> Element[n, Integers]]", "False", 0);
    assert_eval_eq("PossibleZeroQ[Zeta[2 n], Assumptions -> Element[n, Integers] && n > 0]", "False", 0);
    assert_eval_eq("PossibleZeroQ[Sin[n Pi/3], Assumptions -> Element[n, Integers]]", "False", 0);
    assert_eval_eq("PossibleZeroQ[Cos[2 n Pi] - 2, Assumptions -> Element[n, Integers]]", "False", 0);
    assert_eval_eq("PossibleZeroQ[Mod[x, 2], Assumptions -> Element[x, Integers]]", "False", 0);
    assert_eval_eq("PossibleZeroQ[HarmonicNumber[n], Assumptions -> Element[n, Integers] && n > 0]", "False", 0);
    assert_eval_eq("PossibleZeroQ[Log[Exp[z]] - z, Assumptions -> -Pi < Im[z] < Pi]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Log[Exp[z]] - z, Assumptions -> Im[z] == Pi]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Log[Exp[z]] - z, Assumptions -> Im[z] == 3 Pi/2]", "False", 0);
    assert_eval_eq("PossibleZeroQ[Sqrt[x^2] - x, Assumptions -> x == 5]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Sqrt[x^2] - x, Assumptions -> x == -5]", "False", 0);
    assert_eval_eq("PossibleZeroQ[ArcSin[Sin[x]] - x, Assumptions -> x == Pi/3]", "True", 0);
    assert_eval_eq("PossibleZeroQ[BesselJ[-n, x] - (-1)^n BesselJ[n, x], Assumptions -> Element[n, Integers] && x > 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Sqrt[x^2] Sqrt[y^2] - x y, Assumptions -> x >= 0 && y >= 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Sqrt[x^2] Sqrt[y^2] - x y, Assumptions -> x >= 0 && y < 0]", "False", 0);
}

static void test_pzq_stress_06(void) {
    assert_eval_eq("PossibleZeroQ[Log[x^n] - n Log[x], Assumptions -> x > 0 && Element[n, Integers]]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Conjugate[z + w] - (Conjugate[z] + Conjugate[w]), Assumptions -> Element[z, Complexes] && Element[w, Complexes]]", "True", 0);
    assert_eval_eq("PossibleZeroQ[(z^n) Conjugate[z^n] - Abs[z]^(2 n), Assumptions -> Element[z, Complexes] && Element[n, Integers] && n > 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Re[z^2] - (Re[z]^2 - Im[z]^2), Assumptions -> Element[z, Complexes]]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Im[z^2] - 2 Re[z] Im[z], Assumptions -> Element[z, Complexes]]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Re[1/z] - Re[z]/Abs[z]^2, Assumptions -> Element[z, Complexes]]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Abs[x + y] - (x + y), Assumptions -> {x > 0, y > 0}]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Sqrt[x^2] - x, Assumptions -> Element[x, Reals] && x > 0]", "True", 0);
}

/* Known limitations (documented; not asserted):
 *   PossibleZeroQ[Max[x, y] - x, Assumptions -> x >= y]
 *   PossibleZeroQ[Abs[x - y] - (x - y), Assumptions -> x >= y]
 *   PossibleZeroQ[LegendreP[n, x], Assumptions -> Element[n, Integers] && n > 0]
 *   PossibleZeroQ[ChebyshevT[n, x], Assumptions -> Element[n, Integers] && n > 0]
 *   PossibleZeroQ[HermiteH[n, x], Assumptions -> Element[n, Integers] && n > 0]
 *   PossibleZeroQ[Sqrt[z^2] - z, Assumptions -> Re[z^2] > 0]
 *   PossibleZeroQ[x - y, Assumptions -> x == y]
 */

int main(void) {
    symtab_init();
    core_init();

    TEST(test_pzq_stress_00);
    TEST(test_pzq_stress_01);
    TEST(test_pzq_stress_02);
    TEST(test_pzq_stress_03);
    TEST(test_pzq_stress_04);
    TEST(test_pzq_stress_05);
    TEST(test_pzq_stress_06);

    printf("All PossibleZeroQ stress tests passed! (104 cases)\n");
    return 0;
}
