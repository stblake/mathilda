#include "test_utils.h"
#include "symtab.h"
#include "core.h"

/*
 * PossibleZeroQ with assumptions.
 *
 * PossibleZeroQ honours an `Assumptions -> ...` option and an ambient
 * Assuming[] / $Assumptions scope by restricting its Schwartz-Zippel sampler
 * to the assumed region (integer / real / complex domain, sign, finite range,
 * Re/Im-part constraints). An identity that holds only on that region is then
 * recognised as zero; a genuine non-zero on that region is still False. No
 * Simplify is ever invoked -- PossibleZeroQ stays a self-contained
 * numeric/structural test.
 *
 * The assumptions below are written with the `Element[x, dom]` function form
 * (the surface `\[Element]` operator is not part of Mathilda's parser); it is
 * semantically identical.
 *
 * NOTE: assert_eval_eq uses libc assert(), which is a no-op under -DNDEBUG.
 * Run this binary and grep the output for `FAIL:`; a green run prints none.
 */

/* ---- Integer-domain identities ---------------------------------------- */

void test_pzq_integer_domain(void) {
    assert_eval_eq("PossibleZeroQ[Sin[n Pi], Assumptions -> Element[n, Integers]]",
                   "True", 0);
    assert_eval_eq("PossibleZeroQ[(-1)^n - Cos[n Pi], Assumptions -> Element[n, Integers]]",
                   "True", 0);
    assert_eval_eq("PossibleZeroQ[BesselJ[-n, x] - (-1)^n BesselJ[n, x], "
                   "Assumptions -> Element[n, Integers]]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Exp[2 Pi I k] - 1, Assumptions -> Element[k, Integers]]",
                   "True", 0);
    assert_eval_eq("PossibleZeroQ[x - Floor[x], Assumptions -> Element[x, Integers]]",
                   "True", 0);
}

void test_pzq_integer_domain_with_sign(void) {
    assert_eval_eq("PossibleZeroQ[Zeta[-2 n], Assumptions -> Element[n, Integers] && n > 0]",
                   "True", 0);
    assert_eval_eq("PossibleZeroQ[1/Gamma[x], Assumptions -> Element[x, Integers] && x <= 0]",
                   "True", 0);
}

/* ---- Real sign / range identities ------------------------------------- */

void test_pzq_real_sign(void) {
    assert_eval_eq("PossibleZeroQ[Sqrt[x^2] - x, Assumptions -> x >= 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Log[a b] - (Log[a] + Log[b]), Assumptions -> a > 0 && b > 0]",
                   "True", 0);
    assert_eval_eq("PossibleZeroQ[Abs[x + y] - (x + y), Assumptions -> x > 0 && y > 0]",
                   "True", 0);
    assert_eval_eq("PossibleZeroQ[Sign[x] + 1, Assumptions -> x < 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Arg[x], Assumptions -> x > 0]", "True", 0);
}

void test_pzq_real_range(void) {
    assert_eval_eq("PossibleZeroQ[Sqrt[1 - Sin[theta]^2] - Cos[theta], "
                   "Assumptions -> -Pi/2 <= theta <= Pi/2]", "True", 0);
    assert_eval_eq("PossibleZeroQ[ArcSin[Sin[x]] - x, Assumptions -> -Pi/2 <= x <= Pi/2]",
                   "True", 0);
    assert_eval_eq("PossibleZeroQ[ArcTan[z] - I/2 Log[(I + z)/(I - z)], "
                   "Assumptions -> Im[z] == 0 && -1 < Re[z] < 1]", "True", 0);
}

/* ---- Reals via Im[x]==0 / Element[_,Reals] ---------------------------- */

void test_pzq_reals_via_im_zero(void) {
    assert_eval_eq("PossibleZeroQ[Conjugate[x + I y] - (x - I y), "
                   "Assumptions -> Element[x, Reals] && Element[y, Reals]]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Re[Exp[I theta]] - Cos[theta], Assumptions -> Im[theta] == 0]",
                   "True", 0);
    assert_eval_eq("PossibleZeroQ[Im[Sin[x + I y]] - Cos[x] Sinh[y], "
                   "Assumptions -> Im[x] == 0 && Im[y] == 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Re[Cos[I y]] - Cosh[y], Assumptions -> Im[y] == 0]",
                   "True", 0);
    assert_eval_eq("PossibleZeroQ[Im[x^3 - 2 x + 1], Assumptions -> Im[x] == 0]", "True", 0);
}

/* ---- Real value with a Re[x] sign constraint (coordination) ----------- */

void test_pzq_real_with_repart_sign(void) {
    assert_eval_eq("PossibleZeroQ[(x^a)^b - x^(a b), Assumptions -> Im[x] == 0 && Re[x] > 0]",
                   "True", 0);
    assert_eval_eq("PossibleZeroQ[Im[Gamma[z]], Assumptions -> Im[z] == 0 && Re[z] > 0]",
                   "True", 0);
}

/* ---- Complex domain / half-plane / strip ------------------------------ */

void test_pzq_complex_domain(void) {
    assert_eval_eq("PossibleZeroQ[-x + Sqrt[x^2], Assumptions -> Re[x] > 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Re[I z] + Im[z], Assumptions -> Element[z, Complexes]]",
                   "True", 0);
    assert_eval_eq("PossibleZeroQ[Abs[z]^2 - (Re[z]^2 + Im[z]^2), "
                   "Assumptions -> Element[z, Complexes]]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Arg[z] - ArcTan[Re[z], Im[z]], Assumptions -> Re[z] > 0]",
                   "True", 0);
    assert_eval_eq("PossibleZeroQ[Sqrt[z^2] - z, Assumptions -> Re[z] > 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Log[z1 z2] - (Log[z1] + Log[z2]), "
                   "Assumptions -> Re[z1] > 0 && Re[z2] > 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Log[Exp[z]] - z, Assumptions -> -Pi < Im[z] <= Pi]",
                   "True", 0);
    assert_eval_eq("PossibleZeroQ[Sign[z] - z/Abs[z], "
                   "Assumptions -> Re[z] != 0 || Im[z] != 0]", "True", 0);
}

/* ---- Ambient Assuming[] scope ----------------------------------------- */

void test_pzq_assuming_wrapper(void) {
    assert_eval_eq("Assuming[Re[x] > 0, PossibleZeroQ[-x + Sqrt[x^2]]]", "True", 0);
    assert_eval_eq("Assuming[x >= 0, PossibleZeroQ[Sqrt[x^2] - x]]", "True", 0);
    assert_eval_eq("Assuming[Element[n, Integers], PossibleZeroQ[Sin[n Pi]]]", "True", 0);
    /* $Assumptions is restored on exit -> no ambient assumption leaks. */
    assert_eval_eq("$Assumptions", "True", 0);
    assert_eval_eq("PossibleZeroQ[Sin[n Pi]]", "False", 0);
}

/* ---- Soundness: genuine non-identities stay False --------------------- */

void test_pzq_soundness_false(void) {
    /* Non-identities that must NOT be swallowed by the assumed region. */
    assert_eval_eq("PossibleZeroQ[x, Assumptions -> x > 0]", "False", 0);
    assert_eval_eq("PossibleZeroQ[x - 1, Assumptions -> x > 0]", "False", 0);
    assert_eval_eq("PossibleZeroQ[Sin[n Pi] + 1, Assumptions -> Element[n, Integers]]",
                   "False", 0);
    assert_eval_eq("PossibleZeroQ[Cos[n Pi], Assumptions -> Element[n, Integers]]",
                   "False", 0);
    /* Sqrt[x^2]-x is +|x|-x = -2x != 0 on x<0. */
    assert_eval_eq("PossibleZeroQ[Sqrt[x^2] - x, Assumptions -> x < 0]", "False", 0);
    assert_eval_eq("PossibleZeroQ[Arg[x], Assumptions -> x < 0]", "False", 0);
    assert_eval_eq("PossibleZeroQ[Re[z], Assumptions -> Element[z, Complexes]]", "False", 0);
}

/* ---- No-assumption path is unchanged ---------------------------------- */

void test_pzq_no_assumptions_regression(void) {
    assert_eval_eq("PossibleZeroQ[Sin[n Pi]]", "False", 0);
    assert_eval_eq("PossibleZeroQ[Sin[x]^2 + Cos[x]^2 - 1]", "True", 0);
    assert_eval_eq("PossibleZeroQ[x - x]", "True", 0);
    assert_eval_eq("PossibleZeroQ[0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[1]", "False", 0);
    assert_eval_eq("PossibleZeroQ[x]", "False", 0);
}

/* ---- Manual first-argument threading (PossibleZeroQ is not Listable) --- */

void test_pzq_threading(void) {
    assert_eval_eq("PossibleZeroQ[{0, 1, x - x}]", "{True, False, True}", 0);
    /* The Assumptions option is broadcast to each list element as a
     * conjunction, not mis-threaded against the list. */
    assert_eval_eq("PossibleZeroQ[{Sqrt[a^2] - a, a - a}, Assumptions -> a >= 0]",
                   "{True, True}", 0);
    assert_eval_eq("PossibleZeroQ[{Sqrt[a^2] - a, a + 1}, Assumptions -> a >= 0]",
                   "{True, False}", 0);
}

/* ---- Stress-derived hardening: branch cuts, special functions, poles ---- */

void test_pzq_branch_cuts(void) {
    /* Sqrt / Log / inverse-trig branch behaviour, both directions. */
    assert_eval_eq("PossibleZeroQ[Sqrt[z^2] + z, Assumptions -> Re[z] < 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Sqrt[z^2] - z, Assumptions -> Re[z] < 0]", "False", 0);
    assert_eval_eq("PossibleZeroQ[Sqrt[a] Sqrt[b] - Sqrt[a b], Assumptions -> a > 0 && b > 0]",
                   "True", 0);
    assert_eval_eq("PossibleZeroQ[Sqrt[a] Sqrt[b] - Sqrt[a b], Assumptions -> a < 0 && b < 0]",
                   "False", 0);
    assert_eval_eq("PossibleZeroQ[Log[z^2] - 2 Log[z], Assumptions -> Re[z] > 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Log[a b] - Log[a] - Log[b], Assumptions -> a < 0 && b < 0]",
                   "False", 0);
    assert_eval_eq("PossibleZeroQ[Log[Exp[z]] - z, Assumptions -> Im[z] > Pi]", "False", 0);
    assert_eval_eq("PossibleZeroQ[ArcTan[1/x] + ArcTan[x] - Pi/2, Assumptions -> x > 0]", "True", 0);
    assert_eval_eq("PossibleZeroQ[ArcTan[1/x] + ArcTan[x] - Pi/2, Assumptions -> x < 0]", "False", 0);
    assert_eval_eq("PossibleZeroQ[Im[Log[x]] - Pi, Assumptions -> x < 0]", "True", 0);
    /* Nonzero one-sided bound and Equal-pinning of a range channel. */
    assert_eval_eq("PossibleZeroQ[Sqrt[x - 3] Sqrt[x - 3] - (x - 3), Assumptions -> x > 3]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Log[Exp[z]] - z, Assumptions -> Im[z] == Pi]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Sqrt[x^2] - x, Assumptions -> x == 5]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Sqrt[x^2] - x, Assumptions -> x == -5]", "False", 0);
}

void test_pzq_special_function_poles(void) {
    /* A sample that numericalizes to a pole (ComplexInfinity) is definitively
     * non-zero, not "unknown" -> must be False, never collapse to True. */
    assert_eval_eq("PossibleZeroQ[Gamma[x], Assumptions -> Element[x, Integers] && x <= 0]",
                   "False", 0);
    assert_eval_eq("PossibleZeroQ[Gamma[-3]]", "False", 0);
    /* But cancellation of infinities (Indeterminate, not Infinity) stays True. */
    assert_eval_eq("PossibleZeroQ[Gamma[x] - Gamma[x], Assumptions -> Element[x, Integers] && x <= 0]",
                   "True", 0);
    /* 1/pole = 0, but 1/finite != 0. */
    assert_eval_eq("PossibleZeroQ[1/Gamma[x], Assumptions -> Element[x, Integers] && x <= 0]",
                   "True", 0);
    assert_eval_eq("PossibleZeroQ[1/Gamma[x], Assumptions -> Element[x, Integers] && x > 0]",
                   "False", 0);
    /* Special-function identities / guards under integer order. */
    assert_eval_eq("PossibleZeroQ[Fibonacci[n + 2] - Fibonacci[n + 1] - Fibonacci[n], "
                   "Assumptions -> Element[n, Integers]]", "True", 0);
    assert_eval_eq("PossibleZeroQ[Zeta[-2 n], Assumptions -> Element[n, Integers]]", "False", 0);
    assert_eval_eq("PossibleZeroQ[Zeta[2 n], Assumptions -> Element[n, Integers] && n > 0]",
                   "False", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_pzq_integer_domain);
    TEST(test_pzq_integer_domain_with_sign);
    TEST(test_pzq_real_sign);
    TEST(test_pzq_real_range);
    TEST(test_pzq_reals_via_im_zero);
    TEST(test_pzq_real_with_repart_sign);
    TEST(test_pzq_complex_domain);
    TEST(test_pzq_assuming_wrapper);
    TEST(test_pzq_soundness_false);
    TEST(test_pzq_no_assumptions_regression);
    TEST(test_pzq_threading);
    TEST(test_pzq_branch_cuts);
    TEST(test_pzq_special_function_poles);

    printf("All PossibleZeroQ assumptions tests passed!\n");
    return 0;
}
