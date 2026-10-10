#include "test_utils.h"
#include "symtab.h"
#include "core.h"

/*
 * Tests for inverse-trig and inverse-hyperbolic simplification. Three
 * families of identity are exercised:
 *
 *   1. Negation reflection           e.g. ArcCos[-x] -> Pi - ArcCos[x]
 *   2. Imaginary-axis bridge         e.g. ArcSin[I*x] -> I ArcSinh[x]
 *   3. Forward-of-inverse algebraic  e.g. Sinh[ArcCosh[x]] -> Sqrt[x-1]*Sqrt[x+1]
 *
 * Plus two assumption-driven sum / composite identities:
 *
 *   4. ArcTan[x] + ArcTan[1/x] -> +/-Pi/2 under x > 0 / x < 0
 *   5. ArcCosh[2 x^2 - 1] -> 2 ArcCosh[x] under x > 0
 *
 * The auto-evaluation rules in trig.c / hyperbolic.c handle (1)-(3) at
 * eval-step time; the assumption rules in simp.c handle (4)-(5).
 */

/* ---- Negation reflection (auto-eval) ---- */

void test_invtrig_arcsin_neg_odd(void) {
    assert_eval_eq("ArcSin[-x]", "-ArcSin[x]", 0);
}

void test_invtrig_arctan_neg_odd(void) {
    assert_eval_eq("ArcTan[-x]", "-ArcTan[x]", 0);
}

void test_invtrig_arccos_neg_pi_minus(void) {
    assert_eval_eq("ArcCos[-x]", "Pi - ArcCos[x]", 0);
}

void test_invtrig_arccos_neg_sum_simplifies(void) {
    /* The user-reference identity: ArcCos[-x] + ArcCos[x] == Pi. */
    assert_eval_eq("Simplify[ArcCos[-x] + ArcCos[x] - Pi]", "0", 0);
}

void test_invtrig_arcsec_neg_pi_minus(void) {
    assert_eval_eq("ArcSec[-x]", "Pi - ArcSec[x]", 0);
}

void test_invhyp_arcsinh_neg_odd(void) {
    assert_eval_eq("ArcSinh[-x]", "-ArcSinh[x]", 0);
}

void test_invhyp_arctanh_neg_odd(void) {
    assert_eval_eq("ArcTanh[-x]", "-ArcTanh[x]", 0);
}

void test_invhyp_arccoth_neg_odd(void) {
    assert_eval_eq("ArcCoth[-x]", "-ArcCoth[x]", 0);
}

void test_invhyp_arccsch_neg_odd(void) {
    assert_eval_eq("ArcCsch[-x]", "-ArcCsch[x]", 0);
}

/* ---- Imaginary-axis bridge (auto-eval) ---- */

void test_invtrig_arcsin_imag_bridge(void) {
    /* ArcSin[I y] -> I ArcSinh[y]. */
    assert_eval_eq("Simplify[ArcSin[I*x] - I*ArcSinh[x]]", "0", 0);
}

void test_invtrig_arccos_imag_bridge(void) {
    /* ArcCos[I y] -> Pi/2 - I ArcSinh[y]. */
    assert_eval_eq("Simplify[ArcCos[I*x] - (Pi/2 - I*ArcSinh[x])]", "0", 0);
}

void test_invtrig_arctan_imag_bridge(void) {
    /* ArcTan[I y] -> I ArcTanh[y]. */
    assert_eval_eq("Simplify[ArcTan[I*x] - I*ArcTanh[x]]", "0", 0);
}

void test_invtrig_arccot_imag_bridge(void) {
    /* ArcCot[I y] -> -I ArcCoth[y]. */
    assert_eval_eq("ArcCot[I*x]", "-I ArcCoth[x]", 0);
}

void test_invtrig_arccsc_imag_bridge(void) {
    /* ArcCsc[I y] -> -I ArcCsch[y]. */
    assert_eval_eq("ArcCsc[I*x]", "-I ArcCsch[x]", 0);
}

void test_invhyp_arcsinh_imag_bridge(void) {
    /* ArcSinh[I y] -> I ArcSin[y]. */
    assert_eval_eq("ArcSinh[I*x]", "I ArcSin[x]", 0);
}

void test_invhyp_arctanh_imag_bridge(void) {
    /* ArcTanh[I y] -> I ArcTan[y]. */
    assert_eval_eq("ArcTanh[I*x]", "I ArcTan[x]", 0);
}

void test_invhyp_arccoth_imag_bridge(void) {
    assert_eval_eq("ArcCoth[I*x]", "-I ArcCot[x]", 0);
}

void test_invhyp_arccsch_imag_bridge(void) {
    assert_eval_eq("ArcCsch[I*x]", "-I ArcCsc[x]", 0);
}

/* ---- Forward-of-inverse trig ---- */

void test_invtrig_sin_arccos(void) {
    /* Sin[ArcCos[x]] = Sqrt[1 - x^2]. */
    assert_eval_eq("Simplify[Sin[ArcCos[x]] - Sqrt[1 - x^2]]", "0", 0);
}

void test_invtrig_cos_arcsin(void) {
    /* Cos[ArcSin[x]] = Sqrt[1 - x^2]. */
    assert_eval_eq("Simplify[Cos[ArcSin[x]] - Sqrt[1 - x^2]]", "0", 0);
}

void test_invtrig_tan_arcsin(void) {
    /* Tan[ArcSin[x]] = x / Sqrt[1 - x^2]. */
    assert_eval_eq("Simplify[Tan[ArcSin[x]] - x / Sqrt[1 - x^2]]", "0", 0);
}

void test_invtrig_tan_arccos(void) {
    /* Tan[ArcCos[x]] = Sqrt[1 - x^2] / x. */
    assert_eval_eq("Simplify[Tan[ArcCos[x]] - Sqrt[1 - x^2] / x]", "0", 0);
}

void test_invtrig_sin_arctan(void) {
    /* Sin[ArcTan[x]] = x / Sqrt[1 + x^2]. */
    assert_eval_eq("Simplify[Sin[ArcTan[x]] - x / Sqrt[1 + x^2]]", "0", 0);
}

void test_invtrig_cos_arctan(void) {
    /* Cos[ArcTan[x]] = 1 / Sqrt[1 + x^2]. */
    assert_eval_eq("Simplify[Cos[ArcTan[x]] - 1 / Sqrt[1 + x^2]]", "0", 0);
}

void test_invtrig_tan_arccot(void) {
    /* Tan[ArcCot[x]] = 1/x. */
    assert_eval_eq("Tan[ArcCot[x]]", "1/x", 0);
}

void test_invtrig_cot_arctan(void) {
    /* Cot[ArcTan[x]] = 1/x. */
    assert_eval_eq("Cot[ArcTan[x]]", "1/x", 0);
}

/* ---- Forward-of-inverse hyperbolic ---- */

void test_invhyp_sinh_arccosh(void) {
    /* Sinh[ArcCosh[x]] = Sqrt[x-1] Sqrt[x+1] (principal branch). */
    assert_eval_eq("Sinh[ArcCosh[x]]", "Sqrt[-1 + x] Sqrt[1 + x]", 0);
}

void test_invhyp_sinh_arccosh_under_x_gt_1(void) {
    /* Under x > 1 the difference simplifies cleanly to 0. */
    assert_eval_eq("Assuming[x > 1, Simplify[Sinh[ArcCosh[x]] - Sqrt[x-1] Sqrt[x+1]]]",
                   "0", 0);
}

void test_invhyp_cosh_arcsinh(void) {
    /* Cosh[ArcSinh[x]] = Sqrt[1 + x^2]. */
    assert_eval_eq("Cosh[ArcSinh[x]]", "Sqrt[1 + x^2]", 0);
}

void test_invhyp_sinh_arctanh(void) {
    /* Sinh[ArcTanh[x]] = x / Sqrt[1 - x^2]. */
    assert_eval_eq("Simplify[Sinh[ArcTanh[x]] - x / Sqrt[1 - x^2]]", "0", 0);
}

void test_invhyp_cosh_arctanh(void) {
    /* Cosh[ArcTanh[x]] = 1 / Sqrt[1 - x^2]. */
    assert_eval_eq("Simplify[Cosh[ArcTanh[x]] - 1 / Sqrt[1 - x^2]]", "0", 0);
}

void test_invhyp_tanh_arcsinh(void) {
    /* Tanh[ArcSinh[x]] = x / Sqrt[1 + x^2]. */
    assert_eval_eq("Simplify[Tanh[ArcSinh[x]] - x / Sqrt[1 + x^2]]", "0", 0);
}

void test_invhyp_tanh_arccosh(void) {
    /* Tanh[ArcCosh[x]] = Sqrt[x-1] Sqrt[x+1] / x. */
    assert_eval_eq("Simplify[Tanh[ArcCosh[x]] - Sqrt[x-1] Sqrt[x+1] / x]", "0", 0);
}

void test_invhyp_tanh_arccoth(void) {
    /* Tanh[ArcCoth[x]] = 1/x. */
    assert_eval_eq("Tanh[ArcCoth[x]]", "1/x", 0);
}

/* ---- Direct-inverse identities (already exercised by Sin[ArcSin[x]] = x,
 *      sanity-checks on the auto-eval cascade) ---- */

void test_invtrig_sin_arcsin(void) {
    assert_eval_eq("Sin[ArcSin[x]]", "x", 0);
}

void test_invhyp_cosh_arccosh(void) {
    assert_eval_eq("Cosh[ArcCosh[x]]", "x", 0);
}

/* ---- Sum identity: ArcSin[x] + ArcCos[x] == Pi/2 ---- */

void test_invtrig_arcsin_plus_arccos(void) {
    assert_eval_eq("Simplify[ArcSin[x] + ArcCos[x] - Pi/2]", "0", 0);
}

/* ---- Reciprocal sum identity (assumption-driven) ---- */

void test_invtrig_arctan_reciprocal_no_assumption(void) {
    /* Without sign info the sum is not collapsible (the value is
     * Pi/2 for x > 0 but -Pi/2 for x < 0). */
    assert_eval_eq("Simplify[ArcTan[x] + ArcTan[1/x]]",
                   "ArcTan[x] + ArcTan[1/x]", 0);
}

void test_invtrig_arctan_reciprocal_positive(void) {
    /* Under x > 0 the sum reduces to Pi/2. */
    assert_eval_eq("Assuming[x > 0, Simplify[ArcTan[x] + ArcTan[1/x] - Pi/2]]",
                   "0", 0);
    assert_eval_eq("Assuming[x > 0, Simplify[ArcTan[x] + ArcTan[1/x]]]",
                   "1/2 Pi", 0);
}

void test_invtrig_arctan_reciprocal_negative(void) {
    /* Under x < 0 the sum reduces to -Pi/2. */
    assert_eval_eq("Assuming[x < 0, Simplify[ArcTan[x] + ArcTan[1/x] + Pi/2]]",
                   "0", 0);
    assert_eval_eq("Assuming[x < 0, Simplify[ArcTan[x] + ArcTan[1/x]]]",
                   "-1/2 Pi", 0);
}

/* ---- ArcCosh double-angle reduction (assumption-driven) ---- */

void test_invhyp_arccosh_double_angle(void) {
    /* Under x > 1 (or any x > 0): ArcCosh[2x^2 - 1] -> 2 ArcCosh[x]. */
    assert_eval_eq("Assuming[x > 1, Simplify[ArcCosh[2 x^2 - 1] - 2 ArcCosh[x]]]",
                   "0", 0);
    assert_eval_eq("Assuming[x > 0, Simplify[ArcCosh[2 x^2 - 1] - 2 ArcCosh[x]]]",
                   "0", 0);
}

/* ---- Numeric-argument auto-eval through the new folds ---- */

void test_invtrig_arccos_neg_half(void) {
    /* ArcCos[-1/2] should evaluate via the reflection rule + the existing
     * exact-value tables: Pi - ArcCos[1/2] = Pi - Pi/3 = 2 Pi/3. */
    assert_eval_eq("ArcCos[-1/2]", "2/3 Pi", 0);
}

void test_invtrig_arcsin_imag_two(void) {
    /* ArcSin[2 I] -> I ArcSinh[2] (no further numeric collapse since
     * ArcSinh[2] has no exact-value entry). */
    assert_eval_eq("ArcSin[2 I]", "I ArcSinh[2]", 0);
}

/* Complementary-angle identities ArcSin+ArcCos, ArcTan+ArcCot, ArcSec+ArcCsc
 * fold to Pi/2 (and to 0 against -Pi/2). */
void test_invtrig_complementary_pairs(void) {
    assert_eval_eq("Simplify[ArcTan[x] + ArcCot[x] - Pi/2]", "0", 0);
    assert_eval_eq("Simplify[ArcSin[x] + ArcCos[x] - Pi/2]", "0", 0);
    assert_eval_eq("Simplify[ArcSec[x] + ArcCsc[x] - Pi/2]", "0", 0);
    /* The bare sum collapses to the x-free constant Pi/2. */
    assert_eval_eq("FreeQ[Simplify[ArcTan[x] + ArcCot[x]], x]", "True", 0);
    /* Soundness: a non-complementary sum stays. */
    assert_eval_eq("Simplify[2 ArcTan[x] + ArcCot[x]]", "ArcCot[x] + 2 ArcTan[x]", 0);
}

/* Machin-type and ArcSin rational-angle addition reduce exactly via e^{i e}. */
void test_invtrig_rational_angle_addition(void) {
    assert_eval_eq("Simplify[4 ArcTan[1/5] - ArcTan[1/239] - Pi/4]", "0", 0);
    assert_eval_eq("Simplify[ArcSin[3/5] + ArcSin[5/13] - ArcSin[56/65]]", "0", 0);
    assert_eval_eq("Simplify[2 ArcTan[1/2] - ArcTan[4/3]]", "0", 0);
    assert_eval_eq("Simplify[ArcTan[1] - Pi/4]", "0", 0);
    /* Soundness: not an identity (2 ArcTan[1/2] = ArcTan[4/3], not 3/4). */
    assert_eval_eq("Simplify[2 ArcTan[1/2] - ArcTan[3/4]] === 0", "False", 0);
}

/* Conditional identities that reduce ONLY under the right Assumptions, via the
 * general assumption machinery (radical-product combine, inverse-of-forward,
 * ComplexExpand-under-reality, the derivative-constancy region engine, and the
 * pointwise hyperbolic-addition recognizer). */
void test_simplify_under_assumptions(void) {
    /* Radical-product combine (Reduce-backed sign proof). */
    assert_eval_eq("Simplify[Sinh[ArcCosh[x]] - Sqrt[x^2-1], x>1]", "0", 0);
    assert_eval_eq("Simplify[Tanh[ArcCosh[x]] - Sqrt[x^2-1]/x, x>1]", "0", 0);
    /* Inverse-of-forward (hyperbolic family). */
    assert_eval_eq("Simplify[ArcTanh[Tanh[x]] - x, Element[x, Reals]]", "0", 0);
    /* ComplexExpand under reality. */
    assert_eval_eq("Simplify[Abs[Sin[x + I y]]^2 - (Sin[x]^2 + Sinh[y]^2), Element[x|y, Reals]]", "0", 0);
    /* Region engine (derivative-constancy) + hyperbolic-addition recognizer. */
    assert_eval_eq("Simplify[ArcTan[x] + ArcTan[y] - ArcTan[(x+y)/(1-x y)], -1<x<1 && -1<y<1]", "0", 0);
    assert_eval_eq("Simplify[ArcTanh[x] + ArcTanh[y] - ArcTanh[(x+y)/(1+x y)], -1<x<1 && -1<y<1]", "0", 0);
    assert_eval_eq("Simplify[ArcSinh[x] + ArcSinh[y] - ArcSinh[x Sqrt[1+y^2] + y Sqrt[1+x^2]], Element[x|y, Reals]]", "0", 0);
    assert_eval_eq("Simplify[ArcCosh[x] + ArcCosh[y] - ArcCosh[x y + Sqrt[x^2-1] Sqrt[y^2-1]], x>1 && y>1]", "0", 0);
    assert_eval_eq("Simplify[ArcCosh[x] - ArcSinh[Sqrt[x^2-1]], x>1]", "0", 0);
}

/* SOUNDNESS: these must NEVER reduce to 0 (a false zero is a wrong answer). */
void test_simplify_assumptions_soundness(void) {
    /* No assumption. */
    assert_eval_eq("Simplify[Sinh[ArcCosh[x]] - Sqrt[x^2-1]] === 0", "False", 0);
    assert_eval_eq("Simplify[ArcTanh[Tanh[x]] - x] === 0", "False", 0);
    assert_eval_eq("Simplify[ArcTan[x] + ArcTan[y] - ArcTan[(x+y)/(1-x y)]] === 0", "False", 0);
    /* Wrong / too-large region (crosses x y = 1). */
    assert_eval_eq("Simplify[ArcTan[x] + ArcTan[y] - ArcTan[(x+y)/(1-x y)], x>2 && y>2] === 0", "False", 0);
    /* Both bases negative: radical combine unsound, must decline. */
    assert_eval_eq("Simplify[Sqrt[x-1] Sqrt[x+1] - Sqrt[x^2-1], x<-1] === 0", "False", 0);
    /* Genuine non-identities under plausible assumptions. */
    assert_eval_eq("Simplify[ArcTan[x] + ArcTan[y] - ArcTan[x+y], -1<x<1 && -1<y<1] === 0", "False", 0);
    assert_eval_eq("Simplify[ArcSinh[x] + ArcSinh[y] - ArcSinh[x+y], Element[x|y, Reals]] === 0", "False", 0);
    assert_eval_eq("Simplify[ArcCosh[x] + ArcCosh[y] - ArcCosh[x y - Sqrt[x^2-1] Sqrt[y^2-1]], x>1 && y>1] === 0", "False", 0);
}

/* v0.336: ArcF[F[x]] adjacent-branch / sub-interval reduction. */
void test_invtrig_arccos_subinterval(void) {
    /* (Pi/2,Pi) is inside [0,Pi] but not the literal bound fact. */
    assert_eval_eq("Simplify[ArcCos[Cos[x]] - x, Pi/2 < x < Pi]", "0", 0);
}
void test_invtrig_arcsin_reflect(void) {
    /* ArcSin[Sin[x]] = Pi - x on (Pi/2, Pi). */
    assert_eval_eq("Simplify[ArcSin[Sin[x]] - (Pi - x), Pi/2 < x < Pi]", "0", 0);
}
void test_invtrig_arctan_reflect(void) {
    /* ArcTan[Tan[x]] = x - Pi on (Pi/2, Pi). */
    assert_eval_eq("Simplify[ArcTan[Tan[x]] - (x - Pi), Pi/2 < x < Pi]", "0", 0);
}
void test_invtrig_arccos_descending_branch(void) {
    /* ArcCos[Cos[x]] = 2 Pi - x on (Pi, 2 Pi). */
    assert_eval_eq("Simplify[ArcCos[Cos[x]] - (2 Pi - x), Pi < x < 3 Pi/2]", "0", 0);
}
void test_invtrig_arcsin_higher_branch(void) {
    /* ArcSin[Sin[x]] = 3 Pi - x on (5Pi/2, 7Pi/2). */
    assert_eval_eq("Simplify[ArcSin[Sin[x]] - (3 Pi - x), 5 Pi/2 < x < 7 Pi/2]", "0", 0);
}
void test_invtrig_branch_straddle_declines(void) {
    /* Interval spanning a branch boundary must not reduce (soundness). */
    assert_eval_eq("Simplify[ArcCos[Cos[x]] - x, 0 < x < 3 Pi/2] === 0", "False", 0);
    assert_eval_eq("Simplify[ArcSin[Sin[x]] - x, 0 < x < Pi] === 0", "False", 0);
}

/* v0.337: Conjugate Schwarz reflection Conjugate[H[z]] -> H[Conjugate[z]]. */
void test_conjugate_tan_reflection(void) {
    assert_eval_eq("Simplify[Conjugate[Tan[x + I y]] - Tan[x - I y], "
                   "Element[x, Reals] && Element[y, Reals]]", "0", 0);
}
void test_conjugate_cot_sec_exp_reflection(void) {
    assert_eval_eq("Simplify[Conjugate[Cot[x + I y]] - Cot[x - I y], "
                   "Element[x, Reals] && Element[y, Reals]]", "0", 0);
    assert_eval_eq("Simplify[Conjugate[Sec[x + I y]] - Sec[x - I y], "
                   "Element[x, Reals] && Element[y, Reals]]", "0", 0);
    assert_eval_eq("Simplify[Conjugate[Exp[x + I y]] - Exp[x - I y], "
                   "Element[x, Reals] && Element[y, Reals]]", "0", 0);
    assert_eval_eq("Simplify[Conjugate[Sinh[x + I y]] - Sinh[x - I y], "
                   "Element[x, Reals] && Element[y, Reals]]", "0", 0);
}
void test_conjugate_reflection_no_assumption_declines(void) {
    /* Without reality of x, y the reflection must not fire. */
    assert_eval_eq("Simplify[Conjugate[Tan[x + I y]] - Tan[x - I y]] === 0", "False", 0);
}

int main(void) {
    symtab_init();
    core_init();

    /* Negation. */
    TEST(test_invtrig_arcsin_neg_odd);
    TEST(test_invtrig_arctan_neg_odd);
    TEST(test_invtrig_arccos_neg_pi_minus);
    TEST(test_invtrig_arccos_neg_sum_simplifies);
    TEST(test_invtrig_arcsec_neg_pi_minus);
    TEST(test_invhyp_arcsinh_neg_odd);
    TEST(test_invhyp_arctanh_neg_odd);
    TEST(test_invhyp_arccoth_neg_odd);
    TEST(test_invhyp_arccsch_neg_odd);

    /* Imaginary-axis bridge. */
    TEST(test_invtrig_arcsin_imag_bridge);
    TEST(test_invtrig_arccos_imag_bridge);
    TEST(test_invtrig_arctan_imag_bridge);
    TEST(test_invtrig_arccot_imag_bridge);
    TEST(test_invtrig_arccsc_imag_bridge);
    TEST(test_invhyp_arcsinh_imag_bridge);
    TEST(test_invhyp_arctanh_imag_bridge);
    TEST(test_invhyp_arccoth_imag_bridge);
    TEST(test_invhyp_arccsch_imag_bridge);

    /* Forward-of-inverse trig. */
    TEST(test_invtrig_sin_arccos);
    TEST(test_invtrig_cos_arcsin);
    TEST(test_invtrig_tan_arcsin);
    TEST(test_invtrig_tan_arccos);
    TEST(test_invtrig_sin_arctan);
    TEST(test_invtrig_cos_arctan);
    TEST(test_invtrig_tan_arccot);
    TEST(test_invtrig_cot_arctan);

    /* Forward-of-inverse hyperbolic. */
    TEST(test_invhyp_sinh_arccosh);
    TEST(test_invhyp_sinh_arccosh_under_x_gt_1);
    TEST(test_invhyp_cosh_arcsinh);
    TEST(test_invhyp_sinh_arctanh);
    TEST(test_invhyp_cosh_arctanh);
    TEST(test_invhyp_tanh_arcsinh);
    TEST(test_invhyp_tanh_arccosh);
    TEST(test_invhyp_tanh_arccoth);

    /* Direct-inverse cascade sanity. */
    TEST(test_invtrig_sin_arcsin);
    TEST(test_invhyp_cosh_arccosh);

    /* Sum / composite identities. */
    TEST(test_invtrig_arcsin_plus_arccos);
    TEST(test_invtrig_arctan_reciprocal_no_assumption);
    TEST(test_invtrig_arctan_reciprocal_positive);
    TEST(test_invtrig_arctan_reciprocal_negative);
    TEST(test_invhyp_arccosh_double_angle);

    /* Numeric collapses through the new folds. */
    TEST(test_invtrig_arccos_neg_half);
    TEST(test_invtrig_arcsin_imag_two);

    TEST(test_invtrig_complementary_pairs);
    TEST(test_invtrig_rational_angle_addition);
    TEST(test_simplify_under_assumptions);
    TEST(test_simplify_assumptions_soundness);

    /* v0.336 ArcF[F[x]] adjacent-branch / sub-interval reduction */
    TEST(test_invtrig_arccos_subinterval);
    TEST(test_invtrig_arcsin_reflect);
    TEST(test_invtrig_arctan_reflect);
    TEST(test_invtrig_arccos_descending_branch);
    TEST(test_invtrig_arcsin_higher_branch);
    TEST(test_invtrig_branch_straddle_declines);

    /* v0.337 Conjugate Schwarz reflection */
    TEST(test_conjugate_tan_reflection);
    TEST(test_conjugate_cot_sec_exp_reflection);
    TEST(test_conjugate_reflection_no_assumption_declines);

    printf("All inverse-trig Simplify tests passed!\n");
    return 0;
}
