/* Builtin-audit regressions: algebra (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* AlgebraicIntegerQ#1 (high): past to_qqbar's degree cap (120) the decline
 * became a definite False. Ring closure now decides sums, products and
 * positive (rational) powers of integers outright -- 2^(1/k) is a root of the
 * monic x^k - 2 -- and a cap decline is undecided (unevaluated), never False. */
static void test_algebraicintegerq_high_degree(void) {
    audit_true("AlgebraicIntegerQ#1",
               "Table[AlgebraicIntegerQ[2^(1/k)], {k, {50, 100, 200, 300, 500}}] === "
               "{True, True, True, True, True}");
    audit_true("AlgebraicIntegerQ#1",
               "Table[AlgebraicIntegerQ[Sum[Sqrt[Prime[k]], {k, n}]], {n, 4, 8}] === "
               "{True, True, True, True, True}");
    audit_true("AlgebraicIntegerQ#1", "AlgebraicIntegerQ[2^(1/1000) + 1]");
    /* A capped non-closure case is left undecided, not answered False. */
    audit_true("AlgebraicIntegerQ#1",
               "Head[AlgebraicIntegerQ[1/2^(1/200)]] === AlgebraicIntegerQ");
    /* Decisions under the cap are unchanged. */
    audit_chk("AlgebraicIntegerQ#1", "AlgebraicIntegerQ[Sqrt[2]/2]", "False");
    audit_chk("AlgebraicIntegerQ#1", "AlgebraicIntegerQ[(1 + Sqrt[5])/2]", "True");
    audit_chk("AlgebraicIntegerQ#1", "AlgebraicIntegerQ[1/2]", "False");
    audit_chk("AlgebraicIntegerQ#1", "AlgebraicIntegerQ[Pi]", "False");
}

/* IrreduciblePolynomialQ#1 (high): Extension -> All on a multivariate input
 * factored over Q(i) only. Now exact: irreducible over Q, then one Factor over
 * Q(t), t a root of a squarefree specialisation (a smooth point fixes the one
 * absolute factor through it, so that factor is defined over Q(t)). */
static void test_irreduciblepolynomialq_absolute(void) {
    audit_chk("IrreduciblePolynomialQ#1", "IrreduciblePolynomialQ[x^2 - 2 y^2, Extension -> All]", "False");
    audit_chk("IrreduciblePolynomialQ#1", "IrreduciblePolynomialQ[x^2 + x y + y^2, Extension -> All]", "False");
    audit_chk("IrreduciblePolynomialQ#1", "IrreduciblePolynomialQ[x^3 - 2 y^3, Extension -> All]", "False");
    /* x^2 y^2 - 2 = (x y - Sqrt[2])(x y + Sqrt[2]) */
    audit_chk("IrreduciblePolynomialQ#1", "IrreduciblePolynomialQ[x^2 y^2 - 2, Extension -> All]", "False");
    audit_chk("IrreduciblePolynomialQ#1", "IrreduciblePolynomialQ[x^2 + y^2, Extension -> All]", "False");
    /* Absolutely irreducible curves and a rank-3 quadric stay True. */
    audit_chk("IrreduciblePolynomialQ#1", "IrreduciblePolynomialQ[x^2 + y^3, Extension -> All]", "True");
    audit_chk("IrreduciblePolynomialQ#1", "IrreduciblePolynomialQ[x^4 + y^4 + 1, Extension -> All]", "True");
    audit_chk("IrreduciblePolynomialQ#1", "IrreduciblePolynomialQ[x^2 + y^2 + z^2, Extension -> All]", "True");
    audit_chk("IrreduciblePolynomialQ#1", "IrreduciblePolynomialQ[x y - 1, Extension -> All]", "True");
    /* Over Q (no extension) is unchanged. */
    audit_chk("IrreduciblePolynomialQ#1", "IrreduciblePolynomialQ[x^2 - 2 y^2]", "True");
}

/* Wave 3 S12/S13: declared argument counts and option checks
 * (src/arity.h, src/arity_table.c).  The evaluator checks a call
 * against its head's declared arity before the C code runs and
 * prints argx/argr/argrx/argt(u)/argb(u)/argm(u), nonopt or optx;
 * the call stays unevaluated. */
static void test_wave3_arity_options(void) {
    audit_chk("AlgebraicNumberNorm#1", "Quiet[Check[AlgebraicNumberNorm[Extension -> Sqrt[2]], fired, AlgebraicNumberNorm::argx]]", "fired");
    audit_chk("NumberFieldIntegralBasis#3", "Quiet[Check[NumberFieldIntegralBasis[Sqrt[2], 3], fired, NumberFieldIntegralBasis::argx]]", "fired");
    audit_chk("RootReduce#1", "Quiet[Check[RootReduce[Sqrt[2], Foo -> 1], fired, RootReduce::optx]]", "fired");
    audit_chk("RootReduce#1", "Quiet[RootReduce[Sqrt[2], Foo -> 1]]", "RootReduce[Sqrt[2], Foo -> 1]");
    audit_chk("SubresultantPolynomials#2", "Quiet[Check[SubresultantPolynomials[x^2 + 1, x + 1], fired, SubresultantPolynomials::argrx]]", "fired");
    audit_chk("Subresultants#3", "Quiet[Check[Subresultants[x^2 + 1, x + 1], fired, Subresultants::argrx]]", "fired");
    audit_chk("ToNumberField#1", "Quiet[Check[ToNumberField[1, 2, 3], fired, ToNumberField::argt]]", "fired");
    audit_chk("ToNumberField#1", "Quiet[ToNumberField[1, 2, 3]]", "ToNumberField[1, 2, 3]");
    audit_chk("Cancel#2", "Quiet[Check[Cancel[x, y], fired, Cancel::nonopt]]", "fired");
    audit_chk("Cancel#2", "Quiet[Cancel[x, y]]", "Cancel[x, y]");
    audit_chk("Cancel#2", "Quiet[Check[Cancel[(x^2 - 1)/(x - 1), Foo -> 1], fired, Cancel::optx]]", "fired");
}

/* Wave 4 (G5): Apart with transcendental kernels, Together Trig option,
 * Numerator/Denominator of bignum complex numbers, IrreduciblePolynomialQ
 * over GF(p). */
static void test_wave4_algebra(void) {
    audit_chk("Apart#1", "Apart[Sin[x]/(x^2 - 1), x]", "Sin[x]/(2 (-1 + x)) - Sin[x]/(2 (1 + x))");
    audit_chk("Apart#1", "Apart[1/(x^2 - 1) + Sin[x], x]", "1/(2 (-1 + x)) - 1/(2 (1 + x)) + Sin[x]");
    audit_chk("Apart#1", "Apart[E^x/(x^2 - 1), x]", "E^x/(2 (-1 + x)) - E^x/(2 (1 + x))");
    audit_chk("Together#1", "Together[Sin[x]/Cos[x] + 1, Trig -> True]", "Sec[x] (Cos[x] + Sin[x])");
    audit_chk("Together#1", "Together[x + 1/x, Trig -> True]", "(1 + x^2)/x");
    audit_chk("Together#1", "Together[x + 1/x, Trig -> False]", "(1 + x^2)/x");
    audit_chk("Denominator#1", "Denominator[(10^30 + I)/7^25]", "1341068619663964900807");
    audit_chk("Numerator#1", "Numerator[1/10^30 + I]", "1 + 1000000000000000000000000000000*I");
    audit_chk("Numerator#1", "Numerator[1/2 + I/3]", "3 + 2*I");
    audit_chk("Subresultants#1", "Subresultants[y[1]^2 + 1, y[1] + 2, y[1]]", "{5, 1}");
    audit_chk("Subresultants#1", "Subresultants[Sin[x]^2 + 1, Sin[x] + 1, Sin[x]]", "{2, 1}");
    audit_chk("AlgebraicNumberDenominator#1", "AlgebraicNumberDenominator[Root[#^5 - # - 1/7 &, 1]]", "7");
    audit_chk("AlgebraicNumberDenominator#1", "AlgebraicNumberDenominator[Root[#^3 - 1/2 &, 1]]", "2");
    audit_chk("AlgebraicNumberDenominator#1", "AlgebraicNumberDenominator[Root[7 #^5 - 7 # - 1 &, 1]]", "7");
    audit_chk("MinimalPolynomial#2", "MinimalPolynomial[Root[#^2 - y &, 1], x]", "MinimalPolynomial[Root[#1^2 - y &, 1], x]");
    audit_chk("IrreduciblePolynomialQ#2", "IrreduciblePolynomialQ[x^2 + 1, Modulus -> 2]", "False");
    audit_chk("IrreduciblePolynomialQ#2", "IrreduciblePolynomialQ[x^2 + 1, Modulus -> 3]", "True");
}

int main(void) {
    audit_init();
    TEST(test_algebraicintegerq_high_degree);
    TEST(test_irreduciblepolynomialq_absolute);
    TEST(test_wave3_arity_options);
    TEST(test_wave4_algebra);
    return audit_finish("audit_algebra");
}
