#include "poly.h"
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void run_test(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* s = expr_to_string_fullform(res);
    if (strcmp(s, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, s);
    }
    ASSERT_MSG(strcmp(s, expected) == 0, "%s: expected %s, got %s", input, expected, s);
    free(s);
    expr_free(res);
    expr_free(e);
}

void test_polynomialq_basic() {
    run_test("PolynomialQ[x, {x}]", "True");
    run_test("PolynomialQ[a + b, x]", "True");
    run_test("PolynomialQ[x^3 - 2x/y + 3x z, x]", "True");
    run_test("PolynomialQ[x^3 - 2x/y + 3x z, y]", "False");
    run_test("PolynomialQ[x^2 + a x y^2 - b Sin[c], {x, y}]", "True");
}

void test_polynomialq_advanced() {
    run_test("PolynomialQ[f[a] + f[a]^2, f[a]]", "True");
    run_test("PolynomialQ[x + 1/2, x]", "True");
    run_test("PolynomialQ[x + I, x]", "True");
    run_test("PolynomialQ[x^2.0, x]", "False");
    run_test("PolynomialQ[Sqrt[x], x]", "False");
    run_test("PolynomialQ[1/x, x]", "False");
    run_test("PolynomialQ[x, x^2]", "False");
}

void test_variables() {
    run_test("Variables[(x+y)^2+3z^2-y z+7]", "List[x, y, z]");
    run_test("Variables[{x^2-a y,y^2-b z,z^2-c x}]", "List[a, b, c, x, y, z]");
    run_test("Variables[(a-b)/(x+y)-2/z]", "List[a, b, x, y, z]");
    run_test("Variables[Sqrt[x+y-z^2]+(-2 t)^(2/3)]", "List[t, x, y, z]");
    run_test("Variables[Sin[x]+Cos[x]]", "List[Cos[x], Sin[x]]");
    run_test("Variables[E^x]", "List[]");
}

void test_level() {
    run_test("Level[a+f[x,y^n], {-1}]", "List[a, x, y, n]");
    run_test("Level[a+f[x,y^n], 2]", "List[a, x, Power[y, n], f[x, Power[y, n]]]");
    run_test("Level[a+f[x,y^n], {0, Infinity}]", "List[a, x, y, n, Power[y, n], f[x, Power[y, n]], Plus[a, f[x, Power[y, n]]]]");
    run_test("Level[{{{{a}}}}, 3]", "List[List[a], List[List[a]], List[List[List[a]]]]");
    run_test("Level[{{{{a}}}}, -2]", "List[List[a], List[List[a]], List[List[List[a]]]]");
    run_test("Level[h0[h1[h2[h3[a]]]], {0, -1}]", "List[a, h3[a], h2[h3[a]], h1[h2[h3[a]]], h0[h1[h2[h3[a]]]]]");
    run_test("Level[x^2+y^3, {1}, Heads->True]", "List[Plus, Power[x, 2], Power[y, 3]]");
    run_test("Level[x^2+y^3, 3, Heads->True]", "List[Plus, Power, x, 2, Power[x, 2], Power, y, 3, Power[y, 3]]");
}

void test_expand_coeff() {
    run_test("Expand[(x+1)^3]", "Plus[1, Times[3, x], Times[3, Power[x, 2]], Power[x, 3]]");
    run_test("Coefficient[(x+1)^3, x, 2]", "3");
    run_test("Coefficient[a x + b y + c, x]", "a");
    run_test("Coefficient[a x^3 + b x^2 + c x + d, x, 2]", "b");
    run_test("Coefficient[(x+2)^2 + (x+3)^3, x, 0]", "31");
    run_test("Coefficient[x^s x, x^s]", "x");
    run_test("Coefficient[x^s x^t, x^s]", "Power[x, t]");
    run_test("Coefficient[(x+y)^4, x y^3]", "4");
}

void test_polynomialgcd() {
    /* GCDs come back expanded (the FLINT fmpq_mpoly primitive, used unconditionally
     * — no degree gate), which is WL-faithful; integer content is preserved. */
    run_test("PolynomialGCD[(1+x)^2(2+x)(4+x), (1+x)(2+x)(3+x)]", "Plus[2, Times[3, x], Power[x, 2]]");
    run_test("PolynomialGCD[(x+y)^2(y+z), (x+y)(y+z)^3]", "Plus[Times[x, y], Power[y, 2], Times[x, z], Times[y, z]]");
    run_test("PolynomialGCD[x^2+4x+4, x^2+2x+1]", "1");
    run_test("PolynomialGCD[x^4-4, x^4+4 x^2+4]", "Plus[2, Power[x, 2]]");
    run_test("PolynomialGCD[x^2+2 x y+y^2, x^3+y^3]", "Plus[x, y]");
    run_test("PolynomialGCD[x^2-1, x^3-1, x^4-1, x^5-1, x^6-1, x^7-1]", "Plus[-1, x]");
    run_test("PolynomialGCD[x^2+1, x^3+3I x^2-3x-I]", "Plus[Complex[0, 1], x]");
    /* Rational-function coefficients (RISCH_AUDIT_A4.md): a polynomial in the GCD
     * variable whose coefficients are rational functions of the other variables —
     * produced by the recursive Risch tower field integrator.  A denominator is a
     * unit in the fraction field, so it does not change the polynomial GCD; these
     * formerly returned unevaluated. */
    run_test("PolynomialGCD[(2 t^2 - 1)/x, t]", "1");
    run_test("PolynomialGCD[t^2, 2 t/x]", "t");
    run_test("PolynomialGCD[(t^2 - 1)/x, (t - 1)/x^2]", "Plus[-1, t]");
    run_test("PolynomialGCD[(x + y)/z, (x^2 - y^2)/z^2]", "Plus[x, y]");
    run_test("PolynomialGCD[6, 4]", "2");             /* numeric denominators unaffected */
}

void test_polynomiallcm() {
    run_test("PolynomialLCM[(1+x)^2(2+x)(4+x), (1+x)(2+x)(3+x)]", "Times[Plus[2, x], Plus[3, x], Plus[4, x], Power[Plus[1, x], 2]]");
    run_test("PolynomialLCM[(1+x)^2, (1+x)(2+x), (3+x)]", "Times[Plus[2, x], Plus[3, x], Power[Plus[1, x], 2]]");
    run_test("PolynomialLCM[(1+x)^2(2+y), (1+x)(2+y)(x+y)]", "Times[Power[Plus[1, x], 2], Plus[2, y], Plus[x, y]]");
    run_test("PolynomialLCM[x^4-4, x^4+4 x^2+4]", "Times[Plus[-2, Power[x, 2]], Plus[4, Times[4, Power[x, 2]], Power[x, 4]]]");
    run_test("PolynomialLCM[x^2+2 x y+y^2, x^3+y^3]", "Times[Plus[x, y], Plus[Power[x, 3], Power[y, 3]]]");
    run_test("PolynomialLCM[x^2-1, x^3-1, x^4-1, x^5-1, x^6-1, x^7-1]", "Times[Plus[1, x], Plus[-1, Power[x, 3]], Plus[1, Power[x, 2]], Plus[1, x, Power[x, 2], Power[x, 3], Power[x, 4]], Plus[1, x, Power[x, 2], Power[x, 3], Power[x, 4], Power[x, 5], Power[x, 6]], Plus[1, Times[-1, x], Power[x, 2]]]");
    run_test("PolynomialLCM[(x-1)(x-2)/(x-4), (x-1)/((x-4)(x-6))]", "Times[Plus[-2, x], Plus[-1, x], Power[Plus[-4, x], -1]]");
}

void test_polynomial_div_rem() {
    run_test("PolynomialQuotient[x^4+2x+1, x^2+1, x]", "Plus[-1, Power[x, 2]]");
    run_test("PolynomialQuotient[x^2, x+a, x]", "Plus[Times[-1, a], x]");
    run_test("PolynomialQuotient[x^2+2x+1, x^3, x]", "0");
    run_test("PolynomialQuotient[x^2+x+1, 2x+1, x]", "Plus[Rational[1, 4], Times[Rational[1, 2], x]]");
    run_test("PolynomialQuotient[x^2+ b x+1, a x+1, x]", "Plus[Times[-1, Power[a, -2]], Times[Power[a, -1], b], Times[Power[a, -1], x]]");
    
    run_test("PolynomialRemainder[x^4+2x+1, x^2+1, x]", "Plus[2, Times[2, x]]");
    run_test("PolynomialRemainder[x^2+2x+1, x+1, x]", "0");
    run_test("PolynomialRemainder[x^3, a x+b, x]", "Times[-1, Power[a, -3], Power[b, 3]]");
    run_test("PolynomialRemainder[x^2+x+1, 2x+1, x]", "Rational[3, 4]");
    run_test("PolynomialRemainder[x^2+ b x+1, a x+1, x]", "Plus[1, Power[a, -2], Times[-1, Times[Power[a, -1], b]]]");

    /* Radical (algebraic-number) coefficients: the symbolic-q_coeff branch
     * subtracts coefficient-wise so the leading term cancels via the field
     * relation Sqrt[2] Sqrt[3] = Sqrt[6]. The old whole-polynomial
     * Together/Cancel path is blind to such relations over a radical tower
     * and spins forever (the IntRationalLogPart hang). */
    run_test("PolynomialRemainder[x^2 - Sqrt[6], Sqrt[2] x - Sqrt[3], x]",
             "Plus[Rational[3, 2], Times[-1, Power[6, Rational[1, 2]]]]");
    run_test("PolynomialRemainder[x^2 - (Sqrt[2]+Sqrt[3]) x + Sqrt[6], x - Sqrt[2], x]", "0");
}

void test_collect() {
    run_test("Collect[b x^2+5x+7x^2+9a x+2, x]", "Plus[2, Times[Plus[5, Times[9, a]], x], Times[Plus[7, b], Power[x, 2]]]");
    run_test("Collect[a x+b y+c x+d y, x]", "Plus[Times[Plus[a, c], x], Times[b, y], Times[d, y]]");
    run_test("Collect[a x+b y+c x+d y, y]", "Plus[Times[a, x], Times[c, x], Times[Plus[b, d], y]]");
    run_test("Collect[(1+a+x)^3, x]", "Plus[1, Times[3, a], Times[3, Power[a, 2]], Power[a, 3], Times[Plus[3, Times[6, a], Times[3, Power[a, 2]]], x], Times[Plus[3, Times[3, a]], Power[x, 2]], Power[x, 3]]");
    run_test("Collect[a x^4+b x^4+2a^2 x-3b x+x-7, x]", "Plus[-7, Times[Plus[1, Times[2, Power[a, 2]], Times[-3, b]], x], Times[Plus[a, b], Power[x, 4]]]");
    run_test("Collect[a Sqrt[x]+Sqrt[x]+x^(2/3)-c x+3x-2b x^(2/3)+5, x]", "Plus[5, Times[Plus[1, a], Power[x, Rational[1, 2]]], Times[Plus[1, Times[-2, b]], Power[x, Rational[2, 3]]], Times[Plus[3, Times[-1, c]], x]]");
    run_test("Collect[(x y+x z+y z+x+y)^3, {x,y}]", "Plus[Times[Power[y, 3], Plus[1, Times[3, z], Times[3, Power[z, 2]], Power[z, 3]]], Times[x, Plus[Times[Power[y, 2], Plus[3, Times[9, z], Times[9, Power[z, 2]], Times[3, Power[z, 3]]]], Times[Power[y, 3], Plus[3, Times[6, z], Times[3, Power[z, 2]]]]]], Times[Power[x, 2], Plus[Times[y, Plus[3, Times[9, z], Times[9, Power[z, 2]], Times[3, Power[z, 3]]]], Times[Power[y, 2], Plus[6, Times[12, z], Times[6, Power[z, 2]]]], Times[Power[y, 3], Plus[3, Times[3, z]]]]], Times[Power[x, 3], Plus[1, Power[y, 3], Times[3, z], Times[3, Power[z, 2]], Power[z, 3], Times[y, Plus[3, Times[6, z], Times[3, Power[z, 2]]]], Times[Power[y, 2], Plus[3, Times[3, z]]]]]]");
    run_test("Collect[(1-x-(1+a)^2 x^2)^2, x]", "Plus[1, Times[-2, x], Times[Plus[1, Times[-2, Power[Plus[1, a], 2]]], Power[x, 2]], Times[2, Times[Power[Plus[1, a], 2], Power[x, 3]]], Times[Power[Plus[1, a], 4], Power[x, 4]]]");
}

void test_coefficientlist() {
    run_test("CoefficientList[1 + 6 x - x^4, x]", "List[1, 6, 0, 0, -1]");
    run_test("CoefficientList[(1 + x)^10, x]", "List[1, 10, 45, 120, 210, 252, 210, 120, 45, 10, 1]");
    run_test("CoefficientList[1 + a x^2 + b x y + c y^2, {x, y}]", "List[List[1, 0, c], List[0, b, 0], List[a, 0, 0]]");
    run_test("CoefficientList[(2 x + 3)^5, x]", "List[243, 810, 1080, 720, 240, 32]");
    run_test("CoefficientList[a x^4 + b x^3 + c x^2 + d x + e, x]", "List[e, d, c, b, a]");
    run_test("CoefficientList[(3 x + 4 y + 1)^3, {x, y}]", "List[List[1, 12, 48, 64], List[9, 72, 144, 0], List[27, 108, 0, 0], List[27, 0, 0, 0]]");
    run_test("CoefficientList[(x + y + z + 1) (2 x + 3 y^2 + 4 z^3 + 5), {x, y, z}]", "List[List[List[5, 5, 0, 4, 4], List[5, 0, 0, 4, 0], List[3, 3, 0, 0, 0], List[3, 0, 0, 0, 0]], List[List[7, 2, 0, 4, 0], List[2, 0, 0, 0, 0], List[3, 0, 0, 0, 0], List[0, 0, 0, 0, 0]], List[List[2, 0, 0, 0, 0], List[0, 0, 0, 0, 0], List[0, 0, 0, 0, 0], List[0, 0, 0, 0, 0]]]");
}

void test_decompose() {
    run_test("Decompose[x^2 + 1, x]", "List[Plus[1, x], Power[x, 2]]");
    run_test("Decompose[x^4 + x^2, x]", "List[Plus[x, Power[x, 2]], Power[x, 2]]");
    run_test("Decompose[(x + x^2)^2, x]", "List[Power[x, 2], Plus[x, Power[x, 2]]]");
    run_test("Decompose[(x^2 + x)^4 + 1, x]", "List[Plus[1, x], Power[x, 4], Plus[x, Power[x, 2]]]");
    run_test("Decompose[(x^2 + 1)^4 + 3, x]", "List[Plus[4, Times[2, x], Power[x, 2]], Plus[Times[2, x], Power[x, 2]], Power[x, 2]]");
    run_test("Decompose[2 x + 1, x]", "List[Plus[1, Times[2, x]]]");
    run_test("Decompose[(a x^3 + 1)^2 + b, x]", "List[Plus[1, b, Times[2, Times[a, x]], Times[Power[a, 2], Power[x, 2]]], Power[x, 3]]");
}

void test_hornerform() {
    run_test("HornerForm[11 x^3 - 4 x^2 + 7 x + 2]", "Plus[2, Times[x, Plus[7, Times[x, Plus[-4, Times[11, x]]]]]]");
    run_test("HornerForm[a + b x + c x^2, x]", "Plus[a, Times[x, Plus[b, Times[c, x]]]]");
    run_test("HornerForm[(11 x^3 - 4 x^2 + 7 x + 2)/(x^2 - 3 x + 1)]", "Times[Plus[2, Times[x, Plus[7, Times[x, Plus[-4, Times[11, x]]]]]], Power[Plus[1, Times[x, Plus[-3, x]]], -1]]");
    run_test("HornerForm[x y + 2 x^2 y + 2 x y^2 + 4 x^2 y^2, {x, y}]", "Times[x, Plus[Times[y, Plus[1, Times[2, y]]], Times[x, y, Plus[2, Times[4, y]]]]]");
    run_test("HornerForm[x y + 2 x^2 y + 2 x y^2 + 4 x^2 y^2, {y, x}]", "Times[y, Plus[Times[x, Plus[1, Times[2, x]]], Times[x, Plus[2, Times[4, x]], y]]]");
}

void test_resultant() {
    // run_test("Resultant[(x-a)(x-b),(x-c)(x-d)(x-e),x]", ...); 
    // The explicit Mathematica string given was: -(b-c) (-a+c) (b-d) (-a+d) (b-e) (-a+e).
    // Due to Expand/Sort we just test Resultant correctness on expanded output.
    
    // Let's use simple numeric tests to verify Resultant correctness exactly.
    run_test("Resultant[x^2 - 2 x + 7, x^3 - x + 5, x]", "265");
    run_test("Resultant[x^3 - 5 x^2 - 7 x + 3, x^3 - 8 x^2 + 9 x - 11, x]", "-10321");
    run_test("Resultant[x^3 - 5 x^2 - 7 x + 14, x^3 - 8 x^2 + 9 x + 58, x]", "0");
    run_test("Resultant[a x + b, c x + d, x]", "Plus[Times[-1, Times[b, c]], Times[a, d]]");
    run_test("Resultant[a, b, x]", "1");
    run_test("Resultant[a x + b, c, x]", "c");
    run_test("Resultant[a, b x + c, x]", "a");
}

void test_discriminant() {
    run_test("Discriminant[a x^2 + b x + c, x]", "Plus[Power[b, 2], Times[-4, Times[a, c]]]");
    run_test("Discriminant[5 x^4 - 3 x + 9, x]", "23273325");
    run_test("Discriminant[(x-1)(x-2)(x-3), x]", "4");
    run_test("Discriminant[(x-1)(x-2)(x-1), x]", "0");
}

void test_polynomialextendedgcd() {
    run_test("PolynomialExtendedGCD[2x^5-2x, (x^2-1)^2, x]", "List[Plus[-1, Power[x, 2]], List[Times[Rational[1, 4], x], Times[Rational[1, 2], Plus[-2, Times[-1, Power[x, 2]]]]]]");
    run_test("PolynomialExtendedGCD[a (x+b)^2, (x+a)(x+b), x]", "List[Plus[b, x], List[Times[-1, Power[Plus[Power[a, 2], Times[-1, Times[a, b]]], -1]], Power[Plus[a, Times[-1, b]], -1]]]");
    run_test("PolynomialExtendedGCD[(x-1)(x-2)^2, (x-1)(x^2-3), x]", "List[Plus[-1, x], List[Plus[7, Times[4, x]], Plus[9, Times[-4, x]]]]");
    run_test("PolynomialExtendedGCD[x^2+x+1, x+1, x, Modulus -> 2]", "List[1, List[1, x]]");
}

void test_power_simplification() {
    run_test("Power[Rational[1, 4], -1]", "4");
    run_test("Power[Rational[1, 4], -2]", "16");
    run_test("Power[Rational[2, 3], 2]", "Rational[4, 9]");
    run_test("Power[Rational[2, 3], -2]", "Rational[9, 4]");
}

void test_bigint_poly() {
    run_test("PolynomialGCD[104857600000000000000000000 (x^2 - 1), 104857600000000000000000000 (x^3 - 1)]", "Plus[-104857600000000000000000000, Times[104857600000000000000000000, x]]");
    run_test("PolynomialQuotient[104857600000000000000000000 x^2 - 104857600000000000000000000, 104857600000000000000000000 x - 104857600000000000000000000, x]", "Plus[1, x]");
    run_test("PolynomialRemainder[104857600000000000000000000 x^2 - 104857600000000000000000000, x - 1, x]", "0");
}

int main() {
    setbuf(stdout, NULL);
    printf("Starting poly_tests\n");
    symtab_init();
    core_init();

    TEST(test_power_simplification);
    TEST(test_polynomialq_basic);
    TEST(test_polynomialq_advanced);
    TEST(test_variables);
    TEST(test_level);
    TEST(test_expand_coeff);
    TEST(test_polynomialgcd);
    TEST(test_polynomialextendedgcd);
    TEST(test_polynomiallcm);
    TEST(test_polynomial_div_rem);
    TEST(test_collect);
    TEST(test_coefficientlist);
    TEST(test_decompose);
    TEST(test_hornerform);
    TEST(test_resultant);
    TEST(test_bigint_poly);
    TEST(test_discriminant);
    TEST(test_polynomialextendedgcd);
    
    printf("All polynomial tests passed!\n");
    return 0;
}
