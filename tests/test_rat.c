#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void run_test(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* res_str = expr_to_string_fullform(res);
    if (strcmp(res_str, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, res_str);
        ASSERT_STR_EQ(res_str, expected);
    }
    free(res_str);
    expr_free(e);
    expr_free(res);
}

void test_rational_simplification() {
    run_test("Rational[1, 1]", "1");
    run_test("Rational[15, 5]", "3");
    run_test("Rational[15, -5]", "-3");
    run_test("Rational[-15, -5]", "3");
    run_test("Rational[2, 4]", "Rational[1, 2]");
    run_test("Rational[0, 4]", "0");
    run_test("Rational[1, 0]", "ComplexInfinity");
    run_test("Rational[0, 0]", "Indeterminate");
}

void test_numerator() {
    run_test("Numerator[2/3]", "2");
    run_test("Numerator[(x-1)(x-2)/(x-3)^2]", "Times[Plus[-2, x], Plus[-1, x]]");
    run_test("Numerator[Sin[x]^a (Sin[x]-2)^-b]", "Power[Sin[x], a]");
    run_test("Numerator[3/7]", "3");
    run_test("Numerator[3/7 + I/11]", "Complex[33, 7]");
    run_test("Numerator[(x-1)^2/((x-2)(x-3))]", "Power[Plus[-1, x], 2]");
    run_test("Numerator[a x^n y^-m Exp[a - b - 2c + 3d]]", "Times[a, Power[E, Plus[a, Times[3, d]]], Power[x, n]]");
    run_test("Numerator[a^-b/x]", "1");
    run_test("Numerator[2 x^y b^2]", "Times[2, Power[b, 2], Power[x, y]]");
    run_test("Numerator[{1,2,3,4,5,6}/3]", "List[1, 2, 1, 4, 5, 2]");
    run_test("Numerator[5]", "5");
}

void test_denominator() {
    run_test("Denominator[2/3]", "3");
    run_test("Denominator[(x-1)(x-2)/(x-3)^2]", "Power[Plus[-3, x], 2]");
    run_test("Denominator[Sin[x]^a (Sin[x]-2)^-b]", "Power[Plus[-2, Sin[x]], b]");
    run_test("Denominator[3/7]", "7");
    run_test("Denominator[3/7 + I/11]", "77");
    run_test("Denominator[(x-1)^2/((x-2)(x-3))]", "Times[Plus[-3, x], Plus[-2, x]]");
    run_test("Denominator[a x^n y^-m Exp[a - b - 2c + 3d]]", "Times[Power[E, Plus[b, Times[2, c]]], Power[y, m]]");
    run_test("Denominator[a^-b/x]", "Times[Power[a, b], x]");
    run_test("Denominator[2 x^y b^2]", "1");
    run_test("Denominator[{1,2,3,4,5,6}/3]", "List[3, 3, 1, 3, 3, 1]");
    run_test("Denominator[5]", "1");
}

void test_cancel() {
    run_test("Cancel[(x^2 - 1)/(x - 1)]", "Plus[1, x]");
    run_test("Cancel[(x - y)/(x^2 - y^2) + (x^3 - 27)/(x^2 - 9) + (x^3 + 1)/(x^2 - x + 1)]", "Plus[1, x, Times[Power[Plus[3, x], -1], Plus[9, Times[3, x], Power[x, 2]]], Power[Plus[x, y], -1]]");
    run_test("Cancel[(x - a)/(x^2 - a^2) == 0 && (x^2 - 2 x + 1)/(x - 1) >= 0]", "And[Equal[Power[Plus[a, x], -1], 0], GreaterEqual[Plus[-1, x], 0]]");
    run_test("Cancel[(x - 1)/(x^2 - 1) + (x - 2)/(x^2 - 4)]", "Plus[Power[Plus[1, x], -1], Power[Plus[2, x], -1]]");
}

/* Cancel must treat a single symbolic base with fractional rational exponents
 * as an algebraic generator: substitute y -> g^m where m = lcm(denoms),
 * cancel polynomially in g, then substitute back. */
void test_cancel_algebraic_generator() {
    /* (y - 1) / (Sqrt[y] - 1)  =  1 + Sqrt[y]  via y = g^2: (g^2-1)/(g-1) */
    run_test("Cancel[(y - 1)/(Sqrt[y] - 1)]",
             "Plus[1, Power[y, Rational[1, 2]]]");

    /* (y^(2/3) - 1) / (y^(1/3) - 1)  =  1 + y^(1/3)  via y = g^3 */
    run_test("Cancel[(y^(2/3) - 1)/(y^(1/3) - 1)]",
             "Plus[1, Power[y, Rational[1, 3]]]");

    /* (y - 1) / (y^(1/3) - 1)  =  1 + y^(1/3) + y^(2/3) */
    run_test("Cancel[(y - 1)/(y^(1/3) - 1)]",
             "Plus[1, Power[y, Rational[1, 3]], Power[y, Rational[2, 3]]]");

    /* (y^(2/3) - y^(1/3)) / (y^(1/3) - 1)  =  y^(1/3) */
    run_test("Cancel[(y^(2/3) - y^(1/3))/(y^(1/3) - 1)]",
             "Power[y, Rational[1, 3]]");

    /* (Sqrt[y] + 1) / (y - 1)  =  1 / (Sqrt[y] - 1)
       (numerator factor cancels with one factor of difference of squares) */
    run_test("Cancel[(Sqrt[y] + 1)/(y - 1)]",
             "Power[Plus[-1, Power[y, Rational[1, 2]]], -1]");

    /* Sqrt[y] (y - 1) / (Sqrt[y] - 1)  =  y + Sqrt[y]
       (cancellation interleaved with a non-cancelling Sqrt[y] factor) */
    run_test("Cancel[Sqrt[y] (y - 1)/(Sqrt[y] - 1)]",
             "Plus[Power[y, Rational[1, 2]], y]");

    /* (y^(1/2) - y^(1/3)) / (y^(1/6) - 1)  =  y^(1/3)
       (mixed denominators 2,3,6 -> generator g = y^(1/6)) */
    run_test("Cancel[(y^(1/2) - y^(1/3))/(y^(1/6) - 1)]",
             "Power[y, Rational[1, 3]]");

    /* Perfect-square trinomial in algebraic generator:
       (y^(2/3) + 2 y^(1/3) + 1) / (y^(1/3) + 1)  =  1 + y^(1/3) */
    run_test("Cancel[(y^(2/3) + 2 y^(1/3) + 1)/(y^(1/3) + 1)]",
             "Plus[1, Power[y, Rational[1, 3]]]");

    /* Inverted Plus with mixed g-exponents: 1/(y^(2/3) - 1/y^(1/3))
       substitutes to 1/(g^2 - 1/g), whose denominator is a Plus of terms
       with different g-denominators. Cancel must run together-style
       fraction combination on the substituted form (not just GCD on
       num/den extracted directly) to produce y^(1/3)/(y - 1). */
    run_test("Cancel[1/(y^(2/3) - 1/y^(1/3))]",
             "Times[Power[y, Rational[1, 3]], Power[Plus[-1, y], -1]]");

    /* y^(73/24) / (y^(2/3) - 1/y^(1/3)) -- same shape with a non-trivial
       numerator. Substituting y = g^24 gives g^73 / (g^16 - 1/g^8) =
       g^81 / (g^24 - 1), and back-substitution yields
       y^(81/24) / (y - 1) = y^(27/8) / (y - 1). */
    run_test("Cancel[y^(73/24) / (y^(2/3) - 1/y^(1/3))]",
             "Times[Power[y, Rational[27, 8]], Power[Plus[-1, y], -1]]");
}

/* The generator pass also handles non-symbol bases such as Log[r] and
 * Sqrt[Log[r]] -- any sub-expression with at least one fractional rational
 * exponent is treated as the algebraic generator (matched via expr_eq). */
void test_cancel_nonsymbol_base() {
    run_test("Cancel[(Log[r] - 1)/(Sqrt[Log[r]] - 1)]",
             "Plus[1, Power[Log[r], Rational[1, 2]]]");
    run_test("Cancel[(Log[r]^(2/3) - 1)/(Log[r]^(1/3) - 1)]",
             "Plus[1, Power[Log[r], Rational[1, 3]]]");
}

/* A symbolic power x^n whose base x is ALSO a bare polynomial generator makes
 * the generic multivariate PolynomialGCD's variable set algebraically
 * dependent; before the rat_has_dependent_power_generators soundness gate this
 * looped forever (e.g. reached from Integrate[1/(1+x^n), x]).  The gate leaves
 * such fractions uncancelled (correctness-preserving); genuinely-cancellable
 * commensurate generators still cancel via the FLINT path. */
void test_cancel_dependent_power_generators() {
    /* Terminates; nothing to cancel (x^(n-1) and x^n are incommensurate). */
    run_test("Cancel[x^(n-1)/(1+x^n)]",
             "Times[Power[x, Plus[-1, n]], Power[Plus[1, Power[x, n]], -1]]");
    /* Commensurate generator: still cancels (FLINT). */
    run_test("Cancel[(x^(2 n) - 1)/(x^n - 1)]", "Plus[1, Power[x, n]]");
    /* Shared symbolic power with a bare-x cofactor: still cancels (FLINT). */
    run_test("Cancel[x^n (x+1)/(x^n (x+2))]",
             "Times[Plus[1, x], Power[Plus[2, x], -1]]");
}

/* Together must apply the same algebraic-generator pass.
 * The headline regression: y^(5/8) (y^(19/8) - y^(73/24)/(y^(2/3) - 1/y^(1/3)))
 * with denominators {8, 8, 24, 3, 3} yields lcm m = 24, and after substitution
 * the rational expression in g collapses to -g^72/(g^24 - 1) = -y^3/(y - 1). */
void test_together_algebraic_generator() {
    run_test("Together[y^(5/8) (y^(19/8) - y^(73/24)/(y^(2/3) - 1/y^(1/3)))]",
             "Times[-1, Power[y, 3], Power[Plus[-1, y], -1]]");

    /* y^(2/3) - 1/y^(1/3)  =  (y - 1)/y^(1/3) */
    run_test("Together[y^(2/3) - 1/y^(1/3)]",
             "Times[Plus[-1, y], Power[y, Rational[-1, 3]]]");

    /* 1/y^(1/2) + 1/y^(1/3)  =  (1 + y^(1/6))/y^(1/2) */
    run_test("Together[1/y^(1/2) + 1/y^(1/3)]",
             "Times[Power[y, Rational[-1, 2]], Plus[1, Power[y, Rational[1, 6]]]]");

    /* Non-symbol base: Log[r]^(p/q) is treated identically. */
    run_test("Together[1/Log[r]^(1/2) + 1/Log[r]^(1/3)]",
             "Times[Power[Log[r], Rational[-1, 2]], Plus[1, Power[Log[r], Rational[1, 6]]]]");

    /* Constant-radical denominator with a NON-UNIT leading coefficient on the
     * polynomial variable: 1/(2x - Sqrt5) - 1/(2x + Sqrt5) = 2 Sqrt5/(4x^2 - 5).
     * Combining runs the together Plus-branch strict exact-division of the LCM
     * denominator (4x^2 - 5) by each factor (2x +/- Sqrt5); the leading-
     * coefficient recursion in exact_poly_div must divide the algebraic
     * coefficient -2 Sqrt5 by the rational 2. Before the fix that division
     * returned NULL (it required BOTH operands rational), so the whole sum was
     * left uncombined for any coefficient != 1. */
    run_test("Together[1/(2 x - Sqrt[5]) - 1/(2 x + Sqrt[5])]",
             "Times[2, Power[5, Rational[1, 2]], Power[Plus[-5, Times[4, Power[x, 2]]], -1]]");
    run_test("Together[1/(3 x - Sqrt[5]) - 1/(3 x + Sqrt[5])]",
             "Times[2, Power[5, Rational[1, 2]], Power[Plus[-5, Times[9, Power[x, 2]]], -1]]");
}

/* Together / Cancel over transcendental-kernel generators (Log, Exp,
 * inverse-trig): the FLINT plain-rational fast path treats each kernel as an
 * atomic generator instead of falling back to the ~50x slower symbolic GCD.
 * Output must match the classical Together/Cancel exactly. */
void test_together_cancel_kernels() {
    /* Log[x] kernel: combines into one reduced fraction. */
    run_test("Together[1/(1 + Log[x]) + 1/(2 + Log[x])]",
             "Times[Plus[3, Times[2, Log[x]]], Power[Plus[2, Power[Log[x], 2], Times[3, Log[x]]], -1]]");
    /* Inverse-trig kernel: same handling as Log. */
    run_test("Together[1/(1 + ArcTan[x]) + 1/(2 + ArcTan[x])]",
             "Times[Plus[3, Times[2, ArcTan[x]]], Power[Plus[2, Power[ArcTan[x], 2], Times[3, ArcTan[x]]], -1]]");
    /* Exp kernel: E^(k x) is treated as (E^x)^k, so the sum combines correctly
     * (the classical path left this uncombined — a latent bug now fixed). */
    run_test("Together[1/(1 + E^x) + 1/(2 + E^x)]",
             "Times[Plus[3, Times[2, Power[E, x]]], Power[Plus[2, Times[3, Power[E, x]], Power[E, Times[2, x]]], -1]]");
    /* Exp cross-power cancellation depends on E^(2x) = (E^x)^2. */
    run_test("Cancel[(E^(2 x) - 1)/(E^x - 1)]", "Plus[1, Power[E, x]]");
    run_test("Cancel[(Log[x]^2 - 1)/(Log[x] - 1)]", "Plus[1, Log[x]]");
    /* Nested independent Log kernels: Log[x] and Log[1+Log[x]] are distinct
     * generators. */
    run_test("Together[1/(1 + Log[1 + Log[x]]) + 1/Log[x]]",
             "Times[Plus[1, Log[x], Log[Plus[1, Log[x]]]], "
             "Power[Plus[Log[x], Times[Log[x], Log[Plus[1, Log[x]]]]], -1]]");
    /* Trig kernels are declined to the classical path (1/Sin -> Csc form). */
    run_test("Together[1/(1 + Cos[x]) + 1/Sin[x]]",
             "Times[Power[Plus[1, Cos[x]], -1], Plus[1, Cot[x], Csc[x]]]");
}

void test_together() {
    run_test("Together[a/b + c/d]", "Times[Power[b, -1], Power[d, -1], Plus[Times[b, c], Times[a, d]]]");
    run_test("Together[x^2/(x^2 - 1) + x/(x^2 - 1)]", "Times[x, Power[Plus[-1, x], -1]]");
    run_test("Together[1/x + 1/(x + 1) + 1/(x + 2) + 1/(x + 3)]", "Times[Plus[6, Times[22, x], Times[18, Power[x, 2]], Times[4, Power[x, 3]]], Power[Plus[Times[6, x], Times[11, Power[x, 2]], Times[6, Power[x, 3]], Power[x, 4]], -1]]");
    run_test("Together[x^2/(x - y) - x y/(x - y)]", "x");
    run_test("Together[{1/x + 1/(x + 1), 1/(x + 2) + 1/(x + 3)}]", "List[Times[Power[Plus[x, Power[x, 2]], -1], Plus[1, Times[2, x]]], Times[Plus[5, Times[2, x]], Power[Plus[6, Times[5, x], Power[x, 2]], -1]]]");
    run_test("Together[(x - 1)/(x^2 - 1) + (x - 2)/(x^2 - 4)]", "Times[Plus[3, Times[2, x]], Power[Plus[2, Times[3, x], Power[x, 2]], -1]]");
}

int main() {
    symtab_init();
    core_init();

    printf("Running rat tests...\n");
    TEST(test_rational_simplification);
    TEST(test_numerator);
    TEST(test_denominator);
    TEST(test_cancel);
    TEST(test_cancel_algebraic_generator);
    TEST(test_cancel_nonsymbol_base);
    TEST(test_cancel_dependent_power_generators);
    TEST(test_together);
    TEST(test_together_algebraic_generator);
    TEST(test_together_cancel_kernels);
    printf("All rat tests passed!\n");

    return 0;
}
