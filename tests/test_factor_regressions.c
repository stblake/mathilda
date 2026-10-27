/*
 * test_factor_regressions.c
 * -------------------------
 * Regressions from the Factor family review: inputs that hung, crashed,
 * came back unfactored or produced a wrong answer.  Each case is an
 * expression that must evaluate to True (property checks are robust to
 * harmless reorderings), or an exact FullForm where the output convention
 * itself is what is being pinned.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

static int failures = 0;

/* Evaluate `input` and compare its FullForm with `expected`. */
static void check(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    if (!e) {
        printf("FAIL (parse): %s\n", input);
        failures++;
        return;
    }
    Expr* r = evaluate(e);
    char* got = expr_to_string_fullform(r);
    if (strcmp(got, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, got);
        failures++;
    } else {
        printf("PASS: %s\n", input);
    }
    free(got);
    expr_free(r);
    expr_free(e);
}

static void check_true(const char* input) { check(input, "True"); }

/* A univariate numerator or denominator of a rational function went to the
 * int64 Berlekamp-Zassenhaus path instead of FLINT: it hung on x^43 - 1,
 * crashed past degree 10000 and left wide-coefficient sides unfactored. */
static void test_rational_function_sides(void) {
    check_true("Length[FactorList[Denominator[Factor[1/(x^43 - 1)]]]] == 3");
    check_true("Length[FactorList[Denominator[Factor[(x + 1)/(x^43 - 1)]]]] == 3");
    check("Factor[Expand[(586 x - 137) (905 x + 872)]/(x + y)]",
          "Times[Plus[-137, Times[586, x]], Plus[872, Times[905, x]], "
          "Power[Plus[x, y], -1]]");
    check("Factor[(x - 1)/Expand[(7 x - 3 y) (4 + 3 y) (1 + x y)]]",
          "Times[Plus[-1, x], Power[Plus[Times[7, x], Times[-3, y]], -1], "
          "Power[Plus[4, Times[3, y]], -1], Power[Plus[1, Times[x, y]], -1]]");
    check_true("Length[FactorList[Numerator[Factor[(x^84 - 1)/(x + y)]]]] == 13");
    check_true("With[{p = Expand[Product[x^k + k, {k, 1, 12}]]}, "
               "Length[FactorList[Numerator[Factor[p/(x^5 + 7)]]]] == Length[FactorList[p]]]");
    check("Factor[(x^2 - 1)/(2 x + 2)]",
          "Times[Rational[1, 2], Plus[-1, x]]");
}

/* PolynomialGCD / PolynomialLCM took their integer content in int64: a
 * BigInt counted as 1 and a product of contents wrapped. */
static void test_polynomial_content_bigint(void) {
    check("PolynomialLCM[7*10^18 x, 13 y]",
          "Times[91000000000000000000, x, y]");
    check("PolynomialLCM[10^30 x, 3 x]",
          "Times[3000000000000000000000000000000, x]");
    check("PolynomialGCD[10^30 x, 4*10^25 x^2]",
          "Times[40000000000000000000000000, x]");
    check("PolynomialGCD[4, 4 + 4 x]", "4");
}

/* Inexact input goes through an exact round trip.  A coefficient below 1e-4
 * was rationalised to 0 (a wrong answer), and a real of magnitude >= 2^63
 * could not be made exact, which sent the round trip back into itself. */
static void test_inexact_round_trip(void) {
    check("Rationalize[0.000001]", "Rational[1, 1000000]");
    check("Rationalize[1.*^-8]", "Rational[1, 100000000]");
    check_true("Rationalize[1.*^-300] === 1.*^-300");
    check_true("Rationalize[1.*^-12] === 1.*^-12");
    check_true("Abs[(x /. Solve[Factor[x - 0.000001] == 0, x][[1]]) - 1.*^-6] < 1.*^-12");
    check_true("Head[Factor[x + 1.0*^19]] === Plus");
    check_true("Head[Factor[2.0*^19 x]] === Times");
}

/* FactorSquareFree runs through FLINT: monomial factors split out, factors
 * sign-normalised with the unit in front, rational functions decomposed. */
static void test_factorsquarefree_conventions(void) {
    check("FactorSquareFree[x + 7 x^2 + 12 x^3]",
          "Times[x, Plus[1, Times[7, x], Times[12, Power[x, 2]]]]");
    check("FactorSquareFree[Expand[2 (-2 - x)^3 (x^2 + 1)]]",
          "Times[-2, Power[Plus[2, x], 3], Plus[1, Power[x, 2]]]");
    check("FactorSquareFree[(x^2 - 1)^2/(x + 1)^3]",
          "Times[Power[Plus[-1, x], 2], Power[Plus[1, x], -1]]");
    check("FactorSquareFree[Expand[Pi (x + 1)^2]]",
          "Times[Pi, Power[Plus[1, x], 2]]");
    check("FactorSquareFree[Expand[(1 - x)^3 (1 + x)]]",
          "Times[-1, Power[Plus[-1, x], 3], Plus[1, x]]");
}

int main(void) {
    symtab_init();
    core_init();

    printf("Running factor regression tests...\n");
    test_rational_function_sides();
    test_polynomial_content_bigint();
    test_inexact_round_trip();
    test_factorsquarefree_conventions();

    if (failures) {
        printf("%d factor regression test(s) FAILED\n", failures);
        return 1;
    }
    printf("All factor regression tests passed!\n");
    return 0;
}
