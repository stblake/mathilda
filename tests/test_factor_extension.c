/*
 * test_factor_extension.c
 * -----------------------
 * Factoring over algebraic number fields (src/poly/flint_trager.c): the
 * Extension and GaussianIntegers options in any number of variables, the
 * default factorisation over Q(i) of a univariate polynomial with Gaussian
 * coefficients, and the FactorList pairs over an extension.
 *
 * Expected values were computed with the reference engine.  Where its unit
 * placement is reproduced, the answer is compared exactly (SameQ of the
 * evaluated forms).  Where the reference places a unit in a way no stated
 * rule recovers (Gaussian units absorbed into one factor, the sign of a
 * factor free of the first generator), the test checks the property
 * instead: the product expands back to the input, with the reference's
 * number of irreducible factors.
 *
 * Every check runs (no abort on the first failure); the exit status is the
 * failure count's truth value.
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

static int failures = 0;
static int passes = 0;

/* Evaluate `input` and require the printed result to be `True`. */
static void check_true(const char* input) {
    Expr* e = parse_expression(input);
    if (!e) {
        printf("FAIL (parse): %s\n", input);
        failures++;
        return;
    }
    Expr* r = evaluate(e);
    char* got = expr_to_string(r);
    if (strcmp(got, "True") != 0) {
        printf("FAIL: %s\n  got: %s\n", input, got);
        failures++;
    } else {
        passes++;
    }
    free(got);
    expr_free(r);
    expr_free(e);
}

/* Evaluate `input` and compare its printed form with `expected`. */
static void check(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    if (!e) {
        printf("FAIL (parse): %s\n", input);
        failures++;
        return;
    }
    Expr* r = evaluate(e);
    char* got = expr_to_string(r);
    if (strcmp(got, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, got);
        failures++;
    } else {
        passes++;
    }
    free(got);
    expr_free(r);
    expr_free(e);
}

/* The evaluated `call` is SameQ to the evaluated `ref`. */
static void same(const char* call, const char* ref) {
    char buf[2048];
    snprintf(buf, sizeof buf, "(%s) === (%s)", call, ref);
    check_true(buf);
}

/* Factor[p, opts] multiplies back to p and has `n` irreducible factors
 * (counted with multiplicity, constants excluded).  ToRadicals lets a Root
 * generator's powers cancel. */
static void prop(const char* p, const char* opts, int n) {
    char buf[4096];
    snprintf(buf, sizeof buf,
             "Module[{r = Factor[%s, %s]}, "
             "Expand[ToRadicals[r] - (%s)] === 0 && tfcount[r] == %d]",
             p, opts, p, n);
    check_true(buf);
}

static void define_helpers(void) {
    /* The factor count of check_ref (constants and numbers count 0). */
    Expr* e = parse_expression(
        "tfcount[e_] := Which["
        "Head[e] === Times, Total[tfcount /@ List @@ e], "
        "Head[e] === Power && IntegerQ[e[[2]]] && e[[2]] > 0 && !FreeQ[e[[1]], x | y | z | a], "
        "e[[2]] tfcount[e[[1]]], "
        "FreeQ[e, x | y | z | a], 0, True, 1]");
    expr_free(evaluate(e));
    expr_free(e);
}

/* ---- Extension -> radicals, several variables ---- */
static void test_multivariate_radicals(void) {
    same("Factor[x^2 - 2 y^2, Extension -> Sqrt[2]]",
         "-((-x + Sqrt[2]*y)*(x + Sqrt[2]*y))");
    same("Factor[2 x^2 - y^2, Extension -> Sqrt[2]]",
         "-1/2*((-2*x + Sqrt[2]*y)*(2*x + Sqrt[2]*y))");
    same("Factor[y^2 - 2 x^2, Extension -> Sqrt[2]]",
         "((-2*x + Sqrt[2]*y)*(2*x + Sqrt[2]*y))/2");
    same("Factor[x^2 - 2 y^2 z^2, Extension -> Sqrt[2]]",
         "-((-x + Sqrt[2]*y*z)*(x + Sqrt[2]*y*z))");
    same("Factor[x^4 - 4 y^4, Extension -> Sqrt[2]]",
         "-((-x + Sqrt[2]*y)*(x + Sqrt[2]*y)*(x^2 + 2*y^2))");
    same("Factor[x^4 - 4 y^4, Extension -> {Sqrt[2], I}]",
         "-((-x + Sqrt[2]*y)*(-I*x + Sqrt[2]*y)*(I*x + Sqrt[2]*y)*(x + Sqrt[2]*y))");
    same("Factor[x^3 + y^3, Extension -> (-1)^(1/3)]",
         "-((x + y)*(-x + (-1)^(1/3)*y)*(x + (-1 + (-1)^(1/3))*y))");
    same("Factor[x^3 - 2 y^3, Extension -> 2^(1/3)]",
         "-((-x + 2^(1/3)*y)*(x^2 + 2^(1/3)*x*y + 2^(2/3)*y^2))");
    same("Factor[x^2 - x y - y^2, Extension -> Sqrt[5]]",
         "-1/4*((2*x + (-1 + Sqrt[5])*y)*(-2*x + (1 + Sqrt[5])*y))");
    same("Factor[x^2 - 3 y^2 + 2 x - 6 y - 2, Extension -> Sqrt[3]]",
         "-((-1 + Sqrt[3] - x + Sqrt[3]*y)*(1 + Sqrt[3] + x + Sqrt[3]*y))");
    same("Factor[x^2 - 2 y^2 - 2 y - 1/2, Extension -> Sqrt[2]]",
         "-1/4*((Sqrt[2] - 2*x + 2*Sqrt[2]*y)*(Sqrt[2] + 2*x + 2*Sqrt[2]*y))");
    same("Factor[x^4 - 10 x^2 y^2 + y^4, Extension -> {Sqrt[2], Sqrt[3]}]",
         "(-x + (Sqrt[2] - Sqrt[3])*y)*(x + (Sqrt[2] - Sqrt[3])*y)*"
         "(-x + (Sqrt[2] + Sqrt[3])*y)*(x + (Sqrt[2] + Sqrt[3])*y)");
    same("Factor[(x^2 - 2 y^2)^2 (x + y), Extension -> Sqrt[2]]",
         "(x + y)*(-x + Sqrt[2]*y)^2*(x + Sqrt[2]*y)^2");
    same("Factor[3 x^2 - 6 y^2, Extension -> Sqrt[2]]",
         "-3*(-x + Sqrt[2]*y)*(x + Sqrt[2]*y)");
    same("Factor[x^2 y - 2 y^3, Extension -> Sqrt[2]]",
         "-(y*(-x + Sqrt[2]*y)*(x + Sqrt[2]*y))");
    same("Factor[a x^2 - 2 a y^2, Extension -> Sqrt[2]]",
         "-(a*(-x + Sqrt[2]*y)*(x + Sqrt[2]*y))");
    same("Factor[x^2 - Sqrt[2] y^2, Extension -> 2^(1/4)]",
         "-((-x + 2^(1/4)*y)*(x + 2^(1/4)*y))");
    same("Factor[x^2 - 3 y^2, Extension -> Sqrt[2]]", "x^2 - 3*y^2");
    same("Factor[x^2 - 2 x y + y^2 - 2 z^2, Extension -> Sqrt[2]]",
         "-((x - y + Sqrt[2]*z)*(-x + y + Sqrt[2]*z))");
    same("Factor[x^2 y - 2 y z^2, Extension -> Sqrt[2]]",
         "-(y*(-x + Sqrt[2]*z)*(x + Sqrt[2]*z))");
    /* generators the review found unrecognised, spelled other ways */
    same("Factor[x^2 - 2 y^2, Extension -> -Sqrt[2]]",
         "-((-x + Sqrt[2]*y)*(x + Sqrt[2]*y))");
    same("Factor[x^2 - 2 y^2, Extension -> Sqrt[8]]",
         "-((-x + Sqrt[2]*y)*(x + Sqrt[2]*y))");
    same("Factor[x^2 - 2 y^2, Extension -> 1 + Sqrt[2]]",
         "-((-x + Sqrt[2]*y)*(x + Sqrt[2]*y))");
    same("Factor[x^2 - 3 y^2, Extension -> {Sqrt[2], Sqrt[2] + Sqrt[3]}]",
         "-((-x + Sqrt[3]*y)*(x + Sqrt[3]*y))");
    /* the Root generator renders as itself, not in I and Sqrt[3] */
    prop("x^3 - y^3", "Extension -> Root[#^2 + # + 1 &, 1]", 3);
}

/* ---- Extension -> Automatic ---- */
static void test_automatic(void) {
    same("Factor[x^2 - 2 y^2, Extension -> Automatic]", "x^2 - 2*y^2");
    same("Factor[x^2 - 2 Sqrt[2] x y + 2 y^2, Extension -> Automatic]",
         "(-x + Sqrt[2]*y)^2");
    same("Factor[x^2 - 2 Sqrt[2] x y + 2 y^2]", "x^2 - 2*Sqrt[2]*x*y + 2*y^2");
    /* F2: two generators, the whole tower (the reference's sign differs) */
    prop("x^2 + (Sqrt[2] - Sqrt[3]) x y - Sqrt[6] y^2", "Extension -> Automatic", 2);
    same("Factor[x^2 + (Sqrt[2] - Sqrt[3]) x - Sqrt[6], Extension -> Automatic]",
         "-((Sqrt[3] - x)*(Sqrt[2] + x))");
    same("Factor[x^2 - 2 Sqrt[2] x + 2, Extension -> Automatic]", "(Sqrt[2] - x)^2");
    /* Gaussian coefficients make I part of the automatic field */
    prop("Expand[(2 x + 3) ((-2 + 2 I) x^2 + (3 - I) x + (-2 - 2 I))^2 "
         "((3 + 2 I) x^2 + (-1 + 2 I) x + (-1 + I))^2 (3 x^2 + 3 x + 3)^2]",
         "Extension -> Automatic", 7);
}

/* ---- Gaussian integers ---- */
static void test_gaussian(void) {
    same("Factor[x^2 + y^2, Extension -> I]", "(x - I*y)*(x + I*y)");
    same("Factor[x^2 + y^2, Extension -> {I}]", "(x - I*y)*(x + I*y)");
    same("Factor[x^2 + y^2, GaussianIntegers -> True]", "(x - I*y)*(x + I*y)");
    same("Factor[x^2 + 1, GaussianIntegers -> True]", "(-I + x)*(I + x)");
    same("Factor[x^4 - y^4, GaussianIntegers -> True]", "(x - y)*(x - I*y)*(x + I*y)*(x + y)");
    same("Factor[x^2 + 4 y^2, GaussianIntegers -> True]", "(x - (2*I)*y)*(x + (2*I)*y)");
    same("Factor[4 x^2 + y^2, GaussianIntegers -> True]", "(2*x - I*y)*(2*x + I*y)");
    same("Factor[x^2 + y^2 + z^2, GaussianIntegers -> True]", "x^2 + y^2 + z^2");
    same("Factor[x^2 + 2 x y + 2 y^2, GaussianIntegers -> True]",
         "(x + (1 - I)*y)*(x + (1 + I)*y)");
    same("Factor[x^2 - 2 I x y - y^2, GaussianIntegers -> True]", "(x - I*y)^2");
    same("Factor[x^4 + 1, GaussianIntegers -> True]", "(-I + x^2)*(I + x^2)");
    same("Factor[x^2 - 2, Extension -> Sqrt[2], GaussianIntegers -> True]",
         "-((Sqrt[2] - x)*(Sqrt[2] + x))");
    same("Factor[x^2 + y^2, Extension -> Sqrt[2], GaussianIntegers -> True]",
         "(x - I*y)*(x + I*y)");
    same("Factor[x^2 + y^2, GaussianIntegers -> False]", "x^2 + y^2");
    same("Factor[(1 + I) x^2 + (1 + I) y^2, GaussianIntegers -> True]",
         "(1 + I)*(x - I*y)*(x + I*y)");
    same("Factor[2 x^2 + 2 y^2, GaussianIntegers -> True]", "2*(x - I*y)*(x + I*y)");
    same("Factor[x^2 + y^2 + 2 x + 1, GaussianIntegers -> True]",
         "(1 + x - I*y)*(1 + x + I*y)");
    same("Factor[x^2 y^2 + 1, GaussianIntegers -> True]", "(-I + x*y)*(I + x*y)");
    same("Factor[x^6 - y^6, GaussianIntegers -> True]",
         "(x - y)*(x + y)*(x^2 - x*y + y^2)*(x^2 + x*y + y^2)");
    same("Factor[2 x^2 + y^2, GaussianIntegers -> True]", "2*x^2 + y^2");
    same("Factor[x^2/4 + y^2, GaussianIntegers -> True]", "((x - (2*I)*y)*(x + (2*I)*y))/4");
    same("Factor[x^4 + 4 y^4, GaussianIntegers -> True]",
         "(x - (1 + I)*y)*(x - (1 - I)*y)*(x + (1 - I)*y)*(x + (1 + I)*y)");
    same("Factor[-x^2 - y^2, GaussianIntegers -> True]", "-((x - I*y)*(x + I*y))");
    same("Factor[x^3 + x y^2, GaussianIntegers -> True]", "x*(x - I*y)*(x + I*y)");
    same("Factor[x^4 + 4 y^4, Extension -> {Sqrt[2], I}]",
         "(x - (1 + I)*y)*(x - (1 - I)*y)*(x + (1 - I)*y)*(x + (1 + I)*y)");
    same("Factor[x^2 + 2 Sqrt[2] I x y - 2 y^2, Extension -> Automatic]",
         "-((-I)*x + Sqrt[2]*y)^2");
    /* the reference absorbs Gaussian units into a factor in ways no rule
     * here reproduces: the product and the factor count are what is fixed */
    prop("I x^2 + I", "GaussianIntegers -> True", 2);
    prop("3 I x^2 + 3 I y^2", "GaussianIntegers -> True", 2);
    prop("Expand[(x + (1 + I) y)(I x + y)]", "GaussianIntegers -> True", 2);
    prop("Expand[((1 + I) x + 2 y)((2 - I) x + y)]", "GaussianIntegers -> True", 2);
}

/* ---- F1: Gaussian coefficients other than I itself ---- */
static void test_complex_coefficients(void) {
    same("Factor[x^2 + 4, Extension -> I]", "(-2*I + x)*(2*I + x)");
    same("Factor[x^2 + 2 Complex[0, 2] x - 4, Extension -> I]", "(2*I + x)^2");
    same("Factor[x^2 + 2 I x - 1, Extension -> I]", "(I + x)^2");
    same("Factor[x^2 + 2 x + 2, Extension -> I]", "((1 - I) + x)*((1 + I) + x)");
    same("Factor[x^4 + 4, Extension -> I]",
         "((-1 - I) + x)*((-1 + I) + x)*((1 - I) + x)*((1 + I) + x)");
    same("Factor[x^2 + 1/4, Extension -> I]", "((-I + 2*x)*(I + 2*x))/4");
}

/* ---- F1b: default Factor of a univariate Gaussian polynomial is over Q(i) ---- */
static void test_default_gaussian_univariate(void) {
    same("Factor[x^2 + 2 I x - 1]", "(I + x)^2");
    same("Factor[x^2 + 2 I x + 3]", "(-I + x)*(3*I + x)");
    same("Factor[(2 + 2 I) x^2 + 2 I + 2]", "(2 + 2*I)*(-I + x)*(I + x)");
    same("Factor[(1 - I) + (4 - 2 I) x + 4 x^2]", "(1 - I)*(1 + (1 + I)*x)*(1 + 2*x)");
    prop("I x^2 + I", "Extension -> None", 2);
    /* still multivariate as before */
    same("Factor[x^2 - 2 I x y - y^2]", "(x - I*y)^2");
    same("Factor[(2 + I) x + (1 - 2 I) y]", "(2 + I)*(x - I*y)");
    /* a radical coefficient keeps the default over Q */
    same("Factor[x^2 + 2 I Sqrt[2] x - 2]", "-2 + (2*I)*Sqrt[2]*x + x^2");
}

/* ---- F3: coefficients past a machine word ---- */
static void test_bigint_coefficients(void) {
    prop("x^2 - 2 (12345678901234567890123/98765432109876543210987)^2",
         "Extension -> Sqrt[2]", 2);
    prop("x^2 + 2 Sqrt[2] (12345678901234567890123/98765432109876543210987) x + "
         "2 (12345678901234567890123/98765432109876543210987)^2",
         "Extension -> Sqrt[2]", 2);
    prop("x^2 - 2 * 10^40 y^2", "Extension -> Sqrt[2]", 2);
}

/* ---- the univariate forms (also the reference's) ---- */
static void test_univariate(void) {
    same("Factor[x^2 - 2, Extension -> Sqrt[2]]", "-((Sqrt[2] - x)*(Sqrt[2] + x))");
    same("Factor[2 x^2 - 1, Extension -> Sqrt[2]]", "-1/2*((Sqrt[2] - 2*x)*(Sqrt[2] + 2*x))");
    same("Factor[x^4 + 1, Extension -> Sqrt[2]]",
         "-((-1 + Sqrt[2]*x - x^2)*(1 + Sqrt[2]*x + x^2))");
    same("Factor[x^2 + Sqrt[2] x + 1/2, Extension -> Sqrt[2]]", "(Sqrt[2] + 2*x)^2/4");
    same("Factor[x^2 - x - 1, Extension -> Sqrt[5]]",
         "-1/4*((1 + Sqrt[5] - 2*x)*(-1 + Sqrt[5] + 2*x))");
    same("Factor[x^2 - 12, Extension -> Sqrt[3]]", "-((2*Sqrt[3] - x)*(2*Sqrt[3] + x))");
    same("Factor[Sqrt[2] x^2 - 2 Sqrt[2], Extension -> Sqrt[2]]",
         "-(Sqrt[2]*(Sqrt[2] - x)*(Sqrt[2] + x))");
    same("Factor[x^4 - 4, Extension -> {Sqrt[2], I}]",
         "-((Sqrt[2] - x)*(Sqrt[2] - I*x)*(Sqrt[2] + I*x)*(Sqrt[2] + x))");
    same("Factor[x^4 + 1, Extension -> {Sqrt[2], I}]",
         "((Sqrt[2] - (1 + I)*x)*(Sqrt[2] - (1 - I)*x)*(Sqrt[2] + (1 - I)*x)*"
         "(Sqrt[2] + (1 + I)*x))/4");
    same("Factor[x^2 + x + 1, Extension -> Sqrt[-3]]",
         "((-I + Sqrt[3] - (2*I)*x)*(I + Sqrt[3] + (2*I)*x))/4");
    /* degree 40 over Q(2^(1/4)): was 17 s through Euclidean GCDs over the
     * field; the norm route is milliseconds */
    /* (verified in Q[t]/(t^4 - 2): Expand alone does not normalise every
     * product of these radicals to the power basis) */
    check_true("Module[{p = Expand[Product[x^5 + (1 + k 2^(1/4) - 2^(3/4)) x^2 + "
               "(k + Sqrt[2]) x + 2^(1/4) - k, {k, 1, 8}]], r}, "
               "r = Factor[p, Extension -> 2^(1/4)]; tfcount[r] == 8 && "
               "Expand[PolynomialRemainder[Expand[(r - p) /. {2^(3/4) -> t^3, Sqrt[2] -> t^2, "
               "2^(1/4) -> t}], t^4 - 2, t]] === 0]");
}

/* ---- FactorList, rational functions, threading ---- */
static void test_lists_and_threading(void) {
    same("FactorList[x^2 + y^2, GaussianIntegers -> True]",
         "{{1, 1}, {x - I*y, 1}, {x + I*y, 1}}");
    same("FactorList[x^2 - 2 y^2, Extension -> Sqrt[2]]",
         "{{-1, 1}, {-x + Sqrt[2]*y, 1}, {x + Sqrt[2]*y, 1}}");
    same("FactorList[-3 x^2 - 3 y^2, GaussianIntegers -> True]",
         "{{-3, 1}, {x - I*y, 1}, {x + I*y, 1}}");
    same("FactorList[-3 x^2 + 6 y^2, Extension -> Sqrt[2]]",
         "{{3, 1}, {-x + Sqrt[2]*y, 1}, {x + Sqrt[2]*y, 1}}");
    same("FactorList[x^4 - 2, Extension -> Sqrt[2]]",
         "{{-1, 1}, {Sqrt[2] - x^2, 1}, {Sqrt[2] + x^2, 1}}");
    same("FactorList[x^2 - 2 Sqrt[2] x + 2, Extension -> Sqrt[2]]", "{{1, 1}, {Sqrt[2] - x, 2}}");
    same("FactorList[Sqrt[2] x^2 - 2 Sqrt[2], Extension -> Sqrt[2]]",
         "{{-Sqrt[2], 1}, {Sqrt[2] - x, 1}, {Sqrt[2] + x, 1}}");
    same("Factor[(x^2 + y^2)/(x^2 - 2 y^2), Extension -> {I, Sqrt[2]}]",
         "-(((x - I*y)*(x + I*y))/((-x + Sqrt[2]*y)*(x + Sqrt[2]*y)))");
    same("Factor[x^2 + y^2 == 0, GaussianIntegers -> True]", "(x - I*y)*(x + I*y) == 0");
    same("Factor[{x^2 + y^2, x^2 + 1}, GaussianIntegers -> True]",
         "{(x - I*y)*(x + I*y), (-I + x)*(I + x)}");
    check("FactorSquareFree[x^2 + y^2, GaussianIntegers -> True]",
          "FactorSquareFree[x^2 + y^2, GaussianIntegers -> True]");
}

int main(void) {
    symtab_init();
    core_init();
    define_helpers();

    printf("Running factor extension tests...\n");
    test_multivariate_radicals();
    test_automatic();
    test_gaussian();
    test_complex_coefficients();
    test_default_gaussian_univariate();
    test_bigint_coefficients();
    test_univariate();
    test_lists_and_threading();

    printf("%d passed, %d failed\n", passes, failures);
    if (failures) {
        printf("%d factor extension test(s) FAILED\n", failures);
        return 1;
    }
    printf("All factor extension tests passed!\n");
    return 0;
}
