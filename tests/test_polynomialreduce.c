/* test_polynomialreduce.c -- PolynomialReduce
 *
 * PolynomialReduce[poly, {p1,...,pn}, {x1,...,xk}] -> {{a1,...,an}, b}.
 *
 * Two idioms:
 *   run_true(expr)      -- evaluate a Boolean and require "True".  Used for the
 *                          Wolfram-exact `===` checks (both sides canonicalise
 *                          identically) and for the domain-independent
 *                          reconstruction invariant q.p + b == poly.
 *   run_test(in, ff)    -- FullForm string compare (for decline / unevaluated).
 * Uses ASSERT_MSG (never libc assert; Release defines NDEBUG).
 */
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

static void run_true(const char* input) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* s = expr_to_string_fullform(res);
    if (strcmp(s, "True") != 0)
        printf("FAIL(true): %s\n  got: %s\n", input, s);
    ASSERT_MSG(strcmp(s, "True") == 0, "%s: expected True, got %s", input, s);
    free(s); expr_free(res); expr_free(e);
}

static void run_test(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* s = expr_to_string_fullform(res);
    if (strcmp(s, expected) != 0)
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, s);
    ASSERT_MSG(strcmp(s, expected) == 0, "%s: expected %s, got %s", input, expected, s);
    free(s); expr_free(res); expr_free(e);
}

/* ------------------------------------------------------------------ */

static void test_rationals_basic(void) {
    /* Canonical example (Lex, via FLINT fast path where present). */
    run_true("PolynomialReduce[x^3+y^3, {x^2-y^2-1, x+2y-7}, {x,y}] === "
             "{{x, 1+y^2}, 7-2y+7y^2-y^3}");
    /* Non-Groebner divisor set -> nonzero remainder. */
    run_true("PolynomialReduce[2x^3+y^3+3y, {x^2+y^2-1, x y-2}, {x,y}] === "
             "{{2x,-2y}, 2x-y+y^3}");
    /* Single-variable list form. */
    run_true("PolynomialReduce[x^3+2, {x^2-1}, x] === {{x}, 2+x}");
    run_true("PolynomialReduce[x^3+2, {x^2-1}, {x}] === {{x}, 2+x}");
    /* Zero remainder for exact division. */
    run_true("PolynomialReduce[x^2-1, {x-1}, {x}] === {{1+x}, 0}");
    /* Divisor order matters (quotients differ, both valid). */
    run_true("First[PolynomialReduce[x^2, {x, x+1}, {x}]] === {x, 0}");
}

static void test_reconstruction(void) {
    /* q.p + b == poly, the domain-independent correctness oracle. */
    run_true("With[{p=x^7+x^3 y^2+1, d={x^2-y, x y-1}}, "
             "Expand[First[#].d + Last[#] - p] == 0 &[PolynomialReduce[p,d,{x,y}]]]");
    run_true("With[{p=(x+y+z)^4, d={x^2-1, y^2-2, z^2-3}}, "
             "Expand[First[#].d + Last[#] - p] == 0 &[PolynomialReduce[p,d,{x,y,z}]]]");
    /* Remainder reduced -> reducing again changes nothing. */
    run_true("With[{p=x^5+y^5, d=GroebnerBasis[{x^2-y^3-5,y^2-x^3+7},{x,y}]}, "
             "Last[PolynomialReduce[p,d,{x,y}]] === "
             "Last[PolynomialReduce[Last[PolynomialReduce[p,d,{x,y}]],d,{x,y}]]]");
}

static void test_ideal_membership(void) {
    run_true("PolynomialReduce[x^6-14x^3-y^4+49, "
             "GroebnerBasis[{x^2-y^3-5,y^2-x^3+7},{x,y}], {x,y}][[2]] === 0");
    run_true("PolynomialReduce[x^5-y^5, "
             "GroebnerBasis[{x^2-y^3-5,y^2-x^3+7},{x,y}], {x,y}][[2]] === "
             "35+5y^2+7y^3");
}

static void test_variable_order(void) {
    run_true("PolynomialReduce[x^5+(x+y)^2, {x^2-y^3-5,y^3-x^3+7}, {x,y}] === "
             "{{1+5x+x^3+x y^3,0}, 5+25x+2x y+y^2+y^3+10x y^3+x y^6}");
    run_true("PolynomialReduce[x^5+(x+y)^2, {x^2-y^3-5,y^3-x^3+7}, {y,x}] === "
             "{{0,0}, x^2+x^5+2x y+y^2}");
}

static void test_monomial_order(void) {
    /* DegreeReverseLexicographic (gb_divmod path), Wolfram-exact. */
    run_true("PolynomialReduce[x y-3 x y^2+11y+x^3, "
             "GroebnerBasis[{3x^2+y-5x-1,2x+3x y+y^2},{x,y},"
             "MonomialOrder->DegreeReverseLexicographic], {x,y}, "
             "MonomialOrder->DegreeReverseLexicographic] === "
             "{{8/9-y,5/9+x/3,1/3}, 23/9+(8x)/3+(103y)/9-(74y^2)/9}");
    /* DegreeLexicographic reconstruction. */
    run_true("With[{p=x^4+y^4, d={x^3-2x y, x^2 y-2y^2+x}}, "
             "Expand[First[#].d + Last[#] - p] == 0 &["
             "PolynomialReduce[p,d,{x,y},MonomialOrder->DegreeLexicographic]]]");
}

static void test_rational_functions(void) {
    /* Parameters are coefficient-field elements: cofactors carry 1/a etc. */
    run_true("PolynomialReduce[a x+b x, {x}, {x}] === {{a+b}, 0}");
    /* Tag-variable / change of generators (the remainder is in the params). */
    run_true("Last[PolynomialReduce[(x+y)^2, "
             "GroebnerBasis[{a-(x^2+y^2), b-(x y)},{x,y}], {x,y}]] === a+2b");
    run_true("Last[PolynomialReduce[(x+y)^4, "
             "GroebnerBasis[{f1-(x^2+y^2), f2-(x y)},{x,y}], {x,y}]] === "
             "f1^2+4 f1 f2+4 f2^2");
    /* Reconstruction over Q(a): q.gb + r == poly. */
    run_true("With[{polys={a x^2+5y-1, 2x+x y-y^2}, poly=a^2 x-x y+y^2-3}, "
             "Together[First[#].GroebnerBasis[polys,{x,y},"
             "CoefficientDomain->RationalFunctions] + Last[#] - poly] == 0 &["
             "PolynomialReduce[poly, GroebnerBasis[polys,{x,y},"
             "CoefficientDomain->RationalFunctions], {x,y}]]]");
    /* The remainder over Q(a) equals Wolfram's exactly (compared as rational
     * functions -- SameQ would be brittle to grouping of the fraction form). */
    run_true("With[{r=Last[PolynomialReduce[a^2 x-x y+y^2-3, "
             "GroebnerBasis[{a x^2+5y-1, 2x+x y-y^2},{x,y},"
             "CoefficientDomain->RationalFunctions], {x,y}]], "
             "e=(2-6a+a^2)/(2a) - (9(2+a^2)y)/(4a) + "
             "((-10+4a-5a^2+2a^3)y^2)/(4a) + (1/4)(-2-a^2)y^3}, "
             "Together[r-e] == 0]");
}

static void test_modulus(void) {
    run_true("PolynomialReduce[x^3, {x^2+1}, {x}, Modulus->7] === {{x}, 6x}");
    /* Reconstruction mod p. */
    run_true("With[{p=x^4+2x^2 y+y^3, d={x^2-y, x y-2}}, "
             "PolynomialMod[Expand[First[#].d + Last[#] - p], 7] == 0 &["
             "PolynomialReduce[p, d, {x,y}, Modulus->7]]]");
    /* Vars omitted -> Variables[...]; reconstruction against a GF(p) basis. */
    run_true("With[{p=z^9-x^2 y^3-3x y^2 z+11y z^2+x^2 z^2-5, "
             "g=GroebnerBasis[{3x^2+y z-5x-1,2x+3x y+y^2,x-3y+x z-2z^2},"
             "{x,y,z},Modulus->7]}, "
             "PolynomialMod[Expand[First[#].g + Last[#] - p], 7] == 0 &["
             "PolynomialReduce[p, g, Modulus->7]]]");
}

static void test_decline(void) {
    /* Deferred coefficient domains and nonzero tolerance stay unevaluated. */
    run_test("PolynomialReduce[x^2, {x}, {x}, CoefficientDomain->Integers]",
             "PolynomialReduce[Power[x, 2], List[x], List[x], "
             "Rule[CoefficientDomain, Integers]]");
    run_test("PolynomialReduce[x^2, {x}, {x}, CoefficientDomain->InexactNumbers]",
             "PolynomialReduce[Power[x, 2], List[x], List[x], "
             "Rule[CoefficientDomain, InexactNumbers]]");
    /* Non-list divisor argument -> unevaluated. */
    run_test("PolynomialReduce[x^2, x, {x}]",
             "PolynomialReduce[Power[x, 2], x, List[x]]");
    /* Modulus with a free parameter -> unsupported, unevaluated. */
    run_test("PolynomialReduce[a x^2, {x}, {x}, Modulus->7]",
             "PolynomialReduce[Times[a, Power[x, 2]], List[x], List[x], "
             "Rule[Modulus, 7]]");
}

static void test_options_and_attributes(void) {
    /* Default options registered. */
    run_true("MemberQ[Options[PolynomialReduce], "
             "(CoefficientDomain -> RationalFunctions)]");
    run_true("MemberQ[Options[PolynomialReduce], (MonomialOrder -> Lexicographic)]");
    /* Protected. */
    run_true("MemberQ[Attributes[PolynomialReduce], Protected]");
}

static void test_groebnerbasis_order_fix(void) {
    /* DegreeLexicographic is now honoured (differs from the Lex basis) rather
     * than silently downgraded. */
    run_true("GroebnerBasis[{x^3-2x y, x^2 y-2y^2+x},{x,y},"
             "MonomialOrder->DegreeLexicographic] =!= "
             "GroebnerBasis[{x^3-2x y, x^2 y-2y^2+x},{x,y}]");
}

int main(void) {
    setbuf(stdout, NULL);
    printf("Starting polynomialreduce_tests\n");
    symtab_init();
    core_init();

    TEST(test_rationals_basic);
    TEST(test_reconstruction);
    TEST(test_ideal_membership);
    TEST(test_variable_order);
    TEST(test_monomial_order);
    TEST(test_rational_functions);
    TEST(test_modulus);
    TEST(test_decline);
    TEST(test_options_and_attributes);
    TEST(test_groebnerbasis_order_fix);

    printf("All PolynomialReduce tests passed!\n");
    return 0;
}
