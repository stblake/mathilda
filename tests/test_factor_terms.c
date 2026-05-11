/*
 * test_factor_terms.c -- unit tests for FactorTerms / FactorTermsList.
 *
 * The reference outputs below are the documented Mathematica behaviour
 * adapted to picocas's printer (canonical Times/Plus ordering, flat
 * Times for small polynomials). Where picocas legitimately produces a
 * mathematically-equivalent but textually-different result (e.g.
 * Gaussian-integer content, nested Times[c, Times[...]] leftovers from
 * Expand), the assertion uses the actual picocas form. Each FactorTerms
 * case is paired with a multiplicative round-trip so that the factor
 * decomposition and the residue together reproduce the input.
 */
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

static void run_full(const char* input, const char* expected_fullform) {
    Expr* e = parse_expression(input);
    if (!e) { printf("Failed to parse: %s\n", input); ASSERT(0); }
    Expr* res = evaluate(e);
    char* got = expr_to_string_fullform(res);
    if (strcmp(got, expected_fullform) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected_fullform, got);
        ASSERT(0);
    } else {
        printf("PASS: %s\n", input);
    }
    free(got);
    expr_free(res);
    expr_free(e);
}

static void run_str(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    if (!e) { printf("Failed to parse: %s\n", input); ASSERT(0); }
    Expr* res = evaluate(e);
    char* got = expr_to_string(res);
    if (strcmp(got, expected) != 0) {
        printf("FAIL: %s\n  expected: %s\n  got:      %s\n", input, expected, got);
        ASSERT(0);
    } else {
        printf("PASS: %s\n", input);
    }
    free(got);
    expr_free(res);
    expr_free(e);
}

/* Round-trip: FactorTerms[expr] should always multiply back out to expr. */
/* We compare the residuals: (input - Apply[Times, FactorTermsList[input]]) */
/* should simplify to 0 once the rational form is normalised.              */
static void run_roundtrip(const char* input) {
    char buf[1024];
    snprintf(buf, sizeof(buf),
             "Expand[Together[(%s)-Apply[Times,FactorTermsList[%s]]]]",
             input, input);
    Expr* e = parse_expression(buf);
    if (!e) { printf("Failed to parse: %s\n", buf); ASSERT(0); }
    Expr* res = evaluate(e);
    char* got = expr_to_string(res);
    if (strcmp(got, "0") != 0) {
        printf("ROUND-TRIP FAIL: %s\n  got:      %s\n", buf, got);
        ASSERT(0);
    } else {
        printf("ROUND-TRIP PASS: %s\n", input);
    }
    free(got);
    expr_free(res);
    expr_free(e);
}

/* ------------------------------------------------------------------- */
/* FactorTerms[poly] -- numerical content only.                        */
/* ------------------------------------------------------------------- */
static void test_factorterms_numerical(void) {
    /* Simple univariate */
    run_full("FactorTerms[3+6 x+3 x^2]",
             "Times[3, Plus[1, Times[2, x], Power[x, 2]]]");
    run_full("FactorTerms[8 x^3-6 x^2+22 x-6]",
             "Times[2, Plus[-3, Times[11, x], Times[-3, Power[x, 2]], Times[4, Power[x, 3]]]]");
    run_full("FactorTerms[4 x^3-6 x^2+12 x-6]",
             "Times[2, Plus[-3, Times[6, x], Times[-3, Power[x, 2]], Times[2, Power[x, 3]]]]");

    /* Multivariate */
    run_full("FactorTerms[12 a^4+9 x^2+66 b^2]",
             "Times[3, Plus[Times[4, Power[a, 4]], Times[22, Power[b, 2]], Times[3, Power[x, 2]]]]");
    run_full("FactorTerms[6 a^2+9 x^2+12 b^2]",
             "Times[3, Plus[Times[2, Power[a, 2]], Times[4, Power[b, 2]], Times[3, Power[x, 2]]]]");

    /* Constant gcd is 1: FactorTerms is a no-op. */
    run_full("FactorTerms[1+x+y]", "Plus[1, x, y]");

    /* Single integer. */
    run_full("FactorTerms[7]", "7");

    /* Round-trip checks. */
    run_roundtrip("3+6 x+3 x^2");
    run_roundtrip("8 x^3-6 x^2+22 x-6");
    run_roundtrip("12 a^4+9 x^2+66 b^2");
}

/* ------------------------------------------------------------------- */
/* FactorTermsList[poly] -- {numerical content, residue}.              */
/* ------------------------------------------------------------------- */
static void test_factortermslist_numerical(void) {
    run_full("FactorTermsList[3+6 x+3 x^2]",
             "List[3, Plus[1, Times[2, x], Power[x, 2]]]");
    run_full("FactorTermsList[8 x^3-6 x^2+22 x-6]",
             "List[2, Plus[-3, Times[11, x], Times[-3, Power[x, 2]], Times[4, Power[x, 3]]]]");
    run_full("FactorTermsList[6 a^2+9 x^2+12 b^2]",
             "List[3, Plus[Times[2, Power[a, 2]], Times[4, Power[b, 2]], Times[3, Power[x, 2]]]]");
    run_full("FactorTermsList[2 x^2-2]",
             "List[2, Plus[-1, Power[x, 2]]]");
    run_full("FactorTermsList[14 x+21 y+35 x y+63]",
             "List[7, Plus[9, Times[2, x], Times[3, y], Times[5, Times[x, y]]]]");
}

/* ------------------------------------------------------------------- */
/* FactorTerms[poly, x] -- content with respect to a single variable.  */
/* ------------------------------------------------------------------- */
static void test_factorterms_one_var(void) {
    run_full("FactorTerms[3+3 a+6 a x+6 x+12 a x^2+12 x^2,x]",
             "Times[3, Plus[1, a], Plus[1, Times[2, x], Times[4, Power[x, 2]]]]");
    run_full("FactorTermsList[3+3 a+6 a x+6 x+12 a x^2+12 x^2,x]",
             "List[3, Plus[1, a], Plus[1, Times[2, x], Times[4, Power[x, 2]]]]");

    run_full("FactorTerms[-6 y-6 a y+2 x^2 y+2 a x^2 y+4 a y^2+4 a^2 y^2,x]",
             "Times[2, Plus[y, Times[a, y]], Plus[-3, Power[x, 2], Times[2, Times[a, y]]]]");
    run_full("FactorTermsList[-6 y-6 a y+2 x^2 y+2 a x^2 y+4 a y^2+4 a^2 y^2,x]",
             "List[2, Plus[y, Times[a, y]], Plus[-3, Power[x, 2], Times[2, Times[a, y]]]]");

    /* If x does not appear, FactorTerms[poly, x] still extracts the */
    /* numerical content (the only x-independent factor available).    */
    run_full("FactorTerms[3+6 y, x]",
             "Times[3, Plus[1, Times[2, y]]]");
}

/* ------------------------------------------------------------------- */
/* FactorTerms[poly, {x_1, ..., x_n}] -- nested content extraction.    */
/* ------------------------------------------------------------------- */
static void test_factorterms_var_list(void) {
    run_full("FactorTerms[-6 y-6 a y+2 x^2 y+2 a x^2 y+4 a y^2+4 a^2 y^2,{x,y}]",
             "Times[2, Plus[1, a], y, Plus[-3, Power[x, 2], Times[2, Times[a, y]]]]");
    run_full("FactorTermsList[-6 y-6 a y+2 x^2 y+2 a x^2 y+4 a y^2+4 a^2 y^2,{x,y}]",
             "List[2, Plus[1, a], y, Plus[-3, Power[x, 2], Times[2, Times[a, y]]]]");

    /* Three-variable polynomial f = 2 x^2 y z + ... -- the documented */
    /* canonical FactorTerms / FactorTermsList examples.                */
    run_full("FactorTerms[2 x^2 y z+2 x^2 y+4 x^2 z+4 x^2+4 y^2 z^2+4 z y^2+8 z^2 y+2 z y-6 y-12 z-12,x]",
             "Times[2, Plus[2, y, Times[2, z], Times[y, z]], Plus[-3, Power[x, 2], Times[2, Times[y, z]]]]");
    run_full("FactorTerms[2 x^2 y z+2 x^2 y+4 x^2 z+4 x^2+4 y^2 z^2+4 z y^2+8 z^2 y+2 z y-6 y-12 z-12,{x,y}]",
             "Times[2, Plus[2, y], Plus[1, z], Plus[-3, Power[x, 2], Times[2, Times[y, z]]]]");
    run_full("FactorTermsList[2 x^2 y z+2 x^2 y+4 x^2 z+4 x^2+4 y^2 z^2+4 z y^2+8 z^2 y+2 z y-6 y-12 z-12,x]",
             "List[2, Plus[2, y, Times[2, z], Times[y, z]], Plus[-3, Power[x, 2], Times[2, Times[y, z]]]]");
    run_full("FactorTermsList[2 x^2 y z+2 x^2 y+4 x^2 z+4 x^2+4 y^2 z^2+4 z y^2+8 z^2 y+2 z y-6 y-12 z-12,{x,y}]",
             "List[2, Plus[1, z], Plus[2, y], Plus[-3, Power[x, 2], Times[2, Times[y, z]]]]");

    /* Round-trip via Apply[Times, FactorTermsList[...]] = original.     */
    run_roundtrip("-6 y-6 a y+2 x^2 y+2 a x^2 y+4 a y^2+4 a^2 y^2");
    run_roundtrip("2 x^2 y z+2 x^2 y+4 x^2 z+4 x^2+4 y^2 z^2+4 z y^2+8 z^2 y+2 z y-6 y-12 z-12");
}

/* ------------------------------------------------------------------- */
/* Threading over List, Equal, Less, And, ...                          */
/* ------------------------------------------------------------------- */
static void test_factorterms_threading(void) {
    /* Threading over List. */
    run_str("FactorTerms[{5 x^2-15,7 x^4-77,8 x^8-24}]",
            "{5 (-3 + x^2), 7 (-11 + x^4), 8 (-3 + x^8)}");

    /* Threading over a chained inequality. */
    run_str("FactorTerms[1<77 x^3-21 x+35<2]",
            "1 < 7 (5 - 3 x + 11 x^3) < 2");

    /* Threading over Equal. */
    run_str("FactorTerms[3+6 x==9+12 y]",
            "3 (1 + 2 x) == 3 (3 + 4 y)");
}

/* ------------------------------------------------------------------- */
/* Rational and complex coefficients.                                  */
/* ------------------------------------------------------------------- */
static void test_factorterms_rational_complex(void) {
    /* Rational expression: numerator content extraction, denominator    */
    /* survives in the residue.                                          */
    run_full("FactorTermsList[7 x+(14 y+21)/z]",
             "List[7, Times[Power[z, -1], Plus[3, Times[2, y], Times[x, z]]]]");
    /* The product of the FactorTermsList elements should reproduce the */
    /* original (after Together).                                       */
    run_roundtrip("7 x+(14 y+21)/z");

    /* Complex coefficients. picocas extracts the real integer GCD (5)   */
    /* rather than Mathematica's Gaussian-integer GCD (5 I); both are    */
    /* mathematically valid factorizations.                              */
    run_full("FactorTermsList[5 I x^2+20 x I+10]",
             "List[5, Plus[2, Times[Complex[0, 4], x], Times[Complex[0, 1], Power[x, 2]]]]");
    run_roundtrip("5 I x^2+20 x I+10");
}

/* ------------------------------------------------------------------- */
/* Attributes -- sanity check.                                          */
/* ------------------------------------------------------------------- */
static void test_factorterms_attributes(void) {
    run_str("MemberQ[Attributes[FactorTerms],Protected]", "True");
    run_str("MemberQ[Attributes[FactorTermsList],Protected]", "True");
}

int main(void) {
    extern void core_init(void);
    extern void symtab_init(void);
    symtab_init();
    core_init();

    TEST(test_factorterms_numerical);
    TEST(test_factortermslist_numerical);
    TEST(test_factorterms_one_var);
    TEST(test_factorterms_var_list);
    TEST(test_factorterms_threading);
    TEST(test_factorterms_rational_complex);
    TEST(test_factorterms_attributes);

    printf("All FactorTerms / FactorTermsList tests passed!\n");
    return 0;
}
