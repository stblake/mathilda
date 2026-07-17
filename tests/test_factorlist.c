#include "expr.h"
#include "eval.h"
#include "core.h"
#include "symtab.h"
#include "test_utils.h"
#include "parse.h"
#include "print.h"
#include <string.h>
#include <stdlib.h>

/* Shared driver: parse, evaluate, compare the printed form. */
static void check(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    Expr* res = evaluate(e);
    char* res_str = expr_to_string(res);
    if (strcmp(res_str, expected) != 0) {
        printf("FactorList test failed: %s\n  expected: %s\n  got:      %s\n",
               input, expected, res_str);
        ASSERT(0);
    }
    free(res_str);
    expr_free(e);
    expr_free(res);
}

/* ---- Univariate polynomials: leading {1,1} numerical factor ---- */
void test_factorlist_univariate() {
    check("FactorList[x^2-1]", "{{1, 1}, {-1 + x, 1}, {1 + x, 1}}");
    check("FactorList[x^3-6x^2+11x-6]",
          "{{1, 1}, {-3 + x, 1}, {-2 + x, 1}, {-1 + x, 1}}");
    check("FactorList[x^2+1]", "{{1, 1}, {1 + x^2, 1}}");  /* irreducible over Q */
}

/* ---- The overall numerical factor ---- */
void test_factorlist_numeric_factor() {
    check("FactorList[2x^3+2x^2-2x-2]", "{{2, 1}, {-1 + x, 1}, {1 + x, 2}}");
    check("FactorList[3/4 x^2 - 3/4]", "{{3/4, 1}, {-1 + x, 1}, {1 + x, 1}}");
    check("FactorList[-x-1]", "{{-1, 1}, {1 + x, 1}}");   /* sign captured */
    /* Pure numbers: the whole expression is the numerical factor. */
    check("FactorList[6]", "{{6, 1}}");
    check("FactorList[7]", "{{7, 1}}");
}

/* ---- Multiplicities and multivariate factors ---- */
void test_factorlist_multivariate() {
    check("FactorList[x^3+2x^2 y+x y^2+x^2 z-y^2 z-x z^2-2y z^2-z^3]",
          "{{1, 1}, {x - z, 1}, {x + y + z, 2}}");
    /* Bare-symbol factors with multiplicity. */
    check("FactorList[x^2 y]", "{{1, 1}, {x, 2}, {y, 1}}");
}

/* ---- Rational functions: denominator factors carry negative exponents ---- */
void test_factorlist_rational() {
    check("FactorList[(x^3+2x^2)/(x^2-4y^2)-(x+2)/(x^2-4y^2)]",
          "{{1, 1}, {-1 + x, 1}, {1 + x, 1}, {2 + x, 1}, {x - 2 y, -1}, {x + 2 y, -1}}");
}

/* ---- Options are forwarded verbatim to Factor ---- */
void test_factorlist_options() {
    /* Extension -> I splits x^2 + 1 over the Gaussian integers. */
    check("FactorList[x^2+1, Extension->I]", "{{1, 1}, {-I + x, 1}, {I + x, 1}}");
    /* Extension -> Sqrt[2]: algebraic factors (Mathilda's Factor normal form). */
    check("FactorList[x^4-2, Extension->Sqrt[2]]",
          "{{1, 1}, {Sqrt[2] + x^2, 1}, {-Sqrt[2] + x^2, 1}}");
}

/* ---- Non-polynomial expressions (Factor treats kernels as generators) ---- */
void test_factorlist_nonpolynomial() {
    /* Sqrt[x] is a generator; Power[x,1/2] is itself an irreducible factor. */
    check("FactorList[x^2-Sqrt[x]]",
          "{{1, 1}, {Sqrt[x], 1}, {-1 + Sqrt[x], 1}, {1 + Sqrt[x] + x, 1}}");
}

/* ---- Round-trip invariant: Times @@ Power @@@ FactorList[p] == p ----
 * Order- and normalization-independent, so it validates the pairing even where
 * Mathilda's Factor orders factors differently from Mathematica. */
void test_factorlist_roundtrip() {
    check("Expand[(Times@@Power@@@FactorList[x^2-1]) - (x^2-1)]", "0");
    check("Expand[(Times@@Power@@@FactorList[2x^3+2x^2-2x-2]) - (2x^3+2x^2-2x-2)]", "0");
    check("Expand[(Times@@Power@@@FactorList["
          "x^8+11x^7+43x^6+59x^5-35x^4-151x^3-63x^2+81x+54]) - ("
          "x^8+11x^7+43x^6+59x^5-35x^4-151x^3-63x^2+81x+54)]", "0");
    check("Expand[(Times@@Power@@@FactorList[2x^3 y-2a^2 x y-3a^2 x^2+3a^4]) - "
          "(2x^3 y-2a^2 x y-3a^2 x^2+3a^4)]", "0");
    check("Expand[(Times@@Power@@@FactorList[x^4-2, Extension->Sqrt[2]]) - (x^4-2)]", "0");
}

/* ---- Attributes ---- */
void test_factorlist_attributes() {
    check("Attributes[FactorList]", "{Listable, Protected}");
    check("MemberQ[Attributes[FactorList], Listable]", "True");
    check("MemberQ[Attributes[FactorList], Protected]", "True");
}

/* ---- Argument errors: unevaluated + diagnostic ---- */
void test_factorlist_arity() {
    check("FactorList[]", "FactorList[]");                 /* argx */
    check("FactorList[1, 2, 3, 4]", "FactorList[1, 2, 3, 4]"); /* nonopt */
}

int main() {
    symtab_init();
    core_init();

    TEST(test_factorlist_univariate);
    TEST(test_factorlist_numeric_factor);
    TEST(test_factorlist_multivariate);
    TEST(test_factorlist_rational);
    TEST(test_factorlist_options);
    TEST(test_factorlist_nonpolynomial);
    TEST(test_factorlist_roundtrip);
    TEST(test_factorlist_attributes);
    TEST(test_factorlist_arity);

    printf("All FactorList tests passed!\n");
    return 0;
}
