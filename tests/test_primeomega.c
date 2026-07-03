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
        printf("PrimeOmega test failed: %s\n  expected: %s\n  got:      %s\n",
               input, expected, res_str);
        ASSERT(0);
    }
    free(res_str);
    expr_free(e);
    expr_free(res);
}

/* ---- Basic values: Omega(n) = sum of exponents in the factorisation ---- */
void test_primeomega_basic() {
    check("PrimeOmega[1]", "0");    /* unit, empty product */
    check("PrimeOmega[7]", "1");    /* prime */
    check("PrimeOmega[30]", "3");   /* 2 * 3 * 5 */
    check("PrimeOmega[12]", "3");   /* 2^2 * 3 */
    check("PrimeOmega[28]", "3");   /* 2^2 * 7 */
    check("PrimeOmega[8]", "3");    /* 2^3 */
    check("PrimeOmega[2^10]", "10");/* 2^10 */
    check("PrimeOmega[6]", "2");    /* 2 * 3 */
}

/* ---- Sign of n is ignored: Omega(-n) == Omega(n); Omega(+-1) == 0 ---- */
void test_primeomega_negative() {
    check("PrimeOmega[-1]", "0");
    check("PrimeOmega[-12]", "3");
    check("PrimeOmega[-30] == PrimeOmega[30]", "True");
}

/* ---- BigInt path: large arguments ---- */
void test_primeomega_bignum() {
    check("PrimeOmega[30!]", "59");
    /* (2^61 - 1) is a Mersenne prime -> Omega = 1 */
    check("PrimeOmega[2^61 - 1]", "1");
    /* A perfect square of a large prime (10^9 + 7) -> Omega = 2 */
    check("PrimeOmega[(10^9 + 7)^2]", "2");
    /* Prime-power products via the bigint path */
    check("PrimeOmega[2^40 * 3^20]", "60");
    check("PrimeOmega[100!]", "239");
}

/* ---- Listable: threads element-wise over lists ---- */
void test_primeomega_listable() {
    check("PrimeOmega[{4, 12, 24}]", "{2, 3, 4}");
    check("MemberQ[Attributes[PrimeOmega], Listable]", "True");
}

/* ---- Gaussian integers: auto-detected and via GaussianIntegers -> True ---- */
void test_primeomega_gaussian() {
    check("PrimeOmega[5 + 9 I]", "2");                    /* norm 106 = 2 * 53 */
    check("PrimeOmega[1 + I]", "1");                      /* norm 2 prime */
    check("PrimeOmega[12, GaussianIntegers -> True]", "5"); /* (1+I)^4 * 3 */
    check("PrimeOmega[28, GaussianIntegers -> True]", "5"); /* (1+I)^4 * 7 */
    check("PrimeOmega[30, GaussianIntegers -> False]", "3");/* forced integer path */
}

/* ---- Arity / non-integer / zero arguments stay unevaluated ---- */
void test_primeomega_unevaluated() {
    check("PrimeOmega[]", "PrimeOmega[]");
    check("PrimeOmega[x]", "PrimeOmega[x]");
    check("PrimeOmega[5/2]", "PrimeOmega[5/2]");
    check("PrimeOmega[2.5]", "PrimeOmega[2.5]");
    check("PrimeOmega[0]", "PrimeOmega[0]");
}

/* ---- Attributes ---- */
void test_primeomega_attributes() {
    check("Attributes[PrimeOmega]", "{Listable, Protected}");
    check("MemberQ[Attributes[PrimeOmega], Protected]", "True");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_primeomega_basic);
    TEST(test_primeomega_negative);
    TEST(test_primeomega_bignum);
    TEST(test_primeomega_listable);
    TEST(test_primeomega_gaussian);
    TEST(test_primeomega_unevaluated);
    TEST(test_primeomega_attributes);

    printf("All PrimeOmega tests passed!\n");
    return 0;
}
