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
        printf("PrimeNu test failed: %s\n  expected: %s\n  got:      %s\n",
               input, expected, res_str);
        ASSERT(0);
    }
    free(res_str);
    expr_free(e);
    expr_free(res);
}

/* ---- Basic values: nu(n) = number of DISTINCT primes in the factorisation - */
void test_primenu_basic() {
    check("PrimeNu[1]", "0");    /* unit, empty product */
    check("PrimeNu[7]", "1");    /* prime */
    check("PrimeNu[24]", "2");   /* 2^3 * 3 */
    check("PrimeNu[30]", "3");   /* 2 * 3 * 5 */
    check("PrimeNu[105]", "3");  /* 3 * 5 * 7 */
    check("PrimeNu[12]", "2");   /* 2^2 * 3 */
    check("PrimeNu[50]", "2");   /* 2 * 5^2 (matches Length[FactorInteger[50]]) */
    check("PrimeNu[32]", "1");   /* 2^5, a prime power */
    check("PrimeNu[49]", "1");   /* 7^2, a prime power */
}

/* ---- Sign of n is ignored: nu(-n) == nu(n); nu(+-1) == 0 ---- */
void test_primenu_negative() {
    check("PrimeNu[-1]", "0");
    check("PrimeNu[-30]", "3");
    check("PrimeNu[-30] == PrimeNu[30]", "True");
}

/* ---- BigInt path: large arguments ---- */
void test_primenu_bignum() {
    check("PrimeNu[50!]", "15");
    /* 49-digit value from the WL reference: 8 distinct prime factors */
    check("PrimeNu[2491230487120948712093481230948273409812734091238]", "8");
    /* (2^61 - 1) is a Mersenne prime -> nu = 1 */
    check("PrimeNu[2^61 - 1]", "1");
    /* A perfect square of a large prime (10^9 + 7) -> nu = 1 (single prime) */
    check("PrimeNu[(10^9 + 7)^2]", "1");
    /* Product of two distinct prime powers via the bigint path -> nu = 2 */
    check("PrimeNu[2^40 * 3^20]", "2");
}

/* ---- Listable: threads element-wise over lists ---- */
void test_primenu_listable() {
    check("PrimeNu[{4, 28, 180}]", "{1, 2, 3}");
    check("MemberQ[Attributes[PrimeNu], Listable]", "True");
}

/* ---- Gaussian integers: auto-detected and via GaussianIntegers -> True ---- */
void test_primenu_gaussian() {
    check("PrimeNu[3 + I]", "2");                       /* (1+I)(2-I), two primes */
    check("PrimeNu[1 + I]", "1");                       /* norm 2, single prime */
    check("PrimeNu[5 + 9 I]", "2");                     /* norm 106 = 2 * 53 */
    check("PrimeNu[105, GaussianIntegers -> True]", "4"); /* 3 | 5split | 7 */
    check("PrimeNu[30, GaussianIntegers -> True]", "4");  /* (1+I) | 3 | 5split */
    check("PrimeNu[30, GaussianIntegers -> False]", "3"); /* forced integer path */
}

/* ---- Relationships with PrimeOmega / MoebiusMu / LiouvilleLambda ---- */
void test_primenu_relations() {
    /* nu == Omega exactly when n is square-free */
    check("PrimeNu[210] == PrimeOmega[210]", "True");   /* 2*3*5*7, square-free */
    check("SquareFreeQ[42]", "True");
    check("MoebiusMu[42] == (-1)^PrimeNu[42]", "True");
    check("LiouvilleLambda[42] == (-1)^PrimeNu[42]", "True");
    /* Additive on coprime arguments: nu(9*40) == nu(9) + nu(40) */
    check("PrimeNu[9 40] == PrimeNu[9] + PrimeNu[40]", "True");
    /* Prime-power test companion */
    check("PrimeNu[32] == 1", "True");
}

/* ---- Arity / non-integer / zero arguments stay unevaluated ---- */
void test_primenu_unevaluated() {
    check("PrimeNu[]", "PrimeNu[]");
    check("PrimeNu[x]", "PrimeNu[x]");
    check("PrimeNu[5/2]", "PrimeNu[5/2]");
    check("PrimeNu[2.5]", "PrimeNu[2.5]");
    check("PrimeNu[0]", "PrimeNu[0]");
    /* Unknown option is left unevaluated */
    check("PrimeNu[30, Foo -> True]", "PrimeNu[30, Foo -> True]");
}

/* ---- Attributes ---- */
void test_primenu_attributes() {
    check("Attributes[PrimeNu]", "{Listable, Protected}");
    check("MemberQ[Attributes[PrimeNu], Protected]", "True");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_primenu_basic);
    TEST(test_primenu_negative);
    TEST(test_primenu_bignum);
    TEST(test_primenu_listable);
    TEST(test_primenu_gaussian);
    TEST(test_primenu_relations);
    TEST(test_primenu_unevaluated);
    TEST(test_primenu_attributes);

    printf("All PrimeNu tests passed!\n");
    return 0;
}
