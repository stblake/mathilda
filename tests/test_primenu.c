#include <stdio.h>
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

/* ---- BigInt path: large arguments ----
 *
 * The 49-digit WL reference value that used to stand here is FLAKY, and the
 * reference was not the problem -- the factorisation was. That number has a
 * 35-digit composite cofactor which ECM splits only sometimes; when it fails,
 * FactorInteger emits
 *     FactorInteger::nofac: ... is composite but no factor was found within
 *     the search bounds; it is returned unfactored with exponent 1
 * and PrimeNu counts that cofactor as ONE prime, giving 7 instead of 8.
 * Measured on a single unmodified build, three consecutive runs of the old
 * assertion returned 7, 8, 8 -- so it tested ECM's luck, not PrimeNu, and a
 * red here carried no information.
 *
 * Replaced by products of EXPLICIT small primes: above 2^63 so the bigint path
 * is still what runs, every factor inside trial division's reach so the result
 * is deterministic, and the expected count known by CONSTRUCTION rather than
 * from an oracle. */
void test_primenu_bignum() {
    check("PrimeNu[50!]", "15");
    /* 1009 * 2003 * 3001 * 4001 * 5003 * 6007 * 7001 > 2^63: seven distinct */
    check("PrimeNu[5105695083667128769810567]", "7");
    check("PrimeNu[5105695083667128769810567] == "
          "PrimeNu[1009 * 2003 * 3001 * 4001 * 5003 * 6007 * 7001]", "True");
    /* the same times an eighth prime (10007) */
    check("PrimeNu[51092690702256957599494343969]", "8");
    /* repeating a factor must not change nu (1009^2 * the other six) */
    check("PrimeNu[5151646339420132928738862103]", "7");
    /* 1009^3 * 2003^2 * 3001 > 2^63: nu counts 3, omega counts 6 */
    check("PrimeNu[12368054568910624561]", "3");
    check("PrimeOmega[12368054568910624561]", "6");
    /* (2^61 - 1) is a Mersenne prime -> nu = 1 */
    check("PrimeNu[2^61 - 1]", "1");
    /* A perfect square of a large prime (10^9 + 7) -> nu = 1 (single prime) */
    check("PrimeNu[(10^9 + 7)^2]", "1");
    /* Product of two distinct prime powers via the bigint path -> nu = 2 */
    check("PrimeNu[2^40 * 3^20]", "2");
    /* nu <= omega always, with equality exactly on the squarefree numbers */
    check("With[{n = 5105695083667128769810567}, PrimeNu[n] == PrimeOmega[n]]", "True");
    check("With[{n = 5151646339420132928738862103}, PrimeNu[n] < PrimeOmega[n]]", "True");
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

/* A number FactorInteger cannot finish must leave these heads unevaluated, not
 * produce a confident wrong integer.
 *
 * FactorInteger's Automatic method is bounded, and on a hard composite that
 * survives trial division, Pollard rho and ECM it returns the cofactor
 * UNFACTORED with exponent 1 (warning as it goes, but still answering). The
 * shared helper df_factor_mpz read that back as a prime factor, so for the
 * 82-digit semiprime below PrimeNu said 1 and PrimeOmega said 1 where the truth
 * is 2, and MoebiusMu said -1 where the truth is 1 -- while PrimeQ on the very
 * same number correctly said False. facint.c's own facint_factor_complete
 * already states the rule ("the caller must then DECLINE rather than trust a
 * possibly-composite prime"); the helper now enforces it, and every consumer
 * already returned NULL on a failed factorisation. */
void test_primenu_incomplete_factorisation_declines() {
    /* p = NextPrime[10^40], q = NextPrime[10^41] */
    const char* m = "1000000000000000000000000000000000000013190000000000000000000000000000000000013189";
    char buf[512], want[512];

    /* Mathilda knows it is composite ... */
    snprintf(buf, sizeof buf, "PrimeQ[%s]", m);
    check(buf, "False");

    /* ... so it must not claim it has one prime factor. */
    snprintf(buf,  sizeof buf,  "PrimeNu[%s]", m);
    snprintf(want, sizeof want, "PrimeNu[%s]", m);
    check(buf, want);

    snprintf(buf,  sizeof buf,  "PrimeOmega[%s]", m);
    snprintf(want, sizeof want, "PrimeOmega[%s]", m);
    check(buf, want);

    snprintf(buf,  sizeof buf,  "MoebiusMu[%s]", m);
    snprintf(want, sizeof want, "MoebiusMu[%s]", m);
    check(buf, want);

    snprintf(buf,  sizeof buf,  "LiouvilleLambda[%s]", m);
    snprintf(want, sizeof want, "LiouvilleLambda[%s]", m);
    check(buf, want);

    /* A big number that DOES factor must still answer: 2^64 - 1 has 7 distinct
     * primes, so the guard must not turn into a blanket refusal on large input. */
    check("PrimeNu[18446744073709551615]", "7");
    check("PrimeOmega[18446744073709551615]", "7");
    check("MoebiusMu[18446744073709551615]", "-1");
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
    TEST(test_primenu_incomplete_factorisation_declines);
    TEST(test_primenu_attributes);

    printf("All PrimeNu tests passed!\n");
    return 0;
}
