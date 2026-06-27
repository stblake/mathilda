/* Tests for Prime[n] (the nth prime) and PrimePi[x] (prime counting),
 * both implemented in src/numbertheory/prime.c. */

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
        printf("Prime test failed: %s\n  expected: %s\n  got:      %s\n",
               input, expected, res_str);
        ASSERT(0);
    }
    free(res_str);
    expr_free(e);
    expr_free(res);
}

/* ---- Small primes read straight out of the sieve table ---- */
void test_prime_small() {
    check("Prime[1]", "2");
    check("Prime[2]", "3");
    check("Prime[3]", "5");
    check("Prime[4]", "7");
    check("Prime[5]", "11");
    check("Prime[6]", "13");
    check("Prime[10]", "29");
    check("Prime[25]", "97");
    check("Prime[100]", "541");
    check("Prime[1000]", "7919");
    check("Prime[10000]", "104729");
    check("Prime[100000]", "1299709");
}

/* ---- Boundary between the sieve table and the estimate/walk path ---- */
void test_prime_table_boundary() {
    check("Prime[78498]", "999983");      /* largest prime < 10^6 */
    check("Prime[78499]", "1000003");     /* first prime above the table */
    check("Prime[78500]", "1000033");
}

/* ---- Large n via PrimePi inversion + walk ---- */
void test_prime_large() {
    check("Prime[1000000]", "15485863");        /* 10^6 */
    check("Prime[10000000]", "179424673");      /* 10^7 */
    check("Prime[100000000]", "2038074743");    /* 10^8 */
    check("Prime[1000000000]", "22801763489");  /* 10^9 */
    check("Prime[10^10]", "252097800623");      /* 10^10 */
}

/* ---- Listable threads over lists ---- */
void test_prime_listable() {
    check("Prime[{1, 3, 4, 10}]", "{2, 5, 7, 29}");
    check("Prime[{100, 1000, 10000}]", "{541, 7919, 104729}");
    check("Attributes[Prime]", "{Listable, Protected}");
    check("MemberQ[Attributes[Prime], Listable]", "True");
    check("MemberQ[Attributes[Prime], Protected]", "True");
}

/* ---- PrimePi exact values ---- */
void test_primepi_values() {
    check("PrimePi[1]", "0");
    check("PrimePi[2]", "1");
    check("PrimePi[10]", "4");
    check("PrimePi[100]", "25");
    check("PrimePi[1000]", "168");
    check("PrimePi[10000]", "1229");
    check("PrimePi[1000000]", "78498");
    check("PrimePi[1000000000]", "50847534");     /* 10^9  */
    check("PrimePi[10^10]", "455052511");          /* 10^10 */
}

/* ---- Prime and PrimePi are inverse ---- */
void test_prime_roundtrip() {
    check("PrimePi[Prime[1]]", "1");
    check("PrimePi[Prime[50]]", "50");
    check("PrimePi[Prime[12345]]", "12345");
    check("PrimePi[Prime[1000000]]", "1000000");
    check("Prime[PrimePi[7919]]", "7919");
    check("Prime[PrimePi[15485863]]", "15485863");
    /* p_n - p_{n-1} is the prime gap; Prime[1000]=7919, Prime[999]=7907. */
    check("Prime[1000] - Prime[999]", "12");
}

/* ---- Invalid arguments: Prime::intpp / Prime::argx, left unevaluated ---- */
void test_prime_errors() {
    check("Prime[6.5]", "Prime[6.5]");     /* non-integer real */
    check("Prime[-19]", "Prime[-19]");     /* negative          */
    check("Prime[0]", "Prime[0]");         /* zero              */
    check("Prime[3/2]", "Prime[3/2]");     /* rational          */
    check("Prime[]", "Prime[]");           /* argx              */
    check("Prime[3, 4]", "Prime[3, 4]");   /* argx              */
    check("Prime[x]", "Prime[x]");         /* symbolic: deferred */
}

/* pi(10^7)=664579, pi(10^8)=5761455. */
void test_primepi_methods_values() {
    check("PrimePi[10^7, Method->\"Sieve\"]", "664579");
    check("PrimePi[10^7, Method->\"Legendre\"]", "664579");
    check("PrimePi[10^7, Method->\"Meissel\"]", "664579");
    check("PrimePi[10^7, Method->\"Lehmer\"]", "664579");
    check("PrimePi[10^7, Method->\"LMO\"]", "664579");
    check("PrimePi[10^7, Method->\"DelegliseRivat\"]", "664579");
    check("PrimePi[10^7, Method->\"LucyHedgehog\"]", "664579");
    check("PrimePi[10^8, Method->\"Meissel\"]", "5761455");
    check("PrimePi[10^8, Method->\"Lehmer\"]", "5761455");
    check("PrimePi[10^8, Method->\"LMO\"]", "5761455");
}

/* All methods produce identical results (cross-validation). */
void test_primepi_methods_agree() {
    check("SameQ @@ {PrimePi[10^7, Method->\"Sieve\"], "
          "PrimePi[10^7, Method->\"Legendre\"], "
          "PrimePi[10^7, Method->\"Meissel\"], "
          "PrimePi[10^7, Method->\"Lehmer\"], "
          "PrimePi[10^7, Method->\"LMO\"], "
          "PrimePi[10^7, Method->\"DelegliseRivat\"], "
          "PrimePi[10^7, Method->\"LucyHedgehog\"], PrimePi[10^7]}", "True");
    /* non-power-of-ten point */
    check("SameQ @@ {PrimePi[12345678, Method->\"Sieve\"], "
          "PrimePi[12345678, Method->\"Meissel\"], "
          "PrimePi[12345678, Method->\"Lehmer\"], "
          "PrimePi[12345678, Method->\"LMO\"], PrimePi[12345678]}", "True");
}

/* Stress: LMO matches the Lucy oracle on every integer in a multi-segment
 * window (off-by-one detector for the special-leaf sieve). */
void test_primepi_lmo_stress() {
    check("Count[Table[PrimePi[10^9 + k, Method->\"LMO\"] - PrimePi[10^9 + k], "
          "{k, 0, 60}], Except[0]]", "0");
    check("Count[Table[PrimePi[800000000 + 173 k, Method->\"Meissel\"] - "
          "PrimePi[800000000 + 173 k, Method->\"LMO\"], {k, 0, 40}], Except[0]]", "0");
}

/* Invalid / unknown Method leaves the call unevaluated. */
void test_primepi_method_errors() {
    check("Head[PrimePi[100, Method->\"Bogus\"]]", "PrimePi");
    check("PrimePi[1000, Method->Automatic]", "168");
    check("PrimePi[1000, Method->\"LMO\"]", "168");   /* small x via table */
}

/* Default options are registered. */
void test_primepi_options() {
    check("Options[PrimePi]", "{Method -> Automatic}");
}

/* ---- Out-of-range but valid inputs stay unevaluated (no crash) ---- */
void test_prime_out_of_range() {
    check("Head[Prime[10^30]]", "Prime");     /* n far beyond computable range */
    check("Head[PrimePi[10^20]]", "PrimePi"); /* x beyond PI_COUNT_MAX         */
    check("PrimePi[-5]", "0");
    check("PrimePi[1.5]", "0");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_prime_small);
    TEST(test_prime_table_boundary);
    TEST(test_prime_large);
    TEST(test_prime_listable);
    TEST(test_primepi_values);
    TEST(test_primepi_methods_values);
    TEST(test_primepi_methods_agree);
    TEST(test_primepi_lmo_stress);
    TEST(test_primepi_method_errors);
    TEST(test_primepi_options);
    TEST(test_prime_roundtrip);
    TEST(test_prime_errors);
    TEST(test_prime_out_of_range);

    printf("All Prime tests passed!\n");
    return 0;
}
