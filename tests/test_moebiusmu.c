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
        printf("MoebiusMu test failed: %s\n  expected: %s\n  got:      %s\n",
               input, expected, res_str);
        ASSERT(0);
    }
    free(res_str);
    expr_free(e);
    expr_free(res);
}

/* ---- Basic values: 0 / 1 / -1 from the prime factorisation ---- */
void test_moebiusmu_basic() {
    check("MoebiusMu[1]", "1");        /* unit, empty product */
    check("MoebiusMu[2]", "-1");       /* one prime */
    check("MoebiusMu[6]", "1");        /* 2 * 3, two primes */
    check("MoebiusMu[30]", "-1");      /* 2 * 3 * 5, three primes */
    check("MoebiusMu[4]", "0");        /* 2^2 */
    check("MoebiusMu[12]", "0");       /* 2^2 * 3 */
    check("MoebiusMu[11]", "-1");      /* prime */
    check("MoebiusMu[10]", "1");       /* 2 * 5 */
    check("MoebiusMu[32538]", "-1");   /* 2 * 3 * 11 * 17 * 29 (five primes) */
    check("MoebiusMu[1440]", "0");     /* 2^5 * 3^2 * 5 */
}

/* ---- Sign of n is ignored: mu(-n) == mu(n) ---- */
void test_moebiusmu_negative() {
    check("MoebiusMu[-6]", "1");
    check("MoebiusMu[-4]", "0");
    check("MoebiusMu[-11]", "-1");
    check("MoebiusMu[-30] == MoebiusMu[30]", "True");
}

/* ---- BigInt path: large argument ---- */
void test_moebiusmu_bignum() {
    check("MoebiusMu[10^50 + 1]", "-1");
    /* (2^61 - 1) is a Mersenne prime -> mu = -1 */
    check("MoebiusMu[2^61 - 1]", "-1");
    /* A perfect square of a large prime -> mu = 0 */
    check("MoebiusMu[(10^25 + 9)^2]", "0");
}

/* ---- Listable: threads element-wise over lists ---- */
void test_moebiusmu_listable() {
    check("MoebiusMu[{4, 10, 17, 20}]", "{0, 1, -1, 0}");
    check("MemberQ[Attributes[MoebiusMu], Listable]", "True");
}

/* ---- Gaussian integers, auto-detected from non-real Complex input ---- */
void test_moebiusmu_gaussian() {
    check("MoebiusMu[5 + 6 I]", "-1");   /* norm 61 prime -> Gaussian prime */
    check("MoebiusMu[1 + I]", "-1");     /* norm 2 prime -> Gaussian prime */
    check("MoebiusMu[1 + 3 I]", "1");    /* (1+I)(2+I): two distinct primes */
    check("MoebiusMu[3 + 4 I]", "0");    /* (2+I)^2: squared prime factor */
}

/* ---- Arity and non-integer arguments stay unevaluated ---- */
void test_moebiusmu_unevaluated() {
    check("MoebiusMu[]", "MoebiusMu[]");
    check("MoebiusMu[2, 3]", "MoebiusMu[2, 3]");
    check("MoebiusMu[x]", "MoebiusMu[x]");
    check("MoebiusMu[5/2]", "MoebiusMu[5/2]");
    check("MoebiusMu[2.5]", "MoebiusMu[2.5]");
    check("MoebiusMu[0]", "MoebiusMu[0]");
}

/* ---- Attributes ---- */
void test_moebiusmu_attributes() {
    check("Attributes[MoebiusMu]", "{Listable, Protected}");
    check("MemberQ[Attributes[MoebiusMu], Protected]", "True");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_moebiusmu_basic);
    TEST(test_moebiusmu_negative);
    TEST(test_moebiusmu_bignum);
    TEST(test_moebiusmu_listable);
    TEST(test_moebiusmu_gaussian);
    TEST(test_moebiusmu_unevaluated);
    TEST(test_moebiusmu_attributes);

    printf("All MoebiusMu tests passed!\n");
    return 0;
}
