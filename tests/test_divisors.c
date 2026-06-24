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
        printf("Divisors test failed: %s\n  expected: %s\n  got:      %s\n",
               input, expected, res_str);
        ASSERT(0);
    }
    free(res_str);
    expr_free(e);
    expr_free(res);
}

/* ---- Ordinary integer divisors ---- */
void test_divisors_basic() {
    check("Divisors[1729]", "{1, 7, 13, 19, 91, 133, 247, 1729}");
    check("Divisors[6]", "{1, 2, 3, 6}");
    check("Divisors[1]", "{1}");
    check("Divisors[2]", "{1, 2}");      /* prime */
    check("Divisors[12]", "{1, 2, 3, 4, 6, 12}");
    check("Divisors[100]", "{1, 2, 4, 5, 10, 20, 25, 50, 100}");
    check("Divisors[36]", "{1, 2, 3, 4, 6, 9, 12, 18, 36}");
    check("Divisors[605]", "{1, 5, 11, 55, 121, 605}");
}

/* ---- Sign is ignored: Divisors[-n] == Divisors[n] ---- */
void test_divisors_negative() {
    check("Divisors[-12]", "{1, 2, 3, 4, 6, 12}");
    check("Divisors[-1]", "{1}");
    check("Divisors[-7]", "{1, 7}");
    check("Divisors[-12] == Divisors[12]", "True");
}

/* ---- Listable: threads element-wise over lists ---- */
void test_divisors_listable() {
    check("Divisors[{605, 871, 824}]",
          "{{1, 5, 11, 55, 121, 605}, {1, 13, 67, 871}, "
          "{1, 2, 4, 8, 103, 206, 412, 824}}");
    check("Divisors[{1, 6}]", "{{1}, {1, 2, 3, 6}}");
    check("MemberQ[Attributes[Divisors], Listable]", "True");
}

/* ---- BigInt path: large value, small divisor count ---- */
void test_divisors_bignum() {
    check("Divisors[2^60] // Length", "61");
    check("Divisors[2^60] // Last", "1152921504606846976");
    check("Divisors[12!] // Length", "792");   /* 2^10 3^5 5^2 7 11 */
    check("Last[Divisors[10^20]]", "100000000000000000000");
    /* sum of divisors cross-check: sigma(28) = 56 (perfect number) */
    check("Total[Divisors[28]]", "56");
}

/* ---- Gaussian divisors, auto-detected from Complex input ---- */
void test_divisors_gaussian_auto() {
    check("Divisors[6 + 4 I]", "{1, 1 + I, 1 + 5*I, 2, 3 + 2*I, 6 + 4*I}");
    /* 1 + I has norm 2 (a single Gaussian prime). */
    check("Divisors[1 + I]", "{1, 1 + I}");
    /* 2 + 3 I has prime norm 13, so it is a Gaussian prime. */
    check("Divisors[2 + 3 I]", "{1, 2 + 3*I}");
}

/* ---- GaussianIntegers -> True forces Z[i] even for plain integers ---- */
void test_divisors_gaussian_option() {
    check("Divisors[2, GaussianIntegers -> True]", "{1, 1 + I, 2}");
    /* 3 is a Gaussian prime (3 mod 4 == 3). */
    check("Divisors[3, GaussianIntegers -> True]", "{1, 3}");
    check("Divisors[3]", "{1, 3}");
    /* 5 = (2+I)(2-I) splits in Z[i]. */
    check("Divisors[5, GaussianIntegers -> True]", "{1, 1 + 2*I, 2 + I, 5}");
    /* GaussianIntegers -> False keeps the ordinary integer behaviour. */
    check("Divisors[12, GaussianIntegers -> False]", "{1, 2, 3, 4, 6, 12}");
}

/* ---- Divisible cross-check: every Divisors entry divides n ---- */
void test_divisors_divisible_consistency() {
    check("And @@ (Divisible[360, #] & /@ Divisors[360])", "True");
    check("Length[Divisors[360]]", "24");   /* 2^3 3^2 5 -> 4*3*2 */
    /* Gaussian: every divisor of 6 + 4 I divides it. */
    check("Divisible[6 + 4 I, Divisors[6 + 4 I]]",
          "{True, True, True, True, True, True}");
}

/* ---- EulerPhi identity from the spec ---- */
void test_divisors_eulerphi_identity() {
    check("100 Product[1 - 1/k, {k, Select[Divisors[100], PrimeQ]}]", "40");
    check("EulerPhi[100]", "40");
}

/* ---- Arity, zero, and non-integer arguments stay unevaluated ---- */
void test_divisors_unevaluated() {
    check("Divisors[]", "Divisors[]");
    check("Divisors[x]", "Divisors[x]");
    check("Divisors[5/2]", "Divisors[5/2]");
    check("Divisors[2.5]", "Divisors[2.5]");
    check("Divisors[0]", "Divisors[0]");
    /* Astronomically large divisor count: must not crash, left unevaluated. */
    check("Length[Divisors[100!]]", "1");
}

/* ---- Attributes ---- */
void test_divisors_attributes() {
    check("Attributes[Divisors]", "{Listable, Protected}");
    check("MemberQ[Attributes[Divisors], Protected]", "True");
    check("Options[Divisors]", "{GaussianIntegers -> False}");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_divisors_basic);
    TEST(test_divisors_negative);
    TEST(test_divisors_listable);
    TEST(test_divisors_bignum);
    TEST(test_divisors_gaussian_auto);
    TEST(test_divisors_gaussian_option);
    TEST(test_divisors_divisible_consistency);
    TEST(test_divisors_eulerphi_identity);
    TEST(test_divisors_unevaluated);
    TEST(test_divisors_attributes);

    printf("All Divisors tests passed!\n");
    return 0;
}
