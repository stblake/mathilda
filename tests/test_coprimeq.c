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
        printf("CoprimeQ test failed: %s\n  expected: %s\n  got:      %s\n",
               input, expected, res_str);
        ASSERT(0);
    }
    free(res_str);
    expr_free(e);
    expr_free(res);
}

/* ---- Two machine integers ---- */
void test_coprimeq_pairs() {
    check("CoprimeQ[8, 11]", "True");
    check("CoprimeQ[2, 4]", "False");
    check("CoprimeQ[5, 6]", "True");
    check("CoprimeQ[9, 28]", "True");
    check("CoprimeQ[14, 21]", "False");   /* share 7 */
    check("CoprimeQ[1, 1]", "True");
    check("CoprimeQ[1, 100]", "True");    /* 1 is coprime to everything */
    check("CoprimeQ[17, 17]", "False");   /* equal primes share themselves */
}

/* ---- Negative arguments: coprimality ignores sign ---- */
void test_coprimeq_negative() {
    check("CoprimeQ[-3, 5]", "True");
    check("CoprimeQ[3, -5]", "True");
    check("CoprimeQ[-3, -5]", "True");
    check("CoprimeQ[-6, 4]", "False");    /* share 2 */
    check("CoprimeQ[2, 3, -5, 7]", "True");
}

/* ---- Zero edge cases: GCD(0, n) = |n|, so coprime only with units ---- */
void test_coprimeq_zero() {
    check("CoprimeQ[0, 1]", "True");      /* GCD(0,1) = 1 */
    check("CoprimeQ[0, 5]", "False");     /* GCD(0,5) = 5 */
    check("CoprimeQ[0, 0]", "False");     /* GCD(0,0) = 0 */
    check("CoprimeQ[1, 0]", "True");
}

/* ---- More than two arguments: pairwise test ---- */
void test_coprimeq_multiple() {
    check("CoprimeQ[2, 3, 5, 7]", "True");
    check("CoprimeQ[6, 35, 143]", "True");   /* 2*3, 5*7, 11*13 — pairwise coprime */
    check("CoprimeQ[2, 3, 4]", "False");     /* 2 and 4 share 2 */
    check("CoprimeQ[6, 35, 12]", "False");   /* 6 and 12 share factors */
    check("CoprimeQ[3, 5, 7, 11, 13]", "True");
    check("CoprimeQ[15, 28, 11]", "True");
    check("CoprimeQ[15, 28, 22]", "False");  /* 28 and 22 share 2 */
}

/* ---- BigInt integers and large-integer coprimality ---- */
void test_coprimeq_bignum() {
    check("CoprimeQ[2^100 - 1, 3^100 - 1]", "False");  /* both even */
    check("CoprimeQ[2^100, 3^100]", "True");
    check("CoprimeQ[2^127 - 1, 2^61 - 1]", "True");    /* two Mersenne primes */
    check("CoprimeQ[10^50, 10^50 + 1]", "True");       /* consecutive */
    check("CoprimeQ[6 * 10^40, 35 * 10^40]", "False"); /* both divisible by 10^40 */
    check("CoprimeQ[3^200, 5^200, 7^200]", "True");
}

/* ---- Gaussian integers (auto-detected from Complex arguments) ---- */
void test_coprimeq_gaussian_auto() {
    check("CoprimeQ[5 + I, 1 - I]", "False");  /* both divisible by 1 + I */
    check("CoprimeQ[2 + 3 I, 1 + I]", "True"); /* norms 13 (prime) and 2 */
    check("CoprimeQ[3 + 4 I, 1 + 2 I]", "True");
    check("CoprimeQ[1 + I, 1 - I]", "False");  /* associates: 1-I = -I(1+I) */
    check("CoprimeQ[2, 1 + I]", "False");      /* 2 = -I(1+I)^2 */
    check("CoprimeQ[3, 1 + I]", "True");       /* 3 is a Gaussian prime, norm 9 */
}

/* ---- GaussianIntegers -> True forces the ring even for plain integers ---- */
void test_coprimeq_gaussian_option() {
    /* 5 = (2+I)(2-I), 7 is a Gaussian prime: no common factor -> coprime. */
    check("CoprimeQ[5, 7, GaussianIntegers -> True]", "True");
    /* 2 = -I(1+I)^2 and 10 = 2*5 share the (1+I) factor. */
    check("CoprimeQ[2, 10, GaussianIntegers -> True]", "False");
    /* Over Z[i], 2 and 5 are coprime (different primes). */
    check("CoprimeQ[2, 5, GaussianIntegers -> True]", "True");
    /* GaussianIntegers -> False keeps the ordinary integer test. */
    check("CoprimeQ[8, 11, GaussianIntegers -> False]", "True");
    check("CoprimeQ[6, 10, GaussianIntegers -> False]", "False");
    /* Option may sit anywhere thanks to Orderless. */
    check("CoprimeQ[GaussianIntegers -> True, 3, 7]", "True");
}

/* ---- Listable: threads element-wise over lists ---- */
void test_coprimeq_listable() {
    check("CoprimeQ[{1, 2, 3, 4, 5}, 6]",
          "{True, False, False, False, True}");
    check("CoprimeQ[6, {5, 7, 8, 9, 11}]",
          "{True, True, False, False, True}");
    check("CoprimeQ[{4, 9, 25}, {3, 4, 6}]",
          "{True, True, True}");
    /* Threads recursively through nested lists. */
    check("CoprimeQ[{{4, 5}, {6, 7}}, 3]",
          "{{True, True}, {False, True}}");
}

/* ---- Single argument and empty call (documented contract) ---- */
void test_coprimeq_arity() {
    check("CoprimeQ[]", "False");
    check("CoprimeQ[7]", "True");
    check("CoprimeQ[12]", "True");
    check("CoprimeQ[0]", "True");
}

/* ---- Symbolic / non-integer arguments: never coprime -> False ---- */
void test_coprimeq_symbolic_false() {
    check("CoprimeQ[a, b]", "False");
    check("CoprimeQ[x, 5]", "False");
    check("CoprimeQ[2, 3, n]", "False");
    check("CoprimeQ[1/2, 3]", "False");   /* rationals are not handled */
    check("CoprimeQ[2.5, 3]", "False");   /* inexact reals are not handled */
    check("CoprimeQ[Pi, 2]", "False");
    /* An unrecognized option makes the result False. */
    check("CoprimeQ[2, 3, x -> 5]", "False");
}

/* ---- Cross-checks against GCD ---- */
void test_coprimeq_consistency() {
    check("CoprimeQ[5, 6] == (GCD[5, 6] == 1)", "True");
    check("CoprimeQ[4, 6] == (GCD[4, 6] == 1)", "True");
    check("CoprimeQ[14, 15] && (GCD[14, 15] == 1)", "True");
    check("CoprimeQ[2^100 - 1, 3^100 - 1] == (GCD[2^100 - 1, 3^100 - 1] == 1)", "True");
}

/* ---- Attributes: Listable, Orderless, Protected ---- */
void test_coprimeq_attributes() {
    check("Attributes[CoprimeQ]", "{Listable, Orderless, Protected}");
    check("MemberQ[Attributes[CoprimeQ], Listable]", "True");
    check("MemberQ[Attributes[CoprimeQ], Orderless]", "True");
    check("MemberQ[Attributes[CoprimeQ], Protected]", "True");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_coprimeq_pairs);
    TEST(test_coprimeq_negative);
    TEST(test_coprimeq_zero);
    TEST(test_coprimeq_multiple);
    TEST(test_coprimeq_bignum);
    TEST(test_coprimeq_gaussian_auto);
    TEST(test_coprimeq_gaussian_option);
    TEST(test_coprimeq_listable);
    TEST(test_coprimeq_arity);
    TEST(test_coprimeq_symbolic_false);
    TEST(test_coprimeq_consistency);
    TEST(test_coprimeq_attributes);

    printf("All CoprimeQ tests passed!\n");
    return 0;
}
