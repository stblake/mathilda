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
        printf("LiouvilleLambda test failed: %s\n  expected: %s\n  got:      %s\n",
               input, expected, res_str);
        ASSERT(0);
    }
    free(res_str);
    expr_free(e);
    expr_free(res);
}

/* ---- Basic values: (-1)^Omega(n), Omega counted with multiplicity ---- */
void test_liouvillelambda_basic() {
    check("LiouvilleLambda[1]", "1");     /* unit, empty product */
    check("LiouvilleLambda[2]", "-1");    /* one prime, Omega = 1 */
    check("LiouvilleLambda[8]", "-1");    /* 2^3, Omega = 3 */
    check("LiouvilleLambda[9]", "1");     /* 3^2, Omega = 2 */
    check("LiouvilleLambda[6]", "1");     /* 2 * 3, Omega = 2 */
    check("LiouvilleLambda[12]", "-1");   /* 2^2 * 3, Omega = 3 */
    check("LiouvilleLambda[2 3 5]", "-1");   /* Omega = 3 */
    check("LiouvilleLambda[2 2 3 5]", "1");  /* Omega = 4 */
}

/* ---- Sign of n is ignored: lambda(-n) == lambda(n) ---- */
void test_liouvillelambda_negative() {
    check("LiouvilleLambda[-8]", "-1");
    check("LiouvilleLambda[-9]", "1");
    check("LiouvilleLambda[-8] == LiouvilleLambda[8]", "True");
}

/* ---- Completely multiplicative: lambda(m n) == lambda(m) lambda(n) ---- */
void test_liouvillelambda_multiplicative() {
    check("LiouvilleLambda[6 27] == LiouvilleLambda[6] LiouvilleLambda[27]", "True");
    check("LiouvilleLambda[10 21] == LiouvilleLambda[10] LiouvilleLambda[21]", "True");
}

/* ---- BigInt path: large argument ---- */
void test_liouvillelambda_bignum() {
    check("LiouvilleLambda[10^30 + 1]", "-1");
    /* (2^61 - 1) is a Mersenne prime -> Omega = 1 -> -1 */
    check("LiouvilleLambda[2^61 - 1]", "-1");
    /* Square of a large prime -> Omega = 2 -> 1 */
    check("LiouvilleLambda[(10^25 + 9)^2]", "1");
}

/* ---- Listable: threads element-wise over lists ---- */
void test_liouvillelambda_listable() {
    check("LiouvilleLambda[{1, 2, 3, 4, 5, 6}]", "{1, -1, -1, 1, -1, 1}");
    check("MemberQ[Attributes[LiouvilleLambda], Listable]", "True");
}

/* ---- Gaussian integers: auto-detected and via explicit option ---- */
void test_liouvillelambda_gaussian() {
    check("LiouvilleLambda[2 + I]", "-1");    /* norm 5 prime -> Gaussian prime */
    /* Over Z[i], 2 = -i (1 + i)^2, so Omega = 2 -> 1 */
    check("LiouvilleLambda[2, GaussianIntegers -> True]", "1");
    check("LiouvilleLambda[2, GaussianIntegers -> False]", "-1");
    /* Over Z[i], 8 = 2^3 has (1 + i)^6, Omega = 6 -> 1 */
    check("LiouvilleLambda[8, GaussianIntegers -> True]", "1");
    /* Explicit False on a Gaussian input forces the rational-integer path */
    check("LiouvilleLambda[4, GaussianIntegers -> False]", "1");
}

/* ---- Arity, non-integer, and zero arguments stay unevaluated ---- */
void test_liouvillelambda_unevaluated() {
    check("LiouvilleLambda[]", "LiouvilleLambda[]");
    check("LiouvilleLambda[2, 3]", "LiouvilleLambda[2, 3]");
    check("LiouvilleLambda[x]", "LiouvilleLambda[x]");
    check("LiouvilleLambda[5/2]", "LiouvilleLambda[5/2]");
    check("LiouvilleLambda[2.5]", "LiouvilleLambda[2.5]");
    check("LiouvilleLambda[0]", "LiouvilleLambda[0]");
}

/* ---- Attributes ---- */
void test_liouvillelambda_attributes() {
    check("Attributes[LiouvilleLambda]", "{Listable, Protected}");
    check("MemberQ[Attributes[LiouvilleLambda], Protected]", "True");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_liouvillelambda_basic);
    TEST(test_liouvillelambda_negative);
    TEST(test_liouvillelambda_multiplicative);
    TEST(test_liouvillelambda_bignum);
    TEST(test_liouvillelambda_listable);
    TEST(test_liouvillelambda_gaussian);
    TEST(test_liouvillelambda_unevaluated);
    TEST(test_liouvillelambda_attributes);

    printf("All LiouvilleLambda tests passed!\n");
    return 0;
}
