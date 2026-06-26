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
        printf("DivisorSigma test failed: %s\n  expected: %s\n  got:      %s\n",
               input, expected, res_str);
        ASSERT(0);
    }
    free(res_str);
    expr_free(e);
    expr_free(res);
}

/* ---- Basic integer sigma_k ---- */
void test_divisorsigma_basic() {
    check("DivisorSigma[1, 20]", "42");
    check("DivisorSigma[2, 20]", "546");
    check("DivisorSigma[2, 6]", "50");
    check("DivisorSigma[1, 6]", "12");
    check("DivisorSigma[1, 1]", "1");        /* empty product */
    check("DivisorSigma[1, 28]", "56");      /* perfect number */
    check("DivisorSigma[3, 4]", "73");       /* 1 + 8 + 64 */
    check("DivisorSigma[1, 7]", "8");        /* prime */
}

/* ---- sigma_0 is the divisor count ---- */
void test_divisorsigma_count() {
    check("DivisorSigma[0, 12]", "6");
    check("DivisorSigma[0, 1]", "1");
    check("DivisorSigma[0, 7]", "2");
    check("DivisorSigma[0, 60]", "12");
    /* Consistency with Divisors. */
    check("DivisorSigma[0, 360] == Length[Divisors[360]]", "True");
    check("DivisorSigma[1, 360] == Total[Divisors[360]]", "True");
    check("DivisorSigma[2, 84] == Total[Divisors[84]^2]", "True");
}

/* ---- Negative powers give rationals ---- */
void test_divisorsigma_negative() {
    check("DivisorSigma[-2, 10]", "13/10");
    check("DivisorSigma[-1, 6]", "2");       /* 1 + 1/2 + 1/3 + 1/6 */
    check("DivisorSigma[-1, 28] == 2", "True");  /* perfect: sigma_-1 == 2 */
}

/* ---- Sign of n is ignored ---- */
void test_divisorsigma_sign() {
    check("DivisorSigma[1, -20]", "42");
    check("DivisorSigma[2, -6]", "50");
}

/* ---- Big integers ---- */
void test_divisorsigma_bignum() {
    check("DivisorSigma[2, 10^30]",
          "1388888888888888888587721618754026546109040824596043713112451");
    check("DivisorSigma[1, 2^60] == 2^61 - 1", "True");  /* sigma_1(2^60) */
}

/* ---- Listable threading ---- */
void test_divisorsigma_listable() {
    check("DivisorSigma[2, {1, 2, 3, 4, 5}]", "{1, 5, 10, 21, 26}");
    check("DivisorSigma[1, {6, 28, 12}]", "{12, 56, 28}");
    check("MemberQ[Attributes[DivisorSigma], Listable]", "True");
}

/* ---- Attributes ---- */
void test_divisorsigma_attributes() {
    check("Attributes[DivisorSigma]", "{Listable, NHoldAll, Protected}");
    check("MemberQ[Attributes[DivisorSigma], Protected]", "True");
    check("MemberQ[Attributes[DivisorSigma], NHoldAll]", "True");
}

/* ---- Symbolic exponent: multiplicative formula ---- */
void test_divisorsigma_symbolic() {
    check("DivisorSigma[k, 2]", "(-1 + 2^(2 k))/(-1 + 2^k)");
    check("DivisorSigma[k, 1]", "1");
    /* Substitution back to a concrete exponent recovers the numeric value. */
    check("DivisorSigma[k, 20] /. k -> 2", "546");
    check("DivisorSigma[k, 20] /. k -> 1", "42");
    check("DivisorSigma[k, 6] /. k -> 2", "50");
}

/* ---- Rational exponent: matches the sum of rational powers of divisors ---- */
void test_divisorsigma_rational() {
    check("DivisorSigma[1/2, 12]",
          "(2 (-1 + 2 Sqrt[2]))/((-1 + Sqrt[2]) (-1 + Sqrt[3]))");
    /* Algebraically equal to Total[Divisors[12]^(1/2)]. */
    check("PossibleZeroQ[DivisorSigma[1/2, 12] - Total[Divisors[12]^(1/2)]]",
          "True");
    check("PossibleZeroQ[DivisorSigma[1/3, 18] - Total[Divisors[18]^(1/3)]]",
          "True");
}

/* ---- Gaussian integers ---- */
void test_divisorsigma_gaussian() {
    /* Non-real input auto-enables Z[i]. */
    check("DivisorSigma[1, 3 + I]", "2 + 6*I");
    /* GaussianIntegers -> True uses the multiplicative formula over Z[i]. */
    check("DivisorSigma[2, 6, GaussianIntegers -> True]", "-30 + 20*I");
    check("DivisorSigma[1, 6, GaussianIntegers -> True]", "8 + 12*I");
    check("DivisorSigma[2, 100, GaussianIntegers -> True]", "6479 - 6018*I");
    /* sigma_0 over Z[i] counts the divisor classes. */
    check("DivisorSigma[0, 5, GaussianIntegers -> True]", "4");
    check("DivisorSigma[0, 5, GaussianIntegers -> True] == "
          "Length[Divisors[5, GaussianIntegers -> True]]", "True");
}

/* ---- Unevaluated cases ---- */
void test_divisorsigma_unevaluated() {
    check("DivisorSigma[2, x]", "DivisorSigma[2, x]");
    check("DivisorSigma[k, n]", "DivisorSigma[k, n]");
    check("DivisorSigma[2, 0]", "DivisorSigma[2, 0]");
    check("DivisorSigma[2, 5/2]", "DivisorSigma[2, 5/2]");
    check("DivisorSigma[2, 2.5]", "DivisorSigma[2, 2.5]");
    /* Wrong argument counts stay unevaluated (with a DivisorSigma::argrx note). */
    check("DivisorSigma[]", "DivisorSigma[]");
    check("DivisorSigma[5]", "DivisorSigma[5]");
    check("DivisorSigma[1, 2, 3]", "DivisorSigma[1, 2, 3]");
}

int main() {
    symtab_init();
    core_init();

    TEST(test_divisorsigma_basic);
    TEST(test_divisorsigma_count);
    TEST(test_divisorsigma_negative);
    TEST(test_divisorsigma_sign);
    TEST(test_divisorsigma_bignum);
    TEST(test_divisorsigma_listable);
    TEST(test_divisorsigma_attributes);
    TEST(test_divisorsigma_symbolic);
    TEST(test_divisorsigma_rational);
    TEST(test_divisorsigma_gaussian);
    TEST(test_divisorsigma_unevaluated);

    printf("All DivisorSigma tests passed!\n");
    return 0;
}
