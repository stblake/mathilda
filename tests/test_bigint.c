#include "core.h"
#include "expr.h"
#include "symtab.h"
#include "eval.h"
#include "test_utils.h"
#include <stdio.h>
#include <string.h>
#include "print.h"
#include "parse.h"
#include <gmp.h>

void test_large_literal_parsing(void) {
    // A number beyond INT64_MAX should parse as EXPR_BIGINT
    Expr* e = parse_expression("99999999999999999999");
    assert(e != NULL);
    assert(e->type == EXPR_BIGINT);
    char* s = expr_to_string(e);
    assert(strcmp(s, "99999999999999999999") == 0);
    free(s);
    expr_free(e);

    // INT64_MAX should still parse as EXPR_INTEGER
    Expr* e2 = parse_expression("9223372036854775807");
    assert(e2 != NULL);
    assert(e2->type == EXPR_INTEGER);
    expr_free(e2);
}

void test_integer_overflow_promotion(void) {
    // INT64_MAX + 1 should auto-promote to EXPR_BIGINT
    assert_eval_eq("Plus[9223372036854775807, 1]", "9223372036854775808", 0);

    // Verify the result type is bigint
    Expr* e = parse_expression("Plus[9223372036854775807, 1]");
    Expr* res = evaluate(e);
    assert(res->type == EXPR_BIGINT);
    expr_free(e);
    expr_free(res);
}

void test_factorial_correctness(void) {
    // 20! = 2432902008176640000 (fits in int64)
    assert_eval_eq("Factorial[20]", "2432902008176640000", 0);

    // 21! = 51090942171709440000 (overflows int64, should be bigint)
    assert_eval_eq("Factorial[21]", "51090942171709440000", 0);

    // 100! should be a large exact integer
    Expr* e = parse_expression("Factorial[100]");
    Expr* res = evaluate(e);
    assert(res->type == EXPR_BIGINT);
    char* s = expr_to_string(res);
    // 100! has 158 digits
    assert(strlen(s) == 158);
    // Should start with "9332621544"
    assert(strncmp(s, "9332621544", 10) == 0);
    free(s);
    expr_free(e);
    expr_free(res);
}

void test_large_multiplication(void) {
    // Large * large should give exact result
    assert_eval_eq("Times[9999999999999999999, 9999999999999999999]",
                    "99999999999999999980000000000000000001", 0);
}

void test_integerq_bigint(void) {
    // IntegerQ on a factorial result (bigint)
    assert_eval_eq("IntegerQ[Factorial[100]]", "True", 0);
    // IntegerQ on a normal integer
    assert_eval_eq("IntegerQ[42]", "True", 0);
    // IntegerQ on a non-integer
    assert_eval_eq("IntegerQ[3.14]", "False", 0);
}

void test_numberq_bigint(void) {
    assert_eval_eq("NumberQ[Factorial[100]]", "True", 0);
}

void test_numericq_bigint(void) {
    assert_eval_eq("NumericQ[Factorial[100]]", "True", 0);
}

void test_evenq_oddq_bigint(void) {
    // 2^100 is even
    assert_eval_eq("EvenQ[Power[2, 100]]", "True", 0);
    assert_eval_eq("OddQ[Power[2, 100]]", "False", 0);

    // 21! is even (since it contains factor 2)
    assert_eval_eq("EvenQ[Factorial[21]]", "True", 0);

    // 3^67 is odd
    assert_eval_eq("OddQ[Power[3, 67]]", "True", 0);
    assert_eval_eq("EvenQ[Power[3, 67]]", "False", 0);
}

void test_bigint_printing(void) {
    // 2^100 = 1267650600228229401496703205376 (31 digits)
    Expr* e = parse_expression("Power[2, 100]");
    Expr* res = evaluate(e);
    assert(res->type == EXPR_BIGINT);
    char* s = expr_to_string(res);
    assert(strlen(s) == 31);
    assert(strncmp(s, "12676506002", 11) == 0);
    free(s);
    expr_free(e);
    expr_free(res);
}

void test_power_large_exponent(void) {
    // 2^100 should be exact 31-digit number
    assert_eval_eq("Power[2, 100]", "1267650600228229401496703205376", 0);
}

void test_abs_bigint(void) {
    // Abs of a negative bigint literal
    assert_eval_eq("Abs[-99999999999999999999]", "99999999999999999999", 0);
}

/* Regression: Plus/Times must handle Rational with BigInt
 * components without crashing.  Pre-fix this segfaulted in
 * builtin_plus -> is_overflow(NULL) because add_numbers /
 * multiply_numbers fell through to `return NULL` when neither the
 * int64-rational nor the BigInt-typed operand path matched, e.g.
 * Rational[BigInt, Integer] + Rational[BigInt, Integer]. */
void test_plus_rational_with_bigint_parts(void) {
    /* Build huge rationals expressed with BigInt numerators above
     * INT64_MAX, and sum them.  Pre-fix these segfaulted instead
     * of returning a value. */
    assert_eval_eq(
        "(99999999999999999999 / 7) + (1 / 7)",
        "100000000000000000000/7", 0);
    /* Multiplication path: Integer * Rational[BigInt, Integer] */
    assert_eval_eq(
        "3 * (99999999999999999999 / 7)",
        "299999999999999999997/7", 0);
    /* Mixed Integer + Rational[BigInt, BigInt] */
    assert_eval_eq(
        "1 + (99999999999999999999 / 99999999999999999998)",
        "199999999999999999997/99999999999999999998", 0);
}

void test_squfof(void) {
    // Standard test cases for ShanksSquareForms
    assert_eval_eq("FactorInteger[11111, Method -> \"ShanksSquareForms\"]", "{{41, 1}, {271, 1}}", 0);
    assert_eval_eq("FactorInteger[11111111111, Method -> \"ShanksSquareForms\"]", "{{21649, 1}, {513239, 1}}", 0);
    // Should still work if factors are small
    assert_eval_eq("FactorInteger[15, Method -> \"ShanksSquareForms\"]", "{{3, 1}, {5, 1}}", 0);
}

/* Issue #95: a 39-digit semiprime (19- and 21-digit factors) that the default
 * method returned unfactored in most runs, because every ECM "retry" replayed
 * the same curve. Repeated so a single lucky curve cannot mask a regression;
 * the ECM-only method exercises the fresh-curve fix without the p-1 pre-pass. */
void test_ecm_semiprime(void) {
#ifndef NO_ECM
    for (int i = 0; i < 3; i++) {
        assert_eval_eq("FactorInteger[(2^61-1)*(10^20+39)]",
                       "{{2305843009213693951, 1}, {100000000000000000039, 1}}", 0);
        assert_eval_eq("FactorInteger[(2^61-1)*(10^20+39), Method -> \"ECM\"]",
                       "{{2305843009213693951, 1}, {100000000000000000039, 1}}", 0);
        assert_eval_eq("FactorInteger[10000000000000012363*300000000000000000797]",
                       "{{10000000000000012363, 1}, {300000000000000000797, 1}}", 0);
    }
#endif
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_large_literal_parsing);
    TEST(test_integer_overflow_promotion);
    TEST(test_factorial_correctness);
    TEST(test_large_multiplication);
    TEST(test_integerq_bigint);
    TEST(test_numberq_bigint);
    TEST(test_numericq_bigint);
    TEST(test_evenq_oddq_bigint);
    TEST(test_bigint_printing);
    TEST(test_power_large_exponent);
    TEST(test_abs_bigint);
    TEST(test_plus_rational_with_bigint_parts);
    TEST(test_squfof);
    TEST(test_ecm_semiprime);

    // Additional GMP integration tests
    assert_eval_eq("Binomial[100, 50]", "100891344545564193334812497256", 0);
    /* Binomial of big arguments (machine-word indices take mpz_bin_uiui):
     * against the factorial quotient, Pascal's rule, the row sum 2^n, the
     * symmetry, and the extensions to negative arguments */
    assert_eval_eq("Table[Binomial[n, k] == n!/(k! (n - k)!), {n, {255, 256, 355, 511, 1000}}, "
                   "{k, {0, 1, 2, 100, 128, 150, 255}}] // Flatten // Union", "{True}", 0);
    assert_eval_eq("Union[Table[Binomial[400, k] == Binomial[399, k - 1] + Binomial[399, k], {k, 1, 399}]]",
                   "{True}", 0);
    assert_eval_eq("Sum[Binomial[300, k], {k, 0, 300}] == 2^300", "True", 0);
    assert_eval_eq("Binomial[355, 150] == Binomial[355, 205]", "True", 0);
    assert_eval_eq("Binomial[3000, 1500] == 3000!/(1500!)^2", "True", 0);
    assert_eval_eq("{Binomial[300, 301], Binomial[300, -1], Binomial[0, 0]}", "{0, 0, 1}", 0);
    /* the machine-word path: values near and past 2^63 (an overflowing step
     * hands over to GMP), against the factorial quotient */
    assert_eval_eq("{Binomial[19, 9], Binomial[62, 31], Binomial[66, 33], Binomial[67, 33], "
                   "Binomial[1000, 3], Binomial[10^6, 2], Binomial[5, 0], Binomial[5, 5], Binomial[3, 5]}",
                   "{92378, 465428353255261088, 7219428434016265740, 14226520737620288370, "
                   "166167000, 499999500000, 1, 1, 0}", 0);
    assert_eval_eq("Union[Flatten[Table[Binomial[n, k] == n!/(k! (n - k)!), {n, 0, 70}, {k, 0, n}]]]",
                   "{True}", 0);
    assert_eval_eq("Binomial[-5, 3]", "-35", 0);
    assert_eval_eq("Binomial[-300, 4] == (-300) (-301) (-302) (-303)/24", "True", 0);
    assert_eval_eq("{Binomial[-1, -1], Binomial[-2, -5]}", "{1, -4}", 0);
    assert_eval_eq("Binomial[2^70 + 1, 2] == (2^70 + 1) 2^70/2", "True", 0);
    assert_eval_eq("Binomial[2^64 + 7, 3] == (2^64 + 7) (2^64 + 6) (2^64 + 5)/6", "True", 0);
    assert_eval_eq("PowerMod[2, 100, 101]", "1", 0); // Fermat's Little Theorem: 2^(101-1) mod 101 = 1
    assert_eval_eq("PowerMod[2, 10, 3]", "1", 0);
    assert_eval_eq("PowerMod[3, -2, 7]", "4", 0);
    assert_eval_eq("Mod[10^30, 7]", "1", 0);

    printf("\nAll bigint tests passed!\n");
    return 0;
}
