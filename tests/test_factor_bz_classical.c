/*
 * test_factor_bz_classical.c
 * --------------------------
 * Regression tests for the classical univariate factoriser over Z
 * (src/poly/facpoly_bz_uni.inc: squarefree decomposition, Cantor-Zassenhaus
 * mod p, quadratic Hensel lifting, Zassenhaus recombination, all on GMP
 * coefficients) and for the integer-content normalisation of the classical
 * polynomial GCD (poly_gcd_internal).
 *
 * The suite runs in both builds.  With FLINT, Factor of a univariate
 * polynomial goes through FLINT, so the Factor-level checks pin the shared
 * contract; the direct bz_factor_to_expr checks exercise the classical
 * factoriser itself in every build (it also serves the multivariate image
 * factoring and the rational-function path there).
 *
 * Every expected factor count is derived from the construction (products of
 * polynomials irreducible by Eisenstein, a negative discriminant, the
 * rational-root test, or the cyclotomic identity x^n - 1 = prod_{d|n} Phi_d,
 * whose count DivisorSigma[0, n] the kernel computes itself) -- never from a
 * remembered answer.  Correctness is checked by expanding the result back.
 *
 * Previously-failing inputs covered here: int64 overflow in the lift and in
 * c[j]*lc (big leading coefficients), the fixed 10^15 lift target (factors
 * with > 15-digit coefficients), the 1 << count subset mask (>= 31 modular
 * factors), the EDF exponent p^d overflow (x^43 - 1, x^128 - 1), the
 * degree-10000 allocation cap, and Factor[Expand[(x + 2)^40]] whose Yun
 * decomposition inherited a non-primitive GCD with exploding content.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "sym_names.h"
#include "test_utils.h"

#undef ASSERT
#define ASSERT(cond) do { \
    if (!(cond)) { \
        fprintf(stderr, "ASSERTION FAILED: %s\n  at %s:%d\n", \
                #cond, __FILE__, __LINE__); \
        abort(); \
    } \
} while (0)

/* The classical factoriser's public entry (src/poly/facpoly.c). */
extern Expr* bz_factor_to_expr(Expr* P, Expr* var);

/* Evaluate `input` and return its InputForm-ish string.  Caller frees. */
static char* eval_str(const char* input) {
    Expr* parsed = parse_expression(input);
    ASSERT(parsed != NULL);
    Expr* result = evaluate(parsed);
    char* got = expr_to_string(result);
    expr_free(parsed);
    expr_free(result);
    return got;
}

static void expect_true(const char* input) {
    char* got = eval_str(input);
    if (strcmp(got, "True") != 0) {
        fprintf(stderr, "FAIL: %s\n  got: %s\n", input, got);
        free(got);
        ASSERT(0);
    }
    free(got);
}

/* Defines the helpers used by the checks below:
 *   bzCnt[f]   number of irreducible non-constant factors, with multiplicity;
 *   bzOK[p, n] Factor[p] multiplies back to p and has exactly n factors. */
static void define_helpers(void) {
    char* r = eval_str(
        "bzCnt[f_] := Which[Head[f] === Times, Total[bzCnt /@ (List @@ f)], "
        "Head[f] === Power && IntegerQ[f[[2]]], Abs[f[[2]]] bzCnt[f[[1]]], "
        "NumericQ[f], 0, True, 1]");
    free(r);
    r = eval_str(
        "bzOK[p_, n_] := Module[{f = Factor[p]}, "
        "Expand[f - p] === 0 && bzCnt[f] === n]");
    free(r);
}

/* ---------------------------------------------------------------------- */

static void test_repeated_binomial_power(void) {
    /* Yun over a primitive GCD: no content blow-up. */
    expect_true("Factor[Expand[(x + 2)^40]] === (2 + x)^40");
    expect_true("Factor[Expand[-3 (x + 2)^7]] === -3 (2 + x)^7");
}

static void test_content_sign_and_monomial(void) {
    /* x^3 (3), (2x - 3)^2 (2), (x^2 + 1)^3 (3): 8 factors, content -6. */
    expect_true("bzOK[Expand[-6 x^3 (2 x - 3)^2 (x^2 + 1)^3], 8]");
    expect_true("Factor[Expand[-6 x^3 (2 x - 3)^2 (x^2 + 1)^3]] === "
                "-6 x^3 (-3 + 2 x)^2 (1 + x^2)^3");
    expect_true("Factor[1 - x^2] === -((-1 + x) (1 + x))");
}

static void test_big_leading_coefficients(void) {
    /* 1000003 x^3 + 2 x + 2: Eisenstein at 2.  1000005 x^2 - 3 x + 11:
     * negative discriminant.  The old int64 lift overflowed c[j] * lc. */
    expect_true("bzOK[Expand[(1000003 x^3 + 2 x + 2) (1000005 x^2 - 3 x + 11)], 2]");
    expect_true("bzOK[Expand[(10^12 x^3 + 2 x + 2) (10^12 + 7 + 10^12 x^2 + x)], 2]");
}

static void test_big_factor_coefficients(void) {
    /* x^2 + (10^20 + 3) x + 1: the discriminant lies strictly between
     * (10^20 + 2)^2 and (10^20 + 3)^2.  x^3 - 7 x + 5: no root among
     * +-1, +-5.  Factor coefficients far beyond the old 10^15 lift. */
    expect_true("bzOK[Expand[(x^2 + (10^20 + 3) x + 1) (x^3 - 7 x + 5)], 2]");
    expect_true("bzOK[Expand[(x^2 + (10^40 + 3) x + 1) (x^3 - 7 x + 5)], 2]");
}

static void test_many_modular_factors(void) {
    /* 31 linear factors: >= 31 modular factors broke the 1 << count mask. */
    expect_true("bzOK[Expand[Product[x - k, {k, 1, 31}]], 31]");
    expect_true("bzOK[Expand[Product[x - k, {k, -16, 16}]], 33]");
}

static void test_cyclotomic(void) {
    /* x^n - 1 has DivisorSigma[0, n] irreducible factors. */
    expect_true("And @@ Table[bzOK[x^n - 1, DivisorSigma[0, n]], "
                "{n, {12, 43, 79, 84, 105, 128, 240, 360}}]");
    /* x^n + 1 has DivisorSigma[0, 2n] - DivisorSigma[0, n] factors. */
    expect_true("And @@ Table[bzOK[x^n + 1, DivisorSigma[0, 2 n] - DivisorSigma[0, n]], "
                "{n, {32, 43, 105}}]");
    /* A perturbed cyclotomic product still factors (general path). */
    expect_true("bzOK[Expand[(x + 2) Sum[x^k, {k, 0, 42}]], 2]");
}

static void test_swinnerton_dyer(void) {
    /* The Swinnerton-Dyer polynomial for {2, 3, 5, 7} is irreducible over Q
     * but splits into factors of degree <= 2 modulo every prime. */
    expect_true("bzOK[Expand[Times @@ (x + {1, -1}[[#[[1]]]] Sqrt[2] + "
                "{1, -1}[[#[[2]]]] Sqrt[3] + {1, -1}[[#[[3]]]] Sqrt[5] + "
                "{1, -1}[[#[[4]]]] Sqrt[7] & /@ Tuples[{1, 2}, 4])], 1]");
}

static void test_gcd_primitive(void) {
    /* poly_gcd_internal must return a primitive GCD (positive lc), not one
     * carrying the integer content of the remainder sequence. */
    expect_true("Module[{p = Expand[(3 + x - 2 x^2 + x^5)^2 (2 x - 5)]}, "
                "Expand[PolynomialGCD[p, D[p, x]]] === 3 + x - 2 x^2 + x^5]");
    expect_true("Expand[PolynomialGCD[Expand[(10^30 x + 1) (x - 2)], "
                "Expand[(10^30 x + 1) (x + 3) 10^25]]] === 1 + 10^30 x");
}

static void test_rational_function_univariate_numerator(void) {
    expect_true("Module[{r = Together[Expand[Product[x - k, {k, 1, 31}]]/(x + y)]}, "
                "Together[Factor[r] - r] === 0 && "
                "bzCnt[Numerator[Factor[r]]] === 31]");
}

/* Direct calls: these always reach the classical factoriser. */
static void check_bz_direct(const char* poly, int expected_factors) {
    Expr* var = expr_new_symbol("x");
    Expr* parsed = parse_expression(poly);
    ASSERT(parsed != NULL);
    Expr* P = evaluate(parsed);
    expr_free(parsed);
    Expr* F = bz_factor_to_expr(P, var);
    ASSERT(F != NULL);

    /* Multiplies back to P. */
    Expr* diff = expr_new_function(expr_new_symbol("Expand"),
        (Expr*[]){ expr_new_function(expr_new_symbol("Subtract"),
                                     (Expr*[]){ expr_copy(F), expr_copy(P) }, 2) }, 1);
    Expr* d = evaluate(diff);
    expr_free(diff);
    if (!(d->type == EXPR_INTEGER && d->data.integer == 0)) {
        char* s = expr_to_string(F);
        fprintf(stderr, "FAIL: bz_factor_to_expr[%s] = %s does not expand back\n", poly, s);
        free(s);
        ASSERT(0);
    }
    expr_free(d);

    /* Factor count, with multiplicity. */
    Expr* cnt = expr_new_function(expr_new_symbol("bzCnt"), (Expr*[]){ expr_copy(F) }, 1);
    Expr* c = evaluate(cnt);
    expr_free(cnt);
    if (!(c->type == EXPR_INTEGER && c->data.integer == expected_factors)) {
        char* s = expr_to_string(F);
        fprintf(stderr, "FAIL: bz_factor_to_expr[%s] = %s, expected %d factors\n",
                poly, s, expected_factors);
        free(s);
        ASSERT(0);
    }
    expr_free(c);
    expr_free(F);
    expr_free(P);
    expr_free(var);
}

static void test_bz_direct(void) {
    check_bz_direct("Expand[(x + 2)^40]", 40);
    check_bz_direct("Expand[-6 x^3 (2 x - 3)^2 (x^2 + 1)^3]", 8);
    check_bz_direct("Expand[(1000003 x^3 + 2 x + 2) (1000005 x^2 - 3 x + 11)]", 2);
    check_bz_direct("Expand[(x^2 + (10^20 + 3) x + 1) (x^3 - 7 x + 5)]", 2);
    check_bz_direct("Expand[(x^2 + (10^40 + 3) x + 1)^2 (x^3 - 7 x + 5)]", 3);
    check_bz_direct("Expand[Product[x - k, {k, 1, 31}]]", 31);
    check_bz_direct("x^43 - 1", 2);
    check_bz_direct("x^128 - 1", 8);
    check_bz_direct("x^84 - 1", 12);
    check_bz_direct("x^105 + 1", 8);
    check_bz_direct("Expand[(x + 2) Sum[x^k, {k, 0, 42}]]", 2);
    /* Eisenstein at 2 (x^k + 2, k = 3, 5, 7) times a product of linears. */
    check_bz_direct("Expand[(x^3 + 2) (x^5 + 2) (x^7 + 2) (x - 1) (x + 1) (x - 3)]", 6);
    /* Irreducible input comes back unchanged. */
    check_bz_direct("x^10001 + 2", 1);
}

int main(void) {
    symtab_init();
    core_init();
    define_helpers();

    printf("Running classical univariate factoriser tests...\n");
    TEST(test_repeated_binomial_power);
    TEST(test_content_sign_and_monomial);
    TEST(test_big_leading_coefficients);
    TEST(test_big_factor_coefficients);
    TEST(test_many_modular_factors);
    TEST(test_cyclotomic);
    TEST(test_swinnerton_dyer);
    TEST(test_gcd_primitive);
    TEST(test_rational_function_univariate_numerator);
    TEST(test_bz_direct);
    printf("All classical univariate factoriser tests passed!\n");
    return 0;
}
