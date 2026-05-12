/*
 * test_extension_autodetect.c
 * ---------------------------
 * Unit tests for `extension_autodetect` (Phase G9): scan an expression
 * for algebraic-number generators and build the corresponding compositum
 * (a QATower).  Pure unit tests — no Cancel / Together / Simplify wiring
 * involved.  Each test:
 *
 *   1. Parses + evaluates an input string (so the expression is in its
 *      canonical picocas form, the same shape Cancel / Together would
 *      see at the entry to their extension-aware path).
 *   2. Calls `extension_autodetect`.
 *   3. Asserts the generator count and tower degree.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "expr.h"
#include "core.h"
#include "symtab.h"
#include "sym_names.h"
#include "qa.h"
#include "qafactor.h"
#include "test_utils.h"

/* Release builds compile assert() to a no-op via NDEBUG, which silently
 * elides side-effecting expressions inside ASSERT.  See tasks/lessons.md. */
#undef ASSERT
#define ASSERT(cond) do { \
    if (!(cond)) { \
        fprintf(stderr, "ASSERT FAIL at %s:%d: %s\n", __FILE__, __LINE__, #cond); \
        exit(1); \
    } \
} while (0)

/* Parse + evaluate `src`.  Caller owns the result. */
static Expr* parse_eval(const char* src) {
    Expr* p = parse_expression(src);
    ASSERT(p != NULL);
    Expr* e = evaluate(p);
    expr_free(p);
    return e;
}

/* Helper: assert exactly N generators detected and the compositum
 * has the expected degree.  Frees both the expression and the tower. */
static void assert_autodetect(const char* src, int expected_n,
                              size_t expected_deg) {
    Expr* e = parse_eval(src);
    QATower* t = extension_autodetect(e);
    if (!t) {
        fprintf(stderr, "FAIL: extension_autodetect returned NULL for: %s\n", src);
        expr_free(e);
        exit(1);
    }
    if (t->n != expected_n || t->ext->deg != expected_deg) {
        fprintf(stderr,
                "FAIL: %s  →  n=%d deg=%zu  (expected n=%d deg=%zu)\n",
                src, t->n, t->ext->deg, expected_n, expected_deg);
        qa_tower_free(t);
        expr_free(e);
        exit(1);
    }
    qa_tower_free(t);
    expr_free(e);
}

/* Helper: assert no extension is detected (input has no algebraic
 * generators, or has generators in shapes the tier-1 detector
 * doesn't support → NULL). */
static void assert_autodetect_null(const char* src) {
    Expr* e = parse_eval(src);
    QATower* t = extension_autodetect(e);
    if (t) {
        fprintf(stderr,
                "FAIL: %s  →  expected NULL, got n=%d deg=%zu\n",
                src, t->n, t->ext->deg);
        qa_tower_free(t);
        expr_free(e);
        exit(1);
    }
    expr_free(e);
}

/* ================================ Tests ================================ */

/* No algebraic generators at all. */
static void test_no_extension_rational(void) {
    assert_autodetect_null("3*x^2 + 2*x + 1");
}

static void test_no_extension_integer(void) {
    assert_autodetect_null("42");
}

static void test_no_extension_symbol(void) {
    assert_autodetect_null("x");
}

/* Single Sqrt — quadratic extension. */
static void test_single_sqrt2(void) {
    assert_autodetect("Sqrt[2]", 1, 2);
}

static void test_single_sqrt3_in_polynomial(void) {
    assert_autodetect("Sqrt[3] x + 1", 1, 2);
}

/* Sqrt[2] used several times — still one generator. */
static void test_repeated_sqrt2(void) {
    assert_autodetect("Sqrt[2] + 3 Sqrt[2] + Sqrt[2] x", 1, 2);
}

/* Cube root and squared cube root collapse to a single generator. */
static void test_cbrt2_and_sqr(void) {
    /* 2^(1/3) + 2^(2/3) → one generator, deg = 3 (LCM of denominators = 3). */
    assert_autodetect("Power[2, 1/3] + Power[2, 2/3]", 1, 3);
}

/* Sqrt[2] and 2^(1/3): same base, LCM(2,3) = 6 → one generator of degree 6. */
static void test_sqrt2_and_cbrt2_merge(void) {
    assert_autodetect("Sqrt[2] + Power[2, 1/3]", 1, 6);
}

/* Sqrt[2] and Sqrt[3]: distinct bases → two generators, compositum deg 4. */
static void test_sqrt2_and_sqrt3_distinct(void) {
    assert_autodetect("Sqrt[2] + Sqrt[3]", 2, 4);
}

/* 2^(1/3) and 3^(1/2): distinct bases, deg 6. */
static void test_cbrt2_and_sqrt3(void) {
    assert_autodetect("Power[2, 1/3] + Sqrt[3]", 2, 6);
}

/* Power[5, 1/4]: single quartic-root generator, deg 4. */
static void test_quartic_root(void) {
    assert_autodetect("Power[5, 1/4]", 1, 4);
}

/* Negative exponents still index into the same generator class.
 * Power[2, -1/3] is 1/2^(1/3), still a Q(2^(1/3)) element. */
static void test_negative_exponent_cbrt(void) {
    assert_autodetect("1 / Power[2, 1/3]", 1, 3);
}

/* Multiple sub-expressions, all over Q(Sqrt[2]).  Should detect a single
 * generator no matter how deeply nested. */
static void test_sqrt2_nested(void) {
    assert_autodetect("Sin[Sqrt[2] x] + Cos[Sqrt[2] / y]", 1, 2);
}

/* Polynomial radicand (with no atomic radical inside).  Tier-1 G8
 * `qa_resolve_nested_radical` requires at least one atomic radical in
 * the base (n_atoms == 0 → returns NULL), so these bail. */
static void test_sqrt_polynomial_bails(void) {
    assert_autodetect_null("Sqrt[x]");
}

static void test_sqrt_poly_arg_bails(void) {
    assert_autodetect_null("Sqrt[x^2 + 1]");
}

/* Nested-radical surfacing (Phase E): Sqrt[non-integer-base-with-
 * atomic-radicals] is now surfaced as a GEN_NESTED generator and
 * resolved by qa_resolve_nested_radical (Phase G8). */
static void test_sqrt_of_sqrt(void) {
    /* Sqrt[1 + Sqrt[2]] is a nested radical; G8 builds the degree-4
     * compositum minpoly y^4 - 2 y^2 - 1. */
    assert_autodetect("Sqrt[1 + Sqrt[2]]", 1, 4);
}

static void test_sqrt_of_cbrt(void) {
    /* Sqrt[2^(1/3)] = 2^(1/6); single generator of degree 6. */
    assert_autodetect("Sqrt[Power[2, 1/3]]", 1, 6);
}

/* Rational base → tier-1 picocas should have split this into integer
 * components, but be defensive: if a Rational base sneaks through, bail
 * cleanly.  Picocas normalises Sqrt[1/2] to 1/Sqrt[2] so we test the
 * normalised form. */
static void test_sqrt_half_normalised(void) {
    /* Picocas evaluates Sqrt[1/2] → 1/Sqrt[2] = Power[2, -1/2], so the
     * effective generator is Sqrt[2].  Confirms the canonicalisation
     * the auto-detector relies on. */
    assert_autodetect("Sqrt[1/2]", 1, 2);
}

/* The reduced motivating expression: x^3 - 2 carries no algebraic
 * generators (just polynomial coefficients).  But composing with a
 * 2^(1/3) should be detected. */
static void test_polynomial_alone(void) {
    assert_autodetect_null("x^3 - 2");
}

static void test_polynomial_times_cbrt(void) {
    assert_autodetect("(x^3 - 2) * Power[2, 1/3]", 1, 3);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_no_extension_rational);
    TEST(test_no_extension_integer);
    TEST(test_no_extension_symbol);

    TEST(test_single_sqrt2);
    TEST(test_single_sqrt3_in_polynomial);
    TEST(test_repeated_sqrt2);
    TEST(test_cbrt2_and_sqr);
    TEST(test_sqrt2_and_cbrt2_merge);
    TEST(test_sqrt2_and_sqrt3_distinct);
    TEST(test_cbrt2_and_sqrt3);
    TEST(test_quartic_root);
    TEST(test_negative_exponent_cbrt);
    TEST(test_sqrt2_nested);

    TEST(test_sqrt_polynomial_bails);
    TEST(test_sqrt_poly_arg_bails);
    TEST(test_sqrt_half_normalised);
    TEST(test_sqrt_of_sqrt);
    TEST(test_sqrt_of_cbrt);

    TEST(test_polynomial_alone);
    TEST(test_polynomial_times_cbrt);

    printf("All extension_autodetect tests passed!\n");
    return 0;
}
