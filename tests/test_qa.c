/*
 * test_qa.c
 * ---------
 * Unit tests for the QA (Q(α) algebraic-number) substrate, Phase G1 of
 * Trager's algebraic factoring (FACTOR_PLAN.md §14).
 *
 * Coverage:
 *   - QAExt construction (raw + Sqrt[c] convenience).
 *   - QANum construction (zero, one, alpha, from_si, copy).
 *   - Predicates (is_zero, is_one, eq).
 *   - Arithmetic (add, sub, neg, scale).
 *   - Multiplication with reduction modulo P_α (the core invariant).
 *   - Inversion via extended Euclidean over Q[y].
 *   - Division (a / b == a * b^-1, when b is invertible).
 *
 * The substrate is field-agnostic, so we exercise it across three
 * representative extensions:
 *
 *   (a) Q(√2) — degree 2, simplest non-trivial extension.
 *   (b) Q(2^(1/3)) — degree 3, exercises higher-degree reduction.
 *   (c) Q(i) = Q[y]/(y² + 1) — complex extension, same machinery.
 *
 * Each test asserts both predicate outcomes and exact-rational
 * coefficient equality of the result.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include "qa.h"
#include "test_utils.h"

/* Build (num/den) as a freshly-initialised mpq_t. Caller mpq_clear-s. */
static void mq_set(mpq_t out, long num, long den) {
    mpq_set_si(out, num, (unsigned long)den);
    mpq_canonicalize(out);
}

/* Assert that `a` has coefficient `num/den` at position i. */
static void assert_coef_si(const QANum* a, size_t i, long num, long den) {
    mpq_t expected;
    mpq_init(expected);
    mq_set(expected, num, den);
    if (!mpq_equal(a->coef[i], expected)) {
        char* got = mpq_get_str(NULL, 10, a->coef[i]);
        char* exp = mpq_get_str(NULL, 10, expected);
        fprintf(stderr, "coef[%zu] mismatch: expected %s, got %s\n",
                i, exp, got);
        free(got); free(exp);
        ASSERT(0);
    }
    mpq_clear(expected);
}

/* ====================== QAExt construction ============================== */

static void test_qaext_sqrt_si(void) {
    QAExt* ext = qaext_sqrt_si(2);
    ASSERT(ext->deg == 2);
    /* P_α(y) = y² - 2 */
    mpq_t expected;
    mpq_init(expected);
    mpq_set_si(expected, -2, 1); ASSERT(mpq_equal(ext->coef[0], expected));
    mpq_set_si(expected, 0, 1);  ASSERT(mpq_equal(ext->coef[1], expected));
    mpq_set_si(expected, 1, 1);  ASSERT(mpq_equal(ext->coef[2], expected));
    mpq_clear(expected);
    qaext_free(ext);
}

static void test_qaext_root_si(void) {
    QAExt* ext = qaext_root_si(2, 3);
    ASSERT(ext->deg == 3);
    /* P_α(y) = y³ - 2 */
    mpq_t expected;
    mpq_init(expected);
    mpq_set_si(expected, -2, 1); ASSERT(mpq_equal(ext->coef[0], expected));
    mpq_set_si(expected, 0, 1);  ASSERT(mpq_equal(ext->coef[1], expected));
    mpq_set_si(expected, 0, 1);  ASSERT(mpq_equal(ext->coef[2], expected));
    mpq_set_si(expected, 1, 1);  ASSERT(mpq_equal(ext->coef[3], expected));
    mpq_clear(expected);
    qaext_free(ext);
}

/* ====================== QANum construction & predicates ================ */

static void test_qa_zero_one_alpha(void) {
    QAExt* ext = qaext_sqrt_si(2);
    QANum* z = qa_zero(ext);
    QANum* o = qa_one(ext);
    QANum* a = qa_alpha(ext);
    ASSERT(qa_is_zero(z));
    ASSERT(qa_is_one(o));
    ASSERT(!qa_is_zero(a));
    ASSERT(!qa_is_one(a));
    /* a represents √2 = 0 + 1·α */
    assert_coef_si(a, 0, 0, 1);
    assert_coef_si(a, 1, 1, 1);
    qa_free(z); qa_free(o); qa_free(a);
    qaext_free(ext);
}

static void test_qa_from_si(void) {
    QAExt* ext = qaext_sqrt_si(2);
    QANum* a = qa_from_si(ext, 3, 4);
    /* 3/4 + 0·α */
    assert_coef_si(a, 0, 3, 4);
    assert_coef_si(a, 1, 0, 1);
    qa_free(a);
    qaext_free(ext);
}

static void test_qa_copy_eq(void) {
    QAExt* ext = qaext_sqrt_si(3);
    QANum* a = qa_alpha(ext);
    QANum* b = qa_copy(a);
    ASSERT(qa_eq(a, b));
    QANum* c = qa_one(ext);
    ASSERT(!qa_eq(a, c));
    qa_free(a); qa_free(b); qa_free(c);
    qaext_free(ext);
}

/* ============================= Arithmetic ============================== */

static void test_qa_add_sub_neg(void) {
    QAExt* ext = qaext_sqrt_si(2);
    /* a = 1 + 2α, b = 3 - α */
    QANum* a = qa_zero(ext);
    mpq_set_si(a->coef[0], 1, 1);
    mpq_set_si(a->coef[1], 2, 1);
    QANum* b = qa_zero(ext);
    mpq_set_si(b->coef[0], 3, 1);
    mpq_set_si(b->coef[1], -1, 1);

    QANum* sum  = qa_add(a, b);   /* 4 + α */
    QANum* diff = qa_sub(a, b);   /* -2 + 3α */
    QANum* neg  = qa_neg(a);      /* -1 - 2α */

    assert_coef_si(sum,  0,  4, 1); assert_coef_si(sum,  1,  1, 1);
    assert_coef_si(diff, 0, -2, 1); assert_coef_si(diff, 1,  3, 1);
    assert_coef_si(neg,  0, -1, 1); assert_coef_si(neg,  1, -2, 1);

    qa_free(a); qa_free(b); qa_free(sum); qa_free(diff); qa_free(neg);
    qaext_free(ext);
}

static void test_qa_scale(void) {
    QAExt* ext = qaext_sqrt_si(2);
    /* a = 2 + 4α; scale by 3/2 => 3 + 6α */
    QANum* a = qa_zero(ext);
    mpq_set_si(a->coef[0], 2, 1);
    mpq_set_si(a->coef[1], 4, 1);
    QANum* r = qa_scale_si(a, 3, 2);
    assert_coef_si(r, 0, 3, 1);
    assert_coef_si(r, 1, 6, 1);
    qa_free(a); qa_free(r);
    qaext_free(ext);
}

/* ====================== Multiplication with reduction =================== */

static void test_qa_mul_sqrt2_squared(void) {
    /* The headline invariant: in Q(√2), √2 · √2 = 2.
     * α · α = α² ≡ 2 (mod y² - 2). */
    QAExt* ext = qaext_sqrt_si(2);
    QANum* a = qa_alpha(ext);
    QANum* aa = qa_mul(a, a);
    /* Expected: 2 + 0·α */
    assert_coef_si(aa, 0, 2, 1);
    assert_coef_si(aa, 1, 0, 1);
    qa_free(a); qa_free(aa);
    qaext_free(ext);
}

static void test_qa_mul_distributive(void) {
    /* (1 + α)(1 - α) = 1 - α² = 1 - 2 = -1 in Q(√2). */
    QAExt* ext = qaext_sqrt_si(2);
    QANum* p = qa_zero(ext);
    mpq_set_si(p->coef[0], 1, 1);
    mpq_set_si(p->coef[1], 1, 1);
    QANum* q = qa_zero(ext);
    mpq_set_si(q->coef[0], 1, 1);
    mpq_set_si(q->coef[1], -1, 1);
    QANum* pq = qa_mul(p, q);
    assert_coef_si(pq, 0, -1, 1);
    assert_coef_si(pq, 1, 0, 1);
    qa_free(p); qa_free(q); qa_free(pq);
    qaext_free(ext);
}

static void test_qa_mul_cube_root(void) {
    /* In Q(2^(1/3)), α · α · α = α³ ≡ 2.
     * Compute α² first then multiply by α: (α²) · α = α³. */
    QAExt* ext = qaext_root_si(2, 3);
    QANum* a = qa_alpha(ext);          /* α */
    QANum* aa = qa_mul(a, a);          /* α² */
    /* Expected α²: 0 + 0·α + 1·α² */
    assert_coef_si(aa, 0, 0, 1);
    assert_coef_si(aa, 1, 0, 1);
    assert_coef_si(aa, 2, 1, 1);

    QANum* aaa = qa_mul(aa, a);        /* α³ ≡ 2 */
    assert_coef_si(aaa, 0, 2, 1);
    assert_coef_si(aaa, 1, 0, 1);
    assert_coef_si(aaa, 2, 0, 1);
    qa_free(a); qa_free(aa); qa_free(aaa);
    qaext_free(ext);
}

static void test_qa_mul_imaginary(void) {
    /* Q(i) = Q[y]/(y² + 1). i² = -1. */
    QAExt* ext = qaext_sqrt_si(-1);
    QANum* i = qa_alpha(ext);
    QANum* ii = qa_mul(i, i);
    assert_coef_si(ii, 0, -1, 1);
    assert_coef_si(ii, 1, 0, 1);

    /* (3 + 4i)(3 - 4i) = 9 - 16i² = 9 + 16 = 25. */
    QANum* p = qa_zero(ext);
    mpq_set_si(p->coef[0], 3, 1);
    mpq_set_si(p->coef[1], 4, 1);
    QANum* q = qa_zero(ext);
    mpq_set_si(q->coef[0], 3, 1);
    mpq_set_si(q->coef[1], -4, 1);
    QANum* pq = qa_mul(p, q);
    assert_coef_si(pq, 0, 25, 1);
    assert_coef_si(pq, 1, 0, 1);

    qa_free(i); qa_free(ii); qa_free(p); qa_free(q); qa_free(pq);
    qaext_free(ext);
}

/* ============================== Inversion ============================== */

static void test_qa_inverse_sqrt2(void) {
    /* In Q(√2), 1/√2 = √2 / 2 = 0 + (1/2)·α.
     * Verify: α · (α/2) ≡ α²/2 = 2/2 = 1. */
    QAExt* ext = qaext_sqrt_si(2);
    QANum* a = qa_alpha(ext);
    QANum* ainv = qa_inverse(a);
    ASSERT(ainv != NULL);
    /* ainv should equal (1/2) α. */
    assert_coef_si(ainv, 0, 0, 1);
    assert_coef_si(ainv, 1, 1, 2);

    QANum* prod = qa_mul(a, ainv);
    ASSERT(qa_is_one(prod));
    qa_free(a); qa_free(ainv); qa_free(prod);
    qaext_free(ext);
}

static void test_qa_inverse_general(void) {
    /* Inverse of (1 + √2) in Q(√2).
     * (1 + α)(? + ?·α) = 1   with α² = 2.
     * (1 + √2)(√2 - 1) = 2 - 1 = 1, so 1/(1+√2) = √2 - 1. */
    QAExt* ext = qaext_sqrt_si(2);
    QANum* a = qa_zero(ext);
    mpq_set_si(a->coef[0], 1, 1);
    mpq_set_si(a->coef[1], 1, 1);
    QANum* ainv = qa_inverse(a);
    ASSERT(ainv != NULL);
    /* Expected: -1 + 1·α */
    assert_coef_si(ainv, 0, -1, 1);
    assert_coef_si(ainv, 1, 1, 1);

    QANum* prod = qa_mul(a, ainv);
    ASSERT(qa_is_one(prod));

    qa_free(a); qa_free(ainv); qa_free(prod);
    qaext_free(ext);
}

static void test_qa_inverse_cube_root(void) {
    /* Inverse of α in Q(2^(1/3)). α · ? = 1.
     * α · α² / 2 = α³ / 2 = 2/2 = 1, so 1/α = α²/2. */
    QAExt* ext = qaext_root_si(2, 3);
    QANum* a = qa_alpha(ext);
    QANum* ainv = qa_inverse(a);
    ASSERT(ainv != NULL);
    assert_coef_si(ainv, 0, 0, 1);
    assert_coef_si(ainv, 1, 0, 1);
    assert_coef_si(ainv, 2, 1, 2);

    QANum* prod = qa_mul(a, ainv);
    ASSERT(qa_is_one(prod));

    qa_free(a); qa_free(ainv); qa_free(prod);
    qaext_free(ext);
}

static void test_qa_inverse_zero_is_null(void) {
    QAExt* ext = qaext_sqrt_si(2);
    QANum* z = qa_zero(ext);
    ASSERT(qa_inverse(z) == NULL);
    qa_free(z);
    qaext_free(ext);
}

/* ============================== Division =============================== */

static void test_qa_div(void) {
    /* (3 + 4i) / (1 + 2i) = ((3+4i)(1-2i)) / 5 = (3 - 6i + 4i + 8)/5 =
     *                       (11 - 2i)/5 = 11/5 + (-2/5) i. */
    QAExt* ext = qaext_sqrt_si(-1);
    QANum* num = qa_zero(ext);
    mpq_set_si(num->coef[0], 3, 1);
    mpq_set_si(num->coef[1], 4, 1);
    QANum* den = qa_zero(ext);
    mpq_set_si(den->coef[0], 1, 1);
    mpq_set_si(den->coef[1], 2, 1);
    QANum* q = qa_div(num, den);
    ASSERT(q != NULL);
    assert_coef_si(q, 0, 11, 5);
    assert_coef_si(q, 1, -2, 5);

    /* Verify by multiplying back. */
    QANum* check = qa_mul(q, den);
    ASSERT(qa_eq(check, num));

    qa_free(num); qa_free(den); qa_free(q); qa_free(check);
    qaext_free(ext);
}

/* =================== Composition / cumulative tests ===================== */

static void test_qa_pythagoras_in_qsqrt2(void) {
    /* (1 + √2)² + (1 - √2)² should equal 6 (real part) + 0·α. */
    QAExt* ext = qaext_sqrt_si(2);
    QANum* p = qa_zero(ext);
    mpq_set_si(p->coef[0], 1, 1); mpq_set_si(p->coef[1], 1, 1);
    QANum* q = qa_zero(ext);
    mpq_set_si(q->coef[0], 1, 1); mpq_set_si(q->coef[1], -1, 1);
    QANum* pp = qa_mul(p, p);                    /* (1+√2)² = 3 + 2√2 */
    QANum* qq = qa_mul(q, q);                    /* (1-√2)² = 3 - 2√2 */
    QANum* sum = qa_add(pp, qq);                 /* 6 + 0·α */
    assert_coef_si(sum, 0, 6, 1);
    assert_coef_si(sum, 1, 0, 1);
    qa_free(p); qa_free(q); qa_free(pp); qa_free(qq); qa_free(sum);
    qaext_free(ext);
}

static void test_qa_division_round_trip(void) {
    /* For a few non-trivial elements of Q(2^(1/3)), check (a/b)*b == a. */
    QAExt* ext = qaext_root_si(2, 3);
    QANum* a = qa_zero(ext);
    mpq_set_si(a->coef[0], 5, 1);
    mpq_set_si(a->coef[1], 7, 3);
    mpq_set_si(a->coef[2], -1, 2);
    mpq_canonicalize(a->coef[1]);
    mpq_canonicalize(a->coef[2]);

    QANum* b = qa_zero(ext);
    mpq_set_si(b->coef[0], 1, 1);
    mpq_set_si(b->coef[1], 2, 1);
    mpq_set_si(b->coef[2], -3, 1);

    QANum* ab = qa_div(a, b);
    ASSERT(ab != NULL);
    QANum* check = qa_mul(ab, b);
    ASSERT(qa_eq(check, a));
    qa_free(a); qa_free(b); qa_free(ab); qa_free(check);
    qaext_free(ext);
}

/* ================================ Driver ================================ */

int main(void) {
    printf("Running QA (Q(α)) tests...\n");

    TEST(test_qaext_sqrt_si);
    TEST(test_qaext_root_si);

    TEST(test_qa_zero_one_alpha);
    TEST(test_qa_from_si);
    TEST(test_qa_copy_eq);

    TEST(test_qa_add_sub_neg);
    TEST(test_qa_scale);

    TEST(test_qa_mul_sqrt2_squared);
    TEST(test_qa_mul_distributive);
    TEST(test_qa_mul_cube_root);
    TEST(test_qa_mul_imaginary);

    TEST(test_qa_inverse_sqrt2);
    TEST(test_qa_inverse_general);
    TEST(test_qa_inverse_cube_root);
    TEST(test_qa_inverse_zero_is_null);

    TEST(test_qa_div);

    TEST(test_qa_pythagoras_in_qsqrt2);
    TEST(test_qa_division_round_trip);

    printf("All QA tests passed!\n");
    return 0;
}
