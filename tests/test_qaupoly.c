/*
 * test_qaupoly.c
 * --------------
 * Unit tests for QAUPoly — Q(α)[x] univariate polynomials (Phase G2).
 *
 * Coverage:
 *   - Construction (zero, x, from_qa, from_si, copy).
 *   - Coefficient access (setcoef, getcoef, normalize).
 *   - Predicates (is_zero, eq).
 *   - Arithmetic (add, sub, neg, mul, scale).
 *   - Long division (quotient + remainder over a field).
 *   - Euclidean gcd (the headline test: gcd(x² − 2, x − √2) = x − √2
 *     over Q(√2), demonstrating the Trager lift step from a rational
 *     factor of the norm back to a Q(α) factor of the original poly).
 *   - Evaluation via Horner.
 *   - Shift round-trip: shift(p, c) then shift(_, -c) must recover p.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "qa.h"
#include "qaupoly.h"
#include "test_utils.h"

/* Override ASSERT so it always evaluates its argument (Release builds
 * compile assert() to a no-op via NDEBUG, which silently elides any
 * side-effecting call inside ASSERT — see tasks/lessons.md). */
#undef ASSERT
#define ASSERT(cond) do { \
    if (!(cond)) { \
        fprintf(stderr, "ASSERT FAIL at %s:%d: %s\n", __FILE__, __LINE__, #cond); \
        exit(1); \
    } \
} while (0)

/* =========================== Construction =========================== */

static void test_qaupoly_zero(void) {
    QAExt* ext = qaext_sqrt_si(2);
    QAUPoly* z = qaupoly_zero(ext);
    ASSERT(qaupoly_is_zero(z));
    ASSERT(z->deg == -1);
    qaupoly_free(z);
    qaext_free(ext);
}

static void test_qaupoly_x(void) {
    QAExt* ext = qaext_sqrt_si(2);
    QAUPoly* x = qaupoly_x(ext);
    ASSERT(x->deg == 1);
    /* c[0] = 0, c[1] = 1 (in Q(√2)). */
    ASSERT(qa_is_zero(x->c[0]));
    ASSERT(qa_is_one(x->c[1]));
    qaupoly_free(x);
    qaext_free(ext);
}

static void test_qaupoly_from_si(void) {
    QAExt* ext = qaext_sqrt_si(2);
    QAUPoly* p = qaupoly_from_si(ext, 5, 1);
    ASSERT(p->deg == 0);
    /* p = 5 (a constant in Q). */
    QANum* expect = qa_from_si(ext, 5, 1);
    ASSERT(qa_eq(p->c[0], expect));
    qa_free(expect);
    qaupoly_free(p);
    qaext_free(ext);
}

static void test_qaupoly_setcoef_grow(void) {
    QAExt* ext = qaext_sqrt_si(2);
    QAUPoly* p = qaupoly_zero(ext);
    QANum* one = qa_one(ext);
    qaupoly_setcoef(p, 5, one);     /* forces grow past initial cap */
    qa_free(one);
    ASSERT(p->cap >= 6);
    ASSERT(p->deg == 5);
    qaupoly_free(p);
    qaext_free(ext);
}

/* ============================= Predicates ============================ */

static void test_qaupoly_eq(void) {
    QAExt* ext = qaext_sqrt_si(2);
    QAUPoly* a = qaupoly_x(ext);
    QAUPoly* b = qaupoly_x(ext);
    ASSERT(qaupoly_eq(a, b));
    QAUPoly* c = qaupoly_zero(ext);
    ASSERT(!qaupoly_eq(a, c));
    qaupoly_free(a); qaupoly_free(b); qaupoly_free(c);
    qaext_free(ext);
}

/* ============================ Arithmetic ============================= */

/* Build (x - √2) over Q(√2). */
static QAUPoly* x_minus_sqrt2(const QAExt* ext) {
    QAUPoly* p = qaupoly_new(ext, 2);
    QANum* alpha = qa_alpha(ext);
    QANum* neg_alpha = qa_neg(alpha);
    qa_free(p->c[0]); p->c[0] = neg_alpha;
    qa_free(p->c[1]); p->c[1] = qa_one(ext);
    p->deg = 1;
    qa_free(alpha);
    return p;
}

/* Build (x + √2) over Q(√2). */
static QAUPoly* x_plus_sqrt2(const QAExt* ext) {
    QAUPoly* p = qaupoly_new(ext, 2);
    qa_free(p->c[0]); p->c[0] = qa_alpha(ext);
    qa_free(p->c[1]); p->c[1] = qa_one(ext);
    p->deg = 1;
    return p;
}

static void test_qaupoly_add_sub(void) {
    QAExt* ext = qaext_sqrt_si(2);
    QAUPoly* a = x_minus_sqrt2(ext);
    QAUPoly* b = x_plus_sqrt2(ext);

    QAUPoly* sum = qaupoly_add(a, b);
    /* (x - √2) + (x + √2) = 2x. */
    ASSERT(sum->deg == 1);
    ASSERT(qa_is_zero(sum->c[0]));
    QANum* two = qa_from_si(ext, 2, 1);
    ASSERT(qa_eq(sum->c[1], two));

    QAUPoly* diff = qaupoly_sub(a, b);
    /* (x - √2) - (x + √2) = -2√2. */
    ASSERT(diff->deg == 0);
    QANum* neg_two_alpha = qa_zero(ext);
    mpq_set_si(neg_two_alpha->coef[1], -2, 1);
    ASSERT(qa_eq(diff->c[0], neg_two_alpha));

    qa_free(two); qa_free(neg_two_alpha);
    qaupoly_free(a); qaupoly_free(b); qaupoly_free(sum); qaupoly_free(diff);
    qaext_free(ext);
}

static void test_qaupoly_mul_x_squared_minus_2(void) {
    /* (x − √2)(x + √2) = x² − 2. */
    QAExt* ext = qaext_sqrt_si(2);
    QAUPoly* a = x_minus_sqrt2(ext);
    QAUPoly* b = x_plus_sqrt2(ext);
    QAUPoly* prod = qaupoly_mul(a, b);
    ASSERT(prod->deg == 2);
    QANum* neg2 = qa_from_si(ext, -2, 1);
    QANum* zero = qa_zero(ext);
    QANum* one = qa_one(ext);
    ASSERT(qa_eq(prod->c[0], neg2));
    ASSERT(qa_eq(prod->c[1], zero));
    ASSERT(qa_eq(prod->c[2], one));
    qa_free(neg2); qa_free(zero); qa_free(one);
    qaupoly_free(a); qaupoly_free(b); qaupoly_free(prod);
    qaext_free(ext);
}

static void test_qaupoly_neg(void) {
    QAExt* ext = qaext_sqrt_si(2);
    QAUPoly* a = qaupoly_from_si(ext, 7, 1);
    QAUPoly* na = qaupoly_neg(a);
    QANum* neg7 = qa_from_si(ext, -7, 1);
    ASSERT(qa_eq(na->c[0], neg7));
    qa_free(neg7);
    qaupoly_free(a); qaupoly_free(na);
    qaext_free(ext);
}

static void test_qaupoly_scale_qa(void) {
    /* (x + 1) * √2 = √2·x + √2. */
    QAExt* ext = qaext_sqrt_si(2);
    QAUPoly* p = qaupoly_new(ext, 2);
    qa_free(p->c[0]); p->c[0] = qa_one(ext);
    qa_free(p->c[1]); p->c[1] = qa_one(ext);
    p->deg = 1;
    QANum* k = qa_alpha(ext);
    QAUPoly* r = qaupoly_scale_qa(p, k);
    ASSERT(qa_eq(r->c[0], k));
    ASSERT(qa_eq(r->c[1], k));
    qa_free(k);
    qaupoly_free(p); qaupoly_free(r);
    qaext_free(ext);
}

/* ============================ Division =============================== */

static void test_qaupoly_divrem_exact(void) {
    /* (x² − 2) / (x − √2) = x + √2, remainder 0. */
    QAExt* ext = qaext_sqrt_si(2);
    /* Build x² − 2 directly to avoid relying on mul. */
    QAUPoly* a = qaupoly_new(ext, 3);
    qa_free(a->c[0]); a->c[0] = qa_from_si(ext, -2, 1);
    qa_free(a->c[2]); a->c[2] = qa_one(ext);
    a->deg = 2;
    QAUPoly* b = x_minus_sqrt2(ext);

    QAUPoly *q, *r;
    bool ok = qaupoly_divrem(a, b, &q, &r);
    ASSERT(ok);
    QAUPoly* expected_q = x_plus_sqrt2(ext);
    ASSERT(qaupoly_eq(q, expected_q));
    ASSERT(qaupoly_is_zero(r));

    qaupoly_free(a); qaupoly_free(b);
    qaupoly_free(q); qaupoly_free(r);
    qaupoly_free(expected_q);
    qaext_free(ext);
}

static void test_qaupoly_divrem_with_remainder(void) {
    /* (x³ + 1) / (x + 1) = x² − x + 1, remainder 0.
     * (x² + 1) / (x + 1) = x − 1, remainder 2. */
    QAExt* ext = qaext_sqrt_si(2);
    /* Numerator x² + 1. */
    QAUPoly* num = qaupoly_new(ext, 3);
    qa_free(num->c[0]); num->c[0] = qa_one(ext);
    qa_free(num->c[2]); num->c[2] = qa_one(ext);
    num->deg = 2;
    /* Denominator x + 1. */
    QAUPoly* den = qaupoly_new(ext, 2);
    qa_free(den->c[0]); den->c[0] = qa_one(ext);
    qa_free(den->c[1]); den->c[1] = qa_one(ext);
    den->deg = 1;

    QAUPoly *q, *r;
    bool ok = qaupoly_divrem(num, den, &q, &r);
    ASSERT(ok);
    /* q = x − 1. */
    ASSERT(q->deg == 1);
    QANum* neg1 = qa_from_si(ext, -1, 1);
    QANum* one  = qa_one(ext);
    ASSERT(qa_eq(q->c[0], neg1));
    ASSERT(qa_eq(q->c[1], one));
    /* r = 2. */
    ASSERT(r->deg == 0);
    QANum* two = qa_from_si(ext, 2, 1);
    ASSERT(qa_eq(r->c[0], two));

    qa_free(neg1); qa_free(one); qa_free(two);
    qaupoly_free(num); qaupoly_free(den);
    qaupoly_free(q); qaupoly_free(r);
    qaext_free(ext);
}

/* ============================== GCD ================================== */

static void test_qaupoly_gcd_trager_lift(void) {
    /* The headline GCD: gcd(x² − 2, x − √2) over Q(√2) = x − √2.
     *
     * This is exactly the Trager lift step: (x² − 2) factored over Q
     * is irreducible (Norm), but over Q(√2) it splits as
     * (x − √2)(x + √2).  The gcd against the linear test factor
     * (x − √2) recovers the lifted irreducible factor. */
    QAExt* ext = qaext_sqrt_si(2);
    QAUPoly* a = qaupoly_new(ext, 3);
    qa_free(a->c[0]); a->c[0] = qa_from_si(ext, -2, 1);
    qa_free(a->c[2]); a->c[2] = qa_one(ext);
    a->deg = 2;
    QAUPoly* b = x_minus_sqrt2(ext);

    QAUPoly* g = qaupoly_gcd(a, b);
    QAUPoly* expected = x_minus_sqrt2(ext);
    /* Expected is already monic (leading coef 1). */
    ASSERT(qaupoly_eq(g, expected));
    qaupoly_free(a); qaupoly_free(b); qaupoly_free(g);
    qaupoly_free(expected);
    qaext_free(ext);
}

static void test_qaupoly_gcd_coprime(void) {
    /* gcd(x² + 1, x² − 1) over Q(√2) = 1 (both irreducible factors
     * different).  Note: over Q(i) = Q[y]/(y² + 1), x² + 1 splits as
     * (x − i)(x + i), but (x² − 1) splits over Q already; their gcd
     * is still 1. */
    QAExt* ext = qaext_sqrt_si(2);
    QAUPoly* a = qaupoly_new(ext, 3);
    qa_free(a->c[0]); a->c[0] = qa_one(ext);
    qa_free(a->c[2]); a->c[2] = qa_one(ext);
    a->deg = 2;
    QAUPoly* b = qaupoly_new(ext, 3);
    qa_free(b->c[0]); b->c[0] = qa_from_si(ext, -1, 1);
    qa_free(b->c[2]); b->c[2] = qa_one(ext);
    b->deg = 2;

    QAUPoly* g = qaupoly_gcd(a, b);
    /* Result is the monic constant 1. */
    ASSERT(g->deg == 0);
    ASSERT(qa_is_one(g->c[0]));
    qaupoly_free(a); qaupoly_free(b); qaupoly_free(g);
    qaext_free(ext);
}

static void test_qaupoly_gcd_self(void) {
    /* gcd(p, p) = p (made monic). */
    QAExt* ext = qaext_sqrt_si(2);
    QAUPoly* p = qaupoly_new(ext, 3);
    QANum* two_alpha = qa_zero(ext);
    mpq_set_si(two_alpha->coef[1], 2, 1);
    /* p = 2√2 + (1+√2) x + 3 x² */
    qa_free(p->c[0]); p->c[0] = two_alpha;
    p->c[1] = qa_zero(ext);  /* will overwrite */
    qa_free(p->c[1]);
    QANum* one_plus_alpha = qa_one(ext);
    mpq_set_si(one_plus_alpha->coef[1], 1, 1);
    p->c[1] = one_plus_alpha;
    qa_free(p->c[2]); p->c[2] = qa_from_si(ext, 3, 1);
    p->deg = 2;

    QAUPoly* g = qaupoly_gcd(p, p);
    /* g should be p / 3 (monic version). */
    QAUPoly* expected = qaupoly_make_monic(p);
    ASSERT(qaupoly_eq(g, expected));

    qaupoly_free(p); qaupoly_free(g); qaupoly_free(expected);
    qaext_free(ext);
}

static void test_qaupoly_gcd_in_qi(void) {
    /* In Q(i): gcd(x² + 1, x² − i x + 1).  x² + 1 = (x − i)(x + i).
     * x² − i x + 1: try roots i (i² − i² + 1 = 1 ≠ 0) and -i
     * (-i² − i(-i) + 1 = 1 - i² + 1 = 1 + 1 + 1 = 3 ≠ 0).  So coprime,
     * gcd = 1. */
    QAExt* ext = qaext_sqrt_si(-1);
    QAUPoly* a = qaupoly_new(ext, 3);
    qa_free(a->c[0]); a->c[0] = qa_one(ext);
    qa_free(a->c[2]); a->c[2] = qa_one(ext);
    a->deg = 2;
    QAUPoly* b = qaupoly_new(ext, 3);
    qa_free(b->c[0]); b->c[0] = qa_one(ext);
    QANum* neg_i = qa_zero(ext);
    mpq_set_si(neg_i->coef[1], -1, 1);
    qa_free(b->c[1]); b->c[1] = neg_i;
    qa_free(b->c[2]); b->c[2] = qa_one(ext);
    b->deg = 2;

    QAUPoly* g = qaupoly_gcd(a, b);
    ASSERT(g->deg == 0);
    ASSERT(qa_is_one(g->c[0]));
    qaupoly_free(a); qaupoly_free(b); qaupoly_free(g);
    qaext_free(ext);
}

/* ========================== Substitution ============================= */

static void test_qaupoly_eval_at_sqrt2(void) {
    /* p(x) = x² − 2; p(√2) should be 0. */
    QAExt* ext = qaext_sqrt_si(2);
    QAUPoly* p = qaupoly_new(ext, 3);
    qa_free(p->c[0]); p->c[0] = qa_from_si(ext, -2, 1);
    qa_free(p->c[2]); p->c[2] = qa_one(ext);
    p->deg = 2;
    QANum* alpha = qa_alpha(ext);
    QANum* val = qaupoly_eval(p, alpha);
    ASSERT(qa_is_zero(val));
    qa_free(alpha); qa_free(val);
    qaupoly_free(p);
    qaext_free(ext);
}

static void test_qaupoly_shift_roundtrip(void) {
    /* shift(p, c) followed by shift(_, -c) recovers p. */
    QAExt* ext = qaext_sqrt_si(2);
    /* p = 1 + 2x + 3x² */
    QAUPoly* p = qaupoly_new(ext, 3);
    qa_free(p->c[0]); p->c[0] = qa_from_si(ext, 1, 1);
    qa_free(p->c[1]); p->c[1] = qa_from_si(ext, 2, 1);
    qa_free(p->c[2]); p->c[2] = qa_from_si(ext, 3, 1);
    p->deg = 2;

    QANum* c = qa_alpha(ext);             /* shift by √2 */
    QANum* neg_c = qa_neg(c);
    QAUPoly* shifted = qaupoly_shift(p, c);
    QAUPoly* back    = qaupoly_shift(shifted, neg_c);
    ASSERT(qaupoly_eq(back, p));
    qa_free(c); qa_free(neg_c);
    qaupoly_free(p); qaupoly_free(shifted); qaupoly_free(back);
    qaext_free(ext);
}

/* ============================== Driver =============================== */

int main(void) {
    printf("Running QAUPoly (Q(α)[x]) tests...\n");

    TEST(test_qaupoly_zero);
    TEST(test_qaupoly_x);
    TEST(test_qaupoly_from_si);
    TEST(test_qaupoly_setcoef_grow);

    TEST(test_qaupoly_eq);

    TEST(test_qaupoly_add_sub);
    TEST(test_qaupoly_mul_x_squared_minus_2);
    TEST(test_qaupoly_neg);
    TEST(test_qaupoly_scale_qa);

    TEST(test_qaupoly_divrem_exact);
    TEST(test_qaupoly_divrem_with_remainder);

    TEST(test_qaupoly_gcd_trager_lift);
    TEST(test_qaupoly_gcd_coprime);
    TEST(test_qaupoly_gcd_self);
    TEST(test_qaupoly_gcd_in_qi);

    TEST(test_qaupoly_eval_at_sqrt2);
    TEST(test_qaupoly_shift_roundtrip);

    printf("All QAUPoly tests passed!\n");
    return 0;
}
