/*
 * test_zupoly.c
 * -------------
 * Unit tests for the ZUPoly primitive operations.  Each test exercises
 * one operation in isolation so that failures localise immediately.
 *
 * Coverage:
 *   - Construction: zupoly_zero, zupoly_from_int, zupoly_copy
 *   - Coefficient access: setcoef, getcoef, normalize
 *   - Predicates: is_zero, eq, cmp
 *   - Arithmetic: add, sub, mul, neg, scale
 *   - Division: divrem_monic, pseudodivrem, divexact
 *   - GCD, content, primitive_part
 *   - Substitution: eval_si, shift_si
 *   - Conversions: expr_to_zupoly, zupoly_to_expr
 *
 * The tests build polynomials by hand using the primitive setters and
 * compare results coefficient-by-coefficient (so the test does NOT
 * depend on the conversion-to-Expr code, which is itself tested).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include "zupoly.h"
#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "test_utils.h"

/* ====================================================================== */
/*  Helpers                                                               */
/* ====================================================================== */

/* Construct a ZUPoly from a length-N int64 array of coefficients (low
 * to high).  Useful for compact test fixtures. */
static ZUPoly* poly_of(const int64_t* coefs, int n) {
    ZUPoly* p = zupoly_new(n);
    for (int i = 0; i < n; i++) zupoly_setcoef_si(p, i, coefs[i]);
    return p;
}

/* Compare a ZUPoly against an int64 coefficient array. */
static int poly_eq_array(const ZUPoly* p, const int64_t* coefs, int n) {
    /* Trim trailing zeros in the array to find expected degree. */
    int expected_deg = n - 1;
    while (expected_deg >= 0 && coefs[expected_deg] == 0) expected_deg--;
    if (p->deg != expected_deg) return 0;
    for (int i = 0; i <= expected_deg; i++) {
        const mpz_t* ci = zupoly_getcoef(p, i);
        if (!ci) return 0;
        if (mpz_cmp_si(*ci, (long)coefs[i]) != 0) return 0;
    }
    return 1;
}

static void fail_dump(const char* label, const ZUPoly* p,
                      const int64_t* expected, int n) {
    fprintf(stderr, "FAIL: %s\n  expected: [", label);
    for (int i = 0; i < n; i++) fprintf(stderr, "%lld%s", (long long)expected[i], i+1<n?", ":"");
    fprintf(stderr, "]\n  got:      ");
    zupoly_print(p, "x");
    fprintf(stderr, "  (deg=%d)\n", p->deg);
}

/* The default test_utils.h ASSERT macro routes through <assert.h>, which
 * is compiled out under -DNDEBUG (CMake's Release build).  We need real
 * runtime checks here.  Override with a hard-failing variant. */
#undef ASSERT
#define ASSERT(cond) do { \
    if (!(cond)) { \
        fprintf(stderr, "ASSERTION FAILED: %s\n  at %s:%d\n", \
                #cond, __FILE__, __LINE__); \
        abort(); \
    } \
} while (0)


#define EXPECT_POLY(p, ...) do { \
    int64_t _exp[] = {__VA_ARGS__}; \
    int _n = (int)(sizeof(_exp) / sizeof(_exp[0])); \
    if (!poly_eq_array((p), _exp, _n)) { \
        fail_dump(#p " == [" #__VA_ARGS__ "]", (p), _exp, _n); \
        ASSERT(0); \
    } \
} while (0)

/* ====================================================================== */
/*  Construction & coefficient access                                     */
/* ====================================================================== */

static void test_zero_construction(void) {
    ZUPoly* p = zupoly_zero();
    ASSERT(zupoly_is_zero(p));
    ASSERT(p->deg == -1);
    zupoly_free(p);
}

static void test_from_int(void) {
    ZUPoly* p = zupoly_from_int(42);
    ASSERT(!zupoly_is_zero(p));
    ASSERT(p->deg == 0);
    EXPECT_POLY(p, 42);
    zupoly_free(p);

    /* Zero constant. */
    ZUPoly* z = zupoly_from_int(0);
    ASSERT(zupoly_is_zero(z));
    zupoly_free(z);
}

static void test_setcoef_growth(void) {
    /* Build x^5 + 3 by setting coefficients out of order and ensure
     * capacity grows correctly. */
    ZUPoly* p = zupoly_zero();
    zupoly_setcoef_si(p, 5, 1);
    zupoly_setcoef_si(p, 0, 3);
    EXPECT_POLY(p, 3, 0, 0, 0, 0, 1);
    ASSERT(p->deg == 5);
    zupoly_free(p);
}

static void test_setcoef_clear_leading(void) {
    /* Setting the leading coefficient to zero must update deg
     * downward to the next non-zero coefficient. */
    ZUPoly* p = poly_of((int64_t[]){1, 2, 3, 4}, 4);
    ASSERT(p->deg == 3);
    mpz_t zero; mpz_init(zero);
    zupoly_setcoef(p, 3, zero);
    ASSERT(p->deg == 2);
    EXPECT_POLY(p, 1, 2, 3);
    /* Also clearing inner coefficients must leave deg unchanged. */
    zupoly_setcoef(p, 1, zero);
    ASSERT(p->deg == 2);
    EXPECT_POLY(p, 1, 0, 3);
    mpz_clear(zero);
    zupoly_free(p);
}

static void test_copy_independence(void) {
    /* Mutating a copy must not affect the original. */
    ZUPoly* a = poly_of((int64_t[]){1, 2, 3}, 3);
    ZUPoly* b = zupoly_copy(a);
    zupoly_setcoef_si(b, 0, 99);
    EXPECT_POLY(a, 1, 2, 3);
    EXPECT_POLY(b, 99, 2, 3);
    zupoly_free(a); zupoly_free(b);
}

/* ====================================================================== */
/*  Predicates                                                            */
/* ====================================================================== */

static void test_eq(void) {
    ZUPoly* a = poly_of((int64_t[]){1, 2, 3}, 3);
    ZUPoly* b = poly_of((int64_t[]){1, 2, 3}, 3);
    ZUPoly* c = poly_of((int64_t[]){1, 2, 4}, 3);
    ZUPoly* d = poly_of((int64_t[]){1, 2}, 2);
    ASSERT(zupoly_eq(a, b));
    ASSERT(!zupoly_eq(a, c));
    ASSERT(!zupoly_eq(a, d));
    zupoly_free(a); zupoly_free(b); zupoly_free(c); zupoly_free(d);
}

static void test_cmp(void) {
    ZUPoly* a = poly_of((int64_t[]){1, 2, 3}, 3);
    ZUPoly* b = poly_of((int64_t[]){1, 2}, 2);   /* lower deg */
    ZUPoly* c = poly_of((int64_t[]){1, 9, 3}, 3);/* higher leading-second */
    ASSERT(zupoly_cmp(a, b) > 0);
    ASSERT(zupoly_cmp(b, a) < 0);
    ASSERT(zupoly_cmp(a, c) < 0); /* same deg, c has 9>2 in middle */
    ASSERT(zupoly_cmp(a, a) == 0);
    zupoly_free(a); zupoly_free(b); zupoly_free(c);
}

/* ====================================================================== */
/*  Arithmetic                                                            */
/* ====================================================================== */

static void test_add(void) {
    /* (1 + 2x + 3x^2) + (4 + 5x) = 5 + 7x + 3x^2 */
    ZUPoly* a = poly_of((int64_t[]){1, 2, 3}, 3);
    ZUPoly* b = poly_of((int64_t[]){4, 5}, 2);
    ZUPoly* r = zupoly_add(a, b);
    EXPECT_POLY(r, 5, 7, 3);
    zupoly_free(a); zupoly_free(b); zupoly_free(r);
}

static void test_add_with_cancellation(void) {
    /* Highest coefficient cancels to zero -- degree must drop. */
    ZUPoly* a = poly_of((int64_t[]){1, 2, 3}, 3);
    ZUPoly* b = poly_of((int64_t[]){0, 0, -3}, 3);
    ZUPoly* r = zupoly_add(a, b);
    ASSERT(r->deg == 1);
    EXPECT_POLY(r, 1, 2);
    zupoly_free(a); zupoly_free(b); zupoly_free(r);
}

static void test_sub(void) {
    ZUPoly* a = poly_of((int64_t[]){5, 7, 3}, 3);
    ZUPoly* b = poly_of((int64_t[]){4, 5}, 2);
    ZUPoly* r = zupoly_sub(a, b);
    EXPECT_POLY(r, 1, 2, 3);
    zupoly_free(a); zupoly_free(b); zupoly_free(r);
}

static void test_mul(void) {
    /* (x + 1)(x - 1) = x^2 - 1 */
    ZUPoly* a = poly_of((int64_t[]){1, 1}, 2);
    ZUPoly* b = poly_of((int64_t[]){-1, 1}, 2);
    ZUPoly* r = zupoly_mul(a, b);
    EXPECT_POLY(r, -1, 0, 1);
    zupoly_free(a); zupoly_free(b); zupoly_free(r);

    /* (x^2 + x + 1)(x - 1) = x^3 - 1 */
    ZUPoly* p = poly_of((int64_t[]){1, 1, 1}, 3);
    ZUPoly* q = poly_of((int64_t[]){-1, 1}, 2);
    ZUPoly* pq = zupoly_mul(p, q);
    EXPECT_POLY(pq, -1, 0, 0, 1);
    zupoly_free(p); zupoly_free(q); zupoly_free(pq);
}

static void test_mul_zero(void) {
    ZUPoly* a = poly_of((int64_t[]){1, 2, 3}, 3);
    ZUPoly* z = zupoly_zero();
    ZUPoly* r = zupoly_mul(a, z);
    ASSERT(zupoly_is_zero(r));
    zupoly_free(a); zupoly_free(z); zupoly_free(r);
}

static void test_neg(void) {
    ZUPoly* a = poly_of((int64_t[]){1, -2, 3}, 3);
    ZUPoly* r = zupoly_neg(a);
    EXPECT_POLY(r, -1, 2, -3);
    zupoly_free(a); zupoly_free(r);
}

static void test_scale(void) {
    ZUPoly* a = poly_of((int64_t[]){1, 2, 3}, 3);
    ZUPoly* r = zupoly_scale_si(a, 5);
    EXPECT_POLY(r, 5, 10, 15);
    zupoly_free(a); zupoly_free(r);
}

/* ====================================================================== */
/*  Division                                                              */
/* ====================================================================== */

static void test_divrem_monic(void) {
    /* Divide x^2 - 1 by x - 1.  Quotient = x + 1, remainder = 0. */
    ZUPoly* a = poly_of((int64_t[]){-1, 0, 1}, 3);
    ZUPoly* b = poly_of((int64_t[]){-1, 1}, 2);
    ZUPoly *q, *r;
    ASSERT(zupoly_divrem_monic(a, b, &q, &r));
    EXPECT_POLY(q, 1, 1);
    ASSERT(zupoly_is_zero(r));
    zupoly_free(a); zupoly_free(b); zupoly_free(q); zupoly_free(r);
}

static void test_divrem_remainder(void) {
    /* Divide x^2 + 1 by x - 1.  Quotient = x + 1, remainder = 2. */
    ZUPoly* a = poly_of((int64_t[]){1, 0, 1}, 3);
    ZUPoly* b = poly_of((int64_t[]){-1, 1}, 2);
    ZUPoly *q, *r;
    ASSERT(zupoly_divrem_monic(a, b, &q, &r));
    EXPECT_POLY(q, 1, 1);
    EXPECT_POLY(r, 2);
    zupoly_free(a); zupoly_free(b); zupoly_free(q); zupoly_free(r);
}

static void test_divrem_failure_when_inexact_in_Z(void) {
    /* Note: zupoly_divrem_monic uses long division, where the only
     * required exactness check is that lc(b) divides the running
     * remainder's leading coefficient at each step.  For (1 + 2x) / (2x)
     * we have lc(b) = 2 and lc(a) = 2; division proceeds and produces
     * (q=1, r=1).  divrem_monic correctly returns SUCCESS with a non-
     * zero remainder -- exactness at the algebraic level is the
     * caller's responsibility (i.e., zupoly_divexact).  This test
     * verifies that contract: divrem_monic returns true here, with
     * remainder 1 != 0; it is NOT an exactness failure. */
    ZUPoly* a = poly_of((int64_t[]){1, 2}, 2);   /* 1 + 2x */
    ZUPoly* b = poly_of((int64_t[]){0, 2}, 2);   /* 2x */
    ZUPoly *q = NULL, *r = NULL;
    bool ok = zupoly_divrem_monic(a, b, &q, &r);
    ASSERT(ok);
    /* Quotient = 1, remainder = 1. */
    EXPECT_POLY(q, 1);
    EXPECT_POLY(r, 1);
    zupoly_free(a); zupoly_free(b); zupoly_free(q); zupoly_free(r);

    /* For real exactness, divexact must return NULL. */
    a = poly_of((int64_t[]){1, 2}, 2);
    b = poly_of((int64_t[]){0, 2}, 2);
    ZUPoly* div = zupoly_divexact(a, b);
    ASSERT(div == NULL);
    zupoly_free(a); zupoly_free(b);
}

static void test_pseudodivrem(void) {
    /* Pseudo-divide x^2 by 2x.  Recipe: lc(2x)^(2-1+1) * x^2 = 4 x^2;
     * 4 x^2 / 2x = 2 x with remainder 0. */
    ZUPoly* a = poly_of((int64_t[]){0, 0, 1}, 3);
    ZUPoly* b = poly_of((int64_t[]){0, 2}, 2);
    ZUPoly *q, *r;
    ASSERT(zupoly_pseudodivrem(a, b, &q, &r));
    EXPECT_POLY(q, 0, 2);
    ASSERT(zupoly_is_zero(r));
    zupoly_free(a); zupoly_free(b); zupoly_free(q); zupoly_free(r);
}

static void test_divexact_success(void) {
    /* (x + 2)(x + 3) = x^2 + 5x + 6.  Dividing by (x + 2) -> (x + 3). */
    ZUPoly* a = poly_of((int64_t[]){6, 5, 1}, 3);
    ZUPoly* b = poly_of((int64_t[]){2, 1}, 2);
    ZUPoly* q = zupoly_divexact(a, b);
    ASSERT(q != NULL);
    EXPECT_POLY(q, 3, 1);
    zupoly_free(a); zupoly_free(b); zupoly_free(q);
}

static void test_divexact_failure(void) {
    /* (x^2 + 1) is not divisible by (x - 1) over Z. */
    ZUPoly* a = poly_of((int64_t[]){1, 0, 1}, 3);
    ZUPoly* b = poly_of((int64_t[]){-1, 1}, 2);
    ZUPoly* q = zupoly_divexact(a, b);
    ASSERT(q == NULL);
    zupoly_free(a); zupoly_free(b);
}

/* ====================================================================== */
/*  GCD, content, primitive part                                          */
/* ====================================================================== */

static void test_content(void) {
    /* content(6 + 9x + 12x^2) = 3. */
    ZUPoly* p = poly_of((int64_t[]){6, 9, 12}, 3);
    mpz_t c; mpz_init(c);
    zupoly_content(p, c);
    ASSERT(mpz_cmp_ui(c, 3) == 0);
    mpz_clear(c);
    zupoly_free(p);
}

static void test_primitive_part(void) {
    /* primitive(6 + 9x + 12x^2) = 2 + 3x + 4x^2. */
    ZUPoly* p = poly_of((int64_t[]){6, 9, 12}, 3);
    ZUPoly* pp = zupoly_primitive_part(p);
    EXPECT_POLY(pp, 2, 3, 4);
    zupoly_free(p); zupoly_free(pp);
}

static void test_gcd_simple(void) {
    /* gcd((x+1)(x-1), (x+1)(x+2)) = x+1. */
    ZUPoly* a = poly_of((int64_t[]){-1, 0, 1}, 3);    /* x^2 - 1 */
    ZUPoly* b = poly_of((int64_t[]){2, 3, 1}, 3);     /* x^2 + 3x + 2 */
    ZUPoly* g = zupoly_gcd(a, b);
    EXPECT_POLY(g, 1, 1);
    zupoly_free(a); zupoly_free(b); zupoly_free(g);
}

static void test_gcd_with_content(void) {
    /* gcd(6(x+1), 4(x+1)) = 2(x+1) */
    ZUPoly* a = poly_of((int64_t[]){6, 6}, 2);
    ZUPoly* b = poly_of((int64_t[]){4, 4}, 2);
    ZUPoly* g = zupoly_gcd(a, b);
    EXPECT_POLY(g, 2, 2);
    zupoly_free(a); zupoly_free(b); zupoly_free(g);
}

static void test_gcd_coprime(void) {
    /* gcd(x^2+1, x+1) = 1 (over Z). */
    ZUPoly* a = poly_of((int64_t[]){1, 0, 1}, 3);
    ZUPoly* b = poly_of((int64_t[]){1, 1}, 2);
    ZUPoly* g = zupoly_gcd(a, b);
    EXPECT_POLY(g, 1);
    zupoly_free(a); zupoly_free(b); zupoly_free(g);
}

/* ====================================================================== */
/*  Substitution / evaluation                                             */
/* ====================================================================== */

static void test_eval(void) {
    /* p(x) = 1 + 2x + 3x^2; p(2) = 1 + 4 + 12 = 17. */
    ZUPoly* p = poly_of((int64_t[]){1, 2, 3}, 3);
    mpz_t out; mpz_init(out);
    zupoly_eval_si(p, 2, out);
    ASSERT(mpz_cmp_ui(out, 17) == 0);

    /* p(0) = 1, p(1) = 6, p(-1) = 2. */
    zupoly_eval_si(p, 0, out); ASSERT(mpz_cmp_ui(out, 1) == 0);
    zupoly_eval_si(p, 1, out); ASSERT(mpz_cmp_ui(out, 6) == 0);
    zupoly_eval_si(p, -1, out); ASSERT(mpz_cmp_ui(out, 2) == 0);

    mpz_clear(out);
    zupoly_free(p);
}

static void test_shift(void) {
    /* p(x) = x^2.  p(x + 1) = x^2 + 2x + 1. */
    ZUPoly* p = poly_of((int64_t[]){0, 0, 1}, 3);
    ZUPoly* s = zupoly_shift_si(p, 1);
    EXPECT_POLY(s, 1, 2, 1);
    zupoly_free(p); zupoly_free(s);

    /* p(x) = x^3 - x.  p(x - 1) = (x-1)^3 - (x-1).
     *   = x^3 - 3x^2 + 3x - 1 - x + 1 = x^3 - 3x^2 + 2x */
    ZUPoly* q = poly_of((int64_t[]){0, -1, 0, 1}, 4);
    ZUPoly* qs = zupoly_shift_si(q, -1);
    EXPECT_POLY(qs, 0, 2, -3, 1);
    zupoly_free(q); zupoly_free(qs);
}

/* ====================================================================== */
/*  Conversion to / from Expr                                             */
/* ====================================================================== */

static void test_expr_roundtrip(void) {
    /* Parse "1 + 2 x + 3 x^2", convert to ZUPoly, back to Expr,
     * and verify the coefficients match. */
    Expr* var = parse_expression("x");
    Expr* e = parse_expression("1 + 2 x + 3 x^2");
    Expr* ev = evaluate(e);
    expr_free(e);

    ZUPoly* p = expr_to_zupoly(ev, var);
    ASSERT(p != NULL);
    EXPECT_POLY(p, 1, 2, 3);

    Expr* back = zupoly_to_expr(p, var);
    /* Compare full forms by re-parsing the canonical input. */
    Expr* expected = parse_expression("1 + 2 x + 3 x^2");
    Expr* expected_ev = evaluate(expected);
    expr_free(expected);
    ASSERT(expr_eq(back, expected_ev));

    expr_free(back); expr_free(expected_ev);
    expr_free(ev); expr_free(var); zupoly_free(p);
}

static void test_expr_zero(void) {
    Expr* var = parse_expression("x");
    Expr* e = parse_expression("0");
    Expr* ev = evaluate(e);
    expr_free(e);
    ZUPoly* p = expr_to_zupoly(ev, var);
    ASSERT(p != NULL);
    ASSERT(zupoly_is_zero(p));
    zupoly_free(p); expr_free(ev); expr_free(var);
}

static void test_expr_non_integer_returns_null(void) {
    /* x + 1/2 has a rational coefficient; conversion must fail. */
    Expr* var = parse_expression("x");
    Expr* e = parse_expression("x + 1/2");
    Expr* ev = evaluate(e);
    expr_free(e);
    ZUPoly* p = expr_to_zupoly(ev, var);
    ASSERT(p == NULL);
    expr_free(ev); expr_free(var);
}

/* ====================================================================== */
/*  Compositional sanity (multiple ops together)                          */
/* ====================================================================== */

static void test_factor_then_multiply_back(void) {
    /* Take p = (x-2)(x+3)(x-5) = x^3 - 4x^2 - 11x + 30, divide by
     * (x-2) and (x+3), confirm each step is exact, multiply factors
     * back to get p. */
    ZUPoly* p = poly_of((int64_t[]){30, -11, -4, 1}, 4);
    ZUPoly* f1 = poly_of((int64_t[]){-2, 1}, 2);   /* x - 2 */
    ZUPoly* f2 = poly_of((int64_t[]){3, 1}, 2);    /* x + 3 */
    ZUPoly* f3 = poly_of((int64_t[]){-5, 1}, 2);   /* x - 5 */

    ZUPoly* q1 = zupoly_divexact(p, f1);
    ASSERT(q1 != NULL);
    ZUPoly* q2 = zupoly_divexact(q1, f2);
    ASSERT(q2 != NULL);
    ASSERT(zupoly_eq(q2, f3));

    /* Multiply back. */
    ZUPoly* m1 = zupoly_mul(f1, f2);
    ZUPoly* m2 = zupoly_mul(m1, f3);
    ASSERT(zupoly_eq(m2, p));

    zupoly_free(p); zupoly_free(f1); zupoly_free(f2); zupoly_free(f3);
    zupoly_free(q1); zupoly_free(q2); zupoly_free(m1); zupoly_free(m2);
}

static void test_subresultant_gcd_squarefree_check(void) {
    /* For squarefree p, gcd(p, p') = 1.  Test: p = x^3 - x.
     * p' = 3x^2 - 1.  Both have gcd 1 over Z. */
    ZUPoly* p = poly_of((int64_t[]){0, -1, 0, 1}, 4);
    ZUPoly* dp = poly_of((int64_t[]){-1, 0, 3}, 3);
    ZUPoly* g = zupoly_gcd(p, dp);
    EXPECT_POLY(g, 1);
    zupoly_free(p); zupoly_free(dp); zupoly_free(g);
}

static void test_non_squarefree_gcd(void) {
    /* For p = (x-1)^2 = x^2 - 2x + 1, p' = 2x - 2.
     * gcd(p, p') = x - 1 (associate).  Subresultant returns +(x-1) up to sign. */
    ZUPoly* p = poly_of((int64_t[]){1, -2, 1}, 3);
    ZUPoly* dp = poly_of((int64_t[]){-2, 2}, 2);
    ZUPoly* g = zupoly_gcd(p, dp);
    /* Accept x - 1 or 1 - x (sign-canonicalised should give x - 1). */
    EXPECT_POLY(g, -1, 1);
    zupoly_free(p); zupoly_free(dp); zupoly_free(g);
}

/* ====================================================================== */
/*  Diophantine equation                                                  */
/* ====================================================================== */

/* For monic coprime u(x), v(x) and arbitrary e(x), we expect
 *   delta_u * v + delta_v * u = e
 * with deg(delta_u) < deg(u).  These tests verify the contract. */

/* Verify the Diophantine identity by re-multiplying. */
static int verify_diophantine(const ZUPoly* u, const ZUPoly* v, const ZUPoly* e,
                              const ZUPoly* du, const ZUPoly* dv) {
    ZUPoly* duv = zupoly_mul(du, v);
    ZUPoly* dvu = zupoly_mul(dv, u);
    ZUPoly* sum = zupoly_add(duv, dvu);
    int ok = zupoly_eq(sum, e);
    zupoly_free(duv); zupoly_free(dvu); zupoly_free(sum);
    return ok;
}

static void test_diophantine_simple_monic_coprime(void) {
    /* u = x - 1, v = x - 2 (monic, coprime).  e = 1.
     * Expected: delta_u, delta_v ∈ Q[x], degrees < 1.
     * Concretely: 1 = (delta_u)*(x-2) + (delta_v)*(x-1).
     * Setting x = 1: 1 = delta_u(1)*(-1), so delta_u(1) = -1.
     *               Since deg<1, delta_u = -1.  Then -1*(x-2) +
     *               delta_v*(x-1) = 1 => delta_v = (1 + (x-2))/(x-1) =
     *               (x-1)/(x-1) = 1.  So (delta_u, delta_v) = (-1, 1).
     * Both are integer, so the function should succeed. */
    ZUPoly* u = poly_of((int64_t[]){-1, 1}, 2);
    ZUPoly* v = poly_of((int64_t[]){-2, 1}, 2);
    ZUPoly* e = zupoly_from_int(1);
    ZUPoly *du = NULL, *dv = NULL;
    bool ok = zupoly_diophantine(u, v, e, &du, &dv);
    ASSERT(ok);
    ASSERT(verify_diophantine(u, v, e, du, dv));
    /* Check the specific values. */
    EXPECT_POLY(du, -1);
    EXPECT_POLY(dv, 1);
    zupoly_free(u); zupoly_free(v); zupoly_free(e);
    zupoly_free(du); zupoly_free(dv);
}

static void test_diophantine_with_polynomial_e(void) {
    /* u = x, v = x + 1.  e = x^2 + x.
     * Solution: e = x*(x+1).
     * Bezout: 1 = -1*(x) + 1*(x+1).
     * delta_u = T*e mod u = 1 * (x^2 + x) mod x = 0.
     * delta_v = (e - delta_u*v)/u = e/u = e/x = x + 1.
     * So (delta_u, delta_v) = (0, x+1) -- but delta_u = 0 has deg < 1
     * which is fine.  Verify by re-multiplying:
     *   0*(x+1) + (x+1)*x = x^2 + x = e. */
    ZUPoly* u = poly_of((int64_t[]){0, 1}, 2);    /* x */
    ZUPoly* v = poly_of((int64_t[]){1, 1}, 2);    /* x + 1 */
    ZUPoly* e = poly_of((int64_t[]){0, 1, 1}, 3); /* x + x^2 */
    ZUPoly *du = NULL, *dv = NULL;
    bool ok = zupoly_diophantine(u, v, e, &du, &dv);
    ASSERT(ok);
    ASSERT(verify_diophantine(u, v, e, du, dv));
    zupoly_free(u); zupoly_free(v); zupoly_free(e);
    zupoly_free(du); zupoly_free(dv);
}

static void test_diophantine_higher_degree(void) {
    /* u = x^2 - 1 (monic, factors as (x-1)(x+1)),
     * v = x^2 + 1 (monic, irreducible over Z).
     * gcd over Q[x]: 1 (they share no common root).
     * e = x^3.
     * Expect a Diophantine solution with delta_u of degree < 2,
     * delta_v of degree < 2. */
    ZUPoly* u = poly_of((int64_t[]){-1, 0, 1}, 3);
    ZUPoly* v = poly_of((int64_t[]){1, 0, 1}, 3);
    ZUPoly* e = poly_of((int64_t[]){0, 0, 0, 1}, 4);
    ZUPoly *du = NULL, *dv = NULL;
    bool ok = zupoly_diophantine(u, v, e, &du, &dv);
    /* The solution involves rationals (1/2 factors), so the function
     * should return false because the result isn't in Z[x]. */
    if (!ok) {
        ASSERT(du == NULL && dv == NULL);
    } else {
        /* If it succeeds, the math must be consistent. */
        ASSERT(verify_diophantine(u, v, e, du, dv));
        ASSERT(du->deg < u->deg);
        zupoly_free(du); zupoly_free(dv);
    }
    zupoly_free(u); zupoly_free(v); zupoly_free(e);
}

static void test_diophantine_non_coprime_fails(void) {
    /* u = x, v = x (gcd = x, not 1).  Diophantine has no unique
     * solution; the function must signal failure. */
    ZUPoly* u = poly_of((int64_t[]){0, 1}, 2);
    ZUPoly* v = poly_of((int64_t[]){0, 1}, 2);
    ZUPoly* e = zupoly_from_int(1);
    ZUPoly *du = NULL, *dv = NULL;
    bool ok = zupoly_diophantine(u, v, e, &du, &dv);
    ASSERT(!ok);
    ASSERT(du == NULL && dv == NULL);
    zupoly_free(u); zupoly_free(v); zupoly_free(e);
}

static void test_diophantine_non_monic_fails(void) {
    /* u = 2x is not monic; the precondition should be enforced. */
    ZUPoly* u = poly_of((int64_t[]){0, 2}, 2);
    ZUPoly* v = poly_of((int64_t[]){1, 1}, 2);
    ZUPoly* e = zupoly_from_int(1);
    ZUPoly *du = NULL, *dv = NULL;
    bool ok = zupoly_diophantine(u, v, e, &du, &dv);
    ASSERT(!ok);
    ASSERT(du == NULL && dv == NULL);
    zupoly_free(u); zupoly_free(v); zupoly_free(e);
}

static void test_diophantine_residue_after_setup(void) {
    /* This is the form that arises in bivariate Hensel:
     *
     *   u = x - 1, v = x + 1 (both monic, coprime over Z, gcd = 1).
     *   e = 4 x.
     *   We want delta_u(x) * (x+1) + delta_v(x) * (x-1) = 4x
     *   with deg(delta_u) < 1.
     *
     *   Try delta_u = constant a, delta_v = constant b.
     *   Then a*(x+1) + b*(x-1) = (a+b)*x + (a-b) = 4x
     *   => a + b = 4, a - b = 0 => a = b = 2.
     *   So (delta_u, delta_v) = (2, 2). */
    ZUPoly* u = poly_of((int64_t[]){-1, 1}, 2);
    ZUPoly* v = poly_of((int64_t[]){1, 1}, 2);
    ZUPoly* e = poly_of((int64_t[]){0, 4}, 2);
    ZUPoly *du = NULL, *dv = NULL;
    bool ok = zupoly_diophantine(u, v, e, &du, &dv);
    ASSERT(ok);
    EXPECT_POLY(du, 2);
    EXPECT_POLY(dv, 2);
    ASSERT(verify_diophantine(u, v, e, du, dv));
    zupoly_free(u); zupoly_free(v); zupoly_free(e);
    zupoly_free(du); zupoly_free(dv);
}

/* ====================================================================== */
/*  Main                                                                  */
/* ====================================================================== */

int main(void) {
    symtab_init();
    core_init();

    printf("Running ZUPoly tests...\n");

    /* Construction & coefficient access */
    TEST(test_zero_construction);
    TEST(test_from_int);
    TEST(test_setcoef_growth);
    TEST(test_setcoef_clear_leading);
    TEST(test_copy_independence);

    /* Predicates */
    TEST(test_eq);
    TEST(test_cmp);

    /* Arithmetic */
    TEST(test_add);
    TEST(test_add_with_cancellation);
    TEST(test_sub);
    TEST(test_mul);
    TEST(test_mul_zero);
    TEST(test_neg);
    TEST(test_scale);

    /* Division */
    TEST(test_divrem_monic);
    TEST(test_divrem_remainder);
    TEST(test_divrem_failure_when_inexact_in_Z);
    TEST(test_pseudodivrem);
    TEST(test_divexact_success);
    TEST(test_divexact_failure);

    /* GCD, content, primitive part */
    TEST(test_content);
    TEST(test_primitive_part);
    TEST(test_gcd_simple);
    TEST(test_gcd_with_content);
    TEST(test_gcd_coprime);

    /* Substitution */
    TEST(test_eval);
    TEST(test_shift);

    /* Conversion */
    TEST(test_expr_roundtrip);
    TEST(test_expr_zero);
    TEST(test_expr_non_integer_returns_null);

    /* Compositional */
    TEST(test_factor_then_multiply_back);
    TEST(test_subresultant_gcd_squarefree_check);
    TEST(test_non_squarefree_gcd);

    /* Diophantine */
    TEST(test_diophantine_simple_monic_coprime);
    TEST(test_diophantine_with_polynomial_e);
    TEST(test_diophantine_higher_degree);
    TEST(test_diophantine_non_coprime_fails);
    TEST(test_diophantine_non_monic_fails);
    TEST(test_diophantine_residue_after_setup);

    printf("All ZUPoly tests passed!\n");
    return 0;
}
