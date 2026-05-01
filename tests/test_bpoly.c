/*
 * test_bpoly.c
 * ------------
 * Unit tests for the BPoly bivariate polynomial type.  Built on top of
 * the ZUPoly tests -- BPoly is parameterised by ZUPoly so its tests
 * focus on bivariate-specific operations: x-major arithmetic, mixing
 * y-coefficients, truncation modulo y^k, evaluation y -> alpha,
 * shift y -> y + alpha, exact bivariate division.
 *
 * The tests construct BPoly values directly via the primitive setters
 * to avoid coupling to the conversion-to-Expr code (which is itself
 * tested last).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include "bpoly.h"
#include "zupoly.h"
#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "test_utils.h"

#undef ASSERT
#define ASSERT(cond) do { \
    if (!(cond)) { \
        fprintf(stderr, "ASSERTION FAILED: %s\n  at %s:%d\n", \
                #cond, __FILE__, __LINE__); \
        abort(); \
    } \
} while (0)

/* ====================================================================== */
/*  Helpers                                                               */
/* ====================================================================== */

/* Build a ZUPoly from int64 coefficients (low to high). */
static ZUPoly* zu(const int64_t* c, int n) {
    ZUPoly* p = zupoly_new(n);
    for (int i = 0; i < n; i++) zupoly_setcoef_si(p, i, c[i]);
    return p;
}

/* Build a BPoly from a 2D coefficient array indexed [x_deg][y_deg].
 * For each x_deg row, builds a ZUPoly from the row's y coefficients. */
static BPoly* bp(int n_x, int n_y, const int64_t* coefs) {
    /* coefs is row-major: coefs[i * n_y + j] = c_{i,j}. */
    BPoly* p = bpoly_new(n_x);
    for (int i = 0; i < n_x; i++) {
        ZUPoly* row = zupoly_new(n_y);
        for (int j = 0; j < n_y; j++) {
            zupoly_setcoef_si(row, j, coefs[i * n_y + j]);
        }
        bpoly_set_xcoef(p, i, row);
    }
    return p;
}

/* Compare two BPolys for structural equality (degree-by-degree). */
static int bpoly_eq_array(const BPoly* p, int n_x, int n_y,
                          const int64_t* coefs) {
    for (int i = 0; i < n_x; i++) {
        const ZUPoly* row = bpoly_get_xcoef(p, i);
        for (int j = 0; j < n_y; j++) {
            int64_t expected = coefs[i * n_y + j];
            const mpz_t* c = row ? zupoly_getcoef(row, j) : NULL;
            int64_t actual = c ? mpz_get_si(*c) : 0;
            if (actual != expected) return 0;
        }
    }
    return 1;
}

#define EXPECT_BPOLY(p, n_x, n_y, ...) do { \
    int64_t _exp[] = {__VA_ARGS__}; \
    if (!bpoly_eq_array((p), n_x, n_y, _exp)) { \
        fprintf(stderr, "FAIL: bpoly mismatch\n  got: "); \
        bpoly_print((p), "x", "y"); \
        fprintf(stderr, "\n  expected layout: %dx%d coefs\n", (n_x), (n_y)); \
        ASSERT(0); \
    } \
} while (0)

/* ====================================================================== */
/*  Construction                                                          */
/* ====================================================================== */

static void test_zero(void) {
    BPoly* p = bpoly_zero();
    ASSERT(bpoly_is_zero(p));
    ASSERT(bpoly_deg_x(p) == -1);
    ASSERT(bpoly_deg_y(p) == -1);
    ASSERT(bpoly_lc_x(p) == NULL);
    bpoly_free(p);
}

static void test_set_xcoef_simple(void) {
    /* Build P(x, y) = (1 + 2y) + (3 + 4y)*x. */
    BPoly* p = bpoly_new(2);
    bpoly_set_xcoef(p, 0, zu((int64_t[]){1, 2}, 2));
    bpoly_set_xcoef(p, 1, zu((int64_t[]){3, 4}, 2));
    ASSERT(bpoly_deg_x(p) == 1);
    ASSERT(bpoly_deg_y(p) == 1);
    int64_t expected[] = {1, 2, 3, 4};
    ASSERT(bpoly_eq_array(p, 2, 2, expected));
    bpoly_free(p);
}

static void test_set_xcoef_zero_clears_deg(void) {
    BPoly* p = bpoly_new(3);
    bpoly_set_xcoef(p, 0, zu((int64_t[]){1}, 1));
    bpoly_set_xcoef(p, 2, zu((int64_t[]){1}, 1));
    ASSERT(bpoly_deg_x(p) == 2);
    /* Clear the leading slot. */
    bpoly_set_xcoef(p, 2, zupoly_zero());
    ASSERT(bpoly_deg_x(p) == 0);
    bpoly_free(p);
}

static void test_copy_independence(void) {
    BPoly* p = bp(2, 2, (int64_t[]){1, 2, 3, 4});
    BPoly* q = bpoly_copy(p);
    /* Mutate q's leading slot. */
    bpoly_set_xcoef(q, 0, zu((int64_t[]){99}, 1));
    int64_t expected_p[] = {1, 2, 3, 4};
    ASSERT(bpoly_eq_array(p, 2, 2, expected_p));
    /* p[0][0] = 1 unchanged */
    bpoly_free(p); bpoly_free(q);
}

/* ====================================================================== */
/*  Arithmetic                                                            */
/* ====================================================================== */

static void test_add(void) {
    /* P = (1 + 2y) + 3*x.   Q = 4 + 5*x.
     * P + Q = (5 + 2y) + 8x. */
    BPoly* p = bp(2, 2, (int64_t[]){1, 2, 3, 0});
    BPoly* q = bp(2, 2, (int64_t[]){4, 0, 5, 0});
    BPoly* r = bpoly_add(p, q);
    EXPECT_BPOLY(r, 2, 2, 5, 2, 8, 0);
    bpoly_free(p); bpoly_free(q); bpoly_free(r);
}

static void test_sub(void) {
    BPoly* p = bp(2, 2, (int64_t[]){5, 2, 8, 0});
    BPoly* q = bp(2, 2, (int64_t[]){4, 0, 5, 0});
    BPoly* r = bpoly_sub(p, q);
    EXPECT_BPOLY(r, 2, 2, 1, 2, 3, 0);
    bpoly_free(p); bpoly_free(q); bpoly_free(r);
}

static void test_add_with_cancellation(void) {
    /* P = x.  Q = -x.  P + Q = 0. */
    BPoly* p = bpoly_new(2);
    bpoly_set_xcoef(p, 1, zupoly_from_int(1));
    BPoly* q = bpoly_new(2);
    bpoly_set_xcoef(q, 1, zupoly_from_int(-1));
    BPoly* r = bpoly_add(p, q);
    ASSERT(bpoly_is_zero(r));
    bpoly_free(p); bpoly_free(q); bpoly_free(r);
}

static void test_mul(void) {
    /* (x + y)(x - y) = x^2 - y^2.
     * Storage:
     *   factor 1: cx[0] = y,    cx[1] = 1
     *   factor 2: cx[0] = -y,   cx[1] = 1
     *   product:  cx[0] = -y^2, cx[1] = 0, cx[2] = 1
     * coefs row-major (n_x=3, n_y=3): {0,0,-1, 0,0,0, 1,0,0} */
    BPoly* a = bpoly_new(2);
    bpoly_set_xcoef(a, 0, zu((int64_t[]){0, 1}, 2));   /* y */
    bpoly_set_xcoef(a, 1, zu((int64_t[]){1}, 1));      /* 1 */

    BPoly* b = bpoly_new(2);
    bpoly_set_xcoef(b, 0, zu((int64_t[]){0, -1}, 2));  /* -y */
    bpoly_set_xcoef(b, 1, zu((int64_t[]){1}, 1));      /* 1 */

    BPoly* p = bpoly_mul(a, b);
    EXPECT_BPOLY(p, 3, 3, 0, 0, -1,  0, 0, 0,  1, 0, 0);
    bpoly_free(a); bpoly_free(b); bpoly_free(p);
}

static void test_mul_zero(void) {
    BPoly* a = bp(2, 2, (int64_t[]){1, 2, 3, 4});
    BPoly* z = bpoly_zero();
    BPoly* r = bpoly_mul(a, z);
    ASSERT(bpoly_is_zero(r));
    bpoly_free(a); bpoly_free(z); bpoly_free(r);
}

static void test_neg(void) {
    BPoly* a = bp(2, 2, (int64_t[]){1, -2, 3, -4});
    BPoly* r = bpoly_neg(a);
    EXPECT_BPOLY(r, 2, 2, -1, 2, -3, 4);
    bpoly_free(a); bpoly_free(r);
}

static void test_mul_zupoly(void) {
    /* P = x.  s = (1 + y).  P * s = (1 + y)*x. */
    BPoly* p = bpoly_new(2);
    bpoly_set_xcoef(p, 1, zupoly_from_int(1));
    ZUPoly* s = zu((int64_t[]){1, 1}, 2);
    BPoly* r = bpoly_mul_zupoly(p, s);
    EXPECT_BPOLY(r, 2, 2, 0, 0, 1, 1);
    bpoly_free(p); zupoly_free(s); bpoly_free(r);
}

/* ====================================================================== */
/*  Truncation, evaluation, shift                                         */
/* ====================================================================== */

static void test_truncate_y(void) {
    /* P = (1 + 2y + 3y^2) + (4 + 5y + 6y^2)*x.  Truncate mod y^2. */
    BPoly* p = bpoly_new(2);
    bpoly_set_xcoef(p, 0, zu((int64_t[]){1, 2, 3}, 3));
    bpoly_set_xcoef(p, 1, zu((int64_t[]){4, 5, 6}, 3));

    BPoly* tr = bpoly_truncate_y(p, 2);
    /* Each cx[i] now has only y^0 and y^1 terms. */
    EXPECT_BPOLY(tr, 2, 2, 1, 2,  4, 5);
    bpoly_free(p); bpoly_free(tr);
}

static void test_truncate_y_drops_zero_cx(void) {
    /* P = y * x.  Truncate mod y -> the cx[1] becomes zero, so deg_x
     * collapses to -1 (zero polynomial). */
    BPoly* p = bpoly_new(2);
    bpoly_set_xcoef(p, 1, zu((int64_t[]){0, 1}, 2));  /* y */
    BPoly* tr = bpoly_truncate_y(p, 1);
    ASSERT(bpoly_is_zero(tr));
    bpoly_free(p); bpoly_free(tr);
}

static void test_eval_y(void) {
    /* P = (1 + 2y) + (3 + 4y)*x.  P(x, y=5) = 11 + 23x. */
    BPoly* p = bp(2, 2, (int64_t[]){1, 2, 3, 4});
    ZUPoly* u = bpoly_eval_y_si(p, 5);
    /* Expected: c[0] = 1 + 2*5 = 11, c[1] = 3 + 4*5 = 23. */
    const mpz_t* c0 = zupoly_getcoef(u, 0);
    const mpz_t* c1 = zupoly_getcoef(u, 1);
    ASSERT(c0 && mpz_cmp_ui(*c0, 11) == 0);
    ASSERT(c1 && mpz_cmp_ui(*c1, 23) == 0);
    bpoly_free(p); zupoly_free(u);
}

static void test_shift_y(void) {
    /* P = y * x.  Shift y -> y + 1.  Result: (y + 1) * x. */
    BPoly* p = bpoly_new(2);
    bpoly_set_xcoef(p, 1, zu((int64_t[]){0, 1}, 2));
    BPoly* sh = bpoly_shift_y_si(p, 1);
    /* Expected: cx[0] = 0, cx[1] = (1 + y).  As 2x2 row-major:
     *   {0, 0,  1, 1}. */
    EXPECT_BPOLY(sh, 2, 2, 0, 0,  1, 1);
    bpoly_free(p); bpoly_free(sh);
}

/* ====================================================================== */
/*  Exact division                                                        */
/* ====================================================================== */

static void test_divexact_simple(void) {
    /* (x + y)(x - y) = x^2 - y^2.  Divide by (x + y) -> (x - y). */
    BPoly* a = bpoly_new(2);
    bpoly_set_xcoef(a, 0, zu((int64_t[]){0, 1}, 2));
    bpoly_set_xcoef(a, 1, zu((int64_t[]){1}, 1));
    BPoly* b = bpoly_new(2);
    bpoly_set_xcoef(b, 0, zu((int64_t[]){0, -1}, 2));
    bpoly_set_xcoef(b, 1, zu((int64_t[]){1}, 1));
    BPoly* prod = bpoly_mul(a, b);

    BPoly* q = bpoly_divexact(prod, a);
    ASSERT(q != NULL);
    ASSERT(bpoly_eq(q, b));

    bpoly_free(a); bpoly_free(b); bpoly_free(prod); bpoly_free(q);
}

static void test_divexact_failure(void) {
    /* x^2 + 1 over Z[x, y] is not divisible by x + y (the quotient
     * would have y in it -- fine -- but the remainder is y^2 + 1, not
     * zero). */
    BPoly* a = bpoly_new(3);
    bpoly_set_xcoef(a, 0, zu((int64_t[]){1}, 1));
    bpoly_set_xcoef(a, 2, zu((int64_t[]){1}, 1));
    BPoly* b = bpoly_new(2);
    bpoly_set_xcoef(b, 0, zu((int64_t[]){0, 1}, 2));
    bpoly_set_xcoef(b, 1, zu((int64_t[]){1}, 1));

    BPoly* q = bpoly_divexact(a, b);
    ASSERT(q == NULL);
    bpoly_free(a); bpoly_free(b);
}

/* ====================================================================== */
/*  Compositional sanity                                                  */
/* ====================================================================== */

static void test_factor_then_multiply(void) {
    /* a = x + y, b = x - y, c = x^2 + 1.
     * Verify (a*b)*c == ((a*b))*c == a*(b*c). */
    BPoly* a = bpoly_new(2);
    bpoly_set_xcoef(a, 0, zu((int64_t[]){0, 1}, 2));
    bpoly_set_xcoef(a, 1, zu((int64_t[]){1}, 1));
    BPoly* b = bpoly_new(2);
    bpoly_set_xcoef(b, 0, zu((int64_t[]){0, -1}, 2));
    bpoly_set_xcoef(b, 1, zu((int64_t[]){1}, 1));
    BPoly* c = bpoly_new(3);
    bpoly_set_xcoef(c, 0, zu((int64_t[]){1}, 1));
    bpoly_set_xcoef(c, 2, zu((int64_t[]){1}, 1));

    BPoly* ab = bpoly_mul(a, b);
    BPoly* bc = bpoly_mul(b, c);
    BPoly* abc1 = bpoly_mul(ab, c);
    BPoly* abc2 = bpoly_mul(a, bc);
    ASSERT(bpoly_eq(abc1, abc2));

    bpoly_free(a); bpoly_free(b); bpoly_free(c);
    bpoly_free(ab); bpoly_free(bc);
    bpoly_free(abc1); bpoly_free(abc2);
}

static void test_eval_after_shift_invariant(void) {
    /* Shift y -> y+1 then evaluate y=0 should equal evaluating y=1
     * on the original.  Tests shift composition correctness. */
    BPoly* p = bp(2, 3, (int64_t[]){1, 2, 3,   4, 5, 6});
    BPoly* sh = bpoly_shift_y_si(p, 1);
    ZUPoly* via_shift = bpoly_eval_y_si(sh, 0);
    ZUPoly* direct = bpoly_eval_y_si(p, 1);
    ASSERT(zupoly_eq(via_shift, direct));
    bpoly_free(p); bpoly_free(sh);
    zupoly_free(via_shift); zupoly_free(direct);
}

static void test_truncate_after_mul(void) {
    /* (1 + y + y^2) * (1 + y) = 1 + 2y + 2y^2 + y^3.
     * Storage as constant in x: cx[0] is the y-poly.  Truncate mod y^2:
     * 1 + 2y. */
    BPoly* a = bpoly_new(1);
    bpoly_set_xcoef(a, 0, zu((int64_t[]){1, 1, 1}, 3));
    BPoly* b = bpoly_new(1);
    bpoly_set_xcoef(b, 0, zu((int64_t[]){1, 1}, 2));
    BPoly* p = bpoly_mul_truncate_y(a, b, 2);
    /* Expected: 1x1 layout, c[0][0]=1, c[0][1]=2. */
    EXPECT_BPOLY(p, 1, 2, 1, 2);
    bpoly_free(a); bpoly_free(b); bpoly_free(p);
}

/* ====================================================================== */
/*  Conversion to / from Expr                                             */
/* ====================================================================== */

static void test_expr_roundtrip(void) {
    /* x^2 + x*y + y^2 + 1 */
    Expr* xv = parse_expression("x");
    Expr* yv = parse_expression("y");
    Expr* e = parse_expression("x^2 + x*y + y^2 + 1");
    Expr* ev = evaluate(e); expr_free(e);

    BPoly* p = expr_to_bpoly(ev, xv, yv);
    ASSERT(p != NULL);

    /* deg_x = 2, deg_y = 2.  Coefficients:
     *   cx[0] = 1 + y^2
     *   cx[1] = y
     *   cx[2] = 1
     * As 3x3 row-major:
     *   {1, 0, 1,  0, 1, 0,  1, 0, 0}. */
    EXPECT_BPOLY(p, 3, 3,  1, 0, 1,  0, 1, 0,  1, 0, 0);

    Expr* back = bpoly_to_expr(p, xv, yv);
    ASSERT(expr_eq(back, ev));

    expr_free(back); expr_free(ev); expr_free(xv); expr_free(yv);
    bpoly_free(p);
}

static void test_expr_zero(void) {
    Expr* xv = parse_expression("x");
    Expr* yv = parse_expression("y");
    Expr* e = parse_expression("0");
    Expr* ev = evaluate(e); expr_free(e);
    BPoly* p = expr_to_bpoly(ev, xv, yv);
    ASSERT(p != NULL);
    ASSERT(bpoly_is_zero(p));
    bpoly_free(p); expr_free(ev); expr_free(xv); expr_free(yv);
}

static void test_expr_unrelated_var_returns_null(void) {
    /* z is neither x nor y; conversion must fail. */
    Expr* xv = parse_expression("x");
    Expr* yv = parse_expression("y");
    Expr* e = parse_expression("x*y + z");
    Expr* ev = evaluate(e); expr_free(e);
    BPoly* p = expr_to_bpoly(ev, xv, yv);
    ASSERT(p == NULL);
    expr_free(ev); expr_free(xv); expr_free(yv);
}

/* ====================================================================== */
/*  Main                                                                  */
/* ====================================================================== */

int main(void) {
    symtab_init();
    core_init();

    printf("Running BPoly tests...\n");

    /* Construction */
    TEST(test_zero);
    TEST(test_set_xcoef_simple);
    TEST(test_set_xcoef_zero_clears_deg);
    TEST(test_copy_independence);

    /* Arithmetic */
    TEST(test_add);
    TEST(test_sub);
    TEST(test_add_with_cancellation);
    TEST(test_mul);
    TEST(test_mul_zero);
    TEST(test_neg);
    TEST(test_mul_zupoly);

    /* Truncation, eval, shift */
    TEST(test_truncate_y);
    TEST(test_truncate_y_drops_zero_cx);
    TEST(test_eval_y);
    TEST(test_shift_y);

    /* Division */
    TEST(test_divexact_simple);
    TEST(test_divexact_failure);

    /* Compositional */
    TEST(test_factor_then_multiply);
    TEST(test_eval_after_shift_invariant);
    TEST(test_truncate_after_mul);

    /* Conversion */
    TEST(test_expr_roundtrip);
    TEST(test_expr_zero);
    TEST(test_expr_unrelated_var_returns_null);

    printf("All BPoly tests passed!\n");
    return 0;
}
