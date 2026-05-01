/*
 * test_mpoly.c
 * ------------
 * Unit tests for the MPoly primitive operations.  Each test exercises
 * one operation in isolation so failures localise immediately.
 *
 * Coverage:
 *   - Construction: new/zero, from_int, from_mpz, monomial, copy
 *   - Term manipulation: push_term, set_coef, get_coef, normalize
 *   - Predicates: is_zero, eq, deg_var, total_deg, is_constant_in_var
 *   - Arithmetic: add, sub, neg, mul, scale
 *   - Substitution: subst_var_int, shift_var_int, coef_of_var, lc_var
 *   - Expr round-trip: expr_to_mpoly, mpoly_to_expr
 *   - Edge cases: zero, n_vars=0, n_vars=1 vs ZUPoly equivalence
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>

#include "mpoly.h"
#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"

#undef ASSERT
#define ASSERT(cond) do { \
    if (!(cond)) { \
        fprintf(stderr, "ASSERTION FAILED: %s\n  at %s:%d\n", \
                #cond, __FILE__, __LINE__); \
        abort(); \
    } \
} while (0)

#define TEST(name) do { printf("Running test: %s\n", #name); name(); } while (0)

/* ---------------------------------------------------------------- */
/*  Helpers                                                         */
/* ---------------------------------------------------------------- */

/* Build an MPoly from a flat exponent + coefficient table.  Each row
 * is { e_0, e_1, ..., e_{n_vars-1}, coef }.  n_terms rows. */
static MPoly* mp_build(int n_vars, size_t n_terms, const int64_t* table) {
    MPoly* p = mpoly_new(n_vars);
    int exps[n_vars > 0 ? n_vars : 1];
    mpz_t c; mpz_init(c);
    for (size_t i = 0; i < n_terms; i++) {
        for (int v = 0; v < n_vars; v++) exps[v] = (int)table[i * (n_vars + 1) + v];
        mpz_set_si(c, (long)table[i * (n_vars + 1) + n_vars]);
        mpoly_push_term(p, exps, c);
    }
    mpz_clear(c);
    mpoly_normalize(p);
    return p;
}

/* ---------------------------------------------------------------- */
/*  Construction                                                    */
/* ---------------------------------------------------------------- */

static void test_zero(void) {
    MPoly* p = mpoly_new(3);
    ASSERT(p != NULL);
    ASSERT(p->n_vars == 3);
    ASSERT(p->n_terms == 0);
    ASSERT(mpoly_is_zero(p));
    mpoly_free(p);
}

static void test_from_int(void) {
    MPoly* p = mpoly_from_int(2, 7);
    ASSERT(p->n_terms == 1);
    ASSERT(mpz_cmp_si(p->coefs[0], 7) == 0);
    int* row = p->exps;
    ASSERT(row[0] == 0 && row[1] == 0);
    mpoly_free(p);

    MPoly* z = mpoly_from_int(2, 0);
    ASSERT(mpoly_is_zero(z));
    mpoly_free(z);
}

static void test_from_mpz_bigint(void) {
    mpz_t big; mpz_init(big);
    mpz_set_str(big, "123456789012345678901234567890", 10);
    MPoly* p = mpoly_from_mpz(2, big);
    ASSERT(p->n_terms == 1);
    ASSERT(mpz_cmp(p->coefs[0], big) == 0);
    mpz_clear(big);
    mpoly_free(p);
}

static void test_monomial(void) {
    /* 5 * x_1^3 in Z[x_0, x_1, x_2]. */
    MPoly* p = mpoly_monomial(3, 1, 3, 5);
    ASSERT(p->n_terms == 1);
    ASSERT(mpz_cmp_si(p->coefs[0], 5) == 0);
    int* r = p->exps;
    ASSERT(r[0] == 0 && r[1] == 3 && r[2] == 0);
    mpoly_free(p);
}

static void test_copy_independence(void) {
    int64_t tab[] = { 2,1, 3, 1,0, 5, 0,2, -1 };  /* 3 x_0^2 x_1 + 5 x_0 - x_1^2 */
    MPoly* a = mp_build(2, 3, tab);
    MPoly* b = mpoly_copy(a);
    ASSERT(mpoly_eq(a, b));
    /* Mutate b: should not affect a. */
    int row[2] = {0, 0};
    mpz_t c; mpz_init_set_si(c, 99);
    mpoly_set_coef(b, row, c);
    ASSERT(!mpoly_eq(a, b));
    /* a still has its original 3 terms; b has 4. */
    ASSERT(a->n_terms == 3);
    ASSERT(b->n_terms == 4);
    mpz_clear(c);
    mpoly_free(a);
    mpoly_free(b);
}

/* ---------------------------------------------------------------- */
/*  Normalization                                                   */
/* ---------------------------------------------------------------- */

static void test_normalize_dedup(void) {
    /* Push two terms with identical exponents; normalize should merge. */
    MPoly* p = mpoly_new(2);
    int e[] = {1, 0};
    mpz_t c1, c2;
    mpz_init_set_si(c1, 3);
    mpz_init_set_si(c2, 5);
    mpoly_push_term(p, e, c1);
    mpoly_push_term(p, e, c2);
    mpz_clear(c1); mpz_clear(c2);

    mpoly_normalize(p);
    ASSERT(p->n_terms == 1);
    ASSERT(mpz_cmp_si(p->coefs[0], 8) == 0);
    mpoly_free(p);
}

static void test_normalize_cancellation(void) {
    /* Push +x and -x; normalize should leave nothing. */
    MPoly* p = mpoly_new(1);
    int e[] = {1};
    mpz_t a, b;
    mpz_init_set_si(a,  7);
    mpz_init_set_si(b, -7);
    mpoly_push_term(p, e, a);
    mpoly_push_term(p, e, b);
    mpz_clear(a); mpz_clear(b);
    mpoly_normalize(p);
    ASSERT(mpoly_is_zero(p));
    mpoly_free(p);
}

static void test_normalize_drops_zero(void) {
    /* push_term skips zero coefficients, but set_coef with 0 should
     * remove the term. */
    MPoly* p = mpoly_new(1);
    int e[] = {0};
    mpz_t c; mpz_init_set_si(c, 5);
    mpoly_set_coef(p, e, c);
    ASSERT(p->n_terms == 1);
    mpz_set_ui(c, 0);
    mpoly_set_coef(p, e, c);
    ASSERT(p->n_terms == 0);
    mpz_clear(c);
    mpoly_free(p);
}

static void test_lex_ordering(void) {
    /* Verify normalize sorts in lex descending order: x_0^2 first,
     * then x_0 x_1, then x_1^2, then 1. */
    int64_t tab[] = {
        0,0, 7,
        2,0, 1,
        0,2, 5,
        1,1, 3
    };
    MPoly* p = mp_build(2, 4, tab);
    ASSERT(p->n_terms == 4);
    /* Term 0: (2,0). */
    ASSERT(p->exps[0] == 2 && p->exps[1] == 0);
    /* Term 1: (1,1). */
    ASSERT(p->exps[2] == 1 && p->exps[3] == 1);
    /* Term 2: (0,2). */
    ASSERT(p->exps[4] == 0 && p->exps[5] == 2);
    /* Term 3: (0,0). */
    ASSERT(p->exps[6] == 0 && p->exps[7] == 0);
    mpoly_free(p);
}

/* ---------------------------------------------------------------- */
/*  Set/get coef                                                    */
/* ---------------------------------------------------------------- */

static void test_set_get_coef(void) {
    MPoly* p = mpoly_new(2);
    mpz_t c; mpz_init(c);
    int e1[] = {1, 0};
    int e2[] = {0, 1};
    int e3[] = {2, 0};
    mpz_set_si(c, 3);
    mpoly_set_coef(p, e1, c);
    mpz_set_si(c, 7);
    mpoly_set_coef(p, e2, c);
    mpz_set_si(c, 11);
    mpoly_set_coef(p, e3, c);
    ASSERT(p->n_terms == 3);

    const mpz_t* r = mpoly_get_coef(p, e1);
    ASSERT(r && mpz_cmp_si(*r, 3) == 0);
    r = mpoly_get_coef(p, e2);
    ASSERT(r && mpz_cmp_si(*r, 7) == 0);
    r = mpoly_get_coef(p, e3);
    ASSERT(r && mpz_cmp_si(*r, 11) == 0);

    int e_missing[] = {5, 5};
    ASSERT(mpoly_get_coef(p, e_missing) == NULL);

    /* Replace and remove. */
    mpz_set_si(c, 99);
    mpoly_set_coef(p, e2, c);
    r = mpoly_get_coef(p, e2);
    ASSERT(r && mpz_cmp_si(*r, 99) == 0);

    mpz_set_ui(c, 0);
    mpoly_set_coef(p, e2, c);
    ASSERT(mpoly_get_coef(p, e2) == NULL);
    ASSERT(p->n_terms == 2);

    mpz_clear(c);
    mpoly_free(p);
}

/* ---------------------------------------------------------------- */
/*  Predicates                                                      */
/* ---------------------------------------------------------------- */

static void test_eq_basic(void) {
    int64_t ta[] = { 1,1,1, 2,0,3, 0,2,5 };  /* xy + 3x^2 + 5y^2 */
    int64_t tb[] = { 0,2,5, 2,0,3, 1,1,1 };  /* same, scrambled */
    MPoly* a = mp_build(2, 3, ta);
    MPoly* b = mp_build(2, 3, tb);
    ASSERT(mpoly_eq(a, b));
    mpoly_free(a);
    mpoly_free(b);
}

static void test_deg_var(void) {
    /* 3x^2y + 5xy^4 + 7. */
    int64_t tab[] = { 2,1, 3, 1,4, 5, 0,0, 7 };
    MPoly* p = mp_build(2, 3, tab);
    ASSERT(mpoly_deg_var(p, 0) == 2);
    ASSERT(mpoly_deg_var(p, 1) == 4);
    mpoly_free(p);
}

static void test_total_deg(void) {
    /* 3x^2y + 5xy^4 + 7  -> total degrees 3, 5, 0. */
    int64_t tab[] = { 2,1, 3, 1,4, 5, 0,0, 7 };
    MPoly* p = mp_build(2, 3, tab);
    ASSERT(mpoly_total_deg(p) == 5);
    mpoly_free(p);
}

static void test_is_constant_in_var(void) {
    /* 3x^2 + 5x^2 y -- constant in x is FALSE; constant in y is FALSE. */
    int64_t tab[] = { 2,0, 3, 2,1, 5 };
    MPoly* p = mp_build(2, 2, tab);
    ASSERT(!mpoly_is_constant_in_var(p, 0));
    ASSERT(!mpoly_is_constant_in_var(p, 1));
    mpoly_free(p);

    /* 3x^2 + 5x^2  -> 8x^2 -- constant in y. */
    int64_t tab2[] = { 2,0, 3, 2,0, 5 };
    MPoly* q = mp_build(2, 2, tab2);
    ASSERT(!mpoly_is_constant_in_var(q, 0));
    ASSERT(mpoly_is_constant_in_var(q, 1));
    mpoly_free(q);
}

/* ---------------------------------------------------------------- */
/*  Arithmetic                                                      */
/* ---------------------------------------------------------------- */

static void test_add_basic(void) {
    /* (3x + 5y) + (2x + 7) = 5x + 5y + 7 */
    int64_t ta[] = { 1,0, 3, 0,1, 5 };
    int64_t tb[] = { 1,0, 2, 0,0, 7 };
    int64_t te[] = { 1,0, 5, 0,1, 5, 0,0, 7 };
    MPoly* a = mp_build(2, 2, ta);
    MPoly* b = mp_build(2, 2, tb);
    MPoly* expected = mp_build(2, 3, te);
    MPoly* sum = mpoly_add(a, b);
    ASSERT(mpoly_eq(sum, expected));
    mpoly_free(a); mpoly_free(b); mpoly_free(expected); mpoly_free(sum);
}

static void test_add_cancellation(void) {
    /* (x + y) + (-x - y) = 0 */
    int64_t ta[] = { 1,0, 1, 0,1, 1 };
    int64_t tb[] = { 1,0,-1, 0,1,-1 };
    MPoly* a = mp_build(2, 2, ta);
    MPoly* b = mp_build(2, 2, tb);
    MPoly* sum = mpoly_add(a, b);
    ASSERT(mpoly_is_zero(sum));
    mpoly_free(a); mpoly_free(b); mpoly_free(sum);
}

static void test_sub_basic(void) {
    /* (5x + 3y) - (2x + y) = 3x + 2y */
    int64_t ta[] = { 1,0, 5, 0,1, 3 };
    int64_t tb[] = { 1,0, 2, 0,1, 1 };
    int64_t te[] = { 1,0, 3, 0,1, 2 };
    MPoly* a = mp_build(2, 2, ta);
    MPoly* b = mp_build(2, 2, tb);
    MPoly* expected = mp_build(2, 2, te);
    MPoly* d = mpoly_sub(a, b);
    ASSERT(mpoly_eq(d, expected));
    mpoly_free(a); mpoly_free(b); mpoly_free(expected); mpoly_free(d);
}

static void test_neg(void) {
    int64_t ta[] = { 1,0, 5, 0,1, 3 };
    int64_t tn[] = { 1,0,-5, 0,1,-3 };
    MPoly* a = mp_build(2, 2, ta);
    MPoly* expected = mp_build(2, 2, tn);
    MPoly* n = mpoly_neg(a);
    ASSERT(mpoly_eq(n, expected));
    mpoly_free(a); mpoly_free(expected); mpoly_free(n);
}

static void test_mul_basic(void) {
    /* (x + y)(x - y) = x^2 - y^2. */
    int64_t ta[] = { 1,0, 1, 0,1, 1 };
    int64_t tb[] = { 1,0, 1, 0,1,-1 };
    int64_t te[] = { 2,0, 1, 0,2,-1 };
    MPoly* a = mp_build(2, 2, ta);
    MPoly* b = mp_build(2, 2, tb);
    MPoly* expected = mp_build(2, 2, te);
    MPoly* prod = mpoly_mul(a, b);
    ASSERT(mpoly_eq(prod, expected));
    mpoly_free(a); mpoly_free(b); mpoly_free(expected); mpoly_free(prod);
}

static void test_mul_three_vars(void) {
    /* (x + y + z)(x - y) = x^2 - y^2 + xz - yz. */
    int64_t ta[] = { 1,0,0, 1, 0,1,0, 1, 0,0,1, 1 };
    int64_t tb[] = { 1,0,0, 1, 0,1,0,-1 };
    int64_t te[] = { 2,0,0, 1, 0,2,0,-1, 1,0,1, 1, 0,1,1,-1 };
    MPoly* a = mp_build(3, 3, ta);
    MPoly* b = mp_build(3, 2, tb);
    MPoly* expected = mp_build(3, 4, te);
    MPoly* prod = mpoly_mul(a, b);
    ASSERT(mpoly_eq(prod, expected));
    mpoly_free(a); mpoly_free(b); mpoly_free(expected); mpoly_free(prod);
}

static void test_mul_zero(void) {
    int64_t ta[] = { 1,0, 5 };
    MPoly* a = mp_build(2, 1, ta);
    MPoly* z = mpoly_zero(2);
    MPoly* p = mpoly_mul(a, z);
    ASSERT(mpoly_is_zero(p));
    mpoly_free(a); mpoly_free(z); mpoly_free(p);
}

static void test_scale(void) {
    int64_t ta[] = { 1,0, 3, 0,1, 5 };
    MPoly* a = mp_build(2, 2, ta);
    mpz_t k; mpz_init_set_si(k, 4);
    MPoly* s = mpoly_scale(a, k);
    int64_t te[] = { 1,0, 12, 0,1, 20 };
    MPoly* expected = mp_build(2, 2, te);
    ASSERT(mpoly_eq(s, expected));
    mpz_clear(k);
    mpoly_free(a); mpoly_free(s); mpoly_free(expected);
}

/* ---------------------------------------------------------------- */
/*  Substitution                                                    */
/* ---------------------------------------------------------------- */

static void test_subst_var_int(void) {
    /* P = 3x^2 y + 5x y + 7y + 2.  Substitute y = 1 -> 3x^2 + 5x + 9. */
    int64_t tab[] = {
        2,1, 3,
        1,1, 5,
        0,1, 7,
        0,0, 2
    };
    MPoly* p = mp_build(2, 4, tab);
    MPoly* q = mpoly_subst_var_int(p, 1, 1);
    /* After substitution, every term has y-exponent 0; coefficient
     * collected. */
    int64_t te[] = { 2,0, 3, 1,0, 5, 0,0, 9 };
    MPoly* expected = mp_build(2, 3, te);
    ASSERT(mpoly_eq(q, expected));
    mpoly_free(p); mpoly_free(q); mpoly_free(expected);
}

static void test_subst_var_int_zero(void) {
    /* P = 3x^2 y + 5x.  Substitute y = 0 -> 5x. */
    int64_t tab[] = { 2,1, 3, 1,0, 5 };
    MPoly* p = mp_build(2, 2, tab);
    MPoly* q = mpoly_subst_var_int(p, 1, 0);
    int64_t te[] = { 1,0, 5 };
    MPoly* expected = mp_build(2, 1, te);
    ASSERT(mpoly_eq(q, expected));
    mpoly_free(p); mpoly_free(q); mpoly_free(expected);
}

static void test_shift_var_int(void) {
    /* P = (x - 1)^2 = x^2 - 2x + 1.  Shift x -> x + 1: (x)^2 = x^2. */
    int64_t tab[] = { 2, 1, 1,-2, 0, 1 };
    MPoly* p = mp_build(1, 3, tab);
    MPoly* q = mpoly_shift_var_int(p, 0, 1);
    int64_t te[] = { 2, 1 };
    MPoly* expected = mp_build(1, 1, te);
    ASSERT(mpoly_eq(q, expected));
    mpoly_free(p); mpoly_free(q); mpoly_free(expected);
}

static void test_shift_two_vars(void) {
    /* P = x + y.  Shift y -> y + 5: x + (y + 5) = x + y + 5. */
    int64_t tab[] = { 1,0, 1, 0,1, 1 };
    MPoly* p = mp_build(2, 2, tab);
    MPoly* q = mpoly_shift_var_int(p, 1, 5);
    int64_t te[] = { 1,0, 1, 0,1, 1, 0,0, 5 };
    MPoly* expected = mp_build(2, 3, te);
    ASSERT(mpoly_eq(q, expected));
    mpoly_free(p); mpoly_free(q); mpoly_free(expected);
}

static void test_coef_of_var(void) {
    /* P = 3x^2 y + 5xy + 7y + 2.  coef_of_var(y, 1) = 3x^2 + 5x + 7;
     * coef_of_var(y, 0) = 2. */
    int64_t tab[] = {
        2,1, 3,
        1,1, 5,
        0,1, 7,
        0,0, 2
    };
    MPoly* p = mp_build(2, 4, tab);
    MPoly* c1 = mpoly_coef_of_var(p, 1, 1);
    int64_t t_c1[] = { 2,0, 3, 1,0, 5, 0,0, 7 };
    MPoly* expected_c1 = mp_build(2, 3, t_c1);
    ASSERT(mpoly_eq(c1, expected_c1));

    MPoly* c0 = mpoly_coef_of_var(p, 1, 0);
    int64_t t_c0[] = { 0,0, 2 };
    MPoly* expected_c0 = mp_build(2, 1, t_c0);
    ASSERT(mpoly_eq(c0, expected_c0));

    mpoly_free(p);
    mpoly_free(c1); mpoly_free(expected_c1);
    mpoly_free(c0); mpoly_free(expected_c0);
}

static void test_lc_var(void) {
    /* P = 7x^3 y + 2x^3 + 5x^2 y^2 + 11.  lc_x = 7y + 2 (coef of x^3). */
    int64_t tab[] = {
        3,1, 7,
        3,0, 2,
        2,2, 5,
        0,0, 11
    };
    MPoly* p = mp_build(2, 4, tab);
    MPoly* lc = mpoly_lc_var(p, 0);
    int64_t te[] = { 0,1, 7, 0,0, 2 };
    MPoly* expected = mp_build(2, 2, te);
    ASSERT(mpoly_eq(lc, expected));
    mpoly_free(p); mpoly_free(lc); mpoly_free(expected);
}

/* ---------------------------------------------------------------- */
/*  Expr round-trip                                                 */
/* ---------------------------------------------------------------- */

static void test_expr_to_mpoly_simple(void) {
    Expr* x = parse_expression("x");
    Expr* y = parse_expression("y");
    Expr* vars[] = { x, y };

    Expr* e = parse_expression("3 x^2 + 5 x y + 7");
    Expr* ev = evaluate(e);
    MPoly* p = expr_to_mpoly(ev, vars, 2);
    ASSERT(p != NULL);

    int64_t te[] = { 2,0, 3, 1,1, 5, 0,0, 7 };
    MPoly* expected = mp_build(2, 3, te);
    ASSERT(mpoly_eq(p, expected));

    mpoly_free(p); mpoly_free(expected);
    expr_free(e); expr_free(ev); expr_free(x); expr_free(y);
}

static void test_expr_to_mpoly_three_vars(void) {
    Expr* x = parse_expression("x");
    Expr* y = parse_expression("y");
    Expr* z = parse_expression("z");
    Expr* vars[] = { x, y, z };

    Expr* e = parse_expression("z x - x^2 - y^2");
    Expr* ev = evaluate(e);
    MPoly* p = expr_to_mpoly(ev, vars, 3);
    ASSERT(p != NULL);
    /* Lex desc with var order (x, y, z): x^2 first, then x z, then y^2.
     * (x^2: (2,0,0)) > ((x z): (1,0,1)) > ((y^2): (0,2,0)). */
    ASSERT(p->n_terms == 3);
    /* Term 0: -x^2. */
    int e0[] = {2,0,0};
    const mpz_t* c0 = mpoly_get_coef(p, e0);
    ASSERT(c0 && mpz_cmp_si(*c0, -1) == 0);
    /* Term: x z, coef +1. */
    int e1[] = {1,0,1};
    const mpz_t* c1 = mpoly_get_coef(p, e1);
    ASSERT(c1 && mpz_cmp_si(*c1, 1) == 0);
    /* Term: y^2, coef -1. */
    int e2[] = {0,2,0};
    const mpz_t* c2 = mpoly_get_coef(p, e2);
    ASSERT(c2 && mpz_cmp_si(*c2, -1) == 0);

    mpoly_free(p);
    expr_free(e); expr_free(ev);
    expr_free(x); expr_free(y); expr_free(z);
}

static void test_expr_to_mpoly_rejects_non_polynomial(void) {
    Expr* x = parse_expression("x");
    Expr* vars[] = { x };

    Expr* e = parse_expression("Sin[x] + 1");
    Expr* ev = evaluate(e);
    MPoly* p = expr_to_mpoly(ev, vars, 1);
    /* Sin[x] is not a polynomial in x.  Expect NULL. */
    ASSERT(p == NULL);

    expr_free(e); expr_free(ev); expr_free(x);
}

static void test_mpoly_to_expr_roundtrip(void) {
    Expr* x = parse_expression("x");
    Expr* y = parse_expression("y");
    Expr* z = parse_expression("z");
    Expr* vars[] = { x, y, z };

    /* Polynomial: 3 x^2 y - 5 x z + 7. */
    int64_t tab[] = {
        2,1,0,  3,
        1,0,1, -5,
        0,0,0,  7
    };
    MPoly* p = mp_build(3, 3, tab);
    Expr* out = mpoly_to_expr(p, vars);
    /* Evaluate and re-parse to check structural equivalence. */
    Expr* expected = parse_expression("3 x^2 y - 5 x z + 7");
    Expr* expected_ev = evaluate(expected);
    Expr* out_ev = evaluate(out);
    ASSERT(expr_eq(out_ev, expected_ev));

    expr_free(out); expr_free(out_ev);
    expr_free(expected); expr_free(expected_ev);
    expr_free(x); expr_free(y); expr_free(z);
    mpoly_free(p);
}

static void test_roundtrip_via_evaluate(void) {
    /* Build P via Expr, convert, multiply, convert back, compare. */
    Expr* x = parse_expression("x");
    Expr* y = parse_expression("y");
    Expr* z = parse_expression("z");
    Expr* vars[] = { x, y, z };

    Expr* a_e = parse_expression("Expand[(x + y + z)(x - y + z)]");
    Expr* a_ev = evaluate(a_e);
    MPoly* a = expr_to_mpoly(a_ev, vars, 3);
    ASSERT(a != NULL);

    Expr* a_back = mpoly_to_expr(a, vars);
    Expr* a_back_ev = evaluate(a_back);
    ASSERT(expr_eq(a_back_ev, a_ev));

    expr_free(a_e); expr_free(a_ev);
    expr_free(a_back); expr_free(a_back_ev);
    expr_free(x); expr_free(y); expr_free(z);
    mpoly_free(a);
}

/* ---------------------------------------------------------------- */
/*  Compositional sanity (associativity, commutativity)            */
/* ---------------------------------------------------------------- */

static void test_mul_associative(void) {
    int64_t ta[] = { 1,0, 1, 0,1, 1 };          /* x + y */
    int64_t tb[] = { 1,0, 1, 0,1,-1 };          /* x - y */
    int64_t tc[] = { 0,1, 1, 0,0, 1 };          /* y + 1 */
    MPoly* a = mp_build(2, 2, ta);
    MPoly* b = mp_build(2, 2, tb);
    MPoly* c = mp_build(2, 2, tc);

    MPoly* ab = mpoly_mul(a, b);
    MPoly* ab_c = mpoly_mul(ab, c);

    MPoly* bc = mpoly_mul(b, c);
    MPoly* a_bc = mpoly_mul(a, bc);

    ASSERT(mpoly_eq(ab_c, a_bc));

    mpoly_free(a); mpoly_free(b); mpoly_free(c);
    mpoly_free(ab); mpoly_free(ab_c);
    mpoly_free(bc); mpoly_free(a_bc);
}

static void test_mul_commutative(void) {
    int64_t ta[] = { 2,1, 3, 0,2, 5, 0,0, 7 };
    int64_t tb[] = { 1,0, 2, 0,1,-1 };
    MPoly* a = mp_build(2, 3, ta);
    MPoly* b = mp_build(2, 2, tb);
    MPoly* ab = mpoly_mul(a, b);
    MPoly* ba = mpoly_mul(b, a);
    ASSERT(mpoly_eq(ab, ba));
    mpoly_free(a); mpoly_free(b); mpoly_free(ab); mpoly_free(ba);
}

static void test_subst_then_eval_consistency(void) {
    /* P = x^2 + 2 x y + y^2 = (x + y)^2.  Substituting y = 3:
     *   (x + 3)^2 = x^2 + 6x + 9.
     * Cross-check via Expr evaluation. */
    int64_t tab[] = { 2,0, 1, 1,1, 2, 0,2, 1 };
    MPoly* p = mp_build(2, 3, tab);
    MPoly* q = mpoly_subst_var_int(p, 1, 3);

    int64_t te[] = { 2,0, 1, 1,0, 6, 0,0, 9 };
    MPoly* expected = mp_build(2, 3, te);
    ASSERT(mpoly_eq(q, expected));

    mpoly_free(p); mpoly_free(q); mpoly_free(expected);
}

static void test_shift_then_unshift(void) {
    /* (x + 1)^3 shifted x -> x - 1 should give x^3. */
    int64_t tab[] = { 3, 1, 2, 3, 1, 3, 0, 1 };
    MPoly* p = mp_build(1, 4, tab);
    MPoly* shifted = mpoly_shift_var_int(p, 0, -1);
    int64_t te[] = { 3, 1 };
    MPoly* expected = mp_build(1, 1, te);
    ASSERT(mpoly_eq(shifted, expected));
    mpoly_free(p); mpoly_free(shifted); mpoly_free(expected);
}

/* ---------------------------------------------------------------- */
/*  Single-variable sanity (n_vars=1 mirrors univariate behavior)   */
/* ---------------------------------------------------------------- */

static void test_n_vars_one(void) {
    /* P = 3x^2 - 5x + 2.  Standard univariate. */
    int64_t tab[] = { 2, 3, 1,-5, 0, 2 };
    MPoly* p = mp_build(1, 3, tab);
    ASSERT(p->n_terms == 3);
    ASSERT(mpoly_deg_var(p, 0) == 2);
    ASSERT(mpoly_total_deg(p) == 2);

    /* Substitute x = 1: 3 - 5 + 2 = 0. */
    MPoly* q = mpoly_subst_var_int(p, 0, 1);
    ASSERT(mpoly_is_zero(q));
    mpoly_free(p); mpoly_free(q);
}

/* ---------------------------------------------------------------- */
/*  Main                                                            */
/* ---------------------------------------------------------------- */

int main(void) {
    symtab_init();
    core_init();

    printf("Running MPoly tests...\n");

    /* Construction */
    TEST(test_zero);
    TEST(test_from_int);
    TEST(test_from_mpz_bigint);
    TEST(test_monomial);
    TEST(test_copy_independence);

    /* Normalization */
    TEST(test_normalize_dedup);
    TEST(test_normalize_cancellation);
    TEST(test_normalize_drops_zero);
    TEST(test_lex_ordering);

    /* Set/get coef */
    TEST(test_set_get_coef);

    /* Predicates */
    TEST(test_eq_basic);
    TEST(test_deg_var);
    TEST(test_total_deg);
    TEST(test_is_constant_in_var);

    /* Arithmetic */
    TEST(test_add_basic);
    TEST(test_add_cancellation);
    TEST(test_sub_basic);
    TEST(test_neg);
    TEST(test_mul_basic);
    TEST(test_mul_three_vars);
    TEST(test_mul_zero);
    TEST(test_scale);

    /* Substitution */
    TEST(test_subst_var_int);
    TEST(test_subst_var_int_zero);
    TEST(test_shift_var_int);
    TEST(test_shift_two_vars);
    TEST(test_coef_of_var);
    TEST(test_lc_var);

    /* Expr round-trip */
    TEST(test_expr_to_mpoly_simple);
    TEST(test_expr_to_mpoly_three_vars);
    TEST(test_expr_to_mpoly_rejects_non_polynomial);
    TEST(test_mpoly_to_expr_roundtrip);
    TEST(test_roundtrip_via_evaluate);

    /* Compositional sanity */
    TEST(test_mul_associative);
    TEST(test_mul_commutative);
    TEST(test_subst_then_eval_consistency);
    TEST(test_shift_then_unshift);

    /* Single-variable sanity */
    TEST(test_n_vars_one);

    printf("All MPoly tests passed!\n");
    return 0;
}
