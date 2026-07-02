/*
 * test_rootofunity.c
 * ------------------
 * Unit tests for the cyclotomic-extension primitives used by the Q(α)
 * extension recognizer (src/poly/rootofunity.{c,h}):
 *
 *   - qaext_cyclotomic(n): the field Q(ζ_n) with minimal polynomial
 *     Φ_n(y), verified against the known low-order cyclotomic polynomials.
 *   - expr_is_root_of_unity_pow(e): recognition of Mathilda's canonical
 *     root-of-unity surface forms Power[-1, Rational[p, q]] and
 *     Complex[0, ±1].
 *
 * These are the foundation for Phase 1 of SIMPLIFY_IMPROVEMENT_PLAN.md:
 * making Together/Cancel/Simplify fast over cyclotomic Q(α).
 */

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "rootofunity.h"
#include "qa.h"
#include "expr.h"
#include "sym_names.h"
#include "symtab.h"
#include "core.h"
#include "test_utils.h"

/* Assert Φ_n has the given degree and integer coefficient list
 * coef[0..deg] (coef[i] = coefficient of y^i). */
static void assert_cyclotomic(unsigned long n, size_t deg, const long* coef) {
    QAExt* ext = qaext_cyclotomic(n);
    ASSERT(ext != NULL);
    if (ext->deg != deg) {
        fprintf(stderr, "Phi_%lu: expected degree %zu, got %zu\n",
                n, deg, ext->deg);
        ASSERT(0);
    }
    mpq_t e;
    mpq_init(e);
    for (size_t i = 0; i <= deg; i++) {
        mpq_set_si(e, coef[i], 1);
        mpq_canonicalize(e);
        if (!mpq_equal(ext->coef[i], e)) {
            char* got = mpq_get_str(NULL, 10, ext->coef[i]);
            fprintf(stderr, "Phi_%lu coef[%zu]: expected %ld, got %s\n",
                    n, i, coef[i], got);
            free(got);
            ASSERT(0);
        }
    }
    mpq_clear(e);
    qaext_free(ext);
}

static void test_cyclotomic_polynomials(void) {
    assert_cyclotomic(1,  1, (long[]){ -1, 1 });            /* y − 1        */
    assert_cyclotomic(2,  1, (long[]){  1, 1 });            /* y + 1        */
    assert_cyclotomic(3,  2, (long[]){  1, 1, 1 });         /* y² + y + 1   */
    assert_cyclotomic(4,  2, (long[]){  1, 0, 1 });         /* y² + 1       */
    assert_cyclotomic(5,  4, (long[]){  1, 1, 1, 1, 1 });   /* Φ_5          */
    assert_cyclotomic(6,  2, (long[]){  1, -1, 1 });        /* y² − y + 1   */
    assert_cyclotomic(8,  4, (long[]){  1, 0, 0, 0, 1 });   /* y⁴ + 1       */
    assert_cyclotomic(10, 4, (long[]){  1, -1, 1, -1, 1 }); /* Φ_10         */
    assert_cyclotomic(12, 4, (long[]){  1, 0, -1, 0, 1 });  /* y⁴ − y² + 1  */
}

/* Φ_105 is the first cyclotomic polynomial with a coefficient of
 * magnitude > 1 (the famous −2 at y^7 and y^41). Spot-check degree and
 * one of the −2 coefficients to confirm the exact integer arithmetic. */
static void test_cyclotomic_105(void) {
    QAExt* ext = qaext_cyclotomic(105);
    ASSERT(ext != NULL);
    ASSERT(ext->deg == 48);   /* φ(105) = φ(3·5·7) = 2·4·6 = 48 */
    mpq_t m2;
    mpq_init(m2);
    mpq_set_si(m2, -2, 1);
    ASSERT(mpq_equal(ext->coef[7], m2));
    ASSERT(mpq_equal(ext->coef[41], m2));
    mpq_clear(m2);
    qaext_free(ext);
}

/* Build Power[-1, Rational[p, q]] directly (the canonical form Mathilda
 * produces for (-1)^(p/q)). Ownership of the integers/Rational is adopted
 * by the constructed nodes. */
static Expr* mk_neg1_pow(long p, long q) {
    Expr* rat = expr_new_function(expr_new_symbol(SYM_Rational),
        (Expr*[]){ expr_new_integer(p), expr_new_integer(q) }, 2);
    return expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ expr_new_integer(-1), rat }, 2);
}

static Expr* mk_complex(long re, long im) {
    return expr_new_function(expr_new_symbol(SYM_Complex),
        (Expr*[]){ expr_new_integer(re), expr_new_integer(im) }, 2);
}

static void assert_rou(Expr* e, long want_p, long want_q) {
    int64_t p = 0, q = 0;
    if (!expr_is_root_of_unity_pow(e, &p, &q)) {
        fprintf(stderr, "expected root-of-unity match (%ld/%ld), got none\n",
                want_p, want_q);
        ASSERT(0);
    }
    ASSERT(p == want_p);
    ASSERT(q == want_q);
    expr_free(e);
}

static void assert_not_rou(Expr* e) {
    int64_t p = 0, q = 0;
    ASSERT(!expr_is_root_of_unity_pow(e, &p, &q));
    expr_free(e);
}

static void test_recognize_roots_of_unity(void) {
    assert_rou(mk_neg1_pow(2, 3), 2, 3);   /* Exp[2 Pi I/3] = (-1)^(2/3) */
    assert_rou(mk_neg1_pow(1, 3), 1, 3);   /* (-1)^(1/3) = primitive ζ_6 */
    assert_rou(mk_neg1_pow(1, 5), 1, 5);   /* Exp[I Pi/5] = (-1)^(1/5)   */
    assert_rou(mk_complex(0, 1),  1, 2);   /* I  = (-1)^(1/2)            */
    assert_rou(mk_complex(0, -1), 3, 2);   /* -I = (-1)^(3/2)            */
    assert_rou(expr_new_symbol(SYM_I), 1, 2);
}

static void test_reject_non_roots(void) {
    /* Sqrt[2] = Power[2, 1/2] — a radical, not a root of unity. */
    Expr* rat = expr_new_function(expr_new_symbol(SYM_Rational),
        (Expr*[]){ expr_new_integer(1), expr_new_integer(2) }, 2);
    Expr* sqrt2 = expr_new_function(expr_new_symbol(SYM_Power),
        (Expr*[]){ expr_new_integer(2), rat }, 2);
    assert_not_rou(sqrt2);

    assert_not_rou(expr_new_symbol("x"));
    assert_not_rou(expr_new_integer(7));
    assert_not_rou(mk_complex(1, 1));   /* 1 + I — not ±I */
}

int main(void) {
    symtab_init();
    core_init();

    printf("Running root-of-unity / cyclotomic tests...\n");

    TEST(test_cyclotomic_polynomials);
    TEST(test_cyclotomic_105);
    TEST(test_recognize_roots_of_unity);
    TEST(test_reject_non_roots);

    printf("All root-of-unity tests passed!\n");
    return 0;
}
