/*
 * test_qafactor.c
 * ---------------
 * Unit tests for qafactor (Phase G3) — norm of a polynomial in Q(α)[x]
 * computed as Resultant_y(P_α(y), g(x, y)) ∈ Q[x].
 *
 * Headline tests:
 *   - Norm(x − √2) over Q(√2) = x² − 2 (the conjugate product)
 *   - Norm(x − ∛2) over Q(∛2) = x³ − 2
 *   - Norm(x² + 1) over Q(i) = (x² + 1)² (Q-rational input squared)
 *   - Norm(x² − √2) over Q(√2) = x⁴ − 2
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "expr.h"
#include "core.h"
#include "symtab.h"
#include "qa.h"
#include "qaupoly.h"
#include "qafactor.h"
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

/* Compare an Expr to the FullForm of an expected source string.
 * Both sides are evaluated to canonical form (Plus is Orderless, so
 * argument order is deterministic after evaluate). */
static void assert_expr_eq_str(Expr* actual, const char* expected_src) {
    Expr* parsed   = parse_expression(expected_src);
    ASSERT(parsed != NULL);
    Expr* expected = evaluate(parsed);
    expr_free(parsed);

    /* Re-evaluate `actual` so both sides share canonical Plus order.
     * (qaupoly_norm calls expr_expand, which itself produces canonical
     * form, but the safe move is to run through evaluate too.) */
    Expr* actual_eval = evaluate(actual);

    char* s_actual   = expr_to_string_fullform(actual_eval);
    char* s_expected = expr_to_string_fullform(expected);
    if (strcmp(s_actual, s_expected) != 0) {
        fprintf(stderr, "FAIL: Expr mismatch\n  Expected: %s\n  Actual:   %s\n",
                s_expected, s_actual);
    }
    ASSERT(strcmp(s_actual, s_expected) == 0);

    free(s_actual);
    free(s_expected);
    expr_free(actual_eval);
    expr_free(expected);
}

/* ===================== Construction-helper tests ===================== */

static void test_qaext_to_expr_sqrt2(void) {
    /* P_α(y) = y² − 2 */
    QAExt* ext = qaext_sqrt_si(2);
    Expr* e = qaext_to_expr(ext, "y");
    assert_expr_eq_str(e, "y^2 - 2");
    qaext_free(ext);
}

static void test_qaext_to_expr_qi(void) {
    /* P_α(y) = y² + 1 */
    QAExt* ext = qaext_sqrt_si(-1);
    Expr* e = qaext_to_expr(ext, "y");
    assert_expr_eq_str(e, "y^2 + 1");
    qaext_free(ext);
}

static void test_qaext_to_expr_cbrt2(void) {
    /* P_α(y) = y³ − 2 */
    QAExt* ext = qaext_root_si(2, 3);
    Expr* e = qaext_to_expr(ext, "y");
    assert_expr_eq_str(e, "y^3 - 2");
    qaext_free(ext);
}

static void test_qaupoly_to_expr_x_minus_alpha(void) {
    /* f(x, α) = x − α  in Q(√2)[x]  ⇒  g(x, y) = x − y */
    QAExt* ext = qaext_sqrt_si(2);
    QAUPoly* p = qaupoly_x(ext);                 /* p = x */
    QANum* alpha = qa_alpha(ext);
    QANum* neg_alpha = qa_neg(alpha);
    qaupoly_setcoef(p, 0, neg_alpha);            /* c[0] = −α */
    qa_free(alpha); qa_free(neg_alpha);
    qaupoly_normalize(p);

    Expr* g = qaupoly_to_expr(p, "x", "y");
    assert_expr_eq_str(g, "x - y");
    qaupoly_free(p);
    qaext_free(ext);
}

static void test_qaupoly_to_expr_constant_in_q(void) {
    /* f(x, α) = x² + 1 (no α dependency)  ⇒  g(x, y) = x² + 1 */
    QAExt* ext = qaext_sqrt_si(-1);  /* Q(i), but f is rational */
    QAUPoly* p = qaupoly_new(ext, 3);
    qa_free(p->c[0]); p->c[0] = qa_from_si(ext, 1, 1);
    qa_free(p->c[2]); p->c[2] = qa_from_si(ext, 1, 1);
    p->deg = 2;

    Expr* g = qaupoly_to_expr(p, "x", "y");
    assert_expr_eq_str(g, "x^2 + 1");
    qaupoly_free(p);
    qaext_free(ext);
}

/* ============================== Norm tests =============================== */

/* Norm(x − √2) over Q(√2)
 * = Resultant_y(y² − 2, x − y)
 * = (x − √2)(x + √2) = x² − 2. */
static void test_norm_x_minus_sqrt2(void) {
    QAExt* ext = qaext_sqrt_si(2);
    QAUPoly* p = qaupoly_x(ext);
    QANum* alpha = qa_alpha(ext);
    QANum* neg_alpha = qa_neg(alpha);
    qaupoly_setcoef(p, 0, neg_alpha);
    qa_free(alpha); qa_free(neg_alpha);
    qaupoly_normalize(p);

    Expr* n = qaupoly_norm(p, "x", "y");
    assert_expr_eq_str(n, "x^2 - 2");

    qaupoly_free(p);
    qaext_free(ext);
}

/* Norm(x − ∛2) over Q(∛2)
 * = Resultant_y(y³ − 2, x − y)
 * = (x − ω⁰∛2)(x − ω¹∛2)(x − ω²∛2) = x³ − 2. */
static void test_norm_x_minus_cbrt2(void) {
    QAExt* ext = qaext_root_si(2, 3);
    QAUPoly* p = qaupoly_x(ext);
    QANum* alpha = qa_alpha(ext);
    QANum* neg_alpha = qa_neg(alpha);
    qaupoly_setcoef(p, 0, neg_alpha);
    qa_free(alpha); qa_free(neg_alpha);
    qaupoly_normalize(p);

    Expr* n = qaupoly_norm(p, "x", "y");
    assert_expr_eq_str(n, "x^3 - 2");

    qaupoly_free(p);
    qaext_free(ext);
}

/* Norm(x² + 1) over Q(i) where x² + 1 is rational (no α).
 * For a degree-d extension, Norm of a Q-rational element is its dth power.
 * So Norm(x² + 1) = (x² + 1)² = x⁴ + 2x² + 1. */
static void test_norm_x2_plus_1_in_qi(void) {
    QAExt* ext = qaext_sqrt_si(-1);
    QAUPoly* p = qaupoly_new(ext, 3);
    qa_free(p->c[0]); p->c[0] = qa_from_si(ext, 1, 1);
    qa_free(p->c[2]); p->c[2] = qa_from_si(ext, 1, 1);
    p->deg = 2;

    Expr* n = qaupoly_norm(p, "x", "y");
    assert_expr_eq_str(n, "x^4 + 2 x^2 + 1");

    qaupoly_free(p);
    qaext_free(ext);
}

/* Norm(x² − √2) over Q(√2)
 * = Resultant_y(y² − 2, x² − y)
 * = (x² − √2)(x² + √2) = x⁴ − 2. */
static void test_norm_x2_minus_sqrt2(void) {
    QAExt* ext = qaext_sqrt_si(2);
    QAUPoly* p = qaupoly_new(ext, 3);
    /* c[0] = −α */
    QANum* alpha = qa_alpha(ext);
    QANum* neg_alpha = qa_neg(alpha);
    qa_free(p->c[0]); p->c[0] = neg_alpha;
    qa_free(alpha);   /* note: neg_alpha is owned by p->c[0] now */
    /* c[2] = 1 */
    qa_free(p->c[2]); p->c[2] = qa_from_si(ext, 1, 1);
    p->deg = 2;

    Expr* n = qaupoly_norm(p, "x", "y");
    assert_expr_eq_str(n, "x^4 - 2");

    qaupoly_free(p);
    qaext_free(ext);
}

/* Norm(x² + α x + 1) over Q(√2) where α = √2.
 * = Resultant_y(y² − 2, x² + y x + 1)
 *   The resultant of a degree-2 and degree-2 in y is the determinant of
 *   the 4×4 Sylvester matrix.  Hand calculation:
 *     (x² + √2 x + 1)(x² − √2 x + 1) = (x² + 1)² − 2x²
 *                                     = x⁴ + 2x² + 1 − 2x²
 *                                     = x⁴ + 1.
 *   So Norm = x⁴ + 1 (the cyclotomic Φ_8). */
static void test_norm_x2_plus_alphax_plus_1(void) {
    QAExt* ext = qaext_sqrt_si(2);
    QAUPoly* p = qaupoly_new(ext, 3);
    qa_free(p->c[0]); p->c[0] = qa_from_si(ext, 1, 1);
    qa_free(p->c[1]); p->c[1] = qa_alpha(ext);
    qa_free(p->c[2]); p->c[2] = qa_from_si(ext, 1, 1);
    p->deg = 2;

    Expr* n = qaupoly_norm(p, "x", "y");
    assert_expr_eq_str(n, "x^4 + 1");

    qaupoly_free(p);
    qaext_free(ext);
}

/* Norm of a constant c ∈ Q(α): N(c) = c^deg(P_α). */
static void test_norm_rational_constant(void) {
    /* p = 3 (a Q-rational constant) over Q(√2): Norm = 3² = 9. */
    QAExt* ext = qaext_sqrt_si(2);
    QAUPoly* p = qaupoly_from_si(ext, 3, 1);

    Expr* n = qaupoly_norm(p, "x", "y");
    assert_expr_eq_str(n, "9");

    qaupoly_free(p);
    qaext_free(ext);
}

/* ============================== Driver =============================== */

int main(void) {
    symtab_init();
    core_init();

    printf("Running qafactor (Phase G3 — norm via resultant) tests...\n");

    TEST(test_qaext_to_expr_sqrt2);
    TEST(test_qaext_to_expr_qi);
    TEST(test_qaext_to_expr_cbrt2);
    TEST(test_qaupoly_to_expr_x_minus_alpha);
    TEST(test_qaupoly_to_expr_constant_in_q);

    TEST(test_norm_x_minus_sqrt2);
    TEST(test_norm_x_minus_cbrt2);
    TEST(test_norm_x2_plus_1_in_qi);
    TEST(test_norm_x2_minus_sqrt2);
    TEST(test_norm_x2_plus_alphax_plus_1);
    TEST(test_norm_rational_constant);

    printf("All qafactor tests passed!\n");
    return 0;
}
