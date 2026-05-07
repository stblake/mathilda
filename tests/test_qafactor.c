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

/* ===================== Phase G4 — sqfr_norm + alg_factor ================ */

/* Helper: build f(x) ∈ Q(α)[x] from an array of α-degree-0 rational
 * coefficients (constant Q-rational coefficients).  Caller frees
 * the returned QAUPoly. */
static QAUPoly* qaupoly_from_si_coefs(const QAExt* ext,
                                      const long* nums,
                                      int deg_plus_1) {
    QAUPoly* p = qaupoly_new(ext, deg_plus_1);
    for (int i = 0; i < deg_plus_1; i++) {
        QANum* qn = qa_from_si(ext, nums[i], 1);
        qaupoly_setcoef(p, i, qn);
        qa_free(qn);
    }
    qaupoly_normalize(p);
    return p;
}

/* Multiply out an array of monic QAUPoly factors and compare to a
 * given input.  Returns true iff prod == expected (over Q(α)[x]). */
static bool product_equals(QAUPoly** factors, int n, const QAUPoly* expected) {
    if (n <= 0) return false;
    QAUPoly* prod = qaupoly_copy(factors[0]);
    for (int i = 1; i < n; i++) {
        QAUPoly* tmp = qaupoly_mul(prod, factors[i]);
        qaupoly_free(prod);
        prod = tmp;
    }
    /* expected may not be monic — make both monic for comparison. */
    QAUPoly* expected_monic = qaupoly_make_monic(expected);
    bool eq = qaupoly_eq(prod, expected_monic);
    qaupoly_free(prod);
    qaupoly_free(expected_monic);
    return eq;
}

/* Convert a QAUPoly to an Expr via qaupoly_to_expr (with α as the
 * second free variable named "y") and compare to a string source.
 * The string source is parsed and evaluated with the symbol "y"
 * standing for α (so e.g. expected == "x - y" is the textual form
 * of x - α). */
static void assert_qaupoly_eq_str(const QAUPoly* f, const char* expected_src) {
    Expr* e = qaupoly_to_expr(f, "x", "y");
    assert_expr_eq_str(e, expected_src);
}

/* sqfr_norm of an already-squarefree-norm input: shift count s == 0. */
static void test_sqfr_norm_s_zero_x_minus_sqrt2(void) {
    /* f(x) = x − √2 over Q(√2) — Norm = x² − 2 (squarefree). */
    QAExt* ext = qaext_sqrt_si(2);
    QAUPoly* p = qaupoly_x(ext);
    QANum* alpha = qa_alpha(ext);
    QANum* neg_alpha = qa_neg(alpha);
    qaupoly_setcoef(p, 0, neg_alpha);
    qa_free(alpha); qa_free(neg_alpha);
    qaupoly_normalize(p);

    QASqfrNormResult r = qa_sqfr_norm(p, "x", "y");
    ASSERT(r.s == 0);
    ASSERT(r.g != NULL);
    ASSERT(r.R != NULL);
    /* g unchanged at s = 0, equals f. */
    ASSERT(qaupoly_eq(r.g, p));
    assert_expr_eq_str(expr_copy(r.R), "x^2 - 2");

    qa_sqfr_norm_result_free(&r);
    qaupoly_free(p);
    qaext_free(ext);
}

/* sqfr_norm of a Q-rational polynomial: needs s ≥ 1. */
static void test_sqfr_norm_s_positive_for_q_input(void) {
    /* f(x) = x² − 2 (no α dependency) over Q(√2).  At s = 0 the
     * norm is (x² − 2)² which is NOT squarefree.  sqfr_norm must
     * pick s ≥ 1. */
    QAExt* ext = qaext_sqrt_si(2);
    long coefs[3] = { -2, 0, 1 };
    QAUPoly* p = qaupoly_from_si_coefs(ext, coefs, 3);

    QASqfrNormResult r = qa_sqfr_norm(p, "x", "y");
    ASSERT(r.s >= 1);
    ASSERT(r.s <= 4);   /* in practice s ∈ {1, 2}, give some headroom */
    ASSERT(r.g != NULL);
    ASSERT(r.R != NULL);

    qa_sqfr_norm_result_free(&r);
    qaupoly_free(p);
    qaext_free(ext);
}

/* alg_factor: x² − 2 over Q(√2) ⇒ (x − √2)(x + √2). */
static void test_alg_factor_x2_minus_2_in_qsqrt2(void) {
    QAExt* ext = qaext_sqrt_si(2);
    long coefs[3] = { -2, 0, 1 };
    QAUPoly* f = qaupoly_from_si_coefs(ext, coefs, 3);

    int n;
    QAUPoly** facs = qa_alg_factor(f, "x", "y", &n);
    ASSERT(facs != NULL);
    ASSERT(n == 2);

    /* Each factor monic of degree 1.  Verify by product == f. */
    for (int i = 0; i < n; i++) {
        ASSERT(facs[i]->deg == 1);
    }
    ASSERT(product_equals(facs, n, f));

    /* Spot-check: each factor is x ± α (i.e. x ± √2). */
    /* Factor identity:  facs is unordered.  Build the textual set. */
    bool saw_minus = false, saw_plus = false;
    for (int i = 0; i < n; i++) {
        Expr* e = qaupoly_to_expr(facs[i], "x", "y");
        char* s = expr_to_string_fullform(e);
        Expr* parsed_minus = parse_expression("x - y");
        Expr* exp_minus = evaluate(parsed_minus); expr_free(parsed_minus);
        Expr* parsed_plus  = parse_expression("x + y");
        Expr* exp_plus  = evaluate(parsed_plus);  expr_free(parsed_plus);
        Expr* eval_e = evaluate(expr_copy(e));
        char* sm = expr_to_string_fullform(exp_minus);
        char* sp = expr_to_string_fullform(exp_plus);
        char* se = expr_to_string_fullform(eval_e);
        if (strcmp(se, sm) == 0) saw_minus = true;
        if (strcmp(se, sp) == 0) saw_plus  = true;
        free(sm); free(sp); free(se);
        free(s); expr_free(e); expr_free(eval_e);
        expr_free(exp_minus); expr_free(exp_plus);
    }
    ASSERT(saw_minus);
    ASSERT(saw_plus);

    for (int i = 0; i < n; i++) qaupoly_free(facs[i]);
    free(facs);
    qaupoly_free(f);
    qaext_free(ext);
}

/* alg_factor: x⁴ + 1 over Q(√2) ⇒ (x² − √2 x + 1)(x² + √2 x + 1). */
static void test_alg_factor_x4_plus_1_in_qsqrt2(void) {
    QAExt* ext = qaext_sqrt_si(2);
    long coefs[5] = { 1, 0, 0, 0, 1 };
    QAUPoly* f = qaupoly_from_si_coefs(ext, coefs, 5);

    int n;
    QAUPoly** facs = qa_alg_factor(f, "x", "y", &n);
    ASSERT(facs != NULL);
    ASSERT(n == 2);
    for (int i = 0; i < n; i++) {
        ASSERT(facs[i]->deg == 2);
    }
    ASSERT(product_equals(facs, n, f));

    for (int i = 0; i < n; i++) qaupoly_free(facs[i]);
    free(facs);
    qaupoly_free(f);
    qaext_free(ext);
}

/* alg_factor: x⁴ − 5x² + 6 = (x² − 2)(x² − 3) over Q(√2)
 *            ⇒ (x − √2)(x + √2)(x² − 3). */
static void test_alg_factor_x4_5x2_6_in_qsqrt2(void) {
    QAExt* ext = qaext_sqrt_si(2);
    long coefs[5] = { 6, 0, -5, 0, 1 };
    QAUPoly* f = qaupoly_from_si_coefs(ext, coefs, 5);

    int n;
    QAUPoly** facs = qa_alg_factor(f, "x", "y", &n);
    ASSERT(facs != NULL);
    ASSERT(n == 3);

    /* Two linear, one quadratic. */
    int deg1_count = 0, deg2_count = 0;
    for (int i = 0; i < n; i++) {
        if (facs[i]->deg == 1) deg1_count++;
        else if (facs[i]->deg == 2) deg2_count++;
    }
    ASSERT(deg1_count == 2);
    ASSERT(deg2_count == 1);
    ASSERT(product_equals(facs, n, f));

    for (int i = 0; i < n; i++) qaupoly_free(facs[i]);
    free(facs);
    qaupoly_free(f);
    qaext_free(ext);
}

/* alg_factor: x³ − 2 over Q(∛2)
 *            ⇒ (x − ∛2)(x² + ∛2 x + ∛4). */
static void test_alg_factor_x3_minus_2_in_qcbrt2(void) {
    QAExt* ext = qaext_root_si(2, 3);
    long coefs[4] = { -2, 0, 0, 1 };
    QAUPoly* f = qaupoly_from_si_coefs(ext, coefs, 4);

    int n;
    QAUPoly** facs = qa_alg_factor(f, "x", "y", &n);
    ASSERT(facs != NULL);
    ASSERT(n == 2);

    int deg1_count = 0, deg2_count = 0;
    for (int i = 0; i < n; i++) {
        if (facs[i]->deg == 1) deg1_count++;
        else if (facs[i]->deg == 2) deg2_count++;
    }
    ASSERT(deg1_count == 1);
    ASSERT(deg2_count == 1);
    ASSERT(product_equals(facs, n, f));

    for (int i = 0; i < n; i++) qaupoly_free(facs[i]);
    free(facs);
    qaupoly_free(f);
    qaext_free(ext);
}

/* alg_factor: x² + 1 over Q(i) ⇒ (x − i)(x + i). */
static void test_alg_factor_x2_plus_1_in_qi(void) {
    QAExt* ext = qaext_sqrt_si(-1);
    long coefs[3] = { 1, 0, 1 };
    QAUPoly* f = qaupoly_from_si_coefs(ext, coefs, 3);

    int n;
    QAUPoly** facs = qa_alg_factor(f, "x", "y", &n);
    ASSERT(facs != NULL);
    ASSERT(n == 2);
    for (int i = 0; i < n; i++) {
        ASSERT(facs[i]->deg == 1);
    }
    ASSERT(product_equals(facs, n, f));

    for (int i = 0; i < n; i++) qaupoly_free(facs[i]);
    free(facs);
    qaupoly_free(f);
    qaext_free(ext);
}

/* alg_factor: x² + x + 1 over Q(√−3) (the Eisenstein extension)
 *            ⇒ (x − (−1 + √−3)/2)(x − (−1 − √−3)/2). */
static void test_alg_factor_x2_x_1_in_qsqrtm3(void) {
    QAExt* ext = qaext_sqrt_si(-3);
    long coefs[3] = { 1, 1, 1 };
    QAUPoly* f = qaupoly_from_si_coefs(ext, coefs, 3);

    int n;
    QAUPoly** facs = qa_alg_factor(f, "x", "y", &n);
    ASSERT(facs != NULL);
    ASSERT(n == 2);
    for (int i = 0; i < n; i++) {
        ASSERT(facs[i]->deg == 1);
    }
    ASSERT(product_equals(facs, n, f));

    for (int i = 0; i < n; i++) qaupoly_free(facs[i]);
    free(facs);
    qaupoly_free(f);
    qaext_free(ext);
}

/* Irreducible-over-Q(α): the input is already irreducible.  Expect a
 * single-factor output equal to the (made-monic) input. */
static void test_alg_factor_irreducible_input(void) {
    /* f(x) = x² − 3 over Q(√2) is irreducible (√3 ∉ Q(√2)). */
    QAExt* ext = qaext_sqrt_si(2);
    long coefs[3] = { -3, 0, 1 };
    QAUPoly* f = qaupoly_from_si_coefs(ext, coefs, 3);

    int n;
    QAUPoly** facs = qa_alg_factor(f, "x", "y", &n);
    ASSERT(facs != NULL);
    ASSERT(n == 1);
    ASSERT(facs[0]->deg == 2);
    ASSERT(product_equals(facs, n, f));

    qaupoly_free(facs[0]);
    free(facs);
    qaupoly_free(f);
    qaext_free(ext);
}

/* Suppress unused-warning for the helper assert_qaupoly_eq_str
 * (kept exposed in case future tests want textual factor checks). */
static void _unused_keep_helpers(void) {
    (void)assert_qaupoly_eq_str;
}

/* ============================ Phase G5 ============================ */
/* End-to-end picocas-API tests: parse `Factor[poly, Extension -> α]`,
 * evaluate, and compare the result (as FullForm) to the expected
 * factored form parsed and evaluated through the same pipeline. */

static void assert_factor_eq_str(const char* input_src,
                                 const char* expected_src) {
    Expr* in_parsed = parse_expression(input_src);
    ASSERT(in_parsed != NULL);
    Expr* actual = evaluate(in_parsed);
    expr_free(in_parsed);
    assert_expr_eq_str(actual, expected_src);
}

/* x² − 2 over Q(√2) → (x − √2)(x + √2). */
static void test_g5_x2_minus_2_qsqrt2(void) {
    assert_factor_eq_str("Factor[x^2 - 2, Extension -> Sqrt[2]]",
                         "(x - Sqrt[2]) (x + Sqrt[2])");
}

/* x⁴ + 1 over Q(√2) → (x² − √2 x + 1)(x² + √2 x + 1). */
static void test_g5_x4_plus_1_qsqrt2(void) {
    assert_factor_eq_str("Factor[x^4 + 1, Extension -> Sqrt[2]]",
                         "(x^2 - Sqrt[2] x + 1) (x^2 + Sqrt[2] x + 1)");
}

/* x⁴ − 5x² + 6 over Q(√2) → (x − √2)(x + √2)(x² − 3). */
static void test_g5_x4_5x2_6_qsqrt2(void) {
    assert_factor_eq_str("Factor[x^4 - 5 x^2 + 6, Extension -> Sqrt[2]]",
                         "(x - Sqrt[2]) (x + Sqrt[2]) (x^2 - 3)");
}

/* x³ − 2 over Q(∛2) → (x − ∛2)(x² + ∛2 x + ∛4). */
static void test_g5_x3_minus_2_qcbrt2(void) {
    assert_factor_eq_str("Factor[x^3 - 2, Extension -> 2^(1/3)]",
                         "(x - 2^(1/3)) (x^2 + 2^(1/3) x + 2^(2/3))");
}

/* x² + 1 over Q(i) → (x − I)(x + I). */
static void test_g5_x2_plus_1_qi(void) {
    assert_factor_eq_str("Factor[x^2 + 1, Extension -> I]",
                         "(x - I) (x + I)");
}

/* x² + 4 over Q(i) → (x − 2I)(x + 2I) — coefficient handling. */
static void test_g5_x2_plus_4_qi(void) {
    assert_factor_eq_str("Factor[x^2 + 4, Extension -> I]",
                         "(x - 2 I) (x + 2 I)");
}

/* x⁴ + 4 over Q(i) — degree-2 extension, Sophie-Germain identity:
 * x⁴ + 4 = (x² + 2x + 2)(x² − 2x + 2) over Q, but each of those
 * further factors over Q(i): (x² + 2x + 2) = (x + 1 − i)(x + 1 + i),
 * etc.  Final form has 4 linear factors. */
static void test_g5_x4_plus_4_qi(void) {
    assert_factor_eq_str(
        "Factor[x^4 + 4, Extension -> I]",
        "(x - 1 - I) (x - 1 + I) (x + 1 - I) (x + 1 + I)");
}

/* x² + x + 1 over Q(√−3) → (x − ω)(x − ω̄) where ω = (−1 + I√3)/2.
 * Picocas renders this as (x + 1/2 + I√3/2)(x + 1/2 − I√3/2). */
static void test_g5_x2_x_1_qsqrtm3(void) {
    assert_factor_eq_str(
        "Factor[x^2 + x + 1, Extension -> Sqrt[-3]]",
        "(x + 1/2 - I Sqrt[3]/2) (x + 1/2 + I Sqrt[3]/2)");
}

/* α-bearing input: x² − 2√2 x + 2 = (x − √2)² over Q(√2).  Tests the
 * recognition of `Sqrt[2]` as α inside the input polynomial and the
 * multiplicity-detection trial-division loop. */
static void test_g5_alpha_in_input(void) {
    assert_factor_eq_str(
        "Factor[x^2 - 2 Sqrt[2] x + 2, Extension -> Sqrt[2]]",
        "(x - Sqrt[2])^2");
}

/* Repeated-factor input: (x² − 2)² over Q(√2) → (x − √2)²(x + √2)². */
static void test_g5_repeated_factor_qsqrt2(void) {
    assert_factor_eq_str(
        "Factor[(x^2 - 2)^2, Extension -> Sqrt[2]]",
        "(x - Sqrt[2])^2 (x + Sqrt[2])^2");
}

/* Irreducible over the supplied extension: x² − 3 over Q(√2)
 * stays unfactored as x² − 3 (the extension does not contain √3). */
static void test_g5_irreducible_in_qsqrt2(void) {
    assert_factor_eq_str(
        "Factor[x^2 - 3, Extension -> Sqrt[2]]",
        "x^2 - 3");
}

/* ============================== Driver =============================== */

int main(void) {
    symtab_init();
    core_init();

    printf("Running qafactor (Phases G3 + G4) tests...\n");

    /* Phase G3 */
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

    /* Phase G4 — sqfr_norm */
    TEST(test_sqfr_norm_s_zero_x_minus_sqrt2);
    TEST(test_sqfr_norm_s_positive_for_q_input);

    /* Phase G4 — alg_factor headlines (§14.2 rows 1-6) */
    TEST(test_alg_factor_x2_minus_2_in_qsqrt2);
    TEST(test_alg_factor_x4_plus_1_in_qsqrt2);
    TEST(test_alg_factor_x4_5x2_6_in_qsqrt2);
    TEST(test_alg_factor_x3_minus_2_in_qcbrt2);
    TEST(test_alg_factor_x2_plus_1_in_qi);
    TEST(test_alg_factor_x2_x_1_in_qsqrtm3);

    /* Phase G4 — irreducible-input edge case */
    TEST(test_alg_factor_irreducible_input);

    /* Phase G5 — picocas-level Factor[..., Extension -> α] API */
    TEST(test_g5_x2_minus_2_qsqrt2);
    TEST(test_g5_x4_plus_1_qsqrt2);
    TEST(test_g5_x4_5x2_6_qsqrt2);
    TEST(test_g5_x3_minus_2_qcbrt2);
    TEST(test_g5_x2_plus_1_qi);
    TEST(test_g5_x2_plus_4_qi);
    TEST(test_g5_x4_plus_4_qi);
    TEST(test_g5_x2_x_1_qsqrtm3);
    TEST(test_g5_alpha_in_input);
    TEST(test_g5_repeated_factor_qsqrt2);
    TEST(test_g5_irreducible_in_qsqrt2);

    _unused_keep_helpers();

    printf("All qafactor tests passed!\n");
    return 0;
}
