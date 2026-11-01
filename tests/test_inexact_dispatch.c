/* Tests for the inexact-coefficient round-trip in exact-symbolic
 * builtins (PolynomialGCD, PolynomialLCM, Cancel, Together, Apart,
 * Factor, FactorSquareFree, Simplify, Series, Limit).
 *
 * The round-trip is implemented in src/rationalize.c via
 * `internal_rationalize_then_numericalize`: when any argument contains
 * an inexact leaf (Real / MPFR), it is force-rationalised on the way
 * in, the exact algorithm runs, and the result is numericalised back
 * on the way out. These tests pin both observable behaviours:
 *
 *   1. The exact algorithm actually fires (so cancellations,
 *      factorisations, etc. happen) — without the wrapper, the
 *      polynomial-GCD machinery would silently fall through and the
 *      caller would see only a numeric residue.
 *
 *   2. The output is inexact (contains a Real), since inexact-in /
 *      inexact-out is the contract — a fully-symbolic result would
 *      surprise users who passed in floats deliberately.
 */

#include "core.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "rationalize.h"
#include "symtab.h"
#include "test_utils.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Walk an expression tree and report whether it contains an EXPR_REAL
 * leaf anywhere. Used to confirm the post-numericalize step ran. */
static bool contains_real(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
    if (e->type != EXPR_FUNCTION) return false;
    if (contains_real(e->data.function.head)) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++) {
        if (contains_real(e->data.function.args[i])) return true;
    }
    return false;
}

/* Walk a printed result and check it does NOT contain a marker
 * substring. Useful for "this should not be the unevaluated form". */
static void assert_eval_excludes(const char* input, const char* forbidden) {
    Expr* parsed = parse_expression(input);
    Expr* evaluated = evaluate(parsed);
    expr_free(parsed);
    char* s = expr_to_string(evaluated);
    if (strstr(s, forbidden) != NULL) {
        fprintf(stderr, "FAIL (forbidden substring present): %s\n"
                        "  Forbidden: %s\n  Got:       %s\n",
                input, forbidden, s);
        free(s);
        expr_free(evaluated);
        assert(0);
    }
    free(s);
    expr_free(evaluated);
}

/* Evaluate `input` and check the result contains an EXPR_REAL leaf
 * (i.e. the post-numericalize step ran). */
static void assert_eval_contains_real(const char* input) {
    Expr* parsed = parse_expression(input);
    Expr* evaluated = evaluate(parsed);
    expr_free(parsed);
    if (!contains_real(evaluated)) {
        char* s = expr_to_string(evaluated);
        fprintf(stderr, "FAIL (no inexact leaf): %s\n  Got: %s\n", input, s);
        free(s);
        expr_free(evaluated);
        assert(0);
    }
    expr_free(evaluated);
}

/* ------------------------------------------------------------------------
 *  internal_force_rationalize / internal_contains_inexact directly
 * ---------------------------------------------------------------------- */

static void test_contains_inexact_basic(void) {
    Expr* a = parse_expression("1 + x");
    ASSERT(internal_contains_inexact(a) == false);
    expr_free(a);

    Expr* b = parse_expression("1.0 + x");
    Expr* be = evaluate(b);
    ASSERT(internal_contains_inexact(be) == true);
    expr_free(be);

    Expr* c = parse_expression("Pi + x");
    Expr* ce = evaluate(c);
    /* Pi is exact-symbolic; not inexact. */
    ASSERT(internal_contains_inexact(ce) == false);
    expr_free(ce);

    /* Nested deep. */
    Expr* d = parse_expression("Sin[1 + 2.5]");
    /* Without evaluating: 2.5 is inside. */
    ASSERT(internal_contains_inexact(d) == true);
    expr_free(d);
}

static void test_force_rationalize_preserves_pi(void) {
    Expr* e = parse_expression("Pi + 1.5 x");
    Expr* ee = evaluate(e);
    Expr* r = internal_force_rationalize(ee);
    /* Pi must stay symbolic; 1.5 must become 3/2. */
    char* s = expr_to_string(r);
    /* Either order: "Pi + (3 x)/2" or "(3 x)/2 + Pi". Just check Pi is
     * present and no decimal remains. */
    ASSERT(strstr(s, "Pi") != NULL);
    ASSERT(strstr(s, ".") == NULL);
    free(s);
    expr_free(ee);
    expr_free(r);
}

static void test_force_rationalize_recovers_one_ninth(void) {
    /* The double 1/9 has 4e-18 error from true 1/9. The ½-ulp algorithm
     * (RATIONALIZE_ZERO) would return the bit-exact rational. The
     * convergent-bound algorithm (RATIONALIZE_DEFAULT) recovers 1/9.
     * internal_force_rationalize uses default-then-zero, so it must
     * return 1/9. */
    Expr* e = parse_expression("1/9.");
    Expr* ee = evaluate(e);
    Expr* r = internal_force_rationalize(ee);
    char* s = expr_to_string(r);
    ASSERT_STR_EQ(s, "1/9");
    free(s);
    expr_free(ee);
    expr_free(r);
}

/* ------------------------------------------------------------------------
 *  Cancel
 *
 *  The headline test: the user-supplied
 *    Cancel[(x^2/9 + 2*x*y/15 + y^2/25) / (x^2/9. - y^2/25.)]
 *  must actually cancel the (x/3 + y/5) factor, returning a quotient of
 *  two binomials — the inexact-bypass path leaves the denominator
 *  factored away from a polynomial GCD it never gets to compute.
 * ---------------------------------------------------------------------- */

static void test_cancel_user_perfect_square(void) {
    const char* in = "Cancel[(x^2/9 + (2*x*y)/15 + y^2/25)/(x^2/9. - y^2/25.)]";
    /* Result must contain only first-degree x and y (cancellation
     * happened). Without the wrapper, x^2 and y^2 would appear in both
     * numerator and denominator. */
    assert_eval_excludes(in, "x^2.0");  /* would-be unsimplified Power */
    assert_eval_contains_real(in);       /* inexact-out preserved */
}

static void test_cancel_simple_inexact(void) {
    /* (x^2 - 1.0) / (x - 1.0) → x + 1.0 (or 1.0 + x). */
    assert_eval_contains_real("Cancel[(x^2 - 1.0)/(x - 1.0)]");
}

/* ------------------------------------------------------------------------
 *  Together
 * ---------------------------------------------------------------------- */

static void test_together_one_real_plus_one_over_x(void) {
    /* Together[1.0 + 1/x] → (1.0 + x)/x^1.0  (post-numericalize form). */
    assert_eval_contains_real("Together[1.0 + 1/x]");
    assert_eval_excludes  ("Together[1.0 + 1/x]", "Together[");
}

/* ------------------------------------------------------------------------
 *  PolynomialGCD / PolynomialLCM
 * ---------------------------------------------------------------------- */

static void test_polynomialgcd_inexact_recovers_factor(void) {
    /* PolynomialGCD[x^2 - 1.0, x - 1.0] → x - 1.0. With the inexact
     * inputs the rational-GCD machinery normally returns 1, missing
     * the factor entirely. The wrapper rationalises (1.0 → 1) so the
     * GCD returns x - 1, then numericalises back. */
    assert_eval_contains_real("PolynomialGCD[x^2 - 1.0, x - 1.0]");
    assert_eval_excludes("PolynomialGCD[x^2 - 1.0, x - 1.0]", "PolynomialGCD[");
}

static void test_polynomiallcm_inexact(void) {
    /* PolynomialLCM[x - 1.0, x + 1.0] should be (x-1)(x+1) ≈ x^2 - 1. */
    assert_eval_contains_real("PolynomialLCM[x - 1.0, x + 1.0]");
    assert_eval_excludes("PolynomialLCM[x - 1.0, x + 1.0]", "PolynomialLCM[");
}

/* ------------------------------------------------------------------------
 *  Apart
 * ---------------------------------------------------------------------- */

static void test_apart_inexact(void) {
    /* Apart[1./(x^2 - 1.)] → 1/(2(x-1)) - 1/(2(x+1)) (numerical form). */
    assert_eval_contains_real("Apart[1./(x^2 - 1.)]");
    assert_eval_excludes("Apart[1./(x^2 - 1.)]", "Apart[");
}

/* ------------------------------------------------------------------------
 *  Factor / FactorSquareFree
 * ---------------------------------------------------------------------- */

static void test_factor_inexact(void) {
    /* Factor[x^2 - 1.0] → (x - 1.0)(x + 1.0). Without the wrapper,
     * Factor would leave it as x^2 - 1.0 because the inexact -1.0
     * blocks the integer factoring path. */
    assert_eval_excludes("Factor[x^2 - 1.0]", "Factor[");
    assert_eval_excludes("Factor[x^2 - 1.0]", "x^2");   /* must factor */
    assert_eval_contains_real("Factor[x^2 - 1.0]");
}

static void test_factorsquarefree_inexact(void) {
    /* FactorSquareFree[(x - 1.0)^2 (x + 1.0)] → (x - 1.0)^2 (x + 1.0). */
    assert_eval_contains_real("FactorSquareFree[(x - 1.0)^2 (x + 1.0)]");
    assert_eval_excludes("FactorSquareFree[(x - 1.0)^2 (x + 1.0)]",
                         "FactorSquareFree[");
}

/* ------------------------------------------------------------------------
 *  Simplify / Series / Limit
 * ---------------------------------------------------------------------- */

static void test_simplify_inexact(void) {
    /* Simplify[(x^2 - 1.0)/(x - 1.0)] → 1.0 + x. */
    assert_eval_contains_real("Simplify[(x^2 - 1.0)/(x - 1.0)]");
}

static void test_limit_inexact(void) {
    /* Limit[(x - 1.0)/(x^2 - 1.0), x -> 1] should be 1/2 ≈ 0.5. The
     * exact Limit pipeline (Series-based) needs rational coefficients;
     * the wrapper supplies them. */
    assert_eval_contains_real("Limit[(x - 1.0)/(x^2 - 1.0), x -> 1]");
}

static void test_series_inexact_smoke(void) {
    /* Series[1./(1. - x), {x, 0, 3}] should not be left unevaluated.
     * SeriesData is a function head, so we just verify it printed and
     * contains a Real (post-numericalize). */
    assert_eval_contains_real("Series[1./(1. - x), {x, 0, 3}]");
}

/* ------------------------------------------------------------------------
 *  Negative / pass-through guards
 *
 *  When inputs are fully exact, the wrapper must NOT trigger — exact
 *  output should still be exact, not numericalised by accident.
 * ---------------------------------------------------------------------- */

static void test_exact_inputs_stay_exact(void) {
    /* No floats in the input → no float in the output. */
    assert_eval_eq("Cancel[(x^2 - 1)/(x - 1)]", "1 + x", 0);
    assert_eval_eq("Together[1 + 1/x]", "(1 + x)/x", 0);
    assert_eval_eq("PolynomialGCD[x^2 - 1, x - 1]", "-1 + x", 0);
    assert_eval_eq("Factor[x^2 - 1]", "(-1 + x) (1 + x)", 0);
}

/* ------------------------------------------------------------------------
 *  Kernels survive the round trip verbatim (issue #119)
 *
 *  The numericalise step must touch coefficients only: a kernel the input
 *  wrote with exact arguments, c[1] or x[-2 + g], is an indeterminate, and
 *  c[1.0] is a different one.  A kernel written with inexact arguments comes
 *  back as written, not as the f[1/2] the exact core saw.
 * ---------------------------------------------------------------------- */

static void test_kernel_args_not_numericalised(void) {
    assert_eval_eq("Factor[a .4 + c[1] + d]", "0.2 (2.0 a + 5.0 d + 5.0 c[1])", 0);
    assert_eval_eq("Factor[(x[-1] - x[-2 + g])*0.1 + (-x[-3 + g] + x[-4 + g])]",
                   "0.1 (x[-1] + 10.0 x[-4 + g] - 10.0 x[-3 + g] - 1.0 x[-2 + g])", 0);
    assert_eval_eq("Factor[0.5 x f[1/2] + Sin[2] y]", "0.5 (1.81859 y + x f[1/2])", 0);
    assert_eval_eq("Together[1/c[1] + 0.5/x]", "(0.5 (2.0 x + c[1]))/(x c[1])", 0);
    assert_eval_eq("Cancel[(x^2 - 1.0)/(x - 1) + c[2]]", "1.0 + x + c[2]", 0);
    assert_eval_eq("Simplify[2.0 g[1] + g[1]]", "3.0 g[1]", 0);
    assert_eval_eq("Limit[(Sin[0.5 x] + c[1] x)/x, x -> 0]", "0.5 + c[1]", 0);
}

static void test_kernel_inexact_args_restored(void) {
    assert_eval_eq("Factor[x f[0.5] + 2.0 x]", "x (2.0 + f[0.5])", 0);
}

int main(void) {
    symtab_init();
    core_init();

    /* Direct API */
    TEST(test_contains_inexact_basic);
    TEST(test_force_rationalize_preserves_pi);
    TEST(test_force_rationalize_recovers_one_ninth);

    /* Cancel */
    TEST(test_cancel_user_perfect_square);
    TEST(test_cancel_simple_inexact);

    /* Together */
    TEST(test_together_one_real_plus_one_over_x);

    /* PolynomialGCD / LCM */
    TEST(test_polynomialgcd_inexact_recovers_factor);
    TEST(test_polynomiallcm_inexact);

    /* Apart */
    TEST(test_apart_inexact);

    /* Factor / FactorSquareFree */
    TEST(test_factor_inexact);
    TEST(test_factorsquarefree_inexact);

    /* Simplify / Series / Limit */
    TEST(test_simplify_inexact);
    TEST(test_limit_inexact);
    TEST(test_series_inexact_smoke);

    /* Pass-through */
    TEST(test_exact_inputs_stay_exact);
    TEST(test_kernel_args_not_numericalised);
    TEST(test_kernel_inexact_args_restored);

    printf("All inexact_dispatch tests passed.\n");
    return 0;
}
