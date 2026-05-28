/*
 * test_factor_baseline.c
 * ----------------------
 * Performance and correctness baseline for the multivariate factoring
 * pipeline.  Each case runs a Factor or Simplify, verifies the result
 * matches the expected form, and asserts a generous wall-clock budget.
 *
 * This file is the integration point for the cumulative Factor work
 * (Phases 0-6 + tactical fixes).  When the budgets here become too
 * tight, treat that as a regression signal: either the algorithm
 * changed shape (in which case update the budget after auditing) or
 * something genuinely slowed down (investigate).
 *
 * Budgets are intentionally generous (3-5x typical observed timings)
 * to remain stable across machines.  Reference timings on a 2020 x86
 * Mac (gcc -O3): see the comment by each test.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

#undef ASSERT
#define ASSERT(cond) do { \
    if (!(cond)) { \
        fprintf(stderr, "ASSERTION FAILED: %s\n  at %s:%d\n", \
                #cond, __FILE__, __LINE__); \
        abort(); \
    } \
} while (0)

/* Evaluate `input`, time it, and check that the FullForm of the result
 * matches `expected`.  Returns the elapsed milliseconds. */
static double eval_check(const char* input, const char* expected,
                         double budget_ms) {
    Expr* parsed = parse_expression(input);
    ASSERT(parsed != NULL);
    clock_t t0 = clock();
    Expr* result = evaluate(parsed);
    clock_t t1 = clock();
    double ms = 1000.0 * (double)(t1 - t0) / (double)CLOCKS_PER_SEC;

    char* got = expr_to_string_fullform(result);
    if (strcmp(got, expected) != 0) {
        fprintf(stderr, "FAIL [correctness]: %s\n"
                        "  expected: %s\n"
                        "  got:      %s\n  (%.2f ms)\n",
                input, expected, got, ms);
        free(got); expr_free(result); expr_free(parsed);
        ASSERT(0);
    }
    if (ms > budget_ms) {
        fprintf(stderr, "FAIL [perf]: %s\n"
                        "  budget: %.2f ms, actual: %.2f ms\n",
                input, budget_ms, ms);
        free(got); expr_free(result); expr_free(parsed);
        ASSERT(0);
    }
    free(got); expr_free(result); expr_free(parsed);
    return ms;
}

/* ====================================================================== */
/*  (1) Direct Factor calls -- correctness + perf                         */
/* ====================================================================== */

/* Univariate Factor: should be sub-millisecond. */
static void test_univariate_factors(void) {
    eval_check("Factor[x^2 - 1]",
               "Times[Plus[-1, x], Plus[1, x]]", 50.0);
    eval_check("Factor[x^3 - 6 x^2 + 11 x - 6]",
               "Times[Plus[-3, x], Plus[-2, x], Plus[-1, x]]", 50.0);
    eval_check("Factor[x^10 - 1]",
               "Times[Plus[-1, x], Plus[1, x], "
               "Plus[1, x, Power[x, 2], Power[x, 3], Power[x, 4]], "
               "Plus[1, Times[-1, x], Power[x, 2], Times[-1, Power[x, 3]], Power[x, 4]]]", 100.0);
}

/* Bivariate via factor_binomial / monomial-content (Phase 0). */
static void test_bivariate_simple(void) {
    eval_check("Factor[x^2 - 4 y^2]",
               "Times[Plus[x, Times[-2, y]], Plus[x, Times[2, y]]]", 50.0);
    eval_check("Factor[x^3 + y^3]",
               "Times[Plus[x, y], Plus[Power[x, 2], Times[-1, Times[x, y]], Power[y, 2]]]", 50.0);
    eval_check("Factor[x^2 y + x y^2]",
               "Times[x, y, Plus[x, y]]", 50.0);
}

/* Bivariate via Hensel lift (Phase 5a) -- previously stalled. */
static void test_bivariate_hensel(void) {
    eval_check("Factor[(x + y - 1)(x + y + 1)]",
               "Times[Plus[-1, x, y], Plus[1, x, y]]", 100.0);
    eval_check("Factor[(x^2 + y - 1)(x + y + 2)]",
               "Times[Plus[2, x, y], Plus[-1, Power[x, 2], y]]", 100.0);
    eval_check("Factor[(x - 1)(x + 1)(x + y)]",
               "Times[Plus[-1, x], Plus[1, x], Plus[x, y]]", 100.0);
}

/* Multivariate via factor_roots fallback (still works after Phase 0+1). */
static void test_multivariate_via_factor_roots(void) {
    eval_check("Factor[2 x^3 y - 2 a^2 x y - 3 a^2 x^2 + 3 a^4]",
               "Times[Plus[a, x], Plus[Times[-1, a], x], "
                     "Plus[Times[-3, Power[a, 2]], Times[2, Times[x, y]]]]", 100.0);
}

/* Irreducibility short-circuit (Phase 1) -- 60x speedup vs main. */
static void test_irreducibility_short_circuit(void) {
    eval_check("Factor[x^2 + y^2 + 1]",
               "Plus[1, Power[x, 2], Power[y, 2]]", 50.0);
    eval_check("Factor[x^4 + y^4 + 1]",
               "Plus[1, Power[x, 4], Power[y, 4]]", 50.0);
    eval_check("Factor[3 x^2 - 3 - y^2]",
               "Plus[-3, Times[3, Power[x, 2]], Times[-1, Power[y, 2]]]", 50.0);
}

/* Threading over comparison heads (Phase 5d, 2nd attempt). */
static void test_factor_threads_over_logic_heads(void) {
    /* Cubic inside Less: factor each comparison argument. */
    eval_check("Factor[1 < 1 + 2 x + x^2 + 1/(1+x) < 2]",
               "Inequality[1, Less, "
                    "Times[Plus[2, x], Power[Plus[1, x], -1], Plus[1, x, Power[x, 2]]], "
                    "Less, 2]",
               150.0);
}

/* Denominator factoring with extra variables (Phase 5d, context-aware). */
static void test_factor_denominator_with_extra_vars(void) {
    eval_check("Factor[(x^3+2 x^2)/(x^2-4 y^2) - (x+2)/(x^2-4 y^2)]",
               "Times[Plus[-1, x], Plus[1, x], Plus[2, x], "
                     "Power[Plus[x, Times[-2, y]], -1], "
                     "Power[Plus[x, Times[2, y]], -1]]",
               150.0);
}

/* ====================================================================== */
/*  (2) Simplify integration -- the user-reported workloads               */
/* ====================================================================== */

/* The user's primary case.  Reference: ~478 ms.  Pre-Phase-0 it was
 * ~1158 ms (2.4x improvement). */
static void test_simplify_user_case(void) {
    double ms = eval_check("Simplify[Sin[x]^3 + Sin[3 x] - 3 Sin[x]]",
                           "Times[-3, Power[Sin[x], 3]]",
                           1500.0);
    fprintf(stderr, "  user case: %.0f ms\n", ms);
}

/* Sensitive Tan case: rolled back v1 of 5d caused a 26x slowdown here.
 * After context-aware fix (Phase 5d v2) this stays around 2.5 s,
 * which is faster than pristine main (3.3 s) due to Factor memo. */
static void test_simplify_tan_double_angle(void) {
    double ms = eval_check("Simplify[Tan[2 x] - 2 Tan[x]/(1 - Tan[x]^2)]",
                           "0", 5000.0);
    fprintf(stderr, "  Tan double-angle: %.0f ms\n", ms);
}

/* Pythagorean to-the-fourth.  Reference: ~80 ms. */
static void test_simplify_pythagorean_power(void) {
    double ms = eval_check("Simplify[(Sin[x]^2 + Cos[x]^2)^4]",
                           "1", 500.0);
    fprintf(stderr, "  Pythag^4: %.0f ms\n", ms);
}

/* Hyperbolic Pythagorean.  Should be near-instant. */
static void test_simplify_hyperbolic_pythag(void) {
    double ms = eval_check("Simplify[Sinh[x]^2 - Cosh[x]^2]",
                           "-1", 200.0);
    fprintf(stderr, "  Sinh^2 - Cosh^2: %.0f ms\n", ms);
}

/* Identity that needs Factor working on the cubic 4 Cos^3 - 3 Cos. */
static void test_simplify_chebyshev_cubic(void) {
    double ms = eval_check("Simplify[Cos[3 x]/(4 Cos[x]^3 - 3 Cos[x])]",
                           "1", 500.0);
    fprintf(stderr, "  Cos[3x]/Cheb: %.0f ms\n", ms);
}

/* Compound expansion: (x+y)^3 - (x-y)^3 = 6 x^2 y + 2 y^3.
 * Reference: ~46 ms. */
static void test_simplify_cubic_expansion(void) {
    double ms = eval_check("Simplify[(x + y)^3 - (x - y)^3]",
                           "Times[2, y, Plus[Times[3, Power[x, 2]], Power[y, 2]]]",
                           300.0);
    fprintf(stderr, "  (x+y)^3 - (x-y)^3: %.0f ms\n", ms);
}

/* Sin double angle through division. */
static void test_simplify_sin2x_over_sinx(void) {
    double ms = eval_check("Simplify[Sin[2 x] / Sin[x]]",
                           "Times[2, Cos[x]]", 500.0);
    fprintf(stderr, "  Sin[2x]/Sin[x]: %.0f ms\n", ms);
}

/* ====================================================================== */
/*  (3) Sanity: previously slow / known-tricky cases                     */
/* ====================================================================== */

/* The user's reported pre-Phase-0 case.  Was returning input
 * unfactored AND in 364 ms; now returns the correct factorisation in
 * single-digit ms. */
static void test_factor_monomial_cube(void) {
    double ms = eval_check(
        "Factor[3 a^2 b - 3 b - b^3]",
        "Times[b, Plus[-3, Times[3, Power[a, 2]], Times[-1, Power[b, 2]]]]",
        100.0);
    fprintf(stderr, "  3a^2 b - 3b - b^3: %.0f ms\n", ms);
}

/* Trig form of the same shape, derived from Simplify intermediate. */
static void test_factor_trig_monomial(void) {
    double ms = eval_check(
        "Factor[3 Cos[x]^2 Sin[x] - 3 Sin[x] - Sin[x]^3]",
        "Times[Sin[x], Plus[-3, Times[3, Power[Cos[x], 2]], "
                          "Times[-1, Power[Sin[x], 2]]]]",
        100.0);
    fprintf(stderr, "  trig monomial: %.0f ms\n", ms);
}

/* User-reported case (2026-04): a trivariate 18-term polynomial whose
 * FactorSquareFree pre-check ran a univariate gcd over a degree-31
 * image.  The original Expr-level `poly_gcd_internal` (Knuth-style
 * primitive PRS) suffered exponential coefficient growth and took
 * >120s to terminate.  The fix routes `univariate_squarefree`'s
 * internal gcd through `zupoly_gcd` (subresultant PRS), bringing the
 * full Factor call to under 1.5 s.  Mathematica solves it in 2 ms,
 * Mathics in 100 ms; we're not at parity but we're no longer hung.
 *
 * Expected: x^2 * (z^13 - x^12) * (17 - 5y - z^14) * (z^4 + 3x^9 - y^13). */
static void test_factor_trivariate_high_degree_squarefree(void) {
    double ms = eval_check(
        "Factor[Expand[x^2 (z^13 - x^12) (z^4 + 3 x^9 - y^13) (17 - 5 y - z^14)]]",
        "Times["
            "Power[x, 2], "
            "Plus[Times[-1, Power[x, 12]], Power[z, 13]], "
            "Plus[17, Times[-5, y], Times[-1, Power[z, 14]]], "
            "Plus[Times[3, Power[x, 9]], Times[-1, Power[y, 13]], Power[z, 4]]"
        "]",
        5000.0);
    fprintf(stderr, "  trivariate user case: %.0f ms\n", ms);
}

int main(void) {
    symtab_init();
    core_init();

    printf("Running factor baseline tests...\n");

    /* Direct Factor */
    TEST(test_univariate_factors);
    TEST(test_bivariate_simple);
    TEST(test_bivariate_hensel);
    TEST(test_multivariate_via_factor_roots);
    TEST(test_irreducibility_short_circuit);
    TEST(test_factor_threads_over_logic_heads);
    TEST(test_factor_denominator_with_extra_vars);

    /* Simplify integration */
    TEST(test_simplify_user_case);
    TEST(test_simplify_tan_double_angle);
    TEST(test_simplify_pythagorean_power);
    TEST(test_simplify_hyperbolic_pythag);
    TEST(test_simplify_chebyshev_cubic);
    TEST(test_simplify_cubic_expansion);
    TEST(test_simplify_sin2x_over_sinx);

    /* Direct Factor regression cases */
    TEST(test_factor_monomial_cube);
    TEST(test_factor_trig_monomial);
    TEST(test_factor_trivariate_high_degree_squarefree);

    printf("All baseline tests passed!\n");
    return 0;
}
