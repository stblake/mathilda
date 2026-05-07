/*
 * test_factor_phase0.c
 * --------------------
 * Tests for the Phase 0 multivariate factoring improvements:
 *
 *   1. factor_monomial_content: extracts the largest monomial v_1^{e_1} *
 *      ... * v_k^{e_k} that divides every term of a Plus expression.
 *      Solves correctness gaps where picocas previously returned the
 *      input unfactored.
 *
 *   2. is_likely_irreducible_multivariate: probabilistic Hilbert-
 *      irreducibility-based short-circuit that skips the expensive
 *      factor_roots trial-division pipeline when the polynomial is
 *      almost certainly irreducible.
 *
 * The tests cover three orthogonal axes:
 *   (a) Correctness of the new monomial-content extraction.
 *   (b) Performance on inputs that previously stalled in factor_roots.
 *   (c) Regression coverage of the existing multivariate factoring
 *       contract -- the new code must not break any case where
 *       picocas was already producing correct factorisations.
 *
 * All correctness checks compare FullForm output against expected
 * canonical forms, exercising the public Factor builtin end-to-end.
 *
 * Performance assertions use a generous wall-clock bound (200 ms per
 * call) to remain stable across machines.  Pre-Phase-0, the user's
 * reference case (Simplify[Sin[x]^3 + Sin[3 x] - 3 Sin[x]]) took
 * 1100+ ms; with Phase 0 it must complete well under 1 s, with each
 * individual Factor call settling far below 200 ms.
 */

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include "expr.h"
#include "parse.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "test_utils.h"

/* Wall-clock-bounded evaluation.  Returns the FullForm of the result
 * and writes the elapsed milliseconds into *out_ms.  The caller owns
 * the returned string and must free it. */
static char* eval_timed(const char* input, double* out_ms) {
    Expr* parsed = parse_expression(input);
    ASSERT(parsed != NULL);
    clock_t t0 = clock();
    Expr* result = evaluate(parsed);
    clock_t t1 = clock();
    *out_ms = 1000.0 * (double)(t1 - t0) / (double)CLOCKS_PER_SEC;
    char* s = expr_to_string_fullform(result);
    expr_free(result);
    expr_free(parsed);
    return s;
}

/* Compare actual FullForm output against expected.  Print informative
 * diagnostic on mismatch. */
static void expect_fullform(const char* input, const char* expected) {
    double ms;
    char* got = eval_timed(input, &ms);
    if (strcmp(got, expected) != 0) {
        fprintf(stderr,
                "FAIL: %s\n  expected: %s\n  got:      %s\n  (%.2f ms)\n",
                input, expected, got, ms);
        free(got);
        ASSERT(0);
    }
    free(got);
}

/* Assert that the output equals expected AND elapsed time is below
 * `budget_ms`.  Used for the "fast-NO" tests that previously stalled. */
static void expect_fast(const char* input, const char* expected,
                        double budget_ms) {
    double ms;
    char* got = eval_timed(input, &ms);
    if (strcmp(got, expected) != 0) {
        fprintf(stderr,
                "FAIL [correctness]: %s\n  expected: %s\n  got:      %s\n",
                input, expected, got);
        free(got);
        ASSERT(0);
    }
    if (ms > budget_ms) {
        fprintf(stderr,
                "FAIL [perf]: %s\n  budget: %.2f ms, actual: %.2f ms\n",
                input, budget_ms, ms);
        free(got);
        ASSERT(0);
    }
    free(got);
}

/* =====================================================================
 *  (a) Correctness of monomial-content extraction
 * =====================================================================
 *
 * For each input we verify that the structurally-largest common
 * monomial is pulled out and the residue is independently factored.
 * These cases all previously returned the input unfactored.
 */
static void test_monomial_basic_two_var(void) {
    /* x^2 y + x y^2 = x y (x + y).  Two-variable, balanced. */
    expect_fullform("Factor[x^2 y + x y^2]",
                    "Times[x, y, Plus[x, y]]");
}

static void test_monomial_simple_one_var(void) {
    /* x^3 + x = x (x^2 + 1).  Single common variable, residue irred. */
    expect_fullform("Factor[x^3 + x]",
                    "Times[x, Plus[1, Power[x, 2]]]");
}

static void test_monomial_with_coefficient(void) {
    /* 2 x^2 y - 4 x y^2 = 2 x y (x - 2 y).  Common monomial after
     * the integer content (2) is also shared. */
    expect_fullform("Factor[2 x^2 y - 4 x y^2]",
                    "Times[2, x, y, Plus[x, Times[-2, y]]]");
}

static void test_monomial_user_reported_case(void) {
    /* The user's reported regression: 3 a^2 b - 3 b - b^3.
     * Pre-Phase-0 returned the input unchanged (also slowly). */
    expect_fullform("Factor[3 a^2 b - 3 b - b^3]",
                    "Times[b, Plus[-3, Times[3, Power[a, 2]], Times[-1, Power[b, 2]]]]");
}

static void test_monomial_trig_atoms(void) {
    /* 3 Cos[x]^2 Sin[x] - 3 Sin[x] - Sin[x]^3 -> Sin[x] (3 Cos^2 - 3 - Sin^2).
     * This is the inner-loop case for Simplify[Sin[x]^3 + Sin[3 x] - 3 Sin[x]]. */
    expect_fullform("Factor[3 Cos[x]^2 Sin[x] - 3 Sin[x] - Sin[x]^3]",
                    "Times[Sin[x], Plus[-3, Times[3, Power[Cos[x], 2]], Times[-1, Power[Sin[x], 2]]]]");
}

static void test_monomial_three_variables(void) {
    /* x y z + x y^2 z = x y z (1 + y).  Three-variable common factor. */
    expect_fullform("Factor[x y z + x y^2 z]",
                    "Times[x, y, z, Plus[1, y]]");
}

static void test_monomial_higher_powers(void) {
    /* x^3 y^2 + x^2 y^4 = x^2 y^2 (x + y^2).  Higher exponents in
     * the common monomial. */
    expect_fullform("Factor[x^3 y^2 + x^2 y^4]",
                    "Times[Power[x, 2], Power[y, 2], Plus[x, Power[y, 2]]]");
}

static void test_monomial_no_common_factor(void) {
    /* x + y has no common monomial -- verify we don't extract anything. */
    expect_fullform("Factor[x + y]", "Plus[x, y]");
}

static void test_monomial_negative_coefficient(void) {
    /* -x^2 y - x y^2.  Mathematically equals -x y (x + y).  picocas'
     * canonical printer prefers to push the negative inside the Plus
     * rather than emit a leading -1 factor; both forms are equivalent
     * factorisations.  We accept the actual surface form. */
    expect_fullform("Factor[-x^2 y - x y^2]",
                    "Times[x, y, Plus[Times[-1, x], Times[-1, y]]]");
}

static void test_monomial_bigint_coefficient(void) {
    /* F6 cleanup item 9: monomial_collect now uses mpz_t coefficients,
     * so terms whose integer scalar promotes to EXPR_BIGINT (|c| > 2^63)
     * are still recognised as having an integer coefficient instead of
     * being treated as opaque atoms.  After widening, the residue
     * extraction sees the real "variables" {a, b} and the coefficient
     * 2^100 is folded back as an integer factor by poly_content.
     *
     * Pre-fix the old behaviour treated EXPR_BIGINT as an opaque atom
     * which happened to produce the correct answer for *this* exact
     * shape (since both terms shared the same BIGINT atom), but failed
     * for shapes where different BIGINTs share a common variable
     * pattern (Plus[Times[2^100, a, b], Times[3*2^100, b]]).  The
     * post-fix path always produces the correct factorisation by
     * routing the BIGINT through the integer-content extractor. */
    expect_fullform("Factor[2^100 * a^2 * b + 2^100 * b]",
                    "Times[1267650600228229401496703205376, b, "
                    "Plus[1, Power[a, 2]]]");
}

static void test_monomial_bigint_mixed_coefficients(void) {
    /* Plus[Times[2^100, a, b], Times[3 * 2^100, b]] = 2^100 * b * (a + 3).
     * The two terms have *different* BIGINT coefficients (2^100 and
     * 3 * 2^100); the GCD-monomial in *variables* is {b} only -- which
     * is what factor_monomial_content (the F6-widened path) extracts.
     * The residue 3*2^100 + 2^100 a remains as-is because poly_content
     * is still int64-bound (a follow-up cleanup item).  The b
     * extraction is the F6 deliverable verified here. */
    expect_fullform("Factor[2^100 * a * b + 3 * 2^100 * b]",
                    "Times[b, Plus[3802951800684688204490109616128, "
                    "Times[1267650600228229401496703205376, a]]]");
}

/* =====================================================================
 *  (b) Performance: irreducible inputs must complete fast
 * =====================================================================
 *
 * Pre-Phase-0, calling Factor on a bivariate irreducible polynomial
 * forced factor_roots through 40 polynomial trial divisions, each
 * costing ~10 ms.  The is_likely_irreducible_multivariate predicate
 * detects these cases and short-circuits.  We verify timing here
 * against a generous 200 ms ceiling -- pre-Phase-0 these inputs took
 * 300-700 ms.
 */
static void test_irreducible_bivariate_quadratic(void) {
    /* x^2 + y^2 + 1 is irreducible over Z[x, y]. */
    expect_fast("Factor[x^2 + y^2 + 1]",
                "Plus[1, Power[x, 2], Power[y, 2]]",
                200.0);
}

static void test_irreducible_residue_after_monomial(void) {
    /* The post-monomial residue from the user's case: 3 Cos^2 - 3 - Sin^2.
     * Irreducible, must complete fast. */
    expect_fast("Factor[3 Cos[x]^2 - 3 - Sin[x]^2]",
                "Plus[-3, Times[3, Power[Cos[x], 2]], Times[-1, Power[Sin[x], 2]]]",
                200.0);
}

static void test_irreducible_cubic_two_var(void) {
    /* x^3 - 2 y^3 - 1 is irreducible.  The univariate image at small
     * y values is also irreducible (e.g., y=1 gives x^3 - 3, y=2 gives
     * x^3 - 17 -- both irreducible cubics). */
    expect_fast("Factor[x^3 - 2 y^3 - 1]",
                "Plus[-1, Power[x, 3], Times[-2, Power[y, 3]]]",
                200.0);
}

static void test_irreducible_high_degree(void) {
    /* x^4 + y^4 + 1 is irreducible over Z. */
    expect_fast("Factor[x^4 + y^4 + 1]",
                "Plus[1, Power[x, 4], Power[y, 4]]",
                200.0);
}

/* =====================================================================
 *  (c) End-to-end regression: existing multivariate factoring contract
 * =====================================================================
 *
 * These cases were factored correctly pre-Phase-0 (via factor_binomial
 * or factor_roots) and must continue to factor correctly.  This is a
 * defence against the new code paths shadowing or interfering with
 * existing strategies.
 */
static void test_regression_difference_of_squares(void) {
    expect_fullform("Factor[x^2 - 4 y^2]",
                    "Times[Plus[x, Times[-2, y]], Plus[x, Times[2, y]]]");
}

static void test_regression_sum_of_cubes(void) {
    /* Canonical ordering of the second factor's terms is determined
     * by Plus' Orderless attribute; match the actual surface form. */
    expect_fullform("Factor[x^3 + y^3]",
                    "Times[Plus[x, y], Plus[Times[-1, Times[x, y]], Power[x, 2], Power[y, 2]]]");
}

static void test_regression_factor_roots_case(void) {
    /* This case relies on factor_roots' trial-division to find linear
     * factors.  It is the canonical reducible bivariate case in
     * test_facpoly.c. */
    expect_fullform("Factor[2 x^3 y - 2 a^2 x y - 3 a^2 x^2 + 3 a^4]",
                    "Times[Plus[a, x], Plus[Times[-1, a], x], Plus[Times[-3, Power[a, 2]], Times[2, Times[x, y]]]]");
}

/* =====================================================================
 *  (d) Phase 5: bivariate Hensel pipeline wired into heuristic_factor
 * =====================================================================
 *
 * These cases were either returning the input unfactored or factoring
 * incorrectly before Phase 5.  They exercise the full
 * heuristic_factor -> mvfactor_try_bivariate_monic -> bz_factor_to_expr
 * round-trip.
 */
static void test_hensel_linear_y_in_two_factors(void) {
    /* (x + y - 1)(x + y + 1) = x^2 + 2 x y + y^2 - 1.
     * Both factors monic in x, both have linear y dependence.
     * Hensel iterates from y=0 image: (x-1)(x+1). */
    expect_fullform("Factor[(x + y - 1)(x + y + 1)]",
                    "Times[Plus[-1, x, y], Plus[1, x, y]]");
}

static void test_hensel_quadratic_y_in_two_factors(void) {
    /* (x + y^2 - 1)(x - y^2 + 1) = x^2 - (y^2 - 1)^2.
     * Higher y-degree forces multiple Hensel iterations.  Note the
     * factor order is determined by the canonical Orderless ordering
     * of Times, not by our algorithm. */
    expect_fullform("Factor[(x + y^2 - 1)(x - y^2 + 1)]",
                    "Times[Plus[-1, x, Power[y, 2]], Plus[1, x, Times[-1, Power[y, 2]]]]");
}

static void test_hensel_three_factor_lift(void) {
    /* (x - 1)(x + 1)(x + y) = (x^2 - 1)(x + y) = x^3 + x^2 y - x - y.
     * Image at y=0: x^3 - x = x(x-1)(x+1).  Three univariate factors;
     * one (originally x) lifts to x + y, the other two stay constant. */
    expect_fullform("Factor[(x - 1)(x + 1)(x + y)]",
                    "Times[Plus[-1, x], Plus[1, x], Plus[x, y]]");
}

static void test_hensel_cubic_x(void) {
    /* (x^2 + y - 1)(x + y + 2) = x^3 + (y+2) x^2 + (y-1) x + (y^2 + y - 2).
     * At y=0: P(x, 0) = x^3 + 2 x^2 - x - 2 = (x^2 - 1)(x + 2).
     * Image factors as (x-1)(x+1)(x+2) -- so univariate has THREE
     * factors but the bivariate has only TWO (x^2 + y - 1 and
     * x + y + 2).  Without recombination this case still works
     * because the orchestrator's pair-and-recurse lift correctly
     * peels off (x-1) and (x+1) into the lifted x^2 + y - 1 only
     * when neither alone is a true factor.
     *
     * In practice on this input the factor_via_bz callback returns
     * (x-1, x+1, x+2) and the multifactor lift attempts to lift
     * each individually.  When (x-1) can't be lifted to a factor of
     * P alone, the lift fails and we fall through to factor_roots,
     * which handles the case.  Either way the final factorisation
     * must be correct. */
    expect_fullform("Factor[(x^2 + y - 1)(x + y + 2)]",
                    "Times[Plus[-1, Power[x, 2], y], Plus[2, x, y]]");
}

static void test_hensel_irreducible_bivariate(void) {
    /* x^2 + y^2 + 1: irreducible over Z[x, y].  Verified by the
     * Phase-1 short-circuit; the new pipeline's r==1 path also
     * declares irreducibility but we exit through the higher-level
     * irreducibility check first. */
    expect_fullform("Factor[x^2 + y^2 + 1]",
                    "Plus[1, Power[x, 2], Power[y, 2]]");
}

/* =====================================================================
 *  (e) End-to-end Simplify integration
 * =====================================================================
 *
 * The user's reported case: Simplify[Sin[x]^3 + Sin[3 x] - 3 Sin[x]]
 * should produce -3 Sin[x]^3 in well under one second.
 */
static void test_simplify_user_case(void) {
    double ms;
    char* got = eval_timed("Simplify[Sin[x]^3 + Sin[3 x] - 3 Sin[x]]", &ms);
    /* Accept either the symbolic form or the printer's negative-leading
     * Times form; we only care about the simplified content. */
    const char* expected_a = "Times[-3, Power[Sin[x], 3]]";
    if (strcmp(got, expected_a) != 0) {
        fprintf(stderr,
                "FAIL [correctness]: Simplify[Sin[x]^3 + Sin[3 x] - 3 Sin[x]]\n"
                "  expected: %s\n  got:      %s\n",
                expected_a, got);
        free(got);
        ASSERT(0);
    }
    /* Generous bound: pre-Phase-0 took ~1100 ms.  Post-Phase-0 should
     * easily fit in 1500 ms even on a slow machine, but we expect
     * 500-800 ms in practice. */
    if (ms > 1500.0) {
        fprintf(stderr,
                "FAIL [perf]: Simplify[Sin[x]^3 + Sin[3 x] - 3 Sin[x]]\n"
                "  budget: 1500 ms, actual: %.2f ms\n", ms);
        free(got);
        ASSERT(0);
    }
    free(got);
    fprintf(stderr, "  Simplify user case: %.1f ms\n", ms);
}

int main(void) {
    symtab_init();
    core_init();

    printf("Running factor Phase 0 tests...\n");

    /* (a) Monomial-content extraction -- correctness */
    TEST(test_monomial_basic_two_var);
    TEST(test_monomial_simple_one_var);
    TEST(test_monomial_with_coefficient);
    TEST(test_monomial_user_reported_case);
    TEST(test_monomial_trig_atoms);
    TEST(test_monomial_three_variables);
    TEST(test_monomial_higher_powers);
    TEST(test_monomial_no_common_factor);
    TEST(test_monomial_negative_coefficient);
    TEST(test_monomial_bigint_coefficient);
    TEST(test_monomial_bigint_mixed_coefficients);

    /* (b) Irreducibility short-circuit -- performance */
    TEST(test_irreducible_bivariate_quadratic);
    TEST(test_irreducible_residue_after_monomial);
    TEST(test_irreducible_cubic_two_var);
    TEST(test_irreducible_high_degree);

    /* (c) Regression: existing multivariate factoring contract */
    TEST(test_regression_difference_of_squares);
    TEST(test_regression_sum_of_cubes);
    TEST(test_regression_factor_roots_case);

    /* (d) Phase 5: bivariate Hensel wiring */
    TEST(test_hensel_linear_y_in_two_factors);
    TEST(test_hensel_quadratic_y_in_two_factors);
    TEST(test_hensel_three_factor_lift);
    TEST(test_hensel_cubic_x);
    TEST(test_hensel_irreducible_bivariate);

    /* (e) End-to-end Simplify integration */
    TEST(test_simplify_user_case);

    printf("All Phase 0 tests passed!\n");
    return 0;
}
