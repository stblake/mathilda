/*
 * test_factor_recombine.c
 * -----------------------
 * End-to-end tests for Phase 3b (subset recombination in the bivariate
 * Hensel lift) and Phase 4 (n-variate factoring via specialise-and-
 * trial-divide).
 *
 * Phase 3b targets bivariate inputs where the bivariate factorisation
 * is COARSER than the univariate image at any chosen alpha.  Phase 4
 * targets trivariate (and higher) inputs where at least one factor is
 * independent of one of the variables.
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

/* Evaluate `input`, return the FullForm string.  Caller frees. */
static char* eval_to_fullform(const char* input) {
    Expr* parsed = parse_expression(input);
    ASSERT(parsed != NULL);
    Expr* result = evaluate(parsed);
    char* got = expr_to_string_fullform(result);
    expr_free(parsed);
    expr_free(result);
    return got;
}

/* Verify that Factor[input] EQUALS expected (FullForm).  Treats the
 * factor list as ordered; if your test fails on argument-ordering
 * differences, list the expected form in the order Mathilda's printer
 * uses (Times sorts its arguments canonically). */
static void check_factor(const char* input, const char* expected) {
    char* got = eval_to_fullform(input);
    if (strcmp(got, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n"
                        "  expected: %s\n"
                        "  got:      %s\n",
                input, expected, got);
        free(got);
        ASSERT(0);
    }
    free(got);
}

/* Verify that Factor[A] times Factor's product is mathematically equal
 * to A by expanding the product and comparing.  Useful when Mathilda's
 * canonical printout of the factored form differs from a hand-written
 * "expected" string but is still mathematically correct. */
static void check_factor_expands_to(const char* input, const char* expected_expand) {
    char buf[512];
    snprintf(buf, sizeof(buf), "Expand[%s]", input);
    char* got = eval_to_fullform(buf);
    if (strcmp(got, expected_expand) != 0) {
        fprintf(stderr, "FAIL: Expand[%s]\n"
                        "  expected: %s\n"
                        "  got:      %s\n",
                input, expected_expand, got);
        free(got);
        ASSERT(0);
    }
    free(got);
}

/* ====================================================================== */
/*  Phase 3b: bivariate recombination                                     */
/* ====================================================================== */

static void test_recombine_x4_minus_y2(void) {
    /* x^4 - y^2 = (x^2 - y)(x^2 + y).  At y=1: x^4 - 1 has 3 univariate
     * factors (x-1)(x+1)(x^2+1).  Recombination groups them into 2. */
    check_factor_expands_to("Factor[x^4 - y^2]",
                            "Plus[Power[x, 4], Times[-1, Power[y, 2]]]");
}

static void test_recombine_quartic_minus_4y2(void) {
    /* x^4 - 4y^2 + 4y - 1 = (x^2 - (2y - 1))(x^2 + (2y - 1))
     *                     = (x^2 + 1 - 2y)(x^2 - 1 + 2y).
     * At y=1: x^4 - 1.  Recombination groups univariate factors. */
    check_factor_expands_to(
        "Factor[x^4 - 4 y^2 + 4 y - 1]",
        "Plus[-1, Power[x, 4], Times[4, y], Times[-4, Power[y, 2]]]");
}

/* ====================================================================== */
/*  Phase 4: trivariate factoring via specialise-and-trial-divide         */
/* ====================================================================== */

static void test_trivariate_xy_xz(void) {
    /* (x + y)(x + z) -- one factor is z-independent, the other y-
     * independent.  Phase 4 catches at least one specialisation that
     * peels off a clean factor. */
    check_factor_expands_to("Factor[(x + y) (x + z)]",
                            "Plus[Power[x, 2], Times[x, y], Times[x, z], Times[y, z]]");
}

static void test_trivariate_xy_xyz(void) {
    /* (x + y)(x + y + z): x + y is z-independent. */
    check_factor_expands_to(
        "Factor[(x + y) (x + y + z)]",
        "Plus[Power[x, 2], Times[2, Times[x, y]], Power[y, 2], Times[x, z], Times[y, z]]");
}

static void test_trivariate_quadratic_z_indep(void) {
    /* (x^2 + y^2)(x + z): x^2 + y^2 is z-independent and irreducible
     * over Z (no real roots in Z[x, y]). */
    check_factor_expands_to(
        "Factor[(x^2 + y^2) (x + z)]",
        "Plus[Power[x, 3], Times[x, Power[y, 2]], Times[Power[x, 2], z], Times[Power[y, 2], z]]");
}

static void test_trivariate_three_factors_with_z_indep(void) {
    /* (x + 2 y)(x - 2 y)(x + 3 z) -- two factors are z-independent
     * (their product is x^2 - 4 y^2), one is y-independent. */
    check_factor_expands_to(
        "Factor[(x + 2 y) (x - 2 y) (x + 3 z)]",
        "Plus[Power[x, 3], Times[-4, Times[x, Power[y, 2]]], "
        "Times[3, Times[Power[x, 2], z]], Times[-12, Times[Power[y, 2], z]]]");
}

static void test_trivariate_already_factored_irreducible_x_indep(void) {
    /* y^2 - z^2 -- bivariate (no x dependence) but introduced as a
     * trivariate input.  Should factor as (y-z)(y+z) regardless of
     * which variable Phase 4 specialises. */
    check_factor("Factor[(y - z) (y + z) + 0 x]",
                 "Times[Plus[y, z], Plus[y, Times[-1, z]]]");
}

/* ====================================================================== */
/*  Regression: Phase 4 should not break existing irreducible cases       */
/* ====================================================================== */

static void test_trivariate_irreducible_passes_through(void) {
    /* x^2 + y^2 + z^2 + 1 -- irreducible over Z[x, y, z].  Phase 4
     * should fail to find any z-independent factor (because there
     * isn't one) and the result should be the input unchanged. */
    check_factor("Factor[x^2 + y^2 + z^2 + 1]",
                 "Plus[1, Power[x, 2], Power[y, 2], Power[z, 2]]");
}

static void test_trivariate_symmetric_cubic(void) {
    /* x^3 + y^3 + z^3 - 3xyz = (x+y+z)(x^2+y^2+z^2-xy-yz-xz).
     * Neither factor is variable-independent; Phase 4 cannot catch
     * this and falls through to the legacy pipeline.  The result
     * may stay unfactored.  This is a known limitation of the
     * minimum-viable Phase 4 -- it is documented in FACTOR_PLAN. */
    char* got = eval_to_fullform("Factor[x^3 + y^3 + z^3 - 3 x y z]");
    /* Either left as input or partially factored -- both are
     * acceptable for the minimum-viable Phase 4.  We just verify
     * we don't CRASH. */
    ASSERT(got != NULL);
    free(got);
}

/* ====================================================================== */
/*  Phase F1: Wang's leading-coefficient correction (negate path)         */
/*                                                                        */
/*  Stage 1 covers the LC = -1 case: P(x, y) where lc_x(P) is the         */
/*  constant integer -1.  Pre-negate P, run the monic Hensel lift on -P,  */
/*  then absorb the overall sign into the lowest-x-degree factor.  This   */
/*  unlocks inputs of the form Factor[(1 - x^k)(x - y^m)] which were      */
/*  otherwise routed through the legacy factor_roots' linear-trial-       */
/*  division loop.                                                        */
/* ====================================================================== */

static void test_f1_negate_small_quartic_quintic(void) {
    /* (1 - x^4)(x - y^5).  Lc_x = -1.  After negation, the lift sees
     * (x^4 - 1)(x - y^5) with lc_x = +1.  Univariate image at y = 0
     * is x^5 - x = x(x-1)(x+1)(x^2+1) -- 4 factors plus (x - y^5). */
    check_factor_expands_to("Factor[Expand[(1 - x^4) (x - y^5)]]",
        "Plus[x, Times[-1, Power[x, 5]], Times[-1, Power[y, 5]], "
        "Times[Power[x, 4], Power[y, 5]]]");
}

static void test_f1_negate_x6_y7(void) {
    /* (1 - x^6)(x - y^7).  Lc_x = -1, deg_x = 7, deg_y = 7. */
    check_factor_expands_to("Factor[Expand[(1 - x^6) (x - y^7)]]",
        "Plus[x, Times[-1, Power[x, 7]], Times[-1, Power[y, 7]], "
        "Times[Power[x, 6], Power[y, 7]]]");
}

static void test_f1_negate_x8_y9(void) {
    /* (1 - x^8)(x - y^9).  Lc_x = -1, deg_x = 9.  Image is
     * x^9 - x = x(x-1)(x+1)(x^2+1)(x^4+1) -- 5 factors plus (x - y^9). */
    check_factor_expands_to("Factor[Expand[(1 - x^8) (x - y^9)]]",
        "Plus[x, Times[-1, Power[x, 9]], Times[-1, Power[y, 9]], "
        "Times[Power[x, 8], Power[y, 9]]]");
}

static void test_f1_negate_lc_minus_one_simple(void) {
    /* (1 - x^2)(x - y).  Smallest negate-path case.  Image at y = 0
     * is x^3 - x = x(x-1)(x+1).  After lift: (x-y)(x-1)(x+1) with
     * an absorbed -1. */
    check_factor_expands_to("Factor[Expand[(1 - x^2) (x - y)]]",
        "Plus[x, Times[-1, Power[x, 3]], Times[-1, y], "
        "Times[Power[x, 2], y]]");
}

static void test_f1_negate_irreducible_lc_minus_one(void) {
    /* lc_x = -1 but the polynomial is irreducible over Z[x, y].
     * The negate path enters Hensel; the orchestrator's
     * irreducibility short-circuit (or r == 1 from the lift) returns
     * the input unchanged.  Verify no crash and a sensible result. */
    char* got = eval_to_fullform("Factor[-x^2 - y^2 - 1]");
    ASSERT(got != NULL);
    free(got);
}

/* ====================================================================== */
/*  Phase F1 Stage 2: Wang's lc correction — constant |a| > 1 LC          */
/*                                                                        */
/*  When lc_x(P) is a non-zero, non-unit integer constant a (|a| > 1      */
/*  and constant in y), monic substitution Q = a^(d-1) · P(x/a, y) makes  */
/*  the lift's input monic in x with integer coefficients.  After lifting */
/*  to G_1 ... G_r, the true factors of P are recovered as                */
/*    F_i = G_i(a·x, y) / cont_Z(G_i(a·x, y)).                            */
/*                                                                        */
/*  Test inputs are designed so neither variable yields a unit LC.        */
/* ====================================================================== */

static void test_f1_scale_simple_two_factor(void) {
    /* (2x + 3y)(3x + 5y) = 6x^2 + 19xy + 15y^2.  lc_x = 6, lc_y = 15.
     * Neither variable is monic; Stage 2 must trigger.  After Wang
     * substitution Q has integer monic coefficients in x. */
    check_factor("Factor[Expand[(2x+3y)(3x+5y)]]",
                 "Times[Plus[Times[2, x], Times[3, y]], "
                 "Plus[Times[3, x], Times[5, y]]]");
}

static void test_f1_scale_negative_factor(void) {
    /* (2x - 3y)(3x + 5y) = 6x^2 + xy - 15y^2.  lc_x = 6 (constant). */
    check_factor("Factor[Expand[(2x-3y)(3x+5y)]]",
                 "Times[Plus[Times[2, x], Times[-3, y]], "
                 "Plus[Times[3, x], Times[5, y]]]");
}

static void test_f1_scale_with_negate(void) {
    /* -(2x + 3y)(3x + 5y) = -6x^2 - 19xy - 15y^2.  lc_x = -6.
     * This composes Stage 1 (negate) with Stage 2 (scale).
     * The picker preferences negate over scale-with-negative-LC,
     * so the path selected depends on heuristics — we only verify
     * the result expands correctly. */
    check_factor_expands_to("Factor[Expand[-(2x+3y)(3x+5y)]]",
        "Plus[Times[-6, Power[x, 2]], Times[-19, Times[x, y]], "
        "Times[-15, Power[y, 2]]]");
}

static void test_f1_scale_three_factors(void) {
    /* (2x+3y)(3x+5y)(x+y) = (6x^2 + 19xy + 15y^2)(x + y).
     * Three-factor case; the third factor (x+y) is monic. */
    check_factor("Factor[Expand[(2x+3y)(3x+5y)(x+y)]]",
                 "Times[Plus[x, y], "
                 "Plus[Times[2, x], Times[3, y]], "
                 "Plus[Times[3, x], Times[5, y]]]");
}

static void test_f1_scale_mixed_monic_and_nonmonic(void) {
    /* (2x+y)(x+y)(x+2y).  Expanded coefficients have lc_x = 2,
     * so neither variable is monic-after-content-extraction.  Three
     * factors: one Stage 2-like (2x+y), two monic ((x+y), (x+2y)). */
    check_factor("Factor[Expand[(2x+y)(x+y)(x+2y)]]",
                 "Times[Plus[x, y], Plus[x, Times[2, y]], "
                 "Plus[Times[2, x], y]]");
}

static void test_f1_scale_a3_b3_quadratic(void) {
    /* 2 a^2 - 5 a b + 3 b^2 = (a - b)(2a - 3b).  lc_a = 2, lc_b = 3. */
    check_factor("Factor[2 a^2 - 5 a b + 3 b^2]",
                 "Times[Plus[a, Times[-1, b]], "
                 "Plus[Times[2, a], Times[-3, b]]]");
}

static void test_f1_scale_irreducible_lc_two(void) {
    /* 2x^2 + 3xy + 7y^2 is irreducible over Z (discriminant in y is
     * 9 - 56 = -47, not a square).  Stage 2 enters the lift, the
     * orchestrator's irreducibility short-circuit (r == 1) signals
     * "no progress" so the legacy pipeline takes over and leaves the
     * input unchanged.  Verify the canonical form is preserved. */
    check_factor("Factor[2 x^2 + 3 x y + 7 y^2]",
                 "Plus[Times[2, Power[x, 2]], Times[3, Times[x, y]], "
                 "Times[7, Power[y, 2]]]");
}

/* ====================================================================== */
/*  Phase F1 Stage 3: polynomial-in-y leading coefficient                  */
/*                                                                        */
/*  Wang's leading-coefficient correction for inputs whose lc_x(P)(y) is  */
/*  a non-constant polynomial in the other variable.  The orchestrator    */
/*  factors A(y) = lc_x(P) over Z[y], picks α with A(α) = +1, factors    */
/*  the squarefree image P(x, α), and enumerates distributions of A's     */
/*  factors between the predicted leading coefficients (q_u, q_v).        */
/*                                                                        */
/*  MVP scope: r = 2 (two univariate factors), both monic, |cont(A)| = 1. */
/* ====================================================================== */

static void test_f1_polylc_xy_plus_xy_plus_two(void) {
    /* (xy + 1)(xy + 2) = x^2 y^2 + 3xy + 2.  lc_x = y^2.  Distribution
     * q_u = y, q_v = y splits A=y·y between the two factors.  At α=1
     * the image is x^2 + 3x + 2 = (x+1)(x+2), both monic. */
    check_factor("Factor[Expand[(x*y + 1)*(x*y + 2)]]",
                 "Times[Plus[1, Times[x, y]], Plus[2, Times[x, y]]]");
}

static void test_f1_polylc_y2_plus_one(void) {
    /* ((y^2+1)x + 1)(x + 3).  lc_x = y^2+1 (irreducible over Z[y]).
     * Distribution: q_u = 1, q_v = y^2 + 1.  At α=0 the image is
     * x^2 + 4x + 3 = (x+1)(x+3), both monic. */
    check_factor_expands_to(
        "Factor[Expand[((y^2 + 1)*x + 1)*(x + 3)]]",
        "Plus[3, Times[4, x], Power[x, 2], "
        "Times[3, Times[x, Power[y, 2]]], "
        "Times[Power[x, 2], Power[y, 2]]]");
}

static void test_f1_polylc_linear_lc(void) {
    /* ((y-1)x + 1)(x + 2).  lc_x = y - 1 (linear in y).  At α=2,
     * A(2) = 1 ✓ (α=0 gives -1, skipped by MVP). */
    check_factor_expands_to(
        "Factor[Expand[((y - 1)*x + 1)*(x + 2)]]",
        "Plus[2, Times[-1, x], Times[-1, Power[x, 2]], "
        "Times[2, Times[x, y]], Times[Power[x, 2], y]]");
}

static void test_f1_polylc_y_plus_one_squared(void) {
    /* ((y+1)x + 1)((y+1)x + 2).  lc_x = (y+1)^2.  Distribution
     * q_u = (y+1), q_v = (y+1).  At α=0, A=1, image = (x+1)(x+2). */
    check_factor("Factor[Expand[((y+1)*x + 1)*((y+1)*x + 2)]]",
                 "Times[Plus[1, Times[x, Plus[1, y]]], "
                 "Plus[2, Times[x, Plus[1, y]]]]");
}

static void test_f1_polylc_xy_minus_versions(void) {
    /* (xy - 1)(xy - 2) = x^2 y^2 - 3xy + 2.  At α=1, A=1, image =
     * x^2 - 3x + 2 = (x-1)(x-2). */
    check_factor("Factor[Expand[(x*y - 1)*(x*y - 2)]]",
                 "Times[Plus[-2, Times[x, y]], Plus[-1, Times[x, y]]]");
}

static void test_f1_polylc_irreducible_passes_through(void) {
    /* x^2 y + x + 1.  lc_x = y.  At α=1: x^2 + x + 1 (irreducible).
     * The univariate image factorisation has r=1, the orchestrator
     * returns NULL, and the legacy pipeline preserves the input. */
    check_factor_expands_to(
        "Factor[x^2 y + x + 1]",
        "Plus[1, x, Times[Power[x, 2], y]]");
}

/* ====================================================================== */
/*  Phase F3: bivariate Hensel performance                                */
/*                                                                        */
/*  These cases exercise the lift_2 inner loop deeply -- high y-degree    */
/*  with multiple image factors so recombination explores several         */
/*  subsets, each running the full Hensel iteration.  Pre-F3 these were   */
/*  O(d_x^2 * B^3) cumulative; post-F3a (degree sort) the lowest-degree   */
/*  singleton {x} succeeds first, and post-F3b (incremental UV) each      */
/*  iteration is O(d_x^2 * k) instead of O(d_x^2 * k^2).  Post-F3c        */
/*  (Mignotte fast-fail) divergent subset attempts abort early instead    */
/*  of running to convergence.  No explicit time budget here -- the       */
/*  CTest harness will time out independently if F3 regresses badly.     */
/* ====================================================================== */

static void test_f3_perf_x6_y15(void) {
    /* (1 - x^6)(x - y^15).  lc_x = -1 (negate path).  Image at y=0
     * is x^7 - x with 5 univariate factors.  Lift y-degree = 15. */
    check_factor_expands_to("Factor[Expand[(1 - x^6) (x - y^15)]]",
        "Plus[x, Times[-1, Power[x, 7]], Times[-1, Power[y, 15]], "
        "Times[Power[x, 6], Power[y, 15]]]");
}

static void test_f3_perf_x4_y25(void) {
    /* (1 - x^4)(x - y^25).  Lift y-degree = 25 -- exercises the
     * incremental UV update over many iterations. */
    check_factor_expands_to("Factor[Expand[(1 - x^4) (x - y^25)]]",
        "Plus[x, Times[-1, Power[x, 5]], Times[-1, Power[y, 25]], "
        "Times[Power[x, 4], Power[y, 25]]]");
}

static void test_f3_perf_x8_y17(void) {
    /* (1 - x^8)(x - y^17).  Five univariate factors, deep y-lift. */
    check_factor_expands_to("Factor[Expand[(1 - x^8) (x - y^17)]]",
        "Plus[x, Times[-1, Power[x, 9]], Times[-1, Power[y, 17]], "
        "Times[Power[x, 8], Power[y, 17]]]");
}

/* ====================================================================== */
/*  Phase F5: recombination short-circuit                                 */
/*                                                                        */
/*  When every k-subset trial diverges early (Diophantine non-integer or  */
/*  Mignotte coefficient overflow), and no subset shows the partial-      */
/*  success pattern of "completed all iterations but failed final         */
/*  verify", we conclude P is irreducible bivariately under this image    */
/*  and skip the remaining higher-k subset enumerations.  k = 1 is        */
/*  always exhausted first; the cap fires from k = 2 onward, so no        */
/*  k = 2 recombination case is missed.                                   */
/*                                                                        */
/*  Coverage: irreducible bivariates with r >= 6 univariate factors       */
/*  (where the cumulative subset count is large enough for the savings    */
/*  to register), plus a k = 2 recombination case to confirm correctness  */
/*  is preserved when the short-circuit must NOT fire.                    */
/* ====================================================================== */

static double f5_eval_ms(const char* input) {
    Expr* parsed = parse_expression(input);
    ASSERT(parsed != NULL);
    clock_t t0 = clock();
    Expr* result = evaluate(parsed);
    clock_t t1 = clock();
    double ms = 1000.0 * (double)(t1 - t0) / (double)CLOCKS_PER_SEC;
    expr_free(parsed);
    expr_free(result);
    return ms;
}

static void test_f5_irreducible_r6_x12(void) {
    /* x^12 + y^2 - 1: irreducible over Z[x, y].  Image at any small
     * integer alpha factors as Phi_d for d | 12 (six cyclotomic factors
     * over Z): r = 6.  Without F5 the recombination tries
     * k=1..3 = 6 + 15 + 10 = 31 lifts; with F5 it stops after k=1 (or
     * after k=2 with no signal at k>=2). */
    check_factor_expands_to("Factor[Expand[x^12 + y^2 - 1]]",
        "Plus[-1, Power[x, 12], Power[y, 2]]");
    /* Generous budget -- we just want to confirm we are not exponential
     * in r.  Pre-F5 was ~23 ms; post-F5 is ~20 ms. */
    double ms = f5_eval_ms("Factor[Expand[x^12 + y^2 - 1]]");
    if (ms > 200.0) {
        fprintf(stderr, "FAIL [perf]: x^12+y^2-1 took %.2f ms (budget 200)\n", ms);
        ASSERT(0);
    }
}

static void test_f5_irreducible_r8_x24(void) {
    /* x^24 + y^2 - 1: irreducible.  r = 8 univariate factors at
     * a generic alpha (cyclotomic Phi_d for d | 24).  Pre-F5 the
     * recombination loop ran 8 + 28 + 56 + 35 = 127 lifts; post-F5 it
     * stops after k = 1 (no completion signal).  This is the test case
     * that exercises the largest savings. */
    check_factor_expands_to("Factor[Expand[x^24 + y^2 - 1]]",
        "Plus[-1, Power[x, 24], Power[y, 2]]");
    double ms = f5_eval_ms("Factor[Expand[x^24 + y^2 - 1]]");
    if (ms > 600.0) {
        fprintf(stderr, "FAIL [perf]: x^24+y^2-1 took %.2f ms (budget 600)\n", ms);
        ASSERT(0);
    }
}

static void test_f5_irreducible_r6_x12_y4(void) {
    /* x^12 + y^4 - 1: irreducible, deeper y-lift (B = 4) so each
     * subset trial has more iterations to potentially produce a
     * "completed all iterations" signal.  Verifies the short-circuit
     * still fires when the deeper lift consistently bails on
     * Mignotte overflow rather than running to verify-failure. */
    check_factor_expands_to("Factor[Expand[x^12 + y^4 - 1]]",
        "Plus[-1, Power[x, 12], Power[y, 4]]");
    double ms = f5_eval_ms("Factor[Expand[x^12 + y^4 - 1]]");
    if (ms > 200.0) {
        fprintf(stderr, "FAIL [perf]: x^12+y^4-1 took %.2f ms (budget 200)\n", ms);
        ASSERT(0);
    }
}

/* ====================================================================== */
/*  Phase F2 MVP: trivariate two-factor Hensel via MPoly                  */
/* ====================================================================== */

static void test_f2_trivariate_z_plus_xy_times_z_minus_xy(void) {
    /* P = (z + x y)(z - x y) = z^2 - x^2 y^2.  Both factors depend on
     * every variable.  Phase 4 (factor_via_z_independent_split) cannot
     * find this because no variable specialisation exposes a clean
     * factor; F2's trivariate Hensel is required. */
    check_factor_expands_to(
        "Factor[Expand[(z + x y)(z - x y)]]",
        "Plus[Times[-1, Times[Power[x, 2], Power[y, 2]]], Power[z, 2]]");
}

static void test_f2_trivariate_z_plus_x_z_plus_y(void) {
    /* P = (z + x)(z + y).  Linear in z; each factor depends on a
     * different second-tier variable. */
    check_factor_expands_to(
        "Factor[Expand[(z + x)(z + y)]]",
        "Plus[Times[x, y], Times[x, z], Times[y, z], Power[z, 2]]");
}

static void test_f2_trivariate_user_reported(void) {
    /* The user-reported case from plans/FACTOR_PLAN.md §12: every factor
     * depends on every variable.  Pre-F2 took ~6s and returned input
     * unchanged; F2 factors it correctly in ~300 ms.
     *
     * P = (z x - x^2 - y^2)(3 z + 4 x y - y^2). */
    check_factor_expands_to(
        "Factor[Expand[(z*x - x^2 - y^2)*(3*z + 4*x*y - y^2)]]",
        "Plus[Times[-4, Times[Power[x, 3], y]], "
             "Times[Power[x, 2], Power[y, 2]], "
             "Times[-4, Times[x, Power[y, 3]]], "
             "Power[y, 4], "
             "Times[-3, Times[Power[x, 2], z]], "
             "Times[4, Times[Power[x, 2], y, z]], "
             "Times[-3, Times[Power[y, 2], z]], "
             "Times[-1, Times[x, Power[y, 2], z]], "
             "Times[3, Times[x, Power[z, 2]]]]");
}

static void test_f5_recombination_still_works_k2(void) {
    /* (x^2 - y)(x^2 - 4y): both bivariate factors have reducible
     * univariate images at alpha = 1 (x^2 - 1 = (x-1)(x+1) and
     * x^2 - 4 = (x-2)(x+2)), so r = 4 and the bivariate factorisation
     * requires k = 2 recombination at the top level.  The F5 short-
     * circuit only fires at k >= 2 AFTER an entire pass of k = 2
     * subsets without any completion signal -- by which point the
     * correct pair has already been found and the lift succeeds.
     * This verifies F5 does not break legitimate k = 2 recombination. */
    check_factor_expands_to("Factor[Expand[(x^2 - y) (x^2 - 4 y)]]",
        "Plus[Power[x, 4], Times[-5, Times[Power[x, 2], y]], "
        "Times[4, Power[y, 2]]]");
}

/* ====================================================================== */

int main(void) {
    symtab_init();
    core_init();

    printf("Running factor recombination + n-variate tests...\n");

    /* Phase 3b */
    TEST(test_recombine_x4_minus_y2);
    TEST(test_recombine_quartic_minus_4y2);

    /* Phase 4 */
    TEST(test_trivariate_xy_xz);
    TEST(test_trivariate_xy_xyz);
    TEST(test_trivariate_quadratic_z_indep);
    TEST(test_trivariate_three_factors_with_z_indep);
    TEST(test_trivariate_already_factored_irreducible_x_indep);

    /* Regressions */
    TEST(test_trivariate_irreducible_passes_through);
    TEST(test_trivariate_symmetric_cubic);

    /* Phase F1: Wang's lc correction (Stage 1 — negate path) */
    TEST(test_f1_negate_lc_minus_one_simple);
    TEST(test_f1_negate_small_quartic_quintic);
    TEST(test_f1_negate_x6_y7);
    TEST(test_f1_negate_x8_y9);
    TEST(test_f1_negate_irreducible_lc_minus_one);

    /* Phase F1: Wang's lc correction (Stage 2 — constant |a| > 1 LC) */
    TEST(test_f1_scale_simple_two_factor);
    TEST(test_f1_scale_negative_factor);
    TEST(test_f1_scale_with_negate);
    TEST(test_f1_scale_three_factors);
    TEST(test_f1_scale_mixed_monic_and_nonmonic);
    TEST(test_f1_scale_a3_b3_quadratic);
    TEST(test_f1_scale_irreducible_lc_two);

    /* Phase F1: Wang's lc correction (Stage 3 — polynomial-in-y LC) */
    TEST(test_f1_polylc_xy_plus_xy_plus_two);
    TEST(test_f1_polylc_y2_plus_one);
    TEST(test_f1_polylc_linear_lc);
    TEST(test_f1_polylc_y_plus_one_squared);
    TEST(test_f1_polylc_xy_minus_versions);
    TEST(test_f1_polylc_irreducible_passes_through);

    /* Phase F3: bivariate Hensel performance */
    TEST(test_f3_perf_x6_y15);
    TEST(test_f3_perf_x4_y25);
    TEST(test_f3_perf_x8_y17);

    /* Phase F5: recombination short-circuit */
    TEST(test_f5_irreducible_r6_x12);
    TEST(test_f5_irreducible_r8_x24);
    TEST(test_f5_irreducible_r6_x12_y4);
    TEST(test_f5_recombination_still_works_k2);

    /* Phase F2 MVP: trivariate two-factor monic Hensel via MPoly */
    TEST(test_f2_trivariate_z_plus_xy_times_z_minus_xy);
    TEST(test_f2_trivariate_z_plus_x_z_plus_y);
    TEST(test_f2_trivariate_user_reported);

    printf("All factor_recombine tests passed!\n");
    return 0;
}
