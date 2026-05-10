/*
 * test_collect_corpus.c -- Regression corpus for Collect.
 *
 * Captures the user-supplied corpus from the 2026-05-10 audit
 * (14 cases the user marked as "passing" against MMA + Bug 1 fix +
 * Bug 2 current behaviour). Each case pins picocas's canonical
 * output string, so any future change in Collect / Power / Times
 * canonicalisation that perturbs these forms surfaces here as a
 * test failure.
 *
 * Bug 1 -- Collect[a x^(2c) + b x^(2c), x^c]
 *
 *   Pre-fix: returned the input unchanged because decompose_to_bp +
 *   the single-base path could not match Power[x, 2c] (an atomic
 *   bp entry) against the target Power[x, c]. The fix adds a
 *   power-of-power fallback in collect_internal: when the target
 *   is Power[B, e_t] and a term factor decomposes to an atomic
 *   Power[B, e_term], compute k = e_term / e_t (via Cancel) and use
 *   it as the multiplier when k is a positive integer.
 *
 * Bug 2 -- Collect[D[f[Sqrt[x^2+1]], {x,3}], Derivative[_][f][_]]
 *
 *   The target is a pattern. Picocas's Collect does NOT yet do
 *   pattern-based grouping, so each Derivative[k][f][...] factor
 *   is treated as a literal opaque atom. Test pins the current
 *   un-grouped 5-term output so an actual pattern-based Collect
 *   implementation (future work) lands deliberately, not silently.
 */

#include "test_utils.h"
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* --- Bug 1: Symbolic-exponent ratio in target ----------------------- */
static void test_bug1_symbolic_exp_ratio(void) {
    assert_eval_eq("Collect[a*x^(2*c) + b*x^(2*c), x^c]",
                   "x^(2 c) (a + b)", 0);
    /* Generalisation -- ratio 3 should also work. */
    assert_eval_eq("Collect[a*x^(3*c) + b*x^(3*c), x^c]",
                   "x^(3 c) (a + b)", 0);
    /* Mixed ratios (k=2 and k=1) should produce two groups. */
    assert_eval_eq("Collect[a*x^c + b*x^(2*c), x^c]",
                   "a x^c + b x^(2 c)", 0);
    /* Non-integer ratio (3/2) -> stays unevaluated against target x^c. */
    assert_eval_eq("Collect[a*x^(3 c/2), x^c]",
                   "a x^(3 c/2)", 0);
}

/* --- Bug 2: Pattern-target Collect (regression of current state) ---- */
static void test_bug2_pattern_target(void) {
    /* Picocas does not currently pattern-match the Collect target;
     * Derivative[_][f][_] is treated as a literal expression, leaving
     * the higher-order chain rule output un-grouped. Pin the actual
     * 5-term form so the eventual pattern-based Collect lands
     * deliberately. */
    assert_eval_eq(
        "e = D[f[Sqrt[x^2 + 1]], {x, 3}]; Collect[e, Derivative[_][f][_]]",
        "-3 (x Derivative[1][f][Sqrt[1 + x^2]])/(1 + x^2)^(3/2)"
        " + 3 (x Derivative[2][f][Sqrt[1 + x^2]])/(1 + x^2)"
        " + 3 (x^3 Derivative[1][f][Sqrt[1 + x^2]])/(1 + x^2)^(5/2)"
        " - 3 (x^3 Derivative[2][f][Sqrt[1 + x^2]])/(1 + x^2)^2"
        " + (x^3 Derivative[3][f][Sqrt[1 + x^2]])/(1 + x^2)^(3/2)",
        0);
}

/* --- 14 user-supplied "passing" cases ------------------------------- */
static void test_pass_corpus(void) {
    /* P1: Multi-target Times grouping. */
    assert_eval_eq(
        "Collect[x^2*y^4 + z*(x*y^2)^2 + z + a*z, {x*y^2, z}]",
        "x^2 y^4 (1 + z) + z (1 + a)", 0);
    /* P2: Integer-power target. */
    assert_eval_eq("Collect[a^2*(a^2 + 1), a^2]",
                   "a^2 + a^4", 0);
    /* P3: Derivative as exact target. */
    assert_eval_eq("Collect[a*D[f[x], x] + b*D[f[x], x], D[f[x], x]]",
                   "Derivative[1][f][x] (a + b)", 0);
    /* P4: Derivative target with mixed terms. */
    assert_eval_eq(
        "Collect[f[x] + f[x]*D[f[x], x] + x*D[f[x], x]*f[x], D[f[x], x]]",
        "f[x] + Derivative[1][f][x] (f[x] + x f[x])", 0);
    /* P5: Negative-power target (1/f[x]) inside coefficients. */
    assert_eval_eq(
        "Collect[1/f[x] + 1/f[x]*D[f[x], x] + x*D[f[x], x]/f[x], D[f[x], x]]",
        "1/f[x] + Derivative[1][f][x] (1/f[x] + x/f[x])", 0);
    /* P6: Same after Expand. */
    assert_eval_eq(
        "Collect[Expand[1/f[x] + 1/f[x]*D[f[x], x] + x*D[f[x], x]/f[x]], "
        "D[f[x], x]]",
        "1/f[x] + Derivative[1][f][x] (1/f[x] + x/f[x])", 0);
    /* P7: Polynomial expansion. */
    assert_eval_eq("Collect[(x + a + 1)^3, x]",
                   "1 + 3 a + 3 a^2 + a^3 + x^3 + x (3 + 6 a + 3 a^2)"
                   " + x^2 (3 + 3 a)", 0);
    /* P8: Simple polynomial. */
    assert_eval_eq(
        "Collect[a*x^4 + b*x^4 + 2*a^2*x - 3*b*x + x - 7, x]",
        "-7 + x (1 + 2 a^2 - 3 b) + x^4 (a + b)", 0);
    /* P9: Mixed Sqrt[x] + x^(2/3) + x. */
    assert_eval_eq(
        "Collect[a*Sqrt[x] + Sqrt[x] + x^(2/3) - c*x + 3*x"
        " - 2*b*x^(2/3) + 5,  x]",
        "5 + x (3 - c) + Sqrt[x] (1 + a) + x^(2/3) (1 - 2 b)", 0);
    /* P10: Multivariate cubic. */
    assert_eval_eq(
        "Collect[(x*y + x*z + y*z + x + y)^3, {x, y}]",
        "x (y^2 (3 + 9 z + 9 z^2 + 3 z^3) + y^3 (3 + 6 z + 3 z^2))"
        " + x^2 (y (3 + 9 z + 9 z^2 + 3 z^3) + y^2 (6 + 12 z + 6 z^2)"
        " + y^3 (3 + 3 z))"
        " + x^3 (1 + y^3 + 3 z + 3 z^2 + z^3 + y (3 + 6 z + 3 z^2)"
        " + y^2 (3 + 3 z))"
        " + y^3 (1 + 3 z + 3 z^2 + z^3)", 0);
    /* P11: Larger multivariate cubic. */
    assert_eval_eq(
        "Collect[Expand[(y + 2*x + 3*x*y + 4*x*z + 5*y*z)^3], {x, y}]",
        "x (y^2 (6 + 72 z + 270 z^2 + 300 z^3)"
        " + y^3 (9 + 90 z + 225 z^2))"
        " + x^2 (y (12 + 108 z + 288 z^2 + 240 z^3)"
        " + y^2 (36 + 252 z + 360 z^2)"
        " + y^3 (27 + 135 z))"
        " + x^3 (8 + 27 y^3 + 48 z + 96 z^2 + 64 z^3"
        " + y (36 + 144 z + 144 z^2)"
        " + y^2 (54 + 108 z))"
        " + y^3 (1 + 15 z + 75 z^2 + 125 z^3)", 0);
    /* P13: Higher-order derivative as target. */
    assert_eval_eq(
        "Collect[a*D[f[x],{x,2}] + b*D[f[x],{x,2}], D[f[x],{x,2}]]",
        "Derivative[2][f][x] (a + b)", 0);
    /* P14: Nested radical pattern target. */
    assert_eval_eq(
        "Collect[(Sqrt[15 + 5*Sqrt[2]]*x + Sqrt[3 + Sqrt[2]]*y)*2, "
        "Sqrt[_]]",
        "2 Sqrt[3 + Sqrt[2]] y + 2 Sqrt[15 + 5 Sqrt[2]] x", 0);
}

/* P12 uses Series, which the parser/REPL accepts but the inputs in
 * test files cannot easily span the multi-statement form. Inline it
 * here using the post-Normal closed-form expression instead. */
static void test_pass_p12_normal_series(void) {
    /* The post-Normal[Series[...]] form can be written directly. */
    assert_eval_eq(
        "Collect[Normal[Series[Sin[a + b], {b, 0, 10}]],"
        " {Sin[a], Cos[a]}]",
        "Cos[a] (b - 1/6 b^3 + 1/120 b^5 - 1/5040 b^7 + 1/362880 b^9)"
        " + Sin[a] (1 - 1/2 b^2 + 1/24 b^4 - 1/720 b^6 + 1/40320 b^8"
        " - 1/3628800 b^10)", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_bug1_symbolic_exp_ratio);
    TEST(test_bug2_pattern_target);
    TEST(test_pass_corpus);
    TEST(test_pass_p12_normal_series);

    printf("All collect_corpus tests passed!\n");
    return 0;
}
