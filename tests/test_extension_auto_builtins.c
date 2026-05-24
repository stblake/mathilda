/*
 * test_extension_auto_builtins.c
 * -------------------------------
 * Integration tests for `Extension -> Automatic` on the polynomial and
 * rational-function builtins.  Each test exercises a builtin with an
 * input that contains algebraic radicals (`Power[c, 1/n]`, `Sqrt[c]`)
 * and confirms that `Extension -> Automatic` produces the same result
 * as the explicit `Extension -> α` form.
 *
 * Covered builtins (Phase B of tasks/algebraic_auto.md):
 *   - Cancel
 *   - Together
 *   - PolynomialGCD
 *   - PolynomialLCM
 *   - PolynomialQuotient
 *   - PolynomialRemainder
 *   - PolynomialQuotientRemainder
 *   - Factor
 *   - Apart
 *
 * Each builtin is tested for two things:
 *   (a) Without any Extension option, behaviour is unchanged.
 *   (b) `Extension -> Automatic` produces output equal to the explicit
 *       `Extension -> α` form when a single algebraic generator is
 *       present in the input.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "core.h"
#include "symtab.h"
#include "test_utils.h"

/* Assert that two input expressions, after evaluation, produce the same
 * canonical InputForm string.  Used to compare an `Extension ->
 * Automatic` call to its explicit `Extension -> α` counterpart. */
static void assert_eval_same(const char* a, const char* b) {
    struct Expr* pa = parse_expression(a);
    ASSERT(pa != NULL);
    struct Expr* ea = evaluate(pa);
    expr_free(pa);

    struct Expr* pb = parse_expression(b);
    ASSERT(pb != NULL);
    struct Expr* eb = evaluate(pb);
    expr_free(pb);

    char* sa = expr_to_string(ea);
    char* sb = expr_to_string(eb);
    if (strcmp(sa, sb) != 0) {
        fprintf(stderr, "FAIL: results differ\n  A: %s  →  %s\n  B: %s  →  %s\n",
                a, sa, b, sb);
        free(sa); free(sb);
        expr_free(ea); expr_free(eb);
        exit(1);
    }
    free(sa); free(sb);
    expr_free(ea); expr_free(eb);
}

/* =========================== Cancel =========================== */

static void test_cancel_no_extension(void) {
    /* Sanity: rational input, no algebraic generators — auto-detect must
     * be a no-op. */
    assert_eval_eq("Cancel[(x^2 - 1)/(x - 1), Extension -> Automatic]",
                   "1 + x", 0);
}

static void test_cancel_cbrt2_explicit_vs_auto(void) {
    /* (x^3 - 2)/(x - 2^(1/3))  factors as x^2 + 2^(1/3) x + 2^(2/3)
     * over Q(2^(1/3)).  Both Extension -> 2^(1/3) and Extension ->
     * Automatic should produce that result. */
    assert_eval_same(
        "Cancel[(x^3 - 2)/(x - Power[2, 1/3]), Extension -> Automatic]",
        "Cancel[(x^3 - 2)/(x - Power[2, 1/3]), Extension -> Power[2, 1/3]]");
}

static void test_cancel_sqrt2_explicit_vs_auto(void) {
    assert_eval_same(
        "Cancel[(x^2 - 2)/(x - Sqrt[2]), Extension -> Automatic]",
        "Cancel[(x^2 - 2)/(x - Sqrt[2]), Extension -> Sqrt[2]]");
}

/* =========================== Together =========================== */

static void test_together_no_extension(void) {
    /* No algebraic generators — should still produce the standard form. */
    assert_eval_eq("Together[1/x + 1/y, Extension -> Automatic]",
                   "(x + y)/(x y)", 0);
}

static void test_together_cbrt2_explicit_vs_auto(void) {
    assert_eval_same(
        "Together[1/(1 - Power[2, 1/3]) + 1/(1 + Power[2, 1/3]), Extension -> Automatic]",
        "Together[1/(1 - Power[2, 1/3]) + 1/(1 + Power[2, 1/3]), Extension -> Power[2, 1/3]]");
}

/* =========================== PolynomialGCD =========================== */

static void test_pgcd_no_extension(void) {
    assert_eval_eq("PolynomialGCD[x^2 - 1, x - 1, Extension -> Automatic]",
                   "-1 + x", 0);
}

static void test_pgcd_cbrt2_explicit_vs_auto(void) {
    assert_eval_same(
        "PolynomialGCD[x^3 - 2, x - Power[2, 1/3], Extension -> Automatic]",
        "PolynomialGCD[x^3 - 2, x - Power[2, 1/3], Extension -> Power[2, 1/3]]");
}

static void test_pgcd_sqrt2_explicit_vs_auto(void) {
    assert_eval_same(
        "PolynomialGCD[x^2 - 2, x - Sqrt[2], Extension -> Automatic]",
        "PolynomialGCD[x^2 - 2, x - Sqrt[2], Extension -> Sqrt[2]]");
}

/* =========================== PolynomialLCM =========================== */

static void test_plcm_cbrt2_explicit_vs_auto(void) {
    assert_eval_same(
        "PolynomialLCM[x - Power[2, 1/3], x^2 + Power[2, 1/3] x + Power[2, 2/3], Extension -> Automatic]",
        "PolynomialLCM[x - Power[2, 1/3], x^2 + Power[2, 1/3] x + Power[2, 2/3], Extension -> Power[2, 1/3]]");
}

/* =========================== PolynomialQuotient/Remainder =========================== */

static void test_pq_cbrt2_explicit_vs_auto(void) {
    assert_eval_same(
        "PolynomialQuotient[x^3 - 2, x - Power[2, 1/3], x, Extension -> Automatic]",
        "PolynomialQuotient[x^3 - 2, x - Power[2, 1/3], x, Extension -> Power[2, 1/3]]");
}

static void test_pr_cbrt2_explicit_vs_auto(void) {
    assert_eval_same(
        "PolynomialRemainder[x^3 - 2, x - Power[2, 1/3], x, Extension -> Automatic]",
        "PolynomialRemainder[x^3 - 2, x - Power[2, 1/3], x, Extension -> Power[2, 1/3]]");
}

static void test_pqr_cbrt2_explicit_vs_auto(void) {
    assert_eval_same(
        "PolynomialQuotientRemainder[x^3 - 2, x - Power[2, 1/3], x, Extension -> Automatic]",
        "PolynomialQuotientRemainder[x^3 - 2, x - Power[2, 1/3], x, Extension -> Power[2, 1/3]]");
}

/* =========================== Factor =========================== */

static void test_factor_with_radical_in_input(void) {
    /* x^2 + 2 Sqrt[2] x + 2  =  (x + Sqrt[2])^2 over Q(Sqrt[2]).
     * The input contains Sqrt[2] explicitly, so auto-detect surfaces
     * the generator and Factor splits the polynomial.
     *
     * Note: `Factor[x^3 - 2, Extension -> Automatic]` does NOT split
     * for tier-1 because auto-detect only finds generators *present in
     * the input* — discovering that x^3 - 2 has a root at 2^(1/3) is a
     * root-finding problem, out of scope. */
    assert_eval_same(
        "Factor[x^2 + 2 Sqrt[2] x + 2, Extension -> Automatic]",
        "Factor[x^2 + 2 Sqrt[2] x + 2, Extension -> Sqrt[2]]");
}

/* =========================== Apart =========================== */

static void test_apart_no_extension(void) {
    /* Standard partial-fractions case — no algebraic generators. */
    assert_eval_eq("Apart[1/(x (x + 1)), x, Extension -> Automatic]",
                   "1/x - 1/(1 + x)", 0);
}

static void test_apart_cbrt2_explicit_vs_auto(void) {
    /* 1/((x - 2^(1/3))(x^2 + 2^(1/3) x + 2^(2/3))) has algebraic
     * coefficients that Apart should pull apart over Q(2^(1/3)). */
    assert_eval_same(
        "Apart[1/((x - Power[2, 1/3]) (x^2 + Power[2, 1/3] x + Power[2, 2/3])), x, Extension -> Automatic]",
        "Apart[1/((x - Power[2, 1/3]) (x^2 + Power[2, 1/3] x + Power[2, 2/3])), x, Extension -> Power[2, 1/3]]");
}

/* =========================== Bail-out paths =========================== */

static void test_extension_none_blocks_autodetect(void) {
    /* Extension -> None is explicitly "do not use any extension";
     * auto-detect must NOT run.  We compare to the no-Extension form
     * to confirm the two match. */
    assert_eval_same(
        "Cancel[(x^3 - 2)/(x - Power[2, 1/3]), Extension -> None]",
        "Cancel[(x^3 - 2)/(x - Power[2, 1/3])]");
}

static void test_explicit_alpha_wins_over_automatic(void) {
    /* When both Extension -> Automatic and Extension -> α are given,
     * the rightmost-wins rule of extract_extension_option_full should
     * pick the rightmost setting.  We give Automatic first then
     * explicit α; result must equal explicit-α-only. */
    assert_eval_same(
        "Cancel[(x^3 - 2)/(x - Power[2, 1/3]), Extension -> Automatic, Extension -> Power[2, 1/3]]",
        "Cancel[(x^3 - 2)/(x - Power[2, 1/3]), Extension -> Power[2, 1/3]]");
}

/* =========================== Multi-generator (n >= 2) =========================== */

static void test_multigen_sqrt2_sqrt3_collapse(void) {
    /* (Sqrt[2] + Sqrt[3])(Sqrt[2] - Sqrt[3]) = 2 - 3 = -1.  Two
     * generators (Sqrt[2], Sqrt[3]); auto-detect builds the
     * primitive-element tower Q(γ) with γ = Sqrt[2] + Sqrt[3].
     * `qa_cancel_with_tower` lifts to QAUPoly over Q(γ), runs
     * qaupoly_gcd, and renders back.  The result must equal -1 (after
     * normalization) regardless of how γ is rendered. */
    assert_eval_eq(
        "Cancel[(Sqrt[2] + Sqrt[3]) (Sqrt[2] - Sqrt[3]), Extension -> Automatic]",
        "-1", 0);
}

static void test_multigen_together_no_op(void) {
    /* Together over Q(Sqrt[2], Sqrt[3]) on input that's already in
     * combined form should leave it equivalent.  We compare to the
     * no-extension Together output. */
    assert_eval_same(
        "Together[Sqrt[2] + Sqrt[3]]",
        "Together[Sqrt[2] + Sqrt[3], Extension -> Automatic]");
}

static void test_nested_radical_skip_matches_no_ext(void) {
    /* Layer-0 prefilter regression test.  Input has a nested radical
     * (Sqrt[5 + 2 Sqrt[6]] surfaced via Power[Plus[...], -3/2]) and
     * no free polynomial variable, so the predicate skips the
     * extension_autodetect + tower-build cascade.  Result must match
     * what no-extension Together would produce — the structural
     * combine is the same regardless of the Extension option, the
     * difference is just throughput. */
    assert_eval_same(
        "Together[1/(5 + 2 Sqrt[6])^(3/2) + (Sqrt[2] - Sqrt[3])^3]",
        "Together[1/(5 + 2 Sqrt[6])^(3/2) + (Sqrt[2] - Sqrt[3])^3, "
        "Extension -> Automatic]");
}

/* =========================== Phase C: multi-radical GCD/LCM =========================== */

static void test_pgcd_multigen_sqrt2_sqrt3(void) {
    /* (x^2 - 2)(x - 1) shares (x^2 - 2) with itself; auto-detect finds
     * the Sqrt[2] generator (only one) and routes through the single-α
     * path.  Sanity test that the existing path still works. */
    assert_eval_eq(
        "PolynomialGCD[(x^2 - 2) (x - 1), (x^2 - 2) (x + 1), Extension -> Automatic]",
        "-2 + x^2", 0);
}

static void test_pgcd_multigen_x_minus_sqrt2(void) {
    /* `PolynomialGCD[x^2 - 2, x - Sqrt[2], Extension -> Automatic]`
     * routes through the single-α path (Sqrt[2] is the only generator).
     * Confirms baseline. */
    assert_eval_eq(
        "PolynomialGCD[x^2 - 2, x - Sqrt[2], Extension -> Automatic]",
        "-Sqrt[2] + x", 0);
}

static void test_plcm_multigen_sqrt2_collapse(void) {
    /* (x - Sqrt[2])(x + Sqrt[2]) = x^2 - 2.  Single generator → single-α
     * path.  Sanity. */
    assert_eval_eq(
        "PolynomialLCM[x - Sqrt[2], x + Sqrt[2], Extension -> Automatic]",
        "-2 + x^2", 0);
}

static void test_plcm_multigen_sqrt2_sqrt3(void) {
    /* Phase C headline: PolynomialLCM with two independent algebraic
     * generators Sqrt[2] and Sqrt[3].  Auto-detect builds a deg-4 tower
     * over Q(Sqrt[2], Sqrt[3]); the LCM is (x - Sqrt[2])(x - Sqrt[3]) =
     * x^2 - (Sqrt[2] + Sqrt[3]) x + Sqrt[6].
     *
     * Without Phase C the multi-gen branch of `builtin_polynomiallcm`
     * dropped the tower on the floor and the no-extension BPList path
     * returned the literal (x - Sqrt[2])(x - Sqrt[3]) product without
     * canonicalising it. */
    assert_eval_eq(
        "PolynomialLCM[x - Sqrt[2], x - Sqrt[3], Extension -> Automatic]",
        "Sqrt[6] - Sqrt[2] x - Sqrt[3] x + x^2", 0);
}

static void test_plcm_multigen_nested_pair(void) {
    /* (x - Sqrt[2] - Sqrt[3])(x - Sqrt[2] + Sqrt[3]) = (x - Sqrt[2])^2 - 3
     * = x^2 - 2 Sqrt[2] x + 2 - 3 = x^2 - 2 Sqrt[2] x - 1.
     *
     * The two inputs share Sqrt[2] and Sqrt[3] as generators (linear in
     * each, no nesting).  Phase C's tower lift handles the multi-α case
     * and the qaupoly_gcd-driven LCM returns the canonical form. */
    assert_eval_eq(
        "PolynomialLCM[x - Sqrt[2] - Sqrt[3], x - Sqrt[2] + Sqrt[3], "
        "Extension -> Automatic]",
        "-1 - 2 Sqrt[2] x + x^2", 0);
}

static void test_pgcd_multigen_share_linear_factor(void) {
    /* p1 = (x - Sqrt[3] - Sqrt[2])(x - Sqrt[3] + Sqrt[2])
     *    = x^2 - 2 Sqrt[3] x + 1
     * p2 = x - Sqrt[2] - Sqrt[3]
     * gcd(p1, p2) = x - Sqrt[2] - Sqrt[3] over Q(Sqrt[2], Sqrt[3]). */
    assert_eval_eq(
        "PolynomialGCD[x^2 - 2 Sqrt[3] x + 1, x - Sqrt[3] - Sqrt[2], "
        "Extension -> Automatic]",
        "-Sqrt[2] - Sqrt[3] + x", 0);
}

static void test_plcm_multigen_no_automatic_unchanged(void) {
    /* Without `Extension -> Automatic` the multi-generator inputs go
     * through the no-extension BPList path, which returns the symbolic
     * factored product unevaluated as a single expression.  Confirm
     * Phase C did not change this default behaviour. */
    assert_eval_eq(
        "PolynomialLCM[x - Sqrt[2], x - Sqrt[3]]",
        "(-Sqrt[2] + x) (-Sqrt[3] + x)", 0);
}

static void test_composite_sqrt_coalesces_to_primes(void) {
    /* Layer-2 regression test.  The input contains Sqrt[3], Sqrt[5],
     * and Sqrt[15]; Sqrt[15] is algebraically dependent on Sqrt[3] and
     * Sqrt[5] (Sqrt[15] = Sqrt[3]·Sqrt[5]).  Without the coalesce,
     * extension_autodetect builds a deg-8 tower over the three
     * generators; with the coalesce it builds a deg-4 tower over
     * {Sqrt[3], Sqrt[5]} and qa_cancel_with_tower's input
     * pre-substitution rewrites Sqrt[15] to Sqrt[3]·Sqrt[5].
     *
     * The cancellation result must be 0: expand the square
     * (2 Sqrt[3] + 3 Sqrt[5])^2 = 4·3 + 12 Sqrt[15] + 9·5 = 57 + 12 Sqrt[15]
     * so the input collapses additively to 0. */
    assert_eval_eq(
        "Simplify[(2 Sqrt[3] + 3 Sqrt[5])^2 - 12 - 12 Sqrt[15] - 45]",
        "0", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_cancel_no_extension);
    TEST(test_cancel_cbrt2_explicit_vs_auto);
    TEST(test_cancel_sqrt2_explicit_vs_auto);

    TEST(test_together_no_extension);
    TEST(test_together_cbrt2_explicit_vs_auto);

    TEST(test_pgcd_no_extension);
    TEST(test_pgcd_cbrt2_explicit_vs_auto);
    TEST(test_pgcd_sqrt2_explicit_vs_auto);

    TEST(test_plcm_cbrt2_explicit_vs_auto);

    TEST(test_pq_cbrt2_explicit_vs_auto);
    TEST(test_pr_cbrt2_explicit_vs_auto);
    TEST(test_pqr_cbrt2_explicit_vs_auto);

    TEST(test_factor_with_radical_in_input);

    TEST(test_apart_no_extension);
    TEST(test_apart_cbrt2_explicit_vs_auto);

    TEST(test_extension_none_blocks_autodetect);
    TEST(test_explicit_alpha_wins_over_automatic);

    TEST(test_multigen_sqrt2_sqrt3_collapse);
    TEST(test_multigen_together_no_op);
    TEST(test_nested_radical_skip_matches_no_ext);
    TEST(test_composite_sqrt_coalesces_to_primes);

    /* Phase C: multi-radical PolynomialGCD / PolynomialLCM */
    TEST(test_pgcd_multigen_sqrt2_sqrt3);
    TEST(test_pgcd_multigen_x_minus_sqrt2);
    TEST(test_plcm_multigen_sqrt2_collapse);
    TEST(test_plcm_multigen_sqrt2_sqrt3);
    TEST(test_plcm_multigen_nested_pair);
    TEST(test_pgcd_multigen_share_linear_factor);
    TEST(test_plcm_multigen_no_automatic_unchanged);

    printf("All extension_auto_builtins tests passed!\n");
    return 0;
}
