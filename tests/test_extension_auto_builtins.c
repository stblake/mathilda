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

/* =========================== Phase D: multivariate tower GCD/LCM =========================== */

/* Phase D extends the tower path to inputs with more than one polynomial
 * variable beyond γ.  The strategy is to substitute each α_i with its
 * γ-polynomial form and call the no-extension multivariate
 * `PolynomialGCD` / `PolynomialLCM`, then substitute γ back through
 * `t->gamma_render` and Expand + evaluate to canonicalise.
 *
 * Correctness boundary: the result is the Q[γ, x, y, ...]-GCD
 * (γ treated as polynomial variable).  This is a Q(γ)-common-divisor
 * but possibly not the maximal one.  Tests below distinguish:
 *   - cases where Q[γ,...]-GCD equals Q(γ)[...]-GCD (most practical
 *     inputs);
 *   - cases where they differ — those return a non-maximal GCD whose
 *     output is still mathematically valid.
 */

/* --- Multivariate single-α: a Q[γ,x,y] GCD that agrees with Q(γ)[x,y] --- */
static void test_pgcd_multivar_single_alpha_shared_no_gamma(void) {
    /* gcd((x^2 + y^2)(x - Sqrt[2]), (x^2 + y^2)(x + Sqrt[2])) = x^2 + y^2.
     * The shared factor is γ-free, so Q[γ,x,y]-GCD finds it exactly. */
    assert_eval_eq(
        "PolynomialGCD[(x^2 + y^2) (x - Sqrt[2]), (x^2 + y^2) (x + Sqrt[2]), "
        "Extension -> Automatic]",
        "x^2 + y^2", 0);
}

static void test_pgcd_multivar_single_alpha_shared_with_gamma(void) {
    /* gcd((a Sqrt[2] + 1) x, (a Sqrt[2] + 1) y) = a Sqrt[2] + 1.
     * The shared factor is linear in γ (Sqrt[2]); Q[γ,a,x,y]-GCD
     * finds it. */
    assert_eval_eq(
        "PolynomialGCD[(a Sqrt[2] + 1) x, (a Sqrt[2] + 1) y, "
        "Extension -> Automatic]",
        "1 + Sqrt[2] a", 0);
}

static void test_pgcd_multivar_single_alpha_coprime(void) {
    /* gcd(x y - Sqrt[2], x y + Sqrt[2]) = 1 in Q(Sqrt[2])[x,y]; the
     * resultant in (xy) gives a nonzero constant 2 Sqrt[2]. */
    assert_eval_eq(
        "PolynomialGCD[x y - Sqrt[2], x y + Sqrt[2], Extension -> Automatic]",
        "1", 0);
}

/* --- Multivariate multi-α: γ-free GCD --- */
static void test_pgcd_multivar_multi_alpha_gamma_free(void) {
    /* gcd((a + b)(Sqrt[2] + 1), (a + b)(Sqrt[2] - 1)) = a + b. */
    assert_eval_eq(
        "PolynomialGCD[(a + b) (Sqrt[2] + 1), (a + b) (Sqrt[2] - 1), "
        "Extension -> Automatic]",
        "a + b", 0);
}

static void test_pgcd_multivar_multi_alpha_two_vars_unit_factor(void) {
    /* gcd((x^2 - y^2)(Sqrt[2] + Sqrt[3]), (x^2 - y^2)(Sqrt[2] - Sqrt[3])).
     *
     * The canonical Q(Sqrt[2], Sqrt[3])[x,y]-GCD is x^2 - y^2 — the
     * factors (Sqrt[2]+Sqrt[3]) and (Sqrt[2]-Sqrt[3]) are units in
     * Q(γ) (their product is 2 - 3 = -1, both nonzero), so neither
     * contributes to the canonical GCD.
     *
     * Phase D limitation: the substitute-back approach computes the
     * Q[γ, x, y]-GCD (γ treated as polynomial variable).  Under the
     * primitive-element substitution γ = Sqrt[2] + Sqrt[3] one input
     * becomes (x^2 - y^2)·γ and the other (x^2 - y^2)·γ·(γ^2 - 10),
     * so the Q[γ,x,y]-GCD is (x^2 - y^2)·γ — which substitutes back
     * to (x^2 - y^2)(Sqrt[2] + Sqrt[3]).  This is still a Q(γ)-divisor
     * of both inputs (Sqrt[2] + Sqrt[3] is a unit in Q(γ)), just not
     * the canonical monic form.  Documented boundary. */
    assert_eval_eq(
        "PolynomialGCD[(x^2 - y^2) (Sqrt[2] + Sqrt[3]), "
        "(x^2 - y^2) (Sqrt[2] - Sqrt[3]), Extension -> Automatic]",
        "Sqrt[2] x^2 + Sqrt[3] x^2 - Sqrt[2] y^2 - Sqrt[3] y^2", 0);
}

/* --- Multivariate LCM: γ-free common factor --- */
static void test_plcm_multivar_single_alpha_distinct_linears(void) {
    /* lcm(a + Sqrt[2], b + Sqrt[2]) = (a + Sqrt[2])(b + Sqrt[2]) in
     * Q(γ)[a, b], expanded by Phase D's canonicalisation pass into
     * `a b + (a + b) Sqrt[2] + 2`.  Q(γ)-coprime, so LCM equals the
     * product. */
    assert_eval_eq(
        "PolynomialLCM[a + Sqrt[2], b + Sqrt[2], Extension -> Automatic]",
        "2 + Sqrt[2] a + Sqrt[2] b + a b", 0);
}

static void test_plcm_multivar_share_in_extension(void) {
    /* lcm((a + Sqrt[2])(x - 1), (a + Sqrt[2])(x + 1)) =
     * (a + Sqrt[2])(x^2 - 1) in Q(γ)[a, x], expanded by Phase D's
     * canonicalisation into `(a + Sqrt[2])(x^2 - 1)` distributed. */
    assert_eval_eq(
        "PolynomialLCM[(a + Sqrt[2]) (x - 1), (a + Sqrt[2]) (x + 1), "
        "Extension -> Automatic]",
        "-Sqrt[2] - a + Sqrt[2] x^2 + a x^2", 0);
}

static void test_plcm_multivar_multi_alpha_no_share(void) {
    /* lcm(Sqrt[2] + x, Sqrt[3] + y) — fully coprime in
     * Q(Sqrt[2], Sqrt[3])[x, y].  LCM is the product, expanded:
     * (Sqrt[2] + x)(Sqrt[3] + y) = Sqrt[6] + Sqrt[3] x + Sqrt[2] y + x y. */
    assert_eval_eq(
        "PolynomialLCM[Sqrt[2] + x, Sqrt[3] + y, Extension -> Automatic]",
        "Sqrt[6] + Sqrt[3] x + Sqrt[2] y + x y", 0);
}

/* --- Three-argument multivariate --- */
static void test_pgcd_multivar_three_arg_share(void) {
    /* gcd of three polynomials each carrying the (x^2 + Sqrt[2]) factor. */
    assert_eval_eq(
        "PolynomialGCD["
        "(x^2 + Sqrt[2]) (a + 1), "
        "(x^2 + Sqrt[2]) (a + 2), "
        "(x^2 + Sqrt[2]) (a + 3), "
        "Extension -> Automatic]",
        "Sqrt[2] + x^2", 0);
}

/* --- No-op and pass-through cases --- */
static void test_pgcd_multivar_no_radical_unchanged(void) {
    /* Multivariate but no algebraic radicals: tower auto-detect is a
     * no-op, the standard multivariate path runs. */
    assert_eval_eq(
        "PolynomialGCD[(x + y) (a + b), (x + y) (c + d), Extension -> Automatic]",
        "x + y", 0);
}

static void test_plcm_multivar_no_extension_unchanged(void) {
    /* Without Extension -> Automatic the multivariate radical case
     * goes through the no-extension BPList path, which returns the
     * literal product as a single expression.  Phase D must not
     * change this default. */
    assert_eval_eq(
        "PolynomialLCM[a + Sqrt[2], b + Sqrt[2]]",
        "(Sqrt[2] + a) (Sqrt[2] + b)", 0);
}

/* --- Cube-root multi-α --- */
static void test_pgcd_multivar_cbrt_share(void) {
    /* gcd((x + 2^(1/3))(a + 1), (x + 2^(1/3))(a + 2)) = x + 2^(1/3). */
    assert_eval_eq(
        "PolynomialGCD["
        "(x + Power[2, 1/3]) (a + 1), "
        "(x + Power[2, 1/3]) (a + 2), "
        "Extension -> Automatic]",
        "2^(1/3) + x", 0);
}

static void test_plcm_multivar_cbrt_distinct(void) {
    /* lcm(x + 2^(1/3), x + 3^(1/3)) over Q(2^(1/3), 3^(1/3))[x,a,...]:
     * coprime, so LCM = product, expanded into
     * `x^2 + (2^(1/3) + 3^(1/3)) x + 6^(1/3)`. */
    assert_eval_eq(
        "PolynomialLCM[x + Power[2, 1/3], x + Power[3, 1/3], "
        "Extension -> Automatic]",
        "6^(1/3) + 2^(1/3) x + 3^(1/3) x + x^2", 0);
}

/* --- Mixed Sqrt + cube root --- */
static void test_pgcd_multivar_sqrt_cbrt_mix(void) {
    /* Two distinct algebraic generators Sqrt[2] and 3^(1/3) appearing in
     * the same multivariate inputs.  Shared factor is γ-free. */
    assert_eval_eq(
        "PolynomialGCD["
        "(x y + 1) (Sqrt[2] + Power[3, 1/3]), "
        "(x y + 1) (Sqrt[2] - Power[3, 1/3]), "
        "Extension -> Automatic]",
        "1 + x y", 0);
}

/* --- Three poly variables --- */
static void test_pgcd_multivar_three_poly_vars(void) {
    /* gcd((a+b+c)(x+Sqrt[2]), (a+b+c)(x-Sqrt[2])) — three free poly
     * vars besides γ.  Canonical Q(γ)-GCD is (a + b + c); the
     * (x±Sqrt[2]) factors are Q(γ)-coprime. */
    assert_eval_eq(
        "PolynomialGCD[(a + b + c) (x + Sqrt[2]), (a + b + c) (x - Sqrt[2]), "
        "Extension -> Automatic]",
        "a + b + c", 0);
}

static void test_plcm_multivar_three_poly_vars_expanded(void) {
    /* lcm((x + Sqrt[2])(a - 1), (x + Sqrt[2])(b - 1)) =
     * (x + Sqrt[2])(a - 1)(b - 1), Q(γ)-coprime in (a, b).  Phase D
     * canonicalises into the fully-expanded distributed form. */
    assert_eval_eq(
        "PolynomialLCM[(x + Sqrt[2]) (a - 1), (x + Sqrt[2]) (b - 1), "
        "Extension -> Automatic]",
        "Sqrt[2] - Sqrt[2] a - Sqrt[2] b + Sqrt[2] a b + x - a x - b x + a b x", 0);
}

/* --- Phase D boundary: divisibility-by-γ-only-unit cases --- */
static void test_pgcd_multivar_divisible_by_unit(void) {
    /* One input divides the other.  Q(γ)-GCD = (x - 1)(Sqrt[2] + Sqrt[3]) —
     * Phase D recovers this exactly because the "unit" Sqrt[2]+Sqrt[3]
     * appears as a polynomial-in-γ factor in both inputs (one occurrence
     * each), so the Q[γ,x]-GCD agrees with the canonical Q(γ)[x]-GCD. */
    assert_eval_eq(
        "PolynomialGCD["
        "(x - 1) (Sqrt[2] + Sqrt[3]), "
        "(x - 1) (Sqrt[2] + Sqrt[3]) (Sqrt[2] - Sqrt[3]), "
        "Extension -> Automatic]",
        "-Sqrt[2] - Sqrt[3] + Sqrt[2] x + Sqrt[3] x", 0);
}

static void test_pgcd_multivar_high_gamma_power_collapse(void) {
    /* Sqrt[2]^2 = 2 (Mathilda auto-evaluation), so `x - Sqrt[2]^2` and
     * `x - 2` are structurally identical even before Phase D runs.
     * gcd = x - 2. */
    assert_eval_eq(
        "PolynomialGCD[x - Sqrt[2]^2, x - 2, Extension -> Automatic]",
        "-2 + x", 0);
}

/* --- Inexact input must still take the rationalize/numericalize path --- */
static void test_pgcd_multivar_inexact_unchanged(void) {
    /* Phase D MUST NOT swallow inexact-coefficient inputs — those go
     * through `internal_rationalize_then_numericalize` so the user
     * sees consistent inexact-in / inexact-out behaviour. */
    assert_eval_eq(
        "PolynomialGCD[x - 1.5, x^2 - 2.25, Extension -> Automatic]",
        "-1.5 + x", 0);
}

/* --- Regression: ensure no Extension -> Automatic side-effects --- */
static void test_pgcd_multivar_explicit_none_unchanged(void) {
    /* Extension -> None must NOT trigger the tower path; result matches
     * no-extension form (the symbolic product). */
    assert_eval_same(
        "PolynomialGCD[(a + Sqrt[2]) x, (a + Sqrt[2]) y, Extension -> None]",
        "PolynomialGCD[(a + Sqrt[2]) x, (a + Sqrt[2]) y]");
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

    /* Phase D: multivariate tower GCD/LCM */
    TEST(test_pgcd_multivar_single_alpha_shared_no_gamma);
    TEST(test_pgcd_multivar_single_alpha_shared_with_gamma);
    TEST(test_pgcd_multivar_single_alpha_coprime);
    TEST(test_pgcd_multivar_multi_alpha_gamma_free);
    TEST(test_pgcd_multivar_multi_alpha_two_vars_unit_factor);
    TEST(test_plcm_multivar_single_alpha_distinct_linears);
    TEST(test_plcm_multivar_share_in_extension);
    TEST(test_plcm_multivar_multi_alpha_no_share);
    TEST(test_pgcd_multivar_three_arg_share);
    TEST(test_pgcd_multivar_no_radical_unchanged);
    TEST(test_plcm_multivar_no_extension_unchanged);
    TEST(test_pgcd_multivar_cbrt_share);
    TEST(test_plcm_multivar_cbrt_distinct);
    TEST(test_pgcd_multivar_sqrt_cbrt_mix);
    TEST(test_pgcd_multivar_three_poly_vars);
    TEST(test_plcm_multivar_three_poly_vars_expanded);
    TEST(test_pgcd_multivar_divisible_by_unit);
    TEST(test_pgcd_multivar_high_gamma_power_collapse);
    TEST(test_pgcd_multivar_inexact_unchanged);
    TEST(test_pgcd_multivar_explicit_none_unchanged);

    printf("All extension_auto_builtins tests passed!\n");
    return 0;
}
