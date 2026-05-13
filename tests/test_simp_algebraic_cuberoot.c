/*
 * test_simp_algebraic_cuberoot.c
 * -------------------------------
 * Tests for `Simplify` on cube-root and higher-order algebraic
 * expressions.  These exercise the Phase G9 / Phase C path in
 * `simp_algebraic` that routes single-rational-base radical generators
 * (q ≥ 2) through `Together[..., Extension -> α]` and the qaupoly
 * substrate.
 *
 * Before Phase C, `simp_algebraic` bailed on any `Power[c, p/q]` with
 * q != 2 (the q=2 multi-Sqrt rationalisation by sign-flip conjugate
 * doesn't generalise to cube roots), and inputs like
 * `Power[2, 1/3] * Power[2, 2/3]` stayed un-simplified.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "core.h"
#include "symtab.h"
#include "test_utils.h"

/* =========================== Cube roots =========================== */

static void test_cuberoot_product_collapses(void) {
    /* 2^(1/3) * 2^(2/3) = 2^(3/3) = 2.  The picocas evaluator already
     * collapses this at multiplication time (Power[2, 1/3] * Power[2,
     * 2/3] folds via Times-canonicalisation), so this is a sanity
     * check that Simplify doesn't disturb the canonical answer. */
    assert_eval_eq("Simplify[Power[2, 1/3] Power[2, 2/3]]", "2", 0);
}

static void test_cuberoot_minpoly_identity(void) {
    /* α^2 + α + 1)(α - 1) = α^3 - 1 = 2 - 1 = 1   where α = 2^(1/3).
     * Mathematica's Simplify collapses this to 1; the Phase C qaupoly
     * shortcut on simp_algebraic does the same via
     * Together[..., Extension -> 2^(1/3)]. */
    assert_eval_eq(
        "Simplify[(1 + Power[2, 1/3] + Power[2, 2/3]) (Power[2, 1/3] - 1)]",
        "1", 0);
}

static void test_cuberoot_polynomial_factor(void) {
    /* (x^3 - 2)/(x - 2^(1/3)) = x^2 + 2^(1/3) x + 2^(2/3) over Q(2^(1/3)).
     * The cancellation goes through Cancel-with-Extension under the
     * hood (auto-detected). */
    assert_eval_eq(
        "Simplify[(x^3 - 2)/(x - Power[2, 1/3])]",
        "2^(2/3) + 2^(1/3) x + x^2", 0);
}

/* =========================== Quartic / quintic roots =========================== */

static void test_quartic_root_product_collapses(void) {
    /* 5^(1/4) * 5^(3/4) = 5^(4/4) = 5. */
    assert_eval_eq("Simplify[Power[5, 1/4] Power[5, 3/4]]", "5", 0);
}

static void test_quartic_root_squared(void) {
    /* (5^(1/4))^2 = 5^(1/2) = Sqrt[5]. */
    assert_eval_eq("Simplify[Power[Power[5, 1/4], 2]]", "Sqrt[5]", 0);
}

/* =========================== Sqrt regression =========================== */
/* The Phase C shortcut must not break the existing q=2 path. */

static void test_sqrt_product_unchanged(void) {
    /* Sqrt[2] * Sqrt[2] = 2.  picocas's Power[2, 1/2]*Power[2, 1/2]
     * folds at the Times level, so this is a sanity check. */
    assert_eval_eq("Simplify[Sqrt[2] Sqrt[2]]", "2", 0);
}

static void test_sqrt_rational_input(void) {
    /* (x^2 - 2)/(x - Sqrt[2]) over Q(Sqrt[2]).  This was already
     * handled by the q=2 path; confirm the Phase C shortcut produces
     * the same result. */
    assert_eval_eq(
        "Simplify[(x^2 - 2)/(x - Sqrt[2])]",
        "Sqrt[2] + x", 0);
}

/* =========================== Bail paths =========================== */

static void test_no_radical_simplify_unchanged(void) {
    /* No algebraic generator — simp_algebraic should be a no-op. */
    assert_eval_eq("Simplify[x^2 + 2 x + 1]", "(1 + x)^2", 0);
}

static void test_polynomial_radicand_falls_back(void) {
    /* Sqrt[x + 1]: polynomial radicand, the qaupoly shortcut returns
     * NULL (auto-detect bails on non-integer base), but the q=2 multi-
     * Sqrt path still applies.  Result: input unchanged because there's
     * nothing to rationalise. */
    assert_eval_eq("Simplify[Sqrt[x + 1]]", "Sqrt[1 + x]", 0);
}

/* =========================== Nested radicals (Phase E) =========================== */

static void test_nested_sqrt_squared(void) {
    /* Sqrt[1 + Sqrt[2]]^2 = 1 + Sqrt[2].  Auto-detect surfaces the
     * nested radical as a GEN_NESTED generator; Together-with-Extension
     * recognises that Sqrt[base]^2 = base and folds. */
    assert_eval_eq("Simplify[Sqrt[1 + Sqrt[2]]^2]", "1 + Sqrt[2]", 0);
}

static void test_sqrt_of_cbrt_squared(void) {
    /* (Sqrt[2^(1/3)])^2 = 2^(1/3). */
    assert_eval_eq("Simplify[(Sqrt[Power[2, 1/3]])^2]", "2^(1/3)", 0);
}

static void test_sqrt_of_cbrt_quotient_squared(void) {
    /* Sqrt[2^(1/3)/6]^2 = 2^(1/3)/6  (which picocas renders as
     * 1/3 · 2^(-2/3) thanks to its rational-exponent canonicalisation
     * of 2^(1/3)/6 = 2^(1/3 - log_2 6) = ...; verify it's a 2^(-2/3)
     * form proportional to 1/3). */
    assert_eval_eq("Simplify[Sqrt[Power[2, 1/3]/6] Sqrt[Power[2, 1/3]/6]]",
                   "1/3/2^(2/3)", 0);
}

/* =========================== Phase F: canonical-dedup =========================== */

static void test_canon_dedup_equal_radicands_subtract(void) {
    /* Two mathematically-equal nested Sqrts with structurally distinct
     * radicands.  Both equal Sqrt[2^(1/3)/6] but appear as
     *   Sqrt[-1/(9·2^(2/3)) + (2/9)·2^(1/3)]  and
     *   Sqrt[1/(3·2^(2/3))]
     * The canonicalise_post pass surfaces 2 as an INT_BASE generator
     * (with q_lcm=3 from the radicand's 2^(2/3)/2^(1/3) sub-expressions),
     * canonicalises each radicand via Together-with-Extension to
     * 1/3·2^(-2/3), then dedups them.  Simplify uses this to drive the
     * difference to zero. */
    assert_eval_eq(
        "Simplify[Sqrt[-1/9/2^(2/3) + 2/9 2^(1/3)] - Sqrt[1/3/2^(2/3)]]",
        "0", 0);
}

/* =========================== Phase G: Power[c, p/q] general p =========================== */

static void test_power_neg_pq_lift(void) {
    /* Power[2, -2/3] is now accepted by qa_resolve_extension's general-p
     * branch.  Cancel with explicit Extension -> Power[2, 1/3] lifts
     * the user's input through expand_radicals_to_atomic_poly, rewriting
     * Power[2, -2/3] into a polynomial in 2^(1/3) before the QAUPoly
     * lift.  The expected output structure depends on picocas's Times
     * canonicaliser; here we just sanity-check that simplification
     * succeeds (result != input) and contains 2^(1/3) but not
     * 2^(-2/3). */
    assert_eval_eq(
        "Cancel[Power[2, -2/3] x^2 - x/Power[2, 1/3], Extension -> Power[2, 1/3]]",
        "1/2 (-2^(2/3) x + 2^(1/3) x^2)", 0);
}

static void test_multigen_qgamma_constant_collapse(void) {
    /* (Sqrt[2] + Sqrt[3])(Sqrt[2] - Sqrt[3]) = 2 - 3 = -1.  Q(γ)-
     * constant case: qa_cancel_with_tower's no-variable branch lifts to
     * a QANum via PolynomialRemainder + qa_div, renders the canonical
     * coefs against t->gamma_render, and the result reduces to -1.
     * Previously this returned a degree-16 polynomial in γ. */
    assert_eval_eq(
        "Cancel[(Sqrt[2] + Sqrt[3]) (Sqrt[2] - Sqrt[3]), Extension -> Automatic]",
        "-1", 0);
}

static void test_together_equal_nested_radicals(void) {
    /* Two mathematically-equal nested Sqrts via Together with
     * Extension -> Automatic: Sqrt[-1/9/2^(2/3) + 2/9·2^(1/3)] -
     * Sqrt[1/3/2^(2/3)] -> 0.  The Together path canonicalises both
     * radicands to 1/3/2^(2/3) (via Phase F's canonicalise_post +
     * positive-exponent normalisation) and the resulting Sqrt[X] -
     * Sqrt[X] collapses to 0 before any tower lift even runs.  Then
     * the tower's no-variable Q(γ)-arithmetic path renders 0 as 0. */
    assert_eval_eq(
        "Together[Sqrt[-1/9/2^(2/3) + 2/9 2^(1/3)] - Sqrt[1/3/2^(2/3)], Extension -> Automatic]",
        "0", 0);
}

static void test_together_qgamma_constant_canonical_linear_basis(void) {
    /* Q(γ)-constants previously returned an unsimplified γ-polynomial
     * surface form after the qa_cancel_with_tower no-variable branch
     * substituted γ_internal → t->gamma_render and evaluated.  For
     * γ = Sqrt[2] + Sqrt[3], `1/(Sqrt[2]+Sqrt[3]) + 1/(Sqrt[2]-Sqrt[3])`
     * came back as `-(Sqrt[2]+Sqrt[3])^3 + 9(Sqrt[2]+Sqrt[3])`
     * (mathematically equal to `-2 Sqrt[2]` but verbose).
     *
     * Fix: after substitution, run Expand to collapse the γ-polynomial
     * powers, then run Together (no extension) to fold radical
     * denominators via Sqrt-base polynomial GCD.  The Together pass
     * also catches arithmetically-zero spurious nonzero forms (e.g.
     * the (Sqrt[2]+Sqrt[3])^2 - 5 - 2 Sqrt[6] case that Expand alone
     * leaves as `-210969/184/Sqrt[6] + 70323/368 Sqrt[6]`). */
    assert_eval_eq(
        "Together[1/(Sqrt[2] + Sqrt[3]) + 1/(Sqrt[2] - Sqrt[3]), Extension -> Automatic]",
        "-2 Sqrt[2]", 0);
    assert_eval_eq(
        "Together[1/(Sqrt[2] + Sqrt[3]) - 1/(Sqrt[2] - Sqrt[3]), Extension -> Automatic]",
        "2 Sqrt[3]", 0);
    assert_eval_eq(
        "Cancel[(Sqrt[2]+Sqrt[3])^2 - 5 - 2 Sqrt[6], Extension -> Automatic]",
        "0", 0);
    assert_eval_eq(
        "Cancel[1/(Sqrt[2] + Sqrt[3]), Extension -> Automatic]",
        "Sqrt[3] - Sqrt[2]", 0);
    assert_eval_eq(
        "Cancel[(Sqrt[2] + Sqrt[3])/(Sqrt[2] - Sqrt[3]), Extension -> Automatic]",
        "-5 - 2 Sqrt[6]", 0);
    assert_eval_eq(
        "Together[(Sqrt[2] - Sqrt[3])^2 + (Sqrt[2] + Sqrt[3])^2, Extension -> Automatic]",
        "10", 0);
}

static void test_together_headline_d_integrate(void) {
    /* The motivating multi-generator example:
     *
     *   Together[D[Integrate[a x/(x^3 + 2), x], x], Extension -> Automatic]
     *
     * D[Integrate[a x/(x^3 + 2), x], x] produces a 3-term sum involving
     * 2^(1/3), 2^(2/3), Sqrt[3] and a nested radical
     * Sqrt[-1/9/2^(2/3) + 2/9·2^(1/3)] (which equals 1/(2^(1/3) Sqrt[3])).
     * The combined-but-not-cancelled form should collapse mathematically
     * to the integrand a x / (x^3 + 2).
     *
     * Implementation path: extension_autodetect builds an n=2 tower over
     * (2^(1/3), Sqrt[3]); qa_cancel_with_tower's Step 1 substitution
     * blows up the leaf count past the safety gate and returns NULL.
     * The builtin_together fallback then tries each tower generator as
     * a single-α extension via together_recursive_ext, runs a final
     * no-extension Together fold-up, and picks the smallest result. */
    assert_eval_eq(
        "Together[D[Integrate[a x/(x^3+2), x], x], Extension -> Automatic]",
        "(a x)/(2 + x^3)", 0);
    assert_eval_eq(
        "Cancel[D[Integrate[a x/(x^3+2), x], x], Extension -> Automatic]",
        "(a x)/(2 + x^3)", 0);
}

static void test_simplify_headline_d_integrate(void) {
    /* The motivating Simplify case (no Extension option).  Simplify must
     * reach (a x)/(x^3 + 2) on its own by detecting the multi-generator
     * algebraic-number Plus and dispatching the input through
     * Together[..., Extension -> Automatic] as a top-level fast path
     * before simp_bottomup's per-subnode descent.  Without that fast
     * path the search hangs (3-term sum over 2^(1/3), Sqrt[3] and a
     * nested radical never collapses through any single transform that
     * simp_bottomup tries).  Tests the failure mode that previously
     * emitted spurious Power::infy: 1/0 and ComplexInfinity. */
    assert_eval_eq(
        "Simplify[D[Integrate[a x/(x^3+2), x], x]]",
        "(a x)/(2 + x^3)", 0);
}

static void test_autodetect_absorbs_canonical_radical(void) {
    /* extension_autodetect previously bailed (returned NULL) when
     * Phase F's canonicalise_post simplified a GEN_NESTED's surface
     * form into a Times of inverse-powers of existing integer-base
     * generators (e.g. Sqrt[-1/9/2^(2/3) + 2/9·2^(1/3)] canonicalises
     * to 1/(2^(1/3) Sqrt[3]) — no longer a single algebraic
     * generator).  Now: the post-canon step harvests integer-base
     * components from the simplified surface (via re-walking it with
     * autodetect_walk) and marks the original entry as GEN_ABSORBED;
     * dedup drops it.  The tower then builds successfully from the
     * remaining integer-base generators.
     *
     * Smoke test: a Cancel input whose autodetected tower previously
     * failed to build now succeeds and the cancellation runs. */
    assert_eval_eq(
        "Cancel[(Sqrt[-1/9/2^(2/3) + 2/9 2^(1/3)]) (Sqrt[3]) (2^(1/3)) - 1, Extension -> Automatic]",
        "0", 0);
}

/* =================== Cube-root denesting via simp_cuberoot =====================
 *
 * Before the nested-radical-penalty addition, simp_cuberoot's `(a+b Sqrt[c])^(1/3)
 * -> p + q*Sqrt[c]` rewrite was rejected by the simp_search complexity gate
 * whenever the denested form had a higher LeafCount than the input. The penalty
 * (a +3 surcharge for Power[Compound, Rational[_, q]] when the compound base
 * itself contains another root) makes the structurally-denested form preferred
 * by score, matching Mathematica's behaviour.
 *
 * These cases all hit Pattern A: small (a, b, c) with rational (p, q).
 */

static void test_cuberoot_denest_rational_pq(void) {
    /* (2 + Sqrt[5])^(1/3) = (1 + Sqrt[5])/2.  ((1+Sqrt[5])/2)^3 expands
     * to 2 + Sqrt[5].  The denested LeafCount (13) is higher than the
     * nested input's (11), so this case depends on the nested-radical
     * penalty firing. */
    assert_eval_eq(
        "Simplify[(2 + Sqrt[5])^(1/3)]",
        "1/2 (1 + Sqrt[5])", 0);
}

static void test_cuberoot_denest_integer_pq(void) {
    /* (20 + 14 Sqrt[2])^(1/3) = 2 + Sqrt[2]. */
    assert_eval_eq(
        "Simplify[(20 + 14 Sqrt[2])^(1/3)]",
        "2 + Sqrt[2]", 0);
}

static void test_cuberoot_denest_integer_pq_sqrt3(void) {
    /* (10 + 6 Sqrt[3])^(1/3) = 1 + Sqrt[3]. */
    assert_eval_eq(
        "Simplify[(10 + 6 Sqrt[3])^(1/3)]",
        "1 + Sqrt[3]", 0);
}

static void test_cuberoot_sum_of_conjugates(void) {
    /* (2 + Sqrt[5])^(1/3) + (2 - Sqrt[5])^(1/3) = 1 under the real-
     * cube-root convention (Cardano discriminant test). */
    assert_eval_eq(
        "Simplify[(2 + Sqrt[5])^(1/3) + (2 - Sqrt[5])^(1/3)]",
        "1", 0);
}

static void test_cuberoot_recip_via_denest(void) {
    /* 1/(2+Sqrt[5])^(1/3) - 2/(1+Sqrt[5]) == 0.  After denesting,
     * 1/((1+Sqrt[5])/2) = 2/(1+Sqrt[5]). */
    assert_eval_eq(
        "Simplify[1/(2+Sqrt[5])^(1/3) - 2/(1+Sqrt[5])]",
        "0", 0);
}

static void test_cuberoot_inert_when_no_denest(void) {
    /* (3 + Sqrt[2])^(1/3) has no rational denesting (a^2 - b^2 c = 7
     * is not a perfect cube), so the input is returned unchanged. */
    assert_eval_eq(
        "Simplify[(3 + Sqrt[2])^(1/3)]",
        "(3 + Sqrt[2])^(1/3)", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_cuberoot_product_collapses);
    TEST(test_cuberoot_minpoly_identity);
    TEST(test_cuberoot_polynomial_factor);

    TEST(test_cuberoot_denest_rational_pq);
    TEST(test_cuberoot_denest_integer_pq);
    TEST(test_cuberoot_denest_integer_pq_sqrt3);
    TEST(test_cuberoot_sum_of_conjugates);
    TEST(test_cuberoot_recip_via_denest);
    TEST(test_cuberoot_inert_when_no_denest);

    TEST(test_quartic_root_product_collapses);
    TEST(test_quartic_root_squared);

    TEST(test_sqrt_product_unchanged);
    TEST(test_sqrt_rational_input);

    TEST(test_no_radical_simplify_unchanged);
    TEST(test_polynomial_radicand_falls_back);

    TEST(test_nested_sqrt_squared);
    TEST(test_sqrt_of_cbrt_squared);
    TEST(test_sqrt_of_cbrt_quotient_squared);

    TEST(test_canon_dedup_equal_radicands_subtract);

    TEST(test_power_neg_pq_lift);
    TEST(test_multigen_qgamma_constant_collapse);
    TEST(test_together_equal_nested_radicals);
    TEST(test_together_qgamma_constant_canonical_linear_basis);
    TEST(test_together_headline_d_integrate);
    TEST(test_simplify_headline_d_integrate);
    TEST(test_autodetect_absorbs_canonical_radical);

    printf("All simp_algebraic_cuberoot tests passed!\n");
    return 0;
}
