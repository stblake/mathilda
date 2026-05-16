#include "test_utils.h"
#include "symtab.h"
#include "core.h"

/*
 * Tests for reciprocal nested-radical denesting:
 *
 *   1 / Sqrt[A + b Sqrt[C]]   ->   (Sqrt[P] -/+ Sqrt[Q]) / s
 *
 * where the forward half-sum identity gives
 *   Sqrt[A + b Sqrt[C]] = Sqrt[P] +/- Sqrt[Q]
 * with P = (A + s) / 2, Q = (A - s) / 2, s = Sqrt[A^2 - b^2 C], and the
 * conjugate-rationalised denominator P - Q = s.
 *
 * The implementation lives next to try_denest_sqrt_radicand in simp.c.
 * The walker (simp_denest_sqrt_walk) dispatches on both Sqrt-shaped
 * (Power[.., 1/2]) and inv-Sqrt-shaped (Power[.., -1/2]) subtrees.
 *
 * What we verify:
 *
 *   - Standalone denesting fires when the rationalised form strictly
 *     beats the input on Simplify's leaf-count measure.
 *
 *   - Cross-term cancellation works in every shape: forward +/- recip,
 *     two recips with opposing signs, products that collapse to a
 *     constant, mixed multi-extension expressions.
 *
 *   - The original Q(sqrt(2),sqrt(3)) cancellation case from the user
 *     report:
 *         Simplify[Sqrt[2] - Sqrt[3] + 1/Sqrt[5 + 2 Sqrt[6]]] == 0.
 *
 *   - Non-denestable inputs and ones with negative discriminant remain
 *     untouched (no wrong-branch results, no crashes).
 *
 *   - The path is robust to deeper nesting, mixed positive/negative
 *     b, and idempotency under further Simplify.
 */

/* Local strong-assert wrapper. The shared test_utils.h `assert_eval_eq`
 * uses assert(), which the cmake build silences via NDEBUG; so a
 * mismatch only prints FAIL and the run exits 0. We need a hard signal,
 * so this file routes through a wrapper that calls exit(1) on
 * mismatch. */
static int g_failures = 0;
static void check_eval_eq(const char* input, const char* expected) {
    struct Expr* parsed = parse_expression(input);
    if (!parsed) {
        fprintf(stderr, "FAIL: parse failure for: %s\n", input);
        g_failures++;
        return;
    }
    struct Expr* evaluated = evaluate(parsed);
    expr_free(parsed);
    char* str = expr_to_string(evaluated);
    if (strcmp(str, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  Expected: %s\n  Actual:   %s\n",
                input, expected, str);
        g_failures++;
    }
    free(str);
    expr_free(evaluated);
}
#define assert_eval_eq(in, exp, ff) check_eval_eq((in), (exp))

/* ---------------------------------------------------------------- */
/* Block A — the user's reported expression and its close relatives  */
/* ---------------------------------------------------------------- */

void test_user_report_q_sqrt2_sqrt3_cancellation(void) {
    /* The exact case from the user's bug report:
     *   e1 - e2 == Sqrt[2] - Sqrt[3] + 1/Sqrt[5 + 2 Sqrt[6]]
     * Mathematically zero: 1/Sqrt[5+2 Sqrt[6]] = Sqrt[3] - Sqrt[2]. */
    assert_eval_eq(
        "Simplify[Sqrt[2] - Sqrt[3] + 1/Sqrt[5 + 2 Sqrt[6]]]",
        "0", 0);
}

void test_user_report_recip_pair_to_2sqrt3(void) {
    /* 1/Sqrt[5+2 Sqrt[6]] + 1/Sqrt[5-2 Sqrt[6]]
     *   = (Sqrt[3]-Sqrt[2]) + (Sqrt[3]+Sqrt[2]) = 2 Sqrt[3]. */
    assert_eval_eq(
        "Simplify[1/Sqrt[5 + 2 Sqrt[6]] + 1/Sqrt[5 - 2 Sqrt[6]]]",
        "2 Sqrt[3]", 0);
}

void test_user_report_recip_pair_diff_to_neg2sqrt2(void) {
    /* 1/Sqrt[5+2 Sqrt[6]] - 1/Sqrt[5-2 Sqrt[6]]
     *   = (Sqrt[3]-Sqrt[2]) - (Sqrt[3]+Sqrt[2]) = -2 Sqrt[2]. */
    assert_eval_eq(
        "Simplify[1/Sqrt[5 + 2 Sqrt[6]] - 1/Sqrt[5 - 2 Sqrt[6]]]",
        "-2 Sqrt[2]", 0);
}

void test_forward_minus_recip_to_2sqrt2(void) {
    /* Sqrt[5+2 Sqrt[6]] - 1/Sqrt[5+2 Sqrt[6]]
     *   = (Sqrt[3]+Sqrt[2]) - (Sqrt[3]-Sqrt[2]) = 2 Sqrt[2]. */
    assert_eval_eq(
        "Simplify[Sqrt[5 + 2 Sqrt[6]] - 1/Sqrt[5 + 2 Sqrt[6]]]",
        "2 Sqrt[2]", 0);
}

void test_forward_plus_recip_to_2sqrt3(void) {
    /* Sqrt[5+2 Sqrt[6]] + 1/Sqrt[5+2 Sqrt[6]]
     *   = (Sqrt[3]+Sqrt[2]) + (Sqrt[3]-Sqrt[2]) = 2 Sqrt[3]. */
    assert_eval_eq(
        "Simplify[Sqrt[5 + 2 Sqrt[6]] + 1/Sqrt[5 + 2 Sqrt[6]]]",
        "2 Sqrt[3]", 0);
}

void test_recip_times_forward_identity(void) {
    /* Sqrt[x] * (1/Sqrt[x]) = 1 for any positive x; verify it holds when
     * the radicand is itself a denestable nested form. */
    assert_eval_eq(
        "Simplify[Sqrt[5 + 2 Sqrt[6]] * (1/Sqrt[5 + 2 Sqrt[6]])]",
        "1", 0);
}

/* ---------------------------------------------------------------- */
/* Block B — standalone denesting where leaf-count strictly reduces  */
/* ---------------------------------------------------------------- */
/* For these cases the rationalised reciprocal is shorter than the
 * input, so Simplify's leaf-count measure picks it. */

void test_standalone_recip_3_plus_2sqrt2(void) {
    /* 1/Sqrt[3+2 Sqrt[2]] = 1/(1+Sqrt[2]) = Sqrt[2]-1. */
    assert_eval_eq(
        "Simplify[1/Sqrt[3 + 2 Sqrt[2]]]",
        "-1 + Sqrt[2]", 0);
}

void test_standalone_recip_3_minus_2sqrt2(void) {
    /* 1/Sqrt[3-2 Sqrt[2]] = 1/(Sqrt[2]-1) = Sqrt[2]+1. */
    assert_eval_eq(
        "Simplify[1/Sqrt[3 - 2 Sqrt[2]]]",
        "1 + Sqrt[2]", 0);
}

void test_standalone_recip_5_minus_2sqrt6(void) {
    /* 1/Sqrt[5-2 Sqrt[6]] = 1/(Sqrt[3]-Sqrt[2]) = Sqrt[2]+Sqrt[3]. */
    assert_eval_eq(
        "Simplify[1/Sqrt[5 - 2 Sqrt[6]]]",
        "Sqrt[2] + Sqrt[3]", 0);
}

/* ---------------------------------------------------------------- */
/* Block C — non-denestable inputs must remain untouched             */
/* ---------------------------------------------------------------- */

void test_negative_discriminant_no_fire(void) {
    /* 1/Sqrt[1 + Sqrt[2]]: A=1, b=1, C=2, D = A^2 - b^2*C = -1.
     * Not a perfect square (the discriminant is negative), so the
     * denester cannot produce a real rationalised form. The input must
     * survive untouched. */
    assert_eval_eq(
        "Simplify[1/Sqrt[1 + Sqrt[2]]]",
        "1/Sqrt[1 + Sqrt[2]]", 0);
}

void test_symbolic_no_assumption_no_fire(void) {
    /* Without assumptions about x or y, we cannot prove the branch-
     * validity conditions (P, Q both nonneg), so 1/Sqrt[x + Sqrt[y]]
     * must survive untouched. The current convention prints it via
     * the Power[..., -1/2] form. */
    assert_eval_eq(
        "Simplify[1/Sqrt[x + Sqrt[y]]]",
        "1/Sqrt[x + Sqrt[y]]", 0);
}

void test_non_nested_recip_sqrt_untouched(void) {
    /* No inner Sqrt at all: the half-sum identity does not apply.
     * 1/Sqrt[2] should remain in its canonical Power form, not
     * accidentally enter the denester. */
    assert_eval_eq(
        "Simplify[1/Sqrt[2]]",
        "1/Sqrt[2]", 0);
}

/* ---------------------------------------------------------------- */
/* Block D — idempotency: a denested form must be a fixed point      */
/* ---------------------------------------------------------------- */

void test_idempotent_sqrt2_minus_1(void) {
    /* Simplify[Sqrt[2] - 1] should be a fixed point: nothing in the
     * denester should turn the rationalised form back into the nested
     * radical. */
    assert_eval_eq(
        "Simplify[Sqrt[2] - 1]",
        "-1 + Sqrt[2]", 0);
}

void test_idempotent_sqrt3_plus_sqrt2(void) {
    /* Same idempotency check at the multi-extension level. */
    assert_eval_eq(
        "Simplify[Sqrt[3] + Sqrt[2]]",
        "Sqrt[2] + Sqrt[3]", 0);
}

void test_idempotent_sqrt3_minus_sqrt2(void) {
    /* The asymmetric (subtraction) form: P > Q.  Mathilda canonical
     * Plus ordering keeps the positive Sqrt term first. */
    assert_eval_eq(
        "Simplify[Sqrt[3] - Sqrt[2]]",
        "Sqrt[3] - Sqrt[2]", 0);
}

/* ---------------------------------------------------------------- */
/* Block E — large discriminants: stress the conjugate rationaliser  */
/* ---------------------------------------------------------------- */
/* These radicands denest cleanly but the rationalised form has more
 * leaves than the input. The denester still produces a correct seed,
 * which we verify by checking cancellation against the expected form. */

void test_stress_7_plus_2sqrt10_cancel(void) {
    /* A=7, b=2, c=10, s=3, P=5, Q=2.
     * 1/Sqrt[7+2 Sqrt[10]] = (Sqrt[5]-Sqrt[2])/3. */
    assert_eval_eq(
        "Simplify[3/Sqrt[7 + 2 Sqrt[10]] - (Sqrt[5] - Sqrt[2])]",
        "0", 0);
}

void test_stress_11_plus_4sqrt6_cancel(void) {
    /* A=11, b=4, c=6, s=5, P=8, Q=3.
     * 1/Sqrt[11+4 Sqrt[6]] = (Sqrt[8]-Sqrt[3])/5 = (2 Sqrt[2]-Sqrt[3])/5. */
    assert_eval_eq(
        "Simplify[5/Sqrt[11 + 4 Sqrt[6]] - (2 Sqrt[2] - Sqrt[3])]",
        "0", 0);
}

void test_stress_9_plus_4sqrt5_cancel(void) {
    /* A=9, b=4, c=5, s=1, P=5, Q=4.
     * 1/Sqrt[9+4 Sqrt[5]] = Sqrt[5]-2. */
    assert_eval_eq(
        "Simplify[1/Sqrt[9 + 4 Sqrt[5]] - (Sqrt[5] - 2)]",
        "0", 0);
}

void test_stress_14_plus_6sqrt5_cancel(void) {
    /* A=14, b=6, c=5, D = 196 - 180 = 16, s=4, P=9, Q=5.
     * 1/Sqrt[14+6 Sqrt[5]] = (3-Sqrt[5])/4. */
    assert_eval_eq(
        "Simplify[4/Sqrt[14 + 6 Sqrt[5]] - (3 - Sqrt[5])]",
        "0", 0);
}

void test_stress_19_plus_6sqrt10_cancel(void) {
    /* A=19, b=6, c=10, D = 361 - 360 = 1, s=1, P=10, Q=9.
     * 1/Sqrt[19+6 Sqrt[10]] = 3 - Sqrt[10]. Wait: Sqrt[9]=3, Sqrt[10].
     * P > Q so result is Sqrt[10] - 3.  Hmm: 1/(Sqrt[10]+3) =
     * (Sqrt[10]-3)/(10-9) = Sqrt[10]-3.  Verify the sign. */
    assert_eval_eq(
        "Simplify[1/Sqrt[19 + 6 Sqrt[10]] - (Sqrt[10] - 3)]",
        "0", 0);
}

/* ---------------------------------------------------------------- */
/* Block F — Q(sqrt(2), sqrt(3)) multi-extension cancellations       */
/* ---------------------------------------------------------------- */
/* These exercise the chain:
 *   recip denester -> multi-extension qa_substrate cancellation
 * The original user case (test_user_report_q_sqrt2_sqrt3_cancellation)
 * lives here too; below are independent identities in the same
 * extension. */

void test_q_sqrt6_via_recip_2_plus_sqrt3(void) {
    /* 1/Sqrt[2+Sqrt[3]] = Sqrt[2]/Sqrt[2+Sqrt[3]] / Sqrt[2]
     * Direct identity: 1/Sqrt[2+Sqrt[3]] = (Sqrt[6]-Sqrt[2])/2.
     * (A=2, b=1, c=3, D=4-3=1, s=1, P=3/2, Q=1/2.
     *  Then Sqrt[3/2]-Sqrt[1/2] = (Sqrt[6]-Sqrt[2])/2.) */
    assert_eval_eq(
        "Simplify[1/Sqrt[2 + Sqrt[3]] - (Sqrt[6] - Sqrt[2])/2]",
        "0", 0);
}

void test_q_sqrt6_via_recip_2_minus_sqrt3(void) {
    /* 1/Sqrt[2-Sqrt[3]] = (Sqrt[6]+Sqrt[2])/2. */
    assert_eval_eq(
        "Simplify[1/Sqrt[2 - Sqrt[3]] - (Sqrt[6] + Sqrt[2])/2]",
        "0", 0);
}

void test_q_double_recip_sum_to_sqrt6(void) {
    /* 1/Sqrt[2+Sqrt[3]] + 1/Sqrt[2-Sqrt[3]] =
     * (Sqrt[6]-Sqrt[2])/2 + (Sqrt[6]+Sqrt[2])/2 = Sqrt[6]. */
    assert_eval_eq(
        "Simplify[1/Sqrt[2 + Sqrt[3]] + 1/Sqrt[2 - Sqrt[3]]]",
        "Sqrt[6]", 0);
}

/* ---------------------------------------------------------------- */
/* Block G — mixed Plus expressions that exercise the walker         */
/* ---------------------------------------------------------------- */

void test_walker_recip_inside_sum(void) {
    /* The denest must fire inside a Plus subtree, not just at the top
     * level.  This is the user's reported case, embedded in a larger
     * expression that includes a non-denestable term. */
    assert_eval_eq(
        "Simplify[x + Sqrt[2] - Sqrt[3] + 1/Sqrt[5 + 2 Sqrt[6]]]",
        "x", 0);
}

void test_walker_recip_inside_product(void) {
    /* The denest must fire inside a Times subtree.
     *   2 * (1/Sqrt[5+2 Sqrt[6]]) = 2 (Sqrt[3]-Sqrt[2]). */
    assert_eval_eq(
        "Simplify[2 (1/Sqrt[5 + 2 Sqrt[6]]) - (2 Sqrt[3] - 2 Sqrt[2])]",
        "0", 0);
}

void test_walker_recip_at_depth_two(void) {
    /* Walker reaches a denestable recip at depth 2 (under Plus then
     * under Times). */
    assert_eval_eq(
        "Simplify[1 + 3*(1/Sqrt[5 + 2 Sqrt[6]]) - 1 - (3 Sqrt[3] - 3 Sqrt[2])]",
        "0", 0);
}

void test_walker_pair_of_recip_in_difference(void) {
    /* Both recips at depth 2 (inside parenthesised sub-sums). We verify
     * the structural value rather than asserting cancellation against a
     * Q(sqrt(2), sqrt(3)) constant -- the latter trips a pre-existing
     * Simplify pathology unrelated to denesting.
     *
     *   1/Sqrt[5+2 Sqrt[6]] = Sqrt[3] - Sqrt[2]
     *   1/Sqrt[5-2 Sqrt[6]] = Sqrt[3] + Sqrt[2]
     *   sum                  = 2 Sqrt[3]. */
    assert_eval_eq(
        "Simplify[(1/Sqrt[5+2 Sqrt[6]]) + (1/Sqrt[5-2 Sqrt[6]])]",
        "2 Sqrt[3]", 0);
}

/* ---------------------------------------------------------------- */
/* Block H — coefficient combinations that exercise rationalisation  */
/* ---------------------------------------------------------------- */

void test_coeff_rational_recip(void) {
    /* (1/2) * (1/Sqrt[3+2 Sqrt[2]]) = (Sqrt[2]-1)/2. */
    assert_eval_eq(
        "Simplify[(1/2)*(1/Sqrt[3 + 2 Sqrt[2]]) - (Sqrt[2] - 1)/2]",
        "0", 0);
}

void test_coeff_negative_recip(void) {
    /* -1/Sqrt[3-2 Sqrt[2]] = -(1+Sqrt[2]). */
    assert_eval_eq(
        "Simplify[-1/Sqrt[3 - 2 Sqrt[2]] + (1 + Sqrt[2])]",
        "0", 0);
}

void test_coeff_integer_recip_then_cancel(void) {
    /* 5/Sqrt[11+4 Sqrt[6]] - (2 Sqrt[2] - Sqrt[3]) = 0
     * since 1/Sqrt[11+4 Sqrt[6]] = (2 Sqrt[2]-Sqrt[3])/5. */
    assert_eval_eq(
        "Simplify[5/Sqrt[11 + 4 Sqrt[6]] - 2 Sqrt[2] + Sqrt[3]]",
        "0", 0);
}

/* ---------------------------------------------------------------- */
/* Block I — paranoia: deep nesting, symbolic assumptions, ordering  */
/* ---------------------------------------------------------------- */

void test_paranoia_symbolic_recip_with_assumption(void) {
    /* Under x > y > 0, 1/Sqrt[x + Sqrt[x^2-y^2]] has a real rationalised
     * form: (Sqrt[(x+y)/2] - Sqrt[(x-y)/2]) / y. The cancellation form
     * uses the identity y * Sqrt[2] / Sqrt[x + Sqrt[x^2-y^2]] =
     * Sqrt[x+y] - Sqrt[x-y]. Verifying it by composition. */
    assert_eval_eq(
        "Assuming[x > y && y > 0, "
        "  Simplify[1/Sqrt[x + Sqrt[x^2 - y^2]] - "
        "          (Sqrt[(x+y)/2] - Sqrt[(x-y)/2])/y]]",
        "0", 0);
}

void test_paranoia_forward_and_recip_in_same_expr(void) {
    /* A forward and a recip denesting share the same radicand:
     *   Sqrt[5+2 Sqrt[6]] - 1/Sqrt[5+2 Sqrt[6]]
     *     = (Sqrt[3]+Sqrt[2]) - (Sqrt[3]-Sqrt[2]) = 2 Sqrt[2].
     * Tests that both paths cooperate inside one Simplify call. */
    assert_eval_eq(
        "Simplify[Sqrt[5+2 Sqrt[6]] - 1/Sqrt[5+2 Sqrt[6]]]",
        "2 Sqrt[2]", 0);
    /* Same expression, opposite sign:
     *   Sqrt[5+2 Sqrt[6]] + 1/Sqrt[5+2 Sqrt[6]]
     *     = (Sqrt[3]+Sqrt[2]) + (Sqrt[3]-Sqrt[2]) = 2 Sqrt[3]. */
    assert_eval_eq(
        "Simplify[Sqrt[5+2 Sqrt[6]] + 1/Sqrt[5+2 Sqrt[6]]]",
        "2 Sqrt[3]", 0);
}

void test_paranoia_ordering_does_not_break_match(void) {
    /* The walker should denest regardless of where the recip sits in
     * a sum. Same identity, args permuted. */
    assert_eval_eq(
        "Simplify[1/Sqrt[5 + 2 Sqrt[6]] + Sqrt[2] - Sqrt[3]]",
        "0", 0);
    assert_eval_eq(
        "Simplify[-Sqrt[3] + 1/Sqrt[5 + 2 Sqrt[6]] + Sqrt[2]]",
        "0", 0);
    assert_eval_eq(
        "Simplify[Sqrt[2] + 1/Sqrt[5 + 2 Sqrt[6]] - Sqrt[3]]",
        "0", 0);
}

/* Block J: half-integer powers beyond m = +/- 1. The walker handles
 * Power[Plus[..], m/2] for any odd m, so Sqrt[X]^k forms (which the
 * evaluator canonicalises to X^(k/2)) get denested too. For m > 0 we
 * return (Sqrt[P] + sign*Sqrt[Q])^m; for m < 0 we rationalise via
 * the conjugate and divide by s^|m|. */

void test_pow_cube_forward(void) {
    /* (5 + 2 Sqrt[6])^(3/2) - (Sqrt[3] + Sqrt[2])^3 == 0 */
    check_eval_eq(
        "Simplify[(5 + 2 Sqrt[6])^(3/2) - (Sqrt[3] + Sqrt[2])^3]",
        "0");
}

void test_pow_cube_recip(void) {
    /* User's cube case: 1/(5+2 Sqrt[6])^(3/2) - -(-Sqrt[3]+Sqrt[2])^3 = 0 */
    check_eval_eq(
        "Simplify[1/(5+2 Sqrt[6])^(3/2) - -(-Sqrt[3] + Sqrt[2])^3]",
        "0");
}

void test_pow_cube_recip_plus_form(void) {
    /* 1/(5+2 Sqrt[6])^(3/2) + (Sqrt[2] - Sqrt[3])^3 = 0
     * because (Sqrt[2]-Sqrt[3])^3 = -(Sqrt[3]-Sqrt[2])^3 = -1/(5+2Sqrt[6])^(3/2). */
    check_eval_eq(
        "Simplify[1/(5+2 Sqrt[6])^(3/2) + (Sqrt[2] - Sqrt[3])^3]",
        "0");
}

void test_pow_quintic_forward(void) {
    /* (5 + 2 Sqrt[6])^(5/2) - (Sqrt[3] + Sqrt[2])^5 == 0 */
    check_eval_eq(
        "Simplify[(5 + 2 Sqrt[6])^(5/2) - (Sqrt[3] + Sqrt[2])^5]",
        "0");
}

void test_pow_quintic_recip(void) {
    /* 1/(5 + 2 Sqrt[6])^(5/2) - (Sqrt[3] - Sqrt[2])^5 == 0 */
    check_eval_eq(
        "Simplify[1/(5 + 2 Sqrt[6])^(5/2) - (Sqrt[3] - Sqrt[2])^5]",
        "0");
}

void test_pow_cube_neg_b(void) {
    /* (3 - 2 Sqrt[2])^(3/2) - (Sqrt[2] - 1)^3 == 0
     * (3 - 2Sqrt[2]) = (Sqrt[2] - 1)^2 so ^(3/2) = |Sqrt[2]-1|^3 = (Sqrt[2]-1)^3
     * since Sqrt[2] > 1. */
    check_eval_eq(
        "Simplify[(3 - 2 Sqrt[2])^(3/2) - (Sqrt[2] - 1)^3]",
        "0");
}

void test_pow_cube_recip_nonunit_s(void) {
    /* P = 4, Q = 2, s = Sqrt[(7)^2 - 4*4*2] = Sqrt[49 - 32] is not clean,
     * so pick A=6, b=2, C=5 -> s^2 = 36 - 4*5 = 16, s = 4.
     * 6 + 2 Sqrt[5] = (Sqrt[5] + 1)^2 since (Sqrt[5]+1)^2 = 5 + 2Sqrt[5] + 1 = 6 + 2 Sqrt[5].
     * (6+2 Sqrt[5])^(3/2) = (Sqrt[5]+1)^3
     * 1/(6+2 Sqrt[5])^(3/2) = (Sqrt[5]-1)^3 / 4^3 = (Sqrt[5]-1)^3 / 64
     * So Simplify[1/(6+2 Sqrt[5])^(3/2) - (Sqrt[5]-1)^3/64] == 0. */
    check_eval_eq(
        "Simplify[1/(6+2 Sqrt[5])^(3/2) - (Sqrt[5]-1)^3/64]",
        "0");
}

/* ---------------------------------------------------------------- */

int main(void) {
    /* Line-buffered stdout so a hang inside any TEST() leaves a clear
     * trail of which test was running. */
    setvbuf(stdout, NULL, _IOLBF, 0);
    symtab_init();
    core_init();

    /* Block A: the user's reported cases. */
    TEST(test_user_report_q_sqrt2_sqrt3_cancellation);
    TEST(test_user_report_recip_pair_to_2sqrt3);
    TEST(test_user_report_recip_pair_diff_to_neg2sqrt2);
    TEST(test_forward_minus_recip_to_2sqrt2);
    TEST(test_forward_plus_recip_to_2sqrt3);
    TEST(test_recip_times_forward_identity);

    /* Block B: standalone denesting. */
    TEST(test_standalone_recip_3_plus_2sqrt2);
    TEST(test_standalone_recip_3_minus_2sqrt2);
    TEST(test_standalone_recip_5_minus_2sqrt6);

    /* Block C: non-denestable inputs remain untouched. */
    TEST(test_negative_discriminant_no_fire);
    TEST(test_symbolic_no_assumption_no_fire);
    TEST(test_non_nested_recip_sqrt_untouched);

    /* Block D: idempotency. */
    TEST(test_idempotent_sqrt2_minus_1);
    TEST(test_idempotent_sqrt3_plus_sqrt2);
    TEST(test_idempotent_sqrt3_minus_sqrt2);

    /* Block E: stress / larger discriminants. */
    TEST(test_stress_7_plus_2sqrt10_cancel);
    TEST(test_stress_11_plus_4sqrt6_cancel);
    TEST(test_stress_9_plus_4sqrt5_cancel);
    TEST(test_stress_14_plus_6sqrt5_cancel);
    TEST(test_stress_19_plus_6sqrt10_cancel);

    /* Block F: Q(sqrt(2),sqrt(3)) multi-extension. */
    TEST(test_q_sqrt6_via_recip_2_plus_sqrt3);
    TEST(test_q_sqrt6_via_recip_2_minus_sqrt3);
    TEST(test_q_double_recip_sum_to_sqrt6);

    /* Block G: walker placement. */
    TEST(test_walker_recip_inside_sum);
    TEST(test_walker_recip_inside_product);
    TEST(test_walker_recip_at_depth_two);
    TEST(test_walker_pair_of_recip_in_difference);

    /* Block H: coefficient combinations. */
    TEST(test_coeff_rational_recip);
    TEST(test_coeff_negative_recip);
    TEST(test_coeff_integer_recip_then_cancel);

    /* Block I: paranoia. */
    TEST(test_paranoia_symbolic_recip_with_assumption);
    TEST(test_paranoia_forward_and_recip_in_same_expr);
    TEST(test_paranoia_ordering_does_not_break_match);

    /* Block J: half-integer powers (Sqrt[X]^k for odd k). */
    TEST(test_pow_cube_forward);
    TEST(test_pow_cube_recip);
    TEST(test_pow_cube_recip_plus_form);
    TEST(test_pow_quintic_forward);
    TEST(test_pow_quintic_recip);
    TEST(test_pow_cube_neg_b);
    TEST(test_pow_cube_recip_nonunit_s);

    if (g_failures > 0) {
        fprintf(stderr, "\n%d simp_denest_inv_sqrt test(s) FAILED.\n",
                g_failures);
        return 1;
    }
    printf("All simp_denest_inv_sqrt tests passed!\n");
    return 0;
}
