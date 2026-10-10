#include "test_utils.h"
#include "symtab.h"
#include "core.h"

/*
 * Tests for the assumption-driven Log/Power identities applied by
 * Simplify. The identities form the strict-positive cascade:
 *
 *   Log[Times[u1, ..., un]] -> Plus[Log[u1], ..., Log[un]]
 *                                          when every ui is positive
 *   Log[Power[x, p]]       -> p Log[x]    when x positive and p real
 *   Power[Times[u1, ..., un], a]
 *                          -> Times[ui^a]  when every ui positive
 *   Power[Power[x, p], q]  -> Power[x, p q]
 *                                          when x positive and p real
 *
 * The general-real and general-complex cascade branches (with the
 * Boole / Floor / Ceiling phase corrections) are not implemented in v1
 * -- see Mathilda_spec.md for v2 scope. These tests cover only the
 * sound, strict-positive subset.
 *
 * The identities typically INCREASE leaf count (Log[a b] is 4 leaves,
 * Log[a] + Log[b] is 6), so they cannot win the standard complexity
 * tiebreak. Simplify therefore force-biases the logexp rewrite to win
 * whenever it changes the input under the supplied assumptions.
 */

/* ---- Log of products ---- */

void test_log_product_two_positive(void) {
    assert_eval_eq("Simplify[Log[a * b], a > 0 && b > 0]",
                   "Log[a] + Log[b]", 0);
}

void test_log_product_three_positive(void) {
    assert_eval_eq("Simplify[Log[a * b * c], a > 0 && b > 0 && c > 0]",
                   "Log[a] + Log[b] + Log[c]", 0);
}

void test_log_product_with_constant(void) {
    /* 5 is positive by numeric_sign, x is positive by assumption. */
    assert_eval_eq("Simplify[Log[5 * x], x > 0]",
                   "Log[5] + Log[x]", 0);
}

void test_log_product_one_positive_only(void) {
    /* b's sign is unknown -- the rewrite must NOT fire. */
    assert_eval_eq("Simplify[Log[a * b], a > 0]",
                   "Log[a b]", 0);
}

void test_log_product_no_assumption(void) {
    /* Without positivity, no expansion. */
    assert_eval_eq("Simplify[Log[a * b]]", "Log[a b]", 0);
}

/* ---- Log of quotients ---- */

void test_log_quotient_two_positive(void) {
    /* Log[a/b] is Log[Times[a, Power[b, -1]]]. The Times rule expands to
     * Log[a] + Log[Power[b, -1]], then the Power rule (b positive, -1
     * real) collapses Log[1/b] to -Log[b]. */
    assert_eval_eq("Simplify[Log[a / b], a > 0 && b > 0]",
                   "Log[a] - Log[b]", 0);
}

void test_log_inverse_positive(void) {
    assert_eval_eq("Simplify[Log[1 / a], a > 0]",
                   "-Log[a]", 0);
}

/* ---- Log of powers ---- */

void test_log_power_real_exponent(void) {
    /* Log[a^p] -> p Log[a] when a positive and p real. */
    assert_eval_eq("Simplify[Log[a^p], a > 0 && Element[p, Reals]]",
                   "Log[a] p", 0);
}

void test_log_power_integer_exponent_through_real(void) {
    /* Integer exponent is real by lattice. */
    assert_eval_eq("Simplify[Log[a^n], a > 0 && Element[n, Integers]]",
                   "Log[a] n", 0);
}

void test_log_power_no_real_assumption(void) {
    /* Without proving the exponent real, the rewrite must NOT fire. */
    assert_eval_eq("Simplify[Log[a^p], a > 0]",
                   "Log[a] p", 0);
}

void test_log_power_negative_base(void) {
    /* a < 0 does not prove positive, so the rewrite must NOT fire. */
    assert_eval_eq("Simplify[Log[a^p], a < 0 && Element[p, Reals]]",
                   "Log[a^p]", 0);
}

/* ---- Power of products ---- */

void test_power_of_product_two_positive(void) {
    assert_eval_eq("Simplify[(a * b)^c, a > 0 && b > 0]",
                   "a^c b^c", 0);
}

void test_power_of_product_three_positive(void) {
    assert_eval_eq("Simplify[(a * b * d)^c, a > 0 && b > 0 && d > 0]",
                   "a^c b^c d^c", 0);
}

void test_sqrt_of_product_positive(void) {
    /* Sqrt[a b] = (a b)^(1/2). With a, b > 0 the Power rule distributes
     * over Times, giving Sqrt[a] Sqrt[b]. */
    assert_eval_eq("Simplify[Sqrt[a * b], a > 0 && b > 0]",
                   "Sqrt[a] Sqrt[b]", 0);
}

void test_power_of_quotient_positive(void) {
    /* (a/b)^c with both positive expands via the Times rule
     * (a/b is Times[a, b^-1]). */
    assert_eval_eq("Simplify[(a / b)^c, a > 0 && b > 0]",
                   "a^c b^(-c)", 0);
}

/* ---- Tower of powers ---- */

void test_power_of_power_real_inner(void) {
    /* (a^p)^q -> a^(p q) when a positive and p real. */
    assert_eval_eq("Simplify[(a^p)^q, a > 0 && Element[p, Reals]]",
                   "a^(p q)", 0);
}

void test_power_of_power_integer_inner(void) {
    /* Integer inner exponent is real, so the rule fires. */
    assert_eval_eq("Simplify[(a^2)^q, a > 0]",
                   "a^(2 q)", 0);
}

void test_power_of_power_no_real_inner(void) {
    /* Inner exponent's signedness/realness unknown -- rule must not fire.
     * The result is the unchanged input Power[Power[a, p], q], which the
     * printer renders with disambiguating parentheses as "(a^p)^q" (NOT
     * a^(p^q)).  The identity (a^p)^q -> a^(p q) needs p real, so under
     * a > 0 alone WL also leaves this unchanged. */
    assert_eval_eq("Simplify[(a^p)^q, a > 0]",
                   "(a^p)^q", 0);
}

/* ---- Composite cancellation cases enabled by the cascade ---- */

void test_log_pow_plus_log_inv(void) {
    /* Log[x^2] + Log[1/x^2] under x > 0 should fully cancel to 0:
     *   2 Log[x] + (-2 Log[x]) = 0. */
    assert_eval_eq("Simplify[Log[x^2] + Log[1/x^2], x > 0]", "0", 0);
}

void test_log_difference_of_constant_factor(void) {
    /* Log[2 x] - Log[x] -> Log[2] under x > 0. The cascade expands
     * Log[2 x] to Log[2] + Log[x], then -Log[x] cancels. */
    assert_eval_eq("Simplify[Log[2 x] - Log[x], x > 0]", "Log[2]", 0);
}

void test_pow_distrib_collapses(void) {
    /* (a/b)^c * b^c -> a^c after distribution. */
    assert_eval_eq("Simplify[(a/b)^c * b^c, a > 0 && b > 0]", "a^c", 0);
}

/* ---- Log-power symmetry: x^Log[y] = y^Log[x] (PowBaseToExp) ---- */

void test_logpow_symmetry_both_positive(void) {
    /* Both bases positive: base^exp -> Exp[exp Log[base]] exposes the
     * symmetry (both become Exp[Log[x] Log[y]]) and the terms cancel. */
    assert_eval_eq("Simplify[x^Log[y] - y^Log[x], x > 0 && y > 0]", "0", 0);
}

void test_logpow_symmetry_naming(void) {
    /* Firing must not depend on variable naming/order. */
    assert_eval_eq("Simplify[a^Log[b] - b^Log[a], a > 0 && b > 0]", "0", 0);
}

void test_logpow_symmetry_no_assumption(void) {
    /* Without positivity the identity fails (x <= 0), so it must NOT fire. */
    assert_eval_eq("Simplify[x^Log[y] - y^Log[x]]",
                   "x^Log[y] - y^Log[x]", 0);
}

void test_logpow_symmetry_one_sided_y(void) {
    /* Only y > 0: the two-sided gate still needs x > 0 (the reverted base x
     * would otherwise be asserted positive, failing at x = 0 / x < 0). */
    assert_eval_eq("Simplify[x^Log[y] - y^Log[x], y > 0]",
                   "x^Log[y] - y^Log[x]", 0);
}

void test_logpow_symmetry_one_sided_x(void) {
    assert_eval_eq("Simplify[x^Log[y] - y^Log[x], x > 0]",
                   "x^Log[y] - y^Log[x]", 0);
}

void test_logpow_idempotent_standalone(void) {
    /* A standalone x^Log[x] must not be force-rewritten into the larger
     * Exp[Log[x]^2] form: the strict-score seed drops the worse candidate. */
    assert_eval_eq("Simplify[x^Log[x], x > 0]", "x^Log[x]", 0);
}

/* ---- Log[Exp[...]] inverse pair ---- */

void test_log_exp_positive(void) {
    /* Log[Exp[x]] for any x: works through the existing TrigToExp /
     * ExpToTrig roundtrip rather than the new cascade, but listed here
     * as a sanity check that the cascade does not break it. */
    assert_eval_eq("Simplify[Log[Exp[x]], x > 0]", "x", 0);
}

/* ---- Assumption-gated log/exp/trig identities (v0.331 campaign) ---- */

/* T1: range-gated inverse-of-direct trig collapse. */
void test_arcsin_sin_collapse(void) {
    assert_eval_eq("Simplify[ArcSin[Sin[x]] - x, -Pi/2 <= x <= Pi/2]", "0", 0);
}
void test_arcsin_sin_no_assumption(void) {
    /* Soundness: ArcSin[Sin[x]] != x off [-Pi/2, Pi/2]. */
    assert_eval_eq("Simplify[ArcSin[Sin[x]] - x]", "-x + ArcSin[Sin[x]]", 0);
}
void test_arccos_cos_collapse(void) {
    assert_eval_eq("Simplify[ArcCos[Cos[x]] - x, 0 <= x <= Pi]", "0", 0);
}
void test_arccos_cos_no_assumption(void) {
    assert_eval_eq("Simplify[ArcCos[Cos[x]] - x]", "-x + ArcCos[Cos[x]]", 0);
}
void test_arccot_cot_collapse(void) {
    assert_eval_eq("Simplify[ArcCot[Cot[x]] - x, 0 < x < Pi]", "0", 0);
}

/* T3a: ArcCos logarithmic form (unconditional, principal-value general). */
void test_arccos_log_form_assumed(void) {
    assert_eval_eq("Simplify[ArcCos[x] + I Log[x + I Sqrt[1 - x^2]], -1 <= x <= 1]", "0", 0);
}
void test_arccos_log_form_unconditional(void) {
    /* The identity is principal-value general, so it reduces with no assumption. */
    assert_eval_eq("Simplify[Log[x + I Sqrt[1 - x^2]] - I ArcCos[x]]", "0", 0);
}

/* T2: ArcCosh logarithmic form, gated x >= 1 (combined radical). */
void test_arccosh_log_form(void) {
    assert_eval_eq("Simplify[ArcCosh[x] - Log[x + Sqrt[x^2 - 1]], x >= 1]", "0", 0);
}
void test_arccosh_log_no_assumption(void) {
    /* Soundness: Log[x+Sqrt[x^2-1]] != ArcCosh[x] off [1, inf). */
    assert_eval_eq("Simplify[ArcCosh[x] - Log[x + Sqrt[x^2 - 1]]]",
                   "ArcCosh[x] - Log[x + Sqrt[-1 + x^2]]", 0);
}
void test_arccosh_log_wrong_region(void) {
    /* Soundness: residual is nonzero at x <= -1, so it must NOT reduce. */
    assert_eval_eq("Simplify[ArcCosh[x] - Log[x + Sqrt[x^2 - 1]], x <= -1]",
                   "ArcCosh[x] - Log[x + Sqrt[-1 + x^2]]", 0);
}

/* T4: Exp additive periodicity, gated n integer. */
void test_exp_periodicity_sum(void) {
    assert_eval_eq("Simplify[Exp[x + 2 I Pi n] - Exp[x], Element[n, Integers]]", "0", 0);
}
void test_exp_periodicity_no_integer(void) {
    assert_eval_eq("Simplify[Exp[x + 2 I Pi n] - Exp[x]]",
                   "-E^x + E^((2*I) Pi n + x)", 0);
}

/* T5: (-1)^(k+n) integer-exponent split; and the full Cos[nPi]-Exp[I n Pi]. */
void test_neg_one_power_split(void) {
    assert_eval_eq("Simplify[(-1)^(n + 1) + (-1)^n, Element[n, Integers]]", "0", 0);
}
void test_cos_minus_exp_integer(void) {
    assert_eval_eq("Simplify[Cos[n Pi] - Exp[I n Pi], Element[n, Integers]]", "0", 0);
}

/* T6: branch-gated Sqrt[E^w] -> E^(w/2) on the principal strip. */
void test_sqrt_exp_strip(void) {
    assert_eval_eq("Simplify[Sqrt[Exp[2 I x]] - Exp[I x], -Pi/2 < x < Pi/2]", "0", 0);
}
void test_sqrt_exp_no_assumption(void) {
    assert_eval_eq("Simplify[Sqrt[Exp[2 I x]] - Exp[I x]]",
                   "-E^(I x) + Sqrt[E^((2*I) x)]", 0);
}
void test_sqrt_exp_outside_strip(void) {
    /* Soundness: Sqrt[E^(2 I x)] = -E^(I x) off the strip, so no reduction. */
    assert_eval_eq("Simplify[Sqrt[Exp[2 I x]] - Exp[I x], Pi < x < 2 Pi]",
                   "-E^(I x) + Sqrt[E^((2*I) x)]", 0);
}

/* ---- v0.332: deep sign oracle for Abs / Sqrt[_^2] under interval
 * assumptions (Reduce/CAD bridge + pole-bearing trig decomposition) ---- */
void test_abs_sin_quadrant1(void) {
    assert_eval_eq("Simplify[Abs[Sin[x]] - Sin[x], 0 < x < Pi/2]", "0", 0);
}
void test_abs_cos_quadrant2(void) {
    assert_eval_eq("Simplify[Abs[Cos[x]] + Cos[x], Pi/2 < x < Pi]", "0", 0);
}
void test_abs_tan_quadrant2(void) {
    /* Tan is pole-bearing: resolved via sign(Tan)=sign(Sin)*sign(Cos). */
    assert_eval_eq("Simplify[Abs[Tan[x]] + Tan[x], Pi/2 < x < Pi]", "0", 0);
}
void test_abs_sec_quadrant1(void) {
    assert_eval_eq("Simplify[Abs[Sec[x]] - Sec[x], 0 < x < Pi/2]", "0", 0);
}
void test_abs_sum_sin_cos_quadrant1(void) {
    assert_eval_eq("Simplify[Abs[Sin[x] + Cos[x]] - (Sin[x] + Cos[x]), 0 < x < Pi/2]", "0", 0);
}
void test_abs_poly_shifted_interval(void) {
    /* x-2 < 0 on (0,1); the cheap provers can't shift the bound, Reduce can. */
    assert_eval_eq("Simplify[Abs[x - 2] - (2 - x), 0 < x < 1]", "0", 0);
}
void test_sqrt_cos_squared_quadrant3(void) {
    assert_eval_eq("Simplify[Sqrt[Cos[x]^2] + Cos[x], Pi < x < 3 Pi/2]", "0", 0);
}
/* Soundness: the oracle must DECLINE where the sign is not fixed. */
void test_abs_sin_no_assumption_declines(void) {
    assert_eval_eq("Simplify[Abs[Sin[x]]]", "Abs[Sin[x]]", 0);
}
void test_abs_cos_sign_change_declines(void) {
    /* Cos flips sign at Pi/2 inside (0,Pi) -> no collapse. */
    assert_eval_eq("Simplify[Abs[Cos[x]], 0 < x < Pi]", "Abs[Cos[x]]", 0);
}
void test_abs_sin_half_bounded_declines(void) {
    /* x>0 is unbounded above: Sin sign not fixed -> no collapse. */
    assert_eval_eq("Simplify[Abs[Sin[x]], x > 0]", "Abs[Sin[x]]", 0);
}

/* ---- v0.333: generalized Sqrt[c f^2 ...] = Sqrt[c] Abs[f] extraction ---- */
void test_sqrt_coeff_square_trig(void) {
    /* Sqrt[2 Sin[x]^2] -> Sqrt[2] Sin[x] under 0<x<Pi/2. */
    assert_eval_eq("Simplify[Sqrt[1 - Cos[2 x]] - Sqrt[2] Sin[x], 0 < x < Pi/2]", "0", 0);
}
void test_sqrt_coeff_square_trig_quadrant2(void) {
    assert_eval_eq("Simplify[Sqrt[1 - Cos[2 x]] - Sqrt[2] Sin[x], Pi/2 < x < Pi]", "0", 0);
}
void test_sqrt_product_of_squares(void) {
    /* Sqrt[Cos[x]^2 Sin[x]^2] -> Cos[x] Sin[x] on quadrant I. */
    assert_eval_eq("Simplify[Sqrt[Cos[x]^2 Sin[x]^2] - Sin[x] Cos[x], 0 < x < Pi/2]", "0", 0);
}
void test_sqrt_coeff_square_composite_base(void) {
    /* Deep sign on a composite base: x+1 > 0 under x>0. */
    assert_eval_eq("Simplify[Sqrt[2 (x + 1)^2] - Sqrt[2] (x + 1), x > 0]", "0", 0);
}
void test_sqrt_coeff_square_real_abs(void) {
    /* Sign undetermined but real -> Sqrt[2] Abs[Sin[x]]. */
    assert_eval_eq("Simplify[Sqrt[2 Sin[x]^2], Element[x, Reals]]", "Sqrt[2] Abs[Sin[x]]", 0);
}
void test_sqrt_coeff_square_no_assumption_declines(void) {
    /* Sin[x] not provably real -> no extraction. */
    assert_eval_eq("Simplify[Sqrt[2 Sin[x]^2]]", "Sqrt[2 Sin[x]^2]", 0);
}
void test_sqrt_negative_coeff_square(void) {
    /* c<0: Sqrt[-2 x^2] = I Sqrt[2] Abs[x] (split justified by x^2>=0). */
    assert_eval_eq("Simplify[Sqrt[-2 x^2], Element[x, Reals]]", "I Sqrt[2] Abs[x]", 0);
}

/* ---- v0.334: symbolic-integer Pi periodicity Sin[t + k Pi] = (-1)^k Sin[t] ---- */
void test_pi_period_sin_even(void) {
    assert_eval_eq("Simplify[Sin[x + 2 n Pi] - Sin[x], Element[n, Integers]]", "0", 0);
}
void test_pi_period_cos_even(void) {
    assert_eval_eq("Simplify[Cos[x + 2 n Pi] - Cos[x], Element[n, Integers]]", "0", 0);
}
void test_pi_period_sin_odd(void) {
    assert_eval_eq("Simplify[Sin[x + (2 n + 1) Pi] + Sin[x], Element[n, Integers]]", "0", 0);
}
void test_pi_period_cos_odd(void) {
    assert_eval_eq("Simplify[Cos[x + (2 n + 1) Pi] + Cos[x], Element[n, Integers]]", "0", 0);
}
void test_pi_period_sin_4n(void) {
    assert_eval_eq("Simplify[Sin[x + 4 n Pi] - Sin[x], Element[n, Integers]]", "0", 0);
}
void test_pi_period_tan_any(void) {
    /* Tan has period Pi: Tan[x + 3 n Pi] = Tan[x] regardless of parity. */
    assert_eval_eq("Simplify[Tan[x + 3 n Pi] - Tan[x], Element[n, Integers]]", "0", 0);
}
void test_pi_period_no_integer_declines(void) {
    /* Without n integer, 2 n Pi is not an integer multiple of Pi. */
    assert_eval_eq("Simplify[Sin[x + 2 n Pi]]", "Sin[2 Pi n + x]", 0);
}
void test_pi_period_half_multiple_declines(void) {
    /* n Pi/2 is not an integer multiple of Pi. */
    assert_eval_eq("Simplify[Sin[x + n Pi/2], Element[n, Integers]]", "Sin[1/2 Pi n + x]", 0);
}

/* ---- v0.335: Sqrt-local radicand preparation (double/half-angle) ---- */
void test_sqrt_double_angle_plus_q1(void) {
    assert_eval_eq("Simplify[Sqrt[1 + Cos[2 x]] - Sqrt[2] Cos[x], 0 < x < Pi/2]", "0", 0);
}
void test_sqrt_double_angle_plus_q2(void) {
    /* Cos < 0 on (Pi/2, Pi): Sqrt[1+Cos[2x]] = Sqrt[2](-Cos[x]). */
    assert_eval_eq("Simplify[Sqrt[1 + Cos[2 x]] + Sqrt[2] Cos[x], Pi/2 < x < Pi]", "0", 0);
}
void test_sqrt_half_angle_sin(void) {
    assert_eval_eq("Simplify[Sqrt[(1 - Cos[x])/2] - Sin[x/2], 0 < x < Pi/2]", "0", 0);
}
void test_sqrt_half_angle_cos(void) {
    assert_eval_eq("Simplify[Sqrt[(1 + Cos[x])/2] - Cos[x/2], 0 < x < Pi/2]", "0", 0);
}
void test_sqrt_half_angle_wide_interval(void) {
    /* x/2 in (0, Pi/2) over 0<x<Pi, so Sin[x/2] > 0. */
    assert_eval_eq("Simplify[Sqrt[(1 - Cos[x])/2] - Sin[x/2], 0 < x < Pi]", "0", 0);
}
void test_sqrt_double_angle_real_abs(void) {
    assert_eval_eq("Simplify[Sqrt[1 + Cos[2 x]], Element[x, Reals]]", "Sqrt[2] Abs[Cos[x]]", 0);
}
void test_sqrt_double_angle_no_assumption_declines(void) {
    assert_eval_eq("Simplify[Sqrt[1 + Cos[2 x]]]", "Sqrt[1 + Cos[2 x]]", 0);
}
void test_sqrt_prepare_longer_plus_declines(void) {
    /* 3-term Plus: the 2-element rule patterns do not match. */
    assert_eval_eq("Simplify[Sqrt[1 + Cos[2 x] + y], Element[x, Reals] && Element[y, Reals]]",
                   "Sqrt[1 + Cos[2 x] + y]", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_log_product_two_positive);
    TEST(test_log_product_three_positive);
    TEST(test_log_product_with_constant);
    TEST(test_log_product_one_positive_only);
    TEST(test_log_product_no_assumption);

    TEST(test_log_quotient_two_positive);
    TEST(test_log_inverse_positive);

    TEST(test_log_power_real_exponent);
    TEST(test_log_power_integer_exponent_through_real);
    TEST(test_log_power_no_real_assumption);
    TEST(test_log_power_negative_base);

    TEST(test_power_of_product_two_positive);
    TEST(test_power_of_product_three_positive);
    TEST(test_sqrt_of_product_positive);
    TEST(test_power_of_quotient_positive);

    TEST(test_power_of_power_real_inner);
    TEST(test_power_of_power_integer_inner);
    TEST(test_power_of_power_no_real_inner);

    TEST(test_log_pow_plus_log_inv);
    TEST(test_log_difference_of_constant_factor);
    TEST(test_pow_distrib_collapses);

    TEST(test_logpow_symmetry_both_positive);
    TEST(test_logpow_symmetry_naming);
    TEST(test_logpow_symmetry_no_assumption);
    TEST(test_logpow_symmetry_one_sided_y);
    TEST(test_logpow_symmetry_one_sided_x);
    TEST(test_logpow_idempotent_standalone);

    TEST(test_log_exp_positive);

    /* v0.331 assumption-gated campaign */
    TEST(test_arcsin_sin_collapse);
    TEST(test_arcsin_sin_no_assumption);
    TEST(test_arccos_cos_collapse);
    TEST(test_arccos_cos_no_assumption);
    TEST(test_arccot_cot_collapse);
    TEST(test_arccos_log_form_assumed);
    TEST(test_arccos_log_form_unconditional);
    TEST(test_arccosh_log_form);
    TEST(test_arccosh_log_no_assumption);
    TEST(test_arccosh_log_wrong_region);
    TEST(test_exp_periodicity_sum);
    TEST(test_exp_periodicity_no_integer);
    TEST(test_neg_one_power_split);
    TEST(test_cos_minus_exp_integer);
    TEST(test_sqrt_exp_strip);
    TEST(test_sqrt_exp_no_assumption);
    TEST(test_sqrt_exp_outside_strip);

    /* v0.332 deep sign oracle */
    TEST(test_abs_sin_quadrant1);
    TEST(test_abs_cos_quadrant2);
    TEST(test_abs_tan_quadrant2);
    TEST(test_abs_sec_quadrant1);
    TEST(test_abs_sum_sin_cos_quadrant1);
    TEST(test_abs_poly_shifted_interval);
    TEST(test_sqrt_cos_squared_quadrant3);
    TEST(test_abs_sin_no_assumption_declines);
    TEST(test_abs_cos_sign_change_declines);
    TEST(test_abs_sin_half_bounded_declines);

    /* v0.333 generalized Sqrt[c f^2] extraction */
    TEST(test_sqrt_coeff_square_trig);
    TEST(test_sqrt_coeff_square_trig_quadrant2);
    TEST(test_sqrt_product_of_squares);
    TEST(test_sqrt_coeff_square_composite_base);
    TEST(test_sqrt_coeff_square_real_abs);
    TEST(test_sqrt_coeff_square_no_assumption_declines);
    TEST(test_sqrt_negative_coeff_square);

    /* v0.334 symbolic-integer Pi periodicity */
    TEST(test_pi_period_sin_even);
    TEST(test_pi_period_cos_even);
    TEST(test_pi_period_sin_odd);
    TEST(test_pi_period_cos_odd);
    TEST(test_pi_period_sin_4n);
    TEST(test_pi_period_tan_any);
    TEST(test_pi_period_no_integer_declines);
    TEST(test_pi_period_half_multiple_declines);

    /* v0.335 Sqrt-local radicand preparation */
    TEST(test_sqrt_double_angle_plus_q1);
    TEST(test_sqrt_double_angle_plus_q2);
    TEST(test_sqrt_half_angle_sin);
    TEST(test_sqrt_half_angle_cos);
    TEST(test_sqrt_half_angle_wide_interval);
    TEST(test_sqrt_double_angle_real_abs);
    TEST(test_sqrt_double_angle_no_assumption_declines);
    TEST(test_sqrt_prepare_longer_plus_declines);

    printf("All logexp Simplify tests passed!\n");
    return 0;
}
