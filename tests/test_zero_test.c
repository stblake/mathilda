/*
 * test_zero_test.c — unit + stress tests for PossibleZeroQ.
 *
 * Test groups
 *   1.  Stage 0  — structural shortcuts (literal zero, Pi, Complex, …)
 *   2.  Stage 1  — polynomial and rational identities over Q[x_1,…,x_n]
 *   3.  Stage 2  — closed-form numeric expressions
 *                  (including catastrophic-cancellation regression)
 *   4.  Stage 3  — Schwartz–Zippel identity testing for transcendental
 *                  expressions with free symbols
 *   5.  Listable threading and registration smoke tests
 *   6.  Stress tests: large polynomial identities, multivariate
 *       Vieta-style, trig/log identities, complex-number identities,
 *       memory-leak smoke load
 *
 * All Stage-3 tests pre-seed the global PRNG via SeedRandom[42] so test
 * outcomes are reproducible.
 */

#include "expr.h"
#include "eval.h"
#include "core.h"
#include "parse.h"
#include "symtab.h"
#include "test_utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Evaluate a Mathilda input string and assert that the result prints as
 * the given expected symbol ("True" / "False" / "{True, False, ...}"). */
static void assert_pzq(const char* input, const char* expected) {
    Expr* parsed = parse_expression(input);
    ASSERT(parsed != NULL);
    Expr* evald  = evaluate(parsed);
    expr_free(parsed);
    char* s = expr_to_string(evald);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  expected: %s\n  actual:   %s\n",
                input, expected, s);
        free(s);
        expr_free(evald);
        exit(1);
    }
    free(s);
    expr_free(evald);
}

/* Convenience: re-seed the PRNG to a known value for reproducibility. */
static void seed_rng(int64_t seed) {
    char buf[64];
    snprintf(buf, sizeof buf, "SeedRandom[%lld]", (long long)seed);
    Expr* p = parse_expression(buf);
    Expr* r = evaluate(p);
    expr_free(p);
    expr_free(r);
}

/* ============================================================== */
/*  1. Stage 0 — structural shortcuts                             */
/* ============================================================== */

static void test_stage0_integer_zero(void) {
    assert_pzq("PossibleZeroQ[0]", "True");
}
static void test_stage0_real_zero(void) {
    assert_pzq("PossibleZeroQ[0.0]", "True");
}
static void test_stage0_real_negzero(void) {
    assert_pzq("PossibleZeroQ[-0.0]", "True");
}
static void test_stage0_integer_one(void) {
    assert_pzq("PossibleZeroQ[1]", "False");
}
static void test_stage0_integer_negone(void) {
    assert_pzq("PossibleZeroQ[-1]", "False");
}
static void test_stage0_bigint(void) {
    assert_pzq("PossibleZeroQ[10^100]", "False");
}
static void test_stage0_bigint_negative(void) {
    assert_pzq("PossibleZeroQ[-2^200]", "False");
}
static void test_stage0_complex_zero(void) {
    assert_pzq("PossibleZeroQ[Complex[0, 0]]", "True");
}
static void test_stage0_complex_pure_imag(void) {
    assert_pzq("PossibleZeroQ[Complex[0, 1]]", "False");
}
static void test_stage0_pi(void) {
    assert_pzq("PossibleZeroQ[Pi]", "False");
}
static void test_stage0_e(void) {
    assert_pzq("PossibleZeroQ[E]", "False");
}
static void test_stage0_eulergamma(void) {
    assert_pzq("PossibleZeroQ[EulerGamma]", "False");
}
static void test_stage0_free_symbol(void) {
    assert_pzq("PossibleZeroQ[someUnboundSymbol]", "False");
}
static void test_stage0_rational_zero(void) {
    assert_pzq("PossibleZeroQ[Rational[0, 5]]", "True");
}
static void test_stage0_rational_nonzero(void) {
    assert_pzq("PossibleZeroQ[1/7]", "False");
}

/* ============================================================== */
/*  2. Stage 1 — rational identities                              */
/* ============================================================== */

static void test_stage1_diff_of_squares(void) {
    assert_pzq("PossibleZeroQ[(x + 1)(x - 1) - x^2 + 1]", "True");
}
static void test_stage1_cube_expansion(void) {
    assert_pzq("PossibleZeroQ[(x + y)^3 - x^3 - 3 x^2 y - 3 x y^2 - y^3]",
               "True");
}
static void test_stage1_common_denominator(void) {
    assert_pzq("PossibleZeroQ[1/x + 1/y - (x + y)/(x y)]", "True");
}
static void test_stage1_cancel_linear(void) {
    /* (x^2 - 1)/(x - 1) is x + 1 after Cancel. */
    assert_pzq("PossibleZeroQ[(x^2 - 1)/(x - 1) - (x + 1)]", "True");
}
static void test_stage1_perfect_square(void) {
    assert_pzq("PossibleZeroQ[x^2 - 2 x + 1 - (x - 1)^2]", "True");
}
static void test_stage1_obvious_nonzero(void) {
    assert_pzq("PossibleZeroQ[x + y]", "False");
}
static void test_stage1_nonzero_constant(void) {
    assert_pzq("PossibleZeroQ[x^2 + 1]", "False");
}
static void test_stage1_triple_product(void) {
    assert_pzq(
        "PossibleZeroQ[(x - a)(x - b)(x - c) - (x^3 - (a + b + c) x^2 + "
        "(a b + a c + b c) x - a b c)]",
        "True");
}
static void test_stage1_partial_fraction(void) {
    /* 1/((x-1)(x-2)) == 1/(x-2) - 1/(x-1) */
    assert_pzq(
        "PossibleZeroQ[1/((x - 1)(x - 2)) - (1/(x - 2) - 1/(x - 1))]",
        "True");
}
static void test_stage1_three_var_rational(void) {
    assert_pzq(
        "PossibleZeroQ[a/(b c) + b/(a c) + c/(a b) - "
        "(a^2 + b^2 + c^2)/(a b c)]",
        "True");
}

/* ============================================================== */
/*  3. Stage 2 — closed-form numeric                              */
/* ============================================================== */

static void test_stage2_exp_ipi_quarter(void) {
    /* E^(I Pi/4) - (-1)^(1/4) = 0 */
    assert_pzq("PossibleZeroQ[E^(I Pi/4) - (-1)^(1/4)]", "True");
}
static void test_stage2_binomial_e_pi(void) {
    assert_pzq(
        "PossibleZeroQ[(E + Pi)^2 - E^2 - Pi^2 - 2 E Pi]", "True");
}
static void test_stage2_e_pi_vs_pi_e(void) {
    /* E^Pi ≈ 23.14, Pi^E ≈ 22.46.  Definitely non-zero. */
    assert_pzq("PossibleZeroQ[E^Pi - Pi^E]", "False");
}
static void test_stage2_sin_pi(void) {
    assert_pzq("PossibleZeroQ[Sin[Pi]]", "True");
}
static void test_stage2_cos_pi_half(void) {
    assert_pzq("PossibleZeroQ[Cos[Pi/2]]", "True");
}
static void test_stage2_sqrt_two_squared(void) {
    assert_pzq("PossibleZeroQ[Sqrt[2]^2 - 2]", "True");
}
static void test_stage2_pyth_at_one(void) {
    assert_pzq("PossibleZeroQ[Sin[1]^2 + Cos[1]^2 - 1]", "True");
}
static void test_stage2_double_angle_sin_at_one(void) {
    assert_pzq("PossibleZeroQ[Sin[2] - 2 Sin[1] Cos[1]]", "True");
}
static void test_stage2_complex_two_2i(void) {
    /* 2^(2I) - 2^(-2I) - 2 I Sin[Log[4]] = 0 (Euler/De Moivre). */
    assert_pzq(
        "PossibleZeroQ[2^(2 I) - 2^(-2 I) - 2 I Sin[Log[4]]]", "True");
}
static void test_stage2_log_product(void) {
    assert_pzq("PossibleZeroQ[Log[2] + Log[3] - Log[6]]", "True");
}
static void test_stage2_log_quotient(void) {
    assert_pzq("PossibleZeroQ[Log[10] - Log[2] - Log[5]]", "True");
}
static void test_stage2_tiny_nonzero(void) {
    assert_pzq("PossibleZeroQ[10^(-30)]", "False");
}
static void test_stage2_sqrt_irrational_close(void) {
    assert_pzq("PossibleZeroQ[Sqrt[2] - 1.41421356]", "False");
}
static void test_stage2_catastrophic_cancellation(void) {
    /* 10^20 + 1 - 10^20 = 1 in exact arithmetic, but at machine
     * precision the 1 is lost.  The precision-bumping ladder must
     * detect this and answer False. */
    assert_pzq("PossibleZeroQ[10^20 + 1 - 10^20]", "False");
}
static void test_stage2_complex_conjugate(void) {
    /* (1 + I)(1 - I) = 2, not 0. */
    assert_pzq("PossibleZeroQ[(1 + I)(1 - I) - 2]", "True");
}
static void test_stage2_complex_imag_squared(void) {
    assert_pzq("PossibleZeroQ[I^2 + 1]", "True");
}

/* ============================================================== */
/*  4. Stage 3 — Schwartz–Zippel symbolic non-rational            */
/* ============================================================== */

static void test_stage3_pythagorean_identity(void) {
    seed_rng(42);
    assert_pzq("PossibleZeroQ[Sin[x]^2 + Cos[x]^2 - 1]", "True");
}
static void test_stage3_double_angle_sin(void) {
    seed_rng(42);
    assert_pzq("PossibleZeroQ[Sin[2 x] - 2 Sin[x] Cos[x]]", "True");
}
static void test_stage3_sqrt_x_squared_minus_x(void) {
    seed_rng(42);
    /* Sqrt[x^2] != x in general (false at x = -2, complex x, …). */
    assert_pzq("PossibleZeroQ[Sqrt[x^2] - x]", "False");
}
static void test_stage3_exp_log_identity(void) {
    seed_rng(42);
    /* Exp[Log[x]] = x is an identity on the principal branch for all
     * non-zero complex x (Log defined as principal-value inverse). */
    assert_pzq("PossibleZeroQ[Exp[Log[x]] - x]", "True");
}
static void test_stage3_cos_double_angle(void) {
    seed_rng(42);
    assert_pzq("PossibleZeroQ[Cos[2 x] - (Cos[x]^2 - Sin[x]^2)]", "True");
}
static void test_stage3_tan_in_terms_of_sin_cos(void) {
    seed_rng(42);
    assert_pzq("PossibleZeroQ[Tan[x] - Sin[x]/Cos[x]]", "True");
}
static void test_stage3_sqrt_tan_derivative_identity(void) {
    /* D[Integrate[Sqrt[Tan[x]], x], x] == Sqrt[Tan[x]].  A true identity
     * that was misreported False because the Schwartz-Zippel sampler drew
     * a huge imaginary part (|Im| ~ 450), where Tan saturates to +-I and
     * 1 + Tan^2 (= Sec^2) catastrophically cancels to 0 in IEEE arithmetic
     * -- so the Sec^2-scaled derivative read as 0 while Sqrt[Tan[x]] stayed
     * O(1).  Bounding the sampled imaginary magnitude fixes it. */
    seed_rng(42);
    assert_pzq(
        "PossibleZeroQ[D[ArcTan[-1 + Sqrt[2] Sqrt[Tan[x]]]/Sqrt[2]"
        " + ArcTan[1 + Sqrt[2] Sqrt[Tan[x]]]/Sqrt[2]"
        " - 1/2 Log[1 + Tan[x] + Sqrt[2] Sqrt[Tan[x]]]/Sqrt[2]"
        " + 1/2 Log[1 + Tan[x] - Sqrt[2] Sqrt[Tan[x]]]/Sqrt[2], x]"
        " - Sqrt[Tan[x]]]", "True");
}
static void test_stage3_branch_dependent_nonzero(void) {
    /* D[2 Sqrt[1-Cos[x]], x] - Sqrt[1+Cos[x]] is identically 0 on (0,pi),
     * (2pi,3pi), ... but genuinely non-zero on (pi,2pi), (3pi,4pi), ... (the
     * two Sqrt factors sit on different analytic branches).  It is NOT an
     * identical zero, so PossibleZeroQ must report False.  With only 4
     * Schwartz-Zippel samples it returned True ~6% of the time -- whenever all
     * four random real parts happened to land in a zero interval.  A small
     * imaginary perturbation does not help (by the identity theorem the
     * function is zero on a whole complex neighbourhood of each zero interval);
     * defeating it requires many sample POINTS, which the machine-precision
     * screen now provides.  Sweep several seeds so the screen is exercised. */
    for (int64_t s = 1; s <= 8; ++s) {
        seed_rng(s);
        assert_pzq(
            "PossibleZeroQ[D[2 Sqrt[1 - Cos[x]], x] - Sqrt[1 + Cos[x]]]",
            "False");
    }
}

/* ============================================================== */
/*  5. Listable threading and registration                        */
/* ============================================================== */

static void test_listable_threading(void) {
    assert_pzq("PossibleZeroQ[{0, x - x, 1, Pi - Pi}]",
               "{True, True, False, True}");
}
static void test_listable_nested_list(void) {
    assert_pzq("PossibleZeroQ[{{0, 1}, {x - x, 2}}]",
               "{{True, False}, {True, False}}");
}
static void test_attributes_registered(void) {
    assert_pzq("Attributes[PossibleZeroQ]", "{Listable, Protected}");
}
static void test_arity_mismatch_unevaluated(void) {
    /* Wrong arity should leave the call unevaluated (returns NULL inside
     * the builtin); Listable could also thread but with 0 args there's
     * nothing to thread so we expect the original form back. */
    assert_pzq("PossibleZeroQ[]", "PossibleZeroQ[]");
}

/* ============================================================== */
/*  6. Stress tests — large polynomial identities                 */
/* ============================================================== */

static void test_stress_binomial_5(void) {
    /* (a + b)^5 = sum of 6 terms. */
    assert_pzq(
        "PossibleZeroQ[(a + b)^5 - (a^5 + 5 a^4 b + 10 a^3 b^2 + "
        "10 a^2 b^3 + 5 a b^4 + b^5)]",
        "True");
}
static void test_stress_binomial_10(void) {
    assert_pzq(
        "PossibleZeroQ[(a + b)^10 - "
        "(a^10 + 10 a^9 b + 45 a^8 b^2 + 120 a^7 b^3 + 210 a^6 b^4 + "
        "252 a^5 b^5 + 210 a^4 b^6 + 120 a^3 b^7 + 45 a^2 b^8 + "
        "10 a b^9 + b^10)]",
        "True");
}
static void test_stress_vieta_quartic(void) {
    /* (x - r1)(x - r2)(x - r3)(x - r4) expansion versus elementary
     * symmetric polynomial form. */
    assert_pzq(
        "PossibleZeroQ["
        "(x - r1)(x - r2)(x - r3)(x - r4) - "
        "(x^4 - (r1 + r2 + r3 + r4) x^3 + "
        "(r1 r2 + r1 r3 + r1 r4 + r2 r3 + r2 r4 + r3 r4) x^2 - "
        "(r1 r2 r3 + r1 r2 r4 + r1 r3 r4 + r2 r3 r4) x + "
        "r1 r2 r3 r4)]",
        "True");
}
static void test_stress_sophie_germain(void) {
    /* Sophie Germain identity: a^4 + 4 b^4 = (a^2 + 2 b^2 + 2 a b)(a^2 + 2 b^2 - 2 a b) */
    assert_pzq(
        "PossibleZeroQ[a^4 + 4 b^4 - "
        "(a^2 + 2 b^2 + 2 a b)(a^2 + 2 b^2 - 2 a b)]",
        "True");
}
static void test_stress_difference_of_cubes(void) {
    assert_pzq(
        "PossibleZeroQ[a^3 - b^3 - (a - b)(a^2 + a b + b^2)]", "True");
}
static void test_stress_sum_of_cubes(void) {
    assert_pzq(
        "PossibleZeroQ[a^3 + b^3 - (a + b)(a^2 - a b + b^2)]", "True");
}
static void test_stress_lagrange_identity_4d(void) {
    /* Lagrange's identity in 2 dimensions:
     * (a1^2 + a2^2)(b1^2 + b2^2) - (a1 b1 + a2 b2)^2 - (a1 b2 - a2 b1)^2 == 0 */
    assert_pzq(
        "PossibleZeroQ["
        "(a1^2 + a2^2)(b1^2 + b2^2) - (a1 b1 + a2 b2)^2 - "
        "(a1 b2 - a2 b1)^2]",
        "True");
}
static void test_stress_rational_telescoping(void) {
    /* 1/(k(k+1)) = 1/k - 1/(k+1) */
    assert_pzq(
        "PossibleZeroQ[1/(k(k + 1)) - (1/k - 1/(k + 1))]", "True");
}
static void test_stress_rational_three_var_sum(void) {
    /* Common-denominator collapse: a/(bc) + b/(ac) + c/(ab) - (a+b+c)^2/(abc)
     * + 2 (ab+ac+bc)/(abc) is identically 0. */
    assert_pzq(
        "PossibleZeroQ["
        "a/(b c) + b/(a c) + c/(a b) - (a + b + c)^2/(a b c) + "
        "2 (a b + a c + b c)/(a b c)]",
        "True");
}

/* ============================================================== */
/*  7. Stress tests — complex / trig / log                        */
/* ============================================================== */

static void test_stress_complex_conjugate_product(void) {
    /* (a + b I)(a - b I) - (a^2 + b^2) == 0 for real a, b. */
    seed_rng(42);
    assert_pzq("PossibleZeroQ[(a + b I)(a - b I) - (a^2 + b^2)]", "True");
}
static void test_stress_log_power_rule(void) {
    /* Log[x^2] - 2 Log[x] holds on the principal branch only when x is
     * positive real.  Random complex sampling reveals the branch
     * dependence, so the answer is False. */
    seed_rng(42);
    assert_pzq("PossibleZeroQ[Log[x^2] - 2 Log[x]]", "False");
}
static void test_stress_complex_e_to_ipi(void) {
    /* Euler's identity: E^(I Pi) + 1 == 0 */
    assert_pzq("PossibleZeroQ[E^(I Pi) + 1]", "True");
}
static void test_stress_sin_cos_pi_third(void) {
    /* Sin[Pi/3] - Sqrt[3]/2 == 0 */
    assert_pzq("PossibleZeroQ[Sin[Pi/3] - Sqrt[3]/2]", "True");
}
static void test_stress_cos_pi_quarter(void) {
    assert_pzq("PossibleZeroQ[Cos[Pi/4] - Sqrt[2]/2]", "True");
}
static void test_stress_sum_of_5th_roots_of_unity(void) {
    /* Sum of all five 5th roots of unity is 0. */
    assert_pzq(
        "PossibleZeroQ["
        "Sum[E^(2 Pi I k/5), {k, 0, 4}]]",
        "True");
}
static void test_stress_log_e_squared(void) {
    /* Log[E^2] - 2 == 0 (Log[Exp[2]] simplifies via the symbolic rule). */
    assert_pzq("PossibleZeroQ[Log[E^2] - 2]", "True");
}

/* ============================================================== */
/*  8. Catastrophic-cancellation regressions                      */
/* ============================================================== */

static void test_cancel_machine_overflow_1(void) {
    /* Same as Stage 2 case but at 30 orders of magnitude. */
    assert_pzq("PossibleZeroQ[10^30 + 1 - 10^30]", "False");
}
static void test_cancel_machine_overflow_2(void) {
    assert_pzq("PossibleZeroQ[2^100 + 1 - 2^100]", "False");
}
static void test_cancel_pair_diff(void) {
    /* Sum that is exactly 0 numerically; should *still* be True. */
    assert_pzq("PossibleZeroQ[10^30 - 10^30]", "True");
}
static void test_cancel_trig_near_zero_real(void) {
    /* Sin[Pi] is ~1e-16 at machine precision — must NOT mistake the
     * tiny non-zero in (Sin[Pi] + 1) for actual zero. */
    assert_pzq("PossibleZeroQ[Sin[Pi] + 1]", "False");
}

/* ============================================================== */
/*  9. Memory-leak smoke load                                     */
/* ============================================================== */

static void test_memory_smoke_load(void) {
    /* 500 invocations across all three stages.  Run under
     * `valgrind --leak-check=full` to confirm no per-call leaks. */
    seed_rng(42);
    for (int i = 0; i < 500; ++i) {
        Expr* p; Expr* r;
        const char* inputs[] = {
            "PossibleZeroQ[0]",
            "PossibleZeroQ[(x + 1)(x - 1) - x^2 + 1]",
            "PossibleZeroQ[E^Pi - Pi^E]",
            "PossibleZeroQ[Sin[x]^2 + Cos[x]^2 - 1]",
            "PossibleZeroQ[Sqrt[x^2] - x]"
        };
        const char* in = inputs[i % 5];
        p = parse_expression(in);
        r = evaluate(p);
        expr_free(p);
        expr_free(r);
    }
}

/* ============================================================== */
/*  Main driver                                                   */
/* ============================================================== */

int main(void) {
    symtab_init();
    core_init();

    /* Group 1 — Stage 0 */
    TEST(test_stage0_integer_zero);
    TEST(test_stage0_real_zero);
    TEST(test_stage0_real_negzero);
    TEST(test_stage0_integer_one);
    TEST(test_stage0_integer_negone);
    TEST(test_stage0_bigint);
    TEST(test_stage0_bigint_negative);
    TEST(test_stage0_complex_zero);
    TEST(test_stage0_complex_pure_imag);
    TEST(test_stage0_pi);
    TEST(test_stage0_e);
    TEST(test_stage0_eulergamma);
    TEST(test_stage0_free_symbol);
    TEST(test_stage0_rational_zero);
    TEST(test_stage0_rational_nonzero);

    /* Group 2 — Stage 1 */
    TEST(test_stage1_diff_of_squares);
    TEST(test_stage1_cube_expansion);
    TEST(test_stage1_common_denominator);
    TEST(test_stage1_cancel_linear);
    TEST(test_stage1_perfect_square);
    TEST(test_stage1_obvious_nonzero);
    TEST(test_stage1_nonzero_constant);
    TEST(test_stage1_triple_product);
    TEST(test_stage1_partial_fraction);
    TEST(test_stage1_three_var_rational);

    /* Group 3 — Stage 2 */
    TEST(test_stage2_exp_ipi_quarter);
    TEST(test_stage2_binomial_e_pi);
    TEST(test_stage2_e_pi_vs_pi_e);
    TEST(test_stage2_sin_pi);
    TEST(test_stage2_cos_pi_half);
    TEST(test_stage2_sqrt_two_squared);
    TEST(test_stage2_pyth_at_one);
    TEST(test_stage2_double_angle_sin_at_one);
    TEST(test_stage2_complex_two_2i);
    TEST(test_stage2_log_product);
    TEST(test_stage2_log_quotient);
    TEST(test_stage2_tiny_nonzero);
    TEST(test_stage2_sqrt_irrational_close);
    TEST(test_stage2_catastrophic_cancellation);
    TEST(test_stage2_complex_conjugate);
    TEST(test_stage2_complex_imag_squared);

    /* Group 4 — Stage 3 */
    TEST(test_stage3_pythagorean_identity);
    TEST(test_stage3_double_angle_sin);
    TEST(test_stage3_sqrt_x_squared_minus_x);
    TEST(test_stage3_exp_log_identity);
    TEST(test_stage3_cos_double_angle);
    TEST(test_stage3_tan_in_terms_of_sin_cos);
    TEST(test_stage3_sqrt_tan_derivative_identity);
    TEST(test_stage3_branch_dependent_nonzero);

    /* Group 5 — Listable / registration */
    TEST(test_listable_threading);
    TEST(test_listable_nested_list);
    TEST(test_attributes_registered);
    TEST(test_arity_mismatch_unevaluated);

    /* Group 6 — Large polynomial identities */
    TEST(test_stress_binomial_5);
    TEST(test_stress_binomial_10);
    TEST(test_stress_vieta_quartic);
    TEST(test_stress_sophie_germain);
    TEST(test_stress_difference_of_cubes);
    TEST(test_stress_sum_of_cubes);
    TEST(test_stress_lagrange_identity_4d);
    TEST(test_stress_rational_telescoping);
    TEST(test_stress_rational_three_var_sum);

    /* Group 7 — Complex / trig / log */
    TEST(test_stress_complex_conjugate_product);
    TEST(test_stress_log_power_rule);
    TEST(test_stress_complex_e_to_ipi);
    TEST(test_stress_sin_cos_pi_third);
    TEST(test_stress_cos_pi_quarter);
    TEST(test_stress_sum_of_5th_roots_of_unity);
    TEST(test_stress_log_e_squared);

    /* Group 8 — Catastrophic cancellation regressions */
    TEST(test_cancel_machine_overflow_1);
    TEST(test_cancel_machine_overflow_2);
    TEST(test_cancel_pair_diff);
    TEST(test_cancel_trig_near_zero_real);

    /* Group 9 — Memory smoke */
    TEST(test_memory_smoke_load);

    printf("\nAll PossibleZeroQ tests passed.\n");
    return 0;
}
