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
     * Mathilda's printer emits Power[Power[a, p], q] without disambiguating
     * parentheses, so the result string reads "a^p^q". This is the
     * unchanged input, NOT a^(p^q). */
    assert_eval_eq("Simplify[(a^p)^q, a > 0]",
                   "a^p^q", 0);
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

/* ---- Log[Exp[...]] inverse pair ---- */

void test_log_exp_positive(void) {
    /* Log[Exp[x]] for any x: works through the existing TrigToExp /
     * ExpToTrig roundtrip rather than the new cascade, but listed here
     * as a sanity check that the cascade does not break it. */
    assert_eval_eq("Simplify[Log[Exp[x]], x > 0]", "x", 0);
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

    TEST(test_log_exp_positive);

    printf("All logexp Simplify tests passed!\n");
    return 0;
}
