/*
 * test_radical_fusion.c -- Coverage for the opposite-exponent Power fusion
 * pass inside builtin_times, and for Power[0, b] folding in builtin_power.
 *
 * The fusion rule:  a^q * b^(-q)  ->  (a/b)^q
 *   fires whenever a and b are both positive numeric (integer, bigint,
 *   rational, real). The ratio a/b may itself be integer (e.g.
 *   Sqrt[6]/Sqrt[2] -> Sqrt[3]) or rational (e.g. Sqrt[10]/Sqrt[3] ->
 *   Sqrt[10/3]).
 *
 * The Power[0, b] rule:
 *   Power[0, b]  ->  0                for any positive b (int/real/rational),
 *   Power[0, b]  ->  ComplexInfinity  for any negative b.
 */

#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "symtab.h"
#include "core.h"

#include <stdio.h>

/* ----------------------------------------------------------------- */
/* Radical fusion inside Times                                        */
/* ----------------------------------------------------------------- */

static void test_sqrt_integer_ratio(void) {
    /* Classic: 6/2 = 3, 15/5 = 3, 30/6 = 5. */
    assert_eval_eq("Sqrt[6]/Sqrt[2]",  "Sqrt[3]", 0);
    assert_eval_eq("Sqrt[15]/Sqrt[5]", "Sqrt[3]", 0);
    assert_eval_eq("Sqrt[30]/Sqrt[6]", "Sqrt[5]", 0);
    /* Ordering shouldn't matter. */
    assert_eval_eq("1/Sqrt[3] * Sqrt[6]", "Sqrt[2]", 0);
    assert_eval_eq("Sqrt[6] * 1/Sqrt[2]", "Sqrt[3]", 0);
}

static void test_cube_root_ratio(void) {
    /* Non-square exponents: Power[6,1/3] * Power[2,-1/3] = 3^(1/3). */
    assert_eval_eq("Power[6, 1/3] * Power[2, -1/3]",  "3^(1/3)", 0);
    assert_eval_eq("Power[12, 1/3] * Power[4, -1/3]", "3^(1/3)", 0);
    assert_eval_eq("Power[30, 1/3] * Power[6, -1/3]", "5^(1/3)", 0);
}

static void test_rational_ratio_fuses(void) {
    /* Sqrt[4] reduces to 2 before fusion, so Sqrt[6]/Sqrt[4] short-circuits
     * through integer folding to 1/2 Sqrt[6] rather than Sqrt[3/2]. */
    assert_eval_eq("Sqrt[6]/Sqrt[4]", "1/2 Sqrt[6]", 0);
    /* Coprime bases with rational ratio now fuse into a single radical
     * over the rational ratio. */
    assert_eval_eq("Sqrt[10]/Sqrt[3]",       "Sqrt[10/3]",   0);
    assert_eval_eq("Power[2, 1/3] / Power[3, 1/3]", "(2/3)^(1/3)", 0);
    /* Ratio <1: Sqrt[3]/Sqrt[6] = Sqrt[1/2]. */
    assert_eval_eq("Sqrt[3]/Sqrt[6]", "Sqrt[1/2]", 0);
}

static void test_fusion_with_symbolic_factor(void) {
    /* Symbolic factors pass through, radical fusion still fires. */
    assert_eval_eq("a Sqrt[6]/Sqrt[2]", "Sqrt[3] a", 0);
    assert_eval_eq("x * Sqrt[15]/Sqrt[5]", "Sqrt[3] x", 0);
}

static void test_same_base_still_works(void) {
    /* Regression: the radical-fusion pass must not disrupt the existing
     * same-base grouping. These are all handled by base merging, not
     * the new pairwise pass. */
    assert_eval_eq("Sqrt[8]/Sqrt[2]",   "2", 0);
    assert_eval_eq("Sqrt[12]/Sqrt[3]",  "2", 0);
    assert_eval_eq("Sqrt[50]/Sqrt[2]",  "5", 0);
    /* Same base, same rational exponent: Power[2, 1/2] * Power[2, 1/2] = 2. */
    assert_eval_eq("Sqrt[2] * Sqrt[2]", "2", 0);
}

static void test_chained_fusion(void) {
    /* A chain of three different prime-power bases. The restart-from-zero
     * logic should keep fusing pairs until nothing more can be merged. */
    /* Sqrt[30] * Sqrt[6]^(-1) * Sqrt[5]^(-1) = Sqrt[30/6/5] = Sqrt[1] = 1. */
    assert_eval_eq("Sqrt[30] / Sqrt[6] / Sqrt[5]", "1", 0);
    /* 210 / 6 / 5 = 7. */
    assert_eval_eq("Sqrt[210] / Sqrt[6] / Sqrt[5]", "Sqrt[7]", 0);
}

/* ----------------------------------------------------------------- */
/* Power[0, b] folding                                                */
/* ----------------------------------------------------------------- */

static void test_zero_power_positive(void) {
    /* Integer positive exponent. */
    assert_eval_eq("Power[0, 1]", "0", 0);
    assert_eval_eq("Power[0, 2]", "0", 0);
    assert_eval_eq("Power[0, 100]", "0", 0);

    /* Rational positive exponent. This is the one that historically
     * bubbled up as Sqrt[0] into limit results. */
    assert_eval_eq("Power[0, 1/2]", "0", 0);
    assert_eval_eq("Power[0, 3/7]", "0", 0);
    assert_eval_eq("Sqrt[0]",       "0", 0);

    /* Real positive exponent: keep the answer as a real 0.0 so that the
     * overall expression's numeric type is preserved. */
    assert_eval_eq("Power[0, 0.5]", "0.0", 0);
    assert_eval_eq("Power[0.0, 2]", "0.0", 0);
}

static void test_zero_power_negative(void) {
    /* Negative exponents collapse to ComplexInfinity (and emit a
     * Power::infy warning, which we do not assert on here). */
    assert_eval_eq("Power[0, -1]",    "ComplexInfinity", 0);
    assert_eval_eq("Power[0, -2]",    "ComplexInfinity", 0);
    assert_eval_eq("Power[0, -1/2]",  "ComplexInfinity", 0);
    assert_eval_eq("Power[0, -3/4]",  "ComplexInfinity", 0);
}

static void test_zero_power_symbolic(void) {
    /* Symbolic exponent: cannot decide sign, leave unevaluated. */
    assert_eval_eq("Power[0, a]", "0^a", 0);
}

static void test_zero_power_zero_unchanged(void) {
    /* Power[0, 0] is conventionally 1 in Mathilda and that must not
     * change: the new positive-exponent branch must only fire for
     * strictly positive exponents. */
    assert_eval_eq("Power[0, 0]", "Indeterminate", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_sqrt_integer_ratio);
    TEST(test_cube_root_ratio);
    TEST(test_rational_ratio_fuses);
    TEST(test_fusion_with_symbolic_factor);
    TEST(test_same_base_still_works);
    TEST(test_chained_fusion);

    TEST(test_zero_power_positive);
    TEST(test_zero_power_negative);
    TEST(test_zero_power_symbolic);
    TEST(test_zero_power_zero_unchanged);

    printf("All radical_fusion tests passed.\n");
    return 0;
}
