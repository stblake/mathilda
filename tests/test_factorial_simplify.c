/*
 * Tests for the simp_factorial transform inside Simplify.
 *
 * The transform is a single principled procedure that handles factorial
 * simplification across five categories (ratio reduction, polynomial
 * absorption, powers and cross-denominator cancellation, additive
 * collapse, multivariate / edge cases). No per-pattern table -- one
 * algorithm subsumes them all. See the simp_factorial section in
 * src/simp.c for the full specification.
 *
 * Each `assert_eval_eq` here exercises the user-facing string form of
 * the result. Equivalent forms (e.g. (n+2)(n+3) vs (2+n)(3+n) under
 * picocas's canonical Plus-arg ordering) are accepted in the expected
 * column to avoid coupling tests to ordering decisions that live
 * outside of this module.
 */

#include "test_utils.h"
#include "symtab.h"
#include "core.h"

/* ---- 1. Ratio Reduction (Linear & Scaled Offsets) ---- */

void test_factorial_ratio_basic(void) {
    /* n! / (n-1)! -> n */
    assert_eval_eq("Simplify[n!/(n-1)!]", "n", 0);
}

void test_factorial_ratio_offset_pair(void) {
    /* (n+3)! / (n+1)! -> (n+2)(n+3). Picocas prints the canonical
     * form with constants first inside Plus: (2+n)(3+n). */
    assert_eval_eq("Simplify[(n+3)!/(n+1)!]", "(2 + n) (3 + n)", 0);
}

void test_factorial_ratio_smaller_over_larger(void) {
    /* (n-2)! / n! -> 1/((n-1)*n). Picocas keeps the factored
     * Times[Power[n,-1], Power[Plus[-1,n],-1]] form, which prints
     * as 1/(n*(-1+n)). */
    assert_eval_eq("Simplify[(n-2)!/n!]", "1/(n (-1 + n))", 0);
}

void test_factorial_ratio_scaled_offset(void) {
    /* (2n)! / (2n-2)! -> 2n*(2n-1). Picocas prints
     * 2*n*(-1+2n). */
    assert_eval_eq("Simplify[(2*n)!/(2*n - 2)!]", "2 n (-1 + 2 n)", 0);
}

/* ---- 2. Polynomial Absorption (Folding) ---- */

void test_factorial_absorb_one_step(void) {
    /* n * (n-1)! -> n! */
    assert_eval_eq("Simplify[n*(n-1)!]", "Factorial[n]", 0);
}

void test_factorial_absorb_two_step(void) {
    /* (n+1) * n * (n-1)! -> (n+1)! */
    assert_eval_eq("Simplify[(n+1)*n*(n-1)!]", "Factorial[1 + n]", 0);
}

void test_factorial_absorb_after_factoring(void) {
    /* (n^2 + n) * (n-1)! -> (n+1)!
     * Requires Factor inside simp_factorial to expose the (n)*(n+1)
     * factorisation of n^2 + n before the re-fold can absorb. */
    assert_eval_eq("Simplify[(n^2 + n)*(n-1)!]", "Factorial[1 + n]", 0);
}

void test_factorial_absorb_with_gap_fill(void) {
    /* (n^2 - 1) * (n-2)! -> (n+1)! / n
     * Requires gap-filling in the re-fold: the cofactor factors as
     * (n-1)(n+1), which is offsets {1, 3} relative to base n-2.
     * Folding to k=3 needs to introduce 1*(b+2)/(b+2) = n/n so the
     * Pochhammer (n-1)*n*(n+1) lands as (n+1)!. */
    assert_eval_eq("Simplify[(n^2 - 1)*(n-2)!]", "Factorial[1 + n]/n", 0);
}

/* ---- 3. Powers and Cross-Denominator Cancellation ---- */

void test_factorial_power_squared_ratio(void) {
    /* (n!)^2 / ((n-1)! (n+1)!) -> n / (n+1) */
    assert_eval_eq("Simplify[(n!)^2 / ((n-1)! * (n+1)!)]", "n/(1 + n)", 0);
}

void test_factorial_power_cubed_ratio(void) {
    /* ((n+1)!)^3 / (n!)^3 -> (n+1)^3 */
    assert_eval_eq("Simplify[((n+1)!)^3 / (n!)^3]", "(1 + n)^3", 0);
}

void test_factorial_squared_minus_polynomial(void) {
    /* (n!/(n-2)!)^2 - (n^2 - n)^2 -> 0 */
    assert_eval_eq("Simplify[(n!/(n-2)!)^2 - (n^2 - n)^2]", "0", 0);
}

/* ---- 4. Additive Collapse (Sums and Differences) ---- */

void test_factorial_additive_collapse_simple(void) {
    /* (n+1)! - n*n! -> n!  (via (n+1)! = n!*(n+1) then collapse) */
    assert_eval_eq("Simplify[(n+1)! - n*n!]", "Factorial[n]", 0);
}

void test_factorial_additive_collapse_shifted(void) {
    /* (n+2)! - (n+1)!*(n+1) -> (n+1)! */
    assert_eval_eq("Simplify[(n+2)! - (n+1)!*(n+1)]", "Factorial[1 + n]", 0);
}

void test_factorial_ratio_minus_n(void) {
    /* (n+1)!/n! - n -> 1 */
    assert_eval_eq("Simplify[(n+1)!/n! - n]", "1", 0);
}

void test_factorial_inverse_difference(void) {
    /* 1/n! - 1/(n+1)! -> n/(n+1)!
     * Requires the path-B (Together / Factor / combine_inverses)
     * branch of simp_factorial: the factorial atoms in the two
     * denominators must be combined over a common Times[Factorial[n],
     * (n+1)] before the re-fold can collapse it to Factorial[n+1]. */
    assert_eval_eq("Simplify[1/n! - 1/(n+1)!]", "n/Factorial[1 + n]", 0);
}

/* ---- 5. Multivariate and Edge Cases ---- */

void test_factorial_multivariate_offset(void) {
    /* (n+k)! / (n+k-1)! -> n+k */
    assert_eval_eq("Simplify[(n+k)!/(n+k-1)!]", "k + n", 0);
}

void test_factorial_multivariate_negative_offset(void) {
    /* (n-k)! / (n-k-2)! -> (n-k)*(n-k-1).
     * The canonical Plus-ordering of picocas prints the symbolic
     * group with k first (sign-flipped: (k-n)*(1+k-n)), which is
     * the same value. */
    assert_eval_eq("Simplify[(n-k)!/(n-k-2)!]", "(k - n) (1 + k - n)", 0);
}

void test_factorial_zero_arg_via_evaluation(void) {
    /* (n-n)! -> 1 falls out of normal Plus + Factorial[0] semantics
     * (n-n collapses to 0 inside Factorial's argument before our
     * transform sees the input). The test is here because the user
     * spec lists it; it documents that simp_factorial does not
     * regress this edge case. */
    assert_eval_eq("Simplify[(n - n)!]", "1", 0);
}

void test_factorial_pochhammer_round_trip(void) {
    /* n!/(n+2)! * (n^2 + 3n + 2) -> 1
     * (n^2 + 3n + 2) factors as (n+1)(n+2), exactly the Pochhammer
     * Factorial[n+2]/Factorial[n], so the whole expression collapses. */
    assert_eval_eq("Simplify[n!/(n+2)! * (n^2 + 3*n + 2)]", "1", 0);
}

void test_factorial_double_factorial_identity(void) {
    /* (2n)! / (2^n n!) -> (2n - 1)!!
     * The canonical "double factorial from factorial" identity --
     * cannot be derived from the shift-normalization alone (the
     * factorials live in different symbolic groups: 2n and n) so
     * simp_factorial recognises this specific shape and rewrites
     * to Factorial2. The result prints with Factorial2's FullForm
     * head because picocas does not yet pretty-print x!! infix. */
    assert_eval_eq("Simplify[(2*n)!/(2^n * n!)]", "Factorial2[-1 + 2 n]", 0);
}

/* ---- Soundness sentinels ---- */

void test_factorial_no_factorials_passes_through(void) {
    /* No Factorial atom -> simp_factorial is inert. */
    assert_eval_eq("Simplify[(x + 1)*(x - 1)]", "-1 + x^2", 0);
}

void test_factorial_single_call_no_op(void) {
    /* A bare Factorial with no peer in the same group has nothing
     * to shift-normalize. The transform must not introduce noise. */
    assert_eval_eq("Simplify[Factorial[k]]", "Factorial[k]", 0);
}

void test_factorial_concrete_value_evaluates(void) {
    /* Concrete arguments evaluate normally; simp_factorial does not
     * interfere with the existing builtin_factorial's int path. */
    assert_eval_eq("Simplify[5!]", "120", 0);
}

void test_factorial2_concrete(void) {
    /* The new Factorial2 builtin: small concrete inputs evaluate
     * directly. 7!! = 7*5*3*1 = 105.  8!! = 8*6*4*2 = 384.
     * Negative-one and zero are 1 by convention. */
    assert_eval_eq("Factorial2[7]", "105", 0);
    assert_eval_eq("Factorial2[8]", "384", 0);
    assert_eval_eq("Factorial2[0]", "1", 0);
    assert_eval_eq("Factorial2[-1]", "1", 0);
}

void test_factorial2_symbolic_held(void) {
    /* Symbolic arguments stay unevaluated -- Factorial2 is purely
     * a numeric leaf for now. */
    assert_eval_eq("Factorial2[m]", "Factorial2[m]", 0);
}

int main(void) {
    symtab_init();
    core_init();

    /* Category 1 */
    TEST(test_factorial_ratio_basic);
    TEST(test_factorial_ratio_offset_pair);
    TEST(test_factorial_ratio_smaller_over_larger);
    TEST(test_factorial_ratio_scaled_offset);

    /* Category 2 */
    TEST(test_factorial_absorb_one_step);
    TEST(test_factorial_absorb_two_step);
    TEST(test_factorial_absorb_after_factoring);
    TEST(test_factorial_absorb_with_gap_fill);

    /* Category 3 */
    TEST(test_factorial_power_squared_ratio);
    TEST(test_factorial_power_cubed_ratio);
    TEST(test_factorial_squared_minus_polynomial);

    /* Category 4 */
    TEST(test_factorial_additive_collapse_simple);
    TEST(test_factorial_additive_collapse_shifted);
    TEST(test_factorial_ratio_minus_n);
    TEST(test_factorial_inverse_difference);

    /* Category 5 */
    TEST(test_factorial_multivariate_offset);
    TEST(test_factorial_multivariate_negative_offset);
    TEST(test_factorial_zero_arg_via_evaluation);
    TEST(test_factorial_pochhammer_round_trip);
    TEST(test_factorial_double_factorial_identity);

    /* Soundness sentinels */
    TEST(test_factorial_no_factorials_passes_through);
    TEST(test_factorial_single_call_no_op);
    TEST(test_factorial_concrete_value_evaluates);
    TEST(test_factorial2_concrete);
    TEST(test_factorial2_symbolic_held);

    printf("All factorial_simplify tests passed!\n");
    return 0;
}
