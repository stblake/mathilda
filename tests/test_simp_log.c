#include "test_utils.h"
#include "symtab.h"
#include "core.h"

/*
 * tests/test_simp_log.c -- Simplify cases driven by src/simp/simp_log.c.
 *
 * The suite mirrors the 17-case logTestSuite from the user that motivated
 * the simp_log module. Each test exercises one of the two primitives:
 *
 *   Pass A -- prime decomposition of Log[positive rational], the engine
 *             that turns Log[4] -> 2 Log[2], Log[72] -> 3 Log[2] + 2 Log[3],
 *             Log[2/3] -> Log[2] - Log[3], Log[Sqrt[12]] -> Log[2] + 1/2 Log[3].
 *
 *   Pass B -- the linear-combination-of-logs fuser, which collapses
 *             Sum c_i Log[a_i] into a single Log[Product a_i^c_i] when
 *             the fused form (after Together / Expand / one bounded
 *             Simplify) is strictly smaller.
 *
 * The trailing tests (12, 13) require Element[_, Reals] assumptions to
 * collapse Log[E^a] -> a for symbolic a; the user explicitly noted that
 * "some will require conditions".
 */

/* ---- Pass A: prime decomposition + cancellation ---- */

void test_simp_log_base_case_integer(void) {
    assert_eval_eq("Simplify[-4 Log[4] + 8 Log[2]]", "0", 0);
}

void test_simp_log_expanded_integer(void) {
    assert_eval_eq("Simplify[-2 Log[9] - 4 Log[4] + 4 Log[3] + 8 Log[2]]",
                   "0", 0);
}

void test_simp_log_rational_arguments(void) {
    assert_eval_eq("Simplify[Log[2/3] + Log[3/4] + Log[4/5] - Log[2/5]]",
                   "0", 0);
}

void test_simp_log_prime_decomposition(void) {
    assert_eval_eq("Simplify[Log[72] - 3 Log[2] - 2 Log[3]]", "0", 0);
}

void test_simp_log_radical_argument(void) {
    /* Log[Sqrt[12]] expands as 1/2 Log[12] -> 1/2 (2 Log[2] + Log[3])
     *                          = Log[2] + 1/2 Log[3]. */
    assert_eval_eq("Simplify[Log[Sqrt[12]] - Log[2] - 1/2 Log[3]]", "0", 0);
}

void test_simp_log_change_of_base(void) {
    /* Log[2,5] Log[5,8] = (Log[5]/Log[2])(Log[8]/Log[5]) = Log[8]/Log[2].
     * With Log[8] -> 3 Log[2] the ratio collapses to 3. */
    assert_eval_eq("Simplify[Log[2,5]*Log[5,8] - 3]", "0", 0);
}

/* ---- Pass B: algebraic fuser (positive symbolic args) ---- */

void test_simp_log_difference_of_squares(void) {
    /* The fused product (x^2-y^2)/((x-y)(x+y)) -> 1 via Together,
     * then Log[1] = 0. The fusion is taken under the "small constant"
     * gate; no positivity assumption needed. */
    assert_eval_eq("Simplify[Log[x^2 - y^2] - Log[x - y] - Log[x + y]]",
                   "0", 0);
}

void test_simp_log_perfect_square_trinomial(void) {
    assert_eval_eq("Simplify[Log[x^2 + 2 x + 1] - 2 Log[x + 1]]",
                   "0", 0);
}

void test_simp_log_difference_of_cubes(void) {
    assert_eval_eq("Simplify[Log[x^3 - 1] - Log[x - 1] - Log[x^2 + x + 1]]",
                   "0", 0);
}

void test_simp_log_golden_ratio(void) {
    /* (1+Sqrt[5])/2 and (Sqrt[5]-1)/2 are both positive; product = 1. */
    assert_eval_eq("Simplify[Log[(1 + Sqrt[5])/2] + Log[(Sqrt[5] - 1)/2]]",
                   "0", 0);
}

/* ---- Pass B + LogExp interaction (Log[E^a] feedback) ---- */

void test_simp_log_linear_nesting(void) {
    /* Log[E^(3x)] = 3x for numeric scalar exponents already; the
     * symbolic-coefficient case Log[E^(3x)] - 3 Log[E^x] reduces to 0
     * via SimpLogRules + LogExpRules iterated in simp_pipeline_logexp. */
    assert_eval_eq("Simplify[Log[E^(3 x)] - 3 Log[E^x]]", "0", 0);
}

void test_simp_log_exp_product_with_realness(void) {
    /* Element[_, Reals] assumption is required so the LogExp cascade
     * collapses Log[E^(x+y)] to x+y. */
    assert_eval_eq("Simplify[Log[E^x * E^y] - x - y, Element[x, Reals] && Element[y, Reals]]",
                   "0", 0);
}

void test_simp_log_reciprocal_exp_with_realness(void) {
    /* Fusion gives Log[E^(-x)] - x; LogExp cascade then -x - x = -2x. */
    assert_eval_eq("Simplify[Log[1/(1 + E^x)] + Log[1 + E^(-x)] - x, Element[x, Reals]]",
                   "-2 x", 0);
}

/* ---- Hyperbolic / Euler identities (already pass without simp_log;
 *      kept here as regression guards). ---- */

void test_simp_log_hyperbolic_defang(void) {
    assert_eval_eq("Simplify[Log[E^x + Sqrt[E^(2 x) + 1]] - ArcSinh[E^x]]",
                   "0", 0);
}

void test_simp_log_euler_complex_zero(void) {
    assert_eval_eq("Simplify[Log[E^(I x)] - Log[Cos[x] + I Sin[x]]]",
                   "0", 0);
}

/* ---- Pass B over trig log args ---- */

void test_simp_log_trig_complement(void) {
    /* Sin / (Tan * Cos) -> 1 via the evaluator's Trig auto-simplification;
     * Log[1] = 0. */
    assert_eval_eq("Simplify[Log[Sin[x]] - Log[Tan[x]] - Log[Cos[x]]]",
                   "0", 0);
}

void test_simp_log_pythagorean_conjugates(void) {
    /* (Sec[x]-Tan[x])(Sec[x]+Tan[x]) = Sec[x]^2 - Tan[x]^2 -> 1 via the
     * bounded Simplify pass on the fused argument; Log[1] = 0. */
    assert_eval_eq("Simplify[Log[Sec[x] - Tan[x]] + Log[Sec[x] + Tan[x]]]",
                   "0", 0);
}

/* ---- Regression guards: confirm we don't over-expand ---- */

void test_simp_log_single_log_unchanged(void) {
    /* A single Log[positive rational] decomposes; but a single
     * Log[prime] stays as Log[2]. */
    assert_eval_eq("Simplify[Log[2]]", "Log[2]", 0);
}

void test_simp_log_log_one_zero(void) {
    assert_eval_eq("Simplify[Log[1]]", "0", 0);
}

void test_simp_log_symbol_unchanged(void) {
    /* Pass B must NOT fuse a single Log term in a Plus. */
    assert_eval_eq("Simplify[Log[x] + 1]", "1 + Log[x]", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_simp_log_base_case_integer);
    TEST(test_simp_log_expanded_integer);
    TEST(test_simp_log_rational_arguments);
    TEST(test_simp_log_prime_decomposition);
    TEST(test_simp_log_radical_argument);
    TEST(test_simp_log_change_of_base);

    TEST(test_simp_log_difference_of_squares);
    TEST(test_simp_log_perfect_square_trinomial);
    TEST(test_simp_log_difference_of_cubes);
    TEST(test_simp_log_golden_ratio);

    TEST(test_simp_log_linear_nesting);
    TEST(test_simp_log_exp_product_with_realness);
    TEST(test_simp_log_reciprocal_exp_with_realness);

    TEST(test_simp_log_hyperbolic_defang);
    TEST(test_simp_log_euler_complex_zero);

    TEST(test_simp_log_trig_complement);
    TEST(test_simp_log_pythagorean_conjugates);

    TEST(test_simp_log_single_log_unchanged);
    TEST(test_simp_log_log_one_zero);
    TEST(test_simp_log_symbol_unchanged);

    printf("All simp_log Simplify tests passed!\n");
    return 0;
}
