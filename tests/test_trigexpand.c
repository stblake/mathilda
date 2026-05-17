#include "test_utils.h"
#include "symtab.h"
#include "core.h"

/*
 * Tests for TrigExpand. TrigExpand is a rule-based expansion, so the output
 * ordering follows Mathilda's canonical Orderless sort (which may differ from
 * Mathematica's pretty-printed form). Each expected string here reflects
 * Mathilda's canonical order.
 */

void test_trigexpand_double_angle(void) {
    assert_eval_eq("TrigExpand[Sin[2x]]", "2 Cos[x] Sin[x]", 0);
    assert_eval_eq("TrigExpand[Cos[2x]]", "Cos[x]^2 - Sin[x]^2", 0);
    assert_eval_eq("TrigExpand[Sinh[2x]]", "2 Cosh[x] Sinh[x]", 0);
    assert_eval_eq("TrigExpand[Cosh[2x]]", "Cosh[x]^2 + Sinh[x]^2", 0);
}

void test_trigexpand_triple_angle(void) {
    assert_eval_eq("TrigExpand[Sin[3x]]", "-Sin[x]^3 + 3 Cos[x]^2 Sin[x]", 0);
    assert_eval_eq("TrigExpand[Cos[3x]]", "Cos[x]^3 - 3 Cos[x] Sin[x]^2", 0);
    assert_eval_eq("TrigExpand[Sinh[3x]]", "Sinh[x]^3 + 3 Cosh[x]^2 Sinh[x]", 0);
    assert_eval_eq("TrigExpand[Cosh[3x]]", "Cosh[x]^3 + 3 Cosh[x] Sinh[x]^2", 0);
}

void test_trigexpand_higher_angle(void) {
    /* Sin[4x] = 4 Cos[x]^3 Sin[x] - 4 Cos[x] Sin[x]^3 */
    assert_eval_eq("TrigExpand[Sin[4x]]",
                   "-4 Cos[x] Sin[x]^3 + 4 Cos[x]^3 Sin[x]", 0);
    /* Sinh[4x] = 4 Cosh[x]^3 Sinh[x] + 4 Cosh[x] Sinh[x]^3 */
    assert_eval_eq("TrigExpand[Sinh[4x]]",
                   "4 Cosh[x] Sinh[x]^3 + 4 Cosh[x]^3 Sinh[x]", 0);
}

void test_trigexpand_negative_integer_multiple(void) {
    /* Sin[-2x] auto-simplifies (odd) then expands. */
    assert_eval_eq("TrigExpand[Sin[-2x]]", "-2 Cos[x] Sin[x]", 0);
    /* Cos is even: Cos[-2x] -> Cos[2x] -> Cos[x]^2 - Sin[x]^2. */
    assert_eval_eq("TrigExpand[Cos[-2x]]", "Cos[x]^2 - Sin[x]^2", 0);
    assert_eval_eq("TrigExpand[Sinh[-3x]]",
                   "-Sinh[x]^3 - 3 Cosh[x]^2 Sinh[x]", 0);
    assert_eval_eq("TrigExpand[Cosh[-3x]]",
                   "Cosh[x]^3 + 3 Cosh[x] Sinh[x]^2", 0);
}

void test_trigexpand_angle_addition(void) {
    assert_eval_eq("TrigExpand[Sin[x+y]]", "Cos[x] Sin[y] + Sin[x] Cos[y]", 0);
    assert_eval_eq("TrigExpand[Cos[x+y]]", "Cos[x] Cos[y] - Sin[x] Sin[y]", 0);
    assert_eval_eq("TrigExpand[Sinh[x+y]]",
                   "Cosh[x] Sinh[y] + Sinh[x] Cosh[y]", 0);
    assert_eval_eq("TrigExpand[Cosh[x-y]]",
                   "Cosh[x] Cosh[y] - Sinh[x] Sinh[y]", 0);
}

void test_trigexpand_multi_term_sum(void) {
    assert_eval_eq("TrigExpand[Cos[x+y+z]]",
        "Cos[x] Cos[y] Cos[z] - Cos[x] Sin[y] Sin[z] - Sin[x] Cos[y] Sin[z] - Sin[x] Sin[y] Cos[z]", 0);
    assert_eval_eq("TrigExpand[Sin[x+y+z]]",
        "-Sin[x] Sin[y] Sin[z] + Cos[x] Cos[y] Sin[z] + Cos[x] Sin[y] Cos[z] + Sin[x] Cos[y] Cos[z]", 0);
}

void test_trigexpand_mixed_sum(void) {
    /* Sinh[x - 2y] = Cosh[y]^2 Sinh[x] + Sinh[x] Sinh[y]^2 - 2 Cosh[x] Cosh[y] Sinh[y] */
    assert_eval_eq("TrigExpand[Sinh[x-2y]]",
        "Sinh[x] Sinh[y]^2 + Sinh[x] Cosh[y]^2 - 2 Cosh[x] Cosh[y] Sinh[y]", 0);
}

void test_trigexpand_pythagorean(void) {
    assert_eval_eq("TrigExpand[Sin[x]^2+Cos[x]^2]", "1", 0);
    assert_eval_eq("TrigExpand[Cosh[x]^2-Sinh[x]^2]", "1", 0);
    /* Squares of expanded multiple-angle forms collapse via the
     * Factor-based Pythagorean reduction:
     *   Sin[n x]^2 + Cos[n x]^2  Expand  (Sin[x]^2+Cos[x]^2)^n  Factor  1. */
    assert_eval_eq("TrigExpand[Sin[2x]^2+Cos[2x]^2]", "1", 0);
    assert_eval_eq("TrigExpand[Sin[3x]^2+Cos[3x]^2]", "1", 0);
    assert_eval_eq("TrigExpand[Sin[4x]^2+Cos[4x]^2]", "1", 0);
    /* Negative-integer arguments fold to positive via parity before expansion. */
    assert_eval_eq("TrigExpand[Sin[-4x]^2+Cos[-4x]^2]", "1", 0);
    /* Negated Pythagorean identity collapses to (-1)^n regardless of inner sign. */
    assert_eval_eq("TrigExpand[-Sin[2x]^2-Cos[2x]^2]", "-1", 0);
    assert_eval_eq("TrigExpand[-Sin[4x]^2-Cos[4x]^2]", "-1", 0);
    /* Scalar multiples pass through the collapse untouched. */
    assert_eval_eq("TrigExpand[5 Sin[4x]^2+5 Cos[4x]^2]", "5", 0);
    assert_eval_eq("TrigExpand[-5 Sin[4x]^2-5 Cos[4x]^2]", "-5", 0);
}

void test_trigexpand_hyperbolic_pythagorean(void) {
    /* Hyperbolic Pythagorean identity collapses for any integer multiple;
     * Factor produces (Cosh+Sinh)^n (Cosh-Sinh)^n, which rewrites to 1. */
    assert_eval_eq("TrigExpand[Cosh[2x]^2-Sinh[2x]^2]", "1", 0);
    assert_eval_eq("TrigExpand[Cosh[3x]^2-Sinh[3x]^2]", "1", 0);
    assert_eval_eq("TrigExpand[Cosh[4x]^2-Sinh[4x]^2]", "1", 0);
    /* Sign-flipped form collapses to -1. */
    assert_eval_eq("TrigExpand[Sinh[2x]^2-Cosh[2x]^2]", "-1", 0);
    assert_eval_eq("TrigExpand[Sinh[4x]^2-Cosh[4x]^2]", "-1", 0);
    /* Even Cosh folds through parity, so negative multiples land on the same
     * collapse. */
    assert_eval_eq("TrigExpand[Cosh[-4x]^2-Sinh[-4x]^2]", "1", 0);
    /* Scalar multiples pass through. */
    assert_eval_eq("TrigExpand[7 Cosh[3x]^2-7 Sinh[3x]^2]", "7", 0);
}

void test_trigexpand_tan_cot_sec_csc(void) {
    /* Tan[2x] = Sin[2x]/Cos[2x] = (2 Cos[x] Sin[x])/(Cos[x]^2 - Sin[x]^2) */
    assert_eval_eq("TrigExpand[Tan[2x]]",
                   "2 Cos[x] Sin[x] Sec[2 x]", 0);
    /* Unchanged on simple arguments. */
    assert_eval_eq("TrigExpand[Tan[x]]", "Tan[x]", 0);
    assert_eval_eq("TrigExpand[Cot[x]]", "Cot[x]", 0);
    assert_eval_eq("TrigExpand[Sec[x]]", "Sec[x]", 0);
    assert_eval_eq("TrigExpand[Csc[x]]", "Csc[x]", 0);
}

void test_trigexpand_hyperbolic_reciprocals(void) {
    /* Tanh[x+y] splits into two terms over a common denominator. */
    assert_eval_eq("TrigExpand[Tanh[x+y]]",
        "Cosh[x] Sinh[y] Sech[x + y] + Sinh[x] Cosh[y] Sech[x + y]", 0);
    assert_eval_eq("TrigExpand[Tanh[x]]", "Tanh[x]", 0);
    assert_eval_eq("TrigExpand[Coth[x]]", "Coth[x]", 0);
    assert_eval_eq("TrigExpand[Sech[x]]", "Sech[x]", 0);
    assert_eval_eq("TrigExpand[Csch[x]]", "Csch[x]", 0);
}

void test_trigexpand_atoms(void) {
    /* Non-trig expressions pass through unchanged. */
    assert_eval_eq("TrigExpand[5]", "5", 0);
    assert_eval_eq("TrigExpand[x]", "x", 0);
    assert_eval_eq("TrigExpand[x+y]", "x + y", 0);
    assert_eval_eq("TrigExpand[x y]", "x y", 0);
}

void test_trigexpand_threads_over_list(void) {
    assert_eval_eq("TrigExpand[{Sin[2x], Cos[2x]}]",
                   "{2 Cos[x] Sin[x], Cos[x]^2 - Sin[x]^2}", 0);
    assert_eval_eq("TrigExpand[{Tan[2x], Sinh[x+y]}]",
        "{2 Cos[x] Sin[x] Sec[2 x], Cosh[x] Sinh[y] + Sinh[x] Cosh[y]}", 0);
}

void test_trigexpand_threads_over_equations(void) {
    assert_eval_eq("TrigExpand[Sin[x+y] == 0]",
                   "Cos[x] Sin[y] + Sin[x] Cos[y] == 0", 0);
    assert_eval_eq("TrigExpand[Cos[2x] > 1/2]",
                   "Cos[x]^2 - Sin[x]^2 > 1/2", 0);
}

void test_trigexpand_threads_over_inequalities(void) {
    assert_eval_eq("TrigExpand[1<Cos[x+y]<2]",
                   "1 < Cos[x] Cos[y] - Sin[x] Sin[y] < 2", 0);
}

void test_trigexpand_threads_over_logic(void) {
    assert_eval_eq("TrigExpand[Sin[2x]==0 && Cos[x+y]==1]",
                   "2 Cos[x] Sin[x] == 0 && Cos[x] Cos[y] - Sin[x] Sin[y] == 1", 0);
    assert_eval_eq("TrigExpand[Not[Sin[2x]==0]]",
                   "Not[2 Cos[x] Sin[x] == 0]", 0);
}

void test_trigexpand_nested(void) {
    /* Non-trig outer heads keep their structure; inner trig gets expanded. */
    assert_eval_eq("TrigExpand[f[Sin[2x], Cos[2x]]]",
                   "f[2 Cos[x] Sin[x], Cos[x]^2 - Sin[x]^2]", 0);
}

void test_trigexpand_numeric(void) {
    /* Sin[0] and Cos[0] get folded at evaluation time. */
    assert_eval_eq("TrigExpand[Sin[0]]", "0", 0);
    assert_eval_eq("TrigExpand[Cos[0]]", "1", 0);
}

void test_trigtoexp_sanity(void) {
    /* Ensure TrigToExp / ExpToTrig still work after the module move. */
    assert_eval_eq("TrigToExp[Cos[x]]", "1/2 E^(-I x) + 1/2 E^(I x)", 0);
    assert_eval_eq("ExpToTrig[Exp[I x]]", "Cos[x] + I Sin[x]", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_trigexpand_double_angle);
    TEST(test_trigexpand_triple_angle);
    TEST(test_trigexpand_higher_angle);
    TEST(test_trigexpand_negative_integer_multiple);
    TEST(test_trigexpand_angle_addition);
    TEST(test_trigexpand_multi_term_sum);
    TEST(test_trigexpand_mixed_sum);
    TEST(test_trigexpand_pythagorean);
    TEST(test_trigexpand_hyperbolic_pythagorean);
    TEST(test_trigexpand_tan_cot_sec_csc);
    TEST(test_trigexpand_hyperbolic_reciprocals);
    TEST(test_trigexpand_atoms);
    TEST(test_trigexpand_threads_over_list);
    TEST(test_trigexpand_threads_over_equations);
    TEST(test_trigexpand_threads_over_inequalities);
    TEST(test_trigexpand_threads_over_logic);
    TEST(test_trigexpand_nested);
    TEST(test_trigexpand_numeric);
    TEST(test_trigtoexp_sanity);

    printf("All TrigExpand tests passed!\n");
    return 0;
}
