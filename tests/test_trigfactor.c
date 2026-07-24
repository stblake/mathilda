#include "test_utils.h"
#include "symtab.h"
#include "core.h"

/*
 * Tests for TrigFactor. TrigFactor runs a rule-driven pipeline -- rewriting
 * reciprocals, combining via Together, running Factor, applying Pythagorean
 * and angle-addition identities, and restoring reciprocals -- so the output
 * follows Mathilda's canonical Orderless/Flat sort. Each expected string here
 * reflects Mathilda's canonical form, which may differ cosmetically from
 * Mathematica's pretty-printed form but is mathematically equivalent.
 */

void test_trigfactor_pythagorean_circular(void) {
    assert_eval_eq("TrigFactor[Sin[x]^2+Cos[x]^2]", "1", 0);
    assert_eval_eq("TrigFactor[Cos[x]^2+Sin[x]^2]", "1", 0);
    /* With arbitrary coefficient. */
    assert_eval_eq("TrigFactor[a Sin[x]^2 + a Cos[x]^2]", "a", 0);
    assert_eval_eq("TrigFactor[3 Sin[x]^2 + 3 Cos[x]^2]", "3", 0);
    /* Trailing term survives the collapse. */
    assert_eval_eq("TrigFactor[a Sin[x]^2 + a Cos[x]^2 + b]", "a + b", 0);
    /* Iterated: two independent Pythagorean pairs. */
    assert_eval_eq(
        "TrigFactor[Sin[x]^2 + Cos[x]^2 + Sin[y]^2 + Cos[y]^2]", "2", 0);
}

void test_trigfactor_pythagorean_hyperbolic(void) {
    assert_eval_eq("TrigFactor[Cosh[x]^2-Sinh[x]^2]", "1", 0);
    assert_eval_eq("TrigFactor[Sinh[x]^2-Cosh[x]^2]", "-1", 0);
    assert_eval_eq("TrigFactor[a Cosh[x]^2 - a Sinh[x]^2]", "a", 0);
    assert_eval_eq("TrigFactor[7 Cosh[x]^2 - 7 Sinh[x]^2]", "7", 0);
}

void test_trigfactor_double_angle_circular(void) {
    assert_eval_eq("TrigFactor[2 Sin[x] Cos[x]]", "Sin[2 x]", 0);
    assert_eval_eq("TrigFactor[-2 Sin[x] Cos[x]]", "-Sin[2 x]", 0);
    assert_eval_eq("TrigFactor[Cos[x]^2 - Sin[x]^2]", "Cos[2 x]", 0);
    assert_eval_eq("TrigFactor[Sin[x]^2 - Cos[x]^2]", "-Cos[2 x]", 0);
}

void test_trigfactor_double_angle_hyperbolic(void) {
    assert_eval_eq("TrigFactor[2 Sinh[x] Cosh[x]]", "Sinh[2 x]", 0);
    assert_eval_eq("TrigFactor[-2 Sinh[x] Cosh[x]]", "-Sinh[2 x]", 0);
    assert_eval_eq("TrigFactor[Cosh[x]^2 + Sinh[x]^2]", "Cosh[2 x]", 0);
}

void test_trigfactor_angle_addition_circular(void) {
    assert_eval_eq("TrigFactor[Sin[a] Cos[b] + Cos[a] Sin[b]]",
                   "Sin[a + b]", 0);
    assert_eval_eq("TrigFactor[Sin[a] Cos[b] - Cos[a] Sin[b]]",
                   "Sin[a - b]", 0);
    assert_eval_eq("TrigFactor[Cos[a] Cos[b] - Sin[a] Sin[b]]",
                   "Cos[a + b]", 0);
    assert_eval_eq("TrigFactor[Cos[a] Cos[b] + Sin[a] Sin[b]]",
                   "Cos[a - b]", 0);
}

void test_trigfactor_angle_addition_hyperbolic(void) {
    assert_eval_eq("TrigFactor[Sinh[a] Cosh[b] + Cosh[a] Sinh[b]]",
                   "Sinh[a + b]", 0);
    assert_eval_eq("TrigFactor[Sinh[a] Cosh[b] - Cosh[a] Sinh[b]]",
                   "Sinh[a - b]", 0);
    assert_eval_eq("TrigFactor[Cosh[a] Cosh[b] + Sinh[a] Sinh[b]]",
                   "Cosh[a + b]", 0);
    assert_eval_eq("TrigFactor[Cosh[a] Cosh[b] - Sinh[a] Sinh[b]]",
                   "Cosh[a - b]", 0);
}

void test_trigfactor_cos_plus_sinsin(void) {
    /* Classic inverse of angle-addition expansion: Cos[x+y] + Sin[x]Sin[y]
     * reduces to Cos[x]Cos[y] after expanding the angle-sum (Path B). */
    assert_eval_eq("TrigFactor[Cos[x+y] + Sin[x] Sin[y]]",
                   "Cos[x] Cos[y]", 0);
    /* Sin[x-y] + 2 Sin[y] Cos[x] = Sin[x]Cos[y] + Cos[x]Sin[y] = Sin[x+y]. */
    assert_eval_eq("TrigFactor[Sin[x-y] + 2 Sin[y] Cos[x]]",
                   "Sin[x + y]", 0);
}

void test_trigfactor_tan_expressions(void) {
    /* Sin[x]^2 + Tan[x]^2 factors as Tan[x]^2 (1 + Cos[x]^2) once the
     * Tan -> Sin/Cos rewrite and Factor pass are applied. */
    assert_eval_eq("TrigFactor[Sin[x]^2 + Tan[x]^2]",
                   "(1 + Cos[x]^2) Tan[x]^2", 0);
    /* Hyperbolic analog. */
    assert_eval_eq("TrigFactor[Sinh[x]^2 + Tanh[x]^2]",
                   "(1 + Cosh[x]^2) Tanh[x]^2", 0);
    /* Reciprocal identities: Tan[x] Cos[x] = Sin[x], Sec[x] Cos[x] = 1. */
    assert_eval_eq("TrigFactor[Tan[x] Cos[x]]", "Sin[x]", 0);
    assert_eval_eq("TrigFactor[Tanh[x] Cosh[x]]", "Sinh[x]", 0);
    assert_eval_eq("TrigFactor[Sec[x] Cos[x]]", "1", 0);
    assert_eval_eq("TrigFactor[Sech[x] Cosh[x]]", "1", 0);
}

void test_trigfactor_angle_sum_tan(void) {
    /* Sin[x+y]^2 + Tan[x+y] factors with Tan[x+y] pulled out; the clean
     * factored form preserves the angle-sum structure. */
    assert_eval_eq("TrigFactor[Sin[x+y]^2 + Tan[x+y]]",
                   "Tan[x + y] (1 + Cos[x + y] Sin[x + y])", 0);
}

void test_trigfactor_cosh_power_factorization(void) {
    /* Cosh[x]^2 - Cosh[x]^4 = -Cosh[x]^2 (Cosh[x]-1)(Cosh[x]+1)
     * = -Cosh[x]^2 Sinh[x]^2 via the hyperbolic Pythagorean identity. */
    assert_eval_eq("TrigFactor[Cosh[x]^2 - Cosh[x]^4]",
                   "-Cosh[x]^2 Sinh[x]^2", 0);
    /* Difference of fourth powers collapses: Cos^4 - Sin^4 = Cos[2x]. */
    assert_eval_eq("TrigFactor[Cos[x]^4 - Sin[x]^4]", "Cos[2 x]", 0);
    /* Perfect square: (Sin^2 + Cos^2)^2 = 1. */
    assert_eval_eq(
        "TrigFactor[Sin[x]^4 + 2 Sin[x]^2 Cos[x]^2 + Cos[x]^4]", "1", 0);
}

void test_trigfactor_reciprocals_alone(void) {
    /* Reciprocal trig survives the round-trip through Sin/Cos form. */
    assert_eval_eq("TrigFactor[Tan[x]]", "Tan[x]", 0);
    assert_eval_eq("TrigFactor[Cot[x]]", "Cot[x]", 0);
    assert_eval_eq("TrigFactor[Sec[x]]", "Sec[x]", 0);
    assert_eval_eq("TrigFactor[Csc[x]]", "Csc[x]", 0);
    assert_eval_eq("TrigFactor[Sec[x]^2]", "Sec[x]^2", 0);
    assert_eval_eq("TrigFactor[Csc[x]^2]", "Csc[x]^2", 0);
    assert_eval_eq("TrigFactor[Tanh[x]]", "Tanh[x]", 0);
    assert_eval_eq("TrigFactor[Coth[x]]", "Coth[x]", 0);
    assert_eval_eq("TrigFactor[Sech[x]]", "Sech[x]", 0);
    assert_eval_eq("TrigFactor[Csch[x]]", "Csch[x]", 0);
}

void test_trigfactor_atoms(void) {
    /* Non-trig expressions pass through unchanged. */
    assert_eval_eq("TrigFactor[5]", "5", 0);
    assert_eval_eq("TrigFactor[x]", "x", 0);
    assert_eval_eq("TrigFactor[x+y]", "x + y", 0);
    assert_eval_eq("TrigFactor[x y]", "x y", 0);
}

void test_trigfactor_numeric(void) {
    /* Numeric trig values fold before TrigFactor runs. */
    assert_eval_eq("TrigFactor[Sin[0]]", "0", 0);
    assert_eval_eq("TrigFactor[Cos[0]]", "1", 0);
    assert_eval_eq("TrigFactor[Sin[Pi/2]]", "1", 0);
}

void test_trigfactor_threads_over_list(void) {
    assert_eval_eq(
        "TrigFactor[{Sin[x]^2+Cos[x]^2, 2 Sin[x] Cos[x]}]",
        "{1, Sin[2 x]}", 0);
    assert_eval_eq(
        "TrigFactor[{Sin[a]Cos[b]+Cos[a]Sin[b], Cosh[x]^2-Sinh[x]^2}]",
        "{Sin[a + b], 1}", 0);
}

void test_trigfactor_threads_over_equations(void) {
    assert_eval_eq("TrigFactor[Sin[x]^2 + Cos[x]^2 == 1]", "True", 0);
    assert_eval_eq("TrigFactor[2 Sin[x] Cos[x] == 0]", "Sin[2 x] == 0", 0);
    assert_eval_eq("TrigFactor[Cos[x]^2 - Sin[x]^2 != 0]",
                   "Cos[2 x] != 0", 0);
}

void test_trigfactor_threads_over_inequalities(void) {
    assert_eval_eq(
        "TrigFactor[Sin[a]Cos[b]+Cos[a]Sin[b] > 0]",
        "Sin[a + b] > 0", 0);
}

void test_trigfactor_threads_over_logic(void) {
    assert_eval_eq(
        "TrigFactor[Sin[a]Cos[b]+Cos[a]Sin[b]==0 && Cos[x]^2-Sin[x]^2==1]",
        "Sin[a + b] == 0 && Cos[2 x] == 1", 0);
    assert_eval_eq(
        "TrigFactor[Not[2 Sin[x] Cos[x] == 0]]",
        "Not[Sin[2 x] == 0]", 0);
}

void test_trigfactor_trigexpand_inverse(void) {
    /* TrigFactor is broadly the inverse of TrigExpand for structural
     * identities. Expanded angle-sums recombine, and expanded double-angles
     * re-fold. */
    assert_eval_eq("TrigFactor[TrigExpand[Sin[x+y]]]", "Sin[x + y]", 0);
    assert_eval_eq("TrigFactor[TrigExpand[Cos[x+y]]]", "Cos[x + y]", 0);
    assert_eval_eq("TrigFactor[TrigExpand[Sinh[x+y]]]", "Sinh[x + y]", 0);
    assert_eval_eq("TrigFactor[TrigExpand[Cosh[x+y]]]", "Cosh[x + y]", 0);
    assert_eval_eq("TrigFactor[TrigExpand[Sin[2x]]]", "Sin[2 x]", 0);
    assert_eval_eq("TrigFactor[TrigExpand[Cos[2x]]]", "Cos[2 x]", 0);
}

void test_trigfactor_mixed_forms(void) {
    /* Mixed circular/hyperbolic factors reduce independently. */
    assert_eval_eq(
        "TrigFactor[(Sin[x]^2 + Cos[x]^2)(Cosh[y]^2 - Sinh[y]^2)]",
        "1", 0);
    /* Two independent angle-addition groups collapse independently. */
    assert_eval_eq(
        "TrigFactor[Sin[a]Cos[b] + Cos[a]Sin[b] + Sin[u]Cos[v] + Cos[u]Sin[v]]",
        "Sin[a + b] + Sin[u + v]", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_trigfactor_pythagorean_circular);
    TEST(test_trigfactor_pythagorean_hyperbolic);
    TEST(test_trigfactor_double_angle_circular);
    TEST(test_trigfactor_double_angle_hyperbolic);
    TEST(test_trigfactor_angle_addition_circular);
    TEST(test_trigfactor_angle_addition_hyperbolic);
    TEST(test_trigfactor_cos_plus_sinsin);
    TEST(test_trigfactor_tan_expressions);
    TEST(test_trigfactor_angle_sum_tan);
    TEST(test_trigfactor_cosh_power_factorization);
    TEST(test_trigfactor_reciprocals_alone);
    TEST(test_trigfactor_atoms);
    TEST(test_trigfactor_numeric);
    TEST(test_trigfactor_threads_over_list);
    TEST(test_trigfactor_threads_over_equations);
    TEST(test_trigfactor_threads_over_inequalities);
    TEST(test_trigfactor_threads_over_logic);
    TEST(test_trigfactor_trigexpand_inverse);
    TEST(test_trigfactor_mixed_forms);

    printf("All TrigFactor tests passed!\n");
    return 0;
}
