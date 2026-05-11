#include "test_utils.h"
#include "symtab.h"
#include "core.h"

/*
 * Tests for TrigReduce. TrigReduce is the inverse-of-TrigExpand
 * direction: products and integer powers of single-argument trig calls
 * are rewritten into single trig calls of compound arguments. The
 * pipeline is rule-based (product-to-sum + power reduction, with a
 * Together / collapse / from-sincos cleanup), so the printed result
 * follows picocas's canonical Orderless / Flat sort. Each expected
 * string here reflects picocas's canonical print form, which is
 * mathematically equivalent to (but may be cosmetically different
 * from) Mathematica's pretty-printed form.
 */

/* Product-to-sum identities (circular). */
void test_trigreduce_sin_cos_product(void) {
    assert_eval_eq("TrigReduce[2 Sin[a] Cos[b]]",
                   "Sin[a + b] + Sin[a - b]", 0);
}
void test_trigreduce_cos_cos_product(void) {
    assert_eval_eq("TrigReduce[2 Cos[a] Cos[b]]",
                   "Cos[a + b] + Cos[a - b]", 0);
}
void test_trigreduce_sin_sin_product(void) {
    assert_eval_eq("TrigReduce[2 Sin[a] Sin[b]]",
                   "-Cos[a + b] + Cos[a - b]", 0);
}

/* Product-to-sum identities (hyperbolic). */
void test_trigreduce_sinh_cosh_product(void) {
    assert_eval_eq("TrigReduce[2 Sinh[a] Cosh[b]]",
                   "Sinh[a + b] + Sinh[a - b]", 0);
}
void test_trigreduce_cosh_cosh_product(void) {
    assert_eval_eq("TrigReduce[2 Cosh[a] Cosh[b]]",
                   "Cosh[a + b] + Cosh[a - b]", 0);
}
void test_trigreduce_sinh_sinh_product(void) {
    assert_eval_eq("TrigReduce[2 Sinh[a] Sinh[b]]",
                   "Cosh[a + b] - Cosh[a - b]", 0);
}

/* Power reduction (circular degree 2). */
void test_trigreduce_sin_squared(void) {
    assert_eval_eq("TrigReduce[2 Sin[x]^2]",  "1 - Cos[2 x]", 0);
}
void test_trigreduce_cos_squared(void) {
    assert_eval_eq("TrigReduce[2 Cos[x]^2]",  "1 + Cos[2 x]", 0);
}

/* Power reduction (hyperbolic degree 2). */
void test_trigreduce_sinh_squared(void) {
    assert_eval_eq("TrigReduce[2 Sinh[x]^2]", "-1 + Cosh[2 x]", 0);
}
void test_trigreduce_cosh_squared(void) {
    assert_eval_eq("TrigReduce[2 Cosh[x]^2]", "1 + Cosh[2 x]", 0);
}

/* Higher odd powers reduce to a sum of single-arg trig calls. */
void test_trigreduce_sin_cubed(void) {
    /* Sin[x]^3 = (3 Sin[x] - Sin[3 x]) / 4 */
    assert_eval_eq("TrigReduce[Sin[x]^3]",
                   "1/4 (3 Sin[x] - Sin[3 x])", 0);
}
void test_trigreduce_cos_cubed(void) {
    /* Cos[x]^3 = (3 Cos[x] + Cos[3 x]) / 4 */
    assert_eval_eq("TrigReduce[Cos[x]^3]",
                   "1/4 (3 Cos[x] + Cos[3 x])", 0);
}

/* Higher even powers reduce via repeated halving. */
void test_trigreduce_sin_fourth(void) {
    /* Sin[x]^4 = (3 - 4 Cos[2 x] + Cos[4 x]) / 8 */
    assert_eval_eq("TrigReduce[Sin[x]^4]",
                   "1/8 (3 + Cos[4 x] - 4 Cos[2 x])", 0);
}
void test_trigreduce_cos_fourth(void) {
    /* Cos[x]^4 = (3 + 4 Cos[2 x] + Cos[4 x]) / 8 */
    assert_eval_eq("TrigReduce[Cos[x]^4]",
                   "1/8 (3 + Cos[4 x] + 4 Cos[2 x])", 0);
}

/* Pythagorean: Sin[x]^2 + Cos[x]^2 collapses to 1 via the (1 -+ Cos[2x])/2
 * substitutions. */
void test_trigreduce_pythag_circular(void) {
    assert_eval_eq("TrigReduce[Sin[x]^2 + Cos[x]^2]", "1", 0);
    assert_eval_eq("TrigReduce[Cosh[x]^2 - Sinh[x]^2]", "1", 0);
}

/* Reverse double-angle: 2 Sin[x] Cos[x] = Sin[2 x]. */
void test_trigreduce_double_angle(void) {
    assert_eval_eq("TrigReduce[2 Sin[x] Cos[x]]",   "Sin[2 x]", 0);
    assert_eval_eq("TrigReduce[2 Sinh[x] Cosh[x]]", "Sinh[2 x]", 0);
}

/* Reverse angle-addition: products that combine into Sin[a + b], etc. */
void test_trigreduce_angle_addition_circular(void) {
    assert_eval_eq("TrigReduce[Sin[a] Cos[b] + Cos[a] Sin[b]]",
                   "Sin[a + b]", 0);
    assert_eval_eq("TrigReduce[Cos[a] Cos[b] - Sin[a] Sin[b]]",
                   "Cos[a + b]", 0);
    assert_eval_eq("TrigReduce[Sin[a] Cos[b] - Cos[a] Sin[b]]",
                   "Sin[a - b]", 0);
    assert_eval_eq("TrigReduce[Cos[a] Cos[b] + Sin[a] Sin[b]]",
                   "Cos[a - b]", 0);
}
/* Coefficient-aware angle-addition collapses.  Same identities, but
 * with a shared scalar coefficient on both Plus children -- the rules
 * use `c_.` to bind k uniformly and lift it through the rewrite. */
void test_trigreduce_angle_addition_circular_with_coeff(void) {
    assert_eval_eq("TrigReduce[Sqrt[3] Sin[a] Cos[b] + Sqrt[3] Cos[a] Sin[b]]",
                   "Sqrt[3] Sin[a + b]", 0);
    assert_eval_eq("TrigReduce[k Cos[a] Cos[b] - k Sin[a] Sin[b]]",
                   "Cos[a + b] k", 0);
}
void test_trigreduce_angle_addition_hyperbolic(void) {
    assert_eval_eq("TrigReduce[Sinh[a] Cosh[b] + Cosh[a] Sinh[b]]",
                   "Sinh[a + b]", 0);
    assert_eval_eq("TrigReduce[Cosh[a] Cosh[b] + Sinh[a] Sinh[b]]",
                   "Cosh[a + b]", 0);
}

/* Compound input from the user's docstring -- the canonical reverse-
 * angle-addition test case. */
void test_trigreduce_compound_angle_addition(void) {
    assert_eval_eq(
        "TrigReduce[Sin[a] (Cos[b] - Sin[b]) + Cos[a] (Sin[b] + Cos[b])]",
        "Cos[a + b] + Sin[a + b]", 0);
}

/* Tan/Cot/Sec/Csc: rewrite as Sin/Cos, run the pipeline, restore
 * reciprocals on output. */
void test_trigreduce_tan_plus_tan(void) {
    /* (Sin[x]/Cos[x]) + (Sin[y]/Cos[y]) = Sin[x+y]/(Cos[x] Cos[y]). */
    assert_eval_eq("TrigReduce[Tan[x] + Tan[y]]",
                   "Sec[x] Sec[y] Sin[x + y]", 0);
}
void test_trigreduce_tan_plus_cot(void) {
    /* Sin[x]/Cos[x] + Cos[y]/Sin[y] -> (Cos[a-b] form). */
    assert_eval_eq("TrigReduce[Tan[x] + Cot[y]]",
                   "Sec[x] Csc[y] Cos[x - y]", 0);
}
void test_trigreduce_coth_plus_coth(void) {
    assert_eval_eq("TrigReduce[Coth[x] + Coth[y]]",
                   "Csch[x] Csch[y] Sinh[x + y]", 0);
}
void test_trigreduce_tanh_minus_coth(void) {
    assert_eval_eq("TrigReduce[Tanh[x] - Coth[y]]",
                   "Tanh[x] - Coth[y]", 0);
}

/* TrigReduce is idempotent on already-reduced inputs and a no-op on
 * non-trig expressions. */
void test_trigreduce_no_op_atoms(void) {
    assert_eval_eq("TrigReduce[Sin[a + b]]", "Sin[a + b]", 0);
    assert_eval_eq("TrigReduce[Cos[2 x]]",   "Cos[2 x]", 0);
    assert_eval_eq("TrigReduce[Sin[x]]",     "Sin[x]", 0);
    assert_eval_eq("TrigReduce[Tan[x]]",     "Tan[x]", 0);
    assert_eval_eq("TrigReduce[Sec[x]]",     "Sec[x]", 0);
    assert_eval_eq("TrigReduce[5]",          "5", 0);
    assert_eval_eq("TrigReduce[x]",          "x", 0);
    assert_eval_eq("TrigReduce[x^2 + 1]",    "1 + x^2", 0);
}

/* Sin[0] auto-evaluates inside a result, eliminating Sin[a - a] terms. */
void test_trigreduce_sin_cos_same_arg(void) {
    /* Sin[a] Cos[a] = (Sin[2 a] + Sin[0])/2 = Sin[2 a]/2. */
    assert_eval_eq("TrigReduce[Sin[a] Cos[a]]", "1/2 Sin[2 a]", 0);
}

/* Mixed circular / hyperbolic: no identity links them, expression is
 * left as-is up to canonical sort. */
void test_trigreduce_mixed_circular_hyperbolic(void) {
    assert_eval_eq("TrigReduce[Sin[x] Cosh[y]]", "Sin[x] Cosh[y]", 0);
}

/* Threading. ATTR_LISTABLE delivers the List threading; the explicit
 * threading inside builtin_trigreduce_impl handles the rest. */
void test_trigreduce_threads_over_list(void) {
    assert_eval_eq("TrigReduce[{2 Sin[x]^2, 2 Cos[x]^2}]",
                   "{1 - Cos[2 x], 1 + Cos[2 x]}", 0);
    assert_eval_eq(
        "TrigReduce[{Tan[x] + Cot[y], Tanh[x] - Coth[y]}]",
        "{Sec[x] Csc[y] Cos[x - y], Tanh[x] - Coth[y]}", 0);
}
void test_trigreduce_threads_over_equation(void) {
    assert_eval_eq("TrigReduce[2 Sin[x]^2 == 1]",
                   "1 - Cos[2 x] == 1", 0);
}
void test_trigreduce_threads_over_inequalities(void) {
    assert_eval_eq("TrigReduce[2 Cos[x]^2 >= 1]",
                   "1 + Cos[2 x] >= 1", 0);
}
void test_trigreduce_threads_over_logic(void) {
    /* Combined And of two threading-eligible heads. */
    assert_eval_eq(
        "TrigReduce[4 Sin[x]^4 == 1 && 2 Cos[x]^2 >= 1]",
        "1/2 (3 + Cos[4 x] - 4 Cos[2 x]) == 1 && 1 + Cos[2 x] >= 1", 0);
    assert_eval_eq("TrigReduce[Not[2 Sin[x]^2 == 0]]",
                   "Not[1 - Cos[2 x] == 0]", 0);
}

/* Compound argument: Sin[x + y]^2 reduces with the doubled-angle
 * argument carried through the Expand inside the rule body. */
void test_trigreduce_sin_compound_squared(void) {
    assert_eval_eq("TrigReduce[2 Sin[x + y]^2]", "1 - Cos[2 x + 2 y]", 0);
}

/* Product of compound-angle Sin and Cos. */
void test_trigreduce_compound_product(void) {
    /* 2 Sin[x+y] Cos[x-y] = Sin[(x+y)+(x-y)] + Sin[(x+y)-(x-y)]
     *                     = Sin[2 x] + Sin[2 y]. */
    assert_eval_eq("TrigReduce[2 Sin[x + y] Cos[x - y]]",
                   "Sin[2 x] + Sin[2 y]", 0);
}

/* Negative-argument cancellation: Sin[a-b] + Sin[b-a] is mathematically
 * zero. picocas does not auto-canonicalise Sin[Plus[Times[-1, a], b]]
 * into -Sin[a - b], so the collapse pass needs explicit cancellation
 * rules. */
void test_trigreduce_sin_pair_cancellation(void) {
    /* Sin[a-b] + Sin[b-a] -> 0 only when this trick is in place. The
     * test exercises it on a simple two-term input that hits the
     * cancellation rule directly. */
    assert_eval_eq("TrigReduce[Sin[a - b] + Sin[b - a]]", "0", 0);
    assert_eval_eq("TrigReduce[Cos[a - b] + Cos[b - a]]",
                   "2 Cos[a - b]", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_trigreduce_sin_cos_product);
    TEST(test_trigreduce_cos_cos_product);
    TEST(test_trigreduce_sin_sin_product);

    TEST(test_trigreduce_sinh_cosh_product);
    TEST(test_trigreduce_cosh_cosh_product);
    TEST(test_trigreduce_sinh_sinh_product);

    TEST(test_trigreduce_sin_squared);
    TEST(test_trigreduce_cos_squared);
    TEST(test_trigreduce_sinh_squared);
    TEST(test_trigreduce_cosh_squared);

    TEST(test_trigreduce_sin_cubed);
    TEST(test_trigreduce_cos_cubed);
    TEST(test_trigreduce_sin_fourth);
    TEST(test_trigreduce_cos_fourth);

    TEST(test_trigreduce_pythag_circular);
    TEST(test_trigreduce_double_angle);

    TEST(test_trigreduce_angle_addition_circular);
    TEST(test_trigreduce_angle_addition_circular_with_coeff);
    TEST(test_trigreduce_angle_addition_hyperbolic);
    TEST(test_trigreduce_compound_angle_addition);

    TEST(test_trigreduce_tan_plus_tan);
    TEST(test_trigreduce_tan_plus_cot);
    TEST(test_trigreduce_coth_plus_coth);
    TEST(test_trigreduce_tanh_minus_coth);

    TEST(test_trigreduce_no_op_atoms);
    TEST(test_trigreduce_sin_cos_same_arg);
    TEST(test_trigreduce_mixed_circular_hyperbolic);

    TEST(test_trigreduce_threads_over_list);
    TEST(test_trigreduce_threads_over_equation);
    TEST(test_trigreduce_threads_over_inequalities);
    TEST(test_trigreduce_threads_over_logic);

    TEST(test_trigreduce_sin_compound_squared);
    TEST(test_trigreduce_compound_product);

    TEST(test_trigreduce_sin_pair_cancellation);

    printf("All TrigReduce tests passed!\n");
    return 0;
}
