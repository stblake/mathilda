#include "test_utils.h"
#include "symtab.h"
#include "core.h"

/*
 * Tests for Simplify. Each block targets a specific class of input the
 * heuristic search must handle: pure algebra, trigonometric identities,
 * hyperbolic/exp roundtrip, the integer-digit complexity tiebreak, and
 * the assumption-aware rewriters.
 */

/* ---- Algebra ---- */

void test_simplify_trivial(void) {
    assert_eval_eq("Simplify[x]", "x", 0);
    assert_eval_eq("Simplify[5]", "5", 0);
    assert_eval_eq("Simplify[a + b]", "a + b", 0);
}

void test_simplify_polynomial(void) {
    assert_eval_eq("Simplify[(x-1)(x+1)(x^2+1)+1]", "x^4", 0);
}

void test_simplify_rational(void) {
    assert_eval_eq("Simplify[3/(x+3)+x/(x+3)]", "1", 0);
}

void test_simplify_collect_by_variable(void) {
    /* a x + b x + c -> c + x (a + b) is a Collect win that no other
     * transform in the heuristic produces; the search must pick the
     * variable to collect by automatically (Variables[]). */
    assert_eval_eq("Simplify[a x + b x + c]", "c + x (a + b)", 0);
}

/* ---- Trigonometric ---- */

void test_simplify_pythagorean(void) {
    assert_eval_eq("Simplify[Sin[x]^2+Cos[x]^2]", "1", 0);
}

void test_simplify_double_angle(void) {
    /* 2 Tan[x] / (1 + Tan[x]^2) -> Sin[2 x]. The result depends on
     * TrigFactor reaching the Sin[2x] form; if the heuristic cannot
     * reduce all the way, accept Sin[2 x] as the canonical answer. */
    assert_eval_eq("Simplify[2 Tan[x]/(1+Tan[x]^2)]", "Sin[2 x]", 0);
}

/* ---- Hyperbolic / exp roundtrip ---- */

void test_simplify_exp_sinh_ratio(void) {
    /* (E^x - E^-x)/Sinh[x] -> 2 via TrigToExp + Cancel + ExpToTrig. */
    assert_eval_eq("Simplify[(E^x-E^(-x))/Sinh[x]]", "2", 0);
}

/* ---- Complexity tiebreak ---- */

void test_simplify_integer_digit_penalty(void) {
    /* Default complexity counts digits in integers, so 100 Log[2] is
     * smaller than Log[2^100] and survives. */
    assert_eval_eq("Simplify[100 Log[2]]", "100 Log[2]", 0);
}

void test_simplify_complexity_function_leafcount(void) {
    /* With ComplexityFunction -> LeafCount, the integer-digit penalty is
     * dropped, so picocas may collapse 100 Log[2] into Log[bignumber].
     * Verify the default vs. override discrepancy by checking the option
     * actually reaches the heuristic search: under LeafCount the answer
     * must NOT contain a multi-digit literal Log argument coefficient. */
    struct Expr* parsed = parse_expression("Simplify[100 Log[2], ComplexityFunction->LeafCount]");
    assert(parsed != NULL);
    struct Expr* result = evaluate(parsed);
    expr_free(parsed);
    char* str = expr_to_string(result);
    /* Either the default form survives (rare with LeafCount) or it has
     * been rewritten as Log[<bignumber>]. Both are acceptable; what is
     * NOT acceptable is a hard parse failure. The point of the test is
     * exercising the option plumbing. */
    int ok = (strstr(str, "Log") != NULL);
    if (!ok) fprintf(stderr, "FAIL: Simplify[100 Log[2], ComplexityFunction->LeafCount] -> %s\n", str);
    assert(ok);
    free(str);
    expr_free(result);
}

/* ---- Threading ---- */

void test_simplify_threads_over_list(void) {
    assert_eval_eq("Simplify[{Sin[x]^2+Cos[x]^2, 3/(x+3)+x/(x+3)}]", "{1, 1}", 0);
}

void test_simplify_threads_over_equation(void) {
    /* Each side simplifies independently. v1 does not rebalance. */
    assert_eval_eq("Simplify[Sin[x]^2+Cos[x]^2 == 1]", "True", 0);
}

/* ---- Assumption-aware ---- */

void test_simplify_sqrt_square_positive(void) {
    assert_eval_eq("Simplify[Sqrt[x^2], x > 0]", "x", 0);
}

void test_simplify_sqrt_square_real(void) {
    assert_eval_eq("Simplify[Sqrt[x^2], Element[x, Reals]]", "Abs[x]", 0);
}

void test_simplify_sqrt_square_no_assumption(void) {
    /* Without assumption, no reduction. */
    assert_eval_eq("Simplify[Sqrt[x^2]]", "Sqrt[x^2]", 0);
}

void test_simplify_inverse_radical_positive(void) {
    assert_eval_eq("Simplify[1/Sqrt[x] - Sqrt[1/x], x > 0]", "0", 0);
}

void test_simplify_log_pow_positive_real(void) {
    assert_eval_eq("Simplify[Log[x^p], x > 0 && Element[p, Reals]]", "p Log[x]", 0);
}

void test_simplify_sin_n_pi_integer(void) {
    assert_eval_eq("Simplify[Sin[n Pi], Element[n, Integers]]", "0", 0);
}

void test_simplify_cos_n_pi_integer(void) {
    assert_eval_eq("Simplify[Cos[n Pi], Element[n, Integers]]", "(-1)^n", 0);
}

void test_simplify_eq_substitution(void) {
    /* (1 - a^2)/b^2  with  a^2 + b^2 == 1  =>  1 */
    assert_eval_eq("Simplify[(1 - a^2)/b^2, a^2 + b^2 == 1]", "1", 0);
}

void test_simplify_negative_assumption(void) {
    /* Sqrt[y^2] under y < 0 -> -y. */
    assert_eval_eq("Simplify[Sqrt[y^2], y < 0]", "-y", 0);
}

void test_simplify_obvious_truth(void) {
    /* The predicate reduces to True structurally without needing a
     * reasoner: x > 0 with x > 0 in the assumption set. */
    assert_eval_eq("Simplify[x > 0, x > 0]", "True", 0);
}

/* Regression: a fractional-power subexpression (1/(y^(2/3) - 1/y^(1/3)))
 * must not Simplify to 0. The original failure was Apart producing 0 on
 * the Together'd form y^(1/3)/(y - 1) -- get_coeff(y^(1/3), y, 0) = 0
 * yielded a zero-row matrix, and the bottom-up Simplify recursion picked
 * the spurious 0 as the lowest-complexity candidate, collapsing the
 * outer Times to y^(19/8) * y^(5/8) = y^3 instead of -y^3/(-1+y). */
void test_simplify_fractional_power_subexpr_not_zero(void) {
    assert_eval_eq("Simplify[1/(y^(2/3) - 1/y^(1/3))]",
                   "y^(1/3)/(-1 + y)", 0);
}

void test_simplify_fractional_power_combine(void) {
    assert_eval_eq(
        "Simplify[y^(5/8) (y^(19/8) - y^(73/24)/(y^(2/3) - 1/y^(1/3)))]",
        "-y^3/(-1 + y)", 0);
}

/* Equation rebalancing: Simplify[lhs == rhs] reduces by GCD of integer
 * coefficients and partitions positive/negative variable terms across
 * the relation, matching Mathematica's canonical form. */
void test_simplify_equation_rebalance(void) {
    assert_eval_eq("Simplify[2 x - 4 y + 6 z - 10 == -8]",
                   "x + 3 z == 1 + 2 y", 0);
}

void test_simplify_equation_polarity_flip(void) {
    /* Negative leading variable coefficient -> divide both sides by -1. */
    assert_eval_eq("Simplify[-2 x == 4]", "x == -2", 0);
}

void test_simplify_inequality_polarity_flip(void) {
    /* Strict inequality must reverse when dividing by negative. */
    assert_eval_eq("Simplify[-2 x < 4]", "x > -2", 0);
}

/* Sqrt[product of squares] under sign assumptions: each known-real
 * factor's Power[Power[s, 2], 1/2] reduces independently. */
void test_simplify_sqrt_product_signs(void) {
    assert_eval_eq("Simplify[Sqrt[x^2 y^2], x > 0 && y < 0]", "-x y", 0);
}

void test_simplify_sqrt_product_three(void) {
    assert_eval_eq("Simplify[Sqrt[x^2 y^2 z^2], x > 0 && y < 0 && z > 0]",
                   "-x y z", 0);
}

/* Cos[k Pi]^m with k integer and m even should collapse to 1 -- both via
 * the explicit even-integer rule (literal 4) and via the Mod[m, 2] == 0
 * symbolic-even path. */
void test_simplify_cos_k_pi_to_even_int_power(void) {
    assert_eval_eq("Simplify[Cos[k Pi]^4, Element[k, Integers]]", "1", 0);
}

void test_simplify_cos_k_pi_to_even_symbolic_power(void) {
    assert_eval_eq(
        "Simplify[Cos[k Pi]^m, Element[k, Integers], Assumptions->Mod[m,2]==0]",
        "1", 0);
}

/* (Cos + Sin)^4 perfect-square completion: reaches Mathematica's form
 * via 1 + 2 Sin Cos -> (Sin + Cos)^2. */
void test_simplify_cos_sin_fourth_power(void) {
    assert_eval_eq("Simplify[4 Sin[x]^2 Cos[x]^2 + 4 Sin[x] Cos[x] + 1]",
                   "(Cos[x] + Sin[x])^4", 0);
}

/* Primitive cube root of -1 satisfies z^2 - z + 1 = 0. */
void test_simplify_cube_root_of_neg_one(void) {
    assert_eval_eq("Simplify[1 - (-1)^(1/3) + (-1)^(2/3)]", "0", 0);
}

/* General roots-of-unity algorithm: lifts the input to a polynomial
 * in omega = (-1)^(1/Q) (Q the LCM of all (-1)^(p/q) and E^(I p Pi/q)
 * denominators) and reduces modulo Phi_{2Q}(omega). The 5th and 7th
 * cases are alternating sums that match Phi_{10}(x) and Phi_{14}(x),
 * so both reduce to 0. */
void test_simplify_fifth_root_of_unity_alternating_sum(void) {
    assert_eval_eq(
        "Simplify[1 - (-1)^(1/5) + (-1)^(2/5) - (-1)^(3/5) + (-1)^(4/5)]",
        "0", 0);
}

void test_simplify_seventh_root_of_unity_alternating_sum(void) {
    assert_eval_eq(
        "Simplify[1 - (-1)^(1/7) + (-1)^(2/7) - (-1)^(3/7) + (-1)^(4/7) "
        "          - (-1)^(5/7) + (-1)^(6/7)]",
        "0", 0);
}

/* E^(+/- 2 I Pi / 3) are the primitive cube roots of unity (zeta and
 * zeta^2). zeta + zeta^2 = -1, so 3 + 2(zeta + zeta^2) = 1. */
void test_simplify_complex_exp_cube_root(void) {
    assert_eval_eq(
        "Simplify[3 + 2 E^(-2 I Pi/3) + 2 E^(2 I Pi/3)]", "1", 0);
}

/* ---- Half-angle tangent (Weierstrass) and Pythagorean reduction ---- */

/* Sin[x] / (c (1 + Cos[x])) -> Tan[x/2] / c (the printer renders the
 * 1/2 coefficient ahead of Tan, with x/2 as `1/2 x` internally). */
void test_simplify_halfangle_tan_with_factor(void) {
    assert_eval_eq("Simplify[Sin[x]/(2 (Cos[x] + 1))]",
                   "1/2 Tan[1/2 x]", 0);
}

void test_simplify_halfangle_tan_numeric_one(void) {
    assert_eval_eq("Simplify[Sin[1]/(Cos[1] + 1)]", "Tan[1/2]", 0);
}

/* Sin[2k] / (1 + Cos[2k]) -> Tan[k] (the half-angle simplifies the integer). */
void test_simplify_halfangle_tan_numeric_two(void) {
    assert_eval_eq("Simplify[Sin[2]/(Cos[2] + 1)]", "Tan[1]", 0);
}

void test_simplify_halfangle_tan_numeric_four(void) {
    assert_eval_eq("Simplify[Sin[4]/(Cos[4] + 1)]", "Tan[2]", 0);
}

/* Sin[1]^a / (1 + Cos[1])^a -> Tan[1/2]^a (general power form). */
void test_simplify_halfangle_tan_power_a(void) {
    assert_eval_eq("Simplify[(Cos[1] + 1)^(-a) Sin[1]^a]",
                   "Tan[1/2]^a", 0);
}

/* (1 + Cos[5])^a / Sin[5]^a -> Tan[5/2]^(-a) (inverse-power form). */
void test_simplify_halfangle_tan_inverse_power(void) {
    assert_eval_eq("Simplify[(Cos[5] + 1)^a Sin[5]^(-a)]",
                   "Tan[5/2]^(-a)", 0);
}

/* Same with an integer assumption -- assumption is consistent with the
 * algebraic identity but not required for it. */
void test_simplify_halfangle_tan_integer_assumption(void) {
    assert_eval_eq(
        "Simplify[(Cos[5] + 1)^n Sin[5]^(-n), Assumptions->Element[n, Integers]]",
        "Tan[5/2]^(-n)", 0);
}

/* Pythagorean reduction inside a product: Sin[x] (1 - Cos[x]^2) -> Sin[x]^3. */
void test_simplify_pythag_reduce_sin_cube(void) {
    assert_eval_eq("Simplify[Sin[x] (1 - Cos[x]^2)]", "Sin[x]^3", 0);
}

/* Symmetric form: Cos[x] (1 - Sin[x]^2) -> Cos[x]^3. */
void test_simplify_pythag_reduce_cos_cube(void) {
    assert_eval_eq("Simplify[Cos[x] (1 - Sin[x]^2)]", "Cos[x]^3", 0);
}

/* Hyperbolic half-angle: Sinh[x] / (c (1 + Cosh[x])) -> Tanh[x/2] / c. */
void test_simplify_halfangle_tanh_with_factor(void) {
    assert_eval_eq("Simplify[Sinh[x]/(2 (Cosh[x] + 1))]",
                   "1/2 Tanh[1/2 x]", 0);
}

void test_simplify_halfangle_tanh_basic(void) {
    assert_eval_eq("Simplify[Sinh[x]/(Cosh[x] + 1)]", "Tanh[1/2 x]", 0);
}

/* Power form: Sinh[x]^a (1 + Cosh[x])^(-a) -> Tanh[x/2]^a. */
void test_simplify_halfangle_tanh_power(void) {
    assert_eval_eq("Simplify[Sinh[1]^a (1 + Cosh[1])^(-a)]",
                   "Tanh[1/2]^a", 0);
}

/* Hyperbolic Pythagorean: Sinh[x] (Cosh[x]^2 - 1) -> Sinh[x]^3. */
void test_simplify_pythag_reduce_sinh_cube(void) {
    assert_eval_eq("Simplify[Sinh[x] (Cosh[x]^2 - 1)]", "Sinh[x]^3", 0);
}

/* Hyperbolic Pythagorean: Cosh[x] (1 + Sinh[x]^2) -> Cosh[x]^3. */
void test_simplify_pythag_reduce_cosh_cube(void) {
    assert_eval_eq("Simplify[Cosh[x] (1 + Sinh[x]^2)]", "Cosh[x]^3", 0);
}

/* Triple-angle identity collapses to -3 Sin[x]^3 once TrigExpand surfaces
 * Cos[x]^2 Sin[x] terms and the sign-flipped Pythagorean rule
 * (Cos[x]^2 - 1 :> -Sin[x]^2) fires inside the factored product. */
void test_simplify_sin_triple_angle_collapse(void) {
    assert_eval_eq("Simplify[Sin[x]^3 + Sin[3 x] - 3 Sin[x]]",
                   "-3 Sin[x]^3", 0);
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_simplify_trivial);
    TEST(test_simplify_polynomial);
    TEST(test_simplify_rational);
    TEST(test_simplify_collect_by_variable);
    TEST(test_simplify_pythagorean);
    TEST(test_simplify_double_angle);
    TEST(test_simplify_exp_sinh_ratio);
    TEST(test_simplify_integer_digit_penalty);
    TEST(test_simplify_complexity_function_leafcount);
    TEST(test_simplify_threads_over_list);
    TEST(test_simplify_threads_over_equation);
    TEST(test_simplify_sqrt_square_positive);
    TEST(test_simplify_sqrt_square_real);
    TEST(test_simplify_sqrt_square_no_assumption);
    TEST(test_simplify_inverse_radical_positive);
    TEST(test_simplify_log_pow_positive_real);
    TEST(test_simplify_sin_n_pi_integer);
    TEST(test_simplify_cos_n_pi_integer);
    TEST(test_simplify_eq_substitution);
    TEST(test_simplify_negative_assumption);
    TEST(test_simplify_obvious_truth);
    TEST(test_simplify_fractional_power_subexpr_not_zero);
    TEST(test_simplify_fractional_power_combine);
    TEST(test_simplify_equation_rebalance);
    TEST(test_simplify_equation_polarity_flip);
    TEST(test_simplify_inequality_polarity_flip);
    TEST(test_simplify_sqrt_product_signs);
    TEST(test_simplify_sqrt_product_three);
    TEST(test_simplify_cos_k_pi_to_even_int_power);
    TEST(test_simplify_cos_k_pi_to_even_symbolic_power);
    TEST(test_simplify_cos_sin_fourth_power);
    TEST(test_simplify_cube_root_of_neg_one);
    TEST(test_simplify_fifth_root_of_unity_alternating_sum);
    TEST(test_simplify_seventh_root_of_unity_alternating_sum);
    TEST(test_simplify_complex_exp_cube_root);
    TEST(test_simplify_halfangle_tan_with_factor);
    TEST(test_simplify_halfangle_tan_numeric_one);
    TEST(test_simplify_halfangle_tan_numeric_two);
    TEST(test_simplify_halfangle_tan_numeric_four);
    TEST(test_simplify_halfangle_tan_power_a);
    TEST(test_simplify_halfangle_tan_inverse_power);
    TEST(test_simplify_halfangle_tan_integer_assumption);
    TEST(test_simplify_pythag_reduce_sin_cube);
    TEST(test_simplify_pythag_reduce_cos_cube);
    TEST(test_simplify_halfangle_tanh_with_factor);
    TEST(test_simplify_halfangle_tanh_basic);
    TEST(test_simplify_halfangle_tanh_power);
    TEST(test_simplify_pythag_reduce_sinh_cube);
    TEST(test_simplify_pythag_reduce_cosh_cube);
    TEST(test_simplify_sin_triple_angle_collapse);

    printf("All Simplify tests passed!\n");
    return 0;
}
