/* test_parametric_solve.c
 *
 * Kernel gaps behind FeynCalc's FCLoopFindMomentumShifts and
 * FCLoopFindTopologyMappings (48 tests).  Expected values computed by the
 * reference kernel.
 *
 *  1. Solve[eqns] and Solve[eqns, {}] solve for every variable the system
 *     mentions (f[a] is a variable, Sin[x] contributes x); a system without
 *     variables is {{}} or {}.
 *  2. GroebnerBasis[..., CoefficientDomain -> RationalFunctions] reduced the
 *     ring basis with a pseudo-division that scaled only half of each
 *     polynomial, leaving a parameter-only "generator" (a false unit ideal);
 *     the reduction now runs over Q(params) on FLINT rational functions, the
 *     ring basis uses a main >> params block order, and each generator is
 *     divided by its parameter content.
 *  3. Solve's nonlinear-system path declined above three parameters; the
 *     cap is gone, and back-substituted roots are cancelled.
 *  4. A univariate polynomial with symbolic coefficients is factored before
 *     the radical formulas, at any degree (not only quadratics).
 *  5. PossibleZeroQ[f[a]] was True: an inert call was never sampled.
 *  6. AssociationThread[keys -> scalar] gives the scalar to every key.
 */

#include "test_utils.h"
#include "core.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <time.h>

static void test_solve_implicit_variables(void) {
    assert_eval_eq("Solve[{}, {}]", "{{}}", 0);
    assert_eval_eq("Solve[{}, x]", "{{}}", 0);
    assert_eval_eq("Solve[{x == 1}, {}]", "{{x -> 1}}", 0);
    assert_eval_eq("Solve[x == 1]", "{{x -> 1}}", 0);
    assert_eval_eq("Solve[{x + y == 1, x - y == 0}, {}]", "{{x -> 1/2, y -> 1/2}}", 0);
    assert_eval_eq("Solve[{f[a] + b == 1, b == 2}, {}]", "{{b -> 2, f[a] -> -1}}", 0);
    assert_eval_eq("Solve[{x^2 == 1}]", "{{x -> -1}, {x -> 1}}", 0);
    assert_eval_eq("Solve[{2 == 2}, {}]", "{{}}", 0);
    assert_eval_eq("Solve[{False}, {}]", "{}", 0);
    assert_eval_eq("Solve[{x == 1, x == 2}, {}]", "{}", 0);
}

static void test_rational_functions_groebner(void) {
    /* No parameter-only generator: the system has the generic solution
     * x1 = p2, x2 = p1, x3 = p3. */
    assert_eval_eq("Select[GroebnerBasis[Expand[{x3, x1 + x2 - Q, x1 + x2 + x3 - Q, x2, x2 - Q}^2"
                   " - {p3, p1 + p2 - Q, p1 + p2 + p3 - Q, p1, p1 - Q}^2], {x1, x2, x3},"
                   " CoefficientDomain -> RationalFunctions], FreeQ[#, x1 | x2 | x3] &]", "{}", 0);
    /* parameter content divided out (the reference: {-p3 + x3, p1 - x2, -p2 + x1}
     * up to the sign of each generator) */
    assert_eval_eq("GroebnerBasis[Expand[{x3, x1 + x2 - Q, x1 + x2 + x3 - Q, x2, x2 - Q, x2 + x3,"
                   " x2 + x3 - Q, x1 - Q, x1}^2 - {p3, p1 + p2 - Q, p1 + p2 + p3 - Q, p1, p1 - Q,"
                   " p1 + p3, p1 + p3 - Q, p2 - Q, p2}^2], {x1, x2, x3},"
                   " CoefficientDomain -> RationalFunctions]",
                   "{-p3 + x3, -p1 + x2, -p2 + x1}", 0);
}

static void test_parametric_systems(void) {
    assert_eval_eq("Solve[Thread[{x3, x1 + x2 - Q, x1 + x2 + x3 - Q, x2, x2 - Q, x2 + x3, x2 + x3 - Q,"
                   " x1 - Q, x1}^2 == {p3, p1 + p2 - Q, p1 + p2 + p3 - Q, p1, p1 - Q, p1 + p3,"
                   " p1 + p3 - Q, p2 - Q, p2}^2], {x1, x2, x3}]",
                   "{{x1 -> p2, x2 -> p1, x3 -> p3}}", 0);
    assert_eval_eq("Solve[{x3^2 == p3^2, x2^2 == p1^2, x1^2 == p2^2, (x2 - Q)^2 == (p1 - Q)^2},"
                   " {x1, x2, x3}]",
                   "{{x1 -> -p2, x2 -> p1, x3 -> -p3}, {x1 -> -p2, x2 -> p1, x3 -> p3},"
                   " {x1 -> p2, x2 -> p1, x3 -> -p3}, {x1 -> p2, x2 -> p1, x3 -> p3}}", 0);
    clock_t t0 = clock();
    assert_eval_eq("Solve[{(a + c)^2 == (l1 - l2)^2, (a - b)^2 == (l1 - q2)^2, (b + c)^2 == (l2 - q2)^2,"
                   " (b - d)^2 == (l2 + q1)^2, b^2 == l2^2}, {a, b, c, d}]",
                   "{{a -> l1 - l2 - q2, b -> -l2, c -> q2, d -> -2 l2 - q1},"
                   " {a -> l1 - l2 - q2, b -> -l2, c -> q2, d -> q1},"
                   " {a -> l1 + l2 - q2, b -> l2, c -> -2 l2 + q2, d -> -q1},"
                   " {a -> l1 + l2 - q2, b -> l2, c -> -2 l2 + q2, d -> 2 l2 + q1},"
                   " {a -> -l1 - l2 + q2, b -> -l2, c -> 2 l2 - q2, d -> -2 l2 - q1},"
                   " {a -> -l1 - l2 + q2, b -> -l2, c -> 2 l2 - q2, d -> q1},"
                   " {a -> -l1 + l2 + q2, b -> l2, c -> -q2, d -> -q1},"
                   " {a -> -l1 + l2 + q2, b -> l2, c -> -q2, d -> 2 l2 + q1}}", 0);
    double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
    ASSERT_MSG(secs < 3.0, "5-equation parametric Solve took %.2f s", secs);
}

static void test_parametric_univariate_factor(void) {
    assert_eval_eq("Solve[d^4 + 4 l2^2 q1^2 + 4 l2 q1^3 + q1^4 + d^2 (-4 l2^2 - 4 l2 q1 - 2 q1^2) == 0, d]",
                   "{{d -> -2 l2 - q1}, {d -> -q1}, {d -> q1}, {d -> 2 l2 + q1}}", 0);
}

static void test_possible_zero_q_inert_calls(void) {
    assert_eval_eq("PossibleZeroQ[f[a]]", "False", 0);
    assert_eval_eq("PossibleZeroQ[f[a, b]]", "False", 0);
    assert_eval_eq("PossibleZeroQ[f[a] - f[a]]", "True", 0);
    assert_eval_eq("PossibleZeroQ[f[a] - f[b]]", "False", 0);
    assert_eval_eq("PossibleZeroQ[Sin[f[a]]^2 + Cos[f[a]]^2 - 1]", "True", 0);
    assert_eval_eq("PossibleZeroQ[C[1] - 1]", "False", 0);
    assert_eval_eq("PossibleZeroQ[Derivative[1][f][x] - Derivative[1][f][x]]", "True", 0);
    assert_eval_eq("PossibleZeroQ[1/(g[l, l] + g[l, n]) - 1/(g[l, l] - g[l, n])]", "False", 0);
}

static void test_association_thread_scalar(void) {
    assert_eval_eq("AssociationThread[{a, b, c} -> True]", "<|a -> True, b -> True, c -> True|>", 0);
    assert_eval_eq("AssociationThread[{a, b}, 1]", "<|a -> 1, b -> 1|>", 0);
    assert_eval_eq("AssociationThread[{} -> x]", "<||>", 0);
    assert_eval_eq("AssociationThread[{a, a} -> 1]", "<|a -> 1|>", 0);
    assert_eval_eq("Quiet[AssociationThread[{a, b} -> {1}]]", "AssociationThread[{a, b} -> {1}]", 0);
}

int main(void) {
    symtab_init();
    core_init();
    TEST(test_solve_implicit_variables);
    TEST(test_rational_functions_groebner);
    TEST(test_parametric_systems);
    TEST(test_parametric_univariate_factor);
    TEST(test_possible_zero_q_inert_calls);
    TEST(test_association_thread_scalar);
    printf("All parametric-Solve tests passed.\n");
    return 0;
}
