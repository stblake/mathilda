/* test_residue.c
 *
 * Tests for the symbolic residue builtin (src/calculus/residue.c):
 *   Residue[f, {z, z0}]  -- coefficient of (z - z0)^-1 in the Laurent expansion.
 *
 * Coverage: simple poles, higher-order poles (zero and nonzero residue),
 * trigonometric poles, unknown-function numerators (Derivative coefficients,
 * exercising the adaptive-order retry), complex and algebraic pole locations,
 * a high-order pole yielding a StieltjesGamma, analytic points (residue 0),
 * branch points (unevaluated), argument-count validation, and numeric
 * cross-checks against NResidue (Cauchy's theorem).
 */

#include "core.h"
#include "test_utils.h"
#include "expr.h"
#include "eval.h"
#include "parse.h"
#include "symtab.h"

#include <stdio.h>
#include <string.h>

static void check_eq(const char* input, const char* expected) {
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* r = evaluate(p);
    char* s = expr_to_string(r);
    if (strcmp(s, expected) != 0) {
        fprintf(stderr, "FAIL: %s\n  expected: %s\n  actual:   %s\n",
                input, expected, s);
    }
    ASSERT_STR_EQ(s, expected);
    free(s);
    expr_free(p);
    expr_free(r);
}

/* Simple poles: residue = coefficient of (z - z0)^-1. */
static void test_simple_poles(void) {
    check_eq("Residue[1/z, {z, 0}]", "1");
    check_eq("Residue[1/(z - 1), {z, 1}]", "1");
    check_eq("Residue[1/(z^3 - 1), {z, 1}]", "1/3");
    /* Residue of a rational at a complex simple pole: 1/(2 z0). */
    check_eq("Residue[1/(z^2 + 1), {z, I}]", "-1/2*I");
    check_eq("Residue[1/(z^2 + 1), {z, -I}]", "1/2*I");
    /* Algebraic (radical) pole location -- Series expands about 2^(1/4). */
    check_eq("Residue[x^3/(x^4 - 2), {x, 2^(1/4)}]", "1/4");
}

/* Higher-order poles. A pure (z - z0)^-n term (n >= 2) has zero residue; a
 * nontrivial numerator gives a nonzero one via the Laurent coefficient. */
static void test_higher_order_poles(void) {
    check_eq("Residue[1/z^2, {z, 0}]", "0");
    check_eq("Residue[1/(z - a)^3, {z, a}]", "0");
    /* Order-2 pole at 0: d/dz[(z+1)/(z-2)] at 0 = -3/4. */
    check_eq("Residue[(z + 1)/(z^2 (z - 2)), {z, 0}]", "-3/4");
    /* Exp[z]/z^3: residue = coefficient of z^2 in Exp[z] = 1/2. */
    check_eq("Residue[Exp[z]/z^3, {z, 0}]", "1/2");
}

/* Poles of transcendental functions, resolved through the series engine. */
static void test_transcendental_poles(void) {
    check_eq("Residue[Cot[z], {z, 0}]", "1");
    check_eq("Residue[Cot[z]^2, {z, 0}]", "0");
    check_eq("Residue[Tan[z], {z, Pi/2}]", "-1");
    check_eq("Residue[1/Sin[z]^5, {z, 0}]", "3/8");
    check_eq("Residue[1/Sin[z]^7, {z, 0}]", "5/16");
    /* High-order pole of Zeta at 1: coefficient carries a StieltjesGamma. */
    check_eq("Residue[Zeta[z]/(z - 1)^10, {z, 1}]",
             "-1/362880 StieltjesGamma[9]");
}

/* Unknown-function numerators. Series carries symbolic Derivative[k][f][z0]
 * coefficients; f[z]/z^(n+1) has residue f^(n)[z0]/n!. Order 0 truncates
 * relative to f's depth, so this exercises the adaptive-order retry. */
static void test_unknown_function(void) {
    check_eq("Residue[f[z]/z^5, {z, 0}]", "1/24 Derivative[4][f][0]");
    check_eq("Residue[g[z]/z^8, {z, 0}]", "1/5040 Derivative[7][g][0]");
    check_eq("Residue[f[z]/z, {z, 0}]", "f[0]");
}

/* Analytic points have no principal part -> residue 0. */
static void test_analytic_points(void) {
    check_eq("Residue[z^2, {z, 0}]", "0");
    check_eq("Residue[Sin[z], {z, 0}]", "0");
    check_eq("Residue[1/(z^2 + 1), {z, 0}]", "0");
    check_eq("Residue[Exp[z], {z, 0}]", "0");
}

/* Branch points: the Laurent expansion is fractional (Puiseux), so the residue
 * is undefined and the call is left unevaluated. */
static void test_branch_points(void) {
    check_eq("Residue[1/Sqrt[z], {z, 0}]", "Residue[1/Sqrt[z], {z, 0}]");
    check_eq("Residue[Sqrt[z]/z^2, {z, 0}]", "Residue[1/z^(3/2), {z, 0}]");
}

/* Argument-count and shape validation: leave the call unevaluated. */
static void test_bad_arguments(void) {
    /* Residue[] and Residue[x] emit Residue::argm on stderr (not checked here)
     * and stay unevaluated. */
    check_eq("Residue[]", "Residue[]");
    check_eq("Residue[x]", "Residue[x]");
    /* A malformed location spec is left alone. */
    check_eq("Residue[1/z, {z}]", "Residue[1/z, {z}]");
    check_eq("Residue[1/z, z]", "Residue[1/z, z]");
}

/* Cauchy's theorem: the symbolic residue must agree with the numerical one. */
static void test_numeric_crosscheck(void) {
    check_eq("N[Abs[Residue[1/Sin[z]^5, {z, 0}] "
             "- NResidue[1/Sin[z]^5, {z, 0}]]] < 1/1000", "True");
    check_eq("N[Abs[Residue[1/Sin[z]^7, {z, 0}] "
             "- NResidue[1/Sin[z]^7, {z, 0}]]] < 1/1000", "True");
    check_eq("N[Abs[Residue[1/(z^2 + 1), {z, I}] "
             "- NResidue[1/(z^2 + 1), {z, I}]]] < 1/1000", "True");
    check_eq("N[Abs[Residue[(z + 1)/(z^2 (z - 2)), {z, 0}] "
             "- NResidue[(z + 1)/(z^2 (z - 2)), {z, 0}]]] < 1/1000", "True");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_simple_poles);
    TEST(test_higher_order_poles);
    TEST(test_transcendental_poles);
    TEST(test_unknown_function);
    TEST(test_analytic_points);
    TEST(test_branch_points);
    TEST(test_bad_arguments);
    TEST(test_numeric_crosscheck);

    printf("All Residue tests passed!\n");
    return 0;
}
