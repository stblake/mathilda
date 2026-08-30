/*
 * test_dsolve.c — correctness tests for symbolic DSolve (Phase 1, milestone M0).
 *
 * Covers the shared substrate (parse / verify / fit / assemble, both output
 * forms) and the first-order methods DSolve`Quadrature, DSolve`LinearFirstOrder
 * and DSolve`Separable, plus initial-value fitting.  Solutions are checked by
 * back-substitution (PossibleZeroQ of the residual / of the difference from the
 * known closed form) so printer-ordering does not make the tests brittle.
 */
#include "test_utils.h"
#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Evaluate `input` and assert it reduces to the symbol True. */
static void check_true(const char* input) {
    Expr* e = parse_expression(input);
    ASSERT_MSG(e != NULL, "Failed to parse: %s", input);
    Expr* v = evaluate(e);
    char* got = expr_to_string_fullform(v);
    ASSERT_MSG(strcmp(got, "True") == 0,
               "Expected True for %s\n    got: %s", input, got);
    free(got);
    expr_free(v);
    expr_free(e);
}

/* Evaluate `input` and assert its FullForm equals `expected`. */
static void check_form(const char* input, const char* expected) {
    Expr* e = parse_expression(input);
    ASSERT_MSG(e != NULL, "Failed to parse: %s", input);
    Expr* v = evaluate(e);
    char* got = expr_to_string_fullform(v);
    ASSERT_MSG(strcmp(got, expected) == 0,
               "Mismatch for %s\n    expected: %s\n    got:      %s", input, expected, got);
    free(got);
    expr_free(v);
    expr_free(e);
}

/* ---- general solutions satisfy the ODE (residual back-substitutes to 0) ---- */
static void t_linear_first_order_general(void) {
    check_true("PossibleZeroQ[(y'[x] + y[x] - a Sin[x]) /. "
               "DSolve[y'[x] + y[x] == a Sin[x], y, x][[1]]]");
}
static void t_linear_homogeneous_general(void) {
    check_true("PossibleZeroQ[(y'[x] - 3 y[x]) /. "
               "DSolve[y'[x] == 3 y[x], y, x][[1]]]");
}
static void t_separable_general(void) {
    check_true("PossibleZeroQ[(y'[x] + 3 y[x]^2) /. "
               "DSolve[y'[x] == -3 y[x]^2, y, x][[1]]]");
}
static void t_quadrature_second_order_general(void) {
    check_true("PossibleZeroQ[(y''[x] - 7) /. DSolve[y''[x] == 7, y, x][[1]]]");
}
static void t_quadrature_first_order_general(void) {
    check_true("PossibleZeroQ[(y'[x] - Cos[x]) /. DSolve[y'[x] == Cos[x], y, x][[1]]]");
}

/* ---- initial-value problems fit the constants ---- */
static void t_ivp_linear_homogeneous(void) {
    check_true("PossibleZeroQ[(y[x] /. DSolve[{y'[x] == 3 y[x], y[0] == 5}, y[x], x][[1]]) "
               "- 5 E^(3 x)]");
}
static void t_ivp_separable(void) {
    check_true("PossibleZeroQ[(y[x] /. DSolve[{y'[x] == -3 y[x]^2, y[1] == 2}, y[x], x][[1]]) "
               "- 2/(6 x - 5)]");
}
static void t_ivp_linear_sin(void) {
    /* y'+y==a Sin[x], y[0]==0 : solution is 0 at x=0 */
    check_true("PossibleZeroQ[(y[0] /. DSolve[{y'[x] + y[x] == a Sin[x], y[0] == 0}, y, x][[1]])]");
}

/* ---- output forms ---- */
static void t_pure_function_form(void) {
    check_true("MatchQ[DSolve[y'[x] == 3 y[x], y, x][[1, 1]], y -> _Function]");
}
static void t_applied_form(void) {
    check_true("MatchQ[DSolve[y'[x] == 3 y[x], y[x], x][[1, 1]], (y[x] -> _)]");
}

/* ---- named methods (DSolve`Method) match the cascade ---- */
static void t_method_quadrature(void) {
    check_true("PossibleZeroQ[(y'[x] - Cos[x]) /. DSolve`Quadrature[y'[x] == Cos[x], y, x][[1]]]");
}
static void t_method_linear(void) {
    check_true("PossibleZeroQ[(y'[x] + y[x] - x) /. "
               "DSolve`LinearFirstOrder[y'[x] + y[x] == x, y, x][[1]]]");
}
static void t_method_separable(void) {
    check_true("PossibleZeroQ[(y'[x] + 3 y[x]^2) /. "
               "DSolve`Separable[y'[x] == -3 y[x]^2, y, x][[1]]]");
}

/* ---- GeneratedParameters option renames the constant ---- */
static void t_generated_parameters(void) {
    check_true("FreeQ[DSolve[y'[x] == 3 y[x], y, x, GeneratedParameters -> k], C]");
}

/* ---- M1: Bernoulli / Homogeneous / Exact / Clairaut ---- */
static void t_bernoulli(void) {
    check_true("PossibleZeroQ[(y'[x] - y[x]^2 + y[x]) /. "
               "DSolve[y'[x] == y[x]^2 - y[x], y, x][[1]]]");
}
static void t_bernoulli_negative_n(void) {
    /* n = -1; also confirms robust exponent detection over a rational RHS */
    check_true("PossibleZeroQ[(y'[x] - (x^2 + y[x]^2)/(x y[x])) /. "
               "DSolve[y'[x] == (x^2 + y[x]^2)/(x y[x]), y, x][[1]]]");
}
static void t_homogeneous(void) {
    check_true("PossibleZeroQ[(y'[x] - (x - y[x])/(x + y[x])) /. "
               "DSolve[y'[x] == (x - y[x])/(x + y[x]), y, x][[1]]]");
}
static void t_exact(void) {
    check_true("PossibleZeroQ[(2 x y[x] + 1 + x^2 y'[x]) /. "
               "DSolve[2 x y[x] + 1 + x^2 y'[x] == 0, y, x][[1]]]");
}
static void t_exact_value(void) {
    check_true("PossibleZeroQ[(y[x] /. DSolve[2 x y[x] + 1 + x^2 y'[x] == 0, y[x], x][[1]]) "
               "- (C[1] - x)/x^2]");
}
static void t_clairaut_general(void) {
    check_true("PossibleZeroQ[(y[x] - x y'[x] - y'[x]^2) /. "
               "DSolve[y[x] == x y'[x] + y'[x]^2, y, x][[1]]]");
}
static void t_clairaut_singular(void) {
    /* with the option there are two branches; the singular one is -x^2/4 */
    check_form("Length[DSolve[y[x] == x y'[x] + y'[x]^2, y[x], x, "
               "IncludeSingularSolutions -> True]]", "2");
    check_true("PossibleZeroQ[(y[x] - x y'[x] - y'[x]^2) /. "
               "DSolve[y[x] == x y'[x] + y'[x]^2, y, x, IncludeSingularSolutions -> True][[2]]]");
}
static void t_method_bernoulli(void) {
    check_true("PossibleZeroQ[(y'[x] - y[x]^2 + y[x]) /. "
               "DSolve`Bernoulli[y'[x] == y[x]^2 - y[x], y, x][[1]]]");
}
static void t_method_exact(void) {
    check_true("PossibleZeroQ[(2 x y[x] + 1 + x^2 y'[x]) /. "
               "DSolve`Exact[2 x y[x] + 1 + x^2 y'[x] == 0, y, x][[1]]]");
}

/* ---- M2: linear constant-coefficient (homogeneous, inhomogeneous, BVP) ---- */
static void t_cc_inhomogeneous(void) {
    check_true("PossibleZeroQ[(y''[x] + 4 y[x] - 7) /. DSolve[y''[x] + 4 y[x] == 7, y, x][[1]]]");
}
static void t_cc_real_roots(void) {
    check_true("PossibleZeroQ[(y''[x] - 4 y[x]) /. DSolve[y''[x] - 4 y[x] == 0, y, x][[1]]]");
}
static void t_cc_complex_roots(void) {
    check_true("PossibleZeroQ[(y''[x] + 4 y'[x] + 5 y[x]) /. "
               "DSolve[y''[x] + 4 y'[x] + 5 y[x] == 0, y, x][[1]]]");
}
static void t_cc_third_order(void) {
    check_true("PossibleZeroQ[(y'''[x] + 4 y'[x] - 5 y[x]) /. "
               "DSolve[y'''[x] + 4 y'[x] == 5 y[x], y, x][[1]]]");
}
static void t_cc_repeated_root(void) {
    check_true("PossibleZeroQ[(y''[x] - 2 y'[x] + y[x]) /. "
               "DSolve[y''[x] - 2 y'[x] + y[x] == 0, y, x][[1]]]");
    /* a repeated root must give a second independent solution (x E^x) */
    check_true("Not[FreeQ[DSolve[y''[x] - 2 y'[x] + y[x] == 0, y, x][[1]], C[2]]]");
}
static void t_cc_ivp(void) {
    check_true("PossibleZeroQ[(y[x] /. DSolve[{y''[x] + 4 y[x] == 7, y[0] == 1, y'[0] == 2}, y[x], x][[1]]) "
               "- 1/4 (7 - 3 Cos[2 x] + 4 Sin[2 x])]");
}
static void t_cc_bvp(void) {
    check_true("PossibleZeroQ[(y[x] /. DSolve[{y''[x] + y[x] == 0, y[0] == 0, y[Pi/2] == 1}, y[x], x][[1]]) "
               "- Sin[x]]");
}
static void t_method_constcoeff(void) {
    check_true("PossibleZeroQ[(y''[x] - 4 y[x]) /. "
               "DSolve`LinearConstantCoefficients[y''[x] - 4 y[x] == 0, y, x][[1]]]");
}

/* ---- M3: Euler-Cauchy (equidimensional) ---- */
static void t_euler_complex(void) {
    check_true("PossibleZeroQ[(x^2 y''[x] + 4 x y'[x] + 7 y[x]) /. "
               "DSolve[x^2 y''[x] + 4 x y'[x] + 7 y[x] == 0, y, x][[1]]]");
}
static void t_euler_real(void) {
    check_true("PossibleZeroQ[(x^2 y''[x] - 2 y[x]) /. "
               "DSolve[x^2 y''[x] - 2 y[x] == 0, y, x][[1]]]");
}
static void t_euler_repeated(void) {
    check_true("PossibleZeroQ[(x^2 y''[x] - x y'[x] + y[x]) /. "
               "DSolve[x^2 y''[x] - x y'[x] + y[x] == 0, y, x][[1]]]");
    /* repeated root gives a Log[x] term (second constant present) */
    check_true("Not[FreeQ[DSolve[x^2 y''[x] - x y'[x] + y[x] == 0, y, x][[1]], C[2]]]");
}
static void t_euler_inhomogeneous(void) {
    check_true("PossibleZeroQ[(x^2 y''[x] - 2 y[x] - x^2) /. "
               "DSolve[x^2 y''[x] - 2 y[x] == x^2, y, x][[1]]]");
}
static void t_method_euler(void) {
    check_true("PossibleZeroQ[(x^2 y''[x] - 2 y[x]) /. "
               "DSolve`EulerCauchy[x^2 y''[x] - 2 y[x] == 0, y, x][[1]]]");
}

/* ---- M3: special-function recognizers (Airy / Bessel) ---- */
static void t_airy(void) {
    check_true("PossibleZeroQ[(y''[x] - x y[x]) /. DSolve[y''[x] - x y[x] == 0, y, x][[1]]]");
    check_true("Not[FreeQ[DSolve[y''[x] - x y[x] == 0, y[x], x], AiryAi[x]]]");
}
static void t_bessel(void) {
    /* v^2 = 4 -> v = 2; the residual reduces via Bessel recurrences that zero_test
     * cannot decide, so check the recognizer emitted the right heads structurally */
    check_true("Not[FreeQ[DSolve[x^2 y''[x] + x y'[x] + (x^2 - 4) y[x] == 0, y[x], x], BesselJ[2, x]]]");
    check_true("Not[FreeQ[DSolve[x^2 y''[x] + x y'[x] + (x^2 - 4) y[x] == 0, y[x], x], BesselY[2, x]]]");
}
static void t_bessel_modified(void) {
    check_true("Not[FreeQ[DSolve[x^2 y''[x] + x y'[x] - (x^2 + 1) y[x] == 0, y[x], x], BesselI[1, x]]]");
    check_true("Not[FreeQ[DSolve[x^2 y''[x] + x y'[x] - (x^2 + 1) y[x] == 0, y[x], x], BesselK[1, x]]]");
}
static void t_method_specialform(void) {
    check_true("PossibleZeroQ[(y''[x] - x y[x]) /. "
               "DSolve`SpecialFunctionForm[y''[x] - x y[x] == 0, y, x][[1]]]");
}

/* ---- M4: systems of ODEs ---- */
static void t_sys_decoupled(void) {
    check_true("And @@ (PossibleZeroQ /@ ({y'[x] - x^2 y[x], z'[x] - 5 z[x]} /. "
               "DSolve[{y'[x] == x^2 y[x], z'[x] == 5 z[x]}, {y, z}, x][[1]]))");
}
static void t_sys_real_eigenvalues(void) {
    check_true("And @@ (PossibleZeroQ /@ ({y'[x] - z[x], z'[x] - y[x]} /. "
               "DSolve[{y'[x] == z[x], z'[x] == y[x]}, {y, z}, x][[1]]))");
}
static void t_sys_complex_ivp(void) {
    /* eigenvalues +-i -> real Cos/Sin form, plus initial conditions */
    check_true("And @@ (PossibleZeroQ /@ ({y'[x] - (y[x] - 2 z[x]), z'[x] - (y[x] - z[x])} /. "
               "DSolve[{y'[x] == y[x] - 2 z[x], z'[x] == y[x] - z[x], y[0] == 1, z[0] == 4}, {y, z}, x][[1]]))");
    check_true("PossibleZeroQ[((y[x] /. DSolve[{y'[x] == y[x] - 2 z[x], z'[x] == y[x] - z[x], "
               "y[0] == 1, z[0] == 4}, {y[x], z[x]}, x][[1]]) /. x -> 0) - 1]");
}
static void t_sys_constant_forcing(void) {
    check_true("And @@ (PossibleZeroQ /@ ({y'[x] - z[x], z'[x] - (-y[x] + 1)} /. "
               "DSolve[{y'[x] == z[x], z'[x] == -y[x] + 1}, {y, z}, x][[1]]))");
}

/* ---- M4: reduction of order (2nd-order missing y) ---- */
static void t_reduce_order(void) {
    check_true("PossibleZeroQ[(y''[x] - y'[x]^2) /. DSolve[y''[x] == y'[x]^2, y, x][[1]]]");
    /* two independent constants (order preserved) */
    check_true("Not[FreeQ[DSolve[y''[x] == y'[x]^2, y, x][[1]], C[2]]]");
}

/* ---- M6: first-order linear PDEs (method of characteristics) ----
 * Verified with a concrete arbitrary function (C[1][z_] :> Sin[z]) after
 * reducing the Function application, to avoid a pre-existing evaluator crash on
 * D[Function[{x,y}, ...C[1]...]][x,y]. */
static void t_pde_transport(void) {
    check_true("With[{uc = (u[t,x] /. DSolve[D[u[t,x],t] + c D[u[t,x],x] == 0, u, {t,x}][[1]]) "
               "/. C[1][z_] :> Sin[z]}, PossibleZeroQ[D[uc,t] + c D[uc,x]]]");
}
static void t_pde_forcing(void) {
    check_true("With[{uc = (u[x,y] /. DSolve[3 D[u[x,y],x] + 5 D[u[x,y],y] == x, u, {x,y}][[1]]) "
               "/. C[1][z_] :> Sin[z]}, PossibleZeroQ[3 D[uc,x] + 5 D[uc,y] - x]]");
}
static void t_pde_zeroth_order(void) {
    check_true("With[{uc = (u[x,y] /. DSolve[D[u[x,y],x] + 3 D[u[x,y],y] + u[x,y] == 1, u, {x,y}][[1]]) "
               "/. C[1][z_] :> Sin[z]}, PossibleZeroQ[D[uc,x] + 3 D[uc,y] + uc - 1]]");
}

/* ---- unsupported equations stay symbolic (declined, not wrong) ---- */
static void t_declines_unsupported(void) {
    /* a genuinely unrecognised variable-coefficient 2nd-order ODE stays symbolic */
    check_form("Head[DSolve[y''[x] + Sin[x] y[x] == 0, y[x], x]]", "DSolve");
}

int main(void) {
    symtab_init();
    core_init();

    TEST(t_linear_first_order_general);
    TEST(t_linear_homogeneous_general);
    TEST(t_separable_general);
    TEST(t_quadrature_second_order_general);
    TEST(t_quadrature_first_order_general);
    TEST(t_ivp_linear_homogeneous);
    TEST(t_ivp_separable);
    TEST(t_ivp_linear_sin);
    TEST(t_pure_function_form);
    TEST(t_applied_form);
    TEST(t_method_quadrature);
    TEST(t_method_linear);
    TEST(t_method_separable);
    TEST(t_generated_parameters);
    TEST(t_bernoulli);
    TEST(t_bernoulli_negative_n);
    TEST(t_homogeneous);
    TEST(t_exact);
    TEST(t_exact_value);
    TEST(t_clairaut_general);
    TEST(t_clairaut_singular);
    TEST(t_method_bernoulli);
    TEST(t_method_exact);
    TEST(t_cc_inhomogeneous);
    TEST(t_cc_real_roots);
    TEST(t_cc_complex_roots);
    TEST(t_cc_third_order);
    TEST(t_cc_repeated_root);
    TEST(t_cc_ivp);
    TEST(t_cc_bvp);
    TEST(t_method_constcoeff);
    TEST(t_euler_complex);
    TEST(t_euler_real);
    TEST(t_euler_repeated);
    TEST(t_euler_inhomogeneous);
    TEST(t_method_euler);
    TEST(t_airy);
    TEST(t_bessel);
    TEST(t_bessel_modified);
    TEST(t_method_specialform);
    TEST(t_sys_decoupled);
    TEST(t_sys_real_eigenvalues);
    TEST(t_sys_complex_ivp);
    TEST(t_sys_constant_forcing);
    TEST(t_reduce_order);
    TEST(t_pde_transport);
    TEST(t_pde_forcing);
    TEST(t_pde_zeroth_order);
    TEST(t_declines_unsupported);

    printf("\nAll DSolve tests passed.\n");
    return 0;
}
