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

/* ---- M11: formal Linear BVP soundness (inconsistent -> {}, not general) ---- */
static void t_bvp_overdetermined(void) {
    /* y[0]==1 forces C1==1; y[Pi]==1 forces -C1==1: inconsistent -> {} (no solution),
     * NOT the silent unfitted general solution it used to return. */
    check_form("DSolve[{y''[x] + y[x] == 0, y[0] == 1, y[Pi] == 1}, y, x]", "List[]");
    check_form("Length[DSolve[{y''[x] + y[x] == 0, y[0] == 1, y[Pi] == 1}, y, x]]", "0");
    /* and the unfitted constants must be gone (would be present in the old behavior) */
    check_true("FreeQ[DSolve[{y''[x] + y[x] == 0, y[0] == 1, y[Pi] == 1}, y, x], C]");
}
static void t_bvp_underdetermined(void) {
    /* y[0]==0 and y[Pi]==0 both give C1==0; C2 stays free -> one branch with a free
     * amplitude constant (not {}, not a decline). */
    check_form("Head[DSolve[{y''[x] + y[x] == 0, y[0] == 0, y[Pi] == 0}, y, x]]", "List");
    check_form("Length[DSolve[{y''[x] + y[x] == 0, y[0] == 0, y[Pi] == 0}, y, x]]", "1");
    check_true("Not[FreeQ[DSolve[{y''[x] + y[x] == 0, y[0] == 0, y[Pi] == 0}, y, x], C[2]]]");
}
static void t_bvp_system_overdetermined(void) {
    /* coupled harmonic system with 3 conditions on a 2-constant general solution,
     * inconsistent (y[0]==1 & y[Pi]==1 clash) -> {}. */
    check_form("DSolve[{y'[x] == z[x], z'[x] == -y[x], y[0] == 1, y[Pi] == 1, z[0] == 0}, {y, z}, x]",
               "List[]");
}
static void t_bvp_undecided_keeps_general(void) {
    /* an IVP whose fit Solve can decide stays fitted; a *no-condition* solve stays
     * general (constants retained) -- guards against an over-eager {} on any solve
     * the fitter cannot decide.  IVP still fits, general solve keeps constants. */
    check_true("Not[FreeQ[DSolve[y''[x] + y[x] == 0, y, x], C[1]]]");
    check_true("PossibleZeroQ[(y[x] /. DSolve[{y'[x] + y[x] == 0, y[0] == 5}, y[x], x][[1]]) - 5 Exp[-x]]");
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
/* Complex-root inhomogeneous Euler: the x = e^t reduction to constant
 * coefficients solves it; the older x-domain variation of parameters hung on
 * the trig-of-Log products (regression for the reported x^2 y'' + y == x^2). */
static void t_euler_inhomogeneous_complex(void) {
    check_true("PossibleZeroQ[(x^2 y''[x] + y[x] - x^2) /. "
               "DSolve[x^2 y''[x] + y[x] == x^2, y, x][[1]]]");
    /* both fundamental (complex-root) modes are present */
    check_true("Not[FreeQ[DSolve[x^2 y''[x] + y[x] == x^2, y, x][[1]], C[2]]]");
    /* non-power forcing g(x) = Log[x] flows through the transformed solve */
    check_true("PossibleZeroQ[(x^2 y''[x] + 4 x y'[x] + 2 y[x] - Log[x]) /. "
               "DSolve[x^2 y''[x] + 4 x y'[x] + 2 y[x] == Log[x], y, x][[1]]]");
    /* an IVP fits both constants through the transformed general solution */
    check_true("PossibleZeroQ[(x^2 y''[x] + y[x] - x^2) /. "
               "DSolve[{x^2 y''[x] + y[x] == x^2, y[1] == 0, y'[1] == 0}, "
               "y, x][[1]]]");
}
/* Reported Euler-Cauchy corpus: forcing of every kind (power, resonant power,
 * double-resonant, 1/x, Log, x^2 e^x -> ExpIntegralEi, Sin[Log], x Log, third
 * order) and a shifted centre.  All checked by residual, since constant
 * absorption / phase combination make the surface form differ from Mathematica
 * while the solution is equivalent. */
static void t_euler_regression_corpus(void) {
    /* In[2]: resonant power r=2 (Log[x] term) */
    check_true("PossibleZeroQ[(x^2 y''[x] - 2 y[x] - x^2) /. "
               "DSolve[x^2 y''[x] - 2 y[x] == x^2, y[x], x][[1]]]");
    /* In[3]: non-resonant power forcing */
    check_true("PossibleZeroQ[(x^2 y''[x] - 2 x y'[x] + 2 y[x] - x^3) /. "
               "DSolve[x^2 y''[x] - 2 x y'[x] + 2 y[x] == x^3, y[x], x][[1]]]");
    /* In[4]: double resonance (repeated root r=2, forcing x^2) -> Log[x]^2 */
    check_true("PossibleZeroQ[(x^2 y''[x] - 3 x y'[x] + 4 y[x] - x^2) /. "
               "DSolve[x^2 y''[x] - 3 x y'[x] + 4 y[x] == x^2, y[x], x][[1]]]");
    check_true("Not[FreeQ[DSolve[x^2 y''[x] - 3 x y'[x] + 4 y[x] == x^2, y[x], x][[1]], "
               "Log[x]^2]]");
    /* In[5]: complex roots +-i, forcing 1/x */
    check_true("PossibleZeroQ[(x^2 y''[x] + x y'[x] + y[x] - 1/x) /. "
               "DSolve[x^2 y''[x] + x y'[x] + y[x] == 1/x, y[x], x][[1]]]");
    /* In[6]: complex roots, Log[x] forcing */
    check_true("PossibleZeroQ[(x^2 y''[x] + y[x] - Log[x]) /. "
               "DSolve[x^2 y''[x] + y[x] == Log[x], y[x], x][[1]]]");
    /* In[7]: real roots, x^2 e^x forcing -> ExpIntegralEi (x-domain VoP) */
    check_true("PossibleZeroQ[(x^2 y''[x] - 2 x y'[x] + 2 y[x] - x^2 Exp[x]) /. "
               "DSolve[x^2 y''[x] - 2 x y'[x] + 2 y[x] == x^2 Exp[x], y[x], x][[1]]]");
    check_true("Not[FreeQ[DSolve[x^2 y''[x] - 2 x y'[x] + 2 y[x] == x^2 Exp[x], y[x], x][[1]], "
               "ExpIntegralEi]]");
    /* In[8]: shifted Euler, centre -1 */
    check_true("PossibleZeroQ[((x+1)^2 y''[x] - 3 (x+1) y'[x] + 3 y[x] - x^2) /. "
               "DSolve[(x+1)^2 y''[x] - 3 (x+1) y'[x] + 3 y[x] == x^2, y[x], x][[1]]]");
    check_true("Not[FreeQ[DSolve[(x+1)^2 y''[x] - 3 (x+1) y'[x] + 3 y[x] == x^2, y[x], x][[1]], "
               "Log[1 + x]]]");
    /* In[9]: real roots, Sin[Log[x]] forcing */
    check_true("PossibleZeroQ[(x^2 y''[x] + 2 x y'[x] - 2 y[x] - Sin[Log[x]]) /. "
               "DSolve[x^2 y''[x] + 2 x y'[x] - 2 y[x] == Sin[Log[x]], y[x], x][[1]]]");
    /* In[10]: real roots, resonant x Log[x] forcing */
    check_true("PossibleZeroQ[(x^2 y''[x] + x y'[x] - y[x] - x Log[x]) /. "
               "DSolve[x^2 y''[x] + x y'[x] - y[x] == x Log[x], y[x], x][[1]]]");
    /* In[11]: third-order Euler, x^4 forcing (three constants) */
    check_true("PossibleZeroQ[(x^3 y'''[x] - 3 x^2 y''[x] + 6 x y'[x] - 6 y[x] - x^4) /. "
               "DSolve[x^3 y'''[x] - 3 x^2 y''[x] + 6 x y'[x] - 6 y[x] == x^4, y[x], x][[1]]]");
    check_true("Not[FreeQ[DSolve[x^3 y'''[x] - 3 x^2 y''[x] + 6 x y'[x] - 6 y[x] == x^4, y[x], x][[1]], "
               "C[3]]]");
}
static void t_method_euler(void) {
    check_true("PossibleZeroQ[(x^2 y''[x] - 2 y[x]) /. "
               "DSolve`EulerCauchy[x^2 y''[x] - 2 y[x] == 0, y, x][[1]]]");
}

/* ---- ExactODE: higher-order exact linear equations (total derivative) ---- */
static void t_method_exactode(void) {
    check_form("Head[DSolve`ExactODE[x y''[x] + y'[x] == 0, y, x]]", "List");
    check_true("PossibleZeroQ[(x y''[x] + y'[x]) /. "
               "DSolve`ExactODE[x y''[x] + y'[x] == 0, y, x][[1]]]");
}
static void t_exactode_more(void) {
    /* another exact 2nd-order homogeneous form */
    check_true("PossibleZeroQ[(x y''[x] + 3 y'[x]) /. "
               "DSolve[x y''[x] + 3 y'[x] == 0, y, x][[1]]]");
    /* inhomogeneous: the forcing is integrated into the first integral */
    check_true("PossibleZeroQ[(x y''[x] + y'[x] - x) /. "
               "DSolve[x y''[x] + y'[x] == x, y, x][[1]]]");
    /* 3rd-order doubly-exact: the recursion reduces order twice (three constants) */
    check_true("PossibleZeroQ[(x y'''[x] + y''[x]) /. "
               "DSolve[x y'''[x] + y''[x] == 0, y, x][[1]]]");
    check_true("Not[FreeQ[DSolve[x y'''[x] + y''[x] == 0, y, x][[1]], C[3]]]");
}
static void t_exactode_declines(void) {
    /* non-exact (Airy): the pinned method stays symbolic */
    check_form("Head[DSolve`ExactODE[y''[x] - x y[x] == 0, y[x], x]]", "DSolve`ExactODE");
    /* first-order is DSolve`Exact — out of scope (the order >= 2 guard) */
    check_form("Head[DSolve`ExactODE[y'[x] + y[x] == 0, y[x], x]]", "DSolve`ExactODE");
}
static void t_exactode_auto(void) {
    /* automatic cascade solves it (via ExactODE, after Euler declines) */
    check_form("Head[DSolve[x y''[x] + y'[x] == 0, y, x]]", "List");
    check_true("PossibleZeroQ[(x y''[x] + y'[x]) /. "
               "DSolve[x y''[x] + y'[x] == 0, y, x][[1]]]");
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
/* Bessel-reducible pure-power potential y'' + A x^m y == 0 -> Sqrt[x] Z(...).
 * y'' - x^4 y == 0 maps to the modified Bessel of order 1/6, argument x^3/3.
 * Verified numerically (the BesselI/BesselK residual is a recurrence identity that
 * zero_test cannot symbolically decide, so we check the emitted heads + a point). */
static void t_bessel_reducible(void) {
    check_form("Head[DSolve[y''[x] - x^4 y[x] == 0, y[x], x]]", "List");
    check_true("Not[FreeQ[DSolve[y''[x] - x^4 y[x] == 0, y[x], x], BesselI[1/6, 1/3 x^3]]]");
    check_true("Not[FreeQ[DSolve[y''[x] - x^4 y[x] == 0, y[x], x], BesselK[1/6, 1/3 x^3]]]");
    /* the oscillatory sign branch (A > 0) uses J/Y */
    check_true("Not[FreeQ[DSolve[y''[x] + x^6 y[x] == 0, y[x], x], BesselJ[1/8, 1/4 x^4]]]");
    check_true("Abs[N[(D[Sqrt[x] BesselI[1/6, x^3/3], {x,2}] - x^4 Sqrt[x] BesselI[1/6, x^3/3]) "
               "/. x -> 13/10]] < 1/1000000");
}
/* Hang guard: a high-degree rational potential must not drive Kovacic into an
 * unbounded Factor/Solve/Integrate; it declines fast to the series fallback.  If
 * the degree gate regresses this loops forever, so reaching the assertion at all
 * is the real test. */
static void t_kovacic_highdegree_no_hang(void) {
    check_form("Head[DSolve[y''[x] + ((x^10 - 1)/(x^12 + 1)) y[x] == 0, y[x], x]]", "List");
}

/* ---- hypergeometric recognizers: Kummer (1F1) / Gauss (2F1) ---- */
static void t_hypergeometric_kummer(void) {
    /* x y'' + (3/2 - x) y' - 2 y == 0 -> Hypergeometric1F1[2, 3/2, x] basis */
    check_true("PossibleZeroQ[(x y''[x] + (3/2 - x) y'[x] - 2 y[x]) /. "
               "DSolve[x y''[x] + (3/2 - x) y'[x] - 2 y[x] == 0, y[x], x][[1]]]");
    check_true("Not[FreeQ[DSolve[x y''[x] + (3/2 - x) y'[x] - 2 y[x] == 0, y[x], x], "
               "HypergeometricPFQ]]");
}
static void t_hypergeometric_gauss(void) {
    /* x(1-x) y'' + (1/2 - 6 x) y' - 6 y == 0 -> Hypergeometric2F1[2, 3, 1/2, x] basis */
    check_true("PossibleZeroQ[(x (1 - x) y''[x] + (1/2 - 6 x) y'[x] - 6 y[x]) /. "
               "DSolve[x (1 - x) y''[x] + (1/2 - 6 x) y'[x] - 6 y[x] == 0, y[x], x][[1]]]");
    check_true("Not[FreeQ[DSolve[x (1 - x) y''[x] + (1/2 - 6 x) y'[x] - 6 y[x] == 0, y[x], x], "
               "HypergeometricPFQ]]");
}
static void t_method_hypergeometric_kummer(void) {
    check_true("PossibleZeroQ[(x y''[x] + (3/2 - x) y'[x] - 2 y[x]) /. "
               "DSolve`SpecialFunctionForm[x y''[x] + (3/2 - x) y'[x] - 2 y[x] == 0, y, x][[1]]]");
}
static void t_hypergeometric_symbolic_a(void) {
    /* symbolic a is permitted because the exponent parameter b = 3/2 is numeric */
    check_true("Not[FreeQ[DSolve[x y''[x] + (3/2 - x) y'[x] - a y[x] == 0, y[x], x], "
               "HypergeometricPFQ]]");
}
static void t_hypergeometric_gauss_symbolic_ab(void) {
    /* symbolic a, b permitted; exponent parameter c = 1/2 is numeric */
    check_true("Not[FreeQ[DSolve[x (1 - x) y''[x] + (1/2 - (a + b + 1) x) y'[x] - a b y[x] == 0, "
               "y[x], x], HypergeometricPFQ]]");
}
static void t_hypergeometric_integer_declines(void) {
    /* integer exponent b = 2 -> SpecialFunctionForm declines (head left unevaluated),
     * so no singular pFq lower-parameter branch is emitted */
    check_true("Head[DSolve`SpecialFunctionForm[x y''[x] + (2 - x) y'[x] - 3 y[x] == 0, y, x]] "
               "=== DSolve`SpecialFunctionForm");
}

/* ---- M17: affine -> Gauss 2F1 for non-canonical regular singular points, and
 *      the Liouville normal-form pre-pass for Bessel/Airy with a y' term. ---- */

/* A Gauss-class equation whose two finite RSPs are {0, 2} (not the canonical
 * {0,1}): mapped affinely and solved as Hypergeometric2F1.  Numeric residual
 * check at an interior point (the 2F1 residual is a contiguous-relation identity
 * zero_test cannot discharge). */
static void t_m17_affine_gauss(void) {
    check_form("Head[DSolve`SpecialFunctionForm[x (x - 2) y''[x] + (2 x - 3) y'[x] + y[x] == 0, y[x], x]]", "List");
    check_true("Not[FreeQ[DSolve`SpecialFunctionForm[x (x - 2) y''[x] + (2 x - 3) y'[x] + y[x] == 0, y[x], x], "
               "HypergeometricPFQ]]");
    check_true("Abs[N[(x (x - 2) y''[x] + (2 x - 3) y'[x] + y[x]) /. "
               "DSolve[x (x - 2) y''[x] + (2 x - 3) y'[x] + y[x] == 0, y, x][[1]] "
               "/. {C[1] -> 13/10, C[2] -> 7/10, x -> 3/5}]] < 1/1000000");
}
/* Gegenbauer at symbolic degree n and symbolic lambda: the F-homotopy pulls the
 * nonzero local exponents at x = +-1 to reach 2F1.  Kovacic owns the integer-n
 * cases; symbolic n is exactly the genuinely-hypergeometric residue. */
static void t_m17_gegenbauer_symbolic(void) {
    check_true("Not[FreeQ[DSolve[(1 - x^2) y''[x] - (2 lam + 1) x y'[x] + n (n + 2 lam) y[x] == 0, y[x], x], "
               "HypergeometricPFQ]]");
    /* numeric residual on the ORIGINAL equation, all parameters instantiated */
    check_true("Abs[N[((1 - x^2) y''[x] - (2 lam + 1) x y'[x] + n (n + 2 lam) y[x]) /. "
               "DSolve[(1 - x^2) y''[x] - (2 lam + 1) x y'[x] + n (n + 2 lam) y[x] == 0, y, x][[1]] "
               "/. {C[1] -> 6/5, C[2] -> 4/5, n -> 23/10, lam -> 7/10, x -> 3/10}]] < 1/100000");
}
/* Associated Legendre at symbolic degree n and order m: declined by the ordinary-
 * Legendre row (mu != 0) and picked up as verifiable 2F1 by the affine/F-homotopy
 * row.  Ordinary Legendre (mu == 0) must still emit LegendreP. */
static void t_m17_associated_legendre(void) {
    check_true("Not[FreeQ[DSolve[(1 - x^2) y''[x] - 2 x y'[x] + (n (n + 1) - m^2/(1 - x^2)) y[x] == 0, y[x], x], "
               "HypergeometricPFQ]]");
    check_true("Not[FreeQ[DSolve[(1 - x^2) y''[x] - 2 x y'[x] + n (n + 1) y[x] == 0, y[x], x], LegendreP]]");
}
/* Liouville normal-form pre-pass: y'' + (2/x) y' + y == 0 has no y'-free form
 * directly, but its normal form is z'' + z == 0 -> spherical Bessel; the recovered
 * solution is Sin[x]/x, Cos[x]/x (verified numerically). */
static void t_m17_normalform_bessel(void) {
    check_form("Head[DSolve`SpecialFunctionForm[y''[x] + (2/x) y'[x] + y[x] == 0, y[x], x]]", "List");
    check_true("Abs[N[(y''[x] + (2/x) y'[x] + y[x]) /. "
               "DSolve`SpecialFunctionForm[y''[x] + (2/x) y'[x] + y[x] == 0, y, x][[1]] "
               "/. {C[1] -> 13/10, C[2] -> 7/10, x -> 6/5}]] < 1/1000000");
}
/* A non-hypergeometric rational equation (only one finite RSP -> Bessel/confluent,
 * deg L != 2): the affine row declines rather than inventing a spurious 2F1. */
static void t_m17_affine_declines_confluent(void) {
    check_true("FreeQ[DSolve`SpecialFunctionForm[x^2 y''[x] + x y'[x] + (x^2 - 4) y[x] == 0, y[x], x], "
               "HypergeometricPFQ]");
}

/* ---- M19: confluent / Whittaker -> 1F1 recogniser (single finite double pole +
 * rank-1 irregular point at infinity, emitted as Exp z^mu Hypergeometric1F1). ---- */

/* corpus 2.1.2-102: x^2 y'' + (c x^2 + b x + a) y == 0 (P == 0), symbolic a,b,c ->
 * verifiable 1F1; numeric residual with the parameters instantiated. */
static void t_m19_whittaker_confluent(void) {
    check_form("Head[DSolve`SpecialFunctionForm[x^2 y''[x] + (c x^2 + b x + a) y[x] == 0, y[x], x]]", "List");
    check_true("Not[FreeQ[DSolve[x^2 y''[x] + (c x^2 + b x + a) y[x] == 0, y[x], x], HypergeometricPFQ]]");
    check_true("Abs[N[(x^2 y''[x] + (c x^2 + b x + a) y[x]) /. "
               "DSolve[x^2 y''[x] + (c x^2 + b x + a) y[x] == 0, y, x][[1]] "
               "/. {C[1] -> 6/5, C[2] -> 4/5, a -> 2/5, b -> 1/3, c -> -3/5, x -> 6/5}]] < 1/1000000");
}

/* Integer 2 mu (mu == 1 here) makes the two 1F1 partners dependent / a lower
 * parameter a non-positive integer -> a correct decline to the series fallback. */
static void t_m19_declines_integer_2mu(void) {
    check_true("FreeQ[DSolve`SpecialFunctionForm[y''[x] + (-1/4 + 1/x + (1/4 - 1)/x^2) y[x] == 0, y[x], x], "
               "HypergeometricPFQ]");
}

/* ---- M20: PolynomialShiftSubstitution (x-dependent-shift substitution ->
 * separable, implicit first integral).  Verified by the implicit-function rule
 * D[G,x] + D[G,y[x]]*RHS == 0 (branch-safe). ---- */

/* corpus 2.1.2-402: y' == x^2/2 + (1+x^2+x^3) Sqrt[x^3 - 6 y]; u = x^3 - 6 y. */
static void t_m20_polyshift_402(void) {
    check_form("Head[DSolve[y'[x] == x^2/2 + (1 + x^2 + x^3) Sqrt[x^3 - 6 y[x]], y, x][[1, 1]]]", "Equal");
    check_true("With[{g = (DSolve[y'[x] == x^2/2 + (1 + x^2 + x^3) Sqrt[x^3 - 6 y[x]], y, x][[1, 1, 1]] "
               "- C[1])}, Abs[N[(D[g, x] /. y'[x] -> (x^2/2 + (1 + x^2 + x^3) Sqrt[x^3 - 6 y[x]])) "
               "/. {x -> 6/5, y[x] -> 1/10}, 20]] < 1/1000000]");
}

/* corpus 2.1.2-371: symbolic parameter a, base Sqrt[a x^4 + 8 y]. */
static void t_m20_polyshift_371(void) {
    check_true("Head[DSolve[y'[x] == -(1/2) Sqrt[a] x^3 (Sqrt[a] + Sqrt[a] x - 2 Sqrt[a x^4 + 8 y[x]])/(1 + x), "
               "y, x][[1, 1]]] === Equal");
}

/* A radical base NONLINEAR in y (x + y^2) is outside the method's class -> not
 * claimed by PolynomialShiftSubstitution. */
static void t_m20_declines_nonlinear_base(void) {
    check_true("Head[DSolve`PolynomialShiftSubstitution[y'[x] == Sqrt[x + y[x]^2], y, x]] =!= List");
}

/* ---- M23: exact ODE with a radical potential -> implicit first integral ---- */
/* corpus 2.2.3-204: an EXACT ODE (M_y == N_x) whose potential
 *   F = 6 x^(3/2) y^(4/3) - 10 x^(6/5) y^(3/2)
 * carries fractional powers of y, so Solve[F == C[1], y] does not terminate.  The
 * explicit Exact entry is gated to a rational-in-y potential (ds_is_rational_in),
 * so this falls through to the implicit entry, which returns F(x, y[x]) == C[1]
 * verbatim (as Maple/Mathematica do).  Verified here by the implicit-function
 * rule: d/dx[F - C[1]] with y' -> -M/N must vanish. */
static void t_m23_exact_radical(void) {
    const char* ode = "9 Sqrt[x] y[x]^(4/3) - 12 x^(1/5) y[x]^(3/2) + "
                      "(8 x^(3/2) y[x]^(1/3) - 15 x^(6/5) Sqrt[y[x]]) y'[x] == 0";
    /* the branch is implicit (an Equal relation), not an explicit Rule */
    check_form("Head[DSolve[9 Sqrt[x] y[x]^(4/3) - 12 x^(1/5) y[x]^(3/2) + "
               "(8 x^(3/2) y[x]^(1/3) - 15 x^(6/5) Sqrt[y[x]]) y'[x] == 0, y, x][[1, 1]]]", "Equal");
    (void)ode;
    check_true("With[{g = DSolve[9 Sqrt[x] y[x]^(4/3) - 12 x^(1/5) y[x]^(3/2) + "
               "(8 x^(3/2) y[x]^(1/3) - 15 x^(6/5) Sqrt[y[x]]) y'[x] == 0, y, x][[1, 1, 1]] - C[1]}, "
               "Abs[N[(D[g, x] /. y'[x] -> "
               "-(9 Sqrt[x] y[x]^(4/3) - 12 x^(1/5) y[x]^(3/2))/"
               "(8 x^(3/2) y[x]^(1/3) - 15 x^(6/5) Sqrt[y[x]])) "
               "/. {x -> 13/10, y[x] -> 7/10}, 20]] < 1/1000000]");
    /* over-restriction guard: a rational-in-y exact ODE still solves EXPLICITLY
     * (an explicit Rule branch), so the gate never demotes an invertible case. */
    check_form("Head[DSolve[2 x y[x] + 1 + x^2 y'[x] == 0, y, x][[1, 1]]]", "Rule");
}

/* ---- M24: trig-power/product forcing via TrigReduce in UndeterminedCoefficients ---- */
/* corpus 2.2.4-326/362/363/365: a constant-coefficient linear ODE whose forcing is a
 * trig POWER or PRODUCT (Sin[x]^2, Cos[x]^3, Sin[3x] Sin[x], x Cos[x]^3) is not itself
 * a UC function, so it previously declined to the (hanging) variation-of-parameters
 * fallback.  UndeterminedCoefficients now TrigReduce-linearises the forcing into a sum
 * of first-harmonic sinusoids, each a UC function, so it solves tidily. */
static void t_m24_trig_power_forcing(void) {
    /* solved (an explicit Rule branch), not declined */
    check_form("Head[DSolve[y''[x]+y'[x]+y[x] == Sin[x]^2, y, x][[1, 1]]]", "Rule");
    /* residual back-substitutes to zero (Sin^2 forcing, and the resonant x Cos^3) */
    check_true("PossibleZeroQ[(y''[x]+y'[x]+y[x] - Sin[x]^2) /. "
               "DSolve[y''[x]+y'[x]+y[x] == Sin[x]^2, y, x][[1]]]");
    check_true("PossibleZeroQ[(y''[x]+y[x] - x Cos[x]^3) /. "
               "DSolve[y''[x]+y[x] == x Cos[x]^3, y, x][[1]]]");
}

/* ---- M24: numeric complex roots concretized so an IVP fits ---- */
/* corpus 2.2.4-312: y'''==y has complex cube roots the const-coeff basis emitted as
 * Re[-(-1)^(1/3)]/Im[-(-1)^(1/3)] (Re/Im do not auto-evaluate on a radical power),
 * which blocked the IVP constant-fit (the general solution leaked C[k]).
 * dsolve_homog_basis now ComplexExpands the real/imag parts of a NUMERIC complex root,
 * so the fit produces a fully-determined concrete solution. */
static void t_m24_complex_cuberoot_ivp(void) {
    /* one explicit branch, fully fitted -> no generated constant leaks */
    check_true("With[{s = DSolve[{y'''[x] == y[x], y[0] == 1, y'[0] == 0, y''[0] == 0}, y, x]}, "
               "MatchQ[s, {{y -> _Function}}] && FreeQ[s, C[_]]]");
    /* satisfies the ODE and the initial value y(0) == 1 */
    check_true("PossibleZeroQ[(y'''[x] - y[x]) /. "
               "DSolve[{y'''[x] == y[x], y[0] == 1, y'[0] == 0, y''[0] == 0}, y, x][[1]]]");
    check_true("With[{s = DSolve[{y'''[x] == y[x], y[0] == 1, y'[0] == 0, y''[0] == 0}, y, x][[1]]}, "
               "Abs[N[(y[0] /. s) - 1, 20]] < 1/1000000]");
}

/* ---- M25: exact ODE -> Erf closed form (verify no longer spins) ---- */
/* corpus 2.2.5-428: y''+x y'+y==0 is exact, reducing to the first-order linear
 * y'+x y==C[2] whose integrating-factor solution carries Erf[-I x/Sqrt[2]].  The
 * closed form was produced correctly, but dsolve_verify_body's zero_test spun on
 * the Gaussian x Erf residual (E^(-x^2/2) . E^(x^2/2) products that never combine
 * defeat the numeric precision ladder, POSSIBLE_ZEROQ_IMPROVEMENTS.md #1), so the
 * whole solve timed out.  The verify now ExpandAll-normalises an Erf/Erfi residual
 * first, so the exact closed form returns.  (PZQ in this test likewise ExpandAlls
 * to sidestep the same standing zero_test limitation.) */
static void t_m25_exact_erf(void) {
    /* returns an explicit closed-form Rule branch (not a timeout / decline) */
    check_form("Head[DSolve[y''[x]+x y'[x]+y[x]==0, y, x][[1, 1]]]", "Rule");
    /* the Erf solution back-substitutes to zero */
    check_true("PossibleZeroQ[ExpandAll[(y''[x]+x y'[x]+y[x]) /. "
               "DSolve[y''[x]+x y'[x]+y[x]==0, y, x][[1]]]]");
}

/* ---- M25: Kovacic must return a fundamental set (independence guard) ---- */
/* corpus 2.2.5-482: 2x y''+(1-2x^2)y'-4x y==0 has ONE Liouvillian solution
 * Sqrt[x] E^(x^2/2); its reduction-of-order second solution is non-elementary.
 * Kovacic's coincident-exponent path collapsed to the rank-deficient
 * (C[1]+C[2]) Sqrt[x] E^(x^2/2) -- it verifies but is not a general solution.
 * The new independence guard rejects a degenerate basis, so the cascade falls
 * through to Frobenius, which returns the correct two-parameter series. */
static void t_m25_kovacic_fundamental_set(void) {
    /* a genuine two-parameter solution: the two basis solutions (coeffs of C[1],
     * C[2]) are linearly independent -- their ratio is non-constant. */
    check_true("With[{b = y[x] /. DSolve[2 x y''[x]+(1-2 x^2) y'[x]-4 x y[x]==0, y, x][[1]]}, "
               "!PossibleZeroQ[D[Normal[D[b, C[1]]]/Normal[D[b, C[2]]], x]]]");
    /* it is a (Frobenius) series solution */
    check_true("With[{s = DSolve[2 x y''[x]+(1-2 x^2) y'[x]-4 x y[x]==0, y, x]}, "
               "MatchQ[s, {{y -> _Function}}] && !FreeQ[s, SeriesData]]");
}

/* ---- M25: transcendental-coefficient Frobenius at a regular singular point ---- */
/* corpus 2.2.5-463/490: x^2 y''+6 Sin[x] y'+6 y==0 and 2x^2 y''+Sin[x] y'-Cos[x] y==0
 * are regular singular at x=0, but forming xP = x*P leaves a removable singularity
 * (6 Sin[x]/x is 6 at 0 but substitutes to 6 Sin[0]/0 = Indeterminate), so the
 * indicial roots came out garbage and Frobenius declined.  frobenius_regsing now
 * Taylor-normalises xP and x^2 Q first, so an analytic transcendental coefficient
 * is handled and a verified Frobenius series is returned. */
static void t_m25_transcendental_frobenius(void) {
    /* both return a Frobenius series (they previously declined) */
    check_true("With[{s = DSolve[x^2 y''[x]+6 Sin[x] y'[x]+6 y[x]==0, y, x]}, "
               "MatchQ[s, {{y -> _Function}}] && !FreeQ[s, SeriesData]]");
    check_true("With[{s = DSolve[2 x^2 y''[x]+Sin[x] y'[x]-Cos[x] y[x]==0, y, x]}, "
               "MatchQ[s, {{y -> _Function}}] && !FreeQ[s, SeriesData]]");
    /* correctness anchor: the 490 series back-substitutes to ~0 near x=0 */
    check_true("With[{b = Normal[y[x] /. DSolve[2 x^2 y''[x]+Sin[x] y'[x]-Cos[x] y[x]==0, y, x][[1]]] "
               "/. {C[1] -> 13/10, C[2] -> 7/10}}, "
               "Abs[N[(2 x^2 D[b, {x,2}] + Sin[x] D[b, x] - Cos[x] b) /. x -> 1/10, 30]] < 1/100000]");
}

/* ---- M26: distribution value rules + DiracDelta sifting + Leibniz FTC ---- */
static void t_m26_distributions(void) {
    /* HeavisideTheta is left-continuous on definite-sign numeric arguments */
    check_form("HeavisideTheta[0]", "0");
    check_form("HeavisideTheta[3]", "1");
    check_form("HeavisideTheta[-2]", "0");
    check_form("HeavisideTheta[3 Pi]", "1");
    check_form("DiracDelta[3]", "0");
    check_form("DiracDelta[-2]", "0");
    /* HeavisideTheta' = DiracDelta */
    check_form("D[HeavisideTheta[t], t]", "DiracDelta[t]");
    /* DiracDelta sifting under a definite integral (interior + endpoint) */
    check_form("Integrate[DiracDelta[t - 1] f[t], {t, 0, 3}]", "f[1]");
    check_true("PossibleZeroQ[Integrate[DiracDelta[t - 5] f[t], {t, 0, 3}]]");   /* outside -> 0 */
    /* variable-limit Leibniz rule: D[Integrate[e,{u,0,t}],t] = e|_{u=t} */
    check_form("D[Integrate[h[u], {u, 0, t}], t]", "h[t]");
    check_true("PossibleZeroQ[D[Integrate[g[u] Sin[t - u], {u, 0, t}], t] "
               "- Integrate[g[u] Cos[t - u], {u, 0, t}]]");
}

/* ---- M26: constant-coefficient ODE with a DiracDelta impulse forcing ---- */
/* corpus 2.2.6-564: x''+4x==DiracDelta[t], x(0)=x'(0)=0 -> the causal impulse
 * response (1/2) Sin[2t] HeavisideTheta[t] (Green's function via the definite
 * variation-of-parameters convolution + the sifting property). */
static void t_m26_impulse_forcing(void) {
    check_form("Head[DSolve[{y''[t]+4 y[t]==DiracDelta[t], y[0]==0, y'[0]==0}, y, t][[1,1]]]",
               "Rule");
    /* value pinned numerically at t=1 (post-impulse): (1/2) Sin[2] */
    check_true("Abs[N[(y[t] /. DSolve[{y''[t]+4 y[t]==DiracDelta[t], y[0]==0, y'[0]==0}, "
               "y, t][[1]]) /. t -> 1] - 1/2 Sin[2]] < 1/1000000");
    /* shifted impulse response starts at t=Pi (0 before, nonzero after) */
    check_true("PossibleZeroQ[(y[t] /. DSolve[{y''[t]+4 y[t]==DiracDelta[t - Pi], y[0]==0, "
               "y'[0]==0}, y, t][[1]]) /. t -> 1]");   /* t=1 < Pi -> 0 */
}

/* ---- M26: constant-coefficient ODE with an ARBITRARY forcing f(t) ---- */
/* corpus 2.2.6-561..563/572..575: x''+..==f(t), zero ICs -> the Duhamel
 * convolution Integrate[G(t,s) f(s), {s,0,t}].  Verified by a concrete probe:
 * with f -> Cos[3 #] the convolution closes and back-substitutes to 0. */
static void t_m26_general_forcing(void) {
    /* an explicit branch carrying the (unevaluated) convolution integral */
    check_true("With[{b = y[t] /. DSolve[{y''[t]+4 y[t]==f[t], y[0]==0, y'[0]==0}, y, t][[1]]}, "
               "MatchQ[b, _] && !FreeQ[b, Integrate] && FreeQ[b, C]]");
    /* probe verify (distinct complex roots): residual reduces to exactly 0 */
    check_true("With[{b = (y[t] /. DSolve[{y''[t]+4 y[t]==f[t], y[0]==0, y'[0]==0}, y, t][[1]]) "
               "/. f -> (Cos[3 #] &)}, PossibleZeroQ[Simplify[D[b, {t,2}] + 4 b - Cos[3 t]]]]");
    /* real-root kernel must not hang and must probe-verify (f -> t) */
    check_true("With[{b = (y[t] /. DSolve[{y''[t]+6 y'[t]+8 y[t]==f[t], y[0]==0, y'[0]==0}, y, t][[1]]) "
               "/. f -> (# &)}, PossibleZeroQ[Simplify[D[b, {t,2}] + 6 D[b, t] + 8 b - t]]]");
}

/* ---- M27: corpus harness now VERIFIES systems (back-substitution) ---- */
/* A constant-coefficient system solves and both equations back-substitute to ~0
 * -- exactly what the M27 corpus harness now checks for a system.  The second
 * case (eigenvalues -10, -100) previously HUNG: dsolve_linsys_tidy Simplify-ed a
 * body that is a sum of exponentials with widely-separated real decay rates, and
 * Simplify's zero-test spins on Exp[-10 t] against Exp[-100 t].  tidy now Expands
 * any exponential body instead of Simplifying it. */
static void t_m27_system_verify(void) {
    check_true("MatchQ[DSolve[{x'[t]==3 x[t]-2 y[t], y'[t]==2 x[t]+y[t]}, {x,y}, t], "
               "{{x -> _Function, y -> _Function}}]");
    check_true("With[{s = DSolve[{x'[t]==-50 x[t]+20 y[t], y'[t]==100 x[t]-60 y[t]}, {x,y}, t][[1]]}, "
               "Max[Abs[N[{(x'[t]-(-50 x[t]+20 y[t])), (y'[t]-(100 x[t]-60 y[t]))} /. s "
               "/. {C[1]->7/10, C[2]->13/10, t->3/10}, 20]]] < 1/1000000]");
}

/* ---- M27: Solve must mint a FRESH periodicity index (no C[1] collision) ---- */
/* corpus 2.2.7-684: y'=2x Sec[y].  DSolve`Separable feeds Solve an equation that
 * already carries the integration constant C[1]; Solve reused C[1] as the 2 Pi k
 * inverse-trig periodicity index, so once its Element[C[1],Integers] constraint
 * was dropped the shared C[1] shifted y by a non-multiple of 2 Pi -- a WRONG
 * answer (masked as UNEVAL by the verifier's leaked->UNFIT rule).  solveinv now
 * SEEDS its mint counter past every C[k] already in the equation, and
 * dsolve_extract_solutions collapses the integer family (Element[C[k],Integers])
 * to its principal branch, so the solution is correct at a CONTINUOUS constant. */
static void t_m27_separable_inverse_constant(void) {
    check_true("With[{s = DSolve[y'[x]==2 x Sec[y[x]], y, x][[1]]}, "
               "Abs[N[(y'[x]-2 x Sec[y[x]]) /. s /. {C[1]->7/10, x->13/10}, 30]] < 1/1000000]");
}

/* ---- M27: the family collapse must NOT eat the integration constant ---- */
/* corpus 2.2.7-695: y'=3x^2(1+y^2), y(0)=1 -> Tan[x^3+Pi/4].  The Tan/ArcTan
 * inversion condition is a RANGE on x^3+C[1] (not Element[_,Integers]); collapsing
 * every C[k] in a condition would zero the integration constant and break the IVP
 * fit.  The collapse is scoped to Element[C[k],Integers] only, so the IC holds. */
static void t_m27_ivp_family_intact(void) {
    check_true("With[{s = DSolve[{y'[x]==3 x^2 (1+y[x]^2), y[0]==1}, y, x][[1]]}, "
               "Abs[N[(y[0] /. s) - 1, 20]] < 1/1000000 && "
               "Abs[N[(y'[x]-3 x^2 (1+y[x]^2)) /. s /. x->2/10, 20]] < 1/1000000]");
    check_true("FreeQ[DSolve[{y'[x]==3 x^2 (1+y[x]^2), y[0]==1}, y, x], C[_]]");
}

/* ---- M28: Bernoulli must decline FAST on a transcendental-in-y RHS ---- */
/* corpus 2.2.8-757: 2 x Sin[y] Cos[y] y' == 4 x^2 + Sin[y]^2 reduces (u = Sin[y]^2)
 * to the linear u' - u/x == 4 x, which Linearizable solves instantly.  But the RHS,
 * solved for y', is a rational function of Sin[y]/Cos[y] -- transcendental in y, not
 * the Bernoulli form A(x) y + B(x) y^n -- and the Bernoulli exponent detector spun on
 * it for seconds (derivative + Cancel + free-of on a trig mess), timing out the whole
 * cascade before Linearizable was reached ($Aborted).  dsolve_bernoulli.c now declines
 * immediately when y appears inside a non-Power function or a power exponent. */
static void t_m28_bernoulli_hang_trig_substitution(void) {
    /* solves explicitly now (the cascade reaches Linearizable instead of $Aborted) */
    check_form("Head[DSolve[2 x Sin[y[x]] Cos[y[x]] y'[x] == 4 x^2 + Sin[y[x]]^2, y, x][[1, 1]]]", "Rule");
    /* the solution back-substitutes to zero (numeric, C[1] and x instantiated) */
    check_true("With[{s = DSolve[2 x Sin[y[x]] Cos[y[x]] y'[x] == 4 x^2 + Sin[y[x]]^2, y, x][[1]] /. C[1] -> 1}, "
               "Abs[N[(2 x Sin[y[x]] Cos[y[x]] y'[x] - (4 x^2 + Sin[y[x]]^2)) /. s /. x -> 13/10, 20]] < 1/1000000]");
    /* forward-generator: the general a-coefficient family solves too (not overfit to a==4) */
    check_form("Head[DSolve[2 x Sin[y[x]] Cos[y[x]] y'[x] == 9 x^2 + Sin[y[x]]^2, y, x][[1, 1]]]", "Rule");
    /* over-restriction guard: a GENUINE Bernoulli (2.2.8-752, F = y - E^(-2x)/(2x) y^3)
     * is algebraic in y, so the new gate never fires and it still solves. */
    check_form("Head[DSolve[2 y'[x] x + y[x]^3 E^(-2 x) == 2 y[x] x, y, x][[1, 1]]]", "Rule");
}

/* ---- M29: Sec-forcing VoP solution verifies via piecewise Floor derivative ---- */
/* corpus 2.2.9-898: y'' + 9y == 2 Sec[3x] is solved by variation of parameters and
 * its solution carries a Floor[...] branch-tracking term.  Before M29 D[Floor[u],x]
 * returned an inert Derivative[1][Floor][u], so the ODE residual never numericized
 * and the corpus harness passed 898 only under the "non-numericizable => trust
 * DSolve" path.  With piecewise rounding-function derivatives (deriv.c) the residual
 * reduces to a genuine numeric ~0. */
static void t_m29_sec_floor_verifies(void) {
    /* the Floor derivative numericizes to 0 off the integer jump set (bare + chain) */
    check_true("PossibleZeroQ[N[D[Floor[x], x] /. x -> 11/10]]");
    check_true("PossibleZeroQ[N[D[Floor[3 x], x] /. x -> 7/10]]");
    /* 898 solves and its residual back-substitutes to a numeric zero (not UNK) */
    check_form("Head[DSolve[y''[x] + 9 y[x] == 2 Sec[3 x], y, x]]", "List");
    check_true("With[{s = DSolve[y''[x] + 9 y[x] == 2 Sec[3 x], y, x][[1]] /. {C[1] -> 13/10, C[2] -> 7/10}}, "
               "Abs[N[(y''[x] + 9 y[x] - 2 Sec[3 x]) /. s /. x -> 11/10, 20]] < 1/1000000]");
    /* anti-overfit: a sibling Sec-forced equation solves and verifies too */
    check_form("Head[DSolve[y''[x] + 4 y[x] == Sec[2 x], y, x]]", "List");
    check_true("With[{s = DSolve[y''[x] + 4 y[x] == Sec[2 x], y, x][[1]] /. {C[1] -> 13/10, C[2] -> 7/10}}, "
               "Abs[N[(y''[x] + 4 y[x] - Sec[2 x]) /. s /. x -> 11/10, 20]] < 1/1000000]");
}

/* ---- M30: forced linear system with an irrational spectrum (no Integrate blow-up) ---- */
/* corpus 2.2.10-924: {x'==2x+4y+3E^t, y'==5x-y-t^2} has eigenvalues (1+-Sqrt[89])/2.
 * The variation-of-parameters integral Integrate[E^{-lambda t} t^m] rationalised the
 * 1/lambda^k coefficient into a hundreds-of-digit integer and ran >90 s.  dsolve_linsys
 * now abstracts a REAL-irrational eigenvalue to a fresh symbol before the integral (and
 * Simplifies the algebraic coefficients, safe once the exponents are symbolic), then
 * substitutes the eigenvalue back -- fast and exact.  Complex spectra stay concrete. */
static void t_m30_linsys_irrational_forcing(void) {
    /* 924 itself: solved, and both residuals back-substitute to ~0 */
    check_true("MatchQ[DSolve[{x'[t]==2 x[t]+4 y[t]+3 E^t, y'[t]==5 x[t]-y[t]-t^2}, {x,y}, t], "
               "{{x -> _Function, y -> _Function}}]");
    check_true("With[{s = DSolve[{x'[t]==2 x[t]+4 y[t]+3 E^t, y'[t]==5 x[t]-y[t]-t^2}, {x,y}, t][[1]]}, "
               "Max[Abs[N[{(x'[t]-(2 x[t]+4 y[t]+3 E^t)), (y'[t]-(5 x[t]-y[t]-t^2))} /. s "
               "/. {C[1]->7/10, C[2]->13/10, t->3/10}, 20]]] < 1/1000000]");
    /* anti-overfit: a DIFFERENT 2x2 forced system with an irrational spectrum
     * (+-Sqrt[5]) and mixed exp/polynomial forcing solves and verifies too */
    check_true("With[{s = DSolve[{x'[t]==x[t]+4 y[t]+E^t, y'[t]==x[t]-y[t]+t}, {x,y}, t][[1]]}, "
               "Max[Abs[N[{(x'[t]-(x[t]+4 y[t]+E^t)), (y'[t]-(x[t]-y[t]+t))} /. s "
               "/. {C[1]->7/10, C[2]->13/10, t->3/10}, 20]]] < 1/1000000]");
    /* regression guard: a COMPLEX-spectrum forced system (924's sibling in the
     * corpus, 2.2.10-927, spectrum {2,-1+-I,0}) must stay concrete and still solve */
    check_true("MatchQ[DSolve[{x1'[t]==x2[t]+x3[t]+1, x2'[t]==x3[t]+x4[t]+t, "
               "x3'[t]==x1[t]+x4[t]+t^2, x4'[t]==x1[t]+x2[t]+t^3}, {x1,x2,x3,x4}, t], "
               "{{x1 -> _Function, x2 -> _Function, x3 -> _Function, x4 -> _Function}}]");
}

/* ---- M30: Kovacic closes an INHOMOGENEOUS variable-coefficient ODE ---- */
/* corpus 2.2.10-907: (x^2-1)y''-2x y'+2y == x^2-1 (Legendre-type) declined -- the
 * Kovacic solver handled only the homogeneous equation.  It now accepts a forcing
 * (dsolve_second_order_PQ_forced), de-obfuscates the fundamental set it recovers
 * (Sqrt[-1+x^2]E^(-1/2 Log[1+x]+3/2 Log[-1+x]) -> the clean {x,(x-1)^2}), and adds the
 * particular solution by variation of parameters, re-verifying the full solution. */
static void t_m30_kovacic_inhomogeneous(void) {
    /* 907: an explicit closed-form branch (not a decline), residual ~0 */
    check_form("Head[DSolve[(x^2-1) y''[x]-2 x y'[x]+2 y[x] == x^2-1, y, x][[1, 1]]]", "Rule");
    check_true("With[{s = DSolve[(x^2-1) y''[x]-2 x y'[x]+2 y[x] == x^2-1, y, x][[1]] "
               "/. {C[1]->6/5, C[2]->3/5}}, "
               "Abs[N[((x^2-1) y''[x]-2 x y'[x]+2 y[x] - (x^2-1)) /. s /. x->7/10, 20]] < 1/1000000]");
    /* anti-overfit: the SAME operator with a different forcing (== x) also closes */
    check_form("Head[DSolve[(x^2-1) y''[x]-2 x y'[x]+2 y[x] == x, y, x][[1, 1]]]", "Rule");
    check_true("With[{s = DSolve[(x^2-1) y''[x]-2 x y'[x]+2 y[x] == x, y, x][[1]] "
               "/. {C[1]->6/5, C[2]->3/5}}, "
               "Abs[N[((x^2-1) y''[x]-2 x y'[x]+2 y[x] - x) /. s /. x->7/10, 20]] < 1/1000000]");
    /* regression: the homogeneous equation still solves (unchanged path) */
    check_form("Head[DSolve[(x^2-1) y''[x]-2 x y'[x]+2 y[x] == 0, y, x]]", "List");
}

/* ---- M31: coupled DAG constant system — multi-term exponential forcing ---- */
/* corpus 2.2.11-1014: {x1'=2x1, x2'=-7x1+9x2+7x3, x3'=2x3} is a DAG (x1,x3 sources,
 * x2 sink) solved by TriangularSystem, which peels x1,x3 and hands the scalar engine
 * x2' = 9 x2 + 7(C[k]-C[j])E^{2x}.  The integrating-factor integrand mu*q =
 * E^{-9x}(7 C[k]E^{2x} - 7 C[j]E^{2x}) was passed to Integrate UN-combined, driving
 * the exponential-substitution path to a 55 s, branch-WRONG (-1)^(1/9) antiderivative
 * (kept because its residual is zero-test-undecidable).  Root fix in Integrate:
 * integrate.c:try_linearity now distributes a product over a sum factor (c(g+h)->cg+ch),
 * so the exponentials collapse to E^{-7x} (the clean single-product path) -- fast and
 * correct.  (The direct Integrate bug is guarded in test_integrate_dispatch.c.) */
static void t_m31_triangular_exp_forcing(void) {
    /* 1014 itself: solves; both the wrong branch and the slowness are gone (residual ~0) */
    check_true("And @@ (PossibleZeroQ /@ ({x1'[x]-2 x1[x], x2'[x]-(-7 x1[x]+9 x2[x]+7 x3[x]), x3'[x]-2 x3[x]} /. "
               "DSolve[{x1'[x]==2 x1[x], x2'[x]==-7 x1[x]+9 x2[x]+7 x3[x], x3'[x]==2 x3[x]}, {x1,x2,x3}, x][[1]]))");
    /* the peeled scalar equation in residual form with multi-term exponential forcing
     * (the exact shape TriangularSystem feeds) must solve clean -- a branch-wrong
     * antiderivative would make this residual nonzero */
    check_true("With[{s = DSolve[x2'[x] - (-7 a E^(2 x) + 9 x2[x] + 7 b E^(2 x)) == 0, x2, x][[1]]}, "
               "Abs[N[(x2'[x] - (-7 a E^(2 x) + 9 x2[x] + 7 b E^(2 x))) /. s "
               "/. {C[1]->7/10, a->3/10, b->9/10, x->6/5}, 20]] < 1/1000000]");
    /* anti-overfit: a DIFFERENT DAG constant system (other eigenvalues/forcing) */
    check_true("And @@ (PossibleZeroQ /@ ({x1'[x]-3 x1[x], x2'[x]-(2 x1[x]+5 x2[x]-4 x3[x]), x3'[x]-3 x3[x]} /. "
               "DSolve[{x1'[x]==3 x1[x], x2'[x]==2 x1[x]+5 x2[x]-4 x3[x], x3'[x]==3 x3[x]}, {x1,x2,x3}, x][[1]]))");
}

/* ---- M31: 4x4 constant system with large eigenvalues — answer is correct ---- */
/* corpus 2.2.11-1001: spectrum {16,32,48,64}; the matrix-exponential solution is
 * CORRECT (Simplify[residual]==0) but back-substitutes to a difference of E^{64x}-scale
 * terms that the corpus prelude's 20-digit numeric sweep at x~1.1..3 misread as nonzero
 * (catastrophic cancellation -> UNEVAL).  The solver was never wrong; the shared
 * verifier now re-checks a not-small sample at 200-digit precision.  This guards that
 * the solution is genuinely correct: residual ~0 at a small point (no cancellation noise). */
static void t_m31_linsys_large_eigenvalue(void) {
    check_true("MatchQ[DSolve[{x1'[x]==47 x1[x]-8 x2[x]+5 x3[x]-5 x4[x], x2'[x]==-10 x1[x]+32 x2[x]+18 x3[x]-2 x4[x], "
               "x3'[x]==139 x1[x]-40 x2[x]-167 x3[x]-121 x4[x], x4'[x]==-232 x1[x]+64 x2[x]+360 x3[x]+248 x4[x]}, {x1,x2,x3,x4}, x], "
               "{{x1 -> _Function, x2 -> _Function, x3 -> _Function, x4 -> _Function}}]");
    check_true("With[{s = DSolve[{x1'[x]==47 x1[x]-8 x2[x]+5 x3[x]-5 x4[x], x2'[x]==-10 x1[x]+32 x2[x]+18 x3[x]-2 x4[x], "
               "x3'[x]==139 x1[x]-40 x2[x]-167 x3[x]-121 x4[x], x4'[x]==-232 x1[x]+64 x2[x]+360 x3[x]+248 x4[x]}, {x1,x2,x3,x4}, x][[1]]}, "
               "Max[Abs[N[{x1'[x]-(47 x1[x]-8 x2[x]+5 x3[x]-5 x4[x]), x2'[x]-(-10 x1[x]+32 x2[x]+18 x3[x]-2 x4[x]), "
               "x3'[x]-(139 x1[x]-40 x2[x]-167 x3[x]-121 x4[x]), x4'[x]-(-232 x1[x]+64 x2[x]+360 x3[x]+248 x4[x])} /. s "
               "/. {C[1]->7/10, C[2]->9/10, C[3]->11/10, C[4]->13/10, x->1/10}, 30]]] < 1/1000000]");
}

/* ---- M32: §2.2.12 fixes (IVP fitter / Separable implicit / Integrate Erf variable) ---- */
static void t_m32_ivp_unsatisfiable_branch_dropped(void) {
    /* 1147 (was a WRONG answer): Sin[2x]+Cos[3y]y'==0, y[Pi/2]==0.  The "+" inverse branch
     * cannot pass through the initial point; it must be DROPPED, not emitted with an Undefined
     * constant.  Every returned branch is free of Undefined and C[_] and back-substitutes ~0. */
    check_true("With[{s = DSolve[{Sin[2 x] + Cos[3 y[x]] y'[x] == 0, y[Pi/2] == 0}, y, x]}, "
               "MatchQ[s, {__List}] && FreeQ[s, Undefined | C[_]] && "
               "AllTrue[s, Abs[N[(Sin[2 x] + Cos[3 y[x]] y'[x]) /. # /. x -> 13/10, 20]] < 1/10^6 &]]");
    /* 1143: the unsatisfiable -Sqrt branch drops, the +Sqrt fits to (-1+Sqrt[4x^2-15])/2 */
    check_true("With[{s = DSolve[{y'[x] == (2 x)/(1 + 2 y[x]), y[2] == 0}, y, x]}, "
               "MatchQ[s, {__List}] && FreeQ[s, C[_]] && "
               "AllTrue[s, PossibleZeroQ[(y[x] /. #) - (-1 + Sqrt[4 x^2 - 15])/2] &]]");
    /* 1138: a negative initial value needs the -Sqrt sign Bernoulli now also emits */
    check_true("With[{s = DSolve[{y'[x] == (1 - 2 x)/y[x], y[1] == -2}, y, x]}, "
               "MatchQ[s, {__List}] && FreeQ[s, C[_]] && "
               "AllTrue[s, (N[y[1] /. #, 20] < 0) && "
               "Abs[N[(y'[x] - (1 - 2 x)/y[x]) /. # /. x -> 6/5, 20]] < 1/10^6 &]]");
}
static void t_m32_bvp_underdetermined_keeps_constant(void) {
    /* REGRESSION GUARD: an under-determined BVP (the two BCs are dependent) must KEEP the free
     * constant -- the fitter drop must not discard y == C[2] Sin[x]. */
    check_true("Not[FreeQ[DSolve[{y''[x] + y[x] == 0, y[0] == 0, y[Pi] == 0}, y, x], C[2]]]");
}
static void t_m32_root_form_ivp(void) {
    /* 1149: cubic-in-y separable IVP -> the explicit Root branch cannot be C-fit, so the implicit
     * first integral fits C = G(x0,y0) with no inversion (fully fitted, no leftover constant). */
    check_true("With[{s = DSolve[{y'[x] == (3 x^2 + 1)/(-6 y[x] + 3 y[x]^2), y[0] == 1}, y, x]}, "
               "MatchQ[s, {__List}] && FreeQ[s, C[_]]]");
}
static void t_m32_separable_implicit_and_widened(void) {
    /* 1173: Cot[t]y/(1+y) -- elementary but non-invertible -> implicit first integral */
    check_true("MatchQ[DSolve[y'[t] == (Cot[t] y[t])/(1 + y[t]), y, t], {__List}]");
    /* 1186: autonomous, non-elementary y-integral -> implicit relation carrying the unevaluated
     * Integrate (verified by the implicit-function rule). */
    check_true("MatchQ[DSolve[y'[x] == -((2 ArcTan[y[x]])/(1 + y[x]^2)), y, x], {__List}]");
}
static void t_m32_integrate_erf_variable(void) {
    /* Integrate's Gaussian->Erf recognizer emitted a literal x for EVERY variable; now it threads
     * the real integration variable (fixed the non-x-variable Bernoulli solves 1182/1190). */
    check_true("FreeQ[Integrate[E^(a^2/2), a], x] && Not[FreeQ[Integrate[E^(a^2/2), a], a]]");
    check_true("PossibleZeroQ[D[Integrate[E^(t^2/2), t], t] - E^(t^2/2)]");
    /* regression: the x-variable case is unchanged */
    check_true("Not[FreeQ[Integrate[E^(x^2/2), x], x]]");
}

/* ---- M33: §2.2.13 exact-method robustness + IVP fall-through ----
 * Implicit first integrals are verified by the implicit-function rule: for an
 * exact M + N y' == 0 the returned relation G(x,y) == C has G_x = mu M, G_y = mu N,
 * so M + N(-G_x/G_y) == 0 identically -- checked here with PossibleZeroQ. */
static void t_m33_exact_transcendental(void) {
    /* 1201: exact with an E^(x y) coefficient -- Integrate[M,x] comes back in a Tan
     * half-angle form, so the potential is built from the OTHER coefficient (Path 2). */
    check_true("With[{G = DSolve[2 x - 2 E^(y[x] x) Sin[2 x] + E^(y[x] x) Cos[2 x] y[x] "
               "+ (-3 + E^(y[x] x) x Cos[2 x]) y'[x] == 0, y, x][[1,1,1]] /. y[x] -> Y}, "
               "PossibleZeroQ[(2 x - 2 E^(Y x) Sin[2 x] + E^(Y x) Cos[2 x] Y) "
               "+ (-3 + E^(Y x) x Cos[2 x]) (-D[G, x]/D[G, Y])]]");
    /* 1233: y' == P/Q with a NEGATIVE exponential E^(-x); cleared by the SYNTACTIC
     * denominator (Together would mis-factor E^(-x) as a spurious E^x). */
    check_true("With[{G = DSolve[y'[x] == (-E^(2 y[x]) Cos[x] + Cos[y[x]] E^(-x))/"
               "(2 E^(2 y[x]) Sin[x] - Sin[y[x]] E^(-x)), y, x][[1,1,1]] /. y[x] -> Y}, "
               "PossibleZeroQ[(E^(2 Y) Cos[x] - Cos[Y] E^(-x)) "
               "+ (2 E^(2 Y) Sin[x] - Sin[Y] E^(-x)) (-D[G, x]/D[G, Y])]]");
}
static void t_m33_exact_rational_clear(void) {
    /* 1216 (1/x,1/y poles, cleared via Together) and 1238 (y'==-P/Q polynomial):
     * cleared to an exact polynomial form -> Root-form explicit branches.  Mathilda
     * cannot differentiate a Root w.r.t. its parameter, so verify each branch lies on
     * the exact first integral F(x,y) == C (with C fixed to a generic value). */
    check_true("With[{s = DSolve[3 x + 6/y[x] + (x^2/y[x] + 3 y[x]/x) y'[x] == 0, y, x] "
               "/. C[1] -> 5}, MatchQ[s, {__List}] && AllTrue[s, Abs[N[(x^3 y[x] + 3 x^2 "
               "+ y[x]^3 - 5) /. # /. x -> 7/5, 20]] < 1/10^6 &]]");
    check_true("With[{s = DSolve[(-4 + 6 y[x] x + 2 y[x]^2)/(3 x^2 + 4 y[x] x + 3 y[x]^2) "
               "+ y'[x] == 0, y, x] /. C[1] -> 5}, MatchQ[s, {__List}] && AllTrue[s, "
               "Abs[N[(-4 x + 3 x^2 y[x] + 2 x y[x]^2 + y[x]^3 - 5) /. # /. x -> 7/5, 20]] "
               "< 1/10^6 &]]");
}
static void t_m33_exact_mu_trig(void) {
    /* 1214: exact via mu(y) = Sin y -- mu(y) is now tried even when mu(x)'s free-of
     * test mis-decides; and the implicit verify combines over a common denominator so
     * the telescoping Csc/Cot residual is not false-NEGATIVED by the numeric zero-test. */
    check_true("With[{G = DSolve[E^x + (E^x Cot[y[x]] + 2 Csc[y[x]] y[x]) y'[x] == 0, y, x]"
               "[[1,1,1]] /. y[x] -> Y}, PossibleZeroQ[E^x + "
               "(E^x Cot[Y] + 2 Csc[Y] Y) (-D[G, x]/D[G, Y])]]");
}
static void t_m33_exact_homogeneous_ivp_fallthrough(void) {
    /* 1205/1231: exact-AND-homogeneous IVPs.  Homogeneous runs first and returns a
     * transcendental log-form whose constant does not inverse-fit (Solve bubbles back);
     * dsolve_run now DECLINES that undecided IVP fit (FIT_UNDECIDED) so the cascade
     * reaches Exact, whose polynomial first integral fits the IC cleanly. */
    check_true("With[{s = DSolve[{2 x - y[x] + (-x + 2 y[x]) y'[x] == 0, y[1] == 3}, y, x]}, "
               "MatchQ[s, {__List}] && FreeQ[s, C[_]] && AllTrue[s, (Abs[N[y[1] /. #, 20] - 3] "
               "< 1/10^6) && Abs[N[(2 x - y[x] + (-x + 2 y[x]) y'[x]) /. # /. x -> 8/5, 20]] "
               "< 1/10^6 &]]");
    check_true("With[{s = DSolve[{x + y[x] + (x + 2 y[x]) y'[x] == 0, y[2] == 3}, y, x]}, "
               "MatchQ[s, {__List}] && FreeQ[s, C[_]] && AllTrue[s, Abs[N[y[2] /. #, 20] - 3] "
               "< 1/10^6 &]]");
    /* regression guard: a genuinely UNDER-determined BVP still keeps its free constant
     * (the fit SUCCEEDS with a residual C -- distinct from an undecided fit). */
    check_true("MatchQ[DSolve[{y''[x] + y[x] == 0, y[0] == 0, y[Pi] == 0}, y, x], "
               "{{y -> Function[{x}, _. C[2] Sin[x]]}}]");
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

/* ---- M8: general linear systems (defective / singular / triangular) ---- */
static void t_sys_defective_singular(void) {
    /* the reported case: coupled, A={{0,0},{-1,0}} is defective (one Jordan
     * block, eigenvalue 0 doubled) AND singular.  Old eigen-only linsys and
     * DecoupleSystem both declined; matrix exponential / triangular solve it. */
    check_form("Head[DSolve[{y'[t] + 1 == 1, x'[t] + y[t] == 0}, {y[t], x[t]}, t]]", "List");
    check_true("And @@ (PossibleZeroQ /@ ({y'[x], x'[x] + y[x]} /. "
               "DSolve[{y'[x] == 0, x'[x] + y[x] == 0}, {y, x}, x][[1]]))");
}
static void t_sys_defective_nontriangular(void) {
    /* eigenvalue 2 doubled, defective, NOT triangular -> x^k e^{2x} via Jordan */
    check_true("And @@ (PossibleZeroQ /@ ({u'[x] - (u[x] - v[x]), v'[x] - (u[x] + 3 v[x])} /. "
               "DSolve[{u'[x] == u[x] - v[x], v'[x] == u[x] + 3 v[x]}, {u, v}, x][[1]]))");
}
static void t_sys_triangular_varcoeff(void) {
    /* coupled-but-triangular at variable coefficient (matrix exponential cannot
     * reach this; TriangularSystem forward-substitution does). */
    check_form("Head[DSolve[{y'[x] == y[x]/x, z'[x] == y[x]}, {y, z}, x]]", "List");
    check_true("And @@ (PossibleZeroQ /@ ({y'[x] - y[x]/x, z'[x] - y[x]} /. "
               "DSolve[{y'[x] == y[x]/x, z'[x] == y[x]}, {y, z}, x][[1]]))");
}
static void t_sys_singular_forcing(void) {
    /* singular A with forcing -> variation of parameters (subsumes -A^{-1}b,
     * which does not exist for singular A). */
    check_true("And @@ (PossibleZeroQ /@ ({y'[x] - 1, x'[x] + y[x]} /. "
               "DSolve[{y'[x] == 1, x'[x] + y[x] == 0}, {y, x}, x][[1]]))");
}
static void t_sys_triangular_ivp(void) {
    /* triangular variable-coefficient IVP: y=2x, z=x^2-1 at x=1 */
    check_true("PossibleZeroQ[(z[x] /. DSolve[{y'[x] == y[x]/x, z'[x] == y[x], "
               "y[1] == 2, z[1] == 0}, {y[x], z[x]}, x][[1]]) - (x^2 - 1)]");
}

/* ---- M11: LinearSystemVarCoeff (scalar-factor A(x)=f(x)B, genuinely coupled) ---- */
static void t_sys_varcoeff_coupled(void) {
    /* A = (1/x){{2,1},{1,2}}, eigenvalues 1,3 -> x^1, x^3 modes.  Genuinely coupled,
     * non-triangular, variable-coefficient: neither DecoupleSystem, TriangularSystem,
     * nor the constant-A LinearFirstOrderSystem reaches it. */
    check_form("Head[DSolve[{y'[x] == (2 y[x] + z[x])/x, z'[x] == (y[x] + 2 z[x])/x}, {y, z}, x]]", "List");
    check_true("And @@ (PossibleZeroQ /@ ({y'[x] - (2 y[x] + z[x])/x, z'[x] - (y[x] + 2 z[x])/x} /. "
               "DSolve[{y'[x] == (2 y[x] + z[x])/x, z'[x] == (y[x] + 2 z[x])/x}, {y, z}, x][[1]]))");
}
static void t_sys_varcoeff_complex(void) {
    /* f = x, A = x{{0,1},{-1,0}}, complex spectrum -> real Cos/Sin of tau = x^2/2 */
    check_true("And @@ (PossibleZeroQ /@ ({y'[x] - x z[x], z'[x] + x y[x]} /. "
               "DSolve[{y'[x] == x z[x], z'[x] == -x y[x]}, {y, z}, x][[1]]))");
    /* the realified body is genuinely Cos[x^2/2]/Sin[x^2/2] (no Arg/Abs leak) */
    check_true("FreeQ[DSolve[{y'[x] == x z[x], z'[x] == -x y[x]}, {y, z}, x], Arg | Abs]");
}
static void t_sys_varcoeff_forced(void) {
    /* forced: A=(1/x){{2,1},{1,2}}, b={1,0}; VoP integral produces Log[x] terms,
     * which must NOT get split into Log[Abs[x]]+I Arg[x] by the realifier. */
    check_form("Head[DSolve[{y'[x] == (2 y[x] + z[x])/x + 1, z'[x] == (y[x] + 2 z[x])/x}, {y, z}, x]]", "List");
    check_true("And @@ (PossibleZeroQ /@ ({y'[x] - ((2 y[x] + z[x])/x + 1), z'[x] - (y[x] + 2 z[x])/x} /. "
               "DSolve[{y'[x] == (2 y[x] + z[x])/x + 1, z'[x] == (y[x] + 2 z[x])/x}, {y, z}, x][[1]]))");
}
static void t_sys_varcoeff_pinned_decline(void) {
    /* pinned method entry solves the coupled varcoeff system */
    check_form("Head[DSolve`LinearSystemVarCoeff[{y'[x] == (2 y[x] + z[x])/x, "
               "z'[x] == (y[x] + 2 z[x])/x}, {y, z}, x]]", "List");
    /* declines a CONSTANT-A system (that is LinearFirstOrderSystem's job) */
    check_form("Head[DSolve`LinearSystemVarCoeff[{y'[x] == 2 y[x] + z[x], "
               "z'[x] == y[x] + 2 z[x]}, {y, z}, x]]", "DSolve`LinearSystemVarCoeff");
    /* declines a non-scalar-factor A = {{1/x, 1}, {0, 1/x}} (not f(x)*constant) */
    check_form("Head[DSolve`LinearSystemVarCoeff[{y'[x] == y[x]/x + z[x], "
               "z'[x] == z[x]/x}, {y, z}, x]]", "DSolve`LinearSystemVarCoeff");
}

/* ---- M11: Sturm-Liouville DSolve`EigenvalueProblem (pinned first cut) ---- */
/* Verify by extracting the eigenvalue and eigenfunction, substituting the family
 * index C[1]->3 and amplitude C[2]->1, and checking the ODE + both BC residuals
 * vanish (PossibleZeroQ of exact Sin/Cos at integer multiples of Pi). */
static void t_eig_dirichlet(void) {
    check_form("Head[DSolve`EigenvalueProblem[{y''[x] + w y[x] == 0, y[0] == 0, y[Pi] == 0}, y, x]]", "List");
    check_true("Module[{s = DSolve`EigenvalueProblem[{y''[x] + w y[x] == 0, y[0] == 0, y[Pi] == 0}, y, x][[1]], lam, yf}, "
        "lam = First[w /. s[[1]]] /. C[1] -> 3; yf = (y /. s[[2]]) /. {C[1] -> 3, C[2] -> 1}; "
        "PossibleZeroQ[yf''[x] + lam yf[x]] && PossibleZeroQ[yf[0]] && PossibleZeroQ[yf[Pi]]]");
    /* eigenvalue family is n^2 on [0,Pi] */
    check_true("(First[w /. DSolve`EigenvalueProblem[{y''[x] + w y[x] == 0, y[0] == 0, y[Pi] == 0}, y, x][[1,1]]] "
        "/. C[1] -> 4) === 16");
}
static void t_eig_neumann(void) {
    check_true("Module[{s = DSolve`EigenvalueProblem[{y''[x] + w y[x] == 0, y'[0] == 0, y'[Pi] == 0}, y, x][[1]], lam, yf}, "
        "lam = First[w /. s[[1]]] /. C[1] -> 3; yf = (y /. s[[2]]) /. {C[1] -> 3, C[2] -> 1}; "
        "PossibleZeroQ[yf''[x] + lam yf[x]] && PossibleZeroQ[yf'[0]] && PossibleZeroQ[yf'[Pi]]]");
}
static void t_eig_mixed(void) {
    /* Dirichlet at 0, Neumann at 1: half-integer family, eigenfunction Sin */
    check_true("Module[{s = DSolve`EigenvalueProblem[{y''[x] + w y[x] == 0, y[0] == 0, y'[1] == 0}, y, x][[1]], lam, yf}, "
        "lam = First[w /. s[[1]]] /. C[1] -> 3; yf = (y /. s[[2]]) /. {C[1] -> 3, C[2] -> 1}; "
        "PossibleZeroQ[yf''[x] + lam yf[x]] && PossibleZeroQ[yf[0]] && PossibleZeroQ[yf'[1]]]");
    /* Neumann at 0, Dirichlet at Pi: eigenfunction Cos */
    check_true("Module[{s = DSolve`EigenvalueProblem[{y''[x] + w y[x] == 0, y'[0] == 0, y[Pi] == 0}, y, x][[1]], lam, yf}, "
        "lam = First[w /. s[[1]]] /. C[1] -> 3; yf = (y /. s[[2]]) /. {C[1] -> 3, C[2] -> 1}; "
        "PossibleZeroQ[yf''[x] + lam yf[x]] && PossibleZeroQ[yf'[0]] && PossibleZeroQ[yf[Pi]]]");
}
static void t_eig_no_misfire(void) {
    /* an ordinary IVP (inhomogeneous conditions) is NOT an eigenvalue problem */
    check_form("Head[DSolve`EigenvalueProblem[{y''[x] + w y[x] == 0, y[0] == 1, y'[0] == 0}, y, x]]",
               "DSolve`EigenvalueProblem");
    /* no free eigenparameter (coefficient is the number 1, not a symbol) */
    check_form("Head[DSolve`EigenvalueProblem[{y''[x] + y[x] == 0, y[0] == 0, y[Pi] == 0}, y, x]]",
               "DSolve`EigenvalueProblem");
}

/* ---- M4: reduction of order (2nd-order missing y) ---- */
static void t_reduce_order(void) {
    check_true("PossibleZeroQ[(y''[x] - y'[x]^2) /. DSolve[y''[x] == y'[x]^2, y, x][[1]]]");
    /* two independent constants (order preserved) */
    check_true("Not[FreeQ[DSolve[y''[x] == y'[x]^2, y, x][[1]], C[2]]]");
}

/* stress helper: assert DSolve[eqn, y, x] actually solves (does not decline) and
 * that its solution back-substitutes the residual to zero.  Checking Head is
 * essential: a declined DSolve leaves [[1]] symbolic, so the residual would not
 * substitute and PossibleZeroQ would return True vacuously. */
static void check_solves(const char* eqn, const char* residual) {
    char buf[640];
    snprintf(buf, sizeof(buf), "Head[DSolve[%s, y, x]]", eqn);
    check_form(buf, "List");
    snprintf(buf, sizeof(buf),
             "PossibleZeroQ[(%s) /. DSolve[%s, y, x][[1]]]", residual, eqn);
    check_true(buf);
}

/* Pinned-method solve check for an elementary method (Kovacic): Head===List then
 * the residual back-substitutes to zero.  `residual` is the ODE lhs (== 0). */
static void check_method(const char* method, const char* eqn, const char* residual) {
    char buf[900];
    snprintf(buf, sizeof(buf), "Head[%s[%s, y, x]]", method, eqn);
    check_form(buf, "List");
    snprintf(buf, sizeof(buf), "PossibleZeroQ[(%s) /. %s[%s, y, x][[1]]]", residual, method, eqn);
    check_true(buf);
}

/* Series-solution check: extract the SeriesData body (b = y[x] /. sol) and verify
 * the ODE `ode_of_b` (written with b, D[b,x], D[b,{x,2}]) is O[x]^k.  Going through
 * the body directly is required because Derivative[k][Function[{x}, SeriesData]] is
 * not reduced by the evaluator. */
static void check_series(const char* method, const char* eqn, const char* ode_of_b) {
    char buf[900];
    snprintf(buf, sizeof(buf), "Head[%s[%s, y, x]]", method, eqn);
    check_form(buf, "List");
    snprintf(buf, sizeof(buf),
             "PossibleZeroQ[Module[{b = y[x] /. %s[%s, y, x][[1]]}, %s]]", method, eqn, ode_of_b);
    check_true(buf);
}

/* ---- M5: NormalForm ---- */
static void t_normalform_bessel(void) {
    /* Bessel nu=2: y'' + (1/x) y' + (1 - 4/x^2) y == 0 -> r = -1 + 15/(4 x^2), w = 1/Sqrt[x] */
    check_true("PossibleZeroQ[DSolve`NormalForm[x^2 y''[x] + x y'[x] + (x^2 - 4) y[x] == 0, y, x][[1]] "
               "- (-1 + 15/(4 x^2))]");
    check_true("PossibleZeroQ[DSolve`NormalForm[x^2 y''[x] + x y'[x] + (x^2 - 4) y[x] == 0, y, x][[2]] "
               "- 1/Sqrt[x]]");
}
static void t_normalform_const(void) {
    /* y'' + 2 y' + y == 0 -> r == 0, w == E^(-x) */
    check_true("PossibleZeroQ[DSolve`NormalForm[y''[x] + 2 y'[x] + y[x] == 0, y, x][[1]]]");
    check_true("PossibleZeroQ[DSolve`NormalForm[y''[x] + 2 y'[x] + y[x] == 0, y, x][[2]] - E^(-x)]");
}
static void t_normalform_declines(void) {
    /* not a homogeneous 2nd-order linear ODE: stays symbolic */
    check_form("Head[DSolve`NormalForm[y'[x] + y[x] == 0, y, x]]", "DSolve`NormalForm");
    check_form("Head[DSolve`NormalForm[y''[x] + y[x]^2 == 0, y, x]]", "DSolve`NormalForm");
}

/* ---- M5: Kovacic (Cases 1 & 2) ---- */
static void t_kovacic_case1_exp(void) {
    /* z'' = (1 + x^2) z -> Exp[x^2/2] family (Case 1, polynomial omega) */
    check_method("DSolve`Kovacic", "y''[x] - (1 + x^2) y[x] == 0", "D[y[x],{x,2}] - (1 + x^2) y[x]");
}
static void t_kovacic_apparent_singularity(void) {
    /* z'' = (x^2 + 3) z -> x Exp[x^2/2]: an apparent singularity (polynomial P factor) */
    check_method("DSolve`Kovacic", "y''[x] - (x^2 + 3) y[x] == 0", "D[y[x],{x,2}] - (x^2 + 3) y[x]");
}
static void t_kovacic_case1_pole(void) {
    /* z'' = (2/x^2) z -> x^2, 1/x (Case 1, order-2 pole) */
    check_method("DSolve`Kovacic", "y''[x] - (2/x^2) y[x] == 0", "D[y[x],{x,2}] - (2/x^2) y[x]");
}
static void t_kovacic_case2(void) {
    /* z'' = (x/4 + 5/(16 x^2)) z -> x^(-1/4) Exp[+-x^(3/2)/3] (Case 2, degree-2 algebraic) */
    check_method("DSolve`Kovacic", "y''[x] - (x/4 + 5/(16 x^2)) y[x] == 0",
                 "D[y[x],{x,2}] - (x/4 + 5/(16 x^2)) y[x]");
}
static void t_kovacic_auto_closed_form(void) {
    /* the automatic cascade prefers Kovacic's closed form over a series */
    check_form("FreeQ[DSolve[y''[x] - (x^2 + 3) y[x] == 0, y, x], SeriesData]", "True");
    check_form("Head[DSolve[y''[x] - (x^2 + 3) y[x] == 0, y, x]]", "List");
}
static void t_kovacic_declines(void) {
    /* Bessel is not Liouvillian: Kovacic declines (Case 2 must not false-positive) */
    check_form("Head[DSolve`Kovacic[x^2 y''[x] + x y'[x] + (x^2 - 4) y[x] == 0, y, x]]", "DSolve`Kovacic");
}

/* ---- Case 1c: apparent singularities for rational r (Legendre / Chebyshev /
 *      Gegenbauer): the classical monic-P completion over the pole exponents ---- */
static void t_kovacic_legendre1(void) {
    /* Legendre n=1: (1-x^2)y'' - 2x y' + 2y == 0 -> C[1] x + C[2](1 - x ArcTanh[x]) */
    check_method("DSolve`Kovacic", "(1 - x^2) y''[x] - 2 x y'[x] + 2 y[x] == 0",
                 "(1 - x^2) D[y[x],{x,2}] - 2 x D[y[x],x] + 2 y[x]");
}
static void t_kovacic_legendre2(void) {
    /* Legendre n=2: (1-x^2)y'' - 2x y' + 6y == 0 (first solution ∝ P_2) */
    check_method("DSolve`Kovacic", "(1 - x^2) y''[x] - 2 x y'[x] + 6 y[x] == 0",
                 "(1 - x^2) D[y[x],{x,2}] - 2 x D[y[x],x] + 6 y[x]");
}
static void t_kovacic_chebyshev2(void) {
    /* Chebyshev n=2: (1-x^2)y'' - x y' + 4y == 0 (first solution ∝ T_2) */
    check_method("DSolve`Kovacic", "(1 - x^2) y''[x] - x y'[x] + 4 y[x] == 0",
                 "(1 - x^2) D[y[x],{x,2}] - x D[y[x],x] + 4 y[x]");
}
static void t_kovacic_complex_poles(void) {
    /* poles at ±i: z'' == ((3 + 2 x^2)/(1 + x^2)^2) z -> x Sqrt[1+x^2] & Sqrt[1+x^2](1 + x ArcTan[x]) */
    check_method("DSolve`Kovacic", "y''[x] - ((3 + 2 x^2)/(1 + x^2)^2) y[x] == 0",
                 "D[y[x],{x,2}] - ((3 + 2 x^2)/(1 + x^2)^2) y[x]");
}
static void t_kovacic_legendre_auto_closed_form(void) {
    /* the automatic cascade returns Legendre's closed form, not a truncated series */
    check_form("FreeQ[DSolve[(1 - x^2) y''[x] - 2 x y'[x] + 2 y[x] == 0, y, x], SeriesData]", "True");
    check_form("Head[DSolve[(1 - x^2) y''[x] - 2 x y'[x] + 2 y[x] == 0, y, x]]", "List");
}
static void t_kovacic_case2_complex_pole_no_hang(void) {
    /* (x^3+1)y'' + x y' + y == 0: r has a complex-conjugate pole pair; no Liouvillian
     * solution.  Case 1c declines fast and the Case-2 guard skips the σ-solve that
     * used to hang, so DSolve returns (a series) rather than spinning.  The test
     * passing at all is the no-hang assertion. */
    check_form("Head[DSolve[(x^3 + 1) y''[x] + x y'[x] + y[x] == 0, y, x]]", "List");
}

/* ---- OperatorFactor: higher-order linear-operator factoring + DFactor ---- */
static void t_method_operfactor(void) {
    /* shifted-Euler at x=1 (pole != 0, so EulerCauchy declines): OperatorFactor's
       unique niche.  L = (D-1/(x-1))(D-2/(x-1))(D-4/(x-1)), cleared of (x-1)^3. */
    check_method("DSolve`OperatorFactor",
        "(x-1)^3 y'''[x] - 7 (x-1)^2 y''[x] + 18 (x-1) y'[x] - 18 y[x] == 0",
        "(x-1)^3 y'''[x] - 7 (x-1)^2 y''[x] + 18 (x-1) y'[x] - 18 y[x]");
    /* three arbitrary constants present (a full order-3 general solution) */
    check_true("Not[FreeQ[DSolve`OperatorFactor["
        "(x-1)^3 y'''[x] - 7 (x-1)^2 y''[x] + 18 (x-1) y'[x] - 18 y[x] == 0, y, x], C[3]]]");
}
static void t_operfactor_more(void) {
    /* constant-coefficient order 3: (D-1)(D-2)(D-3) */
    check_method("DSolve`OperatorFactor",
        "y'''[x] - 6 y''[x] + 11 y'[x] - 6 y[x] == 0",
        "y'''[x] - 6 y''[x] + 11 y'[x] - 6 y[x]");
    /* order 4: (D-1)(D-2)(D-3)(D-4) */
    check_method("DSolve`OperatorFactor",
        "y''''[x] - 10 y'''[x] + 35 y''[x] - 50 y'[x] + 24 y[x] == 0",
        "y''''[x] - 10 y'''[x] + 35 y''[x] - 50 y'[x] + 24 y[x]");
    /* resonant repeated shifted-Euler factor (secular Log term in the basis) */
    check_method("DSolve`OperatorFactor",
        "(x-1)^3 y'''[x] - 5 (x-1)^2 y''[x] + 10 (x-1) y'[x] - 10 y[x] == 0",
        "(x-1)^3 y'''[x] - 5 (x-1)^2 y''[x] + 10 (x-1) y'[x] - 10 y[x]");
}
static void t_operfactor_ivp(void) {
    check_form("Head[DSolve[{y'''[x] - 6 y''[x] + 11 y'[x] - 6 y[x] == 0, "
               "y[0]==0, y'[0]==0, y''[0]==2}, y[x], x]]", "List");
}
static void t_operfactor_declines(void) {
    /* order-3 with no rational first-order factor (Airy-type) stays symbolic */
    check_form("Head[DSolve`OperatorFactor[y'''[x] - x y[x] == 0, y[x], x]]", "DSolve`OperatorFactor");
    /* order 2 belongs to Kovacic; the n>=3 guard declines */
    check_form("Head[DSolve`OperatorFactor[y''[x] - x y[x] == 0, y[x], x]]", "DSolve`OperatorFactor");
}
static void t_dfactor(void) {
    /* factor the constant operator into three first-order factors {Dx-1,Dx-2,Dx-3} */
    check_form("Length[DSolve`DFactor[y'''[x] - 6 y''[x] + 11 y'[x] - 6 y[x] == 0, y[x], x]]", "3");
    /* reconstruct: applying the Dx-factors (innermost first, i.e. in list order) to
       a concrete test function must reproduce the operator applied to it. */
    check_true("Module[{fs = DSolve`DFactor[y'''[x] - 6 y''[x] + 11 y'[x] - 6 y[x] == 0, y[x], x], "
               "tf = Exp[x] + x^4, recon, opv}, "
               "recon = Fold[Function[{w, f}, D[w, x] + (f /. Dx -> 0) w], tf, fs]; "
               "opv = D[tf,{x,3}] - 6 D[tf,{x,2}] + 11 D[tf,x] - 6 tf; "
               "PossibleZeroQ[recon - opv]]");
}
static void t_operfactor_auto(void) {
    /* the shifted-Euler flagship solves through the automatic cascade slot */
    check_solves("(x-1)^3 y'''[x] - 7 (x-1)^2 y''[x] + 18 (x-1) y'[x] - 18 y[x] == 0",
                 "(x-1)^3 y'''[x] - 7 (x-1)^2 y''[x] + 18 (x-1) y'[x] - 18 y[x]");
}

/* ---- M5: Frobenius / PowerSeries ---- */
static void t_powerseries_ordinary(void) {
    /* pinned power series about the ordinary point 0 for y'' + y == 0 (cos/sin) */
    check_series("DSolve`PowerSeries", "y''[x] + y[x] == 0", "D[b,{x,2}] + b");
}
static void t_powerseries_auto(void) {
    /* auto cascade falls through to a series for an analytic-but-non-closed-form coeff */
    check_form("Head[DSolve[y''[x] + Sin[x] y[x] == 0, y, x]]", "List");
    check_series("DSolve", "y''[x] + Sin[x] y[x] == 0", "D[b,{x,2}] + Sin[x] b");
}
static void t_frobenius_regsing_distinct(void) {
    /* 4 x y'' + 2 y' + y == 0: regular singular, roots 0 & 1/2 (cos/sin of Sqrt[x]) */
    check_series("DSolve`FrobeniusSeries", "4 x y''[x] + 2 y'[x] + y[x] == 0",
                 "4 x D[b,{x,2}] + 2 D[b,x] + b");
}
static void t_frobenius_regsing_log(void) {
    /* x y'' + y' == 0: regular singular, double root 0 -> {1, Log[x]} (Log solution) */
    check_series("DSolve`FrobeniusSeries", "x y''[x] + y'[x] == 0", "x D[b,{x,2}] + D[b,x]");
    check_true("Not[FreeQ[DSolve`FrobeniusSeries[x y''[x] + y'[x] == 0, y, x], Log]]");
}
static void t_frobenius_declines_irregular(void) {
    /* essential singularity at 0: series fallback declines */
    check_form("Head[DSolve`FrobeniusSeries[y''[x] + Exp[1/x] y[x] == 0, y, x]]", "DSolve`FrobeniusSeries");
}

/* ---- 1a: FirstOrderSubstitution — y'==F(a x + b y + c) ---- */
static void t_fos_quadratic(void) {
    check_true("PossibleZeroQ[(y'[x] - (x + y[x])^2) /. "
               "DSolve[y'[x] == (x + y[x])^2, y, x][[1]]]");
    check_form("Head[DSolve[y'[x] == (x + y[x])^2, y, x]]", "List");
}
static void t_fos_shifted(void) {
    check_true("PossibleZeroQ[(y'[x] - (x + y[x] + 1)^2) /. "
               "DSolve[y'[x] == (x + y[x] + 1)^2, y, x][[1]]]");
}
static void t_fos_distinct_coeff(void) {
    /* combination 2y - x: ratio r = -1/2 recovered after a common factor cancels */
    check_true("PossibleZeroQ[(y'[x] - (2 y[x] - x)^2) /. "
               "DSolve[y'[x] == (2 y[x] - x)^2, y, x][[1]]]");
}
static void t_fos_method(void) {
    check_true("PossibleZeroQ[(y'[x] - (x + y[x])^2) /. "
               "DSolve`FirstOrderSubstitution[y'[x] == (x + y[x])^2, y, x][[1]]]");
}
static void t_lie_abaco2_similar(void) {
    /* Sqrt[x+y]: FirstOrderSubstitution/Separable cannot invert the antiderivative,
     * and the rational Lie heuristics (abaco1_simple/linear/abaco1_product) decline,
     * but abaco2_similar (Cheb-Terrab & Roche §4.3) finds the symmetry
     * [F(x), H(x)] = [1, -1] and returns the verified implicit first integral
     * x - 2 Sqrt[x+y] + 2 Log[1 + Sqrt[x+y]] == C[1] — coverage the rational-only
     * heuristics lack.  Needs init.m (the antiderivative uses a CRC table rule). */
    check_true("Head[DSolve`LieSymmetry[y'[x] == Sqrt[x + y[x]], y, x][[1,1]]] === Equal");
    check_true("Head[DSolve[y'[x] == Sqrt[x + y[x]], y[x], x][[1,1]]] === Equal");
    check_true("PossibleZeroQ[Module[{eq = DSolve[y'[x] == Sqrt[x + y[x]], y[x], x][[1,1]]}, "
               "D[eq[[1]] - eq[[2]], x] /. y'[x] -> Sqrt[x + y[x]]]]");
}
static void t_fos_stress(void) {
    char eqn[256], res[256];
    int ks[] = {1, 2}, as[] = {1, 2}, cs[] = {-1, 0, 1, 2};
    for (size_t ki = 0; ki < 2; ki++)
        for (size_t ai = 0; ai < 2; ai++)
            for (size_t ci = 0; ci < 4; ci++) {
                int k = ks[ki], a = as[ai], c = cs[ci];
                snprintf(eqn, sizeof(eqn), "y'[x] == %d (%d x + %d y[x] + %d)^2", k, a, a, c);
                snprintf(res, sizeof(res), "y'[x] - %d (%d x + %d y[x] + %d)^2", k, a, a, c);
                check_solves(eqn, res);
            }
    int ms[] = {2, 3, 4};                       /* distinct x/y coefficients, r=m */
    for (size_t mi = 0; mi < 3; mi++) {
        snprintf(eqn, sizeof(eqn), "y'[x] == (%d x + y[x])^2", ms[mi]);
        snprintf(res, sizeof(res), "y'[x] - (%d x + y[x])^2", ms[mi]);
        check_solves(eqn, res);
    }
}

/* ---- 1a: Riccati — y'==q0(x)+q1(x) y+q2(x) y^2 (linearise y=-u'/(q2 u)) ---- */
static void t_method_riccati(void) {
    /* from-spectrum r=-1,-2: q1=r1+r2=-3, q0=r1 r2=2 -> elementary (Tanh) */
    check_method("DSolve`Riccati", "y'[x] == 2 - 3 y[x] + y[x]^2",
                 "y'[x] - (2 - 3 y[x] + y[x]^2)");
}
static void t_riccati_more(void) {
    /* r=2,3 (elementary) */
    check_method("DSolve`Riccati", "y'[x] == 6 + 5 y[x] + y[x]^2",
                 "y'[x] - (6 + 5 y[x] + y[x]^2)");
    /* genuine Riccati whose linearisation is the Airy equation u''==-x u: the
     * body is an AiryAi/AiryBi ratio, verified by numeric PossibleZeroQ sampling */
    check_solves("y'[x] == y[x]^2 + x", "y'[x] - (y[x]^2 + x)");
    /* variable-coefficient q2=x: u'' - (1/x) u' + x^2 u == 0 */
    check_form("Head[DSolve[y'[x] == x + x y[x]^2, y, x]]", "List");
    /* q2==0 is linear, not Riccati: the pinned method declines (stays symbolic) */
    check_form("Head[DSolve`Riccati[y'[x] == x + y[x], y[x], x]]", "DSolve`Riccati");
}
static void t_ivp_riccati(void) {
    /* logistic-type IVP: y'==y^2-y, y[0]==1/2 -> y == 1/(1+E^x) */
    check_form("Head[DSolve[{y'[x] == y[x]^2 - y[x], y[0] == 1/2}, y, x]]", "List");
    /* value at x=0 is 1/2 (constant fitted) */
    check_true("PossibleZeroQ[((y[x] /. DSolve[{y'[x] == y[x]^2 - y[x], y[0] == 1/2}, y, x][[1]]) "
               "/. x -> 0) - 1/2]");
    /* residual back-substitutes to zero */
    check_true("PossibleZeroQ[(y'[x] - (y[x]^2 - y[x])) /. "
               "DSolve[{y'[x] == y[x]^2 - y[x], y[0] == 1/2}, y, x][[1]]]");
}

/* ---- 1d: AutonomousReduction — y''==f(y, y') missing x ---- */
static void t_auto_exp(void) {
    check_true("PossibleZeroQ[(y[x] y''[x] - y'[x]^2) /. "
               "DSolve[y[x] y''[x] == y'[x]^2, y, x][[1]]]");
    /* order preserved: a second independent constant is present */
    check_true("Not[FreeQ[DSolve[y[x] y''[x] == y'[x]^2, y, x][[1]], C[2]]]");
}
static void t_auto_power(void) {
    check_true("PossibleZeroQ[(2 y[x] y''[x] - y'[x]^2) /. "
               "DSolve[2 y[x] y''[x] == y'[x]^2, y, x][[1]]]");
}
static void t_auto_reciprocal(void) {
    check_true("PossibleZeroQ[(y[x] y''[x] - 2 y'[x]^2) /. "
               "DSolve[y[x] y''[x] == 2 y'[x]^2, y, x][[1]]]");
}
static void t_auto_method(void) {
    check_true("PossibleZeroQ[(y[x] y''[x] - y'[x]^2) /. "
               "DSolve`AutonomousReduction[y[x] y''[x] == y'[x]^2, y, x][[1]]]");
}
static void t_auto_declines_elliptic(void) {
    /* y''==2y^3 reduces to an elliptic integral: stays symbolic (not wrong) */
    check_form("Head[DSolve[y''[x] == 2 y[x]^3, y[x], x]]", "DSolve");
}
static void t_auto_stress(void) {
    char eqn[256], res[256];
    int ab[][2] = {{1, 1}, {1, 2}, {2, 1}, {1, 3}, {3, 1}};
    for (size_t i = 0; i < 5; i++) {
        int a = ab[i][0], b = ab[i][1];
        snprintf(eqn, sizeof(eqn), "%d y[x] y''[x] == %d y'[x]^2", a, b);
        snprintf(res, sizeof(res), "%d y[x] y''[x] - %d y'[x]^2", a, b);
        check_solves(eqn, res);
    }
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

/* ---- M6 (Phase 2): second-order constant-coefficient linear PDE (operator
 * factoring).  Verified with two distinct concrete arbitrary functions
 * (C[1] -> Sin, C[2] -> Cos/Exp) reducing the Function application. ---- */
static void t_pde2_wave(void) {
    /* wave / d'Alembert: u_tt == c^2 u_xx -> C[1][x - c t] + C[2][x + c t] */
    check_true("With[{uc = (u[t,x] /. DSolve[D[u[t,x],{t,2}] == c^2 D[u[t,x],{x,2}], u, {t,x}][[1]]) "
               "/. {C[1] -> Sin, C[2] -> Cos}}, PossibleZeroQ[D[uc,{t,2}] - c^2 D[uc,{x,2}]]]");
}
static void t_pde2_laplace(void) {
    /* elliptic: complex-characteristic form C[1][y - I x] + C[2][y + I x] */
    check_true("With[{uc = (u[x,y] /. DSolve[D[u[x,y],{x,2}] + D[u[x,y],{y,2}] == 0, u, {x,y}][[1]]) "
               "/. {C[1] -> Sin, C[2] -> Cosh}}, PossibleZeroQ[D[uc,{x,2}] + D[uc,{y,2}]]]");
}
static void t_pde2_mixed(void) {
    /* pure mixed u_xy == 0 -> C[1][x] + C[2][y] */
    check_true("With[{uc = (u[x,y] /. DSolve[D[u[x,y],x,y] == 0, u, {x,y}][[1]]) "
               "/. {C[1] -> Exp, C[2] -> Sin}}, PossibleZeroQ[D[uc,x,y]]]");
}
static void t_pde2_repeated(void) {
    /* repeated root (parabolic): (D_x - D_y)^2 u == 0 -> C[1][w] + x C[2][w] */
    check_true("With[{uc = (u[x,y] /. DSolve[D[u[x,y],{x,2}] - 2 D[u[x,y],x,y] + D[u[x,y],{y,2}] == 0, "
               "u, {x,y}][[1]]) /. {C[1] -> Exp, C[2] -> Sin}}, "
               "PossibleZeroQ[D[uc,{x,2}] - 2 D[uc,x,y] + D[uc,{y,2}]]]");
}
static void t_pde2_distinct_asymmetric(void) {
    /* distinct rational roots (-2, -1/2): 2 u_xx + 5 u_xy + 2 u_yy == 0 */
    check_true("With[{uc = (u[x,y] /. DSolve[2 D[u[x,y],{x,2}] + 5 D[u[x,y],x,y] + 2 D[u[x,y],{y,2}] == 0, "
               "u, {x,y}][[1]]) /. {C[1] -> Sin, C[2] -> Exp}}, "
               "PossibleZeroQ[2 D[uc,{x,2}] + 5 D[uc,x,y] + 2 D[uc,{y,2}]]]");
}
static void t_pde2_pinned_and_declines(void) {
    /* pinned method solves the wave equation */
    check_true("With[{uc = (u[t,x] /. DSolve`PDELinearSecondOrder[D[u[t,x],{t,2}] == 4 D[u[t,x],{x,2}], "
               "u, {t,x}][[1]]) /. {C[1] -> Sin, C[2] -> Cos}}, "
               "PossibleZeroQ[D[uc,{t,2}] - 4 D[uc,{x,2}]]]");
    /* the 2nd-order method declines a first-order PDE (head stays symbolic) */
    check_form("Head[DSolve`PDELinearSecondOrder[D[u[x,y],x] + D[u[x,y],y] == 0, u, {x,y}]]",
               "DSolve`PDELinearSecondOrder");
    /* it declines a scalar ODE too */
    check_form("Head[DSolve`PDELinearSecondOrder[y''[x] == y[x], y, x]]",
               "DSolve`PDELinearSecondOrder");
    /* a first-order PDE still solves via the auto cascade (pde1) */
    check_true("With[{uc = (u[x,y] /. DSolve[D[u[x,y],x] + 3 D[u[x,y],y] == 0, u, {x,y}][[1]]) "
               "/. C[1][z_] :> Sin[z]}, PossibleZeroQ[D[uc,x] + 3 D[uc,y]]]");
}

/* ---- M6 (Phase 2): separation of variables — DSolve`SeparationOfVariables
 * (pinned-only). Product mode u = X(v1) Y(v2), back-substitution verified. ---- */
static void t_pdesep(void) {
    /* heat equation u_t == u_xx -> product mode, residual zero */
    check_true("With[{r = DSolve`SeparationOfVariables[D[u[x,t],t] == D[u[x,t],{x,2}], u, {x,t}]}, "
               "Head[r] === List && r =!= {} && "
               "PossibleZeroQ[(D[#,t] - D[#,{x,2}]) &[u[x,t] /. r[[1]]]]]");
    /* heat with a symbolic diffusivity k */
    check_true("With[{r = DSolve`SeparationOfVariables[D[u[x,t],t] == k D[u[x,t],{x,2}], u, {x,t}]}, "
               "Head[r] === List && PossibleZeroQ[(D[#,t] - k D[#,{x,2}]) &[u[x,t] /. r[[1]]]]]");
    /* Helmholtz-type u_xx + u_yy + u == 0 (both sides 2nd order; has a u term) */
    check_true("With[{r = DSolve`SeparationOfVariables[D[u[x,y],{x,2}] + D[u[x,y],{y,2}] + u[x,y] == 0, "
               "u, {x,y}]}, PossibleZeroQ[(D[#,{x,2}] + D[#,{y,2}] + #) &[u[x,y] /. r[[1]]]]]");
    /* declines a mixed-derivative term (pde2 territory) */
    check_form("Head[DSolve`SeparationOfVariables[D[u[x,y],x,y] == 0, u, {x,y}]]",
               "DSolve`SeparationOfVariables");
    /* declines inhomogeneous forcing */
    check_form("Head[DSolve`SeparationOfVariables[D[u[x,t],t] - D[u[x,t],{x,2}] == 1, u, {x,t}]]",
               "DSolve`SeparationOfVariables");
    /* declines an ODE-in-disguise (no v2 derivative) */
    check_form("Head[DSolve`SeparationOfVariables[D[u[x,t],{x,2}] == 0, u, {x,t}]]",
               "DSolve`SeparationOfVariables");
}

/* ---- M6 (Phase 2): PDEClassify — discriminant B^2 - 4 A C of the principal part. */
static void t_pdeclassify(void) {
    check_form("PDEClassify[D[u[t,x],{t,2}] == 4 D[u[t,x],{x,2}], u, {t,x}]", "\"Hyperbolic\"");
    check_form("PDEClassify[D[u[x,y],{x,2}] + D[u[x,y],{y,2}] == 0, u, {x,y}]", "\"Elliptic\"");
    check_form("PDEClassify[D[u[x,t],t] == D[u[x,t],{x,2}], u, {x,t}]", "\"Parabolic\"");
    /* mixed term u_xy: Δ = 1 > 0 -> Hyperbolic */
    check_form("PDEClassify[D[u[x,y],x,y] == 0, u, {x,y}]", "\"Hyperbolic\"");
    /* only the principal part matters: lower-order terms do not change the type */
    check_form("PDEClassify[D[u[x,y],{x,2}] + D[u[x,y],{y,2}] + D[u[x,y],x] + u[x,y] == 0, u, {x,y}]",
               "\"Elliptic\"");
    /* mixed-type / parameter-dependent discriminant (Tricomi) stays unevaluated */
    check_form("Head[PDEClassify[y D[u[x,y],{x,2}] + D[u[x,y],{y,2}] == 0, u, {x,y}]]", "PDEClassify");
    /* a first-order PDE has no principal 2nd-order part -> unevaluated */
    check_form("Head[PDEClassify[D[u[x,y],x] + D[u[x,y],y] == 0, u, {x,y}]]", "PDEClassify");
}

/* ---- M6 (Phase 2): wave-equation IVP by d'Alembert's formula (DSolve auto +
 * pinned DSolve`WaveDAlembert).  Concrete data is fully back-substitution
 * verified; undefined data is checked at the displacement condition (the velocity
 * integral over a zero-width interval vanishes). ---- */
static void t_wave_dalembert(void) {
    /* f = Sin, g = 0, c = 2: u = (Sin[x-2t] + Sin[x+2t])/2, verified */
    check_true("With[{b = u[x,t] /. DSolve[{D[u[x,t],{t,2}] == 4 D[u[x,t],{x,2}], "
               "u[x,0] == Sin[x], Derivative[0,1][u][x,0] == 0}, u, {x,t}][[1]]}, "
               "PossibleZeroQ[D[b,{t,2}] - 4 D[b,{x,2}]] && PossibleZeroQ[(b /. t->0) - Sin[x]] "
               "&& PossibleZeroQ[D[b,t] /. t->0]]");
    /* f = 0, g = Sin, c = 1: pure velocity, verified */
    check_true("With[{b = u[x,t] /. DSolve[{D[u[x,t],{t,2}] == D[u[x,t],{x,2}], "
               "u[x,0] == 0, Derivative[0,1][u][x,0] == Sin[x]}, u, {x,t}][[1]]}, "
               "PossibleZeroQ[D[b,{t,2}] - D[b,{x,2}]] && PossibleZeroQ[b /. t->0] "
               "&& PossibleZeroQ[(D[b,t] /. t->0) - Sin[x]]]");
    /* undefined f, g (auto-dispatch): solved, and the displacement IC holds */
    check_true("With[{r = DSolve[{D[u[x,t],{t,2}] == c^2 D[u[x,t],{x,2}], u[x,0] == f[x], "
               "Derivative[0,1][u][x,0] == g[x]}, u, {x,t}]}, Head[r] === List && r =!= {} && "
               "PossibleZeroQ[((u[x,t] /. r[[1]]) /. t->0) - f[x]]]");
    /* pinned method solves the same */
    check_true("With[{b = u[x,t] /. DSolve`WaveDAlembert[{D[u[x,t],{t,2}] == D[u[x,t],{x,2}], "
               "u[x,0] == Cos[x], Derivative[0,1][u][x,0] == 0}, u, {x,t}][[1]]}, "
               "PossibleZeroQ[D[b,{t,2}] - D[b,{x,2}]] && PossibleZeroQ[(b /. t->0) - Cos[x]]]");
    /* declines an elliptic equation with 'ICs' (not a wave) */
    check_form("Head[DSolve`WaveDAlembert[{D[u[x,y],{x,2}] + D[u[x,y],{y,2}] == 0, "
               "u[x,0] == f[x], Derivative[0,1][u][x,0] == g[x]}, u, {x,y}]]",
               "DSolve`WaveDAlembert");
    /* declines the bare wave equation (no initial conditions, neq == 1) */
    check_form("Head[DSolve`WaveDAlembert[D[u[x,t],{t,2}] == c^2 D[u[x,t],{x,2}], u, {x,t}]]",
               "DSolve`WaveDAlembert");
}

/* ---- M6 (Phase 2): heat-equation Cauchy problem by the heat kernel (DSolve auto
 * + pinned DSolve`HeatKernel).  The Gaussian convolution is nonelementary, so the
 * solution is the unevaluated heat-kernel integral; the kernel itself is verified
 * to solve the PDE (which the method checks internally). ---- */
static void t_heat_kernel(void) {
    /* auto-dispatch: solved (List), heat-kernel convolution present, no decline */
    check_true("With[{r = DSolve[{D[u[x,t],t] == k D[u[x,t],{x,2}], u[x,0] == f[x]}, u, {x,t}]}, "
               "Head[r] === List && r =!= {} && FreeQ[r, DSolve] && !FreeQ[r, Integrate]]");
    /* pinned method, same */
    check_true("With[{r = DSolve`HeatKernel[{D[u[x,t],t] == D[u[x,t],{x,2}], u[x,0] == f[x]}, u, {x,t}]}, "
               "Head[r] === List && FreeQ[r, DSolve] && !FreeQ[r, Integrate]]");
    /* the heat kernel (verified internally) solves the PDE */
    check_true("With[{g = 1/(2 Sqrt[Pi k t]) Exp[-(x-y)^2/(4 k t)]}, "
               "PossibleZeroQ[D[g,t] - k D[g,{x,2}]]]");
    /* declines the backward heat equation (k < 0, ill-posed) */
    check_form("Head[DSolve`HeatKernel[{D[u[x,t],t] == -D[u[x,t],{x,2}], u[x,0] == f[x]}, u, {x,t}]]",
               "DSolve`HeatKernel");
    /* declines a second-order-in-time (wave) equation with one condition */
    check_form("Head[DSolve`HeatKernel[{D[u[x,t],{t,2}] == D[u[x,t],{x,2}], u[x,0] == f[x]}, u, {x,t}]]",
               "DSolve`HeatKernel");
    /* declines advection-diffusion (a u_x term is present) */
    check_form("Head[DSolve`HeatKernel[{D[u[x,t],t] == D[u[x,t],{x,2}] + D[u[x,t],x], u[x,0] == f[x]}, u, {x,t}]]",
               "DSolve`HeatKernel");
}

/* Pinned system + PDE method builtins: each is REPL-callable as DSolve`<Name>[...]
 * (M8 systems, M6 PDE), verified by back-substitution, and declines a wrong-shape
 * input (no silent wrong answer). */
static void t_sys_pde_pinned_methods(void) {
    /* DecoupleSystem: independent equations */
    check_true("And @@ (PossibleZeroQ /@ ({y'[x] - y[x], z'[x] - 2 z[x]} /. "
               "DSolve`DecoupleSystem[{y'[x] == y[x], z'[x] == 2 z[x]}, {y, z}, x][[1]]))");
    /* TriangularSystem: DAG (z depends on y, y independent) */
    check_true("And @@ (PossibleZeroQ /@ ({y'[x], z'[x] + y[x]} /. "
               "DSolve`TriangularSystem[{y'[x] == 0, z'[x] + y[x] == 0}, {y, z}, x][[1]]))");
    /* LinearFirstOrderSystem: coupled constant matrix -> real Cos/Sin */
    check_true("And @@ (PossibleZeroQ /@ ({y'[x] - z[x], z'[x] + y[x]} /. "
               "DSolve`LinearFirstOrderSystem[{y'[x] == z[x], z'[x] == -y[x]}, {y, z}, x][[1]]))");
    /* PDELinearFirstOrder pinned, verified via a concrete C[1] */
    check_true("With[{uc = (u[x,y] /. DSolve`PDELinearFirstOrder[D[u[x,y],x] + 3 D[u[x,y],y] "
               "+ u[x,y] == 1, u, {x,y}][[1]]) /. C[1][z_] :> Sin[z]}, "
               "PossibleZeroQ[D[uc,x] + 3 D[uc,y] + uc - 1]]");
    /* the PDE method declines a scalar ODE (head stays symbolic) */
    check_form("Head[DSolve`PDELinearFirstOrder[y'[x] == y[x], y, x]]", "DSolve`PDELinearFirstOrder");
}

/* ---- M6: first-order nonlinear PDEs — PDEQuasilinear (Lagrange) ---- */
static void t_pde_quasilinear(void) {
    /* semilinear, variable coefficients: x u_x + y u_y == u -> x C[1][y/x]. */
    check_true("With[{uc = (u[x,y] /. DSolve[x D[u[x,y],x] + y D[u[x,y],y] == u[x,y], "
               "u, {x,y}][[1]]) /. C[1][z_] :> Sin[z]}, "
               "PossibleZeroQ[x D[uc,x] + y D[uc,y] - uc]]");
    /* semilinear needing y expressed along the characteristic: u_x + x u_y == y. */
    check_true("With[{uc = (u[x,y] /. DSolve[D[u[x,y],x] + x D[u[x,y],y] == y, "
               "u, {x,y}][[1]]) /. C[1][z_] :> Sin[z]}, "
               "PossibleZeroQ[D[uc,x] + x D[uc,y] - y]]");
    /* conservation law (inviscid Burgers) is IMPLICIT: {{ G == C[1][u] }}. */
    check_true("Head[DSolve[u[x,y] D[u[x,y],x] + D[u[x,y],y] == 0, u, {x,y}][[1,1]]] === Equal");
    /* verify the implicit relation by implicit differentiation (C[1] pinned to #^2):
     * with Psi(x,y,U) = G - U^2, u_x = -Psi_x/Psi_U, u_y = -Psi_y/Psi_U. */
    check_true("Module[{rel, Psi}, "
               "rel = DSolve[u[x,y] D[u[x,y],x] + D[u[x,y],y] == 0, u, {x,y}][[1,1]] /. C[1] -> (#^2 &); "
               "Psi = (rel[[1]] - rel[[2]]) /. u[x,y] -> U; "
               "PossibleZeroQ[(U (-D[Psi,x]/D[Psi,U]) + (-D[Psi,y]/D[Psi,U])) /. U -> u[x,y]]]");
    /* pinned builtin solves a constant-coefficient linear PDE too (a superset). */
    check_true("With[{uc = (u[x,y] /. DSolve`PDEQuasilinear[D[u[x,y],x] + 2 D[u[x,y],y] == 0, "
               "u, {x,y}][[1]]) /. C[1][z_] :> Cos[z]}, PossibleZeroQ[D[uc,x] + 2 D[uc,y]]]");
    /* declines a scalar ODE and a genuinely quasilinear non-conservation form. */
    check_form("Head[DSolve`PDEQuasilinear[y'[x] == y[x], y, x]]", "DSolve`PDEQuasilinear");
    check_form("Head[DSolve`PDEQuasilinear[u[x,y] D[u[x,y],x] + D[u[x,y],y] == u[x,y], u, {x,y}]]",
               "DSolve`PDEQuasilinear");
}

/* ---- M6: PDE Clairaut — complete integral + singular envelope ---- */
static void t_pde_clairaut(void) {
    /* u == x u_x + y u_y + u_x u_y -> complete integral C[1] C[2] + C[1] x + C[2] y. */
    check_true("With[{uc = u[x,y] /. DSolve[u[x,y] == x D[u[x,y],x] + y D[u[x,y],y] + "
               "D[u[x,y],x] D[u[x,y],y], u, {x,y}][[1]]}, "
               "PossibleZeroQ[uc - (x D[uc,x] + y D[uc,y] + D[uc,x] D[uc,y])]]");
    /* IncludeSingularSolutions adds the envelope (u = -x y); EVERY branch verifies. */
    check_true("And @@ Map[Function[br, With[{uc = u[x,y] /. br}, "
               "PossibleZeroQ[uc - (x D[uc,x] + y D[uc,y] + D[uc,x] D[uc,y])]]], "
               "DSolve[u[x,y] == x D[u[x,y],x] + y D[u[x,y],y] + D[u[x,y],x] D[u[x,y],y], "
               "u, {x,y}, IncludeSingularSolutions -> True]]");
    check_true("Length[DSolve[u[x,y] == x D[u[x,y],x] + y D[u[x,y],y] + D[u[x,y],x] D[u[x,y],y], "
               "u, {x,y}, IncludeSingularSolutions -> True]] == 2");
    /* pinned builtin with a different f (u_x^2), and declines a scalar ODE. */
    check_true("With[{uc = u[x,y] /. DSolve`PDEClairaut[u[x,y] == x D[u[x,y],x] + y D[u[x,y],y] "
               "+ D[u[x,y],x]^2, u, {x,y}][[1]]}, "
               "PossibleZeroQ[uc - (x D[uc,x] + y D[uc,y] + D[uc,x]^2)]]");
    check_form("Head[DSolve`PDEClairaut[y'[x] == y[x], y, x]]", "DSolve`PDEClairaut");
}

/* ---- M6: 2nd-order constant-coeff PDE with lower-order terms (operator
 * factoring → exponential-damped arbitrary functions) ---- */
static void t_pde2_lower_order(void) {
    /* distortionless telegraph u_tt − c² u_xx + a u_t + (a²/4) u == 0
     * → e^{−a t/2}(C[1][x+c t] + C[2][x−c t]) */
    check_true("With[{uc = (u[t,x] /. DSolve[D[u[t,x],{t,2}] == c^2 D[u[t,x],{x,2}] "
               "- a D[u[t,x],t] - (a^2/4) u[t,x], u, {t,x}][[1]]) /. {C[1] -> Sin, C[2] -> Cos}}, "
               "PossibleZeroQ[D[uc,{t,2}] - (c^2 D[uc,{x,2}] - a D[uc,t] - (a^2/4) uc)]]");
    /* pure-mixed with a lower-order term: u_xy + u_x == 0 → C[1][y] + e^{−y} C[2][x] */
    check_true("With[{uc = (u[x,y] /. DSolve[D[u[x,y],x,y] + D[u[x,y],x] == 0, u, {x,y}][[1]]) "
               "/. {C[1] -> Sin, C[2] -> Cos}}, PossibleZeroQ[D[uc,x,y] + D[uc,x]]]");
    /* a damped convection case (distinct λ, equal m): u_tt − u_xx + 2 u_t + u == 0 */
    check_true("With[{uc = (u[t,x] /. DSolve[D[u[t,x],{t,2}] - D[u[t,x],{x,2}] + 2 D[u[t,x],t] "
               "+ u[t,x] == 0, u, {t,x}][[1]]) /. {C[1] -> Sin, C[2] -> Cos}}, "
               "PossibleZeroQ[D[uc,{t,2}] - D[uc,{x,2}] + 2 D[uc,t] + uc]]");
    /* pinned builtin; regression: a principal-part-only equation is unchanged;
     * a non-factorable general telegraph declines (stays symbolic) */
    check_true("MatchQ[DSolve`PDELinearSecondOrder[D[u[t,x],{t,2}] == c^2 D[u[t,x],{x,2}] "
               "- a D[u[t,x],t] - (a^2/4) u[t,x], u, {t,x}], {{u -> _Function}}]");
    check_form("Head[DSolve[D[u[t,x],{t,2}] == D[u[t,x],{x,2}] - 3 D[u[t,x],t] - u[t,x], u, {t,x}]]",
               "DSolve");
}

/* ---- M6: Charpit's method — first-order fully nonlinear PDE, standard forms ---- */
static void t_pde_charpit(void) {
    /* Type I  F(p,q):  p^2 + q^2 == 1  (explicit, two branches, each verified) */
    check_true("And @@ Map[Function[br, With[{uc = u[x,y] /. br}, "
               "PossibleZeroQ[D[uc,x]^2 + D[uc,y]^2 - 1]]], "
               "DSolve[D[u[x,y],x]^2 + D[u[x,y],y]^2 == 1, u, {x,y}]]");
    /* Type III separable  p^2 - q^2 == x - y  (explicit) */
    check_true("With[{uc = u[x,y] /. DSolve[D[u[x,y],x]^2 - D[u[x,y],y]^2 == x - y, u, {x,y}][[1]]}, "
               "PossibleZeroQ[D[uc,x]^2 - D[uc,y]^2 - (x - y)]]");
    /* Type II  F(u,p,q):  p q == u  (implicit relation, verified by implicit diff) */
    check_true("Head[DSolve[D[u[x,y],x] D[u[x,y],y] == u[x,y], u, {x,y}][[1,1]]] === Equal");
    check_true("Module[{rel, Psi}, "
               "rel = DSolve[D[u[x,y],x] D[u[x,y],y] == u[x,y], u, {x,y}][[1,1]]; "
               "Psi = (rel[[1]] - rel[[2]]) /. u[x,y] -> U; "
               "PossibleZeroQ[((-D[Psi,x]/D[Psi,U]) (-D[Psi,y]/D[Psi,U]) - U) /. U -> u[x,y]]]");
    /* pinned builtin; declines a scalar ODE and a linear (non-Charpit) PDE */
    check_true("With[{uc = u[x,y] /. DSolve`PDECharpit[D[u[x,y],x]^2 + D[u[x,y],y]^2 == 4, u, {x,y}][[1]]}, "
               "PossibleZeroQ[D[uc,x]^2 + D[uc,y]^2 - 4]]");
    check_form("Head[DSolve`PDECharpit[y'[x] == y[x], y, x]]", "DSolve`PDECharpit");
    check_form("Head[DSolve`PDECharpit[D[u[x,y],x] + D[u[x,y],y] == 0, u, {x,y}]]",
               "DSolve`PDECharpit");
}

/* ---- DSolve is NOT HoldAll: an equation stored in a variable must solve ---- */
static void t_not_holdall(void) {
    check_true("FreeQ[Attributes[DSolve], HoldAll]");
    check_true("Module[{dseq = y'[x] + y[x] == a Sin[x]}, "
               "Head[DSolve[dseq, y, x]] === List]");
    check_true("Module[{dseq = y'[x] + y[x] == a Sin[x]}, "
               "PossibleZeroQ[(y'[x] + y[x] - a Sin[x]) /. DSolve[dseq, y, x][[1]]]]");
}

/* ---- unsupported equations stay symbolic (declined, not wrong) ---- */
static void t_declines_unsupported(void) {
    /* an irregular singular point (essential singularity at 0) is beyond the
     * series fallback and has no Liouvillian solution: DSolve stays symbolic.
     * (y'' + Sin[x] y == 0 now returns a power series — see t_powerseries_auto.) */
    check_form("Head[DSolve[y''[x] + Exp[1/x] y[x] == 0, y[x], x]]", "DSolve");
}

/* ---- M9: backfill unit + pinned-method coverage for thin methods ---- */

/* Homogeneous: was one auto-DSolve case; add a pinned-method test, two more
 * in-domain forms (y=v x separates and inverts), and an IVP. */
static void t_method_homogeneous(void) {
    check_true("PossibleZeroQ[(y'[x] - (x - y[x])/(x + y[x])) /. "
               "DSolve`Homogeneous[y'[x] == (x - y[x])/(x + y[x]), y, x][[1]]]");
}
static void t_homogeneous_more(void) {
    check_solves("y'[x] == (2 x - y[x])/(x + y[x])", "y'[x] - (2 x - y[x])/(x + y[x])");
    check_solves("y'[x] == y[x]/x + (y[x]/x)^2", "y'[x] - (y[x]/x + (y[x]/x)^2)");
}
static void t_ivp_homogeneous(void) {
    check_form("Head[DSolve[{y'[x] == (x - y[x])/(x + y[x]), y[1] == 1}, y, x]]", "List");
    check_true("PossibleZeroQ[(y'[x] - (x - y[x])/(x + y[x])) /. "
               "DSolve[{y'[x] == (x - y[x])/(x + y[x]), y[1] == 1}, y, x][[1]]]");
}

/* ReductionOfOrder: was one auto case; add a pinned-method test + two more
 * missing-y forms of the shape y'' == f(x) y'. */
static void t_method_reduce_order(void) {
    check_true("PossibleZeroQ[(y''[x] - y'[x]^2) /. "
               "DSolve`ReductionOfOrder[y''[x] == y'[x]^2, y, x][[1]]]");
    check_true("Not[FreeQ[DSolve`ReductionOfOrder[y''[x] == y'[x]^2, y, x][[1]], C[2]]]");
}
static void t_reduce_order_more(void) {
    check_method("DSolve`ReductionOfOrder", "y''[x] == y'[x]/x", "y''[x] - y'[x]/x");
    check_method("DSolve`ReductionOfOrder", "y''[x] == 2 y'[x]", "y''[x] - 2 y'[x]");
}

/* Clairaut: was general + singular; add a pinned-method test and three more
 * f(y') forms. */
static void t_method_clairaut(void) {
    check_true("PossibleZeroQ[(y[x] - x y'[x] - y'[x]^2) /. "
               "DSolve`Clairaut[y[x] == x y'[x] + y'[x]^2, y, x][[1]]]");
}
static void t_clairaut_more(void) {
    check_method("DSolve`Clairaut", "y[x] == x y'[x] + 1/y'[x]", "y[x] - x y'[x] - 1/y'[x]");
    check_method("DSolve`Clairaut", "y[x] == x y'[x] + Sqrt[1 + y'[x]^2]",
                 "y[x] - x y'[x] - Sqrt[1 + y'[x]^2]");
    check_method("DSolve`Clairaut", "y[x] == x y'[x] - Log[y'[x]]", "y[x] - x y'[x] + Log[y'[x]]");
}

/* ---- 1a: Lagrange / d'Alembert (parametric general solution) ---- */
/* Verify a parametric solution sol = {x->Function[{t},X], y->Function[{t},Y]}:
 * substitute x->X(t), y[x]->Y(t), y'[x]->Y'(t)/X'(t) into `resid` (the ODE lhs,
 * written in x, y[x], y'[x]) and require PossibleZeroQ. */
static void check_lagrange(const char* method, const char* eqn, const char* resid) {
    char buf[1200];
    snprintf(buf, sizeof(buf), "Head[%s[%s, y, x]]", method, eqn);
    check_form(buf, "List");
    snprintf(buf, sizeof(buf),
        "Module[{s = %s[%s, y, x][[1]], X, Y, yp}, "
        "X = (x /. s)[t]; Y = (y /. s)[t]; yp = D[Y,t]/D[X,t]; "
        "PossibleZeroQ[(%s) /. {Derivative[1][y][x] -> yp, y[x] -> Y, x -> X}]]",
        method, eqn, resid);
    check_true(buf);
}
static void t_method_lagrange(void) {
    /* y == 2 x y' + (y')^2 : phi=2p, psi=p^2 -> rational parametric (decidable) */
    check_lagrange("DSolve`Lagrange", "y[x] == 2 x y'[x] + (y'[x])^2",
                   "y[x] - (2 x y'[x] + (y'[x])^2)");
}
static void t_lagrange_more(void) {
    /* phi=2p, psi=p^3 */
    check_lagrange("DSolve`Lagrange", "y[x] == 2 x y'[x] + (y'[x])^3",
                   "y[x] - (2 x y'[x] + (y'[x])^3)");
    /* phi constant (=1), psi nonlinear (p^2): still Lagrange, transcendental (Log) */
    check_lagrange("DSolve`Lagrange", "y[x] == x + (y'[x])^2",
                   "y[x] - (x + (y'[x])^2)");
    /* automatic dispatch reaches it too */
    check_form("Head[DSolve[y[x] == 2 x y'[x] + (y'[x])^2, y, x]]", "List");
}
static void t_lagrange_declines(void) {
    /* Clairaut (phi==p) is owned by DSolve`Clairaut, not Lagrange */
    check_form("Head[DSolve`Lagrange[y[x] == x y'[x] + (y'[x])^2, y, x]]", "DSolve`Lagrange");
    /* genuinely linear in y' (phi const, psi affine): owned by LinearFirstOrder */
    check_form("Head[DSolve`Lagrange[y'[x] + y[x] == x, y, x]]", "DSolve`Lagrange");
    /* a parametric IVP is deferred: declines rather than ignoring the condition */
    check_form("Head[DSolve[{y[x] == 2 x y'[x] + (y'[x])^2, y[1] == 0}, y, x]]", "DSolve");
}
static void t_lagrange_singular(void) {
    /* y == x (y')^2 + (y')^3: phi(p)=p^2, roots of phi(p)=p are p=0,1 -> singular
     * lines y=0 and y=x+1, emitted alongside the parametric general branch. */
    check_true("Length[DSolve[y[x] == x (y'[x])^2 + (y'[x])^3, y, x, "
               "IncludeSingularSolutions -> True]] == 3");
    check_true("Module[{sing = Select[DSolve[y[x] == x (y'[x])^2 + (y'[x])^3, y, x, "
               "IncludeSingularSolutions -> True], Length[#] == 1 &]}, "
               "Length[sing] == 2 && And @@ (PossibleZeroQ /@ "
               "((y[x] - (x (y'[x])^2 + (y'[x])^3)) /. sing))]");
    /* default (no option): general parametric branch only */
    check_true("Length[DSolve[y[x] == x (y'[x])^2 + (y'[x])^3, y, x]] == 1");
}

/* PowerSeries: was ordinary + auto; add two more ordinary-point forms. */
static void t_powerseries_more(void) {
    check_series("DSolve`PowerSeries", "y''[x] + x y'[x] + y[x] == 0",
                 "D[b,{x,2}] + x D[b,x] + b");
    check_series("DSolve`PowerSeries", "y''[x] - x y[x] == 0", "D[b,{x,2}] - x b");
}

/* Homogeneous now inverts the pure-log (algebraic) family via exponentiation:
 * the exponentiated relation is algebraic and Solve returns Root branches (see
 * dsolve_homogeneous.c homog_exp_log_invert).  The transcendental (ArcTan) family
 * still declines — it has no explicit inverse. */
static void t_homogeneous_algebraic(void) {
    check_method("DSolve`Homogeneous", "y'[x] == (x + 2 y[x])/(2 x + y[x])",
                 "y'[x] - (x + 2 y[x])/(2 x + y[x])");
    check_method("DSolve`Homogeneous", "y'[x] == (2 x + y[x])/(x + 2 y[x])",
                 "y'[x] - (2 x + y[x])/(x + 2 y[x])");
    /* the transcendental case (x+y)/(x-y) has no explicit inverse; it is returned
     * as an implicit first integral instead of declining — see t_homogeneous_implicit. */
}

/* ReductionOfOrder now accepts a correct-but-unsimplified antiderivative (its
 * D[yint]-p guard falls back to PossibleZeroQ), so the autonomous a+b(y')^2 and
 * the Riccati-in-p c x (y')^2 families solve. */
static void t_reduce_order_riccati(void) {
    check_method("DSolve`ReductionOfOrder", "y''[x] == 1 + y'[x]^2", "y''[x] - 1 - y'[x]^2");
    check_method("DSolve`ReductionOfOrder", "y''[x] == -2 x y'[x]^2", "y''[x] + 2 x y'[x]^2");
}

/* The transcendental (ArcTan log-spiral) homogeneous family has no explicit
 * inverse and is returned as the implicit first integral G(x,y[x]) == C[1] (an
 * Equal, not a y[x] -> rule).  Verify by implicit differentiation: from
 * d/dx[G == C] the ODE forces G_x + G_y y'[x] == 0, so substituting the ODE RHS
 * for y'[x] must vanish. */
static void check_implicit(const char* rhs) {
    char buf[768];
    snprintf(buf, sizeof(buf), "Head[DSolve[y'[x] == %s, y, x][[1,1]]] === Equal", rhs);
    check_true(buf);
    snprintf(buf, sizeof(buf),
        "PossibleZeroQ[Module[{eq = DSolve[y'[x] == %s, y, x][[1,1]]}, "
        "D[eq[[1]] - eq[[2]], x] /. y'[x] -> (%s)]]", rhs, rhs);
    check_true(buf);
}
/* Like check_implicit but for a PINNED method: method[y'==rhs,y,x] returns the
 * implicit first integral {{G==C[1]}}; verify by implicit differentiation. */
static void check_pinned_implicit(const char* method, const char* rhs) {
    char buf[900];
    snprintf(buf, sizeof(buf), "Head[%s[y'[x] == %s, y, x][[1,1]]] === Equal", method, rhs);
    check_true(buf);
    snprintf(buf, sizeof(buf),
        "PossibleZeroQ[Module[{eq = %s[y'[x] == %s, y, x][[1,1]]}, "
        "D[eq[[1]] - eq[[2]], x] /. y'[x] -> (%s)]]", method, rhs, rhs);
    check_true(buf);
}
static void t_homogeneous_implicit(void) {
    check_implicit("(x + y[x])/(x - y[x])");
    check_implicit("(x + y[x])/(2 x + y[x])");
    check_implicit("(3 x + y[x])/(x + 2 y[x])");
    /* IVP fits the constant: the relation passes through (1, 0) */
    check_true("PossibleZeroQ[Module[{eq = DSolve[{y'[x] == (x + y[x])/(x - y[x]), y[1] == 0}, "
               "y, x][[1,1]]}, (eq[[1]] - eq[[2]]) /. {x -> 1, y[x] -> 0}]]");
    /* the pinned method also returns the implicit form when no explicit inverse exists */
    check_true("Head[DSolve`Homogeneous[y'[x] == (x + y[x])/(x - y[x]), y, x][[1,1]]] === Equal");
}

/* ---- 1a: Chini / Abel (implicit first integral, reducible-to-autonomous) ---- */
static void t_method_chini(void) {
    /* f=x^2, n=3, B=0, C=1 -> u'=u^3+1, u = x y ; implicit first integral */
    check_implicit("x^2 y[x]^3 - y[x]/x + 1/x");
    check_true("Head[DSolve`Chini[y'[x] == x^2 y[x]^3 - y[x]/x + 1/x, y, x][[1,1]]] === Equal");
}
static void t_chini_more(void) {
    /* n=4: f=x^3 -> f^(1/3)=x, u'=u^4+1 */
    check_implicit("x^3 y[x]^4 - y[x]/x + 1/x");
    /* radical reduction (f=x, n=3 -> Sqrt[x]) still verifies */
    check_implicit("x y[x]^3 - (1/(2 x)) y[x] + 1/Sqrt[x]");
    /* declines: non-reducible (h not matched), and Riccati n=2 */
    check_form("Head[DSolve`Chini[y'[x] == x y[x]^3 + y[x] + x, y, x]]", "DSolve`Chini");
    check_form("Head[DSolve`Chini[y'[x] == y[x]^2 + x, y, x]]", "DSolve`Chini");
}
static void t_method_abel(void) {
    /* the f=x^2 Chini above, shifted by z = y + 1, introduces the y^2 term */
    check_implicit("x^2 y[x]^3 + 3 x^2 y[x]^2 + (3 x^2 - 1/x) y[x] + x^2");
    check_true("Head[DSolve`Abel[y'[x] == x^2 y[x]^3 + 3 x^2 y[x]^2 + (3 x^2 - 1/x) y[x] + x^2, "
               "y, x][[1,1]]] === Equal");
    /* Abel declines a Chini (f2 == 0) — DSolve`Chini owns that */
    check_form("Head[DSolve`Abel[y'[x] == x^2 y[x]^3 - y[x]/x + 1/x, y, x]]", "DSolve`Abel");
}

/* ---- 1a: Lie point-symmetry (heuristic; M10 L1: abaco1_simple) ---- */
/* Like check_implicit but pins DSolve`LieSymmetry: the abaco1_simple ansatze
 * overlap linear/separable, which the automatic cascade claims first, so the
 * method must be exercised through its own builtin. */
static void check_lie_implicit(const char* rhs) {
    char buf[768];
    snprintf(buf, sizeof(buf),
             "Head[DSolve`LieSymmetry[y'[x] == %s, y, x][[1,1]]] === Equal", rhs);
    check_true(buf);
    snprintf(buf, sizeof(buf),
        "PossibleZeroQ[Module[{eq = DSolve`LieSymmetry[y'[x] == %s, y, x][[1,1]]}, "
        "D[eq[[1]] - eq[[2]], x] /. y'[x] -> (%s)]]", rhs, rhs);
    check_true(buf);
}
static void t_method_lie(void) {
    /* case A (linear): xi=0, eta=Exp[Integrate[omega_y, x]] */
    check_lie_implicit("y[x]/x + x");
    /* case B/C (separable): omega_x/omega or omega_y/omega free of the other var */
    check_lie_implicit("x y[x]^2");
    /* another linear, non-elementary-looking integrating factor Exp[-x^2] */
    check_lie_implicit("2 x y[x] + x");
    /* alias resolves to the same method */
    check_true("Head[DSolve`LieGroup[y'[x] == y[x]/x + x, y, x][[1,1]]] === Equal");
}
static void t_lie_declines(void) {
    /* Riccati y' == y^2 + x has no abaco1_simple / affine symmetry: the method
     * declines (no wrong answer, head stays symbolic). */
    check_form("Head[DSolve`LieSymmetry[y'[x] == y[x]^2 + x, y, x]]", "DSolve`LieSymmetry");
}
/* Regression: an omega carrying an UNDEFINED function of both variables used to hang
 * the quadrature heuristics (the classifier / free-of / zero-test ops balloon on the
 * transcendental derivatives of an arbitrary function).  It must now decline promptly:
 * the rational/algebraic ansatze are skipped, abaco2_unique_* reads the symmetry off
 * R = M_y/M_x but the Lie quadrature is non-elementary, so the whole method declines
 * with no inert head.  (If this regresses it HANGS — a ctest timeout catches it.) */
static void t_lie_undefined_function_declines(void) {
    /* user-reported form: no elementary [G(y),F(x)] symmetry -> declines */
    check_form("Head[DSolve`LieSymmetry[y'[x] == Tan[ArcTan[y[x]] + F[x^2 + y[x]^2]], "
               "y, x]]", "DSolve`LieSymmetry");
    /* Cheb-Terrab & Roche Eq 70: symmetry [y,-x] exists but its quadrature
     * int 1/(a^2 - sin a) da is non-elementary, so the method declines (no inert
     * integral head), matching the documented arbitrary-function policy */
    check_form("Head[DSolve`LieSymmetry[y'[x] == -Tan[ArcTan[x/y[x]] + H[x^2 + y[x]^2]], "
               "y, x]]", "DSolve`LieSymmetry");
    /* full cascade also declines cleanly (head stays DSolve) */
    check_form("Head[DSolve[y'[x] == Tan[ArcTan[y[x]] + F[x^2 + y[x]^2]], y[x], x]]",
               "DSolve");
}
/* L2 `linear` heuristic (affine symmetry): the linear-coefficients class
 * y' == (a1 x + b1 y + c1)/(a2 x + b2 y + c2).  The deterministic
 * DSolve`LinearCoefficients (M9) now claims this class in the AUTOMATIC cascade
 * with an explicit (Root-form) solution, so exercise Lie's `linear` heuristic
 * through its pinned builtin, which returns the implicit first integral. */
static void t_lie_linear_coefficients(void) {
    check_pinned_implicit("DSolve`LieSymmetry", "(x + 2 y[x] - 4)/(2 x + y[x] - 5)");
    check_pinned_implicit("DSolve`LieSymmetry", "(2 x + 3 y[x] - 1)/(3 x + 2 y[x] + 2)");
}
/* L3 `bivariate` heuristic: a genuinely degree-2 polynomial symmetry
 * (xi = x^2, eta = x y) where NO affine (`linear`) or one-variable
 * (`abaco1_simple`) symmetry exists, so `bivariate` is the ONLY heuristic that
 * can solve these — not a vacuous pass through an earlier one.  (Verified in the
 * REPL: the degree-1 determining NullSpace is trivial and every abaco1_simple
 * ratio depends on both x and y.)  Both ODEs are members of the family
 * omega = y/x + A(y/x)/x, whose symmetry is (x^2, x y). */
static void t_lie_bivariate(void) {
    check_pinned_implicit("DSolve`LieSymmetry", "-1/x + y[x]/x + y[x]^2/x^3");
    check_pinned_implicit("DSolve`LieSymmetry", "y[x]/x^2 + y[x]/x + y[x]^2/x^3");
    /* alias resolves to the same method on a bivariate-only ODE */
    check_true("Head[DSolve`LieGroup[y'[x] == -1/x + y[x]/x + y[x]^2/x^3, y, x]"
               "[[1,1]]] === Equal");
}
/* `abaco1_product` (Cheb-Terrab & Roche 1998, §4.1): the symmetry [F(x) G(y), 0]
 * (and its inverse [0, F(x) G(y)]) — a rational-but-non-polynomial infinitesimal
 * (here xi = y/x) that the polynomial `linear`/`bivariate` ansatze miss.  These are
 * members of the invariant family omega = f_x/(g(y) f(x) + J(y)); the inhomogeneous
 * J(y) breaks the scaling symmetry so `abaco1_simple`/`linear`/`bivariate` all
 * decline (verified in the REPL: the full cascade returns unevaluated on each), and
 * abaco1_product — last in the chain — is the only heuristic that solves them. */
static void t_lie_abaco1_product(void) {
    /* direct [F(x) G(y), 0]: xi = y/x, omega = 2 x y/(x^2 + P(y)) from f = x^2/2,
     * g = 1/y (symmetry is independent of the invariant family's J(y)).  The chosen
     * y^4 forms of P admit NO polynomial symmetry of degree <= 3 and no one-variable
     * symmetry (verified in the REPL: deg-1/2/3 determining NullSpace empty, every
     * abaco1_simple ratio depends on both x and y), so `linear`/`bivariate`/
     * `abaco1_simple` all decline and abaco1_product is the only heuristic that
     * solves them. */
    check_pinned_implicit("DSolve`LieSymmetry", "2 x y[x]/(x^2 + 2 y[x]^4 + 2)");
    check_pinned_implicit("DSolve`LieSymmetry", "2 x y[x]/(x^2 + 2 y[x]^4 - 3)");
    /* inverse pattern [0, F(x) G(y)] via the inverse ODE (also fully isolating) */
    check_pinned_implicit("DSolve`LieSymmetry", "(y[x]^2 + 2 x^4 + 2)/(2 x y[x])");
    /* alias resolves to the same method */
    check_true("Head[DSolve`LieGroup[y'[x] == 2 x y[x]/(x^2 + 2 y[x]^4 + 2), y, x]"
               "[[1,1]]] === Equal");
}
/* `function_sum` (Cheb-Terrab & Roche §4.2): the additive symmetry [F(x)+G(y), 0].
 * These omega are members of the §4.2 invariant family (F=1/x, G=y with J=0 and J=1/y),
 * whose 1/omega is transcendental (Log) so every rational Lie heuristic
 * (abaco1_simple/linear/abaco1_product) declines and function_sum is the first (and
 * only) to solve them — verified in the REPL via per-heuristic attribution, so the pass
 * is not vacuous.  The classifying quantity is the *rational* factor
 * omega.d^2/dx^2(1/omega) = F''/(F+G); the leading omega cancels the transcendental
 * part.  Verified by implicit differentiation.  (The F=1/x^2 ArcTan member also solves
 * and verifies in the REPL but its verify is ~7 s, so it is omitted from the suites.) */
static void t_lie_function_sum(void) {
    check_pinned_implicit("DSolve`LieSymmetry",
        "(x y[x]^3)/(-1 + x y[x] + x^2 y[x]^2 - 2 Log[1 + x y[x]] "
        "- 2 x y[x] Log[1 + x y[x]])");
    check_pinned_implicit("DSolve`LieSymmetry",
        "(x y[x]^3)/(-1 + y[x]^2 + x y[x]^3 + x y[x] + x^2 y[x]^2 "
        "- 2 Log[1 + x y[x]] - 2 x y[x] Log[1 + x y[x]])");
    /* alias resolves to the same method */
    check_true("Head[DSolve`LieGroup[y'[x] == (x y[x]^3)/(-1 + x y[x] + x^2 y[x]^2 "
               "- 2 Log[1 + x y[x]] - 2 x y[x] Log[1 + x y[x]]), y, x][[1,1]]] === Equal");
}
/* `abaco2_unique_unknown` (Cheb-Terrab & Roche §4.4.1): the symmetries [F(x),G(y)] /
 * [G(y),F(x)], found from a non-integer power (here (x^2+y^2)^p) of both variables in
 * omega.  For omega = (x/y)(x^2+y^2)^p the mapping M = (x^2+y^2)^p gives R = M_y/M_x =
 * y/x (the p and the power cancel), x-factor 1/x, and the symmetry [1/x, -1/y].  The
 * irrational omega makes every rational heuristic decline; abaco2_similar declines too
 * (verified via attribution: abaco2_unique_unknown fires), so the pass is not vacuous.
 * Verified by implicit differentiation. */
static void t_lie_abaco2_unique_unknown(void) {
    check_pinned_implicit("DSolve`LieSymmetry", "(x/y[x]) (x^2 + y[x]^2)^(1/3)");
    check_pinned_implicit("DSolve`LieSymmetry", "(x/y[x]) Sqrt[x^2 + y[x]^2]");
    check_pinned_implicit("DSolve`LieSymmetry", "(x/y[x]) (2 x^2 + y[x]^2)^(1/3)");
    /* alias resolves to the same method */
    check_true("Head[DSolve`LieGroup[y'[x] == (x/y[x]) (x^2 + y[x]^2)^(1/3), y, x]"
               "[[1,1]]] === Equal");
}
/* `abaco2_unique_unknown` §4.4.1 "differential invariant of order zero" extension
 * (Cheb-Terrab & Roche Eqs 73-81): the order-zero candidates [-R,1] / [1,-R] / [1,-1/R]
 * built directly from R = M_y/M_x without a separability test.  Kamke's first order
 * ODE 433, (x y' + y + 2x)^2 == 4(x y + x^2 + a), has the mapping M = Sqrt[x y + x^2 + a]
 * whose ratio R = x/(2x+y) does NOT separate by product, so only the order-zero
 * candidate [1, -R] finds the symmetry -> first integral x - Sqrt[x^2 + x y + a] == C[1]
 * (verified by implicit differentiation on the isolated y'-branch). */
static void t_lie_abaco2_order_zero(void) {
    check_pinned_implicit("DSolve`LieSymmetry",
        "(-y[x] - 2 x + 2 Sqrt[x y[x] + x^2 + 1])/x");
    check_pinned_implicit("DSolve`LieSymmetry",
        "(-y[x] - 2 x - 2 Sqrt[x y[x] + x^2 + 1])/x");   /* the other branch */
}
/* `chi` (Cheb-Terrab, Duarte & da Mota, CPC 101 1997, 5th algorithm): the
 * eta = xi omega + chi reformulation, with chi from a rich-basis (transcendental atoms
 * of omega) determining system — the one heuristic whose chi may be a genuine
 * transcendental beyond `bivariate`'s polynomial reach.  Kamke's first order ODE 357,
 * x y' ln(x) sin(y) + cos(y)(1 - x cos(y)) == 0, has the symmetry
 * [0, cos(y)^2/(ln(x) sin(y))] -> first integral -x + Log[x] Sec[y[x]] == C[1] (which
 * no earlier heuristic finds: omega is trig, so the rational/algebraic ansatze are
 * skipped and there is no [F(x),G(y)] kernel).  Verified by implicit differentiation. */
static void t_lie_chi(void) {
    check_pinned_implicit("DSolve`LieSymmetry",
        "-Cos[y[x]] (1 - x Cos[y[x]])/(x Log[x] Sin[y[x]])");
    /* alias resolves to the same method */
    check_true("Head[DSolve`LieGroup[y'[x] == -Cos[y[x]] (1 - x Cos[y[x]])/"
               "(x Log[x] Sin[y[x]]), y, x][[1,1]]] === Equal");
}
/* Chini reduction (b): linear-term removal y = e^(int g) w -> separable, for the
 * sub-class where reduction (a) (B,C constant) fails.  y' == x E^(2x) y^3 - y -
 * x E^-x -> w' == x(w^3 - 1) via y = E^-x w; implicit first integral. */
static void t_chini_linremoval(void) {
    check_implicit("-x Exp[-x] - y[x] + x Exp[2 x] y[x]^3");
    check_true("Head[DSolve`Chini[y'[x] == -x Exp[-x] - y[x] + x Exp[2 x] y[x]^3, "
               "y, x][[1,1]]] === Equal");
    /* the pure-Bernoulli variant with the E^(2x) coefficient (previously mis-
     * classified by the pre-simplify guard) now solves */
    check_form("Head[DSolve[y'[x] == x Exp[2 x] y[x]^3 - y[x], y, x]]", "List");
}
/* Exact via mu = x^a y^b integrating factor (constant exponents): the equation
 * (x y - 2 x) y' == y - y^2 + 3 x^2 y^3 is exact under mu = x^-2 y^-3. */
static void t_exact_xayb(void) {
    check_form("Head[DSolve[(x y[x] - 2 x) y'[x] == y[x] - y[x]^2 + 3 x^2 y[x]^3, y, x]]", "List");
    /* Function-form solve so y'[x] is rewritten; residual back-substitutes to 0 */
    check_true("PossibleZeroQ[((x y[x] - 2 x) y'[x] - (y[x] - y[x]^2 + 3 x^2 y[x]^3)) /. "
               "DSolve[(x y[x] - 2 x) y'[x] == y[x] - y[x]^2 + 3 x^2 y[x]^3, y, x][[1]]]");
    check_form("Head[DSolve`Exact[(x y[x] - 2 x) y'[x] == y[x] - y[x]^2 + 3 x^2 y[x]^3, y, x]]", "List");
    /* a non-x^a y^b equation still declines the exact method (no wrong answer) */
    check_form("Head[DSolve`Exact[y'[x] == Sqrt[y[x]^4 + 1], y, x]]", "DSolve`Exact");
}

/* ---- M9: SymPy deterministic parity gaps ---- */

/* Every branch of a (possibly multi-branch) result back-substitutes to zero.
 * The second DSolve argument is `y` (Function form) so y'[x] substitutes too. */
static void check_all_branches(const char* call, const char* resid) {
    char buf[1024];
    snprintf(buf, sizeof(buf), "Head[%s] === List", call);
    check_true(buf);                                /* non-vacuous: it actually solved */
    snprintf(buf, sizeof(buf),
             "And @@ Map[PossibleZeroQ[(%s) /. #] &, %s]", resid, call);
    check_true(buf);
}

/* Factorable: (y'-y)(y'+y)==0 splits into y'==y and y'==-y. */
static void t_method_factorable(void) {
    check_all_branches(
        "DSolve`Factorable[(y'[x] - y[x]) (y'[x] + y[x]) == 0, y, x]",
        "(y'[x] - y[x]) (y'[x] + y[x])");
    /* automatic cascade claims it too (Factorable runs at the front) */
    check_all_branches(
        "DSolve[(y'[x] - y[x]) (y'[x] + y[x]) == 0, y, x]",
        "(y'[x] - y[x]) (y'[x] + y[x])");
}
static void t_factorable_more(void) {
    /* three linear factors, distinct spectra */
    check_all_branches(
        "DSolve`Factorable[(y'[x] - y[x]) (y'[x] - 2 y[x]) (y'[x] + y[x]) == 0, y, x]",
        "(y'[x] - y[x]) (y'[x] - 2 y[x]) (y'[x] + y[x])");
    /* a product mixing a linear-inhomogeneous factor */
    check_all_branches(
        "DSolve`Factorable[(y'[x] - x) (y'[x] + y[x]) == 0, y, x]",
        "(y'[x] - x) (y'[x] + y[x])");
}
static void t_factorable_declines(void) {
    /* irreducible (a single differential factor) — Factorable declines */
    check_form("Head[DSolve`Factorable[y'[x] + y[x] == 0, y, x]]", "DSolve`Factorable");
    /* a Sqrt[y] coefficient is non-polynomial in the funcapps: declines, no hang */
    check_form("Head[DSolve`Factorable[y'[x] == 2 Sqrt[y[x]], y, x]]", "DSolve`Factorable");
    /* a pure-function factor (y itself) is not a differential factor: p p' == p^2/y
     * has one genuine factor, so Factorable declines (AutonomousReduction owns it) */
    check_form("Head[DSolve`Factorable[p[y] p'[y] == p[y]^2/y, p, y]]", "DSolve`Factorable");
}

/* NthAlgebraic: algebraic of degree >= 2 in the top derivative. */
static void t_method_nth_algebraic(void) {
    /* (y')^2 == 4 y  ->  y == (x + C)^2 (two sign branches) */
    check_all_branches(
        "DSolve`NthAlgebraic[(y'[x])^2 == 4 y[x], y, x]",
        "(y'[x])^2 - 4 y[x]");
}
static void t_nth_algebraic_more(void) {
    /* branches free of y integrate directly (Quadrature): (y')^2 == 2
     * -> y' == +/- Sqrt[2] -> y == +/- Sqrt[2] x + C[1] (irreducible over Q, so
     * this is NthAlgebraic's alone, not Factorable's) */
    check_all_branches(
        "DSolve`NthAlgebraic[(y'[x])^2 == 2, y, x]",
        "(y'[x])^2 - 2");
    /* automatic cascade: (y')^2 == 4 y is claimed by NthAlgebraic at the front */
    check_all_branches(
        "DSolve[(y'[x])^2 == 4 y[x], y, x]",
        "(y'[x])^2 - 4 y[x]");
}
static void t_nth_algebraic_declines(void) {
    /* linear in the top derivative (the normal case) -> the specialists own it */
    check_form("Head[DSolve`NthAlgebraic[y''[x] + y[x] == 0, y, x]]", "DSolve`NthAlgebraic");
    check_form("Head[DSolve`NthAlgebraic[y'[x] + y[x] == 0, y, x]]", "DSolve`NthAlgebraic");
}

/* LinearCoefficients: y' == (a1 x+b1 y+c1)/(a2 x+b2 y+c2). */
static void t_method_lincoeff(void) {
    /* det != 0 -> explicit (Root-form) branches; verify against the cleared eqn */
    check_all_branches(
        "DSolve`LinearCoefficients[y'[x] == (x + 2 y[x] - 4)/(2 x + y[x] - 5), y, x]",
        "y'[x] (2 x + y[x] - 5) - (x + 2 y[x] - 4)");
    /* det == 0 (parallel) -> separable, implicit first integral */
    check_pinned_implicit("DSolve`LinearCoefficients",
                          "(x + y[x] + 1)/(2 x + 2 y[x] - 1)");
}
static void t_lincoeff_more(void) {
    /* a second det != 0 example (previously unsolved by any deterministic method) */
    check_all_branches(
        "DSolve`LinearCoefficients[y'[x] == (2 x + 3 y[x] - 1)/(3 x + 2 y[x] + 2), y, x]",
        "y'[x] (3 x + 2 y[x] + 2) - (2 x + 3 y[x] - 1)");
    /* det != 0 log-spiral (no explicit inverse) -> implicit first integral */
    check_pinned_implicit("DSolve`LinearCoefficients",
                          "(x + y[x] + 1)/(x - y[x] + 3)");
    /* the automatic cascade also solves it */
    check_form("Head[DSolve[y'[x] == (x + 2 y[x] - 4)/(2 x + y[x] - 5), y, x]]", "List");
}
static void t_lincoeff_declines(void) {
    /* nonlinear (y^2) -> not a ratio of affine forms */
    check_form("Head[DSolve`LinearCoefficients[y'[x] == y[x]^2 + x, y, x]]",
               "DSolve`LinearCoefficients");
    /* no y-coupling (a pure quadrature) */
    check_form("Head[DSolve`LinearCoefficients[y'[x] == x, y, x]]",
               "DSolve`LinearCoefficients");
}

/* AlmostLinear: f(x)g(y) y' + k(x)l(y) + m(x) == 0.  2 y y' + y^2 - x == 0
 * (u = y^2 -> u' + u == x) is the flagship. */
static void t_method_almostlinear(void) {
    check_all_branches("DSolve`AlmostLinear[2 y[x] y'[x] + y[x]^2 - x == 0, y, x]",
                       "2 y[x] y'[x] + y[x]^2 - x");
    /* automatic cascade also solves it */
    check_all_branches("DSolve[2 y[x] y'[x] + y[x]^2 - x == 0, y, x]",
                       "2 y[x] y'[x] + y[x]^2 - x");
}
static void t_almostlinear_declines(void) {
    /* nonlinear in y' -> not almost-linear */
    check_form("Head[DSolve`AlmostLinear[y'[x]^2 == y[x], y, x]]", "DSolve`AlmostLinear");
}

/* SeparableReduced: x y'/y == G(x^n y).  y' == y^2/(1 + x y) -> w = x y,
 * implicit first integral (verified by implicit differentiation). */
static void t_method_sepreduced(void) {
    check_pinned_implicit("DSolve`SeparableReduced", "y[x]^2/(1 + x y[x])");
}
static void t_sepreduced_declines(void) {
    /* x r_x/(y r_y) is not constant -> not the x^n y form */
    check_form("Head[DSolve`SeparableReduced[y'[x] == x^2 + y[x], y, x]]",
               "DSolve`SeparableReduced");
}

/* Liouville: y'' + g(y)(y')^2 + h(x)y' == 0 (both y and x present). */
static void t_method_liouville(void) {
    /* g = 1/y, h = 1/x  ->  Exp[G]=y, EG=y^2/2; EH=Log[x]; y^2/2 == C[1]Log[x]+C[2] */
    check_all_branches(
        "DSolve`Liouville[y''[x] + (y'[x])^2/y[x] + y'[x]/x == 0, y, x]",
        "y''[x] + (y'[x])^2/y[x] + y'[x]/x");
    /* automatic cascade also solves it (missing-y/missing-x reductions decline) */
    check_all_branches(
        "DSolve[y''[x] + (y'[x])^2/y[x] + y'[x]/x == 0, y, x]",
        "y''[x] + (y'[x])^2/y[x] + y'[x]/x");
}
static void t_liouville_declines(void) {
    /* linear -> not a Liouville nonlinearity */
    check_form("Head[DSolve`Liouville[y''[x] + y[x] == 0, y, x]]", "DSolve`Liouville");
    /* g == 0 (no (y')^2 term) is ReductionOfOrder's, not Liouville's */
    check_form("Head[DSolve`Liouville[y''[x] + y'[x]/x == 0, y, x]]", "DSolve`Liouville");
}

/* UndeterminedCoefficients: tidy particular for UC forcing of a const-coeff ODE. */
static void t_method_undetcoeff(void) {
    check_all_branches("DSolve`UndeterminedCoefficients[y''[x] + y[x] == x^2, y, x]",
                       "y''[x] + y[x] - x^2");
    check_all_branches("DSolve`UndeterminedCoefficients[y''[x] + y[x] == Sin[2 x], y, x]",
                       "y''[x] + y[x] - Sin[2 x]");
    /* resonance: forcing coincides with a homogeneous mode -> x-multiplied trial */
    check_all_branches("DSolve`UndeterminedCoefficients[y''[x] - y[x] == Exp[x], y, x]",
                       "y''[x] - y[x] - Exp[x]");
    check_all_branches("DSolve`UndeterminedCoefficients[y''[x] + y[x] == Cos[x], y, x]",
                       "y''[x] + y[x] - Cos[x]");
    /* sum of terms (superposition) */
    check_all_branches("DSolve`UndeterminedCoefficients[y''[x] + y[x] == x + Exp[2 x], y, x]",
                       "y''[x] + y[x] - x - Exp[2 x]");
    /* automatic cascade claims it too (before constcoeff) */
    check_all_branches("DSolve[y''[x] + y[x] == x^2, y, x]", "y''[x] + y[x] - x^2");
}
static void t_undetcoeff_declines(void) {
    /* variable coefficients -> not this method (Euler/other own it) */
    check_form("Head[DSolve`UndeterminedCoefficients[x y''[x] + y[x] == x, y, x]]",
               "DSolve`UndeterminedCoefficients");
    /* non-UC forcing (Log) -> declines; constcoeff's var-params handles it */
    check_form("Head[DSolve`UndeterminedCoefficients[y''[x] + y[x] == Log[x], y, x]]",
               "DSolve`UndeterminedCoefficients");
    /* homogeneous (g == 0) is left to LinearConstantCoefficients */
    check_form("Head[DSolve`UndeterminedCoefficients[y''[x] + y[x] == 0, y, x]]",
               "DSolve`UndeterminedCoefficients");
}

/* FirstOrderPowerSeries: y' == F(x, y) about x0 = 0, truncated SeriesData.
 * The truncated residual is O[x]^N, so Normal[residual] == 0. */
static void t_method_first_order_series(void) {
    check_true("Head[DSolve`FirstOrderPowerSeries[y'[x] == x + y[x], y, x]] === List");
    check_true("PossibleZeroQ[Normal[(y'[x] - x - y[x]) /. "
               "DSolve`FirstOrderPowerSeries[y'[x] == x + y[x], y, x][[1]]]]");
    /* nonlinear -> genuinely new coverage (no closed form) */
    check_true("PossibleZeroQ[Normal[(y'[x] - x - y[x]^2) /. "
               "DSolve`FirstOrderPowerSeries[y'[x] == x + y[x]^2, y, x][[1]]]]");
    /* pinned-only: the automatic cascade does NOT auto-apply it, so a first-order
     * ODE with no closed form stays unevaluated (matching Mathematica / SymPy).
     * y' == x^2 + y^3 has no elementary/Lie closed form (the whole cascade declines)
     * yet has an ordinary point at 0, so FirstOrderPowerSeries alone closes it.
     * (Sqrt[x+y] used to serve here but now solves via Lie/abaco2_similar.) */
    check_true("Head[DSolve`FirstOrderPowerSeries[y'[x] == x^2 + y[x]^3, y, x]] === List");
    check_form("Head[DSolve[y'[x] == x^2 + y[x]^3, y[x], x]]", "DSolve");
}
static void t_first_order_series_declines(void) {
    /* x0 = 0 is a pole of F (not ordinary) -> declines */
    check_form("Head[DSolve`FirstOrderPowerSeries[y'[x] == y[x]/x, y, x]]",
               "DSolve`FirstOrderPowerSeries");
}

/* Regression: heap double-free in the Risch-Norman monomial enumerator
 * (intrischnorman.c enumerate_monomials cap-cleanup) reached via the Lie
 * first-integral quadrature. Pre-fix these corrupted the heap and crashed.
 * The assertion is that evaluation completes (a crash aborts the binary) and
 * returns either a solution (_List) or a clean decline (_DSolve); a generalized
 * family exercises the same path so the guard is not overfit to one input. */
/* M34 (§2.2.14) — the four fixes' flagship corpus cases; forward-generator grids
 * live in test_dsolve_m34_stress.c.  Each wraps DSolve in TimeConstrained so a
 * regression that re-introduces a hang FAILS (=== $Aborted) rather than stalls. */
static void t_m34_corpus_cases(void) {
    /* 1337 — VoP verify short-circuit: the -Cos Log[Sec+Tan] answer's residual
     * spins zero_test; the numeric-zero verify keeps it. */
    check_true("Abs[N[(y''[x] + y[x] - Tan[x]) /. TimeConstrained[DSolve[y''[x] + y[x] == "
               "Tan[x], y, x], 8, $Aborted][[1]] /. {C[1] -> 13/10, C[2] -> 7/10, x -> 1/2}, 20]] < 10^-6");
    /* 1350 — forced Bessel operator, arbitrary g: inert-Integrate VoP integral form. */
    check_true("Head[TimeConstrained[DSolve[x^2 y''[x] + x y'[x] + (x^2 - 1/4) y[x] == g[x], "
               "y, x], 8, $Aborted]] === List");
    /* 1384 — 2nd-order exact -> Frobenius series IVP, fit at x=0. */
    check_true("FreeQ[TimeConstrained[DSolve[{y''[x] + Sin[x] y'[x] + Cos[x] y[x] == 0, y[0] == 0, "
               "y'[0] == 1}, y, x], 8, $Aborted], C[_]]");
    /* 1385 — IC-point transcendental Taylor series. */
    check_true("FreeQ[TimeConstrained[DSolve[{x^2 y''[x] + (x + 1) y'[x] + 3 Log[x] y[x] == 0, "
               "y[1] == 2, y'[1] == 0}, y, x], 8, $Aborted], C[_]]");
    /* 1392 — cubic-coefficient Heun: Kovacic declines the complex poles -> series. */
    check_true("Head[TimeConstrained[DSolve[(x^3 + 1) y''[x] + 4 x y'[x] + y[x] == 0, y, x], "
               "8, $Aborted]] === List");
}

/* M35 — DSolve`PiecewiseForcing (§2.2.15): linear IVPs with piecewise/step/Heaviside
 * forcing, solved by interval continuation into a verified Piecewise closed form. */
static void t_m35_piecewise_forcing(void) {
    /* 1492 — step forcing {1 on [0,Pi), 0 after}: a genuine Integrate-free Piecewise,
     * verified in the SECOND interval (forcing 0) and at the IC (was an inert
     * Integrate[Piecewise[..]*Cos,t] before M35). */
    check_true("FreeQ[TimeConstrained[DSolve[{y''[t] + 4 y[t] == Piecewise[{{1, 0 <= t < Pi}, "
               "{0, Pi <= t}}, 0], y[0] == 1, y'[0] == 0}, y, t], 15, $Aborted], Integrate]");
    check_true("Module[{b}, b = y[t] /. DSolve[{y''[t] + 4 y[t] == Piecewise[{{1, 0 <= t < Pi}, "
               "{0, Pi <= t}}, 0], y[0] == 1, y'[0] == 0}, y, t][[1]]; "
               "!FreeQ[b, Piecewise] && Abs[N[(D[b, {t, 2}] + 4 b) /. t -> 4]] < 10^-6 && "
               "Abs[N[b /. t -> 0] - 1] < 10^-6]");
    /* 1497 — UnitStep forcing Sin[t] - UnitStep[t-2Pi] Sin[t]: the residual carries a
     * UnitStep, which the distributional-keep verify accepts (zero_test spuriously
     * rejected it before). Verified in the second interval (t > 2Pi, forcing 0). */
    check_true("Module[{b}, b = y[t] /. DSolve[{y''[t] + 4 y[t] == Sin[t] - UnitStep[t - 2 Pi] Sin[t], "
               "y[0] == 0, y'[0] == 0}, y, t][[1]]; FreeQ[b, Integrate] && "
               "Abs[N[(D[b, {t, 2}] + 4 b) /. t -> 7]] < 10^-6]");
    /* 1494 — three-piece triangular pulse: residual verified in the MIDDLE interval
     * [1,2) where the forcing is 2-t (exercises the continuity handoff at both ends). */
    check_true("Module[{b}, b = y[t] /. DSolve[{y''[t] + y[t] == Piecewise[{{t, 0 <= t < 1}, "
               "{2 - t, 1 <= t < 2}, {0, 2 <= t}}, 0], y[0] == 1, y'[0] == 0}, y, t][[1]]; "
               "Abs[N[(D[b, {t, 2}] + b - (2 - t)) /. t -> 3/2]] < 10^-6]");
    /* 1499 — damped ramp (UnitStep) forcing on y''+y'+5/4 y: genuine closed form. */
    check_true("FreeQ[TimeConstrained[DSolve[{y''[t] + y'[t] + 5/4 y[t] == "
               "t - UnitStep[t - Pi/2] (t - Pi/2), y[0] == 0, y'[0] == 0}, y, t], 15, $Aborted], Integrate]");
    /* Gate: a SMOOTH nonhomogeneous IVP is NOT piecewise-solved (falls through to
     * UndeterminedCoefficients) — the method fires only on step/piecewise forcing. */
    check_true("FreeQ[DSolve[{y''[t] + w^2 y[t] == Cos[2 t], y[0] == 1, y'[0] == 0}, y, t], Piecewise]");
}

static void t_m36_separable_cubic_log(void) {
    /* 2.2.16-1590 — cubic-in-y separable: Integrate[1/((1+y)(y-1)(y-2)), y] is a sum
     * of THREE distinct logs, which ds_solve cannot invert to y (it churns).  The
     * explicit Separable path now declines a >=3-log relation and the implicit first-
     * integral twin returns G(x,y) == C[1].  Verified here by the implicit-function
     * rule: the total x-derivative of G, with y' replaced by the ODE, is zero. */
    check_true("Module[{r, g}, r = DSolve[y'[x] == -((1 + y[x]) (-1 + y[x]) (y[x] - 2))/(x + 1), y, x]; "
               "MatchQ[r, {{_Equal}}] && (g = r[[1, 1, 1]] - r[[1, 1, 2]]; "
               "PossibleZeroQ[(D[g, x] /. Derivative[1][y][x] -> "
               "-((1 + y[x]) (-1 + y[x]) (y[x] - 2))/(x + 1))])]");
    /* IVP: the constant is fitted (C = G(x0, y0), no inversion), so no C[k] remains. */
    check_true("FreeQ[DSolve[{y'[x] == -((1 + y[x]) (-1 + y[x]) (y[x] - 2))/(x + 1), "
               "y[1] == 0}, y, x], C[_]]");
    /* Anti-overfit: a DIFFERENT cubic-log separable also solves implicitly. */
    check_true("MatchQ[DSolve[y'[x] == (y[x] (y[x] - 1) (y[x] - 3))/x, y, x], {{_Equal}}]");
    /* Regression: a TWO-log (logistic) separable stays on the EXPLICIT path — the
     * >=3-log gate does not over-fire and swallow invertible relations. */
    check_true("MatchQ[DSolve[y'[x] == (y[x]^2 - 1)/x, y, x], {{_Rule}}]");
}

static void t_m37_fractional_power_branch(void) {
    /* §2.2.17 — fractional-power first-order IVPs whose wrong principal-root / spurious
     * fitted-constant branch used to be shipped as a corpus FAIL.  The verifying-root
     * fitter + the numeric branch filters now select a branch that back-substitutes to
     * zero.  Checked at an in-domain point (the closed forms are real only on the
     * sub-interval containing the IC — √y, (y-1)^(1/3) flip sign past a pole). */
    /* 1638: Bernoulli √y IVP -> y=(2E^x-1)^2, NOT the spurious constant y=1. */
    check_true("PossibleZeroQ[((y'[x]-2 y[x]-2 Sqrt[y[x]]) /. "
               "DSolve[{y'[x]-2 y[x]==2 Sqrt[y[x]], y[0]==1}, y, x][[1]]) /. x->1/2]");
    /* 1641: √y, wrong-sign branch rejected. */
    check_true("PossibleZeroQ[((y'[x]-y[x]-x Sqrt[y[x]]) /. "
               "DSolve[{y'[x]-y[x]==x Sqrt[y[x]], y[0]==4}, y, x][[1]]) /. x->1/3]");
    /* 1636: Bernoulli y^(3/2) IVP (solves to a verifiable implicit ArcCoth relation). */
    check_true("MatchQ[DSolve[{y'[x]-y[x] x==y[x]^(3/2) x, y[1]==4}, y, x], {{_Equal}} | {{_Rule}}]");
    /* 1622: cube-root general solution -> the real branch (1+x^3); the spurious complex
     * twin is dropped, so every returned branch back-substitutes to zero for x>0. */
    check_true("Module[{s = DSolve[y'[x]==3 x (-1+y[x])^(1/3), y, x]}, "
               "Head[s]===List && Length[s]>=1 && "
               "AllTrue[s, PossibleZeroQ[((y'[x]-3 x (-1+y[x])^(1/3)) /. #[[1]]) /. x->3/2] &]]");
}

static void t_m37_abel_air(void) {
    /* §2.2.17-1604/1607/1606/1676 — Abel-2nd-kind / rational-in-y forms solved by
     * DSolve`AbelAIR (scaling u=y/s -> separable), returned as an implicit first integral
     * and verified by the implicit-function rule: the total x-derivative of G, with y'
     * replaced by the ODE, is identically zero. */
    /* 1604 (class B): y' + y == 2x E^-x/(1+E^x y). */
    check_true("Module[{r, g, f}, f = (2 x E^(-x))/(1+E^x y[x]) - y[x]; "
               "r = DSolve[y'[x]+y[x]==(2 x E^(-x))/(1+E^x y[x]), y, x]; "
               "MatchQ[r, {{_Equal}}] && (g = r[[1,1,1]] - r[[1,1,2]]; "
               "PossibleZeroQ[D[g, x] /. Derivative[1][y][x] -> f])]");
    /* 1607 (class A): y' - 2y == x E^2x/(1 - y E^-2x). */
    check_true("Module[{r, g, f}, f = (x E^(2 x))/(1-y[x] E^(-2 x)) + 2 y[x]; "
               "r = DSolve[y'[x]-2 y[x]==(x E^(2 x))/(1-y[x] E^(-2 x)), y, x]; "
               "MatchQ[r, {{_Equal}}] && (g = r[[1,1,1]] - r[[1,1,2]]; "
               "PossibleZeroQ[D[g, x] /. Derivative[1][y][x] -> f])]");
    /* 1606 (squared denominator (E^x+y)^2): the degree-2 linear-factor branch. */
    check_true("MatchQ[DSolve[y'[x]-y[x]==((x+1) E^(4 x))/((E^x+y[x])^2), y, x], {{_Equal}}]");
    /* Anti-overfit: a DIFFERENT scaling-reducible Abel-2nd-kind form, (y+x)y' == y + x^2/(y+x). */
    check_true("MatchQ[DSolve`AbelAIR[(y[x]+x) y'[x]==y[x]+x^2/(y[x]+x), y, x], {{_Equal}}]");
    /* Regression: a genuine Riccati (polynomial-in-y denominator free of y) is NOT claimed
     * by AbelAIR -- it declines so Riccati/Chini keep it. */
    check_true("MatchQ[DSolve`AbelAIR[y'[x]==x^2+y[x]^2, y, x], DSolve`AbelAIR[__]]");
}

static void t_m38_forced_decay_wronskian(void) {
    /* §2.2.18-1763: y'' + 4x y' + (4x^2+2) y == 8 E^(-x(x+2)).  Homogeneous set
     * (C1 + C2 x) E^(-x^2); its Wronskian E^(-2x^2) DECAYS, which zero_test's decay
     * false-positive (PossibleZeroQ[E^(-2x^2)] === True) read as identically zero.
     * That corrupted two gates: the Wronskian nonzero-check in
     * dsolve_variation_of_parameters (-> VoP declined) and Kovacic's forcing
     * detection (a decaying forcing read as homogeneous -> dropped the particular,
     * shipping a homogeneous-only WRONG answer masked as UNEVAL by the corpus
     * prelude's leaked-C[k] leniency).  ds_is_structural_zero (Expand[.]===0) gates
     * both now, so Kovacic adds the VoP particular 2 E^(-x^2-2x).  Verified by
     * numeric back-substitution at x=-1, where the forcing is O(1) (8e): a
     * dropped-forcing answer's residual there is large, not ~0. */
    check_true("Module[{s,r}, s=DSolve[y''[x]+4 y'[x] x+(4 x^2+2) y[x]==8 E^(-x (x+2)), y, x]; "
               "MatchQ[s,{{_Rule}}] && (r=(y''[x]+4 y'[x] x+(4 x^2+2) y[x]-8 E^(-x (x+2)))"
               "/.s[[1]]/.{C[1]->1,C[2]->1}; Abs[N[r /. x->-1]] < 1/1000000)]");
    /* The particular is actually present (not homogeneous-only): with both generated
     * constants zeroed the solution is the nonzero particular 2 E^(-x^2-2x) (=2e at x=-1). */
    check_true("Module[{s,p}, s=DSolve[y''[x]+4 y'[x] x+(4 x^2+2) y[x]==8 E^(-x (x+2)), y, x]; "
               "p=(y[x]/.s[[1]]/.{C[1]->0,C[2]->0}); Abs[N[p /. x->-1]] > 1]");
    /* Anti-overfit: a DIFFERENT decaying-Wronskian forced 2nd-order equation --
     * y'' + 2x y' + (x^2+1) y == E^(-x^2/2-x), homogeneous (C1+C2 x)E^(-x^2/2),
     * Wronskian E^(-x^2) (also decaying). */
    check_true("Module[{s,r}, s=DSolve[y''[x]+2 y'[x] x+(x^2+1) y[x]==E^(-x^2/2-x), y, x]; "
               "MatchQ[s,{{_Rule}}] && (r=(y''[x]+2 y'[x] x+(x^2+1) y[x]-E^(-x^2/2-x))"
               "/.s[[1]]/.{C[1]->1,C[2]->1}; Abs[N[r /. x->-1]] < 1/1000000)]");
    /* Regression: the homogeneous equation still solves (fast, correct) -- the fix
     * only changes forcing/Wronskian gating, not the homogeneous Kovacic path. */
    check_true("MatchQ[DSolve[y''[x]+4 y'[x] x+(4 x^2+2) y[x]==0, y, x], {{_Rule}}]");
}

static void t_rischnorman_enum_cap_no_crash(void) {
    check_true("MatchQ[DSolve[x^2 - 1 + (y[x]^2 x^2 + x^3 + x) y'[x] == 0, "
               "y[x], x], _List | _DSolve]");
    check_true("MatchQ[DSolve[x^3 - 1 + (y[x]^2 x^2 + x^3 + x) y'[x] == 0, "
               "y[x], x], _List | _DSolve]");
    check_true("MatchQ[DSolve[2 x^2 - 3 + (y[x]^2 x^2 + x^3 + 2 x) y'[x] == 0, "
               "y[x], x], _List | _DSolve]");
    check_true("MatchQ[DSolve[x^4 - 1 + (y[x]^2 x^4 + x^5 + x) y'[x] == 0, "
               "y[x], x], _List | _DSolve]");
    /* direct multi-kernel Risch-Norman integrands drive a large basis */
    check_true("MatchQ[Integrate[Exp[x] Log[x] Sin[x], x], _]");
    check_true("MatchQ[Integrate[x^3 Exp[x] Log[x]^2 Sin[x], x], _]");
}

/* First-order linear ODEs with trig/elementary coefficients, which route
 * through the integrating-factor solver mu = PowerExpand[Simplify[Exp[Int p]]].
 * Pre-fix these hung (mu = Exp[messy Int Tan] left an un-collapsed
 * Sqrt[Sec^2]/rational-trig integrand that spun in Integrate).  Verified by
 * back-substitution residual. */
static void t_trig_coeff_linear_first_order(void) {
    check_true("PossibleZeroQ[(y'[x] + Tan[x] y[x] - 3 Cos[x]^2) /. "
               "DSolve[y'[x] + Tan[x] y[x] == 3 Cos[x]^2, y[x], x][[1]]]");
    check_true("PossibleZeroQ[(y'[x] - Tan[x] y[x] - Cos[x]) /. "
               "DSolve[y'[x] - Tan[x] y[x] == Cos[x], y[x], x][[1]]]");
    check_true("PossibleZeroQ[(y'[x] + Cot[x] y[x] - Csc[x]) /. "
               "DSolve[y'[x] + Cot[x] y[x] == Csc[x], y[x], x][[1]]]");
    check_true("PossibleZeroQ[(y'[x] + Tanh[x] y[x] - 1) /. "
               "DSolve[y'[x] + Tanh[x] y[x] == 1, y[x], x][[1]]]");
}

/* First-order ODEs linearizable by u = phi(y) (DSolve`Linearizable): reduce to
 * a linear/Bernoulli equation in u, solve, map back y = phi^{-1}(u).  Verified
 * by back-substitution (y -> Function so y'[x] resolves). */
static void t_linearizable_first_order(void) {
    /* Log-family -> linear in Log[y] */
    check_true("PossibleZeroQ[(y'[x] - y[x](Exp[x] + Log[y[x]])) /. "
               "DSolve[y'[x] == y[x](Exp[x] + Log[y[x]]), y, x][[1]]]");
    /* Exp-family -> linear in Exp[y] */
    check_true("PossibleZeroQ[(y'[x] - Exp[x - y[x]](Exp[x] - Exp[y[x]])) /. "
               "DSolve[y'[x] == Exp[x - y[x]](Exp[x] - Exp[y[x]]), y, x][[1]]]");
    /* Sin-family -> linear in Sin[y] */
    check_true("PossibleZeroQ[(y'[x] - Tan[y[x]]/(x+1) - (x+1) Exp[x] Sec[y[x]]) /. "
               "DSolve[y'[x] - Tan[y[x]]/(x+1) == (x+1) Exp[x] Sec[y[x]], y, x][[1]]]");
    /* Sin-family -> Bernoulli in Sin[y] */
    check_true("PossibleZeroQ[(y'[x] Cos[y[x]] - Cos[x] Sin[y[x]]^2 - Sin[y[x]]) /. "
               "DSolve[y'[x] Cos[y[x]] - Cos[x] Sin[y[x]]^2 - Sin[y[x]] == 0, y, x][[1]]]");
    /* Cos-family -> Bernoulli in Cos[y], with a Log[x] coefficient */
    check_true("PossibleZeroQ[(y'[x] - (-2 Cos[y[x]] + x^3 Cos[2 y[x]] Log[x] + x^3 Log[x])"
               "/(2 Sin[y[x]] Log[x] x)) /. "
               "DSolve[y'[x] == (-2 Cos[y[x]] + x^3 Cos[2 y[x]] Log[x] + x^3 Log[x])"
               "/(2 Sin[y[x]] Log[x] x), y, x][[1]]]");
}

/* M12: DSolve`SecondOrderSymmetry — nonlinear second-order Lie point symmetry.
 * Verified NUMERICALLY (the solutions carry logs/radicals PossibleZeroQ cannot
 * decide; lie2 itself relies on a numeric back-substitution guard). */
static void t_second_order_symmetry(void) {
    /* projective symmetry [x,y],[x^2,xy] — Cheb-Terrab et al. Eq.3 */
    check_true("With[{s=DSolve[x^3 y''[x]==(y[x]-x y'[x])^2, y, x]}, "
               "Head[s]===List && Abs[N[(x^3 y''[x]-(y[x]-x y'[x])^2) /. s[[1]] "
               "/. {C[1]->7/5,C[2]->3/4} /. x->13/10, 14]] < 10^-6]");
    /* separable-scaling symmetry [x,0],[0,y] (N5) */
    check_true("With[{s=DSolve[2 x^2 y''[x] y[x]+y[x]^2==x^2 y'[x]^2, y, x]}, "
               "Head[s]===List && Abs[N[(2 x^2 y''[x] y[x]+y[x]^2-x^2 y'[x]^2) /. s[[1]] "
               "/. {C[1]->9/7,C[2]->2/5} /. x->8/5, 14]] < 10^-6]");
    /* pinned method solves; declines a linear ODE (its domain is nonlinear) */
    check_form("Head[DSolve`SecondOrderSymmetry[x^3 y''[x]==(y[x]-x y'[x])^2, y, x]]", "List");
    check_true("Head[DSolve`SecondOrderSymmetry[y''[x]+y[x]==0, y, x]] =!= List");
}

/* M14: DSolve`ChangeOfVariable — transcendental-coefficient 2nd-order linear ODE
 * rationalized by t = phi(x) (Cos/Sin/Tan), e.g. Legendre via t = Cos[x]. */
static void t_change_of_variable(void) {
    /* y'' + Cot[x] y' + k(k+1) y == 0  --(t=Cos[x])-->  Legendre; check C-coeffs vanish */
    check_true("With[{s=DSolve[y''[x]+Cot[x] y'[x]+6 y[x]==0, y, x]}, "
               "Head[s]===List && Module[{r=Simplify[(y''[x]+Cot[x] y'[x]+6 y[x]) /. s[[1]]]}, "
               "Abs[N[Coefficient[r,C[1]] /. x->13/10, 12]] < 10^-5 && "
               "Abs[N[Coefficient[r,C[2]] /. x->13/10, 12]] < 10^-5]]");
    check_form("Head[DSolve`ChangeOfVariable[y''[x]+Cot[x] y'[x]+6 y[x]==0, y, x]]", "List");
    /* declines an already-rational ODE (owned by the direct methods) */
    check_true("Head[DSolve`ChangeOfVariable[y''[x]+y[x]==0, y, x]] =!= List");
}

/* Higher-order linear ODEs whose input form hides their class, fixed by
 * dsolve_linear_normalize: a rational RHS solved for the top derivative clears to
 * Euler form; a common coefficient factor divides out to constant-coefficient. */
static void t_linear_coeff_normalization(void) {
    /* #95: y''' == (24x+24y)/x^3  -> Euler x^3 y''' - 24 y == 24 x */
    check_true("PossibleZeroQ[(y'''[x] - (24 x + 24 y[x])/x^3) /. "
               "DSolve[y'''[x] == (24 x + 24 y[x])/x^3, y, x][[1]]]");
    /* #96: x(y'''+2y''-y'-2y) == 1 -> const-coeff y'''+2y''-y'-2y == 1/x */
    check_true("PossibleZeroQ[(x y'''[x]+2 y''[x] x-y'[x] x-2 y[x] x-1) /. "
               "DSolve[x y'''[x]+2 y''[x] x-y'[x] x-2 y[x] x==1, y, x][[1]]]");
    /* common factor -> const-coeff (2nd order) */
    check_true("PossibleZeroQ[(x y''[x]+3 x y'[x]+2 x y[x]-1) /. "
               "DSolve[x y''[x]+3 x y'[x]+2 x y[x]==1, y, x][[1]]]");
}

/* Variable-coefficient first-order linear systems: the scalar-factor class
 * A(t) = f(t) B (extract_Ab now admits a t-dependent leading coefficient) and
 * the 2x2 commutative class A = a(t) I + b(t) K0 (DSolve`LinearSystemCommutative:
 * rotation / hyperbolic / nilpotent).  Verified by back-substitution residual. */
static void t_system_varcoeff(void) {
    /* scalar-factor A = (1/t) B */
    check_true("PossibleZeroQ[(t x'[t] + y[t]) /. "
               "DSolve[{t x'[t]+y[t]==0, t y'[t]+x[t]==0}, {x[t],y[t]}, t][[1]]]");
    /* scalar-factor with forcing */
    check_true("PossibleZeroQ[(t x'[t]+2 x[t]-2 y[t]-t) /. "
               "DSolve[{t x'[t]+2 x[t]-2 y[t]==t, t y'[t]+x[t]+5 y[t]==t^2}, "
               "{x[t],y[t]}, t][[1]]]");
    /* commutative hyperbolic K0^2 = I */
    check_true("PossibleZeroQ[(x'[t]+x[t]-t y[t]) /. "
               "DSolve[{x'[t]==-x[t]+t y[t], y'[t]==t x[t]-y[t]}, {x[t],y[t]}, t][[1]]]");
    /* commutative rotation K0^2 = -I */
    check_true("PossibleZeroQ[(x'[t]-x[t] Cos[t]+Sin[t] y[t]) /. "
               "DSolve[{x'[t]==x[t] Cos[t]-Sin[t] y[t], y'[t]==x[t] Sin[t]+y[t] Cos[t]}, "
               "{x[t],y[t]}, t][[1]]]");
    /* commutative rotation with a = 1/t */
    check_true("PossibleZeroQ[(x'[t]-x[t]/t-y[t]) /. "
               "DSolve[{x'[t]==x[t]/t+y[t], y'[t]==-x[t]+y[t]/t}, {x[t],y[t]}, t][[1]]]");
}

/* 2-D autonomous (nonlinear) systems via the phase-plane reduction
 * (DSolve`AutonomousSystem): solve the orbit dy/dx==g/f, reconstruct x(t).
 * Verified by back-substitution; radical-orbit systems must decline (not crash). */
static void t_system_autonomous(void) {
    /* orbit y = C x -> exponential */
    check_true("PossibleZeroQ[(x'[t] - y[t]) /. "
               "DSolve[{x'[t]==y[t], y'[t]==y[t]^2/x[t]}, {x[t],y[t]}, t][[1]]]");
    check_true("PossibleZeroQ[(y'[t] - y[t]^2/x[t]) /. "
               "DSolve[{x'[t]==y[t], y'[t]==y[t]^2/x[t]}, {x[t],y[t]}, t][[1]]]");
    /* orbit x y = C */
    check_true("PossibleZeroQ[(x'[t] + 1/y[t]) /. "
               "DSolve[{x'[t]==-1/y[t], y'[t]==1/x[t]}, {x[t],y[t]}, t][[1]]]");
    /* orbit y = C x, Sqrt reconstruction */
    check_true("PossibleZeroQ[(x'[t] - 1/y[t]) /. "
               "DSolve[{x'[t]==1/y[t], y'[t]==1/x[t]}, {x[t],y[t]}, t][[1]]]");
    /* radical orbit -> declines cleanly, no crash (returns unevaluated) */
    check_true("MatchQ[DSolve[{x'[t]==y[t]/(x[t]-y[t]), y'[t]==x[t]/(x[t]-y[t])}, "
               "{x[t],y[t]}, t], _DSolve | _List]");
}

/* Higher-order coupled linear ODE systems via state augmentation
 * (DSolve`SystemReduce).  Verified by back-substitution into the ORIGINAL
 * higher-order equations.  Covers: pure 2nd-order coupled, mixed order (one
 * function 2nd, one 1st), a complex/mixed spectrum, and a forced system
 * (variation of parameters). */
static void t_system_higher_order(void) {
    /* pure 2nd order, x''=4y, y''=4x (spectrum +/-2, +/-2i) */
    check_true("PossibleZeroQ[(x''[t] - 4 y[t]) /. "
               "DSolve[{x''[t]==4 y[t], y''[t]==4 x[t]}, {x[t],y[t]}, t][[1]]]");
    check_true("PossibleZeroQ[(y''[t] - 4 x[t]) /. "
               "DSolve[{x''[t]==4 y[t], y''[t]==4 x[t]}, {x[t],y[t]}, t][[1]]]");
    /* mixed order: x is 2nd order, y is 1st order */
    check_true("PossibleZeroQ[(x''[t]+x'[t]+y'[t]-2 y[t]) /. "
               "DSolve[{x''[t]+x'[t]+y'[t]-2 y[t]==0, x'[t]+x[t]-y'[t]==0}, "
               "{x[t],y[t]}, t][[1]]]");
    check_true("PossibleZeroQ[(x'[t]+x[t]-y'[t]) /. "
               "DSolve[{x''[t]+x'[t]+y'[t]-2 y[t]==0, x'[t]+x[t]-y'[t]==0}, "
               "{x[t],y[t]}, t][[1]]]");
    /* forced 2nd-order system (variation of parameters) */
    check_true("PossibleZeroQ[(x''[t]-4 y[t]-Exp[t]) /. "
               "DSolve[{x''[t]==4 y[t]+Exp[t], y''[t]==4 x[t]-Exp[t]}, "
               "{x[t],y[t]}, t][[1]]]");
    /* the number of arbitrary constants equals the total order (2+2 = 4) */
    check_true("Length[Union[Cases[DSolve[{x''[t]==2 x[t]-3 y[t], "
               "y''[t]==x[t]-2 y[t]}, {x[t],y[t]}, t], C[_], Infinity]]] == 4");
}

/* ---- M18: DSolve`ReducibleIntegratingFactor (2nd-order integrating factor) ---- */

/* mu(x,y) Case A (Section 2.1): Phi degree-2 poly in y', closed-form mu.
 * y y' + y'' == 1  has mu = 1 -> first integral y' - x + y^2/2 == C[1] (Riccati)
 * -> Airy.  Verified numerically on the original 2nd-order residual. */
static void t_m18_mu_xy_caseA(void) {
    check_form("Head[DSolve`ReducibleIntegratingFactor[y[x] y'[x] + y''[x] == 1, y[x], x]]",
               "List");
    check_true("Abs[N[(y[x] y'[x] + y''[x] - 1) /. "
               "DSolve`ReducibleIntegratingFactor[y[x] y'[x] + y''[x] == 1, y, x][[1]] "
               "/. {C[1] -> 13/10, C[2] -> 7/10, x -> 6/5}]] < 1/1000000");
    /* an arbitrary-coefficient generalization y'' - k y y' == c also closes */
    check_true("Abs[N[(y''[x] - y[x] y'[x] - 6) /. "
               "DSolve`ReducibleIntegratingFactor[y''[x] - y[x] y'[x] == 6, y, x][[1]] "
               "/. {C[1] -> 11/10, C[2] -> 3/5, x -> 7/5}]] < 1/1000000");
}

/* mu(x,y) Case B (Section 2.1, 2.18-2.21): the linear-nu-ODE subcase, with
 * transcendental (Coth) coefficients.  -(y'^2/y^2)+y''/y+2Coth[2x]y'/y == 2.
 * The correct answer needs TrigToExp[Coth] (fixed) and a numeric verify since the
 * residual carries ArcTanh/Log terms zero_test cannot discharge. */
static void t_m18_mu_xy_caseB(void) {
    check_true("With[{sol = DSolve`ReducibleIntegratingFactor["
               "-(y'[x]^2/y[x]^2)+(y''[x]/y[x])+(2 Coth[2 x] y'[x]/y[x]) == 2, y, x]}, "
               "Head[sol] === List && Length[sol] >= 1 && "
               "Abs[N[(-(y'[x]^2/y[x]^2)+(y''[x]/y[x])+(2 Coth[2 x] y'[x]/y[x]) - 2) /. "
               "sol[[1]] /. {C[1] -> 13/10, C[2] -> 7/10, x -> 3/5}, 20]] < 10^-6]");
}

/* The method is auto-dispatched (not only pinned): DSolve solves case 14. */
static void t_m18_auto_dispatch(void) {
    check_true("With[{sol = DSolve["
               "-(y'[x]^2/y[x]^2)+(y''[x]/y[x])+(2 Coth[2 x] y'[x]/y[x]) == 2, y, x]}, "
               "Head[sol] === List && Length[sol] >= 1]");
}

/* Linearity gate: a LINEAR 2nd-order ODE is not the method's domain -> the pinned
 * method declines (returns unevaluated), so the linear specialists own it and the
 * Case-B nu-ODE recursion terminates. */
static void t_m18_declines_linear(void) {
    check_form("Head[DSolve`ReducibleIntegratingFactor[y''[x] + y[x] == 0, y[x], x]]",
               "DSolve`ReducibleIntegratingFactor");
    check_form("Head[DSolve`ReducibleIntegratingFactor[y''[x] + (2/x) y'[x] + y[x] == 0, y[x], x]]",
               "DSolve`ReducibleIntegratingFactor");
}

/* Regression for the TrigToExp[Coth] sign bug this wave surfaced (denominator was
 * E^-x - E^x, i.e. -Coth).  Guard it directly so it cannot silently return. */
static void t_m18_trigtoexp_coth(void) {
    check_true("PossibleZeroQ[ExpToTrig[TrigToExp[Coth[z]]] - Coth[z]]");
    check_true("Abs[N[TrigToExp[Coth[2 z]] /. z -> 7/10] - N[Coth[7/5]]] < 1/1000000");
}

int main(void) {
    symtab_init();
    core_init();
    test_load_init_m();   /* match production: deriv.m rules + CRC integral tables */

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
    TEST(t_bvp_overdetermined);
    TEST(t_bvp_underdetermined);
    TEST(t_bvp_system_overdetermined);
    TEST(t_bvp_undecided_keeps_general);
    TEST(t_method_constcoeff);
    TEST(t_euler_complex);
    TEST(t_euler_real);
    TEST(t_euler_repeated);
    TEST(t_euler_inhomogeneous);
    TEST(t_euler_inhomogeneous_complex);
    TEST(t_euler_regression_corpus);
    TEST(t_m34_corpus_cases);
    TEST(t_m35_piecewise_forcing);
    TEST(t_m36_separable_cubic_log);
    TEST(t_m37_fractional_power_branch);
    TEST(t_m37_abel_air);
    TEST(t_m38_forced_decay_wronskian);
    TEST(t_rischnorman_enum_cap_no_crash);
    TEST(t_trig_coeff_linear_first_order);
    TEST(t_linearizable_first_order);
    TEST(t_second_order_symmetry);
    TEST(t_change_of_variable);
    TEST(t_linear_coeff_normalization);
    TEST(t_system_varcoeff);
    TEST(t_system_autonomous);
    TEST(t_system_higher_order);
    TEST(t_method_euler);
    TEST(t_method_exactode);
    TEST(t_exactode_more);
    TEST(t_exactode_declines);
    TEST(t_exactode_auto);
    TEST(t_airy);
    TEST(t_bessel);
    TEST(t_bessel_modified);
    TEST(t_method_specialform);
    TEST(t_bessel_reducible);
    TEST(t_kovacic_highdegree_no_hang);
    TEST(t_hypergeometric_kummer);
    TEST(t_hypergeometric_gauss);
    TEST(t_method_hypergeometric_kummer);
    TEST(t_hypergeometric_symbolic_a);
    TEST(t_hypergeometric_gauss_symbolic_ab);
    TEST(t_hypergeometric_integer_declines);
    /* M17: affine -> Gauss 2F1 + normal-form pre-pass */
    TEST(t_m17_affine_gauss);
    TEST(t_m17_gegenbauer_symbolic);
    TEST(t_m17_associated_legendre);
    TEST(t_m17_normalform_bessel);
    TEST(t_m17_affine_declines_confluent);
    TEST(t_m19_whittaker_confluent);
    TEST(t_m19_declines_integer_2mu);
    TEST(t_m20_polyshift_402);
    TEST(t_m20_polyshift_371);
    TEST(t_m20_declines_nonlinear_base);
    TEST(t_m23_exact_radical);
    TEST(t_m24_trig_power_forcing);
    TEST(t_m24_complex_cuberoot_ivp);
    TEST(t_m25_exact_erf);
    TEST(t_m25_kovacic_fundamental_set);
    TEST(t_m25_transcendental_frobenius);
    TEST(t_m26_distributions);
    TEST(t_m26_impulse_forcing);
    TEST(t_m26_general_forcing);
    TEST(t_m27_system_verify);
    TEST(t_m27_separable_inverse_constant);
    TEST(t_m27_ivp_family_intact);
    TEST(t_m28_bernoulli_hang_trig_substitution);
    TEST(t_m29_sec_floor_verifies);
    TEST(t_m30_linsys_irrational_forcing);
    TEST(t_m30_kovacic_inhomogeneous);
    TEST(t_m31_triangular_exp_forcing);
    TEST(t_m31_linsys_large_eigenvalue);
    /* M32: §2.2.12 IVP-fitter / Separable-implicit / Integrate-Erf-variable */
    TEST(t_m32_ivp_unsatisfiable_branch_dropped);
    TEST(t_m32_bvp_underdetermined_keeps_constant);
    TEST(t_m32_root_form_ivp);
    TEST(t_m32_separable_implicit_and_widened);
    TEST(t_m32_integrate_erf_variable);
    /* M33: §2.2.13 exact-method robustness (transcendental / denominator-clearing /
     * mu(y) trig + implicit-verify Together) and exact-vs-homogeneous IVP fall-through */
    TEST(t_m33_exact_transcendental);
    TEST(t_m33_exact_rational_clear);
    TEST(t_m33_exact_mu_trig);
    TEST(t_m33_exact_homogeneous_ivp_fallthrough);
    /* M5: NormalForm + Kovacic + Frobenius/PowerSeries */
    TEST(t_normalform_bessel);
    TEST(t_normalform_const);
    TEST(t_normalform_declines);
    TEST(t_kovacic_case1_exp);
    TEST(t_kovacic_apparent_singularity);
    TEST(t_kovacic_case1_pole);
    TEST(t_kovacic_case2);
    TEST(t_kovacic_auto_closed_form);
    TEST(t_kovacic_declines);
    TEST(t_kovacic_legendre1);
    TEST(t_kovacic_legendre2);
    TEST(t_kovacic_chebyshev2);
    TEST(t_kovacic_complex_poles);
    TEST(t_kovacic_legendre_auto_closed_form);
    TEST(t_kovacic_case2_complex_pole_no_hang);
    TEST(t_method_operfactor);
    TEST(t_operfactor_more);
    TEST(t_operfactor_ivp);
    TEST(t_operfactor_declines);
    TEST(t_dfactor);
    TEST(t_operfactor_auto);
    TEST(t_powerseries_ordinary);
    TEST(t_powerseries_auto);
    TEST(t_frobenius_regsing_distinct);
    TEST(t_frobenius_regsing_log);
    TEST(t_frobenius_declines_irregular);
    TEST(t_sys_decoupled);
    TEST(t_sys_real_eigenvalues);
    TEST(t_sys_complex_ivp);
    TEST(t_sys_constant_forcing);
    TEST(t_sys_defective_singular);
    TEST(t_sys_defective_nontriangular);
    TEST(t_sys_triangular_varcoeff);
    TEST(t_sys_singular_forcing);
    TEST(t_sys_triangular_ivp);
    TEST(t_sys_varcoeff_coupled);
    TEST(t_sys_varcoeff_complex);
    TEST(t_sys_varcoeff_forced);
    TEST(t_sys_varcoeff_pinned_decline);
    TEST(t_eig_dirichlet);
    TEST(t_eig_neumann);
    TEST(t_eig_mixed);
    TEST(t_eig_no_misfire);
    TEST(t_reduce_order);
    TEST(t_fos_quadratic);
    TEST(t_fos_shifted);
    TEST(t_fos_distinct_coeff);
    TEST(t_fos_method);
    TEST(t_fos_stress);
    TEST(t_method_riccati);
    TEST(t_riccati_more);
    TEST(t_ivp_riccati);
    TEST(t_method_chini);
    TEST(t_chini_more);
    TEST(t_chini_linremoval);
    TEST(t_method_abel);
    TEST(t_method_lie);
    TEST(t_lie_declines);
    TEST(t_lie_undefined_function_declines);
    TEST(t_lie_linear_coefficients);
    TEST(t_lie_bivariate);
    TEST(t_lie_abaco1_product);
    TEST(t_lie_abaco2_similar);
    TEST(t_lie_function_sum);
    TEST(t_lie_abaco2_unique_unknown);
    TEST(t_lie_abaco2_order_zero);
    TEST(t_lie_chi);
    TEST(t_exact_xayb);
    TEST(t_auto_exp);
    TEST(t_auto_power);
    TEST(t_auto_reciprocal);
    TEST(t_auto_method);
    TEST(t_auto_declines_elliptic);
    TEST(t_auto_stress);
    TEST(t_not_holdall);
    TEST(t_pde_transport);
    TEST(t_pde_forcing);
    TEST(t_pde_zeroth_order);
    TEST(t_pde2_wave);
    TEST(t_pde2_laplace);
    TEST(t_pde2_mixed);
    TEST(t_pde2_repeated);
    TEST(t_pde2_distinct_asymmetric);
    TEST(t_pde2_pinned_and_declines);
    TEST(t_pdesep);
    TEST(t_pdeclassify);
    TEST(t_wave_dalembert);
    TEST(t_heat_kernel);
    TEST(t_sys_pde_pinned_methods);
    TEST(t_pde_quasilinear);
    TEST(t_pde_clairaut);
    TEST(t_pde2_lower_order);
    TEST(t_pde_charpit);
    TEST(t_declines_unsupported);
    /* M9: backfill for thin methods */
    TEST(t_method_homogeneous);
    TEST(t_homogeneous_more);
    TEST(t_ivp_homogeneous);
    TEST(t_method_reduce_order);
    TEST(t_reduce_order_more);
    TEST(t_method_clairaut);
    TEST(t_clairaut_more);
    TEST(t_method_lagrange);
    TEST(t_lagrange_more);
    TEST(t_lagrange_declines);
    TEST(t_lagrange_singular);
    TEST(t_powerseries_more);
    TEST(t_homogeneous_algebraic);
    TEST(t_reduce_order_riccati);
    TEST(t_homogeneous_implicit);
    /* M9: SymPy deterministic parity gaps */
    TEST(t_method_factorable);
    TEST(t_factorable_more);
    TEST(t_factorable_declines);
    TEST(t_method_nth_algebraic);
    TEST(t_nth_algebraic_more);
    TEST(t_nth_algebraic_declines);
    TEST(t_method_lincoeff);
    TEST(t_lincoeff_more);
    TEST(t_lincoeff_declines);
    TEST(t_method_almostlinear);
    TEST(t_almostlinear_declines);
    TEST(t_method_sepreduced);
    TEST(t_sepreduced_declines);
    TEST(t_method_liouville);
    TEST(t_liouville_declines);
    TEST(t_method_undetcoeff);
    TEST(t_undetcoeff_declines);
    TEST(t_method_first_order_series);
    TEST(t_first_order_series_declines);
    /* M18: reducible-mu integrating factor (2nd-order) */
    TEST(t_m18_mu_xy_caseA);
    TEST(t_m18_mu_xy_caseB);
    TEST(t_m18_auto_dispatch);
    TEST(t_m18_declines_linear);
    TEST(t_m18_trigtoexp_coth);

    printf("\nAll DSolve tests passed.\n");
    return 0;
}
