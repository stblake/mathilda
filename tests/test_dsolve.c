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
    TEST(t_sys_pde_pinned_methods);
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

    printf("\nAll DSolve tests passed.\n");
    return 0;
}
