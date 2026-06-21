/* Tests for NIntegrate — numerical integration.
 *
 * Phase 1 covers one-dimensional integrals over a finite real interval at
 * machine precision (globally-adaptive Gauss-Kronrod): smooth integrands,
 * polynomial exactness, option handling (PrecisionGoal / AccuracyGoal /
 * MaxRecursion / MaxPoints / Method), the HoldAll/Protected attributes, and
 * memory hygiene.  Later phases extend this file with endpoint singularities,
 * infinite ranges, complex contours, arbitrary precision, multidimensional
 * iteration, oscillatory and Monte-Carlo methods, Exclusions and principal
 * values.
 *
 * Numerical results are compared *inside* the language
 * (N[Abs[result - expected]] < tol) rather than by parsing the printed form.
 */

#include "core.h"
#include "eval.h"
#include "expr.h"
#include "parse.h"
#include "print.h"
#include "symtab.h"
#include "test_utils.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char* eval_str(const char* input) {
    Expr* p = parse_expression(input);
    ASSERT(p != NULL);
    Expr* e = evaluate(p);
    expr_free(p);
    char* s = expr_to_string(e);
    expr_free(e);
    return s;
}

/* |input - expected| < tol, evaluated at full internal precision; the L1
 * magnitude keeps a spurious tiny imaginary residual from leaving it complex. */
static bool close_to(const char* input, const char* expected, double tol) {
    char buf[4096];
    snprintf(buf, sizeof buf,
             "N[Abs[Re[(%s) - (%s)]] + Abs[Im[(%s) - (%s)]]] < %.17g",
             input, expected, input, expected, tol);
    char* s = eval_str(buf);
    bool ok = (strcmp(s, "True") == 0);
    if (!ok) fprintf(stderr, "  close_to FAIL: |%s - (%s)| < %g  =>  %s\n",
                     input, expected, tol, s);
    free(s);
    return ok;
}

#define ASSERT_CLOSE(input, expected, tol)                                   \
    ASSERT_MSG(close_to((input), (expected), (tol)),                         \
               "%s ~= %s (tol %g)", (input), (expected), (tol))

#define ASSERT_EQ_STR(input, expected)                                       \
    do {                                                                     \
        char* _s = eval_str(input);                                          \
        ASSERT_MSG(strcmp(_s, (expected)) == 0, "%s => %s (want %s)",        \
                   (input), _s, (expected));                                 \
        free(_s);                                                            \
    } while (0)

/* ---------------------------------------------------------------------- */

static void test_smooth_finite(void) {
    ASSERT_CLOSE("NIntegrate[Sin[Sin[x]],{x,0,2}]", "1.247056", 1e-5);
    ASSERT_CLOSE("NIntegrate[Exp[Exp[-x]],{x,0,5}]", "6.311156", 1e-5);
    ASSERT_CLOSE("NIntegrate[x^x,{x,1,5}]", "1241.0332840", 1e-4);
    /* closed forms */
    ASSERT_CLOSE("NIntegrate[Sin[x],{x,0,Pi}]", "2", 1e-9);
    ASSERT_CLOSE("NIntegrate[Exp[x],{x,0,1}]", "E - 1", 1e-10);
    ASSERT_CLOSE("NIntegrate[1/(1+x^2),{x,0,1}]", "Pi/4", 1e-10);
    ASSERT_CLOSE("NIntegrate[Cos[x],{x,0,Pi/2}]", "1", 1e-10);
}

static void test_polynomial_exact(void) {
    /* G7-K15 integrates degree <= 21 exactly on a single panel. */
    ASSERT_CLOSE("NIntegrate[1+x+x^3,{x,0,10}]", "2560", 1e-6);
    ASSERT_CLOSE("NIntegrate[x^5,{x,0,2}]", "32/3", 1e-9);
    ASSERT_CLOSE("NIntegrate[3 x^2 - 2 x + 1,{x,-1,2}]", "9", 1e-9);
}

static void test_negative_and_shifted(void) {
    ASSERT_CLOSE("NIntegrate[x^2,{x,-2,2}]", "16/3", 1e-9);
    ASSERT_CLOSE("NIntegrate[Sin[x]^2,{x,0,2 Pi}]", "Pi", 1e-9);
    /* reversed limits negate */
    ASSERT_CLOSE("NIntegrate[x,{x,3,1}]", "-4", 1e-9);
}

static void test_endpoint_singularities(void) {
    /* Algebraic singularities at the lower endpoint. */
    ASSERT_CLOSE("NIntegrate[1/Sqrt[x],{x,0,1}]", "2", 1e-6);
    ASSERT_CLOSE("NIntegrate[1/x^(1/5),{x,0,1}]", "5/4", 1e-6);
    ASSERT_CLOSE("NIntegrate[1/x^(4/5),{x,0,1}]", "5", 1e-5);
    /* Algebraic singularity at the upper endpoint. */
    ASSERT_CLOSE("NIntegrate[1/Sqrt[1-x],{x,0,1}]", "2", 1e-6);
    /* Logarithmic singularity. */
    ASSERT_CLOSE("NIntegrate[Log[x],{x,0,1}]", "-1", 1e-7);
    ASSERT_CLOSE("NIntegrate[Sqrt[x] Log[x],{x,0,1}]", "-4/9", 1e-7);
    /* Singularities at BOTH endpoints (Beta-type). */
    ASSERT_CLOSE("NIntegrate[1/Sqrt[x(1-x)],{x,0,1}]", "Pi", 1e-6);
    /* Forced DoubleExponential on a smooth integrand still works. */
    ASSERT_CLOSE("NIntegrate[Exp[Cos[x]],{x,0,10},Method->\"DoubleExponential\"]",
                 "12.15611", 1e-4);
}

static void test_infinite_ranges(void) {
    /* Semi-infinite (a, ∞). */
    ASSERT_CLOSE("NIntegrate[1/(1+x^2),{x,0,Infinity}]", "Pi/2", 1e-7);
    ASSERT_CLOSE("NIntegrate[Exp[-x],{x,0,Infinity}]", "1", 1e-8);
    ASSERT_CLOSE("NIntegrate[x Exp[-x],{x,0,Infinity}]", "1", 1e-7);
    ASSERT_CLOSE("NIntegrate[Exp[-x^2],{x,0,Infinity}]", "Sqrt[Pi]/2", 1e-7);
    /* Lower limit -∞. */
    ASSERT_CLOSE("NIntegrate[1/(1+x^2),{x,-Infinity,0}]", "Pi/2", 1e-7);
    ASSERT_CLOSE("NIntegrate[Exp[x],{x,-Infinity,0}]", "1", 1e-8);
    /* Doubly-infinite (-∞, ∞). */
    ASSERT_CLOSE("NIntegrate[Exp[-x^2],{x,-Infinity,Infinity}]", "Sqrt[Pi]", 1e-7);
    ASSERT_CLOSE("NIntegrate[1/(1+x^2),{x,-Infinity,Infinity}]", "Pi", 1e-6);
    /* Reversed orientation negates. */
    ASSERT_CLOSE("NIntegrate[Exp[-x],{x,Infinity,0}]", "-1", 1e-8);
}

static void test_complex_contours(void) {
    /* Straight-line segment: ∫ of an analytic function = antiderivative diff. */
    ASSERT_CLOSE("NIntegrate[2 z,{z,0,1+I}]", "2 I", 1e-9);          /* z^2 | = (1+I)^2 */
    ASSERT_CLOSE("NIntegrate[1,{z,I,3-I}]", "3 - 2 I", 1e-9);        /* endpoint diff */
    ASSERT_CLOSE("NIntegrate[Sqrt[x],{x,I,3-I}]", "3.79214 - 2.21131 I", 1e-3);
    /* Real-variable parametrisation of a circular contour: 2 Pi I (residue 1). */
    ASSERT_CLOSE("NIntegrate[(I Exp[I w])/(Exp[I w]+1/2),{w,0,2 Pi}]", "2 Pi I", 1e-6);
    /* Closed piecewise-linear polygon around -1/2: encloses one simple pole. */
    ASSERT_CLOSE("NIntegrate[1/(z+1/2),{z,1,E^(I Pi/3),E^(2 I Pi/3),-1,"
                 "E^(-2 I Pi/3),E^(-I Pi/3),1}]", "2 Pi I", 1e-5);
}

static void test_multidim(void) {
    ASSERT_CLOSE("NIntegrate[x y,{x,0,1},{y,0,1}]", "1/4", 1e-8);
    ASSERT_CLOSE("NIntegrate[x^2+y^2,{x,0,1},{y,0,1}]", "2/3", 1e-8);
    ASSERT_CLOSE("NIntegrate[x y z,{x,0,1},{y,0,1},{z,0,1}]", "1/8", 1e-7);
    /* Mixed finite / semi-infinite: (∫_0^∞ e^-x)(∫_0^1 e^y) = E - 1. */
    ASSERT_CLOSE("NIntegrate[Exp[-x+y],{x,0,Infinity},{y,0,1}]", "E - 1", 1e-6);
    /* Variable-dependent inner bound: area of the unit disk. */
    ASSERT_CLOSE("NIntegrate[1,{x,-1,1},{y,-Sqrt[1-x^2],Sqrt[1-x^2]}]", "Pi", 1e-5);
    /* Adaptive Genz-Malik cubature over a constant box: a corner singularity
     * (integrable; never sampled) that hangs naive iterated quadrature. */
    ASSERT_CLOSE("NIntegrate[1/Sqrt[x+y^2+z^3],{x,0,1},{y,0,1},{z,0,1}]",
                 "1.088537", 1e-4);
    /* Cancelling integrand whose signed value is 0 — the magnitude-scaled
     * tolerance must converge rather than chase an unreachable relative goal. */
    ASSERT_CLOSE("NIntegrate[Sin[x+y],{x,0,Pi},{y,0,Pi}]", "0", 1e-6);
    ASSERT_CLOSE("NIntegrate[Cos[x+y],{x,0,Pi},{y,0,Pi}]", "-4", 1e-6);
}

static void test_vector_integrand(void) {
    ASSERT_EQ_STR("Head[NIntegrate[{x, Sin[x]},{x,0,1}]]", "List");
    /* {∫x, ∫1/Sqrt[x], ∫Sin[x]} over (0,5) = {25/2, 2 Sqrt[5], 1-Cos[5]}. */
    ASSERT_EQ_STR("N[Total[Abs[NIntegrate[{x, 1/Sqrt[x], Sin[x]},{x,0,5}] "
                  "- {25/2, 2 Sqrt[5], 1 - Cos[5]}]]] < 1/100000", "True");
    /* Matrix-valued integrand threads element-wise. */
    ASSERT_EQ_STR("Dimensions[NIntegrate[{{x, x^2},{x^3, x}},{x,0,1}]]", "{2, 2}");
}

static void test_oscillatory(void) {
    /* Semi-infinite oscillatory (integrate between zeros + Wynn epsilon):
     * ∫_1^∞ Sin[x]/x dx = Pi/2 - SinIntegral[1] = 0.62471325... */
    ASSERT_CLOSE("NIntegrate[Sin[x]/x,{x,1,Infinity}]", "0.6247132564277136", 1e-5);
    ASSERT_CLOSE("NIntegrate[Sin[x]/x,{x,1,Infinity},Method->\"LevinRule\"]",
                 "0.6247132564277136", 1e-5);
    /* Finite oscillatory over many periods (half-period panel summation). */
    ASSERT_EQ_STR("N[Total[Abs[NIntegrate[{Cos[x],Sin[x]},{x,1,100}] "
                  "- {Sin[100]-Sin[1], Cos[1]-Cos[100]}]]] < 1/100000", "True");
    ASSERT_CLOSE("NIntegrate[Cos[x],{x,0,100}]", "Sin[100]", 1e-5);
    /* Doubly-infinite, non-decaying oscillation (Fresnel): sinh-sinh diverges,
     * so the line is split at 0 and each half integrated between the zeros.
     * ∫_{-∞}^{∞} Exp[I x^2] dx = Sqrt[Pi/2] (1 + I); the real and imaginary
     * Fresnel integrals each give Sqrt[Pi/2]. */
    ASSERT_CLOSE("NIntegrate[Cos[x^2],{x,-Infinity,Infinity}]", "Sqrt[Pi/2]", 1e-5);
    ASSERT_CLOSE("NIntegrate[Sin[x^2],{x,-Infinity,Infinity}]", "Sqrt[Pi/2]", 1e-5);
    ASSERT_EQ_STR("Abs[NIntegrate[Exp[I x^2],{x,-Infinity,Infinity}] "
                  "- Sqrt[Pi/2](1+I)] < 1/100000", "True");
}

static void test_oscillatory_singularity(void) {
    /* Oscillatory endpoint singularity: the integrand oscillates infinitely
     * fast under a 1/x envelope as x -> 0.  Adaptive Gauss-Kronrod and
     * tanh-sinh both fail; the exponential endpoint map x = e^{-t} carries it to
     * the non-decaying oscillatory half-line integral Cos[t Exp[t]], integrated
     * between the zeros with Wynn extrapolation.  (Mathematica: 0.3233674...) */
    ASSERT_CLOSE("NIntegrate[Cos[Log[x]/x]/x,{x,0,1}]", "0.32336743", 1e-5);
    /* Same, requested explicitly by name. */
    ASSERT_CLOSE("NIntegrate[Cos[Log[x]/x]/x,{x,0,1},"
                 "Method->\"OscillatorySingularity\"]", "0.32336743", 1e-5);
    /* Mirror singularity at the upper endpoint (exercises the b-side map). */
    ASSERT_CLOSE("NIntegrate[Cos[Log[1-x]/(1-x)]/(1-x),{x,0,1}]", "0.32336743", 1e-5);
    /* The named method also handles a plain algebraic endpoint singularity. */
    ASSERT_CLOSE("NIntegrate[1/Sqrt[x],{x,0,1},"
                 "Method->\"OscillatorySingularity\"]", "2", 1e-5);
}

static void test_montecarlo(void) {
    /* Region (Boole) integrands route to Monte-Carlo automatically. */
    ASSERT_CLOSE("NIntegrate[Boole[x^2+y^2<1],{x,-1,1},{y,-1,1}]", "Pi", 0.02);
    ASSERT_CLOSE("NIntegrate[Boole[1/4<=x^2+y^2<=1],{x,-1,1},{y,-1,1}]", "3 Pi/4", 0.02);
    /* Explicit Monte-Carlo / quasi-Monte-Carlo methods. */
    ASSERT_CLOSE("NIntegrate[Exp[Cos[x]],{x,0,10},Method->\"MonteCarlo\"]", "12.15611", 0.05);
    ASSERT_CLOSE("NIntegrate[Exp[Cos[x]],{x,0,10},Method->\"QuasiMonteCarlo\"]", "12.15611", 0.02);
    /* High-dimensional box (>=6 vars) routes to Monte-Carlo automatically:
     * ∫_[0,1]^6 (x1+...+x6) dV = 3. */
    ASSERT_CLOSE("NIntegrate[x1+x2+x3+x4+x5+x6,{x1,0,1},{x2,0,1},{x3,0,1},"
                 "{x4,0,1},{x5,0,1},{x6,0,1}]", "3", 0.05);
}

static void test_exclusions_pv(void) {
    /* Interior singularity handled by splitting (iterator form and Exclusions). */
    ASSERT_CLOSE("NIntegrate[Log[(1-x)^2],{x,0,1,2}]", "-4", 1e-4);
    ASSERT_CLOSE("NIntegrate[Log[(1-x)^2],{x,0,2},Exclusions->{1}]", "-4", 1e-4);
    /* Cauchy principal values about the Exclusions poles. */
    ASSERT_CLOSE("NIntegrate[1/(x-x^2),{x,-1,2},Method->\"PrincipalValue\","
                 "Exclusions->x-x^2==0]", "Log[4]", 1e-4);
    ASSERT_CLOSE("NIntegrate[1/Log[x],{x,0,10},Method->\"PrincipalValue\",Exclusions->1]",
                 "6.16560", 1e-3);
}

static void test_high_precision(void) {
    /* Finite, smooth: ∫_0^12 e^{-t^2} ≈ Sqrt[Pi]/2 (tail erfc(12)/2 ~ 1e-64). */
    ASSERT_CLOSE("NIntegrate[Exp[-t^2],{t,0,12},WorkingPrecision->100]",
                 "Sqrt[Pi]/2", 1e-60);
    /* Endpoint singular at high precision. */
    ASSERT_CLOSE("NIntegrate[1/(1+Sqrt[x]),{x,0,1},WorkingPrecision->60,PrecisionGoal->40]",
                 "2 - 2 Log[2]", 1e-38);
    ASSERT_CLOSE("NIntegrate[1/Sqrt[x],{x,0,1},WorkingPrecision->40]", "2", 1e-28);
    /* Semi-infinite and doubly-infinite at high precision. */
    ASSERT_CLOSE("NIntegrate[1/(1+x^2),{x,0,Infinity},WorkingPrecision->50]", "Pi/2", 1e-40);
    ASSERT_CLOSE("NIntegrate[Exp[-x^2],{x,-Infinity,Infinity},WorkingPrecision->40]",
                 "Sqrt[Pi]", 1e-30);
    /* The result carries the requested precision (not machine). */
    ASSERT_EQ_STR("Precision[NIntegrate[Exp[-x^2],{x,0,Infinity},WorkingPrecision->50]] > 40",
                  "True");
}

static void test_accuracygoal(void) {
    /* Integrand ~1e-5 scale; AccuracyGoal->8 should still nail it. */
    ASSERT_CLOSE("NIntegrate[10^(-5)/(1+10^2 (x-1/2)^2),{x,0,1},AccuracyGoal->8]",
                 "ArcTan[5]/500000", 1e-11);
}

static void test_methods_accepted(void) {
    /* Implemented strategies/rules give the (method-independent) smooth result. */
    ASSERT_CLOSE("NIntegrate[Exp[Cos[x]],{x,0,10},Method->\"GlobalAdaptive\"]",
                 "12.15611", 1e-4);
    ASSERT_CLOSE("NIntegrate[Exp[Cos[x]],{x,0,10},Method->\"GaussKronrodRule\"]",
                 "12.15611", 1e-4);
}

/* ---- Equally-spaced composite rules (Riemann / Trapezoidal / Newton–Cotes) ---- */

static void test_riemann_rule(void) {
    /* Midpoint rule is exact for linear integrands at any panel count. */
    ASSERT_CLOSE("NIntegrate[x,{x,0,1},Method->{\"RiemannRule\",\"Type\"->\"Midpoint\"}]",
                 "1/2", 1e-12);
    ASSERT_CLOSE("NIntegrate[2 x + 3,{x,-1,2},Method->{\"RiemannRule\",\"Type\"->\"Midpoint\"}]",
                 "12", 1e-12);
    /* Midpoint is second order: a modest PrecisionGoal converges on smooth f. */
    ASSERT_CLOSE("NIntegrate[Exp[x],{x,0,1},Method->{\"RiemannRule\",\"Type\"->\"Midpoint\"},"
                 "PrecisionGoal->6]", "E - 1", 1e-4);
    ASSERT_CLOSE("NIntegrate[Sin[x],{x,0,Pi},Method->{\"RiemannRule\",\"Type\"->\"Midpoint\"},"
                 "PrecisionGoal->5]", "2", 1e-3);
    /* Left / Right are first order: they bracket a monotone integrand. */
    ASSERT_CLOSE("NIntegrate[Exp[x],{x,0,1},Method->{\"RiemannRule\",\"Type\"->\"Left\"},"
                 "PrecisionGoal->4]", "E - 1", 5e-4);
    ASSERT_CLOSE("NIntegrate[Exp[x],{x,0,1},Method->{\"RiemannRule\",\"Type\"->\"Right\"},"
                 "PrecisionGoal->4]", "E - 1", 5e-4);
    /* Left underestimates, Right overestimates a strictly increasing integrand. */
    ASSERT_EQ_STR("NIntegrate[Exp[x],{x,0,1},Method->{\"RiemannRule\",\"Type\"->\"Left\"},"
                  "PrecisionGoal->4] < E - 1", "True");
    ASSERT_EQ_STR("NIntegrate[Exp[x],{x,0,1},Method->{\"RiemannRule\",\"Type\"->\"Right\"},"
                  "PrecisionGoal->4] > E - 1", "True");
    /* Default Type is Left when the sub-option is omitted. */
    ASSERT_CLOSE("NIntegrate[Exp[x],{x,0,1},Method->\"RiemannRule\",PrecisionGoal->4]",
                 "E - 1", 5e-4);
    ASSERT_EQ_STR("NIntegrate[Exp[x],{x,0,1},Method->\"RiemannRule\",PrecisionGoal->4] =="
                  "NIntegrate[Exp[x],{x,0,1},Method->{\"RiemannRule\",\"Type\"->\"Left\"},"
                  "PrecisionGoal->4]", "True");
    /* The documented examples (uniform-schedule values; Left > true > Right). */
    ASSERT_CLOSE("NIntegrate[Exp[Cos[x]],{x,0,10},Method->{\"RiemannRule\",\"Type\"->\"Left\"},"
                 "PrecisionGoal->2]", "12.2455", 1e-3);
    ASSERT_CLOSE("NIntegrate[Exp[Cos[x]],{x,0,10},Method->{\"RiemannRule\",\"Type\"->\"Right\"},"
                 "PrecisionGoal->2]", "12.0669", 1e-3);
}

static void test_trapezoidal_rule(void) {
    /* Trapezoidal is exact for linear integrands. */
    ASSERT_CLOSE("NIntegrate[3 x - 1,{x,0,2},Method->\"TrapezoidalRule\"]", "4", 1e-12);
    /* Romberg (default) — smooth integrands to high accuracy. */
    ASSERT_CLOSE("NIntegrate[Exp[Cos[x]],{x,0,10},Method->\"TrapezoidalRule\"]", "12.15611", 1e-4);
    ASSERT_CLOSE("NIntegrate[4/(1+x^2),{x,0,1},Method->\"TrapezoidalRule\"]", "Pi", 1e-9);
    ASSERT_CLOSE("NIntegrate[Sin[x],{x,0,Pi},Method->\"TrapezoidalRule\"]", "2", 1e-9);
    ASSERT_CLOSE("NIntegrate[x^2,{x,-2,2},Method->\"TrapezoidalRule\"]", "16/3", 1e-9);
    /* Romberg requested explicitly is identical to the default. */
    ASSERT_CLOSE("NIntegrate[Exp[Cos[x]],{x,0,10},"
                 "Method->{\"TrapezoidalRule\",\"RombergQuadrature\"->True}]", "12.15611", 1e-4);
    /* Plain piecewise-linear (no extrapolation) still converges, more slowly. */
    ASSERT_CLOSE("NIntegrate[Exp[Cos[x]],{x,0,10},"
                 "Method->{\"TrapezoidalRule\",\"RombergQuadrature\"->False}]", "12.15611", 1e-3);
    ASSERT_CLOSE("NIntegrate[4/(1+x^2),{x,0,1},"
                 "Method->{\"TrapezoidalRule\",\"RombergQuadrature\"->False}]", "Pi", 1e-5);
}

static void test_newtoncotes_rule(void) {
    /* Default (Points->3, Simpson) on the documented example. */
    ASSERT_CLOSE("NIntegrate[Exp[Cos[x]],{x,0,10},Method->\"NewtonCotesRule\"]", "12.15611", 1e-4);
    ASSERT_CLOSE("NIntegrate[4/(1+x^2),{x,0,1},Method->\"NewtonCotesRule\"]", "Pi", 1e-9);
    /* Simpson (Points->3) is exact for cubics. */
    ASSERT_CLOSE("NIntegrate[x^3,{x,0,2},Method->{\"NewtonCotesRule\",\"Points\"->3}]", "4", 1e-9);
    ASSERT_CLOSE("NIntegrate[x^2,{x,0,2},Method->{\"NewtonCotesRule\",\"Points\"->3}]", "8/3", 1e-9);
    /* Points->2 is the (composite) trapezoid rule. */
    ASSERT_CLOSE("NIntegrate[Sin[x],{x,0,Pi},Method->{\"NewtonCotesRule\",\"Points\"->2}]", "2", 1e-9);
    /* Points->4 (Simpson 3/8) is exact for cubics. */
    ASSERT_CLOSE("NIntegrate[x^3,{x,0,2},Method->{\"NewtonCotesRule\",\"Points\"->4}]", "4", 1e-9);
    ASSERT_CLOSE("NIntegrate[Exp[x],{x,0,1},Method->{\"NewtonCotesRule\",\"Points\"->4}]", "E - 1", 1e-9);
    /* Points->5 (Boole) is exact for quintics. */
    ASSERT_CLOSE("NIntegrate[x^5,{x,0,2},Method->{\"NewtonCotesRule\",\"Points\"->5}]", "32/3", 1e-8);
    ASSERT_CLOSE("NIntegrate[1/(1+x^2),{x,0,1},Method->{\"NewtonCotesRule\",\"Points\"->5}]",
                 "Pi/4", 1e-10);
    /* Out-of-range Points is clamped (1 -> 2, 99 -> 5), still a valid result. */
    ASSERT_CLOSE("NIntegrate[Sin[x],{x,0,Pi},Method->{\"NewtonCotesRule\",\"Points\"->99}]", "2", 1e-8);
}

static void test_fixed_rule_high_precision(void) {
    /* Romberg trapezoidal at arbitrary precision (smooth, no endpoint trouble). */
    ASSERT_CLOSE("NIntegrate[4/(1+x^2),{x,0,1},Method->\"TrapezoidalRule\",WorkingPrecision->30]",
                 "Pi", 1e-25);
    ASSERT_CLOSE("NIntegrate[4/(1+x^2),{x,0,1},Method->\"TrapezoidalRule\",WorkingPrecision->40]",
                 "Pi", 1e-33);
    /* Newton–Cotes at high precision: the weights are exact rationals, so Simpson
     * is exact for a cubic to the full working precision (not capped at 1e-16). */
    ASSERT_CLOSE("NIntegrate[x^3,{x,0,2},Method->{\"NewtonCotesRule\",\"Points\"->3},"
                 "WorkingPrecision->40]", "4", 1e-33);
    ASSERT_CLOSE("NIntegrate[x^2,{x,0,2},Method->{\"NewtonCotesRule\",\"Points\"->3},"
                 "WorkingPrecision->30]", "8/3", 1e-25);
    ASSERT_CLOSE("NIntegrate[4/(1+x^2),{x,0,1},Method->{\"NewtonCotesRule\",\"Points\"->5},"
                 "WorkingPrecision->30]", "Pi", 1e-22);
    /* The high-precision result actually carries the requested precision. */
    ASSERT_EQ_STR("Precision[NIntegrate[4/(1+x^2),{x,0,1},Method->\"TrapezoidalRule\","
                  "WorkingPrecision->30]] > 25", "True");
}

static void test_unimplemented_method(void) {
    /* A recognised-but-unimplemented method must warn and stay unevaluated, not
     * silently approximate, so the missing implementation is visible. */
    char* s = eval_str("NIntegrate[Exp[Cos[x]],{x,0,10},Method->\"ClenshawCurtisRule\"]");
    ASSERT_MSG(strstr(s, "NIntegrate") != NULL, "unimplemented method should stay unevaluated, got %s", s);
    free(s);
    s = eval_str("NIntegrate[Exp[Cos[x]],{x,0,10},Method->\"LobattoKronrodRule\"]");
    ASSERT_MSG(strstr(s, "NIntegrate") != NULL, "unimplemented method should stay unevaluated, got %s", s);
    free(s);
}

static void test_maxrecursion_zero(void) {
    /* MaxRecursion->0 forbids subdivision: a single GK15 panel on a wiggly
     * integrand is inexact but must still return a finite number. */
    char* s = eval_str("NIntegrate[Exp[Cos[x]],{x,0,10},MaxRecursion->0]");
    ASSERT(strstr(s, "NIntegrate") == NULL);   /* evaluated to a number */
    free(s);
}

static void test_attributes(void) {
    ASSERT_EQ_STR("MemberQ[Attributes[NIntegrate], HoldAll]", "True");
    ASSERT_EQ_STR("MemberQ[Attributes[NIntegrate], Protected]", "True");
    /* HoldAll + localisation: the variable must not leak a global value. */
    ASSERT_EQ_STR("NIntegrate[xfresh^2,{xfresh,0,1}]; xfresh", "xfresh");
}

static void test_unevaluated_forms(void) {
    /* Symbolic (non-numeric) bound stays unevaluated. */
    char* s = eval_str("NIntegrate[x,{x,0,a}]");
    ASSERT(strstr(s, "NIntegrate") != NULL);
    free(s);
}

static void test_memory_loop(void) {
    const char* inputs[] = {
        "NIntegrate[Sin[Sin[x]],{x,0,2}]",
        "NIntegrate[x^x,{x,1,5}]",
        "NIntegrate[1/(1+x^2),{x,0,1}]",
        "NIntegrate[Exp[Cos[x]],{x,0,10},Method->\"GaussKronrodRule\"]",
        "NIntegrate[10^(-5)/(1+10^2 (x-1/2)^2),{x,0,1},AccuracyGoal->8]",
        "NIntegrate[1/Sqrt[x],{x,0,1}]",
        "NIntegrate[Log[x],{x,0,1}]",
        "NIntegrate[1/Sqrt[x(1-x)],{x,0,1}]",
        "NIntegrate[Cos[Log[x]/x]/x,{x,0,1}]",
        "NIntegrate[1/(1+x^2),{x,0,Infinity}]",
        "NIntegrate[Exp[-x^2],{x,-Infinity,Infinity}]",
        "NIntegrate[x,{x,0,a}]",
    };
    for (int i = 0; i < 3; i++)
        for (size_t k = 0; k < sizeof(inputs)/sizeof(inputs[0]); k++) {
            Expr* p = parse_expression(inputs[k]);
            ASSERT(p != NULL);
            Expr* v = evaluate(p);
            expr_free(p);
            expr_free(v);
        }
}

int main(void) {
    symtab_init();
    core_init();

    TEST(test_smooth_finite);
    TEST(test_polynomial_exact);
    TEST(test_negative_and_shifted);
    TEST(test_endpoint_singularities);
    TEST(test_infinite_ranges);
    TEST(test_complex_contours);
    TEST(test_multidim);
    TEST(test_vector_integrand);
    TEST(test_oscillatory);
    TEST(test_oscillatory_singularity);
    TEST(test_montecarlo);
    TEST(test_exclusions_pv);
    TEST(test_high_precision);
    TEST(test_accuracygoal);
    TEST(test_methods_accepted);
    TEST(test_riemann_rule);
    TEST(test_trapezoidal_rule);
    TEST(test_newtoncotes_rule);
    TEST(test_fixed_rule_high_precision);
    TEST(test_unimplemented_method);
    TEST(test_maxrecursion_zero);
    TEST(test_attributes);
    TEST(test_unevaluated_forms);
    TEST(test_memory_loop);

    printf("All nint_tests passed.\n");
    return 0;
}
