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
    /* Semi-infinite *oscillatory* Fourier integrals at high precision, via the
     * Ooura–Mori double-exponential rule (exp-sinh cannot resolve the slowly
     * decaying oscillatory tail).  Regression for the WorkingPrecision branch
     * that used to return a wrong, non-converged value here. */
    ASSERT_CLOSE("NIntegrate[(x Sin[x])/(x^2+4),{x,0,Infinity},WorkingPrecision->100]",
                 "Pi/(2 E^2)", 1e-90);
    ASSERT_CLOSE("NIntegrate[Cos[x]/(x^2+1),{x,0,Infinity},WorkingPrecision->80]",
                 "Pi/(2 E)", 1e-70);
    ASSERT_CLOSE("NIntegrate[Sin[x]/x,{x,0,Infinity},WorkingPrecision->60]",
                 "Pi/2", 1e-50);
    ASSERT_CLOSE("NIntegrate[(x Sin[2 x])/(x^2+1),{x,0,Infinity},WorkingPrecision->80]",
                 "Pi/(2 E^2)", 1e-70);
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

/* Levin collocation: genuine Levin (regular oscillatory, all phase/kernel
 * forms, ultra-high frequency Gauss-Kronrod cannot resolve), the singular-phase
 * fallback, arbitrary precision, multivariate reduction, and Automatic use. */
static void test_levin_rule(void) {
    /* 1-D: amplitude 1 / x, linear and nonlinear phase, complex kernel. */
    ASSERT_CLOSE("NIntegrate[Cos[1000 x],{x,0,1},Method->\"LevinRule\"]",
                 "Sin[1000]/1000", 1e-9);
    ASSERT_CLOSE("NIntegrate[x Cos[100 x],{x,0,1},Method->\"LevinRule\"]",
                 "(Cos[100]+100 Sin[100])/10000 - 1/10000", 1e-9);
    ASSERT_CLOSE("NIntegrate[Exp[I 50 x],{x,0,1},Method->\"LevinRule\"]",
                 "(Exp[I 50]-1)/(I 50)", 1e-9);
    /* Ultra-high frequency: Gauss-Kronrod cannot resolve this in budget; Levin
     * is exact. */
    ASSERT_CLOSE("NIntegrate[Cos[100000 x],{x,0,1},Method->\"LevinRule\"]",
                 "Sin[100000]/100000", 1e-9);
    /* Automatic picks Levin up when the smooth rules fail. */
    ASSERT_CLOSE("NIntegrate[Cos[100000 x],{x,0,1}]", "Sin[100000]/100000", 1e-9);
    /* Singular-phase fallback (1/x diverges at 0): collocation declines, the
     * cascade still integrates it. */
    ASSERT_CLOSE("NIntegrate[Sin[1/x],{x,0,1},Method->\"LevinRule\"]",
                 "0.50417", 1e-3);
    /* Arbitrary precision. */
    ASSERT_CLOSE("NIntegrate[Cos[1000 x],{x,0,1},Method->\"LevinRule\","
                 "WorkingPrecision->30]", "N[Sin[1000]/1000,30]", 1e-25);
    /* Multivariate reduction over a box. */
    ASSERT_CLOSE("NIntegrate[Cos[50 x] Cos[50 y],{x,0,1},{y,0,1},"
                 "Method->\"LevinRule\"]", "(Sin[50]/50)^2", 1e-9);
    ASSERT_CLOSE("NIntegrate[Cos[3 x] Cos[1000 y],{x,0,1},{y,0,1},"
                 "Method->\"LevinRule\"]", "(Sin[3]/3)(Sin[1000]/1000)", 1e-9);
    /* The motivating separable example with a singular oscillatory outer axis. */
    ASSERT_CLOSE("NIntegrate[Sin[1/x] Cos[1000 y],{x,0,1},{y,0,1},"
                 "Method->\"LevinRule\"]", "0.0004168", 1e-6);
}

/* NIntegrate localises its variables and evaluates the integrand symbolically
 * once before sampling, which is the reference's order. A held wrapper such
 * as f /. b -> 1 used to reach the arbitrary-precision sampler unevaluated
 * and gave 0. */
static void test_symbolic_preevaluation(void) {
    assert_eval_eq("Abs[NIntegrate[y^(b - 1)/(1 + y)^2 /. b -> 1, {y, 0, Infinity},"
                   " WorkingPrecision -> 20] - 1] < 10^-15", "True", 0);
    assert_eval_eq("Abs[NIntegrate[y^(b - 1)/(1 + y)^2 /. b -> 3/2, {y, 0, Infinity},"
                   " WorkingPrecision -> 20] - Pi/2] < 10^-15", "True", 0);
    /* a localised variable keeps its outer value */
    assert_eval_eq("Abs[Block[{y = 5}, NIntegrate[y^2, {y, 0, 1}]] - 1/3] < 10^-12", "True", 0);
}

/* Double-exponential tail truncation must not mistake an interior zero of the
 * integrand for the tail.  A squared Bessel / trig factor has a double zero at
 * every root, so a transformed term landing near one is tiny and the next is
 * large; the drivers used to read that dip as "decayed, then noise" and cut
 * away everything beyond it -- the bulk of the integral -- and report the
 * truncated sum as converged (0.636355 for 0.642638 on [0, oo); 9.7e-15 at
 * WorkingPrecision -> 20 on [0, 60]).  References: closed forms, else 40-digit
 * piecewise Gauss-Legendre. */
static void test_de_interior_zeros(void) {
    const char* k0 = "N[2 EllipticK[4/5]/(Pi Sqrt[5]), 40]";   /* 0.6426376817731... */
    /* Machine: exp-sinh half line, and a finite range. */
    ASSERT_CLOSE("NIntegrate[Exp[-x] BesselJ[0, x]^2, {x, 0, Infinity}]", k0, 1e-10);
    ASSERT_CLOSE("NIntegrate[Exp[x] BesselJ[0, x]^2, {x, -Infinity, 0}]", k0, 1e-10);
    ASSERT_CLOSE("NIntegrate[Exp[-x] BesselJ[0, x]^2, {x, 0, 20}]",
                 "0.642637681745200618177181786944", 1e-10);
    ASSERT_CLOSE("NIntegrate[Exp[-x] BesselJ[0, x] BesselJ[1, x], {x, 0, Infinity}]",
                 "0.178681159113437763681178671163", 1e-10);
    ASSERT_CLOSE("NIntegrate[Exp[-x/3] BesselJ[0, 2 x] BesselJ[1, x], {x, 0, Infinity}]",
                 "0.0267677840697094977663843343592", 1e-11);
    ASSERT_CLOSE("NIntegrate[Exp[-x] Sin[x]^2, {x, 0, Infinity}]", "2/5", 1e-10);
    ASSERT_CLOSE("NIntegrate[(SinIntegral[x] - Pi/2) Exp[-x], {x, 0, Infinity}]",
                 "-Pi/4", 1e-10);
    /* WorkingPrecision -> 20: tanh-sinh (finite) and exp-sinh (half line). */
    ASSERT_CLOSE("NIntegrate[Exp[-x] BesselJ[0, x]^2, {x, 0, Infinity}, "
                 "WorkingPrecision -> 20]", k0, 1e-18);
    /* 18 digits is past the default PrecisionGoal (WorkingPrecision - 10), so
     * ask for it: the graded split of this endpoint-heavy range meets the goal
     * it is given and stops there */
    ASSERT_CLOSE("NIntegrate[Exp[-x] BesselJ[0, x]^2, {x, 0, 60}, "
                 "WorkingPrecision -> 20, PrecisionGoal -> 18]", k0, 1e-18);
    ASSERT_CLOSE("NIntegrate[Exp[-x] BesselJ[0, x]^2, {x, 0, 20}, "
                 "WorkingPrecision -> 20]", "0.642637681745200618177181786944", 1e-18);
    ASSERT_CLOSE("NIntegrate[BesselJ[0, x], {x, 0, 20}, WorkingPrecision -> 20]",
                 "1.05837882142112775848794021926", 1e-18);
    ASSERT_CLOSE("NIntegrate[Exp[-x] BesselJ[0, x], {x, 0, 20}, WorkingPrecision -> 20]",
                 "0.707106781077612320576590651691", 1e-18);
    ASSERT_CLOSE("NIntegrate[Exp[-x] BesselJ[0, x] BesselJ[1, x], {x, 0, Infinity}, "
                 "WorkingPrecision -> 20]", "0.178681159113437763681178671163", 1e-18);
    ASSERT_CLOSE("NIntegrate[(SinIntegral[x] - Pi/2) Exp[-x], {x, 0, Infinity}, "
                 "WorkingPrecision -> 20]", "N[-Pi/4, 40]", 1e-18);
}

/* ---- Integrals outside the double range, and the MPFR driver gauges ----
 * Expected values are independent of NIntegrate: closed forms, an exact
 * coefficient extraction (the Laurent coefficient below was computed with
 * exact integer polynomial arithmetic), and a separate 50-digit quadrature of
 * x^x^x split at its boundary layer. */

/* |result / expected - 1| < tol, the result computed once. */
static bool rel_close(const char* input, const char* expected, const char* tol) {
    char buf[4096];
    snprintf(buf, sizeof buf,
             "With[{r = (%s)}, NumberQ[r] && N[Abs[r/(%s) - 1], 20] < %s]",
             input, expected, tol);
    char* s = eval_str(buf);
    bool ok = (strcmp(s, "True") == 0);
    if (!ok) fprintf(stderr, "  rel_close FAIL: %s vs %s (rel tol %s) => %s\n",
                     input, expected, tol, s);
    free(s);
    return ok;
}
#define ASSERT_REL(input, expected, tol)                                     \
    ASSERT_MSG(rel_close((input), (expected), (tol)),                        \
               "%s ~= %s (rel tol %s)", (input), (expected), (tol))

#define XXX_TRUE "1.10266499936194065557541124033087278915258846709*10^36300"
#define CONTOUR_BODY "(((x+1)^2000 (x^2+x+1)^1000 (x^4+x^3+x^2+x+1)^500/x^3001 " \
                     "/. x -> E^(I phi)) I E^(I phi))"
/* Im[contour integral]/(2 Pi) = coefficient of x^3000 in the numerator. */
#define CONTOUR_TRUE "3.97394226558004303969667626329*10^1426"

static void test_beyond_double_range(void) {
    /* WorkingPrecision: the tail/convergence gauges no longer pass through
     * doubles (an overflowing term was Inf and Inf <= reltol Inf "converged";
     * this returned 6.6e25591). */
    ASSERT_REL("NIntegrate[x^x^x, {x, 1, 6}, WorkingPrecision -> 30]",
               XXX_TRUE, "10^-25");
    ASSERT_REL("NIntegrate[x^x^x, {x, 1, 6}, WorkingPrecision -> 30, "
               "PrecisionGoal -> 7, MaxRecursion -> 20, Method -> \"DoubleExponential\"]",
               XXX_TRUE, "10^-7");
    ASSERT_REL("NIntegrate[x^x^x, {x, 1, 6}, WorkingPrecision -> 50, PrecisionGoal -> 40]",
               XXX_TRUE, "10^-40");
    /* Machine precision: a sample beyond 1.8e308 re-runs the integral on the
     * MPFR path at machine-equivalent precision (was unevaluated). */
    ASSERT_REL("NIntegrate[x^x^x, {x, 1, 6}]", XXX_TRUE, "10^-10");
    ASSERT_EQ_STR("Precision[NIntegrate[x^x^x, {x, 1, 6}]] < 17", "True");
    ASSERT_REL("NIntegrate[Exp[1000 x], {x, 0, 1}]", "(E^1000 - 1)/1000", "10^-10");
    /* The same integral after x^x -> z: a boundary layer turned into a
     * ProductLog integrand over [1, 6^6] (was 0. and 7.6e35391). */
    ASSERT_REL("With[{s = (Log[z]/ProductLog[Log[z]])^z D[Log[z]/ProductLog[Log[z]], z]}, "
               "NIntegrate[s, {z, 1, 6^6}]]", XXX_TRUE, "10^-10");
    ASSERT_REL("With[{s = (Log[z]/ProductLog[Log[z]])^z D[Log[z]/ProductLog[Log[z]], z]}, "
               "NIntegrate[s, {z, 1, 6^6}, WorkingPrecision -> 30, PrecisionGoal -> 20]]",
               XXX_TRUE, "10^-20");
    /* A complex machine integrand whose powers overflow (Overflow[] samples
     * used to be dropped, giving a silent 0.). */
    ASSERT_REL("Im[NIntegrate[" CONTOUR_BODY ", {phi, 0, 2 Pi}]]/(2 Pi)",
               CONTOUR_TRUE, "10^-10");
    /* An endpoint singularity is not an overflow of the bulk: still the
     * oscillatory endpoint transform, at machine precision. */
    ASSERT_CLOSE("NIntegrate[Cos[Log[x]/x]/x, {x, 0, 1}]", "0.32336743", 1e-5);
    ASSERT_EQ_STR("Precision[NIntegrate[Cos[Log[x]/x]/x, {x, 0, 1}]]", "MachinePrecision");
}

static void test_mpfr_quadrature_gauges(void) {
    /* A side that starts at a deep zero and climbs: (1 + Cos[x])^2000 vanishes
     * to order 4000 at the centre of [0, 2 Pi]; the cuts used to read the
     * other side's peak and dropped this side's bulk (half the integral).
     * Closed form: 2 Pi Binomial[2 n, n]/2^n. */
    ASSERT_REL("NIntegrate[(1 + Cos[x])^2000, {x, 0, 2 Pi}, WorkingPrecision -> 30]",
               "2 Pi Binomial[4000, 2000]/2^2000", "10^-25");
    ASSERT_REL("Im[NIntegrate[" CONTOUR_BODY ", {phi, 0, 2 Pi}, WorkingPrecision -> 30]]/(2 Pi)",
               CONTOUR_TRUE, "10^-19");   /* default PrecisionGoal 20 */
    /* Exact endpoints at WorkingPrecision: 2 Pi as a double cut 2.4e-16 off
     * the period (relative error 1e-16 at any precision). */
    ASSERT_REL("NIntegrate[Exp[Cos[x]], {x, 0, 2 Pi}, WorkingPrecision -> 40]",
               "2 Pi BesselI[0, 1]", "10^-35");
    /* Real contour nodes at WorkingPrecision run each segment in MPFR (the
     * segments were computed in doubles). */
    ASSERT_REL("NIntegrate[Exp[Cos[x]], {x, 0, Pi, 2 Pi}, WorkingPrecision -> 40]",
               "2 Pi BesselI[0, 1]", "10^-35");
    ASSERT_EQ_STR("Precision[NIntegrate[Exp[Cos[x]], {x, 0, Pi, 2 Pi}, "
                  "WorkingPrecision -> 40]] >= 40", "True");
    /* The half-line Fourier probe flattened Times into the shared integrand:
     * Exp[-x^2] became E^((-1)^(x^2)) for every later sample. */
    ASSERT_REL("NIntegrate[Exp[-x^2], {x, 0, Infinity}, WorkingPrecision -> 50]",
               "Sqrt[Pi]/2", "10^-45");
    ASSERT_REL("NIntegrate[Exp[-x], {x, 0, Infinity}, WorkingPrecision -> 50]",
               "1", "10^-45");
}

static void test_trapezoidal_periodic(void) {
    /* Coarse levels agree by aliasing: 1 and 2 panels sample Sin[x]^2 only at
     * its zeros.  Refinement starts from 4 panels. */
    ASSERT_CLOSE("NIntegrate[Sin[x]^2, {x, 0, 2 Pi}, Method -> \"Trapezoidal\"]", "Pi", 1e-9);
    ASSERT_CLOSE("NIntegrate[Sin[x]^2, {x, 0, 2 Pi}, Method -> \"TrapezoidalRule\"]", "Pi", 1e-9);
    /* "Trapezoidal" is plain doubling (no Romberg table, which mixes in the
     * aliased coarse levels of a periodic integrand); MaxRecursion -> r allows
     * 2^(r+3) panels. */
    ASSERT_REL("NIntegrate[Exp[Cos[x]], {x, 0, 2 Pi}, Method -> \"Trapezoidal\", "
               "WorkingPrecision -> 40, PrecisionGoal -> 30, MaxRecursion -> 3]",
               "2 Pi BesselI[0, 1]", "10^-35");
    /* Exact once the node count exceeds the Laurent degree 3000: two levels
     * agree to 400 digits (stopped at 2 panels and returned 1.5e1428). */
    ASSERT_REL("Im[NIntegrate[" CONTOUR_BODY ", {phi, 0, 2 Pi}, Method -> \"Trapezoidal\", "
               "WorkingPrecision -> 400, PrecisionGoal -> 16, MaxRecursion -> 10]]/(2 Pi)",
               CONTOUR_TRUE, "10^-28");
}

/* ---- Strong algebraic endpoint singularities, option values, non-numeric
 * integrands ---- */
static void test_power_singularity(void) {
    /* x^(a-1) E^-x on [0, 1], a = 1/1000: half the integral lies below 1e-301.
     * Reference: the closed form Gamma[a] - Gamma[a, 1]. */
    ASSERT_REL("NIntegrate[x^(-999/1000) E^-x, {x, 0, 1}]",
               "Gamma[1/1000] - Gamma[1/1000, 1]", "10^-12");
    ASSERT_REL("NIntegrate[x^(-999/1000) E^-x, {x, 0, 1}, WorkingPrecision -> 40]",
               "Gamma[1/1000] - Gamma[1/1000, 1]", "10^-28");
    /* at the upper endpoint, next to a non-zero end: (1-x)^(-9/10) = t^(-9/10) */
    ASSERT_REL("NIntegrate[(1 - x)^(-9/10), {x, 0, 1}, WorkingPrecision -> 30]", "10", "10^-19");
    ASSERT_REL("NIntegrate[(1 - x)^(-9/10) Cos[x], {x, 0, 1}, WorkingPrecision -> 30]",
               "10 NIntegrate[Cos[1 - u^10], {u, 0, 1}, WorkingPrecision -> 40]", "10^-19");
}

static void test_option_values_and_nonnumeric(void) {
    /* option values are evaluated (NIntegrate is HoldAll) */
    ASSERT_CLOSE("Block[{i = 1}, NIntegrate[Exp[y], {y, 0, 1}, PrecisionGoal -> i]]",
                 "E - 1", 1e-6);
    ASSERT_REL("Block[{w = 20}, NIntegrate[Exp[y], {y, 0, 1}, WorkingPrecision -> w + 5]]",
               "E - 1", "10^-14");
    /* no numeric value: unevaluated (was 0.) */
    ASSERT_EQ_STR("Quiet[NIntegrate[ff[x], {x, 0, 1}]]", "NIntegrate[ff[x], {x, 0, 1}]");
    ASSERT_EQ_STR("Quiet[NIntegrate[ff[x], {x, 0, 1}, WorkingPrecision -> 30]]",
                  "NIntegrate[ff[x], {x, 0, 1}, WorkingPrecision -> 30]");
    /* machine-number values at WorkingPrecision: stops at their accuracy
     * (precw) instead of refining level after level */
    /* an under-resolved oscillation also stalls for a few levels, at
     * differences of the integral's own size: not the noise floor */
    ASSERT_REL("NIntegrate[Cos[1000 Sin[x]], {x, 0, Pi}, WorkingPrecision -> 40]",
               "Pi BesselJ[0, 1000]", "10^-25");
    ASSERT_EQ_STR("Module[{h}, h[x_?NumericQ] := N[Exp[x]]; "
                  "Abs[Quiet[NIntegrate[h[x], {x, 0, 1}, WorkingPrecision -> 25]] - (E - 1)] < 10^-14]",
                  "True");
}

static void test_polylog_integrand(void) {
    /* PolyLog near z = 1 at 60 digits; reference: the closed form
     * of Integrate[x^2 PolyLog[3, 1/(x+1)], {x, 0, 1}] */
    ASSERT_REL("NIntegrate[x^2 PolyLog[3, 1/(x + 1)], {x, 0, 1}, WorkingPrecision -> 60]",
               "-17/108 - Pi^2/36 + 35/27 Log[2] - Pi^2 Log[2]/18 - 4/9 Log[2]^2 + Log[2]^3/9 + Zeta[3]/4",
               "10^-45");
}

/* ---- power-law phases: integrated in the phase variable ---- *
 * Closed forms: Sin[1] - CosIntegral[1]; Gamma[4/3]/2 (= Gamma[1 + 1/n]
 * Sin[Pi/(2n)], n = 3); int_0^1 x Cos[2/x^2] dx = int_2^oo Cos[t]/t^2 dt
 * = Cos[2]/2 - Pi/2 + SinIntegral[2]. */
static void test_phase_substitution(void) {
    ASSERT_REL("NIntegrate[Sin[1/x], {x, 0, 1}]", "Sin[1] - CosIntegral[1]", "10^-8");
    ASSERT_REL("NIntegrate[Sin[1/x], {x, 0, 1}, WorkingPrecision -> 40]",
               "Sin[1] - CosIntegral[1]", "10^-28");
    ASSERT_REL("NIntegrate[Sin[x^3], {x, 0, Infinity}]", "Gamma[4/3]/2", "10^-8");
    ASSERT_REL("NIntegrate[Sin[x^3], {x, 0, Infinity}, WorkingPrecision -> 40]",
               "Gamma[4/3]/2", "10^-28");
    ASSERT_REL("NIntegrate[x Cos[2/x^2], {x, 0, 1}, WorkingPrecision -> 30]",
               "Cos[2]/2 - Pi/2 + SinIntegral[2]", "10^-19");
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
    TEST(test_levin_rule);
    TEST(test_symbolic_preevaluation);
    TEST(test_de_interior_zeros);
    TEST(test_beyond_double_range);
    TEST(test_mpfr_quadrature_gauges);
    TEST(test_trapezoidal_periodic);
    TEST(test_power_singularity);
    TEST(test_option_values_and_nonnumeric);
    TEST(test_polylog_integrand);
    TEST(test_phase_substitution);

    printf("All nint_tests passed.\n");
    return 0;
}
