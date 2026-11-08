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
#include "gkadapt.h"

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
    /* Bounded amplitude: no magnitude blow-up to detect, only the unbounded
     * phase (Gauss-Kronrod used to settle on 0.32397, silently). */
    ASSERT_CLOSE("NIntegrate[Cos[Log[x]/x],{x,0,1}]", "0.32473401260165490375", 1e-7);
    /* Oscillatory singularities at BOTH ends: split at the midpoint, each
     * half mapped at its own end (one map over the whole range gave -0.41956). */
    ASSERT_CLOSE("NIntegrate[Sin[1/(x (1-x))],{x,0,1}]", "-0.42056493771992560235", 1e-7);
}

/* The same integrals at WorkingPrecision: the MPFR counterpart of the
 * exponential endpoint map, lobes between MPFR-refined zeros, Wynn's epsilon.
 * The cascade of non-oscillatory rules returned -0.1186 for C1 at
 * WorkingPrecision -> 60 (true 0.3234), and the shortfall warning was
 * withheld. */
static void test_oscillatory_singularity_wp(void) {
    ASSERT_CLOSE("NIntegrate[Cos[Log[x]/x]/x,{x,0,1},WorkingPrecision->60,PrecisionGoal->50]",
                 "0.32336743167777876139937008795217044665104662572546966168103644343",
                 1e-48);
    /* at the upper endpoint, on another interval: Pi/2 - Si(1/2) */
    ASSERT_CLOSE("NIntegrate[Sin[1/(2-x)]/(2-x),{x,0,2},WorkingPrecision->40]",
                 "Pi/2 - SinIntegral[1/2]", 1e-28);
    /* both endpoints at WorkingPrecision */
    ASSERT_CLOSE("NIntegrate[Sin[1/(x (1-x))],{x,0,1},WorkingPrecision->40]",
                 "-0.420564937719925602353566105440624804841397291901640261734259",
                 1e-28);
    /* an inexact 60-digit parameter (an optimiser's iterate started from a
     * double): the samples carry fewer bits than the working precision.  The
     * zero refinement keeps its best iterate at that noise floor (it fell
     * back to the double zero, flooring the extrapolation at 1e-33 for all
     * 2000 lobes: 28 s), and the geometry samples run on the compiled
     * machine twin of the integrand.  Under a second now. */
    ASSERT_CLOSE("With[{a = SetPrecision[0.64, 60]}, TimeConstrained["
                 "NIntegrate[x^a Sin[a/(2 - x)], {x, 0, 2}, WorkingPrecision -> 60], 20]]",
                 "0.987080936877956809184028586034308987522867390394703303662265", 1e-45);
    /* a fast oscillation of FIXED frequency is no endpoint singularity: the
     * ordinary rules keep it (and agree with the closed form) */
    ASSERT_CLOSE("NIntegrate[Sin[1000 x],{x,0,1},WorkingPrecision->40]",
                 "(1 - Cos[1000])/1000", 1e-30);
    /* a removable singularity whose machine samples cancel to noise next to
     * the endpoint is not mistaken for an oscillation */
    ASSERT_CLOSE("NIntegrate[(Exp[x]-1-x)/x^2,{x,0,1},WorkingPrecision->40]",
                 "2 - E + ExpIntegralEi[1] - EulerGamma", 1e-30);
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

/* ---- Cauchy contour integrals at WorkingPrecision (compiled complex MPFR) ----
 * The trapezoidal rule on N points of a trig polynomial is exact once N exceeds
 * its band, so each is a coefficient to all digits.  The bodies exercise the
 * compiled complex-MPFR integer power (binary powering), a bignum Rational
 * coefficient ((99/100)^k), and Tan / Csc of a complex argument.  Expected
 * values: exact integer / rational arithmetic in Python (sum C(n,j)
 * C(n-j,n-2j); the series of tan from t' = 1 + t^2; [t^20] (1 - t^2)^(-1/2) =
 * C(20,10)/2^20 for the residue of Csc^21 at 0). */
#define NI_CW "cw[f_, k_, r_, wp_] := Re[NIntegrate[Evaluate[(f /. x -> r E^(I t)) E^(-I k t)], " \
              "{t, 0, 2 Pi}, Method -> \"Trapezoidal\", WorkingPrecision -> wp, "             \
              "PrecisionGoal -> 16, MaxRecursion -> 14]]/(2 Pi r^k); "
#define NI_C200 "9168205127076471545205440477960770886979628912652148138101796804930035476232905514720690661207"
static void test_contour_wp_compiled(void) {
    ASSERT_EQ_STR(NI_CW "Abs[cw[(1 + x + x^2)^200, 200, 1, 30]/" NI_C200 " - 1] < 10^-20", "True");
    ASSERT_EQ_STR(NI_CW "Abs[cw[(1 + x + x^2)^200, 200, 99/100, 30]/" NI_C200 " - 1] < 10^-20", "True");
    /* negative powers: [x^300] of the coin-change series = 11080130 */
    ASSERT_EQ_STR(NI_CW "Abs[cw[1/((1 - x) (1 - x^2) (1 - x^3) (1 - x^5) (1 - x^7) (1 - x^11)), "
                  "300, 9/10, 30] - 11080130] < 10^-12", "True");
    ASSERT_EQ_STR(NI_CW "Abs[cw[Tan[x], 41, 3/2, 30]/(1410211493828985228276049834684/"
                  "121699582862361447435141825020548828125) - 1] < 10^-20", "True");
    ASSERT_EQ_STR("Abs[Re[NIntegrate[Evaluate[(1/Sin[x]^21 /. x -> E^(I t)) E^(I t)], {t, 0, 2 Pi}, "
                  "Method -> \"Trapezoidal\", WorkingPrecision -> 30, PrecisionGoal -> 16, "
                  "MaxRecursion -> 14]]/(2 Pi) - 46189/262144] < 10^-20", "True");
}

/* ---- Fourier integrals and smooth WorkingPrecision integrals ----
 * Machine precision: a well-aligned amp(x) {Sin|Cos}[w x] on [0, oo) runs the
 * Ooura-Mori rule (it used to fall to integrate-between-the-zeros, 2e-9 off on
 * Sin[1/x]).  Closed forms: int_0^oo Cos[x]/(1+x)^2 = 1 - Sin[1] Ci(1) - Cos[1]
 * (Pi/2 - Si(1)) (by parts); Sin[1] - Ci(1); Gamma[4/3]/2.  WorkingPrecision:
 * the bounded Clenshaw-Curtis probe on analytic integrands, the weak-oscillation
 * Levin gate (first lobe of Sin[x^3]); references from 50-digit mpmath quad. */
static void test_fourier_and_smooth_wp(void) {
    ASSERT_REL("NIntegrate[Cos[x]/(1 + x)^2, {x, 0, Infinity}]",
               "1 - Sin[1] CosIntegral[1] - Cos[1] (Pi/2 - SinIntegral[1])", "10^-12");
    ASSERT_REL("NIntegrate[Sin[1/x], {x, 0, 1}]", "Sin[1] - CosIntegral[1]", "10^-12");
    ASSERT_REL("NIntegrate[Sin[x^3], {x, 0, Infinity}]", "Gamma[4/3]/2", "10^-12");
    ASSERT_REL("NIntegrate[x PolyLog[3, -x], {x, 0, 1}, WorkingPrecision -> 30]",
               "-0.30765458032881955246584991473604059788450112197684", "10^-27");
    ASSERT_REL("NIntegrate[Sin[x^3], {x, 0, Pi^(1/3)}, WorkingPrecision -> 40]",
               "0.59096435190956840433899567161260526283135974435192", "10^-30");
    /* an endpoint singularity: the probe declines, tanh-sinh as before */
    ASSERT_REL("NIntegrate[Sqrt[x] Log[x], {x, 0, 1}, WorkingPrecision -> 30]", "-4/9", "10^-27");
    ASSERT_REL("NIntegrate[BesselJ[0, x], {x, 0, 20}, WorkingPrecision -> 30]",
               "NIntegrate[BesselJ[0, x], {x, 0, 20}, WorkingPrecision -> 30, "
               "Method -> \"DoubleExponential\"]", "10^-27");
}

/* ---- Narrow features at WorkingPrecision: the recentred DE rule ----
 * Closed forms: int_0^oo Exp[-x^2 - a^2/x^2] = Sqrt[Pi]/2 E^(-2a) (substitute
 * u = x - a/x); int_0^oo x^n E^-x = n!; a Lorentzian integrates to ArcTan;
 * the near-pole 1/(eps^2 + x^2) on [0, 1] to ArcTan[1/eps]/eps.  Exp-sinh and
 * tanh-sinh needed 10^4 samples (or did not converge: the Lorentzian on
 * [0, oo) with eps = 10^-3 returned 0.59 relative error after 20 s); the
 * recentred rule takes a few hundred. */
static void test_narrow_features_wp(void) {
    ASSERT_REL("NIntegrate[E^(-x^2 - 1/x^2), {x, 0, Infinity}, WorkingPrecision -> 40]",
               "Sqrt[Pi]/2 E^-2", "10^-28");
    ASSERT_REL("NIntegrate[E^(-x^2 - 100/x^2), {x, 0, Infinity}, WorkingPrecision -> 40]",
               "Sqrt[Pi]/2 E^-20", "10^-28");
    ASSERT_REL("NIntegrate[E^(-x^2 - 10^6/x^2), {x, 0, Infinity}, WorkingPrecision -> 40]",
               "Sqrt[Pi]/2 E^-2000", "10^-28");
    ASSERT_REL("NIntegrate[E^(-x^2 - 10^8/x^2), {x, 0, Infinity}, WorkingPrecision -> 40]",
               "Sqrt[Pi]/2 E^-20000", "10^-28");
    ASSERT_REL("NIntegrate[x^50 E^-x, {x, 0, Infinity}, WorkingPrecision -> 40]",
               "50!", "10^-28");
    ASSERT_REL("NIntegrate[x^1000 E^-x, {x, 0, Infinity}, WorkingPrecision -> 40]",
               "1000!", "10^-28");
    /* Lorentzian peaks on the half line, one and two of them (the second
     * needs the split at the first) */
    ASSERT_REL("NIntegrate[1/((x - 100)^2 + 1), {x, 0, Infinity}, WorkingPrecision -> 40]",
               "Pi/2 + ArcTan[100]", "10^-28");
    ASSERT_REL("NIntegrate[1/((x - 100)^2 + 10^-6), {x, 0, Infinity}, WorkingPrecision -> 40]",
               "(Pi/2 + ArcTan[10^5]) 10^3", "10^-28");
    ASSERT_REL("NIntegrate[1/((x - 100)^2 + 10^-6) + 1/((x - 300)^2 + 10^-6), "
               "{x, 0, Infinity}, WorkingPrecision -> 40]",
               "(Pi + ArcTan[10^5] + ArcTan[3 10^5]) 10^3", "10^-28");
    /* finite intervals: an interior Lorentzian, two of them, a near-pole at
     * either endpoint */
    ASSERT_REL("NIntegrate[1/((x - 1/3)^2 + 10^-8), {x, 0, 1}, WorkingPrecision -> 40]",
               "(ArcTan[(2/3) 10^4] + ArcTan[(1/3) 10^4]) 10^4", "10^-28");
    ASSERT_REL("NIntegrate[1/((x - 1/3)^2 + 10^-16), {x, 0, 1}, WorkingPrecision -> 40]",
               "(ArcTan[(2/3) 10^8] + ArcTan[(1/3) 10^8]) 10^8", "10^-28");
    ASSERT_REL("NIntegrate[1/((x - 1/3)^2 + 10^-8) + 1/((x - 2/3)^2 + 10^-8), {x, 0, 1}, "
               "WorkingPrecision -> 40]",
               "2 (ArcTan[(2/3) 10^4] + ArcTan[(1/3) 10^4]) 10^4", "10^-28");
    ASSERT_REL("NIntegrate[1/(10^-20 + x^2), {x, 0, 1}, WorkingPrecision -> 40]",
               "ArcTan[10^10] 10^10", "10^-28");
    ASSERT_REL("NIntegrate[1/(10^-40 + x^2), {x, 0, 1}, WorkingPrecision -> 40]",
               "ArcTan[10^20] 10^20", "10^-28");
    ASSERT_REL("NIntegrate[1/(10^-20 + (1 - x)^2), {x, 0, 1}, WorkingPrecision -> 40]",
               "ArcTan[10^10] 10^10", "10^-28");
    /* a narrow Gaussian on a peak of mass, plus a broad one */
    ASSERT_REL("NIntegrate[E^(-(x - 3)^2) + 1/((x - 50)^2 + 10^-6), {x, 0, Infinity}, "
               "WorkingPrecision -> 40]",
               "Sqrt[Pi]/2 (1 + Erf[3]) + (Pi/2 + ArcTan[5 10^4]) 10^3", "10^-28");
    /* oscillatory integrands are not features: the rules built for them */
    ASSERT_REL("NIntegrate[Sin[100 x]^2, {x, 0, 1}, WorkingPrecision -> 40]",
               "1/2 - Sin[200]/400", "10^-28");
    ASSERT_REL("NIntegrate[E^-x BesselJ[0, x]^2, {x, 0, Infinity}, WorkingPrecision -> 30]",
               "2 EllipticK[-4]/Pi", "10^-20");
}

/* A removable singularity at an endpoint cancels catastrophically at the
 * deepest tanh-sinh nodes (relative distance 2^-300 at 30 digits): the
 * samples must be recomputed with the bits the cancellation eats.  Closed
 * forms by termwise integration of the Taylor series:
 *   int_0^1 (e^x - 1 - x)/x^2 = sum_{k>=2} 1/(k! (k-1)),
 *   int_0^1 (e^x - 1 - x - x^2/2)/x^3 = sum_{k>=3} 1/(k! (k-2)),
 * (truncated at k = 60, 70: the remainder is below 10^-80), the mirror image
 * at the upper endpoint, and int_0^2 (1 - cos x)/x^2 = Si(2) - (1 - cos 2)/2
 * (by parts).  Before: -188.488 for 0.59962. */
static void test_endpoint_cancellation_wp(void) {
    ASSERT_CLOSE("NIntegrate[(Exp[x] - 1 - x)/x^2, {x, 0, 1}, WorkingPrecision -> 30]",
                 "Sum[1/(k! (k - 1)), {k, 2, 60}]", 1e-24);
    ASSERT_CLOSE("NIntegrate[(Exp[1 - x] - 1 - (1 - x))/(1 - x)^2, {x, 0, 1}, "
                 "WorkingPrecision -> 30]", "Sum[1/(k! (k - 1)), {k, 2, 60}]", 1e-24);
    ASSERT_CLOSE("NIntegrate[(Exp[x] - 1 - x - x^2/2)/x^3, {x, 0, 1}, WorkingPrecision -> 40]",
                 "Sum[1/(k! (k - 2)), {k, 3, 70}]", 1e-33);
    ASSERT_CLOSE("NIntegrate[(1 - Cos[x])/x^2, {x, 0, 2}, WorkingPrecision -> 30]",
                 "SinIntegral[2] - (1 - Cos[2])/2", 1e-24);
    /* A rational constant inside the cancelling sum (prenumericalized at the
     * working precision for the rules, but re-evaluated exactly by the guard;
     * 0.19066889 for 0.19066925 when it was not):
     *   int_0^1 (sin x - x + x^3/6)/x^5 = sum_{k>=2} (-1)^k/((2k+1)! (2k-3)),
     *   int_0^1 (e^x - 1 - x - x^2/2 - x^3/6)/x^4 = sum_{k>=4} 1/(k! (k-3)). */
    ASSERT_CLOSE("NIntegrate[(Sin[x] - x + x^3/6)/x^5, {x, 0, 1}, WorkingPrecision -> 40]",
                 "Sum[(-1)^k/((2 k + 1)! (2 k - 3)), {k, 2, 40}]", 1e-35);
    ASSERT_CLOSE("NIntegrate[(Exp[x] - 1 - x - x^2/2 - x^3/6)/x^4, {x, 0, 1}, "
                 "WorkingPrecision -> 40]", "Sum[1/(k! (k - 3)), {k, 4, 70}]", 1e-35);
    ASSERT_EQ_STR("Precision[NIntegrate[(Exp[x] - 1 - x - x^2/2)/x^3, {x, 0, 1}, "
                  "WorkingPrecision -> 40]] >= 39", "True");
    /* an inexact constant is widened exactly for the guard (it held every
     * deep re-evaluation to its own 45 or 60 digits, and the rule took the
     * resulting 10^127 for 1/6 at face value) */
    ASSERT_CLOSE("NIntegrate[(Exp[x] - 1 - x - N[1/2, 45] x^2)/x^3, {x, 0, 1}, "
                 "WorkingPrecision -> 40]", "Sum[1/(k! (k - 2)), {k, 3, 70}]", 1e-35);
    ASSERT_CLOSE("NIntegrate[(Exp[x] - 1 - x - N[1/2, 60] x^2)/x^3, {x, 0, 1}, "
                 "WorkingPrecision -> 40]", "Sum[1/(k! (k - 2)), {k, 3, 70}]", 1e-35);
    /* no message: the samples are now accurate (the cancelled zeros used to
     * be reported as a 1-bit integrand precision) */
    ASSERT_EQ_STR("Check[NIntegrate[(Exp[x] - 1 - x)/x^2, {x, 0, 1}, WorkingPrecision -> 30]; ok, "
                  "warned]", "ok");
}

/* Interior points where the integrand is not analytic are read off its
 * Log/Abs/fractional-power arguments and become exact breakpoints at
 * WorkingPrecision.  Closed forms: F(x) = x Log|x^2 - 2| - 2x + Sqrt[2]
 * Log|(x + Sqrt[2])/(x - Sqrt[2])|; int Log|x - c| = (x - c)(Log|x - c| - 1);
 * int (x - c) e^x = (x - c - 1) e^x. */
static void test_interior_singular_points_wp(void) {
    ASSERT_CLOSE("NIntegrate[Log[Abs[x^2 - 2]], {x, 0, 4}, WorkingPrecision -> 30]",
                 "4 Log[14] - 8 + Sqrt[2] Log[(4 + Sqrt[2])/(4 - Sqrt[2])]", 1e-27);
    ASSERT_CLOSE("NIntegrate[Log[Abs[x^2 - 2]], {x, 4, 0}, WorkingPrecision -> 30]",
                 "-(4 Log[14] - 8 + Sqrt[2] Log[(4 + Sqrt[2])/(4 - Sqrt[2])])", 1e-27);
    ASSERT_CLOSE("NIntegrate[Log[Abs[x - Pi/3]], {x, 0, 4}, WorkingPrecision -> 30]",
                 "(Pi/3) (Log[Pi/3] - 1) + (4 - Pi/3) (Log[4 - Pi/3] - 1)", 1e-27);
    ASSERT_CLOSE("NIntegrate[Abs[x - 1/3] Exp[x], {x, 0, 1}, WorkingPrecision -> 30]",
                 "2 E^(1/3) - E/3 - 4/3", 1e-27);
}

/* The default AccuracyGoal floor is relative to the integrand's scale: an
 * integrand of size 1e-40 is integrated to its PrecisionGoal, not stopped at
 * an absolute 1e-16.  Closed forms: int_{-1}^1 dx/(c + x^2) = (2/Sqrt[c])
 * ArcTan[1/Sqrt[c]]; int_0^1 Sqrt[x] Log[x] = -4/9; the residue of Gamma^3 at
 * -20, (1/20!^3)((9/2)(H_20 - EulerGamma)^2 + Pi^2/4 + (3/2) H_20^(2)), as a
 * Cauchy integral on |z + 20| = 1/2 (was 3.1745e-54 for 3.18865e-54). */
static void test_tiny_integrand_scale(void) {
    ASSERT_CLOSE("NIntegrate[10^-40/(10^-4 + x^2), {x, -1, 1}]/(10^-40 2 100 ArcTan[100])",
                 "1", 1e-9);
    ASSERT_CLOSE("NIntegrate[10^-30 Sqrt[x] Log[x], {x, 0, 1}]/(-4/9 10^-30)", "1", 1e-8);
    ASSERT_CLOSE("Re[NIntegrate[Gamma[-20 + Exp[I t]/2]^3 Exp[I t]/2, {t, 0, 2 Pi}]/(2 Pi)]"
                 "/((1/20!^3) ((9/2) (HarmonicNumber[20] - EulerGamma)^2 + Pi^2/4"
                 " + (3/2) HarmonicNumber[20, 2]))", "1", 1e-9);
    /* an O(1) integrand that cancels to 0 keeps the absolute floor */
    ASSERT_CLOSE("NIntegrate[Sin[x], {x, 0, 2 Pi}]", "0", 1e-12);
}

/* Machine-precision NIntegrate delivers the digits the machine can, not just
 * PrecisionGoal (true values: mpmath tanh-sinh at 40 digits, or closed forms).
 *  - QK15 carried its nodes and weights to 15 significant digits; the Kronrod
 *    weights summed to 2 - 6e-15, a 3e-15 relative error in EVERY integral:
 *    Sin[x^3] on [0, 3] had 14 digits (3.3e-15), now 1 ulp.
 *  - An integrable endpoint singularity was bisected until the error estimate
 *    met the 1e-8 goal, leaving ~11 digits: (1-x^4)^(1/4) 4.4e-12, Sqrt[x]
 *    3.4e-12, Pi - 2 ArcSin[1 - Sin[x]] 3.5e-12 (a hidden sqrt at both ends).
 *    The endpoint panel is now settled by tanh-sinh after 4 halvings. */
static void test_machine_last_digits(void) {
    ASSERT_REL("NIntegrate[Sin[x^3], {x, 0, 3}]",
               "0.4564240461870738388407399501352734962993`40", "1*^-15");
    ASSERT_REL("NIntegrate[Exp[x], {x, 0, 1}]", "E - 1", "5*^-16");
    ASSERT_REL("NIntegrate[2 (1 - x^4)^(1/4), {x, -1, 1}]",
               "3.708149354602743836867700694390520092435`40", "1*^-15");
    ASSERT_REL("NIntegrate[Pi - 2 ArcSin[1 - Sin[x]], {x, 0, Pi}]",
               "7.294882384504139948018321633528322796869`40", "1*^-15");
    ASSERT_REL("NIntegrate[Sqrt[x], {x, 0, 1}]", "2/3", "1*^-15");
    ASSERT_REL("NIntegrate[Sqrt[1 - x^2], {x, -1, 1}]", "Pi/2", "1*^-15");
    ASSERT_REL("NIntegrate[Sqrt[x] Log[x], {x, 0, 1}]", "-4/9", "1*^-15");
    ASSERT_REL("NIntegrate[x^(-3/4), {x, 0, 1}]", "4", "1*^-15");
    /* a fast oscillation running into the upper end (SIAM problem 6 of
     * appendix D) also takes the tanh-sinh end panel: 13 digits -> 15 */
    ASSERT_REL("NIntegrate[Exp[x + Sin[Exp[Exp[x + 1/3]]]], {x, -1, 1}]",
               "3.070993788451709547418076196696115317568`40", "1*^-14");
}

/* Oscillatory half lines at WorkingPrecision with several frequencies or a
 * Bessel factor: each linear-phase group is rotated off the real axis.
 * Closed forms (Weber / Sonine): int_0^oo J_0(x) cos(b x) = (1 - b^2)^(-1/2)
 * for b < 1, int_0^oo J_0(x) sin(b x) = (b^2 - 1)^(-1/2) for b > 1,
 * int_0^oo J_0 = int_0^oo J_1(x)/x = 1, int_0^oo sin x cos(b x)/x = Pi/2 for
 * |b| < 1.  J_0(x) Cos[x/2] took 26 s in the Fourier rule (its amplitude
 * oscillates) and more than 100 s between the zeros before. */
static void test_rotated_halfline_wp(void) {
    ASSERT_CLOSE("NIntegrate[BesselJ[0, x] Cos[x/2], {x, 0, Infinity}, WorkingPrecision -> 30]",
                 "2/Sqrt[3]", 1e-24);
    ASSERT_CLOSE("NIntegrate[BesselJ[0, x] Sin[2 x], {x, 0, Infinity}, WorkingPrecision -> 30]",
                 "1/Sqrt[3]", 1e-24);
    ASSERT_CLOSE("NIntegrate[BesselJ[0, x], {x, 0, Infinity}, WorkingPrecision -> 30]",
                 "1", 1e-24);
    ASSERT_CLOSE("NIntegrate[BesselJ[1, x]/x, {x, 0, Infinity}, WorkingPrecision -> 30]",
                 "1", 1e-24);
    ASSERT_CLOSE("NIntegrate[Sin[x] Cos[x/2]/x, {x, 0, Infinity}, WorkingPrecision -> 30]",
                 "Pi/2", 1e-24);
    ASSERT_CLOSE("NIntegrate[BesselJ[0, x] Cos[x/2], {x, 0, Infinity}]", "2/Sqrt[3]", 1e-7);
}

/* Exact constants of the integrand numericalized once (ni_prenumericalize):
 * a Cauchy coefficient whose weight (99/100 E^(I t))^-200 evaluates to the
 * exact rational (100/99)^200 times E^(-200 I t), [z^200] 1/(1 - z) = 1 and
 * [z^20] 1/(1 - z)^2 = 21; constants in sums and products of the variable;
 * a constant under a holding head is left alone. */
static void test_prenumericalized_constants_wp(void) {
    ASSERT_REL("Re[NIntegrate[1/(1 - 99/100 Exp[I t]) (99/100 Exp[I t])^(-200), {t, 0, 2 Pi}, "
               "Method -> \"Trapezoidal\", WorkingPrecision -> 40]]/(2 Pi)", "1", "10^-25");
    ASSERT_REL("Re[NIntegrate[1/(1 - 1/2 Exp[I t])^2 (1/2 Exp[I t])^(-20), {t, 0, 2 Pi}, "
               "WorkingPrecision -> 40]]/(2 Pi)", "21", "10^-25");
    ASSERT_REL("NIntegrate[1/3 + Sqrt[2] x + Pi x^2, {x, 0, 1}, WorkingPrecision -> 40]",
               "1/3 + Sqrt[2]/2 + Pi/3", "10^-30");
    ASSERT_REL("NIntegrate[ReleaseHold[Hold[1/3] + x], {x, 0, 1}, WorkingPrecision -> 40]",
               "1/3 + 1/2", "10^-30");
}

static void test_symbolic_swell(void) {
    /* Symbolic processing that swells the integrand is not used when the held
     * integrand samples numerically.  Expanding LegendreP[200, x]^2 gives a
     * degree-400 polynomial whose terms (~1e117) cancel to O(1): the
     * expanded integrand returned garbage with NIntegrate::accgl.  Norms by
     * orthogonality: int_{-1}^{1} P_n^2 = 2/(2n + 1). */
    ASSERT_REL("NIntegrate[LegendreP[200, x]^2, {x, -1, 1}]", "2/401", "10^-10");
    ASSERT_REL("NIntegrate[LegendreP[400, Cos[t]]^2 Sin[t], {t, 0, Pi}, WorkingPrecision -> 30]",
               "2/801", "10^-20");
    /* Argument principle for the zeros of L_n (all real and positive) inside
     * |z| = r: (1/(2 Pi I)) oint L_n'/L_n dz, L_n' = -L_{n-1}^(1).  The count
     * is checked against the real zeros located by NSolve. */
    assert_eval_eq("Round[Im[NIntegrate[(-(LaguerreL[59, 1, 10 E^(I phi)]/LaguerreL[60, 10 E^(I phi)])) "
                   "I 10 E^(I phi), {phi, 0, Pi, 2 Pi}, WorkingPrecision -> 30]]/(2 Pi)] == "
                   "Count[x /. NSolve[LaguerreL[60, x] == 0, x, WorkingPrecision -> 60], r_ /; Re[r] < 10] == 15",
                   "True", 0);
    /* the held form is the same function: a held Expand / Sum samples to the
     * value of the expanded integrand */
    ASSERT_REL("NIntegrate[Expand[(1 + x)^40], {x, 0, 1}, WorkingPrecision -> 30]",
               "(2^41 - 1)/41", "10^-25");
    ASSERT_REL("NIntegrate[Sum[x^k/k!, {k, 0, 60}], {x, 0, 1}, WorkingPrecision -> 30]",
               "Sum[1/(k + 1)!, {k, 0, 60}]", "10^-25");
    ASSERT_REL("NIntegrate[Expand[(x + y)^30], {x, 0, 1}, {y, 0, 1}]", "(2^32 - 2)/(31 32)", "10^-10");
    /* a held derivative has no numeric value at a number: the evaluated form
     * is used, as before */
    ASSERT_CLOSE("NIntegrate[D[Sin[x]^2 Cos[x]^30, x], {x, 0, 1}]", "N[Sin[1]^2 Cos[1]^30, 30]", 1e-15);
}

/* An oscillatory endpoint singularity whose phase scale is a parameter:
 * x^a Sin[a/(2 - x)] on [0, 2].  Three silent wrong answers, each at a few
 * percent of the parameters a (found by NMaximize, which climbed exactly onto
 * the spurious values):
 *   - the sign-change probe of the endpoint shells saw random signs and missed
 *     its threshold by chance (~3% of a), so Gauss-Kronrod settled on an
 *     under-sampled value (a = 0.785..., 2.3956: 2e-3 and 10% off);
 *   - the first lobe of the exponential-endpoint map starts at x = 0, where
 *     x^a is itself singular, and eight bisections did not resolve it
 *     (a = 0.138: 4e-3 off);
 *   - Wynn's extrapolation of the lobe sums met the tolerance by coincidence
 *     at seven lobes (a = 3.139: 7e-5 off), now confirmed by the next lobe.
 * A real breakpoint list {x, 0, 1, 2} runs each piece on the plain real axis,
 * where the oscillatory-endpoint rule applies.  References: mpmath (t = 1/(2-x),
 * tanh-sinh + quadosc) at 30 digits. */
static void test_parametric_oscillatory_endpoint(void) {
    ASSERT_REL("NIntegrate[x^0.7850486024081051 Sin[0.7850486024081051/(2 - x)], {x, 0, 2}]",
               "1.011232310036737021329", "10^-9");
    ASSERT_REL("NIntegrate[x^2.395615544756091 Sin[2.395615544756091/(2 - x)], {x, 0, 2}]",
               "-0.01866731679544101743162", "10^-8");
    ASSERT_REL("NIntegrate[x^0.13806243266311857 Sin[0.13806243266311857/(2 - x)], {x, 0, 2}]",
               "0.4401784246408458778806", "10^-9");
    ASSERT_REL("NIntegrate[x^3.139242625941186 Sin[3.139242625941186/(2 - x)], {x, 0, 2}]",
               "-0.302859995551580093191", "10^-9");
    ASSERT_REL("NIntegrate[x^0.13806243266311857 Sin[0.13806243266311857/(2 - x)], {x, 0, 1, 2}]",
               "0.4401784246408458778806", "10^-9");
}

/* Integrable point singularities at a vertex of a box.  The cubature bisects
 * toward the vertex until the integrand can no longer be evaluated in doubles
 * (3 - Cos[x] - Cos[y] - Cos[z] cancels to 0 near the origin); that region is
 * then frozen at its estimate instead of abandoning the cubature for the
 * iterated rule, which never finished.  At WorkingPrecision the same cubature
 * runs on MPFR-evaluated samples.  Expected values are closed forms: Watson's
 * simple-cubic integral (the Glasser-Zucker Gamma product), and 1/r on the
 * unit square / cube by polar / spherical integration over the triangle /
 * pyramid pieces. */
static void test_vertex_singularity(void) {
    const char* watson = "Sqrt[6]/96 Gamma[1/24] Gamma[5/24] Gamma[7/24] Gamma[11/24]";
    ASSERT_REL("NIntegrate[1/(3 - Cos[x] - Cos[y] - Cos[z]), {x, 0, Pi}, {y, 0, Pi}, {z, 0, Pi}]",
               watson, "10^-6");
    ASSERT_REL("NIntegrate[1/(3 - Cos[x] - Cos[y] - Cos[z]), {x, 0, Pi}, {y, 0, Pi}, {z, 0, Pi}, "
               "WorkingPrecision -> 20, PrecisionGoal -> 8]", watson, "10^-8");
    ASSERT_REL("NIntegrate[1/Sqrt[x^2 + y^2], {x, 0, 1}, {y, 0, 1}]", "2 Log[1 + Sqrt[2]]", "10^-8");
    ASSERT_REL("NIntegrate[1/Sqrt[x^2 + y^2 + z^2], {x, 0, 1}, {y, 0, 1}, {z, 0, 1}]",
               "3/2 Log[2 + Sqrt[3]] - Pi/4", "10^-7");
    ASSERT_REL("NIntegrate[1/Sqrt[x^2 + y^2 + z^2], {x, 0, 1}, {y, 0, 1}, {z, 0, 1}, "
               "WorkingPrecision -> 20, PrecisionGoal -> 8]", "3/2 Log[2 + Sqrt[3]] - Pi/4", "10^-8");
    /* Lattice Green's functions (SIAM 100-digit problem 6 and relatives):
     * singular at k = 0 inside a parabolic cusp, and the formula cancels
     * there (1 - c is exactly 0 in doubles within 1e-8 of the corner).  The
     * vertex guard re-evaluates those samples in MPFR; the Automatic goal
     * keeps two digits of headroom; 8 digits before, whatever the goal. */
    ASSERT_REL("1 - Pi^2/NIntegrate[With[{c = (Cos[k1] + Cos[k2])/2, s = Sin[k1]/5}, "
               "(1 - c)/((1 - c)^2 + s^2)], {k1, 0, Pi}, {k2, 0, Pi}]",
               "0.412749950715419796784923382261052434492", "10^-10");
    /* singular point on a FACE, at the 0 of [-Pi, Pi]: the box is split
     * there so it becomes a vertex of every piece (4-5 digits before) */
    ASSERT_REL("1 - 2 Pi^2/NIntegrate[With[{c = (Cos[k1] + Cos[k2])/2, "
               "s = Sin[k1]/10 + Sin[k2]/5}, (1 - c)/((1 - c)^2 + s^2)], "
               "{k1, 0, Pi}, {k2, -Pi, Pi}]",
               "0.389798040363126905621467224292308950758", "10^-10");
    /* a WorkingPrecision cubature result keeps the requested precision */
    ASSERT_EQ_STR("Precision[NIntegrate[x y, {x, 0, 1}, {y, 0, 1}, WorkingPrecision -> 20]] >= 19",
                  "True");
}

/* QUADPACK QAGS (gk_qags_machine), the 1-D rule of the iterated fallback:
 * QUADPACK's own test integrals, which its epsilon extrapolation settles in
 * a few hundred samples without sampling the singular endpoint. */
static bool qags_logsqrt(void* ctx, double x, double _Complex* out) {
    (void)ctx; *out = log(x) / sqrt(x); return true;
}
static bool qags_xalpha(void* ctx, double x, double _Complex* out) {
    (void)ctx; *out = pow(x, -0.9); return true;     /* integral 10 */
}
static void test_qags_rule(void) {
    GkResult R = gk_qags_machine(qags_logsqrt, NULL, 0.0, 1.0, 0.0, 1e-12, 50);
    ASSERT(R.status == GK_OK);
    ASSERT_MSG(fabs(creal(R.value) + 4.0) < 1e-12, "QAGS log/sqrt: %.17g", creal(R.value));
    ASSERT_MSG(R.abs_err >= fabs(creal(R.value) + 4.0), "QAGS error estimate %g", R.abs_err);
    ASSERT_MSG(R.n_eval < 600, "QAGS log/sqrt took %ld samples", R.n_eval);
    R = gk_qags_machine(qags_xalpha, NULL, 0.0, 1.0, 0.0, 1e-10, 50);
    ASSERT(R.status == GK_OK);
    ASSERT_MSG(fabs(creal(R.value) - 10.0) < 1e-9, "QAGS x^-0.9: %.17g", creal(R.value));
}

/* A 12-digit goal over a box with an integrable vertex singularity (the
 * lattice Green's functions of SIAM problem 6, Watson's integral).  The
 * cubature's error falls only algebraically there, so it stops as soon as its
 * forecast says it cannot finish, and the iterated QAGS rule (with the vertex
 * guard's zone sized to the goal) gives the digits: 3.4 s and a 12-digit
 * answer at the edge of the goal for the 2-D one before, ~30 ms now; the 3-D
 * one missed the goal (10 digits, with a warning).  Truth: hc6 / c6 truth
 * tables (mpmath), Watson's Gamma product. */
static void test_vertex_singularity_high_goal(void) {
    /* corner singularity, [0, Pi]^2 */
    ASSERT_REL("1 - Pi^2/NIntegrate[With[{c = (Cos[k1] + Cos[k2])/2, s = Sin[k1]/5}, "
               "(1 - c)/((1 - c)^2 + s^2)], {k1, 0, Pi}, {k2, 0, Pi}, PrecisionGoal -> 12]",
               "0.412749950715419796784923382261052434492", "3*10^-12");
    /* the same, near the root of SIAM problem 6 (eps = 0.0619...: G = 2) */
    ASSERT_REL("NIntegrate[With[{c = (Cos[k1] + Cos[k2])/2}, (1 - c)/((1 - c)^2 + "
               "(2 6191395447399094/10^17 Sin[k1])^2)], {k1, 0, Pi}, {k2, 0, Pi}, "
               "PrecisionGoal -> 12]/Pi^2", "2", "3*10^-12");
    /* singular point on a face: the box is split there first */
    ASSERT_REL("1 - 2 Pi^2/NIntegrate[With[{c = (Cos[k1] + Cos[k2])/2, "
               "s = Sin[k1]/10 + Sin[k2]/5}, (1 - c)/((1 - c)^2 + s^2)], "
               "{k1, 0, Pi}, {k2, -Pi, Pi}, PrecisionGoal -> 12]",
               "0.389798040363126905621467224292308950758", "3*10^-12");
    /* 3-D: Watson's integral, converged (no NIntegrate::accgl) */
    ASSERT_EQ_STR("Check[w = NIntegrate[1/(3 - Cos[x] - Cos[y] - Cos[z]), {x, 0, Pi}, "
                  "{y, 0, Pi}, {z, 0, Pi}, PrecisionGoal -> 12]; True, False]", "True");
    ASSERT_REL("w", "Sqrt[6]/96 Gamma[1/24] Gamma[5/24] Gamma[7/24] Gamma[11/24]", "10^-12");
    /* and fast: the cubature no longer burns its 10^6 samples first */
    ASSERT_EQ_STR("First[AbsoluteTiming[NIntegrate[With[{c = (Cos[k1] + Cos[k2])/2}, "
                  "(1 - c)/((1 - c)^2 + (Sin[k1]/5)^2)], {k1, 0, Pi}, {k2, 0, Pi}, "
                  "PrecisionGoal -> 12]]] < 1", "True");
}

/* Endpoint guard at an endpoint known only to the working precision: the
 * ends of each chord interval are FindRoot roots, so a deep tanh-sinh node,
 * carried exactly at raised precision, can fall outside the domain (Surd of a
 * complex number). The probe is declined at its first non-number instead of
 * doubling the precision eight times, and its messages are probe noise: not
 * printed, not seen by Check. True value: mpmath tanh-sinh (SIAM D7). */
static void test_guard_probe_outside_domain(void) {
    ASSERT_REL("gA[y_] := y^2/32 + Exp[Sin[y]]; "
               "wA[y_] := Surd[-y + Sqrt[-Log[gA[y]]], 3] - Surd[-y - Sqrt[-Log[gA[y]]], 3]; "
               "y1 = y /. FindRoot[gA[y] == 1, {y, -28/10}, WorkingPrecision -> 60]; "
               "yc = y /. FindRoot[-Log[gA[y]] == y^2, {y, -8/10}, WorkingPrecision -> 60]; "
               "y2 = y /. FindRoot[gA[y] == 1, {y, 38/10}, WorkingPrecision -> 60]; "
               "y3 = y /. FindRoot[gA[y] == 1, {y, 445/100}, WorkingPrecision -> 60]; "
               "Quiet[NIntegrate[wA[y], {y, y1, yc, 0}, WorkingPrecision -> 60, PrecisionGoal -> 50] + "
               "NIntegrate[wA[y], {y, y2, y3}, WorkingPrecision -> 60, PrecisionGoal -> 50]]",
               "2.11502988094614658580840342658669512126535126194010611388164",
               "10^-48");
}

/* Nested exact endpoints at WorkingPrecision.  The machine evaluation of
 * Log[2 + Sqrt[3]] is not the correctly rounded double of its value (2 +
 * Sqrt[3] is rounded before the Log), so the endpoint lookup that recovers
 * the working-precision value of a range endpoint from its double missed it,
 * and the range was integrated over a machine endpoint: 16 digits labelled 60.
 * Truths: mpmath at 75 digits. */
static void test_nested_exact_endpoints_wp(void) {
    ASSERT_REL("With[{L = Log[2 + Sqrt[3]]}, NIntegrate[1, {x, -L, L}, WorkingPrecision -> 60]]",
               "2.6339157938496334172500926946159368880539639429350329595369445138",
               "10^-56");
    ASSERT_REL("With[{L = Log[2 + Sqrt[3]]}, NIntegrate[Log[4 - Exp[x]] + x, {x, -L, L}, "
               "WorkingPrecision -> 60, PrecisionGoal -> 50]]",
               "2.3301701728850429909428440792847316044229378153969491194146380214",
               "10^-49");
    static const char* const ends[][3] = {
        { "Log[2 + Sqrt[3]]",
          "0.92135342654514780173556525438019796436893051038251230947655283228",
          "0.83080378398920202654098973214708350262868510109094301827506405817" },
        { "Sqrt[2 + Sqrt[2]]",
          "1.0747374656441173711015082201414248872833265455331395831139171884",
          "0.87827585368019038267881716391579222280612426402738659371107121604" },
        { "Pi/Sqrt[3]",
          "1.0669334232178860175140896465385058803133071111753523567879800224",
          "0.87708583378558952082206773773527168208344280438557651945628917106" },
        { "ArcCosh[2]",
          "0.92135342654514780173556525438019796436893051038251230947655283228",
          "0.83080378398920202654098973214708350262868510109094301827506405817" },
        { "E^(1/3)",
          "0.94906147063328688027389299713641144418995169333600690325364279631",
          "0.8433188790813059725333647064625167354293439722331574439061217275" },
    };
    char in[512];
    for (size_t i = 0; i < sizeof ends / sizeof ends[0]; i++) {
        snprintf(in, sizeof in, "NIntegrate[1/(1 + x^2), {x, 0, %s}, WorkingPrecision -> 60]",
                 ends[i][0]);
        ASSERT_REL(in, ends[i][1], "10^-55");
        snprintf(in, sizeof in, "NIntegrate[Exp[-x^2], {x, 0, %s}, WorkingPrecision -> 60]",
                 ends[i][0]);
        ASSERT_REL(in, ends[i][2], "10^-55");
    }
    /* an intermediate breakpoint, a reflected range and an infinite range */
    ASSERT_REL("NIntegrate[1/(1 + x^2), {x, 0, Sqrt[2 + Sqrt[2]], Pi/Sqrt[3]}, "
               "WorkingPrecision -> 60]",
               "1.0669334232178860175140896465385058803133071111753523567879800224", "10^-55");
    ASSERT_REL("NIntegrate[1/(1 + x^2), {x, Log[2 + Sqrt[3]], 0}, WorkingPrecision -> 60]",
               "-0.92135342654514780173556525438019796436893051038251230947655283228",
               "10^-55");
    ASSERT_REL("NIntegrate[1/(1 + x^2), {x, Sqrt[2 + Sqrt[2]], Infinity}, WorkingPrecision -> 60]",
               "0.49605886115077924812981347149832655481525815415441332737355510772", "10^-55");
    ASSERT_REL("NIntegrate[1, {x, 0, Log[2 + Sqrt[3]]}, {y, 0, Pi/Sqrt[3]}, WorkingPrecision -> 60]",
               "N[Log[2 + Sqrt[3]] Pi/Sqrt[3], 80]", "10^-55");
    /* the machine path is unchanged */
    ASSERT_EQ_STR("With[{L = Log[2 + Sqrt[3]]}, r = NIntegrate[Log[4 - Exp[x]] + x, {x, -L, L}]; "
                  "{MachineNumberQ[r], Abs[r - 2.330170172885043] < 10^-13}]",
                  "{True, True}");
    ASSERT_EQ_STR("MachineNumberQ[NIntegrate[1/(1 + x^2), {x, 0, Log[2 + Sqrt[3]]}]]", "True");
    /* an exact endpoint with only a machine value: the result claims 16
     * digits, never the working precision (precw); an inexact endpoint is
     * taken as its binary value, silently, as N takes it */
    ASSERT_EQ_STR("Module[{g}, N[g[x_], ___] := N[Sqrt[x]]; SetAttributes[g, NumericFunction]; "
                  "r = Quiet[NIntegrate[1, {x, 0, g[2]}, WorkingPrecision -> 60]]; "
                  "Precision[r] < 17 && Abs[r - Sqrt[2]] < 10^-15]",
                  "True");
    ASSERT_EQ_STR("Precision[NIntegrate[1, {x, 0, 1.1}, WorkingPrecision -> 60]] == 60", "True");
}

/* An oscillatory endpoint whose phase outgrows the double range: in the phase
 * coordinate v = g(x) the integral is a Fourier integral on a half line
 * (Cos[Exp[1/x]] on [0, 1] is Cos[u]/(u Log[u]^2) on [E, oo)).  Before, the
 * machine pass spent 8 s on samples of cosines of 10^1000-sized numbers and
 * returned 1 digit; WorkingPrecision -> 60 did not finish.  True values:
 * mpmath quadosc of the substituted integrals. */
static void test_phase_coordinate_endpoint(void) {
    ASSERT_REL("NIntegrate[Cos[Exp[1/x]], {x, 0, 1}]",
               "-0.216061756334172478596923174411721720448", "10^-11");
    /* the upper endpoint, an amplitude, a faster phase */
    ASSERT_REL("NIntegrate[Cos[Exp[1/(1 - x)]], {x, 0, 1}]",
               "-0.216061756334172478596923174411721720448", "10^-11");
    ASSERT_REL("NIntegrate[Sin[Exp[1/x]]/x, {x, 0, 1}]",
               "-0.17949998145734447352", "10^-11");
    ASSERT_REL("NIntegrate[x Cos[Exp[1/x^2]], {x, 0, 1}]",
               "-0.108030878167086239298461587205860860224", "10^-11");
    ASSERT_EQ_STR("First[AbsoluteTiming[NIntegrate[Cos[Exp[1/x]], {x, 0, 1}]]] < 1",
                  "True");
    ASSERT_REL("NIntegrate[Cos[Exp[1/x]], {x, 0, 1}, WorkingPrecision -> 60, "
               "PrecisionGoal -> 50]",
               "-0.216061756334172478596923174411721720447946726377356852789162913963",
               "10^-48");
}

/* NIntegrate::accgl states the tolerance the run was asked for: an Automatic
 * PrecisionGoal is 8 digits at machine precision (the goal the rules use),
 * not the ~14 digits the warning used to test against, which both warned on
 * results that had met their goal and printed a tolerance of ~1e-14 |I|. */
static void test_accgl_states_goal(void) {
    const char* path = "/tmp/mathilda_nint_accgl_stderr.log";
    fflush(stderr);
    int saved = dup(STDERR_FILENO);
    FILE* sink = freopen(path, "w+", stderr);
    ASSERT(sink != NULL);
    char* r = eval_str("NIntegrate[Abs[Sin[1000 x]], {x, 0, 1}, MaxRecursion -> 2]");
    fflush(stderr);
    if (saved != -1) { dup2(saved, STDERR_FILENO); close(saved); }
    free(r);
    FILE* f = fopen(path, "r");
    ASSERT(f != NULL);
    char buf[2048];
    size_t n = fread(buf, 1, sizeof buf - 1, f);
    buf[n] = '\0';
    fclose(f);
    remove(path);
    ASSERT_MSG(strstr(buf, "NIntegrate::accgl") != NULL, "no accgl: %s", buf);
    /* 10^-7.98 x 0.6487 */
    ASSERT_MSG(strstr(buf, "tolerance 6.8") != NULL && strstr(buf, "e-09") != NULL,
               "accgl tolerance is not the goal's: %s", buf);
}

/* The machine geometry twin of an iterated integral at WorkingPrecision is
 * its inner integral at machine precision: the outer rule's probes (an
 * oscillation probe alone takes ~46 samples) each used to run a whole inner
 * integral at the working precision -- 27% of the time of a 2-D lattice
 * Green's function at WorkingPrecision -> 30 (81 s -> 59 s).  The values are
 * unchanged.  True value: mpmath. */
static void test_iterated_wp_machine_twin(void) {
    ASSERT_REL("NIntegrate[Exp[-x y] Sin[x + y], {x, 0, 1}, {y, 0, 2}, WorkingPrecision -> 40]",
               "1.05541194957058443667825859745073286821672065", "10^-35");
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
    TEST(test_oscillatory_singularity_wp);
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
    TEST(test_contour_wp_compiled);
    TEST(test_fourier_and_smooth_wp);
    TEST(test_narrow_features_wp);
    TEST(test_symbolic_swell);
    TEST(test_endpoint_cancellation_wp);
    TEST(test_interior_singular_points_wp);
    TEST(test_tiny_integrand_scale);
    TEST(test_rotated_halfline_wp);
    TEST(test_vertex_singularity);
    TEST(test_qags_rule);
    TEST(test_vertex_singularity_high_goal);
    TEST(test_prenumericalized_constants_wp);
    TEST(test_parametric_oscillatory_endpoint);
    TEST(test_machine_last_digits);
    TEST(test_guard_probe_outside_domain);
    TEST(test_nested_exact_endpoints_wp);
    TEST(test_phase_coordinate_endpoint);
    TEST(test_accgl_states_goal);
    TEST(test_iterated_wp_machine_twin);

    printf("All nint_tests passed.\n");
    return 0;
}
