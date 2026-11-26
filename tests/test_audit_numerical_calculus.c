/* Builtin-audit regressions: numerical calculus (review/issues/; tests/test_audit.h). */
#include "test_audit.h"

/* NSeries#1 (critical): an order of 2^31 truncated to a negative int and the
 * sampler overran its buffer.  Orders past the documented maximum are
 * refused with NSeries::ivar and the call stays unevaluated. */
static void test_nseries_huge_order(void) {
    audit_chk("NSeries#1", "Check[NSeries[Exp[x], {x, 0, 2^31}], fired, NSeries::ivar]", "fired");
    audit_chk("NSeries#1", "Quiet[NSeries[Exp[x], {x, 0, 2^31}]]", "NSeries[E^x, {x, 0, 2147483648}]");
    audit_chk("NSeries#1", "Quiet[NSeries[Exp[x], {x, 0, 2^64}]]",
              "NSeries[E^x, {x, 0, 18446744073709551616}]");
    audit_true("NSeries#1",
               "Abs[SeriesCoefficient[NSeries[Exp[x], {x, 0, 4}], 3] - 1/6] < 10^-10");
}

/* NSeries#2: huge orders were narrowed to int (2^40 -> 0, 2^32+1 -> 1); they
 * are refused since the NSeries#1 fix. */
static void test_nseries_order_narrowing(void) {
    audit_chk("NSeries#2", "Check[NSeries[Exp[x], {x, 0, 2^40}], fired, NSeries::ivar]", "fired");
    audit_chk("NSeries#2", "Quiet[NSeries[Exp[x], {x, 0, 2^32 + 1}]]",
              "NSeries[E^x, {x, 0, 4294967297}]");
}

/* NIntegrate#1: Method -> "PrincipalValue" reads its poles from the interior
 * nodes of {x, a, s, b}; the contour split integrated through them. */
static void test_nintegrate_pv_nodes(void) {
    audit_true("NIntegrate#1",
               "Abs[NIntegrate[1/x, {x, -1, 0, 2}, Method -> \"PrincipalValue\"] - Log[2]] < 10^-8");
    audit_true("NIntegrate#1",
               "Abs[NIntegrate[1/(x - 1/2), {x, 0, 1/2, 1}, Method -> \"PrincipalValue\"]] < 10^-8");
    audit_true("NIntegrate#1",   /* reversed range */
               "Abs[NIntegrate[1/x, {x, 2, 0, -1}, Method -> \"PrincipalValue\"] + Log[2]] < 10^-8");
    /* 1/(x (x - 1)) = 1/(x - 1) - 1/x: the PV over [-1, 3] is -Log[3] */
    audit_true("NIntegrate#1",
               "Abs[NIntegrate[1/(x (x - 1)), {x, -1, 0, 1, 3}, Method -> \"PrincipalValue\"] + Log[3]] < 10^-8");
    /* the Exclusions spelling and plain node splits are unchanged */
    audit_true("NIntegrate#1",
               "Abs[NIntegrate[1/x, {x, -1, 2}, Method -> \"PrincipalValue\", Exclusions -> {0}] - Log[2]] < 10^-8");
    audit_true("NIntegrate#1", "Abs[NIntegrate[x^2, {x, 0, 1, 2}] - 8/3] < 10^-10");
}

/* NSum#1: an infinite lower bound was summed as if -Infinity were a number. */
static void test_nsum_doubly_infinite(void) {
    audit_true("NSum#1", "Abs[NSum[1/(n^2 + 1), {n, -Infinity, Infinity}] - Pi Coth[Pi]] < 10^-8");
    audit_true("NSum#1",
               "Abs[NSum[Exp[-n^2], {n, -Infinity, Infinity}] - N[Sum[Exp[-n^2], {n, -40, 40}], 20]] < 10^-8");
    audit_true("NSum#1", "Abs[NSum[1/n^2, {n, -Infinity, -1}] - Pi^2/6] < 10^-8");
    audit_true("NSum#1",
               "Abs[NSum[1/(n^2 + 1), {n, -Infinity, Infinity}, WorkingPrecision -> 20] - Pi Coth[Pi]] < 10^-15");
    audit_true("NSum#1", "Abs[NSum[1/n^2, {n, 1, Infinity}] - Pi^2/6] < 10^-8");
}

/* NMaximize#1 / NMinimize#1 / NMinimize#2: the search region (into which every
 * trial point is clamped) must lie inside the box constraints. */
static void test_nminimize_box_region(void) {
    audit_chk("NMaximize#1", "NMaximize[{x, Element[x, Integers], x <= 7.5}, x]", "{7.0, {x -> 7}}");
    audit_chk("NMaximize#1", "NMaximize[{-x, Element[x, Integers], x >= 2.5}, x]", "{-3.0, {x -> 3}}");
    audit_chk("NMaximize#1", "NMinimize[{-x, Element[x, Integers], x <= 7.5}, x]", "{-7.0, {x -> 7}}");
    audit_chk("NMinimize#1", "NMinimize[{x^2, x >= 3, x <= 3}, x]", "{9.0, {x -> 3.0}}");
    audit_true("NMinimize#1",
               "Module[{r = NMinimize[{x^2 + (y - 1)^2, x >= 3, x <= 3}, {x, y}]}, "
               "Abs[r[[1]] - 9] < 10^-6 && (x /. r[[2]]) == 3.]");
    audit_chk("NMinimize#2", "NMinimize[{x^2, x > 5}, {x, 0, 1}]", "{25.0, {x -> 5.0}}");
    audit_chk("NMinimize#2", "NMinimize[{x^2, x >= 3}, x]", "{9.0, {x -> 3.0}}");
}

/* NMinimize#3: an unbounded integer coordinate pinned to the edge of the
 * default search region grows the region. */
static void test_nminimize_integer_edge(void) {
    audit_chk("NMinimize#3", "NMinimize[{(x - 50)^2, Element[x, Integers]}, x]", "{0.0, {x -> 50}}");
    audit_chk("NMinimize#3", "NMinimize[{(x - 50)^2, Element[x, Integers], x >= 0}, x]",
              "{0.0, {x -> 50}}");
    audit_chk("NMinimize#3", "NMinimize[{(x - 3)^2, Element[x, Integers]}, x]", "{0.0, {x -> 3}}");
    /* still on the edge after the last growth: unbounded, no minimum reported */
    audit_chk("NMinimize#3", "Check[NMinimize[{x, Element[x, Integers]}, x], fired, NMinimize::ubnd]",
              "fired");
    audit_chk("NMinimize#3", "Quiet[NMinimize[{x, Element[x, Integers]}, x]]",
              "NMinimize[{x, Element[x, Integers]}, x]");
    audit_chk("NMinimize#3", "Quiet[NMaximize[{x, Element[x, Integers]}, x]]",
              "NMaximize[{x, Element[x, Integers]}, x]");
}

/* NMaximize#2: the 1e300 "not a number" sentinel was reported as the optimum. */
static void test_nmaximize_nonnumeric(void) {
    audit_chk("NMaximize#2", "Check[NMaximize[a x - x^2, x], fired, NMaximize::nnum]", "fired");
    audit_chk("NMaximize#2", "Quiet[NMaximize[a x - x^2, x]]", "NMaximize[a x - x^2, x]");
    audit_chk("NMaximize#2", "Check[NMinimize[I x, x], fired, NMinimize::nnum]", "fired");
    audit_chk("NMaximize#2", "NMaximize[2 x - x^2, x]", "{1.0, {x -> 1.0}}");
}

/* FindMaximum#1: a Gradient option is the gradient of f; the internal
 * FindMinimum of -f must receive its negation. */
static void test_findmaximum_gradient(void) {
    audit_true("FindMaximum#1",
               "Module[{r = FindMaximum[-x^2 - y^2 + 2 x + 4 y, {{x, 0}, {y, 0}}, "
               "Gradient -> {-2 x + 2, -2 y + 4}]}, Abs[r[[1]] - 5] < 10^-8 && "
               "Abs[(x /. r[[2]]) - 1] < 10^-6 && Abs[(y /. r[[2]]) - 2] < 10^-6]");
    audit_true("FindMaximum#1",
               "Module[{r = FindMaximum[-x^2 + 2 x, {x, 0}, Gradient -> {-2 x + 2}, "
               "Method -> \"Newton\"]}, Abs[r[[1]] - 1] < 10^-8]");
    audit_true("FindMaximum#1",
               "Module[{r = FindMaximum[-x^2 + 2 x, {x, 0}]}, Abs[(x /. r[[2]]) - 1] < 10^-6]");
}

/* NProduct#1: a zero factor makes the product 0 (Log[0] broke the log-sum). */
static void test_nproduct_zero_factor(void) {
    audit_chk("NProduct#1", "NProduct[k - 3, {k, 1, 5}]", "0.0");
    audit_chk("NProduct#1", "NProduct[0, {k, 1, 3}]", "0.0");
    audit_chk("NProduct#1", "NProduct[1 - 1/k^2, {k, 1, Infinity}]", "0.0");
    audit_true("NProduct#1", "Abs[NProduct[1 - 1/k^2, {k, 2, Infinity}] - 1/2] < 10^-8");
    audit_chk("NProduct#1", "NProduct[k, {k, 1, 5}]", "120.0");
}

/* NRoots#1: the clustering/chopping windows were absolute below |z| = 1. */
static void test_nroots_small_roots(void) {
    audit_chk("NRoots#1", "NRoots[10^20 x^2 - 1 == 0, x]", "x == -1.0*^-10 || x == 1.0*^-10");
    audit_chk("NRoots#1", "NRoots[x^2 - 10^-30 == 0, x]", "x == -1.0*^-15 || x == 1.0*^-15");
    audit_chk("NRoots#1", "NRoots[10^30 x^2 + 1 == 0, x]",
              "x == 0.0 - 1.0*^-15*I || x == 0.0 + 1.0*^-15*I");
    audit_chk("NRoots#1", "NRoots[(x^2 - 10^-20) (x - 1) == 0, x]",
              "x == -1.0*^-10 || x == 1.0*^-10 || x == 1.0");
    /* multiple roots still cluster */
    audit_chk("NRoots#1", "NRoots[(x - 1)^2 (x + 2) == 0, x]", "x == -2.0 || x == 1.0 || x == 1.0");
    audit_chk("NRoots#1", "NRoots[(x - 10^-10)^2 == 0, x]", "x == 1.0*^-10 || x == 1.0*^-10");
    /* an inexact multiple root keeps the scale-free window (eps^(1/m) error) */
    audit_chk("NRoots#1", "NRoots[(x - 1.*10^-3)^3 == 0, x]", "x == 0.001 || x == 0.001 || x == 0.001");
}

/* NRoots#2: a non-numeric constant equation is not False. */
static void test_nroots_nonnumeric_constant(void) {
    audit_chk("NRoots#2", "Check[NRoots[y^2 == 1, x], fired, NRoots::nnum]", "fired");
    audit_chk("NRoots#2", "Quiet[NRoots[y^2 == 1, x]]", "NRoots[y^2 == 1, x]");
    audit_chk("NRoots#2", "NRoots[1 == 0, x]", "False");
    audit_chk("NRoots#2", "NRoots[1 == 1, x]", "True");
}

/* NSolve#1: a disjunction is the union of its branches. */
static void test_nsolve_or(void) {
    audit_chk("NSolve#1", "NSolve[x^2 == 1 || x == 3, x]", "{{x -> -1.0}, {x -> 1.0}, {x -> 3.0}}");
    audit_chk("NSolve#1", "NSolve[x^2 == 1 || x == 3, x, Reals]",
              "{{x -> -1.0}, {x -> 1.0}, {x -> 3.0}}");
    audit_true("NSolve#1", "Length[NSolve[x^2 == 1 || x^3 == 8, x]] == 5");
    audit_chk("NSolve#1", "NSolve[x^2 == 1, x]", "{{x -> -1.0}, {x -> 1.0}}");
}

/* WhenEvent#1: a sign change across a jump of lhs - rhs is not an event. */
static void test_whenevent_jump(void) {
    audit_chk("WhenEvent#1",
              "Reap[NDSolve[{y'[t] == 1, y[0] == 0, WhenEvent[Mod[t, 1] == 1/2, Sow[t]]}, y, "
              "{t, 0, 3.5}]][[2]]",
              "{{0.5, 1.5, 2.5, 3.5}}");
    audit_chk("WhenEvent#1",
              "Reap[NDSolve[{y'[t] == 1, y[0] == 0, WhenEvent[y[t] == 1.5, Sow[t]]}, y, {t, 0, 3}]][[2]]",
              "{{1.5}}");
}

/* NSum#1 follow-up (v0.6726): D#1 made D reject a numeric variable, which
 * exposed NSum taking its Euler-Maclaurin derivatives with the summation
 * variable still bound to a sample point after a failed evaluation -- every
 * nested NSum with a finite inner range flooded D::ivar, so Check took the
 * failure branch.  The values are the closed forms zeta(2)(1 + 1/4 + 1/9) and
 * zeta(3)(Sin[1] + Sin[2]). */
static void test_nsum_nested_no_d_ivar(void) {
    audit_true("NSum#1", "Quiet[Abs[Check[NSum[1/(i j)^2, {i, 1, Infinity}, {j, 1, 3}], 0, D::ivar] - "
                         "N[Zeta[2] (1 + 1/4 + 1/9)]] < 10^-8, NSum::nnum]");
    audit_true("NSum#1", "Quiet[Abs[Check[NSum[Sin[j]/i^3, {i, 1, Infinity}, {j, 1, 2}], 0, D::ivar] - "
                         "N[Zeta[3] (Sin[1] + Sin[2])]] < 10^-8, NSum::nnum]");
}

/* Wave 3 S12/S13: declared argument counts and option checks
 * (src/arity.h, src/arity_table.c).  The evaluator checks a call
 * against its head's declared arity before the C code runs and
 * prints argx/argr/argrx/argt(u)/argb(u)/argm(u), nonopt or optx;
 * the call stays unevaluated. */
static void test_wave3_arity_options(void) {
    audit_chk("FindMinimum#6", "Quiet[Check[FindMinimum[x^2, {x, 1}, Foo -> 1], fired, FindMinimum::optx]]", "fired");
    audit_chk("NMinimize#5", "Quiet[Check[NMinimize[(x - 2)^2, x, Foo -> 1], fired, NMinimize::optx]]", "fired");
    audit_chk("Accuracy#3", "Quiet[Check[Accuracy[1., 2.], fired, Accuracy::argx]]", "fired");
    audit_chk("FindMaximum#2", "Quiet[Check[FindMaximum[-x^2 + 2 x, {x, 0}, Foo -> 1], fired, FindMaximum::optx]]", "fired");
    audit_chk("NSum#5", "Quiet[Check[NSum[1/n^2, {n, 1, Infinity}, Foo -> 1], fired, NSum::optx]]", "fired");
    audit_chk("NDSolve#5", "Quiet[Check[NDSolve[{y'[x] == -y[x], y[0] == 1}, y, {x, 0, 1}, Foo -> 3], fired, NDSolve::optx]]", "fired");
    audit_chk("NDSolveValue#3", "Quiet[Check[NDSolveValue[{y'[x] == y[x], y[0] == 1}, y[1], {x, 0, 1}, Foo -> 1], fired, NDSolveValue::optx]]", "fired");
    audit_chk("NSeries#3", "Quiet[Check[NSeries[Exp[x]], fired, NSeries::argr]]", "fired");
    audit_chk("NSeries#3", "Quiet[NSeries[Exp[x]]]", "NSeries[E^x]");
}

/* Wave 4 (G5): region integrands, constant and infeasible optimisation,
 * NMaximize's infeasible sentinel, NSolve's integer parameters, Fourier of
 * exact constants, NDSolveValue of a conjunction. */
static void test_wave4_numerical(void) {
    audit_chk("NIntegrate#2", "Abs[NIntegrate[Boole[x < 1/3], {x, 0, 1}] - 1/3] < 10^-12", "True");
    audit_chk("NIntegrate#3", "Abs[NIntegrate[Boole[x^2 + y^2 < 1], {x, -1, 1}, {y, -Sqrt[1 - x^2], Sqrt[1 - x^2]}] - Pi] < 10^-6", "True");
    audit_chk("FindMinimum#1", "FindMinimum[5, {x, 1}]", "{5.0, {x -> 1.0}}");
    audit_chk("FindMinimum#3", "Check[FindMinimum[{x^2, x >= 1, x <= 0}, {x, 0.5}], fired, FindMinimum::nsol]", "fired");
    audit_chk("NMaximize#3", "Check[NMaximize[{x, x^2 + 1 <= 0}, x], fired, NMaximize::nsol]", "fired");
    audit_chk("NMaximize#3", "Quiet[NMaximize[{x, x^2 + 1 <= 0}, x]]", "{-Infinity, {x -> Indeterminate}}");
    audit_chk("NSolve#2", "FreeQ[Quiet[NSolve[Sin[x] == 0, x]], C[1.]]", "True");
    audit_chk("NSolve#2", "MemberQ[Quiet[NSolve[Sin[x] == 0, x]], {x -> ConditionalExpression[_ C[1], Element[C[1], Integers]]}, Infinity]", "True");
    audit_chk("Fourier#1", "Max[Abs[Fourier[{Pi, E}] - {4.14356, 0.299326}]] < 10^-5", "True");
    audit_chk("Fourier#1", "Head[First[Fourier[{Sqrt[2], 1}]]]", "Real");
    audit_chk("NDSolveValue#1", "Abs[NDSolveValue[y'[x] == y[x] && y[0] == 1, y[1], {x, 0, 1}] - E] < 10^-5", "True");
}







/* ---- Wave 5 G5 BEGIN (calculus-solve audit fixes) ---- */
static void test_wave5_g5(void) {
    audit_chk("Or#1", "StringContainsQ[ToString[Information[Or]], \"logical OR\"]", "True");
    audit_chk("NSum#4", "Abs[Check[NSum[1/(i j)^2, {i, 1, Infinity}, {j, 1, Infinity}], fired, NSum::nnum] - Zeta[2]^2] < 10^-6", "True");
    audit_chk("FindMinimum#4", "Check[FindMinimum[x, {x, 0}], fired, FindMinimum::cvmit]", "fired");
    audit_chk("FindMinimum#5", "Check[FindMinimum[x^2, {{x, 1}, {x, 2}}], fired, FindMinimum::dupv]", "fired");
    audit_chk("FindMinimum#5", "Quiet[FindMinimum[x^2, {{x, 1}, {x, 2}}]]", "FindMinimum[x^2, {{x, 1}, {x, 2}}]");
    audit_chk("FindMinimum#7", "Check[FindMinimum[f[x], {x, 1}], fired, FindMinimum::nrnum]", "fired");
    audit_chk("FindMinimum#7", "Quiet[FindMinimum[f[x], {x, 1}]]", "FindMinimum[f[x], {x, 1}]");
    audit_chk("NSolve#3", "Head[Check[NSolve[x^2 == a, x], fired, NRoots::nnum]]", "List");
    audit_chk("NSolve#4", "Block[{$MessageList = {}}, NSolve[x^x == 2, x, Reals]; Length[$MessageList]]", "0");
    audit_chk("NSolve#4", "Abs[x - 1.5596104694623694 /. First[NSolve[x^x == 2, x, Reals]]] < 10^-9", "True");
    audit_chk("ND#2", "Check[ND[x^2, {x, 4294967298}, 1], fired, ND::ord]", "fired");
    audit_chk("NDEigenvalues#2", "OrderedQ[TimeConstrained[NDEigenvalues[{-Laplacian[u[x, y], {x, y}], DirichletCondition[u[x, y] == 0, True]}, u[x, y], {x, y} \\[Element] Rectangle[{0, 0}, {1, 1}], 6], 120]]", "True");
    audit_chk("NDSolve#2", "Check[NDSolve[{y'[x] == -y[x], y[0] == 1}, y, {x, 0, 100}, MaxSteps -> 5], fired, NDSolve::mxst]", "fired");
    audit_chk("NDSolve#3", "Check[NDSolve[{y'[x] == -y[x], y[0] == a}, y, {x, 0, 1}], fired, NDSolve::ndinnt]", "fired");
    audit_chk("NDSolve#4", "Check[NDSolve[{y'[x] == -y[x], y[0] == 1}, y, {x, 0, a}], fired, NDSolve::ndnl]", "fired");
    audit_chk("NDSolve#6", "Abs[Check[y[1] /. First[NDSolve[{y'[x] == -y[x], y[0] == 1}, y, {x, 0, 1}, StartingStepSize -> 0.25, Method -> \"RK4\"]], fired] - E^-1.] < 10^-6", "True");
    audit_chk("WhenEvent#4", "Check[(y /. First[NDSolve[{y'[t] == 1, y[0] == 0, WhenEvent[t == 1/2, y[t] -> 0]}, y, {t, 0, 1}]])[1], fired] - 0.5 < 10^-9", "True");
    audit_chk("NDSolveValue#2", "Quiet[Check[NDSolveValue[{y'[x] == y[x], y[0] == 1}, y[1], {x, 0, 1}, MaxSteps -> 3], fired, NDSolveValue::mxst]]", "fired");
    audit_chk("NDSolveValue#4", "Check[NDSolveValue[{y'[x] == y[x], y[0] == 1}, y[1], {x, 0, a}], fired, NDSolveValue::ndnl]", "fired");
    audit_chk("NDSolveValue#5", "Quiet[Check[NDSolveValue[{y'[x] == 1/(x - 0.5), y[0] == 1}, y[1], {x, 0, 1}], fired, NDSolveValue::ndsz]]", "fired");
    audit_chk("WhenEvent#4", "(y /. First[NDSolve[{y'[t] == 1, y[0] == 0, WhenEvent[t == 1/2, y[t] -> 0]}, y, {t, 0, 1}]])[\"Domain\"]", "{{0.0, 1.0}}");
    audit_chk("NProduct#3", "NProduct[2, {3}] === 8.", "True");
    audit_chk("NRoots#3", "Check[NRoots[x^2 == 2], fired, NRoots::argr]", "fired");
    audit_chk("NRoots#4", "Check[NRoots[x^2 == 2, x, Foo -> 1], fired, NRoots::optx]", "fired");
    audit_chk("Fourier#2", "{Check[Fourier[{}], fired, Fourier::fftl], Check[Fourier[5], fired, Fourier::fftl], Check[Fourier[{{1, 2}, {3}}], fired, Fourier::fftl], Check[Fourier[{\"a\", \"b\"}], fired, Fourier::fftl]}", "{fired, fired, fired, fired}");
    audit_chk("Fourier#2", "Quiet[Fourier[{\"a\", \"b\"}]]", "Fourier[{\"a\", \"b\"}]");
    audit_chk("NLimit#2", "{Check[NLimit[Sin[x]/x, x], fired, NLimit::lim], Check[NLimit[Sin[x]/x, x -> 0, Scale -> 0], fired, NLimit::badopt], Check[NLimit[Sin[x]/x, x -> 0, Direction -> 0], fired, NLimit::badopt]}", "{fired, fired, fired}");
    audit_chk("NResidue#1", "Abs[Check[NResidue[1/x, {x, 0}, Radius -> 1/10, MaxRecursion -> 0], fired] - 1] < 10^-12", "True");
    audit_chk("ND#3", "Abs[Check[ND[x, {x, -1/2}, 1, Method -> NIntegrate], fired] - 4/(3 Sqrt[Pi])] < 10^-8", "True");
    audit_chk("NSeries#4", "Max[Abs[Normal[CoefficientList[Normal[NSeries[Exp[x], {x, 0, 2}]] x^2, x]] - {0, 0, 1, 1, 1/2}]] < 10^-10", "True");
}
/* ---- Wave 5 G5 END ---- */

/* Wave 6 G2: dependent inner ranges in NSum, a minimum at the edge of the
 * real domain, decaying oscillation in NLimit, periodic and high-precision
 * events, sub-machine WorkingPrecision, default ND goals, Fourier precision. */
static void test_wave6_g2(void) {
    audit_true("NSum#3", "Abs[NSum[1/(i j)^2, {i, 1, Infinity}, {j, 1, i}] - N[(Zeta[2]^2 + Zeta[4])/2]] < 10^-10");
    audit_true("FindMinimum#2", "Module[{r = Quiet[FindMinimum[Sqrt[x], {x, 1}]]}, "
               "(x /. r[[2]]) >= 0 && r[[1]] == (Sqrt[x] /. r[[2]]) && r[[1]] < 10^-10]");
    audit_true("NLimit#1", "Abs[NLimit[Sin[x]/x, x -> Infinity]] < 10^-12");
    audit_true("NLimit#1", "Quiet[Check[NLimit[Sin[x]/x, x -> Infinity], \"m\"]] =!= \"m\"");
    audit_true("NLimit#1", "Quiet[Check[NLimit[x Sin[x], x -> Infinity], \"m\", NLimit::osc]] === \"m\"");
    audit_chk("WhenEvent#2", "Reap[NDSolve[{y'[t] == 1, y[0] == 0, WhenEvent[Mod[t, 1] == 0, Sow[t]]}, y, {t, 0, 3.5}]][[2]]",
              "{{1.0, 2.0, 3.0}}");
    audit_chk("WhenEvent#2", "Reap[NDSolve[{y'[t] == 1, y[0] == 0, WhenEvent[Mod[t, 1] == 1/2, Sow[t]]}, y, {t, 0, 3.5}]][[2]]",
              "{{0.5, 1.5, 2.5, 3.5}}");
    audit_chk("WhenEvent#3", "Precision[Reap[NDSolve[{y'[t] == y[t], y[0] == 1, WhenEvent[y[t] == 2, Sow[t]]}, y, {t, 0, 1}, "
              "WorkingPrecision -> 30]][[2, 1, 1]]]", "30.0");
    audit_chk("NIntegrate#4", "Precision[NIntegrate[Exp[x], {x, 0, 1}, WorkingPrecision -> 10]]", "10.0");
    audit_true("NIntegrate#4", "Abs[NIntegrate[Exp[x], {x, 0, 1}, WorkingPrecision -> 10] - (E - 1)] < 10^-9");
    audit_true("ND#1", "Quiet[Check[ND[Cos[x]^3, {x, 2}, 0], \"m\"]] =!= \"m\"");
    audit_true("ND#1", "Abs[ND[Cos[x]^3, {x, 2}, 0] + 3] < 10^-9");
    audit_chk("Fourier#3", "Precision[Fourier[{2.`30, 1.`10}]]", "10.0");
    audit_true("Fourier#3", "Abs[Fourier[{2.`30, 1.`10}][[1]] - 3/Sqrt[2]] < 10^-9");
    audit_true("NDSolve#1", "Module[{s = NDSolve[{y'[x] == {{0, 1}, {-1, 0}}.y[x], y[0] == {1, 0}}, y, {x, 0, 1}]}, "
               "Norm[(y[1] /. s[[1]]) - {Cos[1], -Sin[1]}] < 10^-6]");
    audit_true("NDSolve#1", "Norm[NDSolveValue[{y'[x] == {{0, 1}, {-1, 0}}.y[x], y[0] == {1, 0}}, y[1], {x, 0, 1}] - "
               "{Cos[1], -Sin[1]}] < 10^-6");
    audit_true("NSeries#4", "Abs[Coefficient[Normal[NSeries[Exp[x], {x, 0, 2}, PrecisionGoal -> Infinity, "
               "AccuracyGoal -> Infinity]], x, -2] - 1/720] < 10^-6");
}

int main(void) {
    audit_init();
    TEST(test_nsum_nested_no_d_ivar);
    TEST(test_nseries_huge_order);
    TEST(test_nseries_order_narrowing);
    TEST(test_nintegrate_pv_nodes);
    TEST(test_nsum_doubly_infinite);
    TEST(test_nminimize_box_region);
    TEST(test_nminimize_integer_edge);
    TEST(test_nmaximize_nonnumeric);
    TEST(test_findmaximum_gradient);
    TEST(test_nproduct_zero_factor);
    TEST(test_nroots_small_roots);
    TEST(test_nroots_nonnumeric_constant);
    TEST(test_nsolve_or);
    TEST(test_whenevent_jump);
    TEST(test_wave3_arity_options);
    TEST(test_wave4_numerical);
    TEST(test_wave5_g5);
    TEST(test_wave6_g2);
    return audit_finish("audit_numerical_calculus");
}
