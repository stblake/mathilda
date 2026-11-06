/* Tests for the special functions behind extreme zeros, and FindRoot on them:
 *
 *   JacobiP              exact polynomials (checked against the ODE, the
 *                        reflection P_n^(a,b)(-x) = (-1)^n P_n^(b,a)(x), the
 *                        Legendre / Chebyshev specialisations and the sum
 *                        form), machine / 40-digit / complex / non-integer
 *                        degree values against mpmath 1.4 (an independent
 *                        implementation), D, surfaces, attributes;
 *   RiemannSiegelZ /     values against mpmath (siegelz, siegeltheta), the
 *   RiemannSiegelTheta   defining identities theta = Im LogGamma(1/4 + i t/2)
 *                        - t log(pi)/2 and Z = Exp(i theta) Zeta(1/2 + i t),
 *                        symmetry, Gram points and published zeta zeros;
 *   LegendreP            an arbitrary-precision x evaluates by the recurrence
 *                        (it used to substitute into the expanded polynomial:
 *                        P_1000(0.999997`50) came out as 2.4e329);
 *   ChebyshevT           high degree through cos(n acos x) / cosh(n acosh x);
 *   FindRoot             the residual test certifies Newton's answer (a start
 *                        within the step tolerance used to return after one
 *                        step), StepMonitor sees the new iterate, and the
 *                        extreme zeros of the ISSAC '97 holdout family 5.
 *
 * Reference values: mpmath 1.4.1 at 50 digits; closed forms where they exist
 * (1 - Cos[Pi/4002] for the largest zero of ChebyshevT[2001, x]).
 */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "test_utils.h"

#include <stdio.h>
#include <stdlib.h>

/* An expression that must evaluate to True. */
static void assert_true(const char* input) { assert_eval_eq(input, "True", 0); }

/* ---- JacobiP ------------------------------------------------------------ */

void test_jacobi_exact() {
    assert_eval_eq("JacobiP[0, a, b, x]", "1", 0);
    assert_eval_eq("JacobiP[n, a, b, 1]", "Binomial[a + n, n]", 0);
    assert_eval_eq("JacobiP[5, 2, 3, 1]", "21", 0);              /* Binomial[7, 5] */
    assert_eval_eq("JacobiP[n, a, b, x]", "JacobiP[n, a, b, x]", 0);
    /* P_1 = (a + 1) + (a + b + 2)(x - 1)/2, the recurrence seed. */
    assert_true("Expand[JacobiP[1, a, b, x] - ((a + 1) + (a + b + 2) (x - 1)/2)] === 0");
    /* Legendre and Chebyshev specialisations (exact). */
    assert_true("Expand[JacobiP[9, 0, 0, x] - LegendreP[9, x]] === 0");
    assert_true("Expand[JacobiP[7, -1/2, -1/2, x] - Binomial[7 - 1/2, 7] ChebyshevT[7, x]] === 0");
    /* The ODE (1-x^2) y'' + (b - a - (a+b+2) x) y' + n(n+a+b+1) y = 0. */
    assert_true("With[{y = JacobiP[8, 1/2, 1/3, x]}, Expand[(1 - x^2) D[y, x, x]"
                " + (1/3 - 1/2 - (1/2 + 1/3 + 2) x) D[y, x] + 8 (8 + 1/2 + 1/3 + 1) y]] === 0");
    assert_true("With[{y = JacobiP[4, a, b, x]}, Expand[(1 - x^2) D[y, x, x]"
                " + (b - a - (a + b + 2) x) D[y, x] + 4 (4 + a + b + 1) y]] === 0");
    /* Reflection P_n^(a,b)(-x) = (-1)^n P_n^(b,a)(x). */
    assert_true("Expand[JacobiP[11, 2/3, -1/4, -x] + JacobiP[11, -1/4, 2/3, x]] === 0");
    /* Exact rational x by Horner == the polynomial evaluated at it. */
    assert_true("JacobiP[30, 3, 1/2, 7/9] === (JacobiP[30, 3, 1/2, x] /. x -> 7/9)");
    /* Negative integer a inside [-n, -1]: still a polynomial (the sum form has
     * no division by a parameter); it satisfies the ODE and P_n(1) = 0. */
    assert_true("With[{y = JacobiP[3, -2, 1/2, x]}, Expand[(1 - x^2) D[y, x, x]"
                " + (1/2 + 2 - (-2 + 1/2 + 2) x) D[y, x] + 3 (3 - 2 + 1/2 + 1) y] === 0 &&"
                " (y /. x -> 1) === 0 && Exponent[y, x] === 3]");
}

void test_jacobi_numeric() {
    /* Machine: the shared recurrence helper, against mpmath. */
    /* Near x = 1 the value moves by |P'/P| eps ~ 1.2e-12 relative when the
     * argument is rounded to a double, so the mpmath value at 0.9998 itself is
     * matched to that, and the value at the double-rounded point to 1e-13. */
    assert_true("Abs[JacobiP[300, 3, 1/2, 0.9998]/"
                "100955.52971661217588429811261857849206544 - 1] < 3*10^-12");
    assert_true("Abs[JacobiP[300, 3, 1/2, 0.9998]/"
                "N[JacobiP[300, 3, 1/2, SetPrecision[0.9998, 40]], 20] - 1] < 10^-13");
    assert_true("Abs[JacobiP[7, 1/3, 2.5, -0.6]/1.6520729335975207031448458060255042422 - 1] < 10^-14");
    assert_true("Abs[JacobiP[40, -0.25, 1.75, 1.3]/1696960108695.6685315698750881346 - 1] < 10^-13");
    assert_true("Abs[JacobiP[1000, 0.5, -0.5, 0.25]/(-0.018627855071485999427926374140588914798) - 1]"
                " < 10^-11");
    /* Machine complex argument. */
    assert_true("Abs[JacobiP[4, 1., 2, 0.3 + 0.2 I] - (0.6345625 - 0.9475 I)] < 10^-14");
    /* Non-integer degree: Binomial[n+a, n] 2F1(-n, n+a+b+1; a+1; (1-x)/2). */
    assert_true("Abs[JacobiP[2.5, 1, 2, 0.3]/(-0.81484244219558440084884398824695282777) - 1] < 10^-13");
    /* Arbitrary precision: the MPFR recurrence against mpmath and against the
     * exact rational evaluation (the sum form). */
    assert_true("Abs[JacobiP[300, 3, 1/2, N[9998/10000, 40]] - "
                "N[JacobiP[300, 3, 1/2, 9998/10000], 40]] < 10^-30");
    assert_true("Precision[JacobiP[300, 3, 1/2, N[9998/10000, 40]]] >= 39");
    assert_true("Abs[JacobiP[1000, 1/2, -1/2, N[1/4, 50]] + "
                "N[18627855071485999427926374140588914798267428947087/10^51, 50]] < 10^-45");
    assert_true("Abs[JacobiP[4, 1, 2, N[3/10 + I/5, 30]] - (6345625/10^7 - 9475/10^4 I)] < 10^-25");
    /* The machine scalar equals the explicit sum at the double-rounded point. */
    assert_true("Abs[JacobiP[25, 2, 3, 0.37] - N[JacobiP[25, 2, 3, x] /. x -> SetPrecision[0.37, 40], 20]]"
                " < 10^-12 Abs[JacobiP[25, 2, 3, 0.37]]");
}

void test_jacobi_derivative() {
    assert_eval_eq("D[JacobiP[n, a, b, x], x]",
                   "1/2 (1 + a + b + n) JacobiP[-1 + n, 1 + a, 1 + b, x]", 0);
    assert_true("Expand[D[JacobiP[6, 2, 1, x], x] - (6 + 2 + 1 + 1)/2 JacobiP[5, 3, 2, x]] === 0");
}

void test_jacobi_surfaces() {
    assert_true("Module[{xs = {0.1, 0.5, 0.9, -0.3, 0.99, 1.4}, sc, cf, cfa},"
                " sc = JacobiP[7, 1/3, 2.5, #] & /@ xs;"
                " cf = Compile[{{x, _Real}}, JacobiP[7, 1/3, 2.5, x]];"
                " cfa = Compile[{{v, _Real, 1}}, JacobiP[7, 1/3, 2.5, v]];"
                " sc === JacobiP[7, 1/3, 2.5, xs] &&"
                " sc === JacobiP[7, 1/3, 2.5, Developer`ToPackedArray[xs]] &&"
                " sc === Normal[JacobiP[7, 1/3, 2.5, NDArray[xs]]] &&"
                " sc === (cf /@ xs) && sc === cfa[xs]]");
    assert_true("Lookup[CompileDiagnostics[{{x, _Real}}, JacobiP[7, 1, 2, x]], \"Compiled\"] &&"
                " Lookup[CompileDiagnostics[{{v, _Real, 1}}, JacobiP[7, 1, 2, v]], \"Compiled\"]");
    assert_true("Developer`PackedArrayQ[JacobiP[3, 1, 2, Developer`ToPackedArray[N[Range[50]/51]]]]");
    assert_eval_eq("Attributes[JacobiP]", "{Listable, NumericFunction, Protected}", 0);
}

/* ---- RiemannSiegelZ / RiemannSiegelTheta ------------------------------ */

void test_riemann_siegel_values() {
    assert_eval_eq("RiemannSiegelTheta[0]", "0", 0);
    assert_eval_eq("RiemannSiegelZ[0]", "Zeta[1/2]", 0);
    assert_eval_eq("RiemannSiegelZ[t]", "RiemannSiegelZ[t]", 0);
    assert_true("Abs[RiemannSiegelZ[20.5] - 0.59932870251475132138720989122873849254] < 10^-14");
    assert_true("Abs[RiemannSiegelZ[1.] + 0.73630546286731773467789982892561467210] < 10^-15");
    assert_true("Abs[RiemannSiegelZ[1000.] - 0.99779463752158661398600268518815709241] < 10^-12");
    assert_true("Abs[RiemannSiegelTheta[10.] + 3.0670743962898952917020135348094859760] < 10^-14");
    assert_true("Abs[RiemannSiegelTheta[1000.5]/2035.8139600703492805571840086616883933 - 1] < 10^-15");
    assert_true("Abs[RiemannSiegelZ[2. + I] - (-0.47014436803007441718693502996321265359"
                " + 0.034599242557948909654678220608299506214 I)] < 10^-14");
    assert_true("Abs[RiemannSiegelTheta[3. + 0.5 I] - (-3.0155383774771599347947547847829829130"
                " - 0.18485654419418169176091552782106416002 I)] < 10^-14");
    /* 30 digits, against mpmath. */
    assert_true("Abs[RiemannSiegelZ[N[41/2, 30]] - "
                "N[59932870251475132138720989122873849254017262464966/10^50, 30]] < 10^-28");
    assert_true("Precision[RiemannSiegelZ[N[41/2, 30]]] >= 29");
    /* A real argument gives a real number. */
    assert_true("Head[RiemannSiegelZ[N[123/4, 25]]] === Real && Head[RiemannSiegelZ[30.75]] === Real");
}

void test_riemann_siegel_identities() {
    /* theta(t) = Im LogGamma(1/4 + i t/2) - t Log(Pi)/2 and
     * Z(t) = Exp(i theta) Zeta(1/2 + i t), at 30 digits, an independent route. */
    assert_true("With[{t = N[173/10, 30]}, Abs[RiemannSiegelTheta[t] -"
                " (Im[LogGamma[1/4 + I t/2]] - t Log[Pi]/2)] < 10^-26]");
    assert_true("With[{t = N[173/10, 30]}, Abs[RiemannSiegelZ[t] -"
                " Re[Exp[I RiemannSiegelTheta[t]] Zeta[1/2 + I t]]] < 10^-26]");
    assert_true("Abs[Abs[RiemannSiegelZ[30.3]] - Abs[Zeta[1/2 + 30.3 I]]] < 10^-14");
    /* theta odd, Z even. */
    assert_true("RiemannSiegelTheta[-7.25] == -RiemannSiegelTheta[7.25]");
    assert_true("Abs[RiemannSiegelZ[-7.25] - RiemannSiegelZ[7.25]] < 10^-15");
    /* Gram points: theta(g_n) = n Pi (g_0, g_1 published). */
    assert_true("Abs[RiemannSiegelTheta[N[17845599540410860816826338412519097/10^33, 35]]] < 10^-30");
    assert_true("Abs[RiemannSiegelTheta[23.170282701246309279] - Pi] < 10^-14");
}

void test_riemann_siegel_zeros() {
    assert_true("(zF[t_?NumericQ] := RiemannSiegelZ[t]; True)");
    /* First and tenth zeros of Zeta on the critical line (published). */
    assert_true("Abs[(t /. FindRoot[zF[t], {t, 14.1}]) - 14.134725141734693790] < 10^-12");
    assert_true("Abs[(t /. FindRoot[zF[t], {t, 4977/100}]) - 49.773832477672302182] < 10^-12");
    assert_true("Abs[(t /. FindRoot[zF[t], {t, 4977/100}, WorkingPrecision -> 50]) - "
                "49773832477672302181916784678563724057723178299677/10^48] < 10^-40");
    /* The packed / Compile surfaces agree with the scalar. */
    assert_true("Module[{ts = {1., 14.134725141734693, 20.5, 49.7, 100.3}, sc, cf, cfa},"
                " sc = RiemannSiegelZ /@ ts;"
                " cf = Compile[{{x, _Real}}, RiemannSiegelZ[x]];"
                " cfa = Compile[{{v, _Real, 1}}, RiemannSiegelZ[v]];"
                " sc === RiemannSiegelZ[ts] && sc === RiemannSiegelZ[Developer`ToPackedArray[ts]] &&"
                " sc === Normal[RiemannSiegelZ[NDArray[ts]]] && sc === (cf /@ ts) && sc === cfa[ts] &&"
                " (RiemannSiegelTheta /@ ts) === RiemannSiegelTheta[Developer`ToPackedArray[ts]]]");
    assert_true("Lookup[CompileDiagnostics[{{x, _Real}}, RiemannSiegelTheta[x]], \"Compiled\"] &&"
                " Lookup[CompileDiagnostics[{{v, _Real, 1}}, RiemannSiegelZ[v]], \"Compiled\"]");
    assert_eval_eq("Attributes[RiemannSiegelZ]", "{Listable, NumericFunction, Protected}", 0);
    assert_eval_eq("Attributes[RiemannSiegelTheta]", "{Listable, NumericFunction, Protected}", 0);
}

/* ---- LegendreP / ChebyshevT at high degree ----------------------------- */

void test_legendre_mpfr_recurrence() {
    /* mpmath legendre(1000, 0.999997) = -0.0235907712313252317571870076602339702585... */
    assert_true("Abs[LegendreP[1000, N[999997/1000000, 50]] + "
                "N[23590771231325231757187007660233970258521559414788/10^51, 50]] < 10^-45");
    assert_true("Precision[LegendreP[1000, N[999997/1000000, 50]]] >= 49");
    assert_true("Abs[LegendreP[50, N[3/10 + I/5, 30]] - N[LegendreP[50, 3/10 + I/5], 30]] < 10^-24");
    assert_true("LegendreP[1, N[1/3, 30]] == N[1/3, 30] && LegendreP[0, N[1/3, 30]] == 1");
    /* Bonnet's recurrence holds between the arbitrary-precision values. */
    assert_true("With[{x = N[7/11, 40]}, Abs[31 LegendreP[31, x] - 61 x LegendreP[30, x]"
                " + 30 LegendreP[29, x]] < 10^-35]");
}

void test_chebyshev_high_degree() {
    /* Machine, degree > 1000: cos(n acos x), within a few ulps of the
     * conditioning |T_n'(x)| eps of the double argument. */
    assert_true("Abs[ChebyshevT[2001, 0.9999997] - "
                "N[ChebyshevT[2001, SetPrecision[0.9999997, 40]], 20]] < 10^-12");
    assert_true("Abs[ChebyshevT[2001, 0.9999997] - 0.020826846965824051968644242578695765] < 10^-9");
    assert_true("Abs[ChebyshevT[2001, 1.0001]/974344123482.34101911098025616821 - 1] < 10^-11");
    assert_true("Abs[ChebyshevT[2001, -1.0001]/974344123482.34101911098025616821 + 1] < 10^-11");
    assert_true("Abs[ChebyshevT[5000, -0.7] - 0.96583540433174708987979841278901005735] < 10^-11");
    /* Arbitrary precision, degree > 64: the same closed form in MPFR. */
    assert_true("Abs[ChebyshevT[1500, N[3/7, 40]] + "
                "N[80193918694834440763482624174775447272660692980456/10^51, 40]] < 10^-34");
    assert_true("Abs[ChebyshevT[2001, N[10001/10000, 40]]/N[ChebyshevT[2001, 10001/10000], 40] - 1]"
                " < 10^-34");
    assert_true("ChebyshevT[2001, N[-10001/10000, 30]] == -ChebyshevT[2001, N[10001/10000, 30]]");
    /* Every surface agrees bit for bit at high degree. */
    assert_true("Module[{cft = Compile[{{n, _Integer}, {x, _Real}}, ChebyshevT[n, x]]},"
                " ChebyshevT[2001, 0.9999997] === cft[2001, 0.9999997] &&"
                " (ChebyshevT[2001, #] & /@ {0.3, -0.99, 1.2}) === "
                " ChebyshevT[2001, Developer`ToPackedArray[{0.3, -0.99, 1.2}]]]");
}

/* ---- FindRoot --------------------------------------------------------- */

void test_findroot_residual_certifies() {
    /* A start within the step tolerance of a root at 1 - 3.08e-7: the first
     * step used to pass the step test with 1 - x right to 3 digits.  Closed
     * form 1 - Cos[Pi/4002] = 3.0811692773857049463438593532e-7. */
    assert_true("(tF[x_?NumericQ] := ChebyshevT[2001, x]; True)");
    assert_true("Abs[(1 - (x /. FindRoot[tF[x], {x, 0.9999997}]))/(1 - Cos[Pi/4002]) - 1] < 10^-6");
    assert_true("Abs[(1 - (x /. FindRoot[tF[x], {x, 9999997/10000000}, WorkingPrecision -> 50]))"
                " - (1 - Cos[Pi/4002])] < 10^-40");
    /* Same root through the analytic derivative (no finite difference). */
    assert_true("Abs[(1 - (x /. FindRoot[ChebyshevT[2001, x], {x, 0.9999997}]))/(1 - Cos[Pi/4002]) - 1]"
                " < 10^-6");
    /* Ordinary roots are unchanged: Sin at Pi, a cubic, a complex root. */
    assert_true("Abs[(x /. FindRoot[Sin[x], {x, 3.}]) - N[Pi]] < 10^-15");
    assert_true("Abs[(x /. FindRoot[x^3 - 2, {x, 1.}]) - 2^(1/3)] < 10^-15");
    assert_true("Abs[(x /. FindRoot[x^2 + 1, {x, 0.5 + 0.5 I}]) - I] < 10^-15");
    /* A residual that can never reach the goal (|f| ~ 1e20 scale): the
     * polish stops in the noise and returns the root. */
    assert_true("Abs[(x /. FindRoot[10^20 (x - 1/3) + 10^20 (x - 1/3)^2, {x, 0.3}]) - 1/3] < 10^-15");
}

void test_findroot_step_monitor_iterate() {
    /* The StepMonitor sees each new iterate, not a finite-difference probe:
     * every reported point is a Newton iterate of h, so h shrinks along them. */
    assert_true("(hh[x_?NumericQ] := Sin[x] - 1/2; True)");
    assert_true("Module[{pts = Reap[FindRoot[hh[x], {x, 0.5}, StepMonitor :> Sow[x]]][[2, 1]]},"
                " Abs[Last[pts] - Pi/6] < 10^-15 && And @@ Thread[Abs[Rest[hh /@ pts]] <"
                " Abs[Most[hh /@ pts]] + 10^-16]]");
}

void test_findroot_extreme_zeros() {
    /* ISSAC '97 holdout family 5 (benchmarks/140): 1 - x at the largest zero. */
    assert_true("(pF[x_?NumericQ] := LegendreP[1000, x]; jF[x_?NumericQ] := JacobiP[300, 3, 1/2, x]; True)");
    assert_true("Abs[(1 - (x /. FindRoot[pF[x], {x, 999997/1000000}, WorkingPrecision -> 50,"
                " MaxIterations -> 200])) - 2888701924489430123709748/10^30] < 10^-26");
    assert_true("Abs[(1 - (x /. FindRoot[jF[x], {x, 9998/10000}]))/2.227861165971667686640295*^-4 - 1]"
                " < 10^-6");
    assert_true("Abs[(1 - (x /. FindRoot[jF[x], {x, 9998/10000}, WorkingPrecision -> 50,"
                " MaxIterations -> 200])) - 2227861165971667686640295/10^28] < 10^-24");
}

int main(void) {
    symtab_init();
    core_init();
    test_mute_messages();

    TEST(test_jacobi_exact);
    TEST(test_jacobi_numeric);
    TEST(test_jacobi_derivative);
    TEST(test_jacobi_surfaces);
    TEST(test_riemann_siegel_values);
    TEST(test_riemann_siegel_identities);
    TEST(test_riemann_siegel_zeros);
    TEST(test_legendre_mpfr_recurrence);
    TEST(test_chebyshev_high_degree);
    TEST(test_findroot_residual_certifies);
    TEST(test_findroot_step_monitor_iterate);
    TEST(test_findroot_extreme_zeros);

    printf("All special-roots tests passed.\n");
    return 0;
}
