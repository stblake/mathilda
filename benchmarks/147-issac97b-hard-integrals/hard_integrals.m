(* bench-timeout: 1200 *)
(* ==========================================================================
   Experiment 147 -- ISSAC '97 holdout B, family 2: hard definite integrals
   (the counterpart of Problem 2, Integrate[x^x^x, {x, 1, 6}])
   ==========================================================================

       cd benchmarks/147-issac97b-hard-integrals
       ../../Mathilda -file hard_integrals.m

   Ten numerically hostile integrals, each hostile in a different way from
   the original's boundary layer and out-of-range answer: infinite
   oscillation at an endpoint, a near-non-integrable power singularity, a
   1e-5-wide interior peak, oscillation accelerating to infinity, a Gaussian
   boundary layer whose answer overflows a double, cancellation in the
   integrand, a 3-D integral with a point singularity, an unannounced
   interior log singularity at an irrational point, 3000 oscillations on
   [0, 1], and a slowly decaying oscillatory Bessel product on [0, oo).
   Every closed form below was derived by hand (the derivation is given) and
   confirmed by an independent quadrature in mpmath (tanh-sinh with the
   singularity split out, or quadosc), to at least 18 digits.

   b2.1  int_0^1 sin(1/x) dx = sin 1 - Ci(1)            (t = 1/x, then parts)
         = 0.504067061906928371989856
   b2.2  int_0^1 x^(-999/1000) cos x dx = sum_k (-1)^k/((2k)! (2k + 1/1000))
         = 999.760310630203964532180          (check: x = u^1000 makes it smooth)
   b2.3  int_{-1}^{1} dx/(10^-10 + x^2) = 2 10^5 ArcTan[10^5]
         = 314157.265358979390512931
   b2.4  int_0^oo sin(x^3) dx = Gamma(4/3) sin(pi/6) = Gamma(4/3)/2
         = 0.446489755784624605609282      (u = x^3: a Mellin transform of sin)
   b2.5  int_0^X exp(x^2) dx = (sqrt(pi)/2) erfi(X): X = 50 gives
         5.448684490011405885572774e1083 (outside the double range, so the
         machine case must escalate), X = 25 (the machine-sized instance)
         gives 5.43787622880010660341099e269.  Second method: the asymptotic
         series e^(X^2)/(2X) sum (2k-1)!!/(2X^2)^k.
   b2.6  int_0^1 (e^x - 1 - x)/x^2 dx = sum_{k>=2} 1/(k! (k-1))
         = 0.599620322995358659499721  (termwise; cancellation near x = 0)
   b2.7  Watson's simple-cubic integral
         int_{[0,pi]^3} dx dy dz/(3 - cos x - cos y - cos z) = pi^3 W/3,
         W = sqrt(6)/(32 pi^3) Gamma(1/24) Gamma(5/24) Gamma(7/24) Gamma(11/24)
         (Glasser-Zucker), = 15.67249523473857324485605.  Second method:
         pi^3 int_0^oo e^(-3t) I_0(t)^3 dt.
   b2.8  int_0^4 log|x^2 - 2| dx, singular at sqrt 2 (not a breakpoint):
         F(x) = x log|x^2-2| - 2x + sqrt2 log|(x+sqrt2)/(x-sqrt2)|, F(4) - F(0)
         = 3.601330233221994217662537
   b2.9  int_0^1 cos(10^4 x^2) dx = sqrt(pi/(2w)) C(sqrt(2w/pi)), w = 10^4
         (C the Fresnel integral) = 0.006251292347636025417761906
   b2.10 int_0^oo J_0(x) cos(x/2) dx = 1/sqrt(1 - 1/4) = 2/sqrt 3
         (int_0^oo J_0(a x) cos(b x) dx = 1/sqrt(a^2 - b^2), b < a)
         = 1.154700538379251529018298

   METHODS.  NIntegrate with default options at machine precision, and with
   WorkingPrecision -> 30 (default PrecisionGoal, i.e. 15 digits); b2.7 at
   WorkingPrecision -> 20, PrecisionGoal -> 8 (a 3-D singular integral).
   CHECKS.  Machine 4-7 digits (default PrecisionGoal is 6-10 digits; 4 for the
   3-D singular one); 12 at 30 digits; 6 for b2.7 at 20 digits.
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"NIntegrate", "BesselJ"}];

mc[d_][v_] := If[NumericQ[v], sig[v, d], "unevaluated"];

run["b2.1 NIntegrate machine", NIntegrate[Sin[1/x], {x, 0, 1}], mc[5], 3];
run["b2.1 NIntegrate wp30", NIntegrate[Sin[1/x], {x, 0, 1}, WorkingPrecision -> 30], mc[12]];

run["b2.2 NIntegrate machine", NIntegrate[x^(-999/1000) Cos[x], {x, 0, 1}], mc[6], 3];
run["b2.2 NIntegrate wp30", NIntegrate[x^(-999/1000) Cos[x], {x, 0, 1}, WorkingPrecision -> 30], mc[12]];

run["b2.3 NIntegrate machine", NIntegrate[1/(10^-10 + x^2), {x, -1, 1}], mc[6], 3];
run["b2.3 NIntegrate wp30", NIntegrate[1/(10^-10 + x^2), {x, -1, 1}, WorkingPrecision -> 30], mc[12]];

run["b2.4 NIntegrate machine", NIntegrate[Sin[x^3], {x, 0, Infinity}], mc[5], 3];
run["b2.4 NIntegrate wp30", NIntegrate[Sin[x^3], {x, 0, Infinity}, WorkingPrecision -> 30], mc[12]];

run["b2.5 NIntegrate machine", NIntegrate[Exp[x^2], {x, 0, 50}], mc[6], 3];
run["b2.5 NIntegrate machine X=25", NIntegrate[Exp[x^2], {x, 0, 25}], mc[6], 3];
run["b2.5 NIntegrate wp30", NIntegrate[Exp[x^2], {x, 0, 50}, WorkingPrecision -> 30], mc[12]];

run["b2.6 NIntegrate machine", NIntegrate[(Exp[x] - 1 - x)/x^2, {x, 0, 1}], mc[7], 3];
run["b2.6 NIntegrate wp30", NIntegrate[(Exp[x] - 1 - x)/x^2, {x, 0, 1}, WorkingPrecision -> 30], mc[12]];

watson = 1/(3 - Cos[x] - Cos[y] - Cos[z]);
run["b2.7 NIntegrate machine", NIntegrate[watson, {x, 0, Pi}, {y, 0, Pi}, {z, 0, Pi}], mc[4]];
run["b2.7 NIntegrate wp20",
    NIntegrate[watson, {x, 0, Pi}, {y, 0, Pi}, {z, 0, Pi}, WorkingPrecision -> 20, PrecisionGoal -> 8], mc[6]];

run["b2.8 NIntegrate machine", NIntegrate[Log[Abs[x^2 - 2]], {x, 0, 4}], mc[6], 3];
run["b2.8 NIntegrate wp30", NIntegrate[Log[Abs[x^2 - 2]], {x, 0, 4}, WorkingPrecision -> 30], mc[12]];

run["b2.9 NIntegrate machine", NIntegrate[Cos[10^4 x^2], {x, 0, 1}], mc[5], 3];
run["b2.9 NIntegrate wp30", NIntegrate[Cos[10^4 x^2], {x, 0, 1}, WorkingPrecision -> 30], mc[12]];

run["b2.10 NIntegrate machine", NIntegrate[BesselJ[0, x] Cos[x/2], {x, 0, Infinity}], mc[5], 3];
run["b2.10 NIntegrate wp30", NIntegrate[BesselJ[0, x] Cos[x/2], {x, 0, Infinity}, WorkingPrecision -> 30], mc[12]];
