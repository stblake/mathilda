(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 153 -- ISSAC '97 holdout B, family 8: definite integrals with
   Zeta / PolyLog / Log / Catalan / EulerGamma closed forms (the counterpart
   of Problem 8, Integrate[x^2 PolyLog[3, 1/(x + 1)], {x, 0, 1}])
   ==========================================================================

       cd benchmarks/153-issac97b-closed-form-integrals
       ../../Mathilda -file closed_form_integrals.m

   The original is a polylogarithm of a Moebius-transformed argument times a
   polynomial.  These integrands are logarithms of algebraic arguments
   (Serret's integral, a cyclotomic quadratic), products of logarithms and
   polylogarithms at weight 4, a trigonometric integral and an arctangent
   square whose answers carry Catalan's constant, a Gaussian times a
   logarithm (EulerGamma), an integrand on [0, oo) with a hyperbolic
   denominator, and LogGamma against sin(2 pi x) (Kummer's Fourier series).
   Every closed form was derived by hand -- the derivation is noted -- and
   confirmed by mpmath tanh-sinh quadrature at 40 digits (30 digits agree;
   Gauss-Legendre agrees to 8+ on the four with smooth integrands):

   b8.1  int_0^1 log(1+x)/(1+x^2) = pi log2/8 = 0.2721982612879502663125861
         (Serret: x = tan t, then t -> pi/4 - t)
   b8.2  int_0^1 log^2 x log(1-x)/x = -2 zeta(4) = -pi^4/45
         = -2.164646467422276383032007   (expand log(1-x), int x^(k-1) log^2 x = 2/k^3)
   b8.3  int_0^{pi/2} x/sin x = 2 G = 1.831931188354438030109207
         (x = 2 arctan u: the integral becomes int_0^1 2 arctan(u)/u du = 2 G)
   b8.4  int_0^oo log(x) e^(-x^2) = -sqrt(pi)(gamma + 2 log 2)/4
         = -0.8700577267283155067346488  (d/ds of Gamma(s/2)/2 at s = 1)
   b8.5  int_0^1 Li_2(x) log(1-x)/x = -Li_2(1)^2/2 = -pi^4/72
         = -1.352904042138922739395005   (Li_2' = -log(1-x)/x)
   b8.6  int_0^1 log(1 - x + x^2)/x = -pi^2/18 = -0.5483113556160754788241384
         (1 - x + x^2 = (1 + x^3)/(1 + x); int log(1+x^k)/x = pi^2/(12 k))
   b8.7  int_0^1 log x log(1+x)/x = -eta(3) = -3 zeta(3)/4
         = -0.9015426773696957140498036
   b8.8  int_0^1 arctan(x)^2/x = pi G/2 - 7 zeta(3)/8
         = 0.3869956005394355761640192   (parts, then x = tan t: a log-sine integral)
   b8.9  int_0^oo x/sinh x = 2 sum_k 1/(2k+1)^2 = pi^2/4 = 2.467401100272339654708623
   b8.10 int_0^1 LogGamma(x) sin(2 pi x) = (gamma + log 2pi)/(2 pi)
         = 0.3843739462134329156026643   (Kummer's Fourier series of log Gamma)

   METHODS.  Integrate (exact), NIntegrate at machine precision and at
   WorkingPrecision -> 40 (default PrecisionGoal, 20 digits).
   CHECKS.  Exact: sig[N[v, 30], 15] (Chop at 1e-20: a closed form evaluated
   through complex intermediates leaves a 0``15 I).  Machine: 8 digits (6 for
   the three with a log-power endpoint singularity, b8.2, b8.4, b8.7).
   WP 40: 18 digits (16 for b8.9, whose 18th digit sits on a rounding edge).
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"Integrate", "NIntegrate", "PolyLog", "LogGamma", "Catalan", "Zeta"}];

exactChk[v_] := If[FreeQ[v, Integrate] && NumericQ[v], sig[Chop[N[v, 30], 10^-20], 15], "unevaluated"];
mc[d_][v_] := If[NumericQ[v], sig[v, d], "unevaluated"];

SetAttributes[three, HoldAll];
three[tag_, f_, iv_, dm_, dw_:18] := (
  run[tag <> " Integrate exact", Integrate[f, iv], exactChk];
  run[tag <> " NIntegrate wp40", NIntegrate[f, iv, WorkingPrecision -> 40], mc[dw], 3];
  run[tag <> " NIntegrate machine", NIntegrate[f, iv], mc[dm], 3]);

three["b8.1", Log[1 + x]/(1 + x^2), {x, 0, 1}, 8];
three["b8.2", Log[x]^2 Log[1 - x]/x, {x, 0, 1}, 6];
three["b8.3", x/Sin[x], {x, 0, Pi/2}, 8];
three["b8.4", Log[x] Exp[-x^2], {x, 0, Infinity}, 6];
three["b8.5", PolyLog[2, x] Log[1 - x]/x, {x, 0, 1}, 8];
three["b8.6", Log[1 - x + x^2]/x, {x, 0, 1}, 8];
three["b8.7", Log[x] Log[1 + x]/x, {x, 0, 1}, 6];
three["b8.8", ArcTan[x]^2/x, {x, 0, 1}, 8];
three["b8.9", x/Sinh[x], {x, 0, Infinity}, 8, 16];
three["b8.10", LogGamma[x] Sin[2 Pi x], {x, 0, 1}, 8];
