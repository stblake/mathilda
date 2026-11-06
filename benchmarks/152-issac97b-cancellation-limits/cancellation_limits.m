(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 152 -- ISSAC '97 holdout B, family 7: limits with deep
   cancellation (the counterpart of Problem 7, the ratio of twelve-term nested
   Sin/Tan/Sinh/Tanh sums)
   ==========================================================================

       cd benchmarks/152-issac97b-cancellation-limits
       ../../Mathilda -file cancellation_limits.m

   The original cancels inside compositions of elementary functions at 0.
   These cancel in other places: inside the Taylor series of LogGamma, among
   exponentials of different growth at infinity (Gruntz), in Arnold's
   tan-sin problem, in the Laurent expansion of zeta at 1 (the answer is a
   Stieltjes constant), in an exponent, against logarithms at 0+, in the
   asymptotics of the harmonic numbers, between two Bessel compositions, in
   the asymptotic series of erfc at infinity, and in a power mean.
   True values, each derived by hand and confirmed by evaluating the
   expression at 10^-10 .. 10^-100 (or 10^10 .. 10^30) in mpmath at up to
   9000 digits (and, where noted, by a sympy series):

   b7.1  lim_{x->0} (LogGamma(1+x) + gamma x - zeta(2) x^2/2 + zeta(3) x^3/3)/x^4
         = zeta(4)/4 = pi^4/360 = 0.2705808084277845478790009
         (log Gamma(1+x) = -gamma x + sum_{k>=2} (-1)^k zeta(k) x^k/k; sympy agrees)
   b7.2  lim_{x->oo} (exp(x e^-x/(e^-x + e^(-2x^2/(x+1)))) - e^x)/x = -e^2
         = -7.389056098930650227230427 (Gruntz's thesis; numerically -7.24418,
         -7.37431, -7.38758 at x = 1e2, 1e3, 1e4, converging as 1/x)
   b7.3  lim_{x->0} (sin tan x - tan sin x)/(arcsin arctan x - arctan arcsin x) = 1
         (Arnold; numerator and denominator are both -x^7/30 + ..., sympy)
   b7.4  lim_{n->oo} n (zeta(1 + 1/n) - n - gamma) = -gamma_1
         = 0.0728158454836767248605863758749547830 (the Laurent series of zeta
         at 1; mpmath stieltjes)
   b7.5  lim_{n->oo} n (((1+1/n)^n/e)^n e^(1/2) - 1) = 1/3
         (n^2 log(1+1/n) - n + 1/2 = 1/(3n) - 1/(4n^2) + ...)
   b7.6  lim_{x->0+} (x^(x^x) - x)/(x^2 log^2 x) = 1
         (x^(x^x) = x exp((x^x - 1) log x), x^x - 1 = x log x + ...)
   b7.7  lim_{n->oo} n^2 (H_n - log n - gamma - 1/(2n)) = -1/12
   b7.8  lim_{x->0} (2 J_1(sin x) - sin(2 J_1(x)))/x^7 = -13/5760
         (both are x - 7x^3/24 + 31 x^5/384 + ...; sympy series to x^11)
   b7.9  lim_{x->oo} x^6 (sqrt(pi) x e^(x^2) erfc x - 1 + 1/(2x^2) - 3/(4x^4)) = -15/8
         (the asymptotic series sum (-1)^m (2m-1)!!/(2x^2)^m)
   b7.10 lim_{x->0} (((2^x + 3^x + 5^x)/3)^(1/x) - 30^(1/3))/x
         = 30^(1/3) Var(log 2, log 3, log 5)/2 = 0.218357912107657334199163911399
         (the cumulant expansion of (1/x) log mean e^(x L))

   METHODS.  Limit (exact); Series (b7.3, b7.8: ratio / coefficient of the
   leading terms); the numeric approach at arbitrary precision -- evaluate at
   10^-k (or 10^k) at a working precision that covers the cancellation, until
   two successive values agree to 1e-15 (b7.2 instead extrapolates 10 values
   at x = 100..550 by Neville's algorithm at 600 digits, since its error
   decays only like 1/x; b7.9 extrapolates x = 20..65 in h = 1/x^2 at 100
   digits, since e^(x^2) at the x ~ 1e8 a 1e-15 approach needs overflows any
   arbitrary-precision exponent range; mpmath: error 8e-24); and at machine precision Neville extrapolation to
   h = 0 of a few samples (b7.6: one sample at x = 1e-8, where the error is
   x log x ~ 1.6e-6).  Machine sample points are chosen so the cancellation
   leaves at least 8 digits (b7.2: 3).
   CHECKS.  Exact: sig[N[v, 30], 15].  Numeric approach and b7.9: 12 digits
   (b7.2: 10).  Machine: the digits Neville's extrapolation of IEEE doubles reaches
   on the same samples, less one (3-8; b7.2, whose samples cancel 13 digits at x = 30: 2).
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"Limit", "Series", "SeriesCoefficient", "LogGamma", "Zeta", "StieltjesGamma", "HarmonicNumber", "BesselJ", "Erfc"}];

symChk[v_] := If[NumericQ[v] && FreeQ[v, Limit | Series | SeriesData], sig[N[v, 30], 15], "unevaluated"];
mc[d_][v_] := If[NumericQ[v], sig[v, d], "unevaluated"];

(* evaluate f at the points xs (at precision wp) until two successive values agree to tol *)
approach[f_, xs_List, wp_, tol_] := Module[{old = 10^100, new = 0, k = 0},
  While[k < Length[xs] && (k = k + 1; new = f[N[xs[[k]], wp]]; Abs[new - old] > tol Abs[new]), old = new];
  new];
(* Neville's algorithm: the interpolating polynomial through (hs, vs) evaluated at h = 0 *)
neville[hs_List, vs_List] := Module[{p = vs, n = Length[hs]},
  Do[p = Table[(hs[[i]] p[[i + 1]] - hs[[i + m]] p[[i]])/(hs[[i]] - hs[[i + m]]), {i, n - m}], {m, 1, n - 1}];
  First[p]];

(* ---- b7.1 LogGamma Taylor cancellation -------------------------------- *)
f1[x_] := (LogGamma[1 + x] + EulerGamma x - Zeta[2] x^2/2 + Zeta[3] x^3/3)/x^4;
run["b7.1 Limit exact", Limit[(LogGamma[1 + x] + EulerGamma x - Zeta[2] x^2/2 + Zeta[3] x^3/3)/x^4, x -> 0], symChk];
run["b7.1 approach wp300", approach[f1, 10^-Range[4, 60, 4], 300, 10^-15], mc[12], 3];
run["b7.1 Neville machine", neville[2.^-Range[2, 7], f1 /@ (2.^-Range[2, 7])], mc[4], 3];

(* ---- b7.2 Gruntz ------------------------------------------------------ *)
f2[x_] := (Exp[x Exp[-x]/(Exp[-x] + Exp[-2 x^2/(x + 1)])] - Exp[x])/x;
run["b7.2 Limit exact", Limit[(Exp[x Exp[-x]/(Exp[-x] + Exp[-2 x^2/(x + 1)])] - Exp[x])/x, x -> Infinity], symChk];
run["b7.2 Neville wp600", Module[{xs = N[Range[100, 550, 50], 600]}, neville[1/xs, f2 /@ xs]], mc[10], 3];
run["b7.2 Neville machine", Module[{xs = {20., 25., 30.}}, neville[1/xs, f2 /@ xs]], mc[2], 3];

(* ---- b7.3 Arnold ------------------------------------------------------ *)
num3 = Sin[Tan[x]] - Tan[Sin[x]]; den3 = ArcSin[ArcTan[x]] - ArcTan[ArcSin[x]];
f3[x_] := (Sin[Tan[x]] - Tan[Sin[x]])/(ArcSin[ArcTan[x]] - ArcTan[ArcSin[x]]);
run["b7.3 Limit exact", Limit[num3/den3, x -> 0], symChk];
run["b7.3 Series exact", SeriesCoefficient[num3, {x, 0, 7}]/SeriesCoefficient[den3, {x, 0, 7}], symChk];
run["b7.3 approach wp200", approach[f3, 10^-Range[2, 30, 2], 200, 10^-15], mc[12], 3];
run["b7.3 Neville machine", neville[(2.^-Range[1, 4])^2, f3 /@ (2.^-Range[1, 4])], mc[5], 3];

(* ---- b7.4 Stieltjes constant ------------------------------------------ *)
f4[n_] := n (Zeta[1 + 1/n] - n - EulerGamma);
run["b7.4 Limit exact", Limit[n (Zeta[1 + 1/n] - n - EulerGamma), n -> Infinity], symChk];
run["b7.4 approach wp100", approach[f4, 10^Range[2, 40, 2], 100, 10^-15], mc[12], 3];
run["b7.4 Neville machine", Module[{ns = 2.^Range[2, 7]}, neville[1/ns, f4 /@ ns]], mc[5], 3];

(* ---- b7.5 cancellation in an exponent --------------------------------- *)
f5[n_] := n (((1 + 1/n)^n/E)^n Exp[1/2] - 1);
run["b7.5 Limit exact", Limit[n (((1 + 1/n)^n/E)^n Exp[1/2] - 1), n -> Infinity], symChk];
run["b7.5 approach wp150", approach[f5, 10^Range[2, 40, 2], 150, 10^-15], mc[12], 3];
run["b7.5 Neville machine", Module[{ns = 2.^Range[3, 8]}, neville[1/ns, f5 /@ ns]], mc[8], 3];

(* ---- b7.6 logarithms at 0+ -------------------------------------------- *)
f6[x_] := (x^(x^x) - x)/(x^2 Log[x]^2);
run["b7.6 Limit exact", Limit[(x^(x^x) - x)/(x^2 Log[x]^2), x -> 0, Direction -> "FromAbove"], symChk];
run["b7.6 approach wp150", approach[f6, 10^-Range[2, 60, 2], 150, 10^-15], mc[12], 3];
run["b7.6 sample machine x=1e-8", f6[1.*^-8], mc[5], 3];

(* ---- b7.7 harmonic numbers -------------------------------------------- *)
f7[n_] := n^2 (HarmonicNumber[n] - Log[n] - EulerGamma - 1/(2 n));
run["b7.7 Limit exact", Limit[n^2 (HarmonicNumber[n] - Log[n] - EulerGamma - 1/(2 n)), n -> Infinity], symChk];
run["b7.7 approach wp150", approach[f7, 10^Range[2, 40, 2], 150, 10^-15], mc[12], 3];
run["b7.7 Neville machine", Module[{ns = 2.^Range[3, 7]}, neville[1/ns^2, f7 /@ ns]], mc[7], 3];

(* ---- b7.8 Bessel compositions ----------------------------------------- *)
num8 = 2 BesselJ[1, Sin[x]] - Sin[2 BesselJ[1, x]];
f8[x_] := (2 BesselJ[1, Sin[x]] - Sin[2 BesselJ[1, x]])/x^7;
run["b7.8 Limit exact", Limit[num8/x^7, x -> 0], symChk];
run["b7.8 Series exact", SeriesCoefficient[num8, {x, 0, 7}], symChk];
run["b7.8 approach wp250", approach[f8, 10^-Range[2, 30, 2], 250, 10^-15], mc[12], 3];
run["b7.8 Neville machine", neville[(2.^-Range[2, 5])^2, f8 /@ (2.^-Range[2, 5])], mc[3], 3];

(* ---- b7.9 erfc asymptotics -------------------------------------------- *)
f9[x_] := x^6 (Sqrt[Pi] x Exp[x^2] Erfc[x] - 1 + 1/(2 x^2) - 3/(4 x^4));
run["b7.9 Limit exact", Limit[x^6 (Sqrt[Pi] x Exp[x^2] Erfc[x] - 1 + 1/(2 x^2) - 3/(4 x^4)), x -> Infinity], symChk];
run["b7.9 Neville wp100", Module[{xs = N[Range[20, 65, 5], 100]}, neville[1/xs^2, f9 /@ xs]], mc[12], 3];
run["b7.9 Neville machine", Module[{xs = {6., 7., 8., 9., 10.}}, neville[1/xs^2, f9 /@ xs]], mc[4], 3];

(* ---- b7.10 power mean ------------------------------------------------- *)
f10[x_] := (((2^x + 3^x + 5^x)/3)^(1/x) - 30^(1/3))/x;
run["b7.10 Limit exact", Limit[(((2^x + 3^x + 5^x)/3)^(1/x) - 30^(1/3))/x, x -> 0], symChk];
run["b7.10 approach wp100", approach[f10, 10^-Range[2, 40, 2], 100, 10^-15], mc[12], 3];
run["b7.10 Neville machine", neville[2.^-Range[4, 9], f10 /@ (2.^-Range[4, 9])], mc[7], 3];
