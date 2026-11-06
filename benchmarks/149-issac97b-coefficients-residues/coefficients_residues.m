(* bench-timeout: 2400 *)
(* ==========================================================================
   Experiment 149 -- ISSAC '97 holdout B, family 4: exact coefficients of
   huge series and residues (the counterpart of Problem 4, the middle
   coefficient of a product of cyclotomic powers)
   ==========================================================================

       cd benchmarks/149-issac97b-coefficients-residues
       ../../Mathilda -file coefficients_residues.m

   The original is a polynomial (a finite product of fixed powers).  These are
   coefficients of an infinite product, a transcendental function, an
   algebraic function, an exponential of an exponential, an orthogonal
   polynomial, and a rational function with poles on the unit circle; and
   residues at a triple pole of a special function, at an essential
   singularity, and at a 41st-order pole of a trigonometric function.  Each
   has an exact case (SeriesCoefficient / Coefficient / Residue), an
   arbitrary-precision case (the Cauchy integral on a circle, trapezoidal
   rule, at a working precision chosen from the digits the contour loses) and
   a machine case (the Cauchy integral at machine precision; a stated reduced
   instance where the full one loses more than a double holds).
   True answers (pure-integer Python unless noted; each with a second method):

   b4.1  [x^1000] prod_{k>=1} 1/(1-x^k) = p(1000)
         = 24061467864032622473692149727991     (pentagonal recurrence;
         second: Rademacher's series).  Machine instance p(100) = 190569292.
         Contour |x| = 9/10 loses 21 digits at n = 1000, 2 at n = 100.
   b4.2  [x^400] sec x = |E_400|/400!  = 4.53899084277565898...e-79, an exact
         rational (zigzag/Entringer integer recurrence; second: mpmath
         eulernum).  Machine instance [x^20] = 370371188237525/2432902008176640000
         = 1.52234322217977e-4.  Contour |x| = 1 loses 79 / 4 digits.
   b4.3  Res_{z=-20} Gamma(z)^3 = (1/20!^3) [ (9/2)(H_20 - gamma)^2
         + pi^2/4 + (3/2) H_20^(2) ] = 3.18864502745108196186819681843e-54
         (derived from Gamma(-n+e) = (-1)^n pi/(sin(pi e) Gamma(n+1-e));
         second: the trapezoidal contour |z + 20| = 1/2 at 80 digits, 30 digits
         agree).
   b4.4  [x^10000] 1/((1-x)(1-x^5)(1-x^10)(1-x^25)(1-x^50)(1-x^100))
         = 139946140451 (ways to make change; two DPs in different coin
         order; second: contour |x| = 99/100, 30000 points, 80 digits).
         Machine instance [x^100] = 293.
   b4.5  Motzkin number M_1000 = [x^1000] (1 - x - sqrt(1-2x-3x^2))/(2x^2),
         473 digits, 6.11327659767718...e472 (three-term recurrence; second:
         sum_k C(n,2k) Catalan(k)).  Machine instance M_30 = 1697385471211.
         Contour |x| = 3/10 (inside the branch point 1/3) loses 52 / 1.5 digits.
   b4.6  Ramanujan tau(1000) = [q^1000] q prod (1-q^k)^24 = -30328412970240000
         (integer power series; second: Hecke multiplicativity
         tau(8) tau(125) = 84480 * (-359001100500)).  Contour |q| = 99/100
         loses 2 digits, so the machine case is the FULL instance.
   b4.7  Res_{z=0} z^20 exp(z + 1/z) = I_21(2) = 2.04822323834387801884284e-20
         (essential singularity: exp(z + 1/z) = sum I_n(2) z^n; second: the
         series sum 1/(k! (k+21)!)).  Machine instance z^3 exp(z + 1/z):
         I_4(2) = 0.050728569979180238238.  Contour |z| = 1 loses 20 digits.
   b4.8  [x^1000] T_2000(x) = (n/2)(-1)^k (n-k-1)!/(k!(n-2k)!) 2^(n-2k), n = 2000,
         k = 500: a 714-digit integer 7.00127388815131...e713 (closed form;
         second: the integer recurrence T_{n+1} = 2x T_n - T_{n-1}).  Machine
         instance [x^20] T_40 = 21002987765760.  Contour |x| = 3/5 loses 2
         digits, but the integrand reaches 1e715 at n = 2000.
   b4.9  Bell number B_500 = 500! [x^500] exp(e^x - 1), 844 digits,
         1.60607260103999...e843 (Bell triangle; second: Dobinski's series at
         1500 digits).  Machine instance B_30 = 846749014511809332450147.
         Saddle-point contour |x| = 5 (n = 500) and 3 (n = 30) loses 5 / 2.5 digits.
   b4.10 Res_{z=0} pi cot(pi z)/z^40 = -2 zeta(40) = -2.00000000000181898956805
         (pi cot(pi z) = 1/z - 2 sum zeta(2k) z^(2k-1); second: contour at 80
         digits).  Machine instance /z^10: -2 zeta(10) = -2.0019891505...

   CHECKS.  Exact integers {sig[v, 15], v mod 10^9}; exact rationals
   {sig[v, 15], numerator mod 10^9, denominator mod 10^9}; symbolic residues
   sig[N[v, 30], 15]; contours sig[v, d] with d = 15 at arbitrary precision
   (PrecisionGoal 20, digits lost subtracted from the working precision) and
   d = 6 at machine precision (default PrecisionGoal ~8; 3 for the 3-digit
   b4.4; 4 and 5 for b4.3 and b4.5, whose next digit sits on a rounding edge).
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"SeriesCoefficient", "Coefficient", "Residue", "NIntegrate", "ChebyshevT"}];

intChk[v_] := If[IntegerQ[v], {sig[v, 15], Mod[v, 10^9]}, "unevaluated"];
ratChk[v_] := If[IntegerQ[v] || Head[v] === Rational,
  {sig[v, 15], Mod[Numerator[v], 10^9], Mod[Denominator[v], 10^9]}, "unevaluated"];
(* No Chop: the residues are as small as 1e-54.  A closed form evaluated through
   complex intermediates may leave a negligible imaginary part; drop it. *)
symChk[v_] := If[NumericQ[v] && FreeQ[v, Residue | SeriesCoefficient],
  Module[{n = N[v, 30]}, sig[If[Abs[Im[n]] < 10^-20 Abs[n], Re[n], n], 15]], "unevaluated"];
mc[d_][v_] := If[NumericQ[v], sig[v, d], "unevaluated"];

(* [w^n] f(c + w) = (1/2 pi) int_0^{2 pi} f(c + r e^{it}) (r e^{it})^-n dt.  f is a
   numeric-only function; a residue is n = -1. *)
coef[f_, c_, n_, r_, opts___] :=
  Re[NIntegrate[f[c + r Exp[I t]] (r Exp[I t])^(-n), {t, 0, 2 Pi}, opts]]/(2 Pi);
trap[wp_] := Sequence[Method -> "Trapezoidal", WorkingPrecision -> wp, PrecisionGoal -> 20, MaxRecursion -> 16];

(* ---- b4.1 partitions -------------------------------------------------- *)
part1000[z_?NumericQ] := 1/Product[1 - z^k, {k, 1, 1000}];
part100[z_?NumericQ] := 1/Product[1 - z^k, {k, 1, 100}];
run["b4.1 SeriesCoefficient exact", SeriesCoefficient[Product[1/(1 - x^k), {k, 1, 1000}], {x, 0, 1000}], intChk];
run["b4.1 contour wp60", coef[part1000, 0, 1000, 9/10, trap[60]], mc[15]];
run["b4.1 contour machine n=100", coef[part100, 0, 100, 9/10], mc[6], 3];

(* ---- b4.2 secant ------------------------------------------------------ *)
secf[z_?NumericQ] := Sec[z];
run["b4.2 SeriesCoefficient exact", SeriesCoefficient[Sec[x], {x, 0, 400}], ratChk];
run["b4.2 contour wp160", coef[secf, 0, 400, 1, trap[160]], mc[15]];
run["b4.2 contour machine n=20", coef[secf, 0, 20, 1], mc[6], 3];

(* ---- b4.3 triple pole of Gamma^3 -------------------------------------- *)
gam3[z_?NumericQ] := Gamma[z]^3;
run["b4.3 Residue exact", Residue[Gamma[z]^3, {z, -20}], symChk];
run["b4.3 contour wp40", coef[gam3, -20, -1, 1/2, trap[40]], mc[15]];
run["b4.3 contour machine", coef[gam3, -20, -1, 1/2], mc[4], 3];

(* ---- b4.4 making change ----------------------------------------------- *)
coins[z_?NumericQ] := 1/((1 - z) (1 - z^5) (1 - z^10) (1 - z^25) (1 - z^50) (1 - z^100));
run["b4.4 SeriesCoefficient exact",
    SeriesCoefficient[1/((1 - x) (1 - x^5) (1 - x^10) (1 - x^25) (1 - x^50) (1 - x^100)), {x, 0, 10000}], intChk];
run["b4.4 contour wp80", coef[coins, 0, 10000, 99/100, trap[80]], mc[11]];
run["b4.4 contour machine n=100", coef[coins, 0, 100, 9/10], mc[3], 3];

(* ---- b4.5 Motzkin ----------------------------------------------------- *)
motz[z_?NumericQ] := (1 - z - Sqrt[1 - 2 z - 3 z^2])/(2 z^2);
run["b4.5 SeriesCoefficient exact", SeriesCoefficient[(1 - x - Sqrt[1 - 2 x - 3 x^2])/(2 x^2), {x, 0, 1000}], intChk];
run["b4.5 contour wp100", coef[motz, 0, 1000, 3/10, trap[100]], mc[15]];
run["b4.5 contour machine n=30", coef[motz, 0, 30, 3/10], mc[5], 3];

(* ---- b4.6 Ramanujan tau ----------------------------------------------- *)
delta[z_?NumericQ] := z Product[1 - z^k, {k, 1, 1000}]^24;
run["b4.6 SeriesCoefficient exact", SeriesCoefficient[q Product[(1 - q^k)^24, {k, 1, 1000}], {q, 0, 1000}], intChk];
run["b4.6 contour wp40", coef[delta, 0, 1000, 99/100, trap[40]], mc[15]];
run["b4.6 contour machine", coef[delta, 0, 1000, 99/100], mc[6]];

(* ---- b4.7 essential singularity --------------------------------------- *)
ess21[z_?NumericQ] := z^20 Exp[z + 1/z];
ess3[z_?NumericQ] := z^3 Exp[z + 1/z];
run["b4.7 Residue exact", Residue[z^20 Exp[z + 1/z], {z, 0}], symChk];
run["b4.7 contour wp50", coef[ess21, 0, -1, 1, trap[50]], mc[15]];
run["b4.7 contour machine z^3", coef[ess3, 0, -1, 1], mc[6], 3];

(* ---- b4.8 Chebyshev coefficient --------------------------------------- *)
cheb2000[z_?NumericQ] := ChebyshevT[2000, z];
cheb40[z_?NumericQ] := ChebyshevT[40, z];
run["b4.8 Coefficient exact", Coefficient[ChebyshevT[2000, x], x, 1000], intChk];
run["b4.8 contour wp30", coef[cheb2000, 0, 1000, 3/5, trap[30]], mc[15]];
run["b4.8 contour machine n=40", coef[cheb40, 0, 20, 3/5], mc[6], 3];

(* ---- b4.9 Bell numbers ------------------------------------------------ *)
bellf[z_?NumericQ] := Exp[Exp[z] - 1];
run["b4.9 SeriesCoefficient exact", 500! SeriesCoefficient[Exp[Exp[x] - 1], {x, 0, 500}], intChk];
run["b4.9 contour wp40", 500! coef[bellf, 0, 500, 5, trap[40]], mc[15]];
run["b4.9 contour machine n=30", 30! coef[bellf, 0, 30, 3], mc[6], 3];

(* ---- b4.10 41st-order pole -------------------------------------------- *)
picot[z_?NumericQ] := Pi Cot[Pi z]/z^40;
picot10[z_?NumericQ] := Pi Cot[Pi z]/z^10;
run["b4.10 Residue exact", Residue[Pi Cot[Pi z]/z^40, {z, 0}], symChk];
run["b4.10 contour wp40", coef[picot, 0, -1, 1/2, trap[40]], mc[15]];
run["b4.10 contour machine z^10", coef[picot10, 0, -1, 1/2], mc[6], 3];
