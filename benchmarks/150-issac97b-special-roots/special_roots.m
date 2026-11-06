(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 150 -- ISSAC '97 holdout B, family 5: zeros of high-degree
   special polynomials and special functions (the counterpart of Problem 5,
   the largest zero of LaguerreL[1000, x])
   ==========================================================================

       cd benchmarks/150-issac97b-special-roots
       ../../Mathilda -file special_roots.m

   Other orthogonal families (Hermite, Legendre, Jacobi) at other ends of the
   spectrum (the largest zero, the smallest positive one, one crowded
   against x = 1), zeros of transcendental functions far out (the 1000th
   Bessel zero, the 1000th Riemann zero, the 500th Airy zero), an inverse of
   Gamma, the real roots of Wilkinson's perturbed polynomial, a Mignotte
   polynomial's pair of roots 1.4e-21 apart, and a complex zero of Erf.
   True values, each established twice without either benchmarked system:

   b5.1  largest zero of H_1000 = 44.20915249799639770159036338329
         (n = 100: 13.40648733814491013849801535865).  (a) Sturm-count
         bisection of the Hermite Jacobi matrix (off-diagonal sqrt(k/2)) at 70
         digits; (b) Newton on the three-term recurrence (mpmath hermite).
   b5.2  smallest positive zero of P_1000 = 0.001570010480083193829005023042
         (a) Sturm bisection of the Legendre Jacobi matrix; (b) Newton.
   b5.3  j_{0,1000} = 3140.807295225078628895545453471
         (a) mpmath besseljzero; (b) McMahon's expansion to beta^-7.
   b5.4  gamma_1000 (the 1000th zero of zeta on the critical line)
         = 1419.4224809459956864659890380799  (a) mpmath zetazero;
         (b) Odlyzko's published table 1419.42248094599568646598903807991.
   b5.5  a_500 (500th zero of Ai) = -177.00925104647920682546615344202
         (a) mpmath airyaizero; (b) the asymptotic series T(3 pi (4k-1)/8)
         to t^-10 (26 digits agree).
   b5.6  the root of Gamma(x) = 10^100: 70.957574457353546136215496627741
         (a) Newton on LogGamma; (b) bisection on Gamma itself.
   b5.7  prod_{k=1}^{20} (x - k) - 2^-23 x^19 (Wilkinson's perturbation):
         10 real roots (exact Sturm count, sympy), the largest
         20.846908101482256914928772892631 (sympy exact isolation refined to
         40 digits; mpmath polyroots at 60 digits).
   b5.8  largest zero of the Jacobi polynomial P^(3/2,-1/2)_400
         = 0.99993721886612628233934595246627 (1 - x = 6.278e-5).
         (a) Sturm bisection of the Jacobi-weight Jacobi matrix; (b) Newton
         on mpmath jacobi.
   b5.9  x^n - 2(10x - 1)^2 (Mignotte): two real roots in (0, 1/5) straddling
         1/10.  Their separation: n = 8: 1.41421371793662116917520670534e-5;
         n = 40: 1.41421356237309504880168872421e-21 (= sqrt2 1e-21 to 20 digits: x = (1 +- x^20/sqrt2)/10).  (a) exact real-root
         isolation (sympy, intervals of width 1e-60); (b) Newton from each
         isolating interval.
   b5.10 first zero of Erf in the first quadrant
         = 1.4506161632436755914236098368742 + 1.8809430001533153719412024587837 I
         (a) mpmath findroot (secant); (b) Newton with the exact derivative
         2/sqrt(pi) e^(-z^2) from another start; agrees with the 8-digit
         table of Fettis, Caslin and Cramer (1973).

   METHODS.  FindRoot at machine precision and at WorkingPrecision 40 or 60
   (default goals = half the working precision); the closed-form zero heads
   (BesselJZero, ZetaZero, AiryAiZero) evaluated to 30 digits where the
   problem has one; CountRoots exactly; NSolve at 40 / 80 digits.  Machine
   instances: H_100 (|H_1000| ~ 1e1400 near its largest zero), Mignotte
   n = 8 (the n = 40 pair is closer than a double resolves); the others run
   the full instance in a double.
   CHECKS.  Machine 8-10 digits (3 for the Mignotte gap, whose roots a
   double carries to ~sqrt(eps)); 18 at 40 digits, 25 at 60, 14 at 30, 14
   for the 80-digit Mignotte gap (its 15th digit sits on a rounding edge).
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"FindRoot", "HermiteH", "LegendreP", "JacobiP", "BesselJ", "BesselJZero", "RiemannSiegelZ",
         "ZetaZero", "AiryAi", "AiryAiZero", "CountRoots", "NSolve"}];

mc[d_][v_] := If[NumericQ[v], sig[v, d], "unevaluated"];
(* a real zero returned with a negligible imaginary part counts as real *)
realOf[v_] := If[NumericQ[v] && Abs[Im[v]] < 10^-8 Abs[v], Re[v], v];
rv[d_][r_] := sig[realOf[ruleVal[r, x]], d];
cpx[d_][v_] := If[NumericQ[v], {sig[Re[v], d], sig[Im[v], d]}, "unevaluated"];

h1000[(x_)?NumericQ] := HermiteH[1000, x];
h100[(x_)?NumericQ] := HermiteH[100, x];
leg[(x_)?NumericQ] := LegendreP[1000, x];
jac[(x_)?NumericQ] := JacobiP[400, 3/2, -1/2, x];

run["b5.1 FindRoot wp60", FindRoot[h1000[x], {x, 45}, WorkingPrecision -> 60], rv[25], 3];
run["b5.1 FindRoot machine n=100", FindRoot[h100[x], {x, 14}], rv[10], 3];

run["b5.2 FindRoot wp40", FindRoot[leg[x], {x, 16/10000}, WorkingPrecision -> 40], rv[18], 3];
run["b5.2 FindRoot machine", FindRoot[leg[x], {x, 0.0016}], rv[10], 3];

run["b5.3 FindRoot wp40", FindRoot[BesselJ[0, x], {x, 31408/10}, WorkingPrecision -> 40], rv[18], 3];
run["b5.3 FindRoot machine", FindRoot[BesselJ[0, x], {x, 3140.8}], rv[10], 3];
run["b5.3 BesselJZero wp30", N[BesselJZero[0, 1000], 30], mc[14], 3];

run["b5.4 FindRoot wp30", FindRoot[RiemannSiegelZ[t], {t, 14194/10}, WorkingPrecision -> 30], sig[realOf[ruleVal[#, t]], 14] &, 3];
run["b5.4 FindRoot machine", FindRoot[RiemannSiegelZ[t], {t, 1419.4}], sig[realOf[ruleVal[#, t]], 10] &, 3];
run["b5.4 ZetaZero wp30", N[Im[ZetaZero[1000]], 30], mc[14], 3];

run["b5.5 FindRoot wp40", FindRoot[AiryAi[x], {x, -177}, WorkingPrecision -> 40], rv[18], 3];
run["b5.5 FindRoot machine", FindRoot[AiryAi[x], {x, -177.}], rv[9], 3];
run["b5.5 AiryAiZero wp30", N[AiryAiZero[500], 30], mc[14], 3];

run["b5.6 FindRoot wp40", FindRoot[Gamma[x] == 10^100, {x, 70}, WorkingPrecision -> 40], rv[18], 3];
run["b5.6 FindRoot machine", FindRoot[Gamma[x] == 10^100, {x, 70.}], rv[10], 3];

wilkP = Product[x - k, {k, 1, 20}] - x^19/2^23;
realMax[sols_] := Module[{r = Select[x /. sols, Abs[Im[#]] < 10^-20 &]}, {Length[r], sig[Max[Re[r]], 18]}];
run["b5.7 CountRoots exact", CountRoots[wilkP, x], # &, 3];
run["b5.7 NSolve wp40", NSolve[wilkP == 0, x, WorkingPrecision -> 40], realMax, 3];
run["b5.7 FindRoot machine", FindRoot[Evaluate[wilkP], {x, 21.}], rv[10], 3];

run["b5.8 FindRoot wp40", FindRoot[jac[x], {x, 1}, WorkingPrecision -> 40], rv[18], 3];
run["b5.8 FindRoot machine", FindRoot[jac[x], {x, 1.}], rv[10], 3];

mig[n_] := x^n - 2 (10 x - 1)^2;
gapNear[sols_, d_] := Module[{r = Sort[Re[Select[x /. sols, Abs[Im[#]] < 10^-30 && 0 < Re[#] < 1/5 &]]]},
  If[Length[r] == 2, sig[r[[2]] - r[[1]], d], {"roots in (0, 1/5)", Length[r]}]];
run["b5.9 CountRoots exact", CountRoots[mig[40], {x, 0, 1/5}], # &, 3];
run["b5.9 NSolve wp80", NSolve[mig[40] == 0, x, WorkingPrecision -> 80], gapNear[#, 14] &, 3];
run["b5.9 NSolve machine n=8", NSolve[mig[8] == 0, x], gapNear[#, 3] &, 3];

run["b5.10 FindRoot wp30", FindRoot[Erf[z], {z, 145/100 + 188/100 I}, WorkingPrecision -> 30], cpx[14][ruleVal[#, z]] &, 3];
run["b5.10 FindRoot machine", FindRoot[Erf[z], {z, 1.45 + 1.88 I}], cpx[8][ruleVal[#, z]] &, 3];
