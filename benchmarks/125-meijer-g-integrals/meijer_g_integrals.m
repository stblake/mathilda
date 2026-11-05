(* ==========================================================================
   Experiment 125 -- definite integrals by Meijer G functions
   ==========================================================================

       cd benchmarks/125-meijer-g-integrals
       ../../Mathilda -file meijer_g_integrals.m

   WHAT IT MEASURES.  The Meijer G machinery: Integrate[..., Method -> "Roach"]
   (src/internal/RoachIntegrate.m -- the function -> G table, the PBM 2.24.1
   Mellin/convolution theorems with their validity conditions, then G -> named
   functions by Roach's Formula algorithm, src/special_functions/hrep_core.c),
   FunctionExpand of a pFq (the same Formula algorithm), and 30-digit
   HypergeometricPFQ / MeijerG numerics.

   WHY THESE CASES.  A representative subset of tests/RoachCorpus.m -- the
   test suites of Kelly Roach's 2008 Quantum paper and his ISSAC '96 / '97
   galleries -- chosen from the rows where a comparable engine exists:

     roach ...        12 integrals over (0, oo), the 2008 tables and the '97
                      gallery at sample parameters; baseline sympy's
                      meijerint (integrate(..., meijerg=True)), the same
                      Meijer G convolution method;
     hyperexpand ...  4 pFq of the ISSAC '96 gallery; baseline
                      sympy.hyperexpand against FunctionExpand;
     N ...            a 3F2 at z = 1 and a G^{1,4}_{4,2} at 30 digits;
                      baseline mpmath.

   The parameters are numbers, not symbols with assumptions: sympy's
   meijerint carries symbolic conditions through every step, and a symbolic
   row would time its condition simplifier rather than the method.

   VALUE CHECKS.  Round[10^6 value] of the closed form (the pFq rows at
   z = 3/10).  A row a system leaves unevaluated checks as "unevaluated",
   never as the number its N would produce by quadrature.

   KNOWN GAPS (2026-11-05), which the job reports as INCOMPLETE:
     roach x^(-1/2) Exp[-x] BesselJ[0,x]   the answer keeps two 2F1 at -1
                                           that FunctionExpand does not close;
     N 3F2(-15/4,1,1; -7/2,2; 1)           a 3F2 at z = 1 has no numeric
                                           value in Mathilda (Saalschutz-type
                                           continuation missing).
   ========================================================================== *)

Get["../harness.m"];

require[{"Integrate", "FunctionExpand", "HypergeometricPFQ", "MeijerG", "MeijerGReduce"}];

chk[v_] := If[!FreeQ[v, Integrate | HypergeometricPFQ | MeijerG | Hypergeometric1F1 | Hypergeometric2F1],
  "unevaluated", Round[10^6 N[v]]];
chkz[v_] := chk[v /. z -> 3/10];

SetAttributes[roach, HoldAll];
roach[label_String, f_] := (
  bench[label, Integrate[f, {x, 0, Infinity}, Method -> "Roach"];];
  check[label, chk[Integrate[f, {x, 0, Infinity}, Method -> "Roach"]]]);

(* ---- Roach 2008, tables 1 and 2, and the introduction ---------------- *)

roach["roach Exp[-3x/2] BesselJ[0,x]^2", E^(-3 x/2) BesselJ[0, x]^2];
roach["roach x^(-1/4)/Sqrt[x^4+1]", x^(-1/4)/Sqrt[x^4 + 1]];
roach["roach x^(-1/2) (1+2x)^(-3/2)", x^(-1/2) (1 + 2 x)^(-3/2)];
roach["roach x^(-1/2) (3+x)^-2", x^(-1/2) (3 + x)^-2];
roach["roach BesselJ[1,x] Exp[-x^2/2]", BesselJ[1, x] E^(-x^2/2)];
roach["roach Exp[-x] Sqrt[x]/(2+x)", E^(-x) Sqrt[x]/(2 + x)];
roach["roach x Erfc[x] Exp[-x^2/2]", x Erfc[x] E^(-x^2/2)];
roach["roach (Si[2x]-Pi/2) Exp[-x]", (SinIntegral[2 x] - Pi/2) E^(-x)];

(* ---- ISSAC '97 gallery, and two Bessel convolutions ------------------- *)

roach["roach Sin[x]/(x^2+4)^(3/4)", Sin[x]/(x^2 + 4)^(3/4)];
roach["roach BesselJ[1,x] BesselJ[1,2x]/x^2", BesselJ[1, x] BesselJ[1, 2 x]/x^2];
roach["roach x^2 Exp[-x^2] BesselJ[0,x]", x^2 E^(-x^2) BesselJ[0, x]];
roach["roach x^(-1/2) Exp[-x] BesselJ[0,x]", x^(-1/2) E^(-x) BesselJ[0, x]];

(* ---- ISSAC '96 gallery: pFq -> named functions ------------------------ *)

SetAttributes[hexp, HoldAll];
hexp[label_String, h_] := (
  bench[label, FunctionExpand[h];];
  check[label, chkz[FunctionExpand[h]]]);

hexp["hyperexpand 1F1(1/2; 9/2; z)", HypergeometricPFQ[{1/2}, {9/2}, z]];
hexp["hyperexpand 2F1(-3/2,-1/2; 2; z)", HypergeometricPFQ[{-3/2, -1/2}, {2}, z]];
hexp["hyperexpand 3F2(1,2,3; 1/2,4; z)", HypergeometricPFQ[{1, 2, 3}, {1/2, 4}, z]];
hexp["hyperexpand 1F2(3/2; 5/2,5; z)", HypergeometricPFQ[{3/2}, {5/2, 5}, z]];

(* ---- 30-digit numerics ------------------------------------------------ *)

bench["N 3F2(-15/4,1,1; -7/2,2; 1) 30 digits", N[HypergeometricPFQ[{-15/4, 1, 1}, {-7/2, 2}, 1], 30];];
check["N 3F2(-15/4,1,1; -7/2,2; 1) 30 digits", chk[N[HypergeometricPFQ[{-15/4, 1, 1}, {-7/2, 2}, 1], 30]]];
bench["N MeijerG G^{1,4}_{4,2}(4) 30 digits", N[MeijerG[{{1/2, 1, 1, 1}, {}}, {{1}, {0}}, 4], 30];];
check["N MeijerG G^{1,4}_{4,2}(4) 30 digits", chk[N[MeijerG[{{1/2, 1, 1, 1}, {}}, {{1}, {0}}, 4], 30]]];
