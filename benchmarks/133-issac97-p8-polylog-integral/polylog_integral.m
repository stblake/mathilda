(* bench-timeout: 900 *)
(* ==========================================================================
   Experiment 133 -- ISSAC '97 Systems Challenge, Problem 8
   ==========================================================================

       cd benchmarks/133-issac97-p8-polylog-integral
       ../../Mathilda -file polylog_integral.m

   PROBLEM.  Integrate[x^2 PolyLog[3, 1/(x + 1)], {x, 0, 1}] in closed form,
   0.21088285956559351177...  (a combination of Zeta[3], Pi^2 Log[2],
   Log[2]^3 and rationals).  The 1997 systems could not do it directly; the
   notebook's trick is to simplify the polylog's argument by substitution.

   METHODS (the notebook's):
     exact direct    Integrate[x^2 PolyLog[3, 1/(x + 1)], {x, 0, 1}];
     exact subst1    1/(x + 1) -> z: Integrate[((1 - z)/z)^2 PolyLog[3, z] (-1/z^2), {z, 1, 1/2}];
     exact subst2    x + 1 -> z:     Integrate[(z - 1)^2 PolyLog[3, 1/z], {z, 1, 2}];
     exact indefinite  both substituted antiderivatives (checked through
                       their values at the endpoints);
     numeric  NIntegrate at machine precision, 30 and 60 digits.
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"Integrate", "PolyLog", "NIntegrate", "FunctionExpand"}];

(* Chop: a closed form evaluated through complex-valued pieces (PolyLog at
   an argument > 1 on the way to the endpoint limit) leaves a 0``15 I. *)
exactChk[v_] := If[FreeQ[v, Integrate] && NumericQ[v], sig[Chop[N[v, 30], 10^-20], 15], "unevaluated"];

run["p8 exact Integrate direct", Integrate[x^2 PolyLog[3, 1/(x + 1)], {x, 0, 1}], exactChk];
run["p8 exact Integrate subst 1/(x+1)->z",
    Integrate[((1 - z)/z)^2 PolyLog[3, z] (-(1/z^2)), {z, 1, 1/2}], exactChk];
run["p8 exact Integrate subst x+1->z",
    Integrate[(z - 1)^2 PolyLog[3, 1/z], {z, 1, 2}], exactChk];

indef[F_, a_, b_] := If[FreeQ[F, Integrate], Limit[F, z -> b] - Limit[F, z -> a], F];
run["p8 exact indefinite subst 1/(x+1)->z",
    indef[Simplify[FunctionExpand[Integrate[((1 - z)/z)^2 PolyLog[3, z] (-(1/z^2)), z]]], 1, 1/2],
    exactChk];
run["p8 exact indefinite subst x+1->z",
    indef[Simplify[FunctionExpand[Integrate[(z - 1)^2 PolyLog[3, 1/z], z]]], 1, 2], exactChk];

run["p8 NIntegrate machine", NIntegrate[x^2 PolyLog[3, 1/(x + 1)], {x, 0, 1}], sig[#, 9] &, 3];   (* 9: default PrecisionGoal is ~6-10 digits *)
run["p8 NIntegrate wp30",
    NIntegrate[x^2 PolyLog[3, 1/(x + 1)], {x, 0, 1}, WorkingPrecision -> 30], sig[#, 25] &, 3];
run["p8 NIntegrate wp60",
    NIntegrate[x^2 PolyLog[3, 1/(x + 1)], {x, 0, 1}, WorkingPrecision -> 60], sig[#, 50] &];
