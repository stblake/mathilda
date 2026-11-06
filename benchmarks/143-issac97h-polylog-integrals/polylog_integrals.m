(* bench-timeout: 900 *)
(* ==========================================================================
   Experiment 143 -- ISSAC '97 holdout, family 8: polylogarithmic integrals
   ==========================================================================

       cd benchmarks/143-issac97h-polylog-integrals
       ../../Mathilda -file polylog_integrals.m

   HOLDOUT.  Ten problems probing what experiment 133
   (Integrate[x^2 PolyLog[3, 1/(x + 1)], {x, 0, 1}]) probes -- definite
   integrals of polynomial/rational x PolyLog/Log/ArcTan whose closed forms
   live in Zeta, PolyLog, Log, Pi and Catalan -- on different integrands.
   See benchmarks/ISSAC97_HOLDOUT.md.

   Each problem has three cases: exact Integrate (the result must be a
   closed form free of Integrate; checked numerically to 20 digits), and its
   numeric twin NIntegrate at machine precision (8 digits) and at
   WorkingPrecision -> 30 (14 digits).

   PROBLEMS AND TRUE VALUES.  Each closed form was derived by hand -- termwise
   integration of the power series (Int_0^1 x^(n-1) Log[x]^k = (-1)^k k!/n^(k+1))
   followed by partial fractions in n, or integration by parts -- and checked
   against a 40-digit adaptive quadrature (mpmath); the two agree to >= 38
   digits:
   h8.1  Int_0^1 Log[x]^2 Log[1 - x]/x = -2 Zeta[4] = -Pi^4/45
         = -2.164646467422276383032007393
   h8.2  Int_0^1 x^3 PolyLog[3, x] = Sum[1/(n^3 (n + 4))]
         = Zeta[3]/4 - Pi^2/96 + 25/768 = 0.2302579299452177524037419258
   h8.3  Int_0^1 Log[1 + x]^3/x = Pi^4/15 + Pi^2 Log[2]^2/4 - Log[2]^4/4
         - 6 PolyLog[4, 1/2] - 21/4 Zeta[3] Log[2] = 0.1425141979357109157087215012
   h8.4  Int_0^1 ArcTan[x] Log[x]/x = -Sum[(-1)^k/(2k + 1)^3] = -Pi^3/32
         = -0.9689461462593693804836348458
   h8.5  Int_0^1 Log[1 + x^2]/(1 + x^2) = Pi Log[2]/2 - Catalan
         = 0.1728274509745820501957409342
   h8.6  Int_0^1 x PolyLog[2, 1/(1 + x)] = Pi^2/12 + 1/4 - Log[2]
         (parts with (x^2/2)' and d/du PolyLog[2, 1/u] = Log[1 - 1/u]/u)
         = 0.3793198528641679088189754619
   h8.7  Int_0^1 Log[x]^3/(1 - x) = -6 Zeta[4] = -Pi^4/15
         = -6.493939402266829149096022179
   h8.8  Int_0^1 Log[1 - x] Log[1 + x]/x = -5 Zeta[3]/8
         = -0.7512855644747464283748363509
   h8.9  Int_0^Infinity Log[x]^2/(1 + x^2) = Pi^3/8 = 3.875784585037477521934539383
   h8.10 Int_0^1 x PolyLog[3, -x] = Sum[(-1)^n/(n^3 (n + 2))]
         = -3 Zeta[3]/8 + Pi^2/48 - 1/16 = -0.3076545803288195524658499147
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"Integrate", "NIntegrate", "PolyLog", "Catalan"}];

(* Chop: a closed form assembled from complex-valued pieces can carry a
   0``k I. *)
exactChk[v_] := If[FreeQ[v, Integrate] && NumericQ[v], sig[Chop[N[v, 40], 10^-30], 20], "unevaluated"];

problems = {
 {"h8.1",  Log[x]^2 Log[1 - x]/x,         {x, 0, 1}},
 {"h8.2",  x^3 PolyLog[3, x],             {x, 0, 1}},
 {"h8.3",  Log[1 + x]^3/x,                {x, 0, 1}},
 {"h8.4",  ArcTan[x] Log[x]/x,            {x, 0, 1}},
 {"h8.5",  Log[1 + x^2]/(1 + x^2),        {x, 0, 1}},
 {"h8.6",  x PolyLog[2, 1/(1 + x)],       {x, 0, 1}},
 {"h8.7",  Log[x]^3/(1 - x),              {x, 0, 1}},
 {"h8.8",  Log[1 - x] Log[1 + x]/x,       {x, 0, 1}},
 {"h8.9",  Log[x]^2/(1 + x^2),            {x, 0, Infinity}},
 {"h8.10", x PolyLog[3, -x],              {x, 0, 1}}};

Do[With[{lab = p[[1]], f = p[[2]], r = p[[3]]},
   run[lab <> " Integrate exact", Integrate[f, r], exactChk];
   run[lab <> " NIntegrate machine", NIntegrate[f, r], sig[#, 8] &, 3];
   run[lab <> " NIntegrate wp30", NIntegrate[f, r, WorkingPrecision -> 30], sig[#, 14] &, 3]],
  {p, problems}];
