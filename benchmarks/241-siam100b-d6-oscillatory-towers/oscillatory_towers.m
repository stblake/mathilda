(* bench-timeout: 1200 *)
(* ==========================================================================
   Experiment 241 -- SIAM 100-Digit Challenge, SEALED HOLDOUT B for Appendix D
   Problem 6 (experiment 181, "quadruple exponential").

   The capability: a smooth (entire or analytic) integrand on a finite
   interval whose phase is a tower of exponentials, so the oscillations crowd
   toward one end -- integrated as written by the system's general-purpose
   quadrature, at machine precision and at 50 digits.  Ten new integrals:

     bd6.1   Int_0^1    exp(sin(e^x))                       (easy, no crowding)
     bd6.2   Int_-1^1   cos(x + e^(e^(x+1/2)))              (~28 sign changes)
     bd6.3   Int_0^2    exp(cos(x^3 + e^x))                 (polynomial+exp phase)
     bd6.4   Int_-1^1   log(2 + sin(e^(e^x)))               (log of a tower)
     bd6.5   Int_0^1    1/(5/4 + cos(e^(4x)))               (near-poles 0.69/(4e^4x) off the axis)
     bd6.6   Int_0^1    sin(e^(e^(2x)))                     (HARD: ~515 sign changes)
     bd6.7   Int_0^Pi   cos(50 sin x - 3x)                  (= Pi BesselJ[3, 50]; periodic)
     bd6.8   Int_-1^1   e^x cos(sinh(e^(x+1/2)))            (sinh tower)
     bd6.9   Int_-1^1/2 cos(e^(e^(e^x)))                    (three exponentials, ~56 sign changes)
     bd6.10  Int_-1^1   exp(x + sin(e^(e^(x+1/2))))         (the original, shifted: phase to e^(e^1.5) = 88)

   Each problem has a machine row (default goals) and a wp50 row
   (WorkingPrecision -> 60, PrecisionGoal -> 50).  No splitting or change of
   variables is applied in any system: the integral as written.
   True values and provenance: truth.tsv (siam100_refs/bd6_ref.{py,gp}).
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NIntegrate"}];

(* {integrand in x, a, b} *)
cases = {
  {Exp[Sin[Exp[x]]],                    0, 1},
  {Cos[x + Exp[Exp[x + 1/2]]],         -1, 1},
  {Exp[Cos[x^3 + Exp[x]]],              0, 2},
  {Log[2 + Sin[Exp[Exp[x]]]],          -1, 1},
  {1/(5/4 + Cos[Exp[4 x]]),             0, 1},
  {Sin[Exp[Exp[2 x]]],                  0, 1},
  {Cos[50 Sin[x] - 3 x],                0, Pi},
  {Exp[x] Cos[Sinh[Exp[x + 1/2]]],     -1, 1},
  {Cos[Exp[Exp[Exp[x]]]],              -1, 1/2},
  {Exp[x + Sin[Exp[Exp[x + 1/2]]]],    -1, 1}};

Do[With[{g = cases[[k, 1]], a = cases[[k, 2]], b = cases[[k, 3]]},
  vrun["bd6." <> ToString[k] <> " machine", NIntegrate[g, {x, a, b}], 10, 3, 60];
  vrun["bd6." <> ToString[k] <> " wp50",
       NIntegrate[g, {x, a, b}, WorkingPrecision -> 60, PrecisionGoal -> 50], 48, 1, 60]],
  {k, 1, 10}];
