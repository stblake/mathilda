(* bench-timeout: 1800 *)
(* ==========================================================================
   Experiment 211 -- SIAM 100-Digit Challenge HOLDOUT A, family D6
   (generalises experiment 181, "quadruple exponential").

   The capability: a smooth (entire) integrand on a FINITE interval whose
   phase grows super-linearly -- polynomially, exponentially, or as a tower
   of exponentials -- so the oscillations crowd toward one end, at machine
   precision and at 50 digits.  Ten instances, varied in kind:

     hd6.1   Int_0^2   cos(e^(x^2)) dx                     ~17 sign changes; the
                                                          value is -0.00187: heavy
                                                          cancellation (harder)
     hd6.2   Int_-1^1  exp(x + cos(e^(e^x))) dx            two-level tower, mild
     hd6.3   Int_0^1   sin(e^(e^(2x))) dx                  ~515 sign changes (hard)
     hd6.4   Int_0^3   sin(x^3) dx                         polynomial phase
     hd6.5   Int_0^pi  exp(sin(e^x)) dx                    one exponential
     hd6.6   Int_-1^1  sin(e^(4x))/(2+x) dx
     hd6.7   Int_0^2   x e^(x^2) sin(e^(x^2)) dx           GROWING amplitude (to 218);
                                                          closed form (cos 1 - cos e^4)/2
     hd6.8   Int_-1^1  exp(x + sin(e^(e^(x+1/2)))) dx      the original, shifted:
                                                          ~28 sign changes
     hd6.9   Int_0^8   cos(x^2 + x) dx                     Fresnel-type
     hd6.10  Int_0^1   cos(1/(x + 1/100)) dx              oscillation crowded at the
                                                          LEFT end (near-singular)

   Idiomatic statement only, as in 181: the integral as written, NIntegrate
   with default goals (machine) and WorkingPrecision -> 60, PrecisionGoal -> 50
   (wp50).  No splitting at the zeros: that is the reference method
   (siam100_refs/hd6_ref.{py,wl}).  Caps: 30 s per machine run (warm-up + 3
   timed), 60 s per wp50 run.  True values and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NIntegrate"}];


vrun["hd6.1 machine", NIntegrate[Cos[Exp[x^2]], {x, 0, 2}], 10, 3, 30];
vrun["hd6.1 wp50", NIntegrate[Cos[Exp[x^2]], {x, 0, 2}, WorkingPrecision -> 60, PrecisionGoal -> 50], 48, 1, 60];
vrun["hd6.2 machine", NIntegrate[Exp[x + Cos[Exp[Exp[x]]]], {x, -1, 1}], 10, 3, 30];
vrun["hd6.2 wp50", NIntegrate[Exp[x + Cos[Exp[Exp[x]]]], {x, -1, 1}, WorkingPrecision -> 60, PrecisionGoal -> 50], 48, 1, 60];
vrun["hd6.3 machine", NIntegrate[Sin[Exp[Exp[2 x]]], {x, 0, 1}], 10, 3, 30];
vrun["hd6.3 wp50", NIntegrate[Sin[Exp[Exp[2 x]]], {x, 0, 1}, WorkingPrecision -> 60, PrecisionGoal -> 50], 48, 1, 60];
vrun["hd6.4 machine", NIntegrate[Sin[x^3], {x, 0, 3}], 10, 3, 30];
vrun["hd6.4 wp50", NIntegrate[Sin[x^3], {x, 0, 3}, WorkingPrecision -> 60, PrecisionGoal -> 50], 48, 1, 60];
vrun["hd6.5 machine", NIntegrate[Exp[Sin[Exp[x]]], {x, 0, Pi}], 10, 3, 30];
vrun["hd6.5 wp50", NIntegrate[Exp[Sin[Exp[x]]], {x, 0, Pi}, WorkingPrecision -> 60, PrecisionGoal -> 50], 48, 1, 60];
vrun["hd6.6 machine", NIntegrate[Sin[Exp[4 x]]/(2 + x), {x, -1, 1}], 10, 3, 30];
vrun["hd6.6 wp50", NIntegrate[Sin[Exp[4 x]]/(2 + x), {x, -1, 1}, WorkingPrecision -> 60, PrecisionGoal -> 50], 48, 1, 60];
vrun["hd6.7 machine", NIntegrate[x Exp[x^2] Sin[Exp[x^2]], {x, 0, 2}], 10, 3, 30];
vrun["hd6.7 wp50", NIntegrate[x Exp[x^2] Sin[Exp[x^2]], {x, 0, 2}, WorkingPrecision -> 60, PrecisionGoal -> 50], 48, 1, 60];
vrun["hd6.8 machine", NIntegrate[Exp[x + Sin[Exp[Exp[x + 1/2]]]], {x, -1, 1}], 10, 3, 30];
vrun["hd6.8 wp50", NIntegrate[Exp[x + Sin[Exp[Exp[x + 1/2]]]], {x, -1, 1}, WorkingPrecision -> 60, PrecisionGoal -> 50], 48, 1, 60];
vrun["hd6.9 machine", NIntegrate[Cos[x^2 + x], {x, 0, 8}], 10, 3, 30];
vrun["hd6.9 wp50", NIntegrate[Cos[x^2 + x], {x, 0, 8}, WorkingPrecision -> 60, PrecisionGoal -> 50], 48, 1, 60];
vrun["hd6.10 machine", NIntegrate[Cos[1/(x + 1/100)], {x, 0, 1}], 10, 3, 30];
vrun["hd6.10 wp50", NIntegrate[Cos[1/(x + 1/100)], {x, 0, 1}, WorkingPrecision -> 60, PrecisionGoal -> 50], 48, 1, 60];
