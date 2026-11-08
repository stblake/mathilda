(* bench-timeout: 1800 *)
(* ==========================================================================
   Experiment 226 -- SIAM 100-Digit Challenge, sealed holdout B for Problem 1
   ("A Twisted Tail"): ten more integrals that oscillate infinitely often at
   an endpoint (x -> 0) or over an infinite range and converge only
   conditionally (or with a slowly decaying envelope).

     bc1.1   Int_0^1 sin(1/x) dx                       = sin 1 - Ci(1)
     bc1.2   Int_0^1 sin(1/x)/x dx                     = pi/2 - Si(1)
     bc1.3   Int_0^1 cos(log(x)/x) dx                  (no closed form)
     bc1.4   Int_0^oo sin(x^2) dx                      = sqrt(pi/8)
     bc1.5   Int_0^oo cos(x)/sqrt(x) dx                = sqrt(pi/2)
     bc1.6   Int_0^1 sin(1/x)/sqrt(x) dx               (Fresnel closed form)
     bc1.7   Int_0^oo sin(x)/(x+1) dx                  = Ci(1) sin 1 + (pi/2 - Si(1)) cos 1
     bc1.8   Int_0^oo J1(x)/sqrt(x) dx                 = Gamma(3/4)/(sqrt 2 Gamma(5/4))
     bc1.9   Int_0^oo cos(x^3) dx                      = Gamma(4/3) cos(pi/6)
     bc1.10  Int_2^oo sin(x)/log(x) dx                 (no closed form; 1/log envelope)

   Idiomatic statement only: the integral as written, NIntegrate at machine
   precision and at WorkingPrecision -> 60, PrecisionGoal -> 50.  Labels
   "bc1.k machine" / "bc1.k wp50".  True values: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NIntegrate", "BesselJ"}];

SetAttributes[c1pair, HoldAll];
c1pair[k_, f_, rng_] := (
  vrun["bc1." <> ToString[k] <> " machine", NIntegrate[f, rng], 10, 3, 60];
  vrun["bc1." <> ToString[k] <> " wp50",
       NIntegrate[f, rng, WorkingPrecision -> 60, PrecisionGoal -> 50], 48, 1, 120]);

c1pair[1, Sin[1/x], {x, 0, 1}];
c1pair[2, Sin[1/x]/x, {x, 0, 1}];
c1pair[3, Cos[Log[x]/x], {x, 0, 1}];
c1pair[4, Sin[x^2], {x, 0, Infinity}];
c1pair[5, Cos[x]/Sqrt[x], {x, 0, Infinity}];
c1pair[6, Sin[1/x]/Sqrt[x], {x, 0, 1}];
c1pair[7, Sin[x]/(x + 1), {x, 0, Infinity}];
c1pair[8, BesselJ[1, x]/Sqrt[x], {x, 0, Infinity}];
c1pair[9, Cos[x^3], {x, 0, Infinity}];
c1pair[10, Sin[x]/Log[x], {x, 2, Infinity}];
