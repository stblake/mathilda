(* bench-timeout: 1800 *)
(* ==========================================================================
   Experiment 196 -- SIAM 100-Digit Challenge, holdout A for Problem 1
   (family hc1: improper integrals with infinitely many oscillations at an
   endpoint).  Ten new instances of the C1 capability; each is asked as
   written, over the stated interval, at machine precision and at 50 digits
   (WorkingPrecision 60, PrecisionGoal 50) -- exactly the C1 statement:

     hc1.1   Int_0^1 sin(1/x) dx                       bounded, abs. integrable (easy)
     hc1.2   Int_0^1 sin(1/x)/x dx                     conditionally convergent
     hc1.3   Int_0^1 cos(x^-2)/x dx                    conditionally convergent, faster phase
     hc1.4   Int_0^1 x^-1/2 cos(1/x) dx                unbounded amplitude, abs. integrable
     hc1.5   Int_0^1 x^-1 sin(x^-1 log x) dx           the C1 twin (sin for cos)
     hc1.6   Int_0^1 cos(x^-1 log x) dx                C1 phase, bounded amplitude
     hc1.7   Int_0^1 cos(e^(1/x)) dx                   super-exponential phase (overflows
                                                       a double for x < 1/709)
     hc1.8   Int_0^1 sin(1/(x(1-x))) dx                oscillation at BOTH endpoints
     hc1.9   Int_0^2 sin(1/(2-x))/(2-x) dx             at the UPPER endpoint, other interval
     hc1.10  Int_0^1 x^-3/2 cos(x^-1 log x) dx         C1 phase, amplitude x^-3/2 (harder)

   True values (two independent methods each): truth.tsv; reference scripts
   siam100_refs/hc1_ref.{py,wl}.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NIntegrate"}];

hc1[1] = {Sin[1/x], {x, 0, 1}};
hc1[2] = {Sin[1/x]/x, {x, 0, 1}};
hc1[3] = {Cos[1/x^2]/x, {x, 0, 1}};
hc1[4] = {Cos[1/x]/Sqrt[x], {x, 0, 1}};
hc1[5] = {Sin[Log[x]/x]/x, {x, 0, 1}};
hc1[6] = {Cos[Log[x]/x], {x, 0, 1}};
hc1[7] = {Cos[Exp[1/x]], {x, 0, 1}};
hc1[8] = {Sin[1/(x (1 - x))], {x, 0, 1}};
hc1[9] = {Sin[1/(2 - x)]/(2 - x), {x, 0, 2}};
hc1[10] = {Cos[Log[x]/x]/x^(3/2), {x, 0, 1}};

Do[
  With[{f = hc1[k][[1]], rng = hc1[k][[2]]},
    vrun["hc1." <> ToString[k] <> " machine", NIntegrate[f, rng], 10, 1, 60];
    vrun["hc1." <> ToString[k] <> " wp50",
         NIntegrate[f, rng, WorkingPrecision -> 60, PrecisionGoal -> 50], 48, 1, 90]],
  {k, 10}];
