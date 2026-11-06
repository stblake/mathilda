(* bench-timeout: 900 *)
(* ==========================================================================
   Experiment 137 -- ISSAC '97 holdout, family 2: hard definite integrals
   ==========================================================================

       cd benchmarks/137-issac97h-hard-integrals
       ../../Mathilda -file hard_integrals.m

   HOLDOUT.  Ten problems probing what experiment 127 (Integrate[x^x^x,
   {x, 1, 6}]) probes -- boundary layers, values far outside the double
   range, endpoint singularities, unbounded oscillation -- on different
   integrands.  See benchmarks/ISSAC97_HOLDOUT.md.

   Each problem has a machine case (default NIntegrate; where the full
   problem leaves the double range, a reduced instance that stays inside it)
   and a WorkingPrecision -> 40 case on the full problem, checked to 15
   digits (at or below the default PrecisionGoal of WP/2).

   PROBLEMS AND TRUE VALUES ("cf" = closed form; "q" = an independent 50-digit
   adaptive Gauss-Legendre quadrature in mpmath, with breakpoints at
   10^-k from the layer; the two agree to every digit shown):
   h2.1  Integrate[x^x^x, {x, 1, 5}] = 1.390195540965462192941711e2180      (q)
         machine: {x, 1, 4} = 1.476250472508459122292813e151               (q)
   h2.2  Integrate[E^E^E^x, {x, 0, 2}] = 4.880478191980581312448124e698    (q)
         machine: {x, 0, 3/2} = 6.207592136617010347600381e35              (q)
   h2.3  Integrate[Sin[1/x], {x, 0, 1}] = Sin[1] - CosIntegral[1]
         = 0.5040670619069283719898561                                    (cf, q)
   h2.4  Integrate[x^-x, {x, 0, 1}] = Sum[n^-n, {n, 1, Infinity}]
         = 1.291285997062663540407283                                     (cf, q)
   h2.5  Integrate[E^(-x^2 - 10^6/x^2), {x, 0, Infinity}] = Sqrt[Pi] E^-2000/2
         = 2.283395465013097564662795e-869                                (cf)
         machine: 100/x^2 = Sqrt[Pi] E^-20/2 = 1.826649837699537925233829e-9 (cf, q)
   h2.6  Integrate[1/(10^-20 + x^2), {x, 0, 1}] = 10^10 ArcTan[10^10]
         = 15707963266.94896619231322                                     (cf, q)
   h2.7  Integrate[Sin[x^3], {x, 0, Infinity}] = Gamma[4/3]/2
         = 0.4464897557846246056092822                                    (cf)
   h2.8  Integrate[x^(-999/1000) E^-x, {x, 0, 1}] = Gamma[1/1000] - Gamma[1/1000, 1]
         = 999.2042906713881633188507                                     (cf)
   h2.9  Integrate[x^1000 E^-x, {x, 0, Infinity}] = 1000!
         = 4.023872600770937735437024e2567                                (cf)
         machine: x^100 = 100! = 9.332621544394415268169924e157           (cf)
   h2.10 Integrate[Cos[1000 Sin[x]], {x, 0, Pi}] = Pi BesselJ[0, 1000]
         = 0.07786967112327907846744718                                   (cf, q)
   Closed forms: h2.5 Int E^(-a x^2 - b/x^2) = Sqrt[Pi/a] E^(-2 Sqrt[a b])/2;
   h2.7 Int Sin[x^n] = Gamma[1 + 1/n] Sin[Pi/(2n)]; h2.8 is the lower
   incomplete gamma function, also summed as Sum[(-1)^k/(k! (k + a))];
   h2.9 Euler's integral for n!.
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"NIntegrate"}];

nm[f_, r_] := NIntegrate[f, r];
nw[f_, r_] := NIntegrate[f, r, WorkingPrecision -> 40];
SetAttributes[{nm, nw}, HoldAll];

run["h2.1 NIntegrate machine [1,4]", nm[x^x^x, {x, 1, 4}], sig[#, 6] &, 3];
run["h2.1 NIntegrate wp40 [1,5]",    nw[x^x^x, {x, 1, 5}], sig[#, 15] &];
run["h2.2 NIntegrate machine [0,3/2]", nm[E^E^E^x, {x, 0, 3/2}], sig[#, 6] &, 3];
run["h2.2 NIntegrate wp40 [0,2]",      nw[E^E^E^x, {x, 0, 2}], sig[#, 15] &];
run["h2.3 NIntegrate machine", nm[Sin[1/x], {x, 0, 1}], sig[#, 6] &, 3];
run["h2.3 NIntegrate wp40",    nw[Sin[1/x], {x, 0, 1}], sig[#, 15] &];
run["h2.4 NIntegrate machine", nm[x^-x, {x, 0, 1}], sig[#, 6] &, 3];
run["h2.4 NIntegrate wp40",    nw[x^-x, {x, 0, 1}], sig[#, 15] &];
run["h2.5 NIntegrate machine b=100",  nm[E^(-x^2 - 100/x^2), {x, 0, Infinity}], sig[#, 6] &, 3];
run["h2.5 NIntegrate wp40 b=10^6",    nw[E^(-x^2 - 10^6/x^2), {x, 0, Infinity}], sig[#, 15] &];
run["h2.6 NIntegrate machine", nm[1/(10^-20 + x^2), {x, 0, 1}], sig[#, 6] &, 3];
run["h2.6 NIntegrate wp40",    nw[1/(10^-20 + x^2), {x, 0, 1}], sig[#, 15] &];
run["h2.7 NIntegrate machine", nm[Sin[x^3], {x, 0, Infinity}], sig[#, 6] &, 3];
run["h2.7 NIntegrate wp40",    nw[Sin[x^3], {x, 0, Infinity}], sig[#, 15] &];
run["h2.8 NIntegrate machine", nm[x^(-999/1000) E^-x, {x, 0, 1}], sig[#, 6] &, 3];
run["h2.8 NIntegrate wp40",    nw[x^(-999/1000) E^-x, {x, 0, 1}], sig[#, 15] &];
run["h2.9 NIntegrate machine x^100", nm[x^100 E^-x, {x, 0, Infinity}], sig[#, 6] &, 3];
run["h2.9 NIntegrate wp40 x^1000",   nw[x^1000 E^-x, {x, 0, Infinity}], sig[#, 15] &];
run["h2.10 NIntegrate machine", nm[Cos[1000 Sin[x]], {x, 0, Pi}], sig[#, 6] &, 3];
run["h2.10 NIntegrate wp40",    nw[Cos[1000 Sin[x]], {x, 0, Pi}], sig[#, 15] &];
