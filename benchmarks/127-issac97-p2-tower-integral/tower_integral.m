(* bench-timeout: 600 *)
(* ==========================================================================
   Experiment 127 -- ISSAC '97 Systems Challenge, Problem 2
   ==========================================================================

       cd benchmarks/127-issac97-p2-tower-integral
       ../../Mathilda -file tower_integral.m

   PROBLEM.  Integrate[x^x^x, {x, 1, 6}] = 1.1026649993619406556e36300.  The
   integrand climbs from 1 to 6^46656 = 2.66e36305 within the last ~1e-5 of
   the interval, so (a) the answer is far outside the double range and (b) all
   of the mass sits in a boundary layer at the right endpoint.

   METHODS (the notebook's two):
     m1  NIntegrate directly: at machine precision (a system must notice the
         integrand leaves the double range and escalate to arbitrary
         precision), and at WorkingPrecision 30 with Method ->
         "DoubleExponential", MaxRecursion -> 20 (verbatim), plus 50 and 100
         digits for the cost-vs-digits curve;
     m2  the substitution x^x -> z, x = Log[z]/ProductLog[Log[z]], which
         removes the boundary layer; machine and 30 digits.
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"NIntegrate", "ProductLog"}];

run["p2 m1 NIntegrate machine",
    NIntegrate[x^x^x, {x, 1, 6}], sig[#, 6] &, 3];
run["p2 m1 NIntegrate DE wp30",
    NIntegrate[x^x^x, {x, 1, 6}, MaxRecursion -> 20, WorkingPrecision -> 30,
               PrecisionGoal -> 7, Method -> "DoubleExponential"], sig[#, 7] &, 3];
run["p2 m1 NIntegrate wp50",
    NIntegrate[x^x^x, {x, 1, 6}, WorkingPrecision -> 50, PrecisionGoal -> 40],
    sig[#, 30] &, 3];
run["p2 m1 NIntegrate wp100",
    NIntegrate[x^x^x, {x, 1, 6}, WorkingPrecision -> 100, PrecisionGoal -> 80],
    sig[#, 60] &];

sub = (Log[z]/ProductLog[Log[z]])^z D[Log[z]/ProductLog[Log[z]], z];
run["p2 m2 ProductLog substitution machine",
    NIntegrate[Evaluate[sub], {z, 1, 6^6}], sig[#, 6] &, 3];
run["p2 m2 ProductLog substitution wp30",
    NIntegrate[Evaluate[sub], {z, 1, 6^6}, WorkingPrecision -> 30, PrecisionGoal -> 20],
    sig[#, 15] &, 3];
