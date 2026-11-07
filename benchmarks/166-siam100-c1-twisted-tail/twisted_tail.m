(* bench-timeout: 600 *)
(* ==========================================================================
   Experiment 166 -- SIAM 100-Digit Challenge, Problem 1 ("A Twisted Tail").

       lim_{eps -> 0}  Int_eps^1  x^-1 cos(x^-1 log x) dx  =  0.3233674316...

   Idiomatic statement only: the integral as written, over (0, 1], at machine
   precision and at 50 digits (WorkingPrecision 60, PrecisionGoal 50).  The
   integrand oscillates infinitely often as x -> 0 with an amplitude that
   GROWS like 1/x; the improper integral converges only conditionally.
   True value and provenance: truth.tsv (and SIAM100.md).
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NIntegrate"}];

vrun["c1.machine", NIntegrate[Cos[Log[x]/x]/x, {x, 0, 1}], 10, 3];
vrun["c1.wp50", NIntegrate[Cos[Log[x]/x]/x, {x, 0, 1},
                           WorkingPrecision -> 60, PrecisionGoal -> 50], 48];
