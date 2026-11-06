(* bench-timeout: 600 *)
(* ==========================================================================
   Experiment 128 -- ISSAC '97 Systems Challenge, Problem 3
   ==========================================================================

       cd benchmarks/128-issac97-p3-slow-series
       ../../Mathilda -file slow_series.m

   PROBLEM.  Sum[1/(n^Pi + n^2 + n^Sqrt[2] + 1)^(1/3), {n, 1, Infinity}]
   = 21.19324037771154094...  The terms decay like n^(-Pi/3) = n^(-1.047): the
   series converges, barely, so direct summation is hopeless and everything
   rests on the tail estimate (Euler-Maclaurin / integral tail).

   METHODS (the notebook's one, at three precisions):
     machine  NSum with default options (the notebook's "first seven digits");
     wp40     NSum[..., WorkingPrecision -> 40], default method;
     wp120    verbatim: VerifyConvergence -> False, WorkingPrecision -> 120,
              PrecisionGoal -> 16, Method -> "NIntegrate" (integral tail).
   The wp40 row can run for many minutes in a system whose default NSum
   method does not scale with WorkingPrecision, so it is capped ($cap).
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"NSum"}];

term = 1/(n^Pi + n^2 + n^Sqrt[2] + 1)^(1/3);

run["p3 NSum machine", NSum[Evaluate[term], {n, 1, Infinity}], sig[#, 8] &, 3];
run["p3 NSum wp40", NSum[Evaluate[term], {n, 1, Infinity}, WorkingPrecision -> 40],
    sig[#, 30] &];
run["p3 NSum Method->NIntegrate wp120",
    NSum[Evaluate[term], {n, 1, Infinity}, VerifyConvergence -> False,
         WorkingPrecision -> 120, PrecisionGoal -> 16, Method -> "NIntegrate"],
    sig[#, 15] &];
