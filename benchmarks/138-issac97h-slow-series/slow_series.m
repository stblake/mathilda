(* bench-timeout: 900 *)
(* ==========================================================================
   Experiment 138 -- ISSAC '97 holdout, family 3: slowly convergent series
   ==========================================================================

       cd benchmarks/138-issac97h-slow-series
       ../../Mathilda -file slow_series.m

   HOLDOUT.  Ten problems probing what experiment 128 (a series with terms
   ~ n^-1.047) probes -- sums whose value rests entirely on the tail
   estimate: exponents barely above 1, logarithmic factors, and
   conditionally convergent alternating and trigonometric series.  See
   benchmarks/ISSAC97_HOLDOUT.md.

   Each problem: NSum at machine precision (checked to 6 digits) and NSum
   with WorkingPrecision -> 30 (checked to 14 digits, below the default
   PrecisionGoal of WP/2).  A double can carry every answer, so the machine
   case is the full problem.

   PROBLEMS AND TRUE VALUES (two independent routes, agreeing to all digits
   shown; "cf" closed form, "em" Euler-Maclaurin: the first N terms summed
   directly plus the tail integral (in closed form) and the Bernoulli
   corrections, at two different N; "hz" the summand expanded in powers of n
   for n >= N and summed term by term as Hurwitz zeta values, at two N):
   h3.1  Sum[n^(-101/100), {n, 1, Infinity}] = Zeta[101/100]
         = 100.5779433384968724902821543                         (cf, em)
   h3.2  Sum[1/(n Log[n]^2), {n, 2, Infinity}]
         = 2.109742801236891974479257198                         (em N = 10^4, 2 10^4)
   h3.3  Sum[(-1)^(n+1)/Sqrt[n], {n, 1, Infinity}] = (1 - Sqrt[2]) Zeta[1/2]
         = 0.6048986434216303702472659142                        (cf, alternating zeta)
   h3.4  Sum[(-1)^n Log[n]/n, {n, 1, Infinity}] = EulerGamma Log[2] - Log[2]^2/2
         = 0.1598689037424309717569478703                        (cf, Levin-type nsum)
   h3.5  Sum[1/(n^(11/10) + Sqrt[n] + 1), {n, 1, Infinity}]
         = 8.957773692316167741131746530                         (hz N = 2000, 4000)
   h3.6  Sum[Sin[n]/n, {n, 1, Infinity}] = (Pi - 1)/2
         = 1.070796326794896619231321692                         (cf, nsum)
   h3.7  Sum[Cos[n]/Sqrt[n], {n, 1, Infinity}] = Re[PolyLog[1/2, E^I]]
         = -0.1941089350921820497391492449                       (cf via the
           Hurwitz-zeta reflection formula for PolyLog, and direct PolyLog)
   h3.8  Sum[1/(n Log[n] Log[Log[n]]^2), {n, 3, Infinity}]   (tail ~ 1/Log[Log[N]])
         = 38.40676809282178631849374770                         (em N = 10^4, 4 10^4)
   h3.9  Sum[(-1)^(n+1)/n^(1/10), {n, 1, Infinity}] = (1 - 2^(9/10)) Zeta[1/10]
         = 0.5222702824645705054936517728                        (cf, alternating zeta)
   h3.10 Sum[1/(n^(21/20) + 1), {n, 1, Infinity}]
         = 1/2 + Sum[(-1)^k (Zeta[21 (k + 1)/20] - 1), {k, 0, Infinity}]
         = 19.64669717802233418224120349                         (cf, hz N = 1000)
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"NSum"}];

terms = {
  {"h3.1",  n^(-101/100), 1},
  {"h3.2",  1/(n Log[n]^2), 2},
  {"h3.3",  (-1)^(n + 1)/Sqrt[n], 1},
  {"h3.4",  (-1)^n Log[n]/n, 1},
  {"h3.5",  1/(n^(11/10) + Sqrt[n] + 1), 1},
  {"h3.6",  Sin[n]/n, 1},
  {"h3.7",  Cos[n]/Sqrt[n], 1},
  {"h3.8",  1/(n Log[n] Log[Log[n]]^2), 3},
  {"h3.9",  (-1)^(n + 1)/n^(1/10), 1},
  {"h3.10", 1/(n^(21/20) + 1), 1}};

Do[With[{lab = t[[1]], f = t[[2]], n0 = t[[3]]},
   run[lab <> " NSum machine", NSum[f, {n, n0, Infinity}], sig[#, 6] &, 3];
   run[lab <> " NSum wp30", NSum[f, {n, n0, Infinity}, WorkingPrecision -> 30], sig[#, 14] &]],
  {t, terms}];
