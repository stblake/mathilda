(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 239 -- SIAM 100-Digit Challenge, sealed holdout B for Appendix D
   problem 4 (largest zero of Riemann's R): ten other roots of functions
   given by special functions or by slowly convergent / cancelling series,
   found with FindRoot.

     bd4.1   LogIntegral[x] == 0                    x = 1.451369234... (Soldner)
     bd4.2   NSum[n^-s, n >= 1] == 2                s = 1.728647238...
     bd4.3   NSum[n^-s/Log[n], n >= 2] == 1         s = 1.576644639...
     bd4.4   RiemannR[Exp[t]] == 10^-3              t = -27.64917570...
     bd4.5   NSum[MoebiusMu[k]/k LogIntegral[Exp[t/k]], k >= 1] == 1/10
                                                    t = -4.569900637...
     bd4.6   NSum[MoebiusMu[k] x^k, k >= 1] == 0, smallest x > 0
                                                    x = 0.5802946238...
     bd4.7   RiemannSiegelZ[t] == 0, first zero     t = 14.13472514...
     bd4.8   PrimeZetaP[s] == 1                     s = 1.399433328...
     bd4.9   RiemannR[Exp[t]] == 10^-9              t = -2798.755341...
     bd4.10  NSum[-Log[n] n^-s, n >= 1] == -1       s = 1.970017855...  (zeta'(s) = -1)

   The answer recorded for bd4.4/5/9 is t = log x itself (x underflows a
   double for bd4.9).  The series rows (bd4.2/3/5/6/10) state the series as
   written inside FindRoot (a ?NumericQ helper around NSum): bd4.2/3/10
   converge like n^(1-s) with s < 2, bd4.5 is the conditionally convergent
   Moebius-li series of experiment 179 at a moderate x, bd4.6 converges
   geometrically.  bd4.4 and bd4.9 are the built-in RiemannR in the region
   where its Gram series cancels from e^u down to R (u = 28: 12 digits lost;
   u = 2799: about 1215 digits lost) -- bd4.9 is the hard one.
   Rows: machine precision; WorkingPrecision -> 60 (FindRoot and NSum).
   True values and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"FindRoot", "NSum", "LogIntegral", "MoebiusMu", "RiemannR", "RiemannSiegelZ",
         "PrimeZetaP"}];

s2[s_?NumericQ] := NSum[n^-s, {n, 1, Infinity}];
s2w[s_?NumericQ] := NSum[n^-s, {n, 1, Infinity}, WorkingPrecision -> 60];
s3[s_?NumericQ] := NSum[n^-s/Log[n], {n, 2, Infinity}];
s3w[s_?NumericQ] := NSum[n^-s/Log[n], {n, 2, Infinity}, WorkingPrecision -> 60];
s5[t_?NumericQ] := NSum[MoebiusMu[k]/k LogIntegral[Exp[t/k]], {k, 1, Infinity}];
s5w[t_?NumericQ] := NSum[MoebiusMu[k]/k LogIntegral[Exp[t/k]], {k, 1, Infinity},
                         WorkingPrecision -> 60];
s6[x_?NumericQ] := NSum[MoebiusMu[k] x^k, {k, 1, Infinity}];
s6w[x_?NumericQ] := NSum[MoebiusMu[k] x^k, {k, 1, Infinity}, WorkingPrecision -> 60];
s10[s_?NumericQ] := NSum[-Log[n] n^-s, {n, 1, Infinity}];
s10w[s_?NumericQ] := NSum[-Log[n] n^-s, {n, 1, Infinity}, WorkingPrecision -> 60];

vrun["bd4.1 machine", ruleVal[FindRoot[LogIntegral[x] == 0, {x, 1.5}], x], 10, 1, 60];
vrun["bd4.1 wp50", ruleVal[FindRoot[LogIntegral[x] == 0, {x, 3/2}, WorkingPrecision -> 60], x], 48, 1, 60];
vrun["bd4.2 machine", ruleVal[FindRoot[s2[s] == 2, {s, 1.7}], s], 10, 1, 60];
vrun["bd4.2 wp50", ruleVal[FindRoot[s2w[s] == 2, {s, 17/10}, WorkingPrecision -> 60], s], 48, 1, 120];
vrun["bd4.3 machine", ruleVal[FindRoot[s3[s] == 1, {s, 1.6}], s], 10, 1, 60];
vrun["bd4.3 wp50", ruleVal[FindRoot[s3w[s] == 1, {s, 8/5}, WorkingPrecision -> 60], s], 48, 1, 120];
vrun["bd4.4 machine", ruleVal[FindRoot[RiemannR[Exp[t]] == 1/1000, {t, -28.}], t], 10, 1, 60];
vrun["bd4.4 wp50", ruleVal[FindRoot[RiemannR[Exp[t]] == 1/1000, {t, -28}, WorkingPrecision -> 60], t], 48, 1, 60];
vrun["bd4.5 machine", ruleVal[FindRoot[s5[t] == 1/10, {t, -5.}], t], 10, 1, 60];
vrun["bd4.5 wp50", ruleVal[FindRoot[s5w[t] == 1/10, {t, -5}, WorkingPrecision -> 60], t], 48, 1, 120];
vrun["bd4.6 machine", ruleVal[FindRoot[s6[x] == 0, {x, 0.6}], x], 10, 1, 60];
vrun["bd4.6 wp50", ruleVal[FindRoot[s6w[x] == 0, {x, 3/5}, WorkingPrecision -> 60], x], 48, 1, 120];
vrun["bd4.7 machine", ruleVal[FindRoot[RiemannSiegelZ[t] == 0, {t, 14.}], t], 10, 1, 60];
vrun["bd4.7 wp50", ruleVal[FindRoot[RiemannSiegelZ[t] == 0, {t, 14}, WorkingPrecision -> 60], t], 48, 1, 60];
vrun["bd4.8 machine", ruleVal[FindRoot[PrimeZetaP[s] == 1, {s, 1.4}], s], 10, 1, 60];
vrun["bd4.8 wp50", ruleVal[FindRoot[PrimeZetaP[s] == 1, {s, 7/5}, WorkingPrecision -> 60], s], 48, 1, 60];
vrun["bd4.9 machine", ruleVal[FindRoot[RiemannR[Exp[t]] == 10^-9, {t, -2800.}], t], 10, 1, 60];
vrun["bd4.9 wp50", ruleVal[FindRoot[RiemannR[Exp[t]] == 10^-9, {t, -2800}, WorkingPrecision -> 60], t], 48, 1, 120];
vrun["bd4.10 machine", ruleVal[FindRoot[s10[s] == -1, {s, 1.9}], s], 10, 1, 60];
vrun["bd4.10 wp50", ruleVal[FindRoot[s10w[s] == -1, {s, 19/10}, WorkingPrecision -> 60], s], 48, 1, 120];
