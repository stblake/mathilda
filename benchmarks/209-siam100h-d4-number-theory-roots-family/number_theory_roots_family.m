(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 209 -- SIAM 100-digit HOLDOUT A, family d4 (roots of functions
   of analytic number theory; original: experiment 179, the largest zero of
   Riemann's R(x) = Sum mu(k)/k li(x^(1/k))).

   R(x) is Riemann's prime-counting function.  Its two classical series are
   the Moebius sum Sum_k mu(k)/k li(x^(1/k)) (conditionally and very slowly
   convergent: li(x^(1/k)) -> log log x - log k + gamma) and the Gram series
   1 + Sum_{k>=1} (log x)^k / (k k! zeta(k+1)) (entire, fast).
     hd4.1   li(x) = 0                     LogIntegral built-in   (Soldner)
     hd4.2   R(x) = 100                    RiemannR built-in
     hd4.3   R(x) = 2                      Gram series by NSum
     hd4.4   R(x) = 1000                   Moebius series by NSum (as exp. 179)
     hd4.5   R(x) = 1/2, 0 < x < 1         Gram series by NSum (alternating)
     hd4.6   Sum mu(k)/k^s = 1/2           Moebius Dirichlet series by NSum
     hd4.7   Sum Lambda(n)/n^s = 1         von Mangoldt Dirichlet series by NSum
     hd4.8   first zero of Z(t)            RiemannSiegelZ built-in
     hd4.9   R(x) = 10^10                  Gram series by NSum
     hd4.10  the SECOND largest zero of R  RiemannR built-in, in t = log x
   The answer is the root x (s, t).  Generalised along the original's axis:
   FindRoot on a number-theoretic function that is either a built-in special
   function or an infinite series stated with NSum; varied in kind: the
   level c of R(x) = c on both sides of 1, the two series of R, the Dirichlet
   series of 1/zeta and -zeta'/zeta (slow, sign-indefinite / slowly
   convergent tails), li itself, the Riemann-Siegel Z function, and the next
   zero of R below the record one (x ~ 2.04 10^-15301, the harder member; started
   at t = -35300, magnitude given as in exp. 179).

   Idiomatic statement only: FindRoot with the series as a ?NumericQ helper
   around NSum, machine precision and WorkingPrecision -> 60 (NSum and
   FindRoot).  True values: mpmath (riemannr, li, zeta, siegelz, the zeta-zero
   expansion of R) == wolframscript (RiemannR, LogIntegral, Zeta, ZetaZero)
   at 90 digits: siam100_refs/hd4_ref.py and .wl; truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"FindRoot", "NSum", "LogIntegral", "RiemannR", "MoebiusMu", "MangoldtLambda",
         "Zeta", "RiemannSiegelZ"}];

gram[x_?NumericQ] := 1 + NSum[Log[x]^k/(k k! Zeta[k + 1]), {k, 1, Infinity}];
gram60[x_?NumericQ] := 1 + NSum[Log[x]^k/(k k! Zeta[k + 1]), {k, 1, Infinity},
                                WorkingPrecision -> 60];
moeb[x_?NumericQ] := NSum[MoebiusMu[k]/k LogIntegral[x^(1/k)], {k, 1, Infinity}];
moeb60[x_?NumericQ] := NSum[MoebiusMu[k]/k LogIntegral[x^(1/k)], {k, 1, Infinity},
                            WorkingPrecision -> 60];
imu[s_?NumericQ] := NSum[MoebiusMu[k]/k^s, {k, 1, Infinity}];
imu60[s_?NumericQ] := NSum[MoebiusMu[k]/k^s, {k, 1, Infinity}, WorkingPrecision -> 60];
lam[s_?NumericQ] := NSum[MangoldtLambda[k]/k^s, {k, 1, Infinity}];
lam60[s_?NumericQ] := NSum[MangoldtLambda[k]/k^s, {k, 1, Infinity}, WorkingPrecision -> 60];

vrun["hd4.1 machine", ruleVal[FindRoot[LogIntegral[x] == 0, {x, 2}], x], 10, 1, 60];
vrun["hd4.1 wp50", ruleVal[FindRoot[LogIntegral[x] == 0, {x, 2}, WorkingPrecision -> 60], x], 48, 1, 60];
vrun["hd4.2 machine", ruleVal[FindRoot[RiemannR[x] == 100, {x, 500}], x], 10, 1, 60];
vrun["hd4.2 wp50", ruleVal[FindRoot[RiemannR[x] == 100, {x, 500}, WorkingPrecision -> 60], x], 48, 1, 60];
vrun["hd4.3 machine", ruleVal[FindRoot[gram[x] == 2, {x, 3}], x], 10, 1, 60];
vrun["hd4.3 wp50", ruleVal[FindRoot[gram60[x] == 2, {x, 3}, WorkingPrecision -> 60], x], 48, 1, 60];
vrun["hd4.4 machine", ruleVal[FindRoot[moeb[x] == 1000, {x, 8000}], x], 10, 1, 60];
vrun["hd4.4 wp50", ruleVal[FindRoot[moeb60[x] == 1000, {x, 8000}, WorkingPrecision -> 60], x], 48, 1, 60];
vrun["hd4.5 machine", ruleVal[FindRoot[gram[x] == 1/2, {x, 1/4}], x], 10, 1, 60];
vrun["hd4.5 wp50", ruleVal[FindRoot[gram60[x] == 1/2, {x, 1/4}, WorkingPrecision -> 60], x], 48, 1, 60];
vrun["hd4.6 machine", ruleVal[FindRoot[imu[s] == 1/2, {s, 2}], s], 10, 1, 60];
vrun["hd4.6 wp50", ruleVal[FindRoot[imu60[s] == 1/2, {s, 2}, WorkingPrecision -> 60], s], 48, 1, 60];
vrun["hd4.7 machine", ruleVal[FindRoot[lam[s] == 1, {s, 2}], s], 10, 1, 60];
vrun["hd4.7 wp50", ruleVal[FindRoot[lam60[s] == 1, {s, 2}, WorkingPrecision -> 60], s], 48, 1, 60];
vrun["hd4.8 machine", ruleVal[FindRoot[RiemannSiegelZ[t] == 0, {t, 14}], t], 10, 1, 60];
vrun["hd4.8 wp50", ruleVal[FindRoot[RiemannSiegelZ[t] == 0, {t, 14}, WorkingPrecision -> 60], t], 48, 1, 60];
vrun["hd4.9 machine", ruleVal[FindRoot[gram[x] == 10^10, {x, 2 10^11}], x], 10, 1, 60];
vrun["hd4.9 wp50", ruleVal[FindRoot[gram60[x] == 10^10, {x, 2 10^11}, WorkingPrecision -> 60], x], 48, 1, 60];
vrun["hd4.10 machine", Exp[ruleVal[FindRoot[RiemannR[Exp[t]] == 0, {t, -35300}], t]], 10, 1, 60];
vrun["hd4.10 wp50", Exp[ruleVal[FindRoot[RiemannR[Exp[t]] == 0, {t, -35300}, WorkingPrecision -> 60], t]], 48, 1, 120];
