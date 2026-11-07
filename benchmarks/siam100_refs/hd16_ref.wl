(* Reference B for holdout family hd16 (exp 221): heat-equation crossing times
   by the METHOD OF IMAGES -- the free-space heat kernel G(s,t) applied to the
   odd (Dirichlet) / even (Neumann) 2L-periodic extension of the initial data,
   in closed form (Erf / Gaussian moments of polynomial pieces), NO Fourier
   coefficients.  Arc lengths: NIntegrate of Sqrt[1 + u_x^2], with u_x from the
   (even) extension of u0' -- Erf pieces for the hat, Gaussian image sums at
   the jumps of the indicator.  Root: FindRoot at WorkingPrecision 90.
   hd16.1 (u0 = Sin[Pi x], u known exactly) is checked by a direct NIntegrate.
   (Reference A: siam100_refs/hd16_ref.py, Fourier series.)
   Run: wolframscript -file hd16_ref.wl *)
nim = 8; wp = 90;
G[s_, t_] := Exp[-s^2/(4 t)]/Sqrt[4 Pi t];
(* Int_a^b (c0 + c1 y + c2 y^2) G(y - s, t) dy in closed form *)
mom[{c0_, c1_, c2_}, a_, b_, s_, t_] := Module[{A = a - s, B = b - s, I0, J1, J2},
  I0 = (Erf[B/(2 Sqrt[t])] - Erf[A/(2 Sqrt[t])])/2;
  J1 = Sqrt[t/Pi] (Exp[-A^2/(4 t)] - Exp[-B^2/(4 t)]);
  J2 = 2 t I0 + Sqrt[t/Pi] (A Exp[-A^2/(4 t)] - B Exp[-B^2/(4 t)]);
  c0 I0 + c1 (s I0 + J1) + c2 (s^2 I0 + 2 s J1 + J2)];
(* u(x,t) for u0 = pieces {{coef list, a, b}} on [0,L]; sgn = -1 Dirichlet, +1 Neumann *)
imgU[pieces_, L_, sgn_, x_, t_] := Sum[Total[(mom[#[[1]], #[[2]], #[[3]], x - 2 n L, t] +
       sgn mom[#[[1]], #[[2]], #[[3]], -x - 2 n L, t]) & /@ pieces], {n, -nim, nim}];
(* note Int u0(y) G(x+y-2nL) dy = Int u0(y) G(y - (-x+2nL)) dy; the n-sum is symmetric *)
root[f_, t0_] := t /. FindRoot[f[t], {t, t0, t0 (1 + 1/50)}, WorkingPrecision -> wp, AccuracyGoal -> 80, PrecisionGoal -> 80];
hat = {{{0, 2, 0}, 0, 1/2}, {{2, -2, 0}, 1/2, 1}};
show[lab_, v_] := (Print[lab, " ", N[v, 70]]; v);
show["hd16.2", root[imgU[hat, 1, -1, 1/2, #] - 1/2 &, 4918/100000]];
show["hd16.3", root[imgU[{{{1, 0, 0}, 0, 1}}, 1, -1, 1/2, #] - 1/2 &, 9468/100000]];
show["hd16.4", root[imgU[{{{0, 1, -1}, 0, 1}}, 1, -1, 1/2, #] - 1/10 &, 9603/100000]];
show["hd16.5", root[imgU[{{{1, 0, 0}, 1/4, 3/4}}, 1, -1, 1/2, #] - 1/2 &, 5988/100000]];
show["hd16.6", root[imgU[{{{0, 1, 0}, 0, 1}}, 1, +1, 0, #] - 1/3 &, 9004/100000]];
show["hd16.9", root[1/4 + imgU[{{{0, -1, 0}, 0, 1}}, 1, -1, 1/4, #] - 1/10 &, 10842/100000]];
show["hd16.10", root[imgU[{{{1, 0, 0}, 0, 1}}, 1, -1, 1/2, #]^2 - 1/2 &, 5927/100000]];
lenSin[s_?NumericQ] := NIntegrate[Sqrt[1 + (Pi Exp[-Pi^2 s] Cos[Pi xx])^2], {xx, 0, 1}, WorkingPrecision -> wp, PrecisionGoal -> 80];
show["hd16.1", root[lenSin[#] - 3/2 &, 6547/100000]];
(* u_x for the hat: u0' = 2 on (0,1/2), -2 on (1/2,1), even extension *)
uxHat[x_, s_] := Sum[mom[{2, 0, 0}, 0, 1/2, x - 2 n, s] + mom[{2, 0, 0}, 0, 1/2, -x - 2 n, s] -
     mom[{2, 0, 0}, 1/2, 1, x - 2 n, s] - mom[{2, 0, 0}, 1/2, 1, -x - 2 n, s], {n, -nim, nim}];
lenHat[s_?NumericQ] := 2 NIntegrate[Sqrt[1 + uxHat[xx, s]^2], {xx, 0, 1/4, 1/2}, WorkingPrecision -> wp, PrecisionGoal -> 80, MaxRecursion -> 30];
show["hd16.7", root[lenHat[#] - 6/5 &, 9910/100000]];
uxInd3[x_, s_] := Sum[G[x - 1 - 6 n, s] - G[x - 2 - 6 n, s] + G[x + 1 - 6 n, s] - G[x + 2 - 6 n, s], {n, -nim, nim}];
lenInd3[s_?NumericQ] := 2 NIntegrate[Sqrt[1 + uxInd3[xx, s]^2], {xx, 0, 1/2, 1, 5/4, 3/2}, WorkingPrecision -> wp, PrecisionGoal -> 80, MaxRecursion -> 30];
show["hd16.8", root[lenInd3[#] - 4 &, 1822/100000]];
