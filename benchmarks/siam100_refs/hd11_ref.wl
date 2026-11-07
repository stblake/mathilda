(* Holdout A, family d11 (experiment 216), method 2: NDSolve (Method ->
   "Extrapolation", WorkingPrecision 80, Precision/AccuracyGoal 70) shot in
   (theta, T) by FindRoot at WorkingPrecision 80, started on the shortest
   branch (hd11_ref.py closed forms / hd11_scan.py).  Same force laws as
   216-siam100h-d11-shooting-targets/shooting_targets.m. *)
grav[ms_] := Function[{x, y}, {Sum[(m[[1]] - x)/((m[[1]] - x)^2 + (m[[2]] - y)^2)^(3/2), {m, ms}],
                               Sum[(m[[2]] - y)/((m[[1]] - x)^2 + (m[[2]] - y)^2)^(3/2), {m, ms}]}];
cs[1] = {{0, -1} &, {0, 0}, 2, {1, 1/2}};
cs[2] = {{-#1, -#2} &, {1, 0}, 1, {0, 1/2}};
cs[3] = {{-#1, -4 #2} &, {0, 0}, 1, {1/2, 1/4}};
cs[4] = {{-#3/2, -#4/2 - 1} &, {0, 0}, 2, {1, 0}};
cs[5] = {-{#1, #2}/(#1^2 + #2^2)^(3/2) &, {1, 0}, 1, {-1/2, 1/2}};
cs[6] = {-{#1, #2}/(#1^2 + #2^2)^(3/2) &, {1, 0}, 6/5, {0, 3/2}};
cs[7] = {grav[{{2, -1}, {2, 1}}][#1, #2] &, {0, 0}, 1, {4, 0}};
cs[8] = {grav[{{1, 1}, {2, -1}, {3, 1}}][#1, #2] &, {0, 0}, 1, {4, 0}};
cs[9] = {{#4, -#3} &, {0, 0}, 1, {1, 1/2}};
cs[10] = {grav[{{2, -1}, {2, 0}, {2, 1}, {3, -1}, {3, 0}, {3, 1}}][#1, #2] &, {0, 0}, 1, {1, 2}};
odes[k_, th_] := Module[{acc = cs[k][[1]][x[t], y[t], x'[t], y'[t]], s0 = cs[k][[2]], v = cs[k][[3]]},
  {x''[t] == acc[[1]], y''[t] == acc[[2]], x[0] == s0[[1]], y[0] == s0[[2]],
   x'[0] == v Cos[th], y'[0] == v Sin[th]}];
pos[k_, th0_?NumericQ, T0_?NumericQ] := Module[{th = SetPrecision[th0, 80], T = SetPrecision[T0, 80]},
  {x[T], y[T]} /. First[NDSolve[odes[k, th], {x, y}, {t, 0, T}, Method -> "Extrapolation",
      WorkingPrecision -> 80, PrecisionGoal -> 70, AccuracyGoal -> 70, MaxSteps -> Infinity]]];
g = {{6/10, 6/10}, {25/10, 9/10}, {55/100, 63/100}, {15/100, 58/100}, {189/100, 150/100}, {154/100, 180/100}, {0, 29/10}, {-35/100, 288/100}, {106/100, 119/100}, {94/100, 203/100}};
Do[Module[{r = FindRoot[pos[k, th, T] == cs[k][[4]], {{th, g[[k, 1]]}, {T, g[[k, 2]]}},
     WorkingPrecision -> 80, AccuracyGoal -> 70, PrecisionGoal -> 70]},
  Print["hd11.", k, " theta ", N[th /. r, 25], " T ", N[T /. r, 70]]], {k, 10}];
