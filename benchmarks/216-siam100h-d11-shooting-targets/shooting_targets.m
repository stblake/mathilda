(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 216 -- SIAM 100-Digit Challenge, holdout A, family d11
   (generalises App. D Problem 11, experiment 186: "shortest time in which a
   unit-speed particle launched from (0,0) among six fixed masses reaches
   (4,1)").

   The capability: TWO-UNKNOWN SHOOTING -- FindRoot in (launch angle theta,
   arrival time T) on {x(T), y(T)} == target, with an NDSolve inside; the
   answer is T on the SHORTEST branch.  Ten instances that vary the force law
   in kind (uniform field, isotropic / anisotropic springs, linear drag,
   Kepler centre, magnetic (velocity-dependent) force, sets of fixed masses):

     hd11.1  r'' = (0,-1), speed 2 from (0,0), target (1, 1/2)       T = Sqrt[7 - 2 Sqrt[11]]
     hd11.2  r'' = -r, speed 1 from (1,0), target (0, 1/2)           Cos[2T] = -1/4
     hd11.3  x'' = -x, y'' = -4y, speed 1 from (0,0), target (1/2, 1/4)
     hd11.4  r'' = -r'/2 + (0,-1), speed 2 from (0,0), target (1, 0)  (linear drag)
     hd11.5  r'' = -r/|r|^3, speed 1 from (1,0), target (-1/2, 1/2)   (Kepler)
     hd11.6  r'' = -r/|r|^3, speed 6/5 from (1,0), target (0, 3/2)    (Kepler)
     hd11.7  unit masses at (2, +-1), speed 1 from (0,0), target (4, 0)
     hd11.8  unit masses at (1,1), (2,-1), (3,1), speed 1 from (0,0), target (4, 0)
     hd11.9  r'' = (y', -x') (unit magnetic field), speed 1 from (0,0), target (1, 1/2)
     hd11.10 exp 186's six masses, speed 1 from (0,0), target (1, 2)
   Masses attract with r^-2 as in exp 186.  Every case starts the root
   finder at the same rounded guess on the shortest branch (branches from
   closed forms, or siam100_refs/hd11_scan.py) and is scored on its T.

   Asked of each tier (the statement of exp 186):
     machine: NDSolve PrecisionGoal/AccuracyGoal 13, FindRoot defaults;
     wp50:    NDSolve and FindRoot at WorkingPrecision -> 60 (NDSolve
              PrecisionGoal/AccuracyGoal 55), exact starting guesses.
   True values and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NDSolve", "FindRoot"}];

grav[ms_] := Function[{x, y}, {Sum[(m[[1]] - x)/((m[[1]] - x)^2 + (m[[2]] - y)^2)^(3/2), {m, ms}],
                               Sum[(m[[2]] - y)/((m[[1]] - x)^2 + (m[[2]] - y)^2)^(3/2), {m, ms}]}];
(* k -> {acceleration as a function of {x, y, x', y'}, start, speed, target, {theta0, T0}} *)
cs[1] = {{0, -1} &, {0, 0}, 2, {1, 1/2}, {6/10, 6/10}};
cs[2] = {{-#1, -#2} &, {1, 0}, 1, {0, 1/2}, {25/10, 9/10}};
cs[3] = {{-#1, -4 #2} &, {0, 0}, 1, {1/2, 1/4}, {55/100, 63/100}};
cs[4] = {{-#3/2, -#4/2 - 1} &, {0, 0}, 2, {1, 0}, {15/100, 58/100}};
cs[5] = {-{#1, #2}/(#1^2 + #2^2)^(3/2) &, {1, 0}, 1, {-1/2, 1/2}, {189/100, 150/100}};
cs[6] = {-{#1, #2}/(#1^2 + #2^2)^(3/2) &, {1, 0}, 6/5, {0, 3/2}, {154/100, 180/100}};
cs[7] = {grav[{{2, -1}, {2, 1}}][#1, #2] &, {0, 0}, 1, {4, 0}, {5/100, 29/10}};
cs[8] = {grav[{{1, 1}, {2, -1}, {3, 1}}][#1, #2] &, {0, 0}, 1, {4, 0}, {-35/100, 288/100}};
cs[9] = {{#4, -#3} &, {0, 0}, 1, {1, 1/2}, {106/100, 119/100}};
cs[10] = {grav[{{2, -1}, {2, 0}, {2, 1}, {3, -1}, {3, 0}, {3, 1}}][#1, #2] &, {0, 0}, 1, {1, 2}, {94/100, 203/100}};

odes[k_, th_] := Module[{acc = cs[k][[1]][x[t], y[t], x'[t], y'[t]], s0 = cs[k][[2]], v = cs[k][[3]]},
  {x''[t] == acc[[1]], y''[t] == acc[[2]], x[0] == s0[[1]], y[0] == s0[[2]],
   x'[0] == v Cos[th], y'[0] == v Sin[th]}];
posM[k_, th_?NumericQ, T_?NumericQ] := {x[T], y[T]} /. First[NDSolve[odes[k, th], {x, y}, {t, 0, T},
    PrecisionGoal -> 13, AccuracyGoal -> 13, MaxSteps -> Infinity]];
pos50[k_, th_?NumericQ, T_?NumericQ] := {x[T], y[T]} /. First[NDSolve[odes[k, th], {x, y}, {t, 0, T},
    WorkingPrecision -> 60, PrecisionGoal -> 55, AccuracyGoal -> 55, MaxSteps -> Infinity]];

Do[With[{k = k, tg = cs[k][[4]], g = cs[k][[5]]},
  vrun["hd11." <> ToString[k] <> " machine",
       T /. FindRoot[posM[k, th, T] == tg, {{th, N[g[[1]]]}, {T, N[g[[2]]]}}], 10, 1, 60];
  vrun["hd11." <> ToString[k] <> " wp50",
       T /. FindRoot[pos50[k, th, T] == tg, {{th, g[[1]]}, {T, g[[2]]}}, WorkingPrecision -> 60], 48, 1, 90]],
  {k, 10}];
