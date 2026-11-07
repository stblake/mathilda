(* Reference for SIAM100 App. D problem 13 (experiment 188), method B:
   NDSolve at WorkingPrecision 90 (Taylor-like "Extrapolation" method),
   FindRoot (secant) on z0 for Norm[u(50)] == 1.  Coefficients exact
   (7/10, 3/100, 1/10). *)
nrm[z0_?NumericQ] := Block[{x, y, z, t, s, n},
  n = Sqrt[x[t]^2 + y[t]^2 + z[t]^2];
  s = NDSolveValue[{x'[t] == -x[t] + 10 y[t] + n (-7/10 y[t] - 3/100 z[t]),
     y'[t] == -y[t] + 10 z[t] + n (7/10 y[t] - 1/10 z[t]),
     z'[t] == -z[t] + n (3/100 x[t] + 1/10 y[t]),
     x[0] == 0, y[0] == 0, z[0] == z0}, {x[50], y[50], z[50]}, {t, 0, 50},
    WorkingPrecision -> 90, PrecisionGoal -> 80, AccuracyGoal -> 80,
    MaxSteps -> Infinity, Method -> "Extrapolation"];
  Norm[s]];
r = FindRoot[nrm[z0] == 1, {z0, 6482396179622983/10^18, 6482396179622990/10^18},
  WorkingPrecision -> 90, PrecisionGoal -> 70, AccuracyGoal -> 80];
Print[N[z0 /. r, 75]];
