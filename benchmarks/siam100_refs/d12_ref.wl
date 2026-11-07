(* Reference for SIAM100 App. D problem 12 (experiment 187), method B:
   limiting ODE on the unit circle (see d12_ref.py for the derivation),
   NDSolve at WorkingPrecision 80 with WhenEvent for the first hit of y = 0,
   FindRoot (secant) on alpha. *)
T[a_?NumericQ] := Block[{th, t},
  NDSolveValue[{th''[t] == -Cos[th[t]] + a Sin[th[t]] Cos[th[t]]/(2 Sqrt[1 + a Cos[th[t]]^2]),
    th[0] == Pi/2, th'[0] == -1, WhenEvent[Sin[th[t]] == 0, "StopIntegration"]},
    th, {t, 0, 20}, WorkingPrecision -> 80, PrecisionGoal -> 70, AccuracyGoal -> 70,
    MaxSteps -> Infinity, Method -> "Extrapolation"]["Domain"][[1, 2]]];
r = FindRoot[T[a] == 10, {a, 106234/10000, 10623/1000}, WorkingPrecision -> 80, PrecisionGoal -> 65, AccuracyGoal -> 65];
Print[N[a /. r, 70]];
