(* Holdout A, family d12 (experiment 217), method 2: NDSolve (Method ->
   "Extrapolation", WorkingPrecision 80) with WhenEvent for the first event,
   FindRoot (secant) on the parameter a at WorkingPrecision 80.  Method 1 is
   hd12_ref.py (closed forms / quadrature / mpmath odefun). *)
hit[eqs_, ic_, ev_, dep_, a0_?NumericQ, tmax_] := Block[{a = SetPrecision[a0, 80]},
  NDSolveValue[Join[eqs /. aa -> a, ic, {WhenEvent[ev, "StopIntegration"]}], dep, {t, 0, tmax},
    WorkingPrecision -> 80, PrecisionGoal -> 70, AccuracyGoal -> 70, MaxSteps -> Infinity,
    Method -> "Extrapolation"]["Domain"][[1, 2]]];
solve[k_, eqs_, ic_, ev_, dep_, T_, g_, tmax_] := Module[{r},
  r = FindRoot[hit[eqs, ic, ev, dep, a, tmax] == T, {a, g, g*(1 + 10^-6)}, WorkingPrecision -> 80,
    PrecisionGoal -> 66, AccuracyGoal -> 66];
  Print["hd12.", k, " ", N[a /. r, 70]]];
solve[1, {y''[t] == -1 - aa y'[t]^2}, {y[0] == 0, y'[0] == 2}, y'[t] == 0, y, 1, 136/100, 10];
solve[2, {th''[t] == -aa Sin[th[t]]}, {th[0] == 0, th'[0] == 1}, th'[t] == 0, th, 2, 75/100, 10];
solve[3, {x''[t] == -x[t] - aa x[t]^3}, {x[0] == 1, x'[0] == 0}, x[t] == 0, x, 6/5, 97/100, 10];
solve[4, {x''[t] + aa x'[t] + x[t] == 0}, {x[0] == 1, x'[0] == 0}, x[t] == 0, x, 2, 64/100, 10];
solve[5, {y'[t] == aa + y[t]^2}, {y[0] == 0}, y[t] == 1, y, 1/2, 17/10, 10];
solve[6, {r''[t] == -aa/r[t]^2}, {r[0] == 1, r'[0] == 0}, r[t] == 1/2, r, 1, 83/100, 10];
solve[7, {y''[t] + aa t y[t] == 0}, {y[0] == 1, y'[0] == 0}, y[t] == 0, y, 1, 78/10, 10];
solve[8, {y''[t] + (aa + Cos[t]) y[t] == 0}, {y[0] == 1, y'[0] == 0}, y[t] == 0, y, 1, 153/100, 10];
solve[9, {th''[t] == -Cos[th[t]] + aa Sin[th[t]] Cos[th[t]]/(2 Sqrt[1 + aa Cos[th[t]]^2])},
  {th[0] == Pi/2, th'[0] == -1}, Sin[th[t]] == 0, th, 5, 10574/1000, 20];
solve[10, {x''[t] - aa (1 - x[t]^2) x'[t] + x[t] == 0}, {x[0] == 2, x'[0] == 0}, x[t] == 0, x, 2, 75/100, 10];
