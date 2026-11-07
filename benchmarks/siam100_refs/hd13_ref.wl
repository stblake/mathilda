(* Holdout A, family d13 (experiment 218), method 2: NDSolve (Method ->
   "Extrapolation", WorkingPrecision 80) for the state at the final time,
   FindRoot (secant) on the initial value at WorkingPrecision 80.  Method 1
   is hd13_ref.py (closed forms; odefun for hd13.6/hd13.9). *)
op = {WorkingPrecision -> 80, PrecisionGoal -> 72, AccuracyGoal -> 72, MaxSteps -> Infinity, Method -> "Extrapolation"};
fr[k_, g_, x0_] := Print["hd13.", k, " ", N[p /. FindRoot[g[p], {p, x0, x0 (1 + 10^-8)}, WorkingPrecision -> 80,
      PrecisionGoal -> 68, AccuracyGoal -> 75], 70]];
g1[p_?NumericQ] := NDSolveValue[{y'[t] == y[t] (1 - y[t]), y[0] == SetPrecision[p, 80]}, y[10], {t, 0, 10}, Evaluate[Sequence @@ op]] - 1/2;
g2[p_?NumericQ] := Norm[NDSolveValue[{x'[t] == -x[t] + 10 y[t], y'[t] == -y[t] + 10 z[t], z'[t] == -z[t],
     x[0] == 0, y[0] == 0, z[0] == SetPrecision[p, 80]}, {x[5], y[5], z[5]}, {t, 0, 5}, Evaluate[Sequence @@ op]]] - 1;
g3[p_?NumericQ] := NDSolveValue[{y'[t] == y[t]^2 - y[t], y[0] == SetPrecision[p, 80]}, y[3], {t, 0, 3}, Evaluate[Sequence @@ op]] - 2;
g4[p_?NumericQ] := NDSolveValue[{th''[t] == -Sin[th[t]], th[0] == 0, th'[0] == SetPrecision[p, 80]}, th[2], {t, 0, 2}, Evaluate[Sequence @@ op]] - 1;
g5[p_?NumericQ] := NDSolveValue[{y''[t] == -Exp[y[t]], y[0] == 0, y'[0] == SetPrecision[p, 80]}, y[1], {t, 0, 1}, Evaluate[Sequence @@ op]];
g6[p_?NumericQ] := NDSolveValue[{f'''[t] + f[t] f''[t]/2 == 0, f[0] == 0, f'[0] == 0, f''[0] == SetPrecision[p, 80]}, f'[10], {t, 0, 10}, Evaluate[Sequence @@ op]] - 1;
g7[p_?NumericQ] := Norm[NDSolveValue[{x''[t] == -x[t]/(x[t]^2 + y[t]^2)^(3/2), y''[t] == -y[t]/(x[t]^2 + y[t]^2)^(3/2),
     x[0] == 1, y[0] == 0, x'[0] == 0, y'[0] == SetPrecision[p, 80]}, {x[2], y[2]}, {t, 0, 2}, Evaluate[Sequence @@ op]]] - 3/2;
g8[p_?NumericQ] := NDSolveValue[{y'[t] == t^2 + y[t]^2, y[0] == SetPrecision[p, 80]}, y[1], {t, 0, 1}, Evaluate[Sequence @@ op]] - 2;
g9[p_?NumericQ] := Block[{n = Sqrt[x[t]^2 + y[t]^2 + z[t]^2]}, Norm[NDSolveValue[{
     x'[t] == -x[t] + 10 y[t] + n (-7/10 y[t] - 3/100 z[t]), y'[t] == -y[t] + 10 z[t] + n (7/10 y[t] - 1/10 z[t]),
     z'[t] == -z[t] + n (3/100 x[t] + 1/10 y[t]), x[0] == 0, y[0] == 0, z[0] == SetPrecision[p, 80]},
     {x[30], y[30], z[30]}, {t, 0, 30}, Evaluate[Sequence @@ op]]] - 1];
g10[p_?NumericQ] := NDSolveValue[{x''[t] + x[t] + x[t]^3 == 0, x[0] == SetPrecision[p, 80], x'[0] == 0}, x[10], {t, 0, 10}, Evaluate[Sequence @@ op]];
fr[1, g1, 454/10^7]; fr[2, g2, 1186/10^4]; fr[3, g3, 10255/10^4]; fr[4, g4, 996/1000]; fr[5, g5, 549/1000];
fr[6, g6, 332/1000]; fr[7, g7, 1184/1000]; fr[8, g8, 570/1000]; fr[9, g9, 6483984378704/10^15]; fr[10, g10, 530/1000];
