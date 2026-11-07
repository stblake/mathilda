(* SIAM100 App. D problem 15 (experiment 190), corroborating method: WL NDSolve, 4th-order FD method of lines on the v = Exp[-u] form, n = 21..61 points per side; add v/|v_t| ~ 1.1e-6 to the stop time. Results: 21 2.1192538, 31 2.1196067, 41 2.1196869, 61 2.1197202 (stop times) -> h^4 convergence to 2.11973 (6 digits). *)
(* v = Exp[-u]: v_t = Lap v - |grad v|^2/v - 1, v = 1 on boundary and at t = 0; blow-up <=> v(center) = 0 *)
tv[n_, ord_] := Module[{v, t, x, y, tt, s},
  s = NDSolveValue[{D[v[t, x, y], t] == D[v[t, x, y], x, x] + D[v[t, x, y], y, y] -
        (D[v[t, x, y], x]^2 + D[v[t, x, y], y]^2)/v[t, x, y] - 1,
     v[0, x, y] == 1, v[t, 0, y] == 1, v[t, 3, y] == 1, v[t, x, 0] == 1, v[t, x, 3] == 1,
     WhenEvent[v[t, 3/2, 3/2] == 10^-6, "StopIntegration"]},
    v, {t, 0, 3}, {x, 0, 3}, {y, 0, 3},
    Method -> {"MethodOfLines", "SpatialDiscretization" -> {"TensorProductGrid",
        "MinPoints" -> n, "MaxPoints" -> n, "DifferenceOrder" -> ord}},
    PrecisionGoal -> 12, AccuracyGoal -> 14, MaxSteps -> Infinity];
  tt = s["Domain"][[1, 2]];
  {tt, s[tt, 3/2, 3/2]}];
Do[Print[n, " ", AbsoluteTiming[tv[n, 4]]], {n, {21, 31, 41, 61}}];
