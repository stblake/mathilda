(* Reference for d11 (six masses), method 1: WL NDSolve (Method -> "Extrapolation";
   the default method fails at WorkingPrecision 60-80 with NDSolve::nderr / ndtol
   "error test failure" at t ~ 0.3-0.7,
   WorkingPrecision 80, Precision/AccuracyGoal 70) shot in (theta, T) by
   FindRoot at WorkingPrecision 80 (arguments lifted to 80 digits by
   SetPrecision so FindRoot's difference steps do not degrade the ODE), started on the shortest branch found by
   siam100_refs/d11_scan.py (theta ~ -0.274881254, T ~ 1.9160969183). *)
ms = {{2, -1}, {2, 0}, {2, 1}, {3, -1}, {3, 0}, {3, 1}};
eqs[th_] := Join[
   Thread[{x''[t], y''[t]} == Total[(# - {x[t], y[t]})/((#[[1]] - x[t])^2 + (#[[2]] - y[t])^2)^(3/2) & /@ ms]],
   {x[0] == 0, y[0] == 0, x'[0] == Cos[th], y'[0] == Sin[th]}];
pos[th0_?NumericQ, T0_?NumericQ] := Module[{th = SetPrecision[th0, 80], T = SetPrecision[T0, 80]}, {x[T], y[T]} /. First[NDSolve[eqs[th], {x, y}, {t, 0, T}, Method -> "Extrapolation",
      WorkingPrecision -> 80, PrecisionGoal -> 70, AccuracyGoal -> 70, MaxSteps -> Infinity]]];
sol = FindRoot[pos[th, T] == {4, 1}, {{th, -274881254/10^9}, {T, 19160969183/10^10}},
   WorkingPrecision -> 80, AccuracyGoal -> 66, PrecisionGoal -> 66];
Print[N[{th, T} /. sol, 70]];
Print["residual ", N[(pos @@ ({th, T} /. sol)) - {4, 1}, 5]];
