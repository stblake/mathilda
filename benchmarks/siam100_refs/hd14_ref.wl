(* Holdout A, family d14 (experiment 219), method B for hd14.7: the exact
   torsion function of the equilateral triangle of side s, u = c y (Sqrt[3] x - y)
   (Sqrt[3] (s - x) - y), with c fixed by -Laplacian u == 1 (checked to be
   constant), integrated symbolically; then Int u == 1/1000 solved for s. *)
u = c y (Sqrt[3] x - y) (Sqrt[3] (s - x) - y);
lap = Simplify[-Laplacian[u, {x, y}]];
cc = c /. First[Solve[(lap /. {x -> 0, y -> 0}) == 1, c]];
Print["-Laplacian u = ", Simplify[lap /. c -> cc]];
int = Simplify[Integrate[u /. c -> cc, {y, 0, s Sqrt[3]/2}, {x, y/Sqrt[3], s - y/Sqrt[3]}]];
Print["Int u = ", int];
Print["hd14.7 ", N[s /. First[Solve[int == 1/1000 && s > 0, s, Reals]], 70]];
(* method B for hd14.3: shooting written out in WL -- NDSolve IVP (Method ->
   "Extrapolation", WorkingPrecision 80) from u(0) = 0, u'(0) = c, an inner
   FindRoot on c for u(1) = 0 (lower branch: c near 1.86), an outer FindRoot
   on lambda for u(1/2) = 1/2. *)
ivp[l_, c_?NumericQ, x_] := NDSolveValue[{w''[t] == -l Exp[w[t]], w[0] == 0, w'[0] == c}, w[x], {t, 0, 1},
   WorkingPrecision -> 80, PrecisionGoal -> 72, AccuracyGoal -> 72, MaxSteps -> Infinity, Method -> "Extrapolation"];
slope[l_?NumericQ] := Module[{ll = SetPrecision[l, 80]},
   c /. FindRoot[ivp[ll, c, 1] == 0, {c, 186/100, 187/100}, WorkingPrecision -> 80, PrecisionGoal -> 72, AccuracyGoal -> 75]];
mid[l_?NumericQ] := ivp[SetPrecision[l, 80], slope[l], 1/2];
Print["hd14.3 ", N[l /. FindRoot[mid[l] == 1/2, {l, 263/100, 264/100}, WorkingPrecision -> 80, PrecisionGoal -> 68, AccuracyGoal -> 70], 70]];
