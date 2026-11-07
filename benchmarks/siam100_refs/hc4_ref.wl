(* hc4 (holdout A, family C4) reference, method B: WL global search
   (NMinimize DifferentialEvolution, 300 search points, three seeds) then a
   Newton polish (FindMinimum at WorkingPrecision 120, Abs[u] -> Sqrt[u^2],
   differentiable at the minimiser) -- for k = 8 the minimiser is on the box
   edge x = 512, so the polish is in y with x = 512.  Method A: hc4_ref.py
   (exhaustive grid + mpmath Newton); closed forms for k = 2, 3. *)
f[1] = (4 - 21/10 x^2 + x^4/3) x^2 + x y + (-4 + 4 y^2) y^2;            b[1] = {{-3, 3}, {-2, 2}};
f[2] = (y - 51/10/(4 Pi^2) x^2 + 5/Pi x - 6)^2 + 10 (1 - 1/(8 Pi)) Cos[x] + 10; b[2] = {{-5, 10}, {0, 15}};
f[3] = Sum[(u^4 - 16 u^2 + 5 u)/2, {u, {x, y, z}}];                     b[3] = {{-5, 5}, {-5, 5}, {-5, 5}};
f[4] = Sum[i Cos[(i + 1) x + i], {i, 5}] Sum[i Cos[(i + 1) y + i], {i, 5}]; b[4] = {{-10, 10}, {-10, 10}};
f[5] = -Abs[Sin[x] Cos[y] Exp[Abs[1 - Sqrt[x^2 + y^2]/Pi]]];            b[5] = {{-10, 10}, {-10, 10}};
f[6] = -(Abs[Sin[x] Sin[y] Exp[Abs[100 - Sqrt[x^2 + y^2]/Pi]]] + 1)^(1/10)/10000; b[6] = {{-10, 10}, {-10, 10}};
f[7] = Exp[Sin[30 x]] + Sin[40 Exp[y]] + Sin[50 Sin[x]] + Sin[Sin[60 y]] - Sin[10 (x + y)] + (x^2 + y^2)/4; b[7] = {{-1, 1}, {-1, 1}};
f[8] = -(y + 47) Sin[Sqrt[Abs[x/2 + y + 47]]] - x Sin[Sqrt[Abs[x - (y + 47)]]]; b[8] = {{-512, 512}, {-512, 512}};
f[9] = (x^2 + y - 11)^2 + (x + y^2 - 7)^2 + x + y;                       b[9] = {{-5, 5}, {-5, 5}};
f[10] = Exp[Sin[70 x]] + Sin[80 Exp[y]] + Sin[90 Sin[x]] + Sin[Sin[100 y]] - Sin[12 (x + y)] + (x^2 + y^2)/4; b[10] = {{-1, 1}, {-1, 1}};
Do[
  vars = If[k == 3, {x, y, z}, {x, y}];
  cons = And @@ MapThread[#2[[1]] <= #1 <= #2[[2]] &, {vars, b[k]}];
  runs = Table[NMinimize[{f[k], cons}, vars, Method -> {"DifferentialEvolution", "SearchPoints" -> 300, "RandomSeed" -> s}], {s, 3}];
  best = First[SortBy[runs, First]];
  fs = f[k] /. Abs[u_] :> Sqrt[u^2];
  If[k == 8,
    pol = FindMinimum[fs /. x -> 512, {y, SetPrecision[y /. best[[2]], 120]}, WorkingPrecision -> 120,
                      AccuracyGoal -> 100, PrecisionGoal -> 100, MaxIterations -> 500],
    pol = FindMinimum[fs, Transpose[{vars, SetPrecision[vars /. best[[2]], 120]}], Method -> "Newton",
                      WorkingPrecision -> 120, AccuracyGoal -> 100, PrecisionGoal -> 100, MaxIterations -> 500]];
  Print[k, " ", ToString[NumberForm[First[pol], 70]], "   DE best ", First[best], " at ", N[vars /. best[[2]], 12]],
  {k, 10}];
Print["closed 2: ", ToString[NumberForm[N[5/(4 Pi), 70], 70]]];
Print["closed 3: ", ToString[NumberForm[N[3 ((t^4 - 16 t^2 + 5 t)/2 /. t -> Root[4 #^3 - 32 # + 5 &, 1]), 70], 70]]];
