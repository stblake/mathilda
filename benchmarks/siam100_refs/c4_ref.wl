(* c4 reference (method 2): FindMinimum (Newton) at WorkingPrecision 120 from
   a seed that NMinimize[..., Method -> "DifferentialEvolution"] finds over the
   box [-1,1]^2 (the minimum must lie there: (x^2+y^2)/4 > 1 + e + 3 outside). *)
f = Exp[Sin[50 x]] + Sin[60 Exp[y]] + Sin[70 Sin[x]] + Sin[Sin[80 y]] - Sin[10 (x + y)] + (x^2 + y^2)/4;
seed = Last[NMinimize[{f, -1 <= x <= 1, -1 <= y <= 1}, {x, y}, Method -> {"DifferentialEvolution", "SearchPoints" -> 400}]];
r = FindMinimum[f, {{x, x /. seed}, {y, y /. seed}}, Method -> "Newton", WorkingPrecision -> 120, PrecisionGoal -> 100, AccuracyGoal -> 100];
Print[seed]; Print[N[First[r], 70]];
