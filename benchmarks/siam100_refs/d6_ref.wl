(* D6 reference: NIntegrate, global adaptive Gauss-Kronrod, WorkingPrecision 90 *)
v = NIntegrate[Exp[x + Sin[Exp[Exp[x + 1/3]]]], {x, -1, 1}, WorkingPrecision -> 90, PrecisionGoal -> 70, MaxRecursion -> 50, Method -> "GaussKronrodRule"];
Print[NumberForm[v, 70]];
