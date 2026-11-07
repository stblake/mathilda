(* D4 reference, method A: the built-in RiemannR (WL evaluates it at tiny x
   with its own high-precision series), FindRoot in t = log x at 110 digits.
   A coarse scan of u^3 R(e^-u) (u = -t) on [10^4, 7 10^4] shows the first
   sign change near u = 34150; R > 0 for all larger x (R ~ 21.9/u^3 > 0
   plus an oscillation of amplitude ~1e-10 u^-1/2 from the first zeta zero). *)
Print[Table[{u, N[u^3 RiemannR[Exp[-u]], 6]}, {u, 34000, 34300, 50}]];
ts = t /. FindRoot[RiemannR[Exp[t]] == 0, {t, -34150, -34100}, WorkingPrecision -> 110,
     AccuracyGoal -> 100, PrecisionGoal -> 100];
Print["t* = ", NumberForm[ts, 80]];
Print["x* = ", NumberForm[N[Exp[ts], 75], 70]];
