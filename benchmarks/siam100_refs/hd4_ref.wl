(* Holdout A, family hd4, method A: WL's own special functions at WP 90
   (independent of mpmath's riemannr/li/zeta/siegelz and of the zeta-zero
   expansion used for problem 10 in hd4_ref.py).  Run: wolframscript -file hd4_ref.wl *)
wp = 90;
r[1] = x /. FindRoot[LogIntegral[x] == 0, {x, 3/2}, WorkingPrecision -> wp, AccuracyGoal -> 80, PrecisionGoal -> 80];
r[2] = x /. FindRoot[RiemannR[x] == 100, {x, 540}, WorkingPrecision -> wp, AccuracyGoal -> 80, PrecisionGoal -> 80];
r[3] = x /. FindRoot[RiemannR[x] == 2, {x, 3}, WorkingPrecision -> wp, AccuracyGoal -> 80, PrecisionGoal -> 80];
r[4] = x /. FindRoot[RiemannR[x] == 1000, {x, 7900}, WorkingPrecision -> wp, AccuracyGoal -> 80, PrecisionGoal -> 80];
r[5] = x /. FindRoot[RiemannR[x] == 1/2, {x, 3/10}, WorkingPrecision -> wp, AccuracyGoal -> 80, PrecisionGoal -> 80];
r[6] = s /. FindRoot[Zeta[s] == 2, {s, 17/10}, WorkingPrecision -> wp, AccuracyGoal -> 80, PrecisionGoal -> 80];
r[7] = s /. FindRoot[-Zeta'[s]/Zeta[s] == 1, {s, 17/10}, WorkingPrecision -> wp, AccuracyGoal -> 80, PrecisionGoal -> 80];
r[8] = Im[N[ZetaZero[1], wp]];
r[9] = x /. FindRoot[RiemannR[x] == 10^10, {x, 25 10^10}, WorkingPrecision -> wp, AccuracyGoal -> 80, PrecisionGoal -> 80];
t10 = t /. FindRoot[RiemannR[Exp[t]] == 0, {t, -35240, -35220}, WorkingPrecision -> wp, AccuracyGoal -> 80, PrecisionGoal -> 80];
r[10] = Exp[t10];
Print["Soldner check: ", N[r[1] - SoldnerConstant, 5]];
Print["scan u^3 R(e^-u) between the first two zeros: ",
  Table[{u, N[u^3 RiemannR[Exp[-u]], 4]}, {u, 34000, 35400, 200}]];
Do[Print[k, " ", NumberForm[r[k], 75]], {k, 1, 10}];
Print["10 t2 = ", NumberForm[t10, 80]];
