(* hc1 (holdout A, family C1) reference, method B: WL at WorkingPrecision 80.
   Closed forms are evaluated with N[, 75]; the others are NIntegrate with
   Method -> "ExtrapolatingOscillatory" on a substitution DIFFERENT from
   hc1_ref.py: for the t e^t phase, u = t e^t (t = W(u)), so the oscillation
   is plain sin u / cos u with a slowly varying Lambert-W amplitude. *)
opt = Sequence[Method -> "ExtrapolatingOscillatory", WorkingPrecision -> 80,
               PrecisionGoal -> 70, MaxRecursion -> 30];
w[u_] := ProductLog[u];
r[1] = NIntegrate[Sin[t]/t^2, {t, 1, Infinity}, Evaluate[opt]];
r[2] = NIntegrate[Sin[t]/t, {t, 1, Infinity}, Evaluate[opt]];
r[3] = NIntegrate[Cos[t]/(2 t), {t, 1, Infinity}, Evaluate[opt]];
r[4] = N[Integrate[Cos[t] t^(-3/2), {t, 1, Infinity}], 75];
(* 5: - int_0^oo sin(u) W/(u (1+W)) du *)
r[5] = -NIntegrate[Sin[u] w[u]/(u (1 + w[u])), {u, 0, Infinity}, Evaluate[opt]];
(* 6: e^-t dt = cos u W^2/(u^2 (1+W)) du *)
r[6] = NIntegrate[Cos[u] w[u]^2/(u^2 (1 + w[u])), {u, 0, Infinity}, Evaluate[opt]];
r[7] = NIntegrate[Cos[u]/(u Log[u]^2), {u, E, Infinity}, Evaluate[opt]];
(* 8: t = 1/(x(1-x)) on (0, 1/2]: 2 int_4^oo sin t/(t^(3/2) sqrt(t - 4)) dt,
   the endpoint square-root singularity on [4, 5] by the default rules *)
r[8] = 2 (NIntegrate[Sin[t]/(t^(3/2) Sqrt[t - 4]), {t, 4, 5}, WorkingPrecision -> 80,
                     PrecisionGoal -> 70, MaxRecursion -> 30] +
          NIntegrate[Sin[t]/(t^(3/2) Sqrt[t - 4]), {t, 5, Infinity}, Evaluate[opt]]);
r[9] = NIntegrate[Sin[t]/t, {t, 1/2, Infinity}, Evaluate[opt]];
(* 10: e^(t/2) dt = cos u sqrt(W/u)/(1+W) du *)
r[10] = NIntegrate[Cos[u] Sqrt[w[u]/u]/(1 + w[u]), {u, 0, Infinity}, Evaluate[opt]];
Do[Print[k, " ", ToString[NumberForm[r[k], 70]]], {k, 10}];
