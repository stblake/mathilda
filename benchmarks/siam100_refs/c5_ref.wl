(* c5 reference (method 2): the discrete minimax problem on the 128th roots of
   unity, independently of the Python active-set solve: no symmetry reduction
   (8 real unknowns: complex c_j = a_j + I b_j), solved as the second-order-cone
   problem min t s.t. |e_k| <= t by ConvexOptimization at machine precision;
   the active set is read off, then the KKT system (8 coefficients, E, and one
   multiplier per active point) is solved by FindRoot at WorkingPrecision 100.
   Convex problem: KKT with positive multipliers and all inactive |e_k| < E
   certify the global optimum. *)
M = 128; zk = Table[Exp[2 Pi I k/M], {k, 0, M - 1}];
fk = N[1/Gamma[zk], 140]; zn = N[zk, 140];
as = {a0, a1, a2, a3}; bs = {b0, b1, b2, b3}; vars = Join[as, bs];
re[k_] := Re[fk[[k]]] - Sum[as[[j + 1]] Re[zn[[k]]^j] - bs[[j + 1]] Im[zn[[k]]^j], {j, 0, 3}];
im[k_] := Im[fk[[k]]] - Sum[as[[j + 1]] Im[zn[[k]]^j] + bs[[j + 1]] Re[zn[[k]]^j], {j, 0, 3}];
e2 = Table[re[k]^2 + im[k]^2, {k, M}];
cons = Table[Norm[N[{re[k], im[k]}]] <= t, {k, M}];
sol = ConvexOptimization[t, cons, Append[vars, t]];
tm = t /. sol; Print["machine: ", tm];
(* The machine SOCP is only good to ~4 digits here (ConvexOptimization::inac),
   so the active set is the one c5_ref.py identifies from its SLSQP solve, the 9
   points k = 28,29,46,47,64,81,82,99,100 (conjugate pairs plus k = 64,
   z = -1); everything else here is independent.  With 9 real unknowns (a, b, E) the equal-modulus conditions
   |e_k| = E on those 9 points are a SQUARE system, solved by FindRoot at
   WorkingPrecision 120 with no symmetry assumption; then every inactive point
   is checked to lie strictly below E.  (The KKT multipliers are certified
   positive in c5_ref.py.) *)
act = 1 + {28, 29, 46, 47, 64, 81, 82, 99, 100};
Print["active k: ", act - 1];
r = FindRoot[Table[e2[[k]] == EE^2, {k, act}], Append[Thread[{vars, SetPrecision[vars /. sol, 120]}], {EE, SetPrecision[tm, 120]}],
   WorkingPrecision -> 120, PrecisionGoal -> 100, AccuracyGoal -> 100, MaxIterations -> 200];
Print["E_128 = ", N[EE /. r, 70]];
Print["coeffs: ", N[vars /. r, 20]];
Print["max inactive / E: ", N[Max[Sqrt[Delete[e2, List /@ act] /. r]]/(EE /. r), 10]];
