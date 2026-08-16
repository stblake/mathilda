(* Experiment 65 -- L-BFGS-B accuracy + efficiency stress (FindMinimum,
   Method -> "LBFGSB").

   MEASURES src/numerical_calculus/findmin.c's limited-memory BFGS on 10 hard
   local-optimization test functions with KNOWN analytic optima, against
   scipy.optimize.minimize(method="L-BFGS-B"). Each case stresses accuracy (a
   curved valley, ill conditioning up to 1e6, a singular Hessian, an active
   bound, a coupled Trid quadratic) and is timed for efficiency. scipy is given
   the analytic Jacobian (jac=), so the race is solver-vs-solver.

   The check is Round[10^6 * f] at the optimum (0 for the f*=0 problems), or the
   exact integer optimum for Trid (-50) and the bound corner (5): both systems
   must reach the true optimum to that precision or the join CHECK-FAILs, so the
   check itself is the accuracy gate. bench times the solve (efficiency). Raw
   achieved errors are reported separately by the accuracy probe. Objectives are
   built with Total[Table[...]] (Sum with a symbolic iterator burns the budget on
   a closed-form attempt). *)

Get["../harness.m"];

require[{"FindMinimum"}];

(* n-variable objectives / specs built out of fresh scalar symbols z1..zn. *)
illObj[n_, c_] := Total[Table[10^(c (i - 1)/(n - 1)) Symbol["z" <> ToString[i]]^2, {i, 1, n}]];
illSpec[n_] := Table[{Symbol["z" <> ToString[i]], 1.0}, {i, 1, n}];
tridObj[n_] := Total[Table[(Symbol["z" <> ToString[i]] - 1)^2, {i, 1, n}]]
             - Total[Table[Symbol["z" <> ToString[i]] Symbol["z" <> ToString[i - 1]], {i, 2, n}]];
tridSpec[n_] := Table[{Symbol["z" <> ToString[i]], 0.0}, {i, 1, n}];

(* Precompute the n-var problems out of the timed region. *)
ill10 = illObj[10, 6];   ill10s = illSpec[10];
ill50 = illObj[50, 4];   ill50s = illSpec[50];
trid6 = tridObj[6];      trid6s = tridSpec[6];

wood = 100 (y - x^2)^2 + (1 - x)^2 + 90 (w - z^2)^2 + (1 - z)^2
     + 10.1 ((y - 1)^2 + (w - 1)^2) + 19.8 (y - 1) (w - 1);
powell = (x + 10 y)^2 + 5 (z - w)^2 + (y - 2 z)^4 + 10 (x - w)^4;
beale = (1.5 - x + x y)^2 + (2.25 - x + x y^2)^2 + (2.625 - x + x y^3)^2;

(* ---- 01 Rosenbrock 2D (curved valley), f* = 0 ---- *)
bench["01 Rosenbrock2D", FindMinimum[(1 - x)^2 + 100 (y - x^2)^2, {{x, -1.2}, {y, 1}}, Method -> "LBFGSB", MaxIterations -> 2000];];
check["01 Rosenbrock2D", Round[10^6 First[FindMinimum[(1 - x)^2 + 100 (y - x^2)^2, {{x, -1.2}, {y, 1}}, Method -> "LBFGSB", MaxIterations -> 2000]]]];

(* ---- 02 ill-conditioned quadratic n=10, condition 1e6, f* = 0 ---- *)
bench["02 illcond n10 c1e6", FindMinimum[Evaluate[ill10], Evaluate[ill10s], Method -> "LBFGSB", MaxIterations -> 5000];];
check["02 illcond n10 c1e6", Round[10^6 First[FindMinimum[Evaluate[ill10], Evaluate[ill10s], Method -> "LBFGSB", MaxIterations -> 5000]]]];

(* ---- 03 Beale, f* = 0 at (3, 0.5) ---- *)
bench["03 Beale", FindMinimum[Evaluate[beale], {{x, 1}, {y, 1}}, Method -> "LBFGSB", MaxIterations -> 2000];];
check["03 Beale", Round[10^6 First[FindMinimum[Evaluate[beale], {{x, 1}, {y, 1}}, Method -> "LBFGSB", MaxIterations -> 2000]]]];

(* ---- 04 Booth, f* = 0 at (1, 3) ---- *)
bench["04 Booth", FindMinimum[(x + 2 y - 7)^2 + (2 x + y - 5)^2, {{x, 0}, {y, 0}}, Method -> "LBFGSB"];];
check["04 Booth", Round[10^6 First[FindMinimum[(x + 2 y - 7)^2 + (2 x + y - 5)^2, {{x, 0}, {y, 0}}, Method -> "LBFGSB"]]]];

(* ---- 05 Matyas, f* = 0 at origin ---- *)
bench["05 Matyas", FindMinimum[0.26 (x^2 + y^2) - 0.48 x y, {{x, 3}, {y, 3}}, Method -> "LBFGSB"];];
check["05 Matyas", Round[10^6 First[FindMinimum[0.26 (x^2 + y^2) - 0.48 x y, {{x, 3}, {y, 3}}, Method -> "LBFGSB"]]]];

(* ---- 06 Wood/Colville 4D, f* = 0 ---- *)
bench["06 Wood", FindMinimum[Evaluate[wood], {{x, -3}, {y, -1}, {z, -3}, {w, -1}}, Method -> "LBFGSB", MaxIterations -> 3000];];
check["06 Wood", Round[10^6 First[FindMinimum[Evaluate[wood], {{x, -3}, {y, -1}, {z, -3}, {w, -1}}, Method -> "LBFGSB", MaxIterations -> 3000]]]];

(* ---- 07 Powell singular 4D (singular Hessian at the optimum), f* = 0 ---- *)
bench["07 Powell-singular", FindMinimum[Evaluate[powell], {{x, 3}, {y, -1}, {z, 0}, {w, 1}}, Method -> "LBFGSB", MaxIterations -> 3000];];
check["07 Powell-singular", Round[10^6 First[FindMinimum[Evaluate[powell], {{x, 3}, {y, -1}, {z, 0}, {w, 1}}, Method -> "LBFGSB", MaxIterations -> 3000]]]];

(* ---- 08 Trid n=6 (coupled quadratic), f* = -50 ---- *)
bench["08 Trid n6", FindMinimum[Evaluate[trid6], Evaluate[trid6s], Method -> "LBFGSB", MaxIterations -> 2000];];
check["08 Trid n6", Round[First[FindMinimum[Evaluate[trid6], Evaluate[trid6s], Method -> "LBFGSB", MaxIterations -> 2000]]]];

(* ---- 09 bound-active quadratic on [0,1]^2 (optimum on the corner), f* = 5 ---- *)
bench["09 Bound-corner", FindMinimum[(x - 2)^2 + (y - 3)^2, {{x, 0, 0, 1}, {y, 0, 0, 1}}, Method -> "LBFGSB"];];
check["09 Bound-corner", Round[First[FindMinimum[(x - 2)^2 + (y - 3)^2, {{x, 0, 0, 1}, {y, 0, 0, 1}}, Method -> "LBFGSB"]]]];

(* ---- 10 ill-conditioned quadratic n=50, condition 1e4, f* = 0 ---- *)
bench["10 illcond n50 c1e4", FindMinimum[Evaluate[ill50], Evaluate[ill50s], Method -> "LBFGSB", MaxIterations -> 3000];];
check["10 illcond n50 c1e4", Round[10^6 First[FindMinimum[Evaluate[ill50], Evaluate[ill50s], Method -> "LBFGSB", MaxIterations -> 3000]]]];
