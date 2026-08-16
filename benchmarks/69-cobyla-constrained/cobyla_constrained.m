(* Experiment 69 -- COBYLA (derivative-free constrained) optimization
   (FindMinimum, Method -> "COBYLA").

   MEASURES src/numerical_calculus/findmin.c's Powell linear-approximation
   trust-region method -- linear models of the objective and every constraint
   (central differences on a coordinate cross), a two-stage trust-region LP
   (feasibility then objective, solved by the same active-set QP as SLSQP), and
   an L-infinity exact-penalty merit -- against
   scipy.optimize.minimize(method="COBYLA"). It is the first derivative-free
   method in Mathilda to accept general (non-box) constraints.

   Both solvers are DERIVATIVE-FREE (no Jacobian to either side), so this is a
   fair derivative-free-vs-derivative-free race (like experiment 66 Powell).
   scipy's COBYLA is inequality-only, so every case here is inequality- (and
   bound-) constrained -- the common ground. The check is Round[10^3 * f] at the
   optimum: both must reach it or the join CHECK-FAILs, so the check is the
   accuracy gate. Objectives/constraints are built with Total[...] over fresh
   global symbols and passed through Evaluate[] past FindMinimum's HoldAll. *)

Get["../harness.m"];

require[{"FindMinimum"}];

(* Separable constrained quadratic: min Total[(x_i - i)^2] s.t. Total[x] <= 10.
   Unconstrained min x_i=i has Total=15>10 -> constraint active at Total=10, so
   x_i = i-1 and f* = Total[1] = n(=5). Fresh global symbols x1..xn. *)
cqVars[n_] := Table[Symbol["x" <> ToString[i]], {i, 1, n}];
cqProb[n_] := With[{v = cqVars[n]},
   {Total[Table[(v[[i]] - i)^2, {i, 1, n}]], Total[v] <= 10}];
cqSpec[n_] := Table[{Symbol["x" <> ToString[i]], 0.0}, {i, 1, n}];
cq5 = cqProb[5]; cq5s = cqSpec[5];

(* ---- 01 QP tutorial: 3 linear inequalities + bounds, f* = 0.8 at (1.4,1.7) ---- *)
qpTut = {(x - 1)^2 + (y - 2.5)^2,
         x - 2 y + 2 >= 0 && -x - 2 y + 6 >= 0 && -x + 2 y + 2 >= 0 && x >= 0 && y >= 0};
bench["01 QP tutorial", FindMinimum[Evaluate[qpTut], {{x, 2}, {y, 0}}, Method -> "COBYLA", MaxIterations -> 3000];];
check["01 QP tutorial", Round[10^3 First[FindMinimum[Evaluate[qpTut], {{x, 2}, {y, 0}}, Method -> "COBYLA", MaxIterations -> 3000]]]];

(* ---- 02 Corner LP: min x+y s.t. 3x+2y>=7, x>=0, y>=0 -> f* = 7/3 at (7/3,0) ---- *)
cornerLP = {x + y, 3 x + 2 y >= 7 && x >= 0 && y >= 0};
bench["02 Corner LP", FindMinimum[Evaluate[cornerLP], {{x, 1}, {y, 1}}, Method -> "COBYLA", MaxIterations -> 3000];];
check["02 Corner LP", Round[10^3 First[FindMinimum[Evaluate[cornerLP], {{x, 1}, {y, 1}}, Method -> "COBYLA", MaxIterations -> 3000]]]];

(* ---- 03 Nonlinear inequality: min x^2+y^2 s.t. x y >= 1 -> f* = 2 at (1,1) ---- *)
nlIneq = {x^2 + y^2, x y >= 1};
bench["03 Nonlinear-ineq xy>=1", FindMinimum[Evaluate[nlIneq], {{x, 2}, {y, 0.5}}, Method -> "COBYLA", MaxIterations -> 3000];];
check["03 Nonlinear-ineq xy>=1", Round[10^3 First[FindMinimum[Evaluate[nlIneq], {{x, 2}, {y, 0.5}}, Method -> "COBYLA", MaxIterations -> 3000]]]];

(* ---- 04 Constrained quadratic n=5, f* = 5 ---- *)
bench["04 Constrained-quad n5", FindMinimum[Evaluate[cq5], Evaluate[cq5s], Method -> "COBYLA", MaxIterations -> 3000];];
check["04 Constrained-quad n5", Round[10^3 First[FindMinimum[Evaluate[cq5], Evaluate[cq5s], Method -> "COBYLA", MaxIterations -> 3000]]]];
