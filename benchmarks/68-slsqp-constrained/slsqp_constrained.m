(* Experiment 68 -- SLSQP (sequential least-squares QP) constrained
   optimization (FindMinimum, Method -> "SLSQP"; alias
   "SequentialQuadraticProgramming").

   MEASURES src/numerical_calculus/findmin.c's Han-Powell SQP -- each outer step
   solves a QP (quadratic model of f with a BFGS Lagrangian-Hessian kept SPD by
   Powell damping; linearized constraints) via a Goldfarb-Idnani dual active-set
   solver, accepted by an L1 exact-penalty merit line search -- against
   scipy.optimize.minimize(method="SLSQP"). Unlike the other gradient methods,
   SLSQP consumes general (non-box) constraints DIRECTLY through the QP rather
   than the augmented-Lagrangian penalty wrapper, so it reaches the true
   constrained optimum with accurate constraint satisfaction.

   Like experiments 65/67, scipy is given analytic Jacobians (objective AND
   constraints), so the race is solver-vs-solver -- both are gradient-based SQP.
   The cases span the constraint kinds SLSQP handles: a linear-inequality QP
   (the scipy SLSQP tutorial), HS71 (equality + inequality + bounds), a
   nonlinear equality, and equal-weight simplex QPs at n=10/50/100 (scaling of
   an equality + nonnegativity feasible set). The check is the objective at the
   optimum, rounded (10^6, or 10^3 for HS71 whose optimum is irrational): both
   systems must reach it or the join CHECK-FAILs, so the check is the accuracy
   gate. Objectives/constraints are built with Total[...] over fresh global
   symbols and passed through Evaluate[] past FindMinimum's HoldAll. *)

Get["../harness.m"];

require[{"FindMinimum"}];

(* Equal-weight simplex QP: min Total[w^2] s.t. Total[w]==1, w_i>=0 -> 1/n at
   the uniform point. Fresh global symbols w1..wn so obj/cons/spec share them. *)
splxVars[n_] := Table[Symbol["w" <> ToString[i]], {i, 1, n}];
splxProb[n_] := With[{w = splxVars[n]}, {Total[w^2], Total[w] == 1 && And @@ (# >= 0 & /@ w)}];
splxSpec[n_] := Table[{Symbol["w" <> ToString[i]], 1./n}, {i, 1, n}];

(* Precompute out of the timed region. *)
splx10 = splxProb[10];   splx10s = splxSpec[10];
splx50 = splxProb[50];   splx50s = splxSpec[50];
splx100 = splxProb[100]; splx100s = splxSpec[100];

(* ---- 01 QP tutorial: 5 linear inequalities + bounds, f* = 0.8 at (1.4,1.7).
        The scipy SLSQP docs example. ---- *)
qpTut = {(x - 1)^2 + (y - 2.5)^2,
         x - 2 y + 2 >= 0 && -x - 2 y + 6 >= 0 && -x + 2 y + 2 >= 0 && x >= 0 && y >= 0};
bench["01 QP tutorial", FindMinimum[Evaluate[qpTut], {{x, 2}, {y, 0}}, Method -> "SLSQP"];];
check["01 QP tutorial", Round[10^6 First[FindMinimum[Evaluate[qpTut], {{x, 2}, {y, 0}}, Method -> "SLSQP"]]]];

(* ---- 02 HS71: equality + inequality + bounds, f* = 17.0140173 at
        (1, 4.743, 3.821, 1.379). The canonical SLSQP test. ---- *)
hs71 = {x1 x4 (x1 + x2 + x3) + x3,
        x1 x2 x3 x4 >= 25 && x1^2 + x2^2 + x3^2 + x4^2 == 40 &&
        1 <= x1 <= 5 && 1 <= x2 <= 5 && 1 <= x3 <= 5 && 1 <= x4 <= 5};
bench["02 HS71", FindMinimum[Evaluate[hs71], {{x1, 1}, {x2, 5}, {x3, 5}, {x4, 1}}, Method -> "SLSQP"];];
check["02 HS71", Round[10^3 First[FindMinimum[Evaluate[hs71], {{x1, 1}, {x2, 5}, {x3, 5}, {x4, 1}}, Method -> "SLSQP"]]]];

(* ---- 03 nonlinear equality: min x^2+y^2 s.t. x y == 1 -> f* = 2 at (1,1)
        (AM-GM: x^2+y^2 >= 2|xy| = 2). ---- *)
nlEq = {x^2 + y^2, x y == 1};
bench["03 Nonlinear-eq xy=1", FindMinimum[Evaluate[nlEq], {{x, 2}, {y, 0.5}}, Method -> "SLSQP"];];
check["03 Nonlinear-eq xy=1", Round[10^6 First[FindMinimum[Evaluate[nlEq], {{x, 2}, {y, 0.5}}, Method -> "SLSQP"]]]];

(* ---- 04 simplex QP n=10, f* = 1/10 = 0.1 ---- *)
bench["04 Simplex QP n10", FindMinimum[Evaluate[splx10], Evaluate[splx10s], Method -> "SLSQP"];];
check["04 Simplex QP n10", Round[10^6 First[FindMinimum[Evaluate[splx10], Evaluate[splx10s], Method -> "SLSQP"]]]];

(* ---- 05 simplex QP n=50, f* = 1/50 = 0.02 ---- *)
bench["05 Simplex QP n50", FindMinimum[Evaluate[splx50], Evaluate[splx50s], Method -> "SLSQP"];];
check["05 Simplex QP n50", Round[10^6 First[FindMinimum[Evaluate[splx50], Evaluate[splx50s], Method -> "SLSQP"]]]];

(* ---- 06 simplex QP n=100 (scaling), f* = 1/100 = 0.01 ---- *)
bench["06 Simplex QP n100", FindMinimum[Evaluate[splx100], Evaluate[splx100s], Method -> "SLSQP"];];
check["06 Simplex QP n100", Round[10^6 First[FindMinimum[Evaluate[splx100], Evaluate[splx100s], Method -> "SLSQP"]]]];
