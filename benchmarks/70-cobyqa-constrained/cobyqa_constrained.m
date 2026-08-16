(* Experiment 70 -- COBYQA (derivative-free quadratic-approximation, constrained)
   optimization (FindMinimum, Method -> "COBYQA").

   MEASURES src/numerical_calculus/findmin.c's derivative-free trust-region SQP
   with QUADRATIC interpolation models -- full quadratic models of f and every
   constraint built by finite differences on a structured stencil, a
   Byrd-Omojokun-style trust-region SQP step (the same active-set QP as SLSQP),
   an L1 penalty merit, and native equality + inequality + bound handling --
   against scipy.optimize.minimize(method="COBYQA") (Ragonneau & Zhang 2023).

   Both solvers are DERIVATIVE-FREE (no Jacobian either side). Unlike COBYLA
   (experiment 69, inequality-only), COBYQA handles equalities natively, so this
   suite mixes eq + ineq + bounds: HS71, a nonlinear equality, the equal-weight
   simplex QP, and constrained Rosenbrock (the curved valley COBYQA's quadratic
   models navigate where COBYLA's linear ones stall). The check is Round[10^3 * f]
   at the optimum (10^2 for HS71, whose optimum is irrational): both must reach it
   or the join CHECK-FAILs. Objectives/constraints are built with Total[...] over
   fresh global symbols and passed through Evaluate[] past FindMinimum's HoldAll. *)

Get["../harness.m"];

require[{"FindMinimum"}];

(* Equal-weight simplex QP: min Total[w^2] s.t. Total[w]==1, w_i>=0 -> 1/n. *)
splxVars[n_] := Table[Symbol["w" <> ToString[i]], {i, 1, n}];
splxProb[n_] := With[{w = splxVars[n]}, {Total[w^2], Total[w] == 1 && And @@ (# >= 0 & /@ w)}];
splxSpec[n_] := Table[{Symbol["w" <> ToString[i]], 1./n}, {i, 1, n}];
splx8 = splxProb[8]; splx8s = splxSpec[8];

(* ---- 01 HS71: eq + ineq + bounds, f* = 17.0140 at (1, 4.743, 3.821, 1.379) ---- *)
hs71 = {x1 x4 (x1 + x2 + x3) + x3,
        x1 x2 x3 x4 >= 25 && x1^2 + x2^2 + x3^2 + x4^2 == 40 &&
        1 <= x1 <= 5 && 1 <= x2 <= 5 && 1 <= x3 <= 5 && 1 <= x4 <= 5};
bench["01 HS71", FindMinimum[Evaluate[hs71], {{x1, 1}, {x2, 5}, {x3, 5}, {x4, 1}}, Method -> "COBYQA", MaxIterations -> 3000];];
check["01 HS71", Round[10^2 First[FindMinimum[Evaluate[hs71], {{x1, 1}, {x2, 5}, {x3, 5}, {x4, 1}}, Method -> "COBYQA", MaxIterations -> 3000]]]];

(* ---- 02 Nonlinear equality: min x^2+y^2 s.t. x y == 1 -> f* = 2 at (1,1).
        Start (1.5, 0.8): scipy COBYQA is start-sensitive on this hyperbola and
        stalls at the (2, 0.5) start; (1.5, 0.8) reaches the optimum in both. ---- *)
nlEq = {x^2 + y^2, x y == 1};
bench["02 Nonlinear-eq xy=1", FindMinimum[Evaluate[nlEq], {{x, 1.5}, {y, 0.8}}, Method -> "COBYQA", MaxIterations -> 3000];];
check["02 Nonlinear-eq xy=1", Round[10^3 First[FindMinimum[Evaluate[nlEq], {{x, 1.5}, {y, 0.8}}, Method -> "COBYQA", MaxIterations -> 3000]]]];

(* ---- 03 Constrained Rosenbrock: min (1-x)^2+100(y-x^2)^2 s.t. x+y<=1.
        On the line x+y==1 (unconstrained min (1,1) has x+y=2>1). The curved
        valley COBYQA's quadratic models navigate; COBYLA's linear ones stall. ---- *)
cros = {(1 - x)^2 + 100 (y - x^2)^2, x + y <= 1};
bench["03 Constrained-Rosenbrock", FindMinimum[Evaluate[cros], {{x, -1}, {y, 1}}, Method -> "COBYQA", MaxIterations -> 3000];];
check["03 Constrained-Rosenbrock", Round[10^3 First[FindMinimum[Evaluate[cros], {{x, -1}, {y, 1}}, Method -> "COBYQA", MaxIterations -> 3000]]]];

(* ---- 04 Simplex QP n=8, f* = 1/8 = 0.125 ---- *)
bench["04 Simplex QP n8", FindMinimum[Evaluate[splx8], Evaluate[splx8s], Method -> "COBYQA", MaxIterations -> 3000];];
check["04 Simplex QP n8", Round[10^3 First[FindMinimum[Evaluate[splx8], Evaluate[splx8s], Method -> "COBYQA", MaxIterations -> 3000]]]];
