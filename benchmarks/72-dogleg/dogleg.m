(* Experiment 72 -- dogleg (Powell dogleg trust region) unconstrained
   optimization (FindMinimum, Method -> "Dogleg"; alias "dogleg").

   MEASURES src/numerical_calculus/findmin.c's dogleg -- it builds the dense
   Hessian and blends the Cauchy (steepest-descent) point with the Newton point
   along the dogleg path, intersected with the trust radius -- against
   scipy.optimize.minimize(method="dogleg"). scipy is given the analytic Jacobian
   AND Hessian, so the race is solver-vs-solver (both dogleg trust region).

   dogleg needs a POSITIVE-DEFINITE model Hessian (scipy raises LinAlgError
   otherwise), so unlike experiments 71/73/74/75 the problem set is PD
   THROUGHOUT: ill-conditioned diagonal quadratics (condition 1, 1e2, 1e4, 1e6,
   constant SPD Hessian) and the Trid function (constant tridiagonal SPD Hessian,
   integer optimum). The check is the objective at the optimum rounded to 10^6:
   0 for the f*=0 quadratics, the exact integer optimum for Trid (n=10 -> -210,
   n=20 -> -1520). Objectives are built with Sum/Total over fresh global symbols
   and passed through Evaluate[] past FindMinimum's HoldAll. *)

Get["../harness.m"];

require[{"FindMinimum"}];


quad[n_, c_] := With[{v = Table[Symbol["qv" <> ToString[i]], {i, 1, n}]},
                     Sum[10^(c (i - 1)/(n - 1)) v[[i]]^2, {i, 1, n}]];
quadS[n_]  := Table[{Symbol["qv" <> ToString[i]], 1.0}, {i, 1, n}];
trid[n_]   := With[{v = Table[Symbol["tv" <> ToString[i]], {i, 1, n}]},
                   Sum[(v[[i]] - 1)^2, {i, 1, n}] - Sum[v[[i]] v[[i - 1]], {i, 2, n}]];
tridS[n_]  := Table[{Symbol["tv" <> ToString[i]], 0.0}, {i, 1, n}];

(* Precompute out of the timed region. *)
q0  = quad[10, 0]; q0s  = quadS[10];
q2  = quad[10, 2]; q2s  = quadS[10];
q4  = quad[10, 4]; q4s  = quadS[10];
q6  = quad[10, 6]; q6s  = quadS[10];
t10 = trid[10];    t10s = tridS[10];
t20 = trid[20];    t20s = tridS[20];

bench["01 Sphere n10", FindMinimum[Evaluate[q0], Evaluate[q0s], Method -> "Dogleg", MaxIterations -> 5000];];
check["01 Sphere n10", Round[10^6 First[FindMinimum[Evaluate[q0], Evaluate[q0s], Method -> "Dogleg", MaxIterations -> 5000]]]];

bench["02 Illcond quad n10 c1e2", FindMinimum[Evaluate[q2], Evaluate[q2s], Method -> "Dogleg", MaxIterations -> 5000];];
check["02 Illcond quad n10 c1e2", Round[10^6 First[FindMinimum[Evaluate[q2], Evaluate[q2s], Method -> "Dogleg", MaxIterations -> 5000]]]];

bench["03 Illcond quad n10 c1e4", FindMinimum[Evaluate[q4], Evaluate[q4s], Method -> "Dogleg", MaxIterations -> 5000];];
check["03 Illcond quad n10 c1e4", Round[10^6 First[FindMinimum[Evaluate[q4], Evaluate[q4s], Method -> "Dogleg", MaxIterations -> 5000]]]];

bench["04 Illcond quad n10 c1e6", FindMinimum[Evaluate[q6], Evaluate[q6s], Method -> "Dogleg", MaxIterations -> 5000];];
check["04 Illcond quad n10 c1e6", Round[10^6 First[FindMinimum[Evaluate[q6], Evaluate[q6s], Method -> "Dogleg", MaxIterations -> 5000]]]];

bench["05 Trid n10", FindMinimum[Evaluate[t10], Evaluate[t10s], Method -> "Dogleg", MaxIterations -> 5000];];
check["05 Trid n10", Round[10^6 First[FindMinimum[Evaluate[t10], Evaluate[t10s], Method -> "Dogleg", MaxIterations -> 5000]]]];

bench["06 Trid n20", FindMinimum[Evaluate[t20], Evaluate[t20s], Method -> "Dogleg", MaxIterations -> 5000];];
check["06 Trid n20", Round[10^6 First[FindMinimum[Evaluate[t20], Evaluate[t20s], Method -> "Dogleg", MaxIterations -> 5000]]]];
