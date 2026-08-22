(* Experiment 75 -- trust-krylov (GLTR / Lanczos) unconstrained optimization
   (FindMinimum, Method -> "TrustKrylov"; alias "trust-krylov").

   MEASURES src/numerical_calculus/findmin.c's trust-krylov -- the GLTR method
   Lanczos-tridiagonalizes the Hessian in the Krylov space of the gradient using
   Hessian-VECTOR products only (finite differences of the exact compiled
   gradient), solves the small tridiagonal trust-region subproblem exactly
   (Moré-Sorensen), and lifts the step back -- against
   scipy.optimize.minimize(method="trust-krylov"). scipy is given the analytic
   Jacobian AND Hessian, so the race is solver-vs-solver (both Hessian-free
   Krylov trust-region; the exact tridiagonal solve handles the indefinite/hard
   case). These are UNCONSTRAINED smooth problems: the far-start Rosenbrock
   family (indefinite Hessian en route), an ill-conditioned quadratic, and the
   Trid function (integer optimum). The check is the objective at the optimum
   rounded to 10^6: 0 for the f*=0 problems, the exact integer optimum for Trid
   (n=10 -> -210, n=20 -> -1520). Objectives are built with Sum/Total over fresh
   global symbols and passed through Evaluate[] past FindMinimum's HoldAll. *)

Get["../harness.m"];

require[{"FindMinimum"}];


rosen[n_]  := With[{v = Table[Symbol["u" <> ToString[i]], {i, 1, n}]},
                   Sum[100 (v[[j + 1]] - v[[j]]^2)^2 + (1 - v[[j]])^2, {j, 1, n - 1}]];
rosenS[n_] := Table[{Symbol["u" <> ToString[i]], -1.2}, {i, 1, n}];
trid[n_]   := With[{v = Table[Symbol["tv" <> ToString[i]], {i, 1, n}]},
                   Sum[(v[[i]] - 1)^2, {i, 1, n}] - Sum[v[[i]] v[[i - 1]], {i, 2, n}]];
tridS[n_]  := Table[{Symbol["tv" <> ToString[i]], 0.0}, {i, 1, n}];
quad[n_, c_] := With[{v = Table[Symbol["qv" <> ToString[i]], {i, 1, n}]},
                     Sum[10^(c (i - 1)/(n - 1)) v[[i]]^2, {i, 1, n}]];
quadS[n_]  := Table[{Symbol["qv" <> ToString[i]], 1.0}, {i, 1, n}];

(* Precompute out of the timed region. *)
r2 = rosen[2];   r2s = rosenS[2];
r10 = rosen[10]; r10s = rosenS[10];
r20 = rosen[20]; r20s = rosenS[20];
q20 = quad[20, 4]; q20s = quadS[20];
t10 = trid[10];  t10s = tridS[10];
t20 = trid[20];  t20s = tridS[20];

bench["01 Rosenbrock 2D",  FindMinimum[Evaluate[r2],  Evaluate[r2s],  Method -> "TrustKrylov", MaxIterations -> 5000];];
check["01 Rosenbrock 2D",  Round[10^6 First[FindMinimum[Evaluate[r2],  Evaluate[r2s],  Method -> "TrustKrylov", MaxIterations -> 5000]]]];

bench["02 Rosenbrock 10D", FindMinimum[Evaluate[r10], Evaluate[r10s], Method -> "TrustKrylov", MaxIterations -> 5000];];
check["02 Rosenbrock 10D", Round[10^6 First[FindMinimum[Evaluate[r10], Evaluate[r10s], Method -> "TrustKrylov", MaxIterations -> 5000]]]];

bench["03 Rosenbrock 20D", FindMinimum[Evaluate[r20], Evaluate[r20s], Method -> "TrustKrylov", MaxIterations -> 5000];];
check["03 Rosenbrock 20D", Round[10^6 First[FindMinimum[Evaluate[r20], Evaluate[r20s], Method -> "TrustKrylov", MaxIterations -> 5000]]]];

bench["04 Illcond quad n20 c1e4", FindMinimum[Evaluate[q20], Evaluate[q20s], Method -> "TrustKrylov", MaxIterations -> 5000];];
check["04 Illcond quad n20 c1e4", Round[10^6 First[FindMinimum[Evaluate[q20], Evaluate[q20s], Method -> "TrustKrylov", MaxIterations -> 5000]]]];

bench["05 Trid n10", FindMinimum[Evaluate[t10], Evaluate[t10s], Method -> "TrustKrylov", MaxIterations -> 5000];];
check["05 Trid n10", Round[10^6 First[FindMinimum[Evaluate[t10], Evaluate[t10s], Method -> "TrustKrylov", MaxIterations -> 5000]]]];

bench["06 Trid n20", FindMinimum[Evaluate[t20], Evaluate[t20s], Method -> "TrustKrylov", MaxIterations -> 5000];];
check["06 Trid n20", Round[10^6 First[FindMinimum[Evaluate[t20], Evaluate[t20s], Method -> "TrustKrylov", MaxIterations -> 5000]]]];
