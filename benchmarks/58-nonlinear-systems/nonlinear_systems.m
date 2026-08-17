(* ==========================================================================
   Experiment 58 -- Nonlinear systems (group E)
   ==========================================================================
   WHAT IT MEASURES.  Multivariate FindRoot -- src/numerical_roots/findroot.c
   and nsolve_system.c -- solving systems of nonlinear equations from a starting
   point, the vector analogue of the scalar root-finding in experiment 14.  Two
   small systems (a polynomial and a transcendental), the Broyden tridiagonal
   test at N=10 and N=40 (a scaling row for the Jacobian solve), and the
   Burden-Faires 3x3 system.

   Baseline is scipy.optimize.fsolve (MINPACK hybrd/hybrj, compiled).  Each
   system has a unique root near the shared starting point, so both solvers land
   on the same solution and the check -- a coordinate or the sum of the solution
   -- agrees.

   Note the plain-symbol requirement: FindRoot rejects indexed variables xb[i],
   so the N-variable systems build distinct symbols v1..vN and inject the
   equations/spec with Evaluate (FindRoot holds its arguments).
   ========================================================================== *)

Get["../harness.m"];
Get["../data.m"];

require[{"FindRoot"}];

(* ---- small systems ----------------------------------------------------- *)

bench["FindRoot 2x2 polynomial", FindRoot[{x^2 + y^2 == 4, x y == 1}, {{x, 1.5}, {y, 0.5}}];];
check["FindRoot 2x2 polynomial",
      Round[10^6 (x /. FindRoot[{x^2 + y^2 == 4, x y == 1}, {{x, 1.5}, {y, 0.5}}])]];

bench["FindRoot 2x2 transcendental",
      FindRoot[{Cos[x] - y == 0, Sin[y] - x/2 == 0}, {{x, 1.0}, {y, 1.0}}];];
check["FindRoot 2x2 transcendental",
      Round[10^6 ({x, y} /. FindRoot[{Cos[x] - y == 0, Sin[y] - x/2 == 0},
                                     {{x, 1.0}, {y, 1.0}}]) . {1, 1}]];

bench["FindRoot Burden-Faires 3x3",
      FindRoot[{3 x - Cos[y z] - 1/2 == 0,
                x^2 - 81 (y + 1/10)^2 + Sin[z] + 106/100 == 0,
                Exp[-x y] + 20 z + (10 Pi - 3)/3 == 0},
               {{x, 1/10}, {y, 1/10}, {z, -1/10}}];];
check["FindRoot Burden-Faires 3x3",
      Round[10^6 ({x, y, z} /. FindRoot[{3 x - Cos[y z] - 1/2 == 0,
                x^2 - 81 (y + 1/10)^2 + Sin[z] + 106/100 == 0,
                Exp[-x y] + 20 z + (10 Pi - 3)/3 == 0},
               {{x, 1/10}, {y, 1/10}, {z, -1/10}}]) . {1, 1, 1}]];

(* ---- Broyden tridiagonal, N=10 and N=40 -------------------------------- *)
(* Built inline with distinct plain symbols p1..p10 / q1..q40 (FindRoot rejects
   indexed variables, and a Module wrapper around Symbol[] recurses). *)

v10 = Table[Symbol["p" <> ToString[i]], {i, 10}];
e10 = Table[(3 - 2 v10[[i]]) v10[[i]] - If[i > 1, v10[[i - 1]], 0]
            - 2 If[i < 10, v10[[i + 1]], 0] + 1 == 0, {i, 10}];
s10 = Transpose[{v10, ConstantArray[-1., 10]}];
bench["FindRoot Broyden tridiagonal N=10", FindRoot[Evaluate[e10], Evaluate[s10]];];
check["FindRoot Broyden tridiagonal N=10",
      Round[10^6 Total[v10 /. FindRoot[Evaluate[e10], Evaluate[s10]]]]];

v40 = Table[Symbol["q" <> ToString[i]], {i, 40}];
e40 = Table[(3 - 2 v40[[i]]) v40[[i]] - If[i > 1, v40[[i - 1]], 0]
            - 2 If[i < 40, v40[[i + 1]], 0] + 1 == 0, {i, 40}];
s40 = Transpose[{v40, ConstantArray[-1., 40]}];
bench["FindRoot Broyden tridiagonal N=40", FindRoot[Evaluate[e40], Evaluate[s40]];];
check["FindRoot Broyden tridiagonal N=40",
      Round[10^6 Total[v40 /. FindRoot[Evaluate[e40], Evaluate[s40]]]]];
