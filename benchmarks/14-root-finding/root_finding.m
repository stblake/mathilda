(* Experiment 14 -- Numerical root finding.
   MEASURES src/numerical_roots/: findroot.c (Newton/secant), nroots_aberth.c and
   nroots_jt.c (all roots of a polynomial), nsolve_system.c (multivariate).
   Checks are the root VALUE rounded, or the root COUNT -- both invariant under
   the iteration path taken to get there. *)

Get["../harness.m"];
Get["../data.m"];

require[{"FindRoot", "NSolve", "NRoots", "Roots", "WorkingPrecision"}];

bench["FindRoot Cos[x]==x", FindRoot[Cos[x] == x, {x, 1}];];
check["FindRoot Cos[x]==x", Round[10^6 (x /. FindRoot[Cos[x] == x, {x, 1}])]];

bench["FindRoot x^3-2x-5", FindRoot[x^3 - 2 x - 5 == 0, {x, 2}];];
check["FindRoot x^3-2x-5",
  Round[10^6 (x /. FindRoot[x^3 - 2 x - 5 == 0, {x, 2}])]];

(* A 2x2 nonlinear system: the Jacobian path rather than the scalar path. *)
bench["FindRoot 2x2 nonlinear system",
  FindRoot[{x^2 + y^2 == 4, x y == 1}, {{x, 1.8}, {y, 0.6}}];];
check["FindRoot 2x2 nonlinear system",
  Round[10^6 (x /. FindRoot[{x^2 + y^2 == 4, x y == 1},
                            {{x, 1.8}, {y, 0.6}}])]];

(* All roots of a degree-50 polynomial: the companion-matrix / Aberth path. *)
p50 = denseUPoly[50, x];
bench["NSolve all roots, degree 50", NSolve[p50 == 0, x];];
check["NSolve all roots, degree 50", Length[NSolve[p50 == 0, x]]];

(* A root of a transcendental with a nearby competing root: tests bracketing. *)
bench["FindRoot Tan[x]==x near 4.5", FindRoot[Tan[x] == x, {x, 4.5}];];
check["FindRoot Tan[x]==x near 4.5",
  Round[10^4 (x /. FindRoot[Tan[x] == x, {x, 4.5}])]];
