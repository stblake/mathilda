(* Experiment 66 -- Powell derivative-free optimization (FindMinimum,
   Method -> "Powell"; alias "PrincipalAxis").

   MEASURES src/numerical_calculus/findmin.c's Powell conjugate-direction method
   -- FindMinimum's only DERIVATIVE-FREE local method -- on 13 local-optimization
   test functions with KNOWN analytic optima, against
   scipy.optimize.minimize(method="Powell"). Cases span a curved valley
   (Rosenbrock), the flat-x-direction Beale (the regression case for the
   fm_bracket_line strict-`<` bracket growth), a singular Hessian (Powell's own
   quartic), a coupled Trid quadratic, an active bound corner, a NON-SMOOTH
   L1 objective (where the gradient methods stall and derivative-free wins), and
   separable / ill-conditioned quadratics for scaling.

   Unlike experiment 65 (L-BFGS-B), scipy is given NO Jacobian -- Powell is
   derivative-free, so the fair race is derivative-free-vs-derivative-free.

   The check is Round[10^6 * f] at the optimum (0 for the f*=0 problems), or the
   exact integer optimum for Trid (-50) and the bound corner (5): both systems
   must reach the true optimum to that precision or the join CHECK-FAILs, so the
   check itself is the accuracy gate. bench times the solve (efficiency).
   Objectives are built with Total[Table[...]] (Sum with a symbolic iterator
   burns the budget on a closed-form attempt). *)

Get["../harness.m"];

require[{"FindMinimum"}];

(* n-variable objectives / specs built out of fresh scalar symbols z1..zn. *)
sphereObj[n_] := Total[Table[Symbol["z" <> ToString[i]]^2, {i, 1, n}]];
sphereSpec[n_] := Table[{Symbol["z" <> ToString[i]], 1.0}, {i, 1, n}];
illObj[n_, c_] := Total[Table[10^(c (i - 1)/(n - 1)) Symbol["z" <> ToString[i]]^2, {i, 1, n}]];
illSpec[n_] := Table[{Symbol["z" <> ToString[i]], 1.0}, {i, 1, n}];
tridObj[n_] := Total[Table[(Symbol["z" <> ToString[i]] - 1)^2, {i, 1, n}]]
             - Total[Table[Symbol["z" <> ToString[i]] Symbol["z" <> ToString[i - 1]], {i, 2, n}]];
tridSpec[n_] := Table[{Symbol["z" <> ToString[i]], 0.0}, {i, 1, n}];
l1Obj[n_] := Total[Table[Abs[Symbol["z" <> ToString[i]] - i], {i, 1, n}]];
l1Spec[n_] := Table[{Symbol["z" <> ToString[i]], 0.0}, {i, 1, n}];

(* Precompute the n-var problems out of the timed region. *)
sph10 = sphereObj[10];   sph10s = sphereSpec[10];
sph20 = sphereObj[20];   sph20s = sphereSpec[20];
ill10c4 = illObj[10, 4]; ill10c4s = illSpec[10];
ill10c6 = illObj[10, 6]; ill10c6s = illSpec[10];
trid6 = tridObj[6];      trid6s = tridSpec[6];
l1n5 = l1Obj[5];         l1n5s = l1Spec[5];

wood = 100 (y - x^2)^2 + (1 - x)^2 + 90 (w - z^2)^2 + (1 - z)^2
     + 10.1 ((y - 1)^2 + (w - 1)^2) + 19.8 (y - 1) (w - 1);
powell = (x + 10 y)^2 + 5 (z - w)^2 + (y - 2 z)^4 + 10 (x - w)^4;
beale = (1.5 - x + x y)^2 + (2.25 - x + x y^2)^2 + (2.625 - x + x y^3)^2;

(* ---- 01 Rosenbrock 2D (curved valley), f* = 0 ---- *)
bench["01 Rosenbrock2D", FindMinimum[(1 - x)^2 + 100 (y - x^2)^2, {{x, -1.2}, {y, 1}}, Method -> "Powell", MaxIterations -> 2000];];
check["01 Rosenbrock2D", Round[10^6 First[FindMinimum[(1 - x)^2 + 100 (y - x^2)^2, {{x, -1.2}, {y, 1}}, Method -> "Powell", MaxIterations -> 2000]]]];

(* ---- 02 Beale, f* = 0 at (3, 0.5) (flat in x at y=1) ---- *)
bench["02 Beale", FindMinimum[Evaluate[beale], {{x, 1}, {y, 1}}, Method -> "Powell", MaxIterations -> 2000];];
check["02 Beale", Round[10^6 First[FindMinimum[Evaluate[beale], {{x, 1}, {y, 1}}, Method -> "Powell", MaxIterations -> 2000]]]];

(* ---- 03 Booth, f* = 0 at (1, 3) ---- *)
bench["03 Booth", FindMinimum[(x + 2 y - 7)^2 + (2 x + y - 5)^2, {{x, 0}, {y, 0}}, Method -> "Powell"];];
check["03 Booth", Round[10^6 First[FindMinimum[(x + 2 y - 7)^2 + (2 x + y - 5)^2, {{x, 0}, {y, 0}}, Method -> "Powell"]]]];

(* ---- 04 Matyas, f* = 0 at origin ---- *)
bench["04 Matyas", FindMinimum[0.26 (x^2 + y^2) - 0.48 x y, {{x, 3}, {y, 3}}, Method -> "Powell"];];
check["04 Matyas", Round[10^6 First[FindMinimum[0.26 (x^2 + y^2) - 0.48 x y, {{x, 3}, {y, 3}}, Method -> "Powell"]]]];

(* ---- 05 Powell singular 4D (singular Hessian at the optimum), f* = 0 ---- *)
bench["05 Powell-singular", FindMinimum[Evaluate[powell], {{x, 3}, {y, -1}, {z, 0}, {w, 1}}, Method -> "Powell", MaxIterations -> 3000];];
check["05 Powell-singular", Round[10^6 First[FindMinimum[Evaluate[powell], {{x, 3}, {y, -1}, {z, 0}, {w, 1}}, Method -> "Powell", MaxIterations -> 3000]]]];

(* ---- 06 Wood/Colville 4D, f* = 0 ---- *)
bench["06 Wood", FindMinimum[Evaluate[wood], {{x, -3}, {y, -1}, {z, -3}, {w, -1}}, Method -> "Powell", MaxIterations -> 3000];];
check["06 Wood", Round[10^6 First[FindMinimum[Evaluate[wood], {{x, -3}, {y, -1}, {z, -3}, {w, -1}}, Method -> "Powell", MaxIterations -> 3000]]]];

(* ---- 07 Trid n=6 (coupled quadratic), f* = -50 ---- *)
bench["07 Trid n6", FindMinimum[Evaluate[trid6], Evaluate[trid6s], Method -> "Powell", MaxIterations -> 2000];];
check["07 Trid n6", Round[First[FindMinimum[Evaluate[trid6], Evaluate[trid6s], Method -> "Powell", MaxIterations -> 2000]]]];

(* ---- 08 bound-active quadratic on [0,1]^2 (optimum on the corner), f* = 5 ---- *)
bench["08 Bound-corner", FindMinimum[(x - 2)^2 + (y - 3)^2, {{x, 0, 0, 1}, {y, 0, 0, 1}}, Method -> "Powell"];];
check["08 Bound-corner", Round[First[FindMinimum[(x - 2)^2 + (y - 3)^2, {{x, 0, 0, 1}, {y, 0, 0, 1}}, Method -> "Powell"]]]];

(* ---- 09 non-smooth L1 sum n=5, f* = 0 at (1,2,3,4,5) (gradient methods stall) ---- *)
bench["09 Nonsmooth-L1 n5", FindMinimum[Evaluate[l1n5], Evaluate[l1n5s], Method -> "Powell", MaxIterations -> 2000];];
check["09 Nonsmooth-L1 n5", Round[10^6 First[FindMinimum[Evaluate[l1n5], Evaluate[l1n5s], Method -> "Powell", MaxIterations -> 2000]]]];

(* ---- 10 separable quadratic n=10, f* = 0 ---- *)
bench["10 Sphere n10", FindMinimum[Evaluate[sph10], Evaluate[sph10s], Method -> "Powell", MaxIterations -> 2000];];
check["10 Sphere n10", Round[10^6 First[FindMinimum[Evaluate[sph10], Evaluate[sph10s], Method -> "Powell", MaxIterations -> 2000]]]];

(* ---- 11 separable quadratic n=20 (scaling), f* = 0 ---- *)
bench["11 Sphere n20", FindMinimum[Evaluate[sph20], Evaluate[sph20s], Method -> "Powell", MaxIterations -> 3000];];
check["11 Sphere n20", Round[10^6 First[FindMinimum[Evaluate[sph20], Evaluate[sph20s], Method -> "Powell", MaxIterations -> 3000]]]];

(* ---- 12 ill-conditioned quadratic n=10, condition 1e4, f* = 0 ---- *)
bench["12 illcond n10 c1e4", FindMinimum[Evaluate[ill10c4], Evaluate[ill10c4s], Method -> "Powell", MaxIterations -> 5000];];
check["12 illcond n10 c1e4", Round[10^6 First[FindMinimum[Evaluate[ill10c4], Evaluate[ill10c4s], Method -> "Powell", MaxIterations -> 5000]]]];

(* ---- 13 ill-conditioned quadratic n=10, condition 1e6, f* = 0 ---- *)
bench["13 illcond n10 c1e6", FindMinimum[Evaluate[ill10c6], Evaluate[ill10c6s], Method -> "Powell", MaxIterations -> 8000];];
check["13 illcond n10 c1e6", Round[10^6 First[FindMinimum[Evaluate[ill10c6], Evaluate[ill10c6s], Method -> "Powell", MaxIterations -> 8000]]]];
