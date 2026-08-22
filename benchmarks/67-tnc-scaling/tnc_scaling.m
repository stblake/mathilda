(* Experiment 67 -- TNC (truncated Newton) scaling + accuracy (FindMinimum,
   Method -> "TNC"; alias "TruncatedNewton").

   MEASURES src/numerical_calculus/findmin.c's Hessian-free truncated Newton --
   an inner conjugate-gradient loop solving H.p=-g via Hessian-vector products
   (finite differences of the exact compiled gradient), with active-set bounds --
   against scipy.optimize.minimize(method="TNC"). Cases stress TNC's niche: large-n
   ill-conditioned quadratics (cond up to 1e6) and extended Rosenbrock valleys,
   where true curvature beats L-BFGS's low-rank model, plus active-bound corners.

   Like experiment 65 (L-BFGS-B) and unlike 66 (Powell), scipy is given the
   analytic Jacobian (jac=), so the race is solver-vs-solver -- both are
   gradient-based. The check is Round[10^6 * f] at the optimum (0 for the f*=0
   problems), or the exact integer optimum for Trid (-50), the bound corner (5),
   and many-active (160): both systems must reach it or the join CHECK-FAILs, so
   the check is the accuracy gate. Objectives are built with Total[Table[...]]. *)

Get["../harness.m"];

require[{"FindMinimum"}];

(* n-variable objectives / specs from fresh scalar symbols z1..zn. *)
illObj[n_, c_] := Total[Table[10^(c (i - 1)/(n - 1)) Symbol["z" <> ToString[i]]^2, {i, 1, n}]];
illSpec[n_] := Table[{Symbol["z" <> ToString[i]], 1.0}, {i, 1, n}];
rosObj[n_] := Total[Table[100 (Symbol["z" <> ToString[i + 1]] - Symbol["z" <> ToString[i]]^2)^2
                        + (1 - Symbol["z" <> ToString[i]])^2, {i, 1, n - 1}]];
rosSpec[n_] := Table[{Symbol["z" <> ToString[i]], -1.2}, {i, 1, n}];
tridObj[n_] := Total[Table[(Symbol["z" <> ToString[i]] - 1)^2, {i, 1, n}]]
             - Total[Table[Symbol["z" <> ToString[i]] Symbol["z" <> ToString[i - 1]], {i, 2, n}]];
tridSpec[n_] := Table[{Symbol["z" <> ToString[i]], 0.0}, {i, 1, n}];
(* many-active: sum (z_i-5)^2 on the box [-50,1] (4-element spec) -> the
   unconstrained min at 5 pins every coordinate to 1, value 10*(5-1)^2=160. *)
maObj[n_] := Total[Table[(Symbol["z" <> ToString[i]] - 5)^2, {i, 1, n}]];
maSpec[n_] := Table[{Symbol["z" <> ToString[i]], 0.0, -50.0, 1.0}, {i, 1, n}];

(* Precompute out of the timed region. *)
ill10c6 = illObj[10, 6];  ill10c6s = illSpec[10];
ill50c4 = illObj[50, 4];  ill50c4s = illSpec[50];
ill100c4 = illObj[100, 4]; ill100c4s = illSpec[100];
ros10 = rosObj[10]; ros10s = rosSpec[10];
ros30 = rosObj[30]; ros30s = rosSpec[30];
trid6 = tridObj[6]; trid6s = tridSpec[6];
ma10 = maObj[10];   ma10s = maSpec[10];

booth = (x + 2 y - 7)^2 + (2 x + y - 5)^2;
wood  = 100 (y - x^2)^2 + (1 - x)^2 + 90 (w - z^2)^2 + (1 - z)^2
      + 10.1 ((y - 1)^2 + (w - 1)^2) + 19.8 (y - 1) (w - 1);

(* ---- 01 Rosenbrock 2D (curved valley), f* = 0 ---- *)
bench["01 Rosenbrock2D", FindMinimum[(1 - x)^2 + 100 (y - x^2)^2, {{x, -1.2}, {y, 1}}, Method -> "TNC", MaxIterations -> 2000];];
check["01 Rosenbrock2D", Round[10^6 First[FindMinimum[(1 - x)^2 + 100 (y - x^2)^2, {{x, -1.2}, {y, 1}}, Method -> "TNC", MaxIterations -> 2000]]]];

(* ---- 02 extended Rosenbrock n=10, f* = 0 ---- *)
bench["02 Rosenbrock n10", FindMinimum[Evaluate[ros10], Evaluate[ros10s], Method -> "TNC", MaxIterations -> 3000];];
check["02 Rosenbrock n10", Round[10^6 First[FindMinimum[Evaluate[ros10], Evaluate[ros10s], Method -> "TNC", MaxIterations -> 3000]]]];

(* ---- 03 extended Rosenbrock n=30 (scaling), f* = 0 ---- *)
bench["03 Rosenbrock n30", FindMinimum[Evaluate[ros30], Evaluate[ros30s], Method -> "TNC", MaxIterations -> 5000];];
check["03 Rosenbrock n30", Round[10^6 First[FindMinimum[Evaluate[ros30], Evaluate[ros30s], Method -> "TNC", MaxIterations -> 5000]]]];

(* ---- 04 ill-conditioned quadratic n=10, cond 1e6, f* = 0 ---- *)
bench["04 illcond n10 c1e6", FindMinimum[Evaluate[ill10c6], Evaluate[ill10c6s], Method -> "TNC", MaxIterations -> 3000];];
check["04 illcond n10 c1e6", Round[10^6 First[FindMinimum[Evaluate[ill10c6], Evaluate[ill10c6s], Method -> "TNC", MaxIterations -> 3000]]]];

(* ---- 05 ill-conditioned quadratic n=50, cond 1e4, f* = 0 ---- *)
bench["05 illcond n50 c1e4", FindMinimum[Evaluate[ill50c4], Evaluate[ill50c4s], Method -> "TNC", MaxIterations -> 3000];];
check["05 illcond n50 c1e4", Round[10^6 First[FindMinimum[Evaluate[ill50c4], Evaluate[ill50c4s], Method -> "TNC", MaxIterations -> 3000]]]];

(* ---- 06 ill-conditioned quadratic n=100, cond 1e4 (scaling), f* = 0 ---- *)
bench["06 illcond n100 c1e4", FindMinimum[Evaluate[ill100c4], Evaluate[ill100c4s], Method -> "TNC", MaxIterations -> 4000];];
check["06 illcond n100 c1e4", Round[10^6 First[FindMinimum[Evaluate[ill100c4], Evaluate[ill100c4s], Method -> "TNC", MaxIterations -> 4000]]]];

(* ---- 07 Trid n=6 (coupled quadratic), f* = -50 ---- *)
bench["07 Trid n6", FindMinimum[Evaluate[trid6], Evaluate[trid6s], Method -> "TNC", MaxIterations -> 2000];];
check["07 Trid n6", Round[First[FindMinimum[Evaluate[trid6], Evaluate[trid6s], Method -> "TNC", MaxIterations -> 2000]]]];

(* ---- 08 bound-active quadratic on [0,1]^2 (corner), f* = 5 ---- *)
bench["08 Bound-corner", FindMinimum[(x - 2)^2 + (y - 3)^2, {{x, 0, 0, 1}, {y, 0, 0, 1}}, Method -> "TNC"];];
check["08 Bound-corner", Round[First[FindMinimum[(x - 2)^2 + (y - 3)^2, {{x, 0, 0, 1}, {y, 0, 0, 1}}, Method -> "TNC"]]]];

(* ---- 09 many-active bounds n=10 (all z_i <= 1), f* = 160 ---- *)
bench["09 Many-active n10", FindMinimum[Evaluate[ma10], Evaluate[ma10s], Method -> "TNC", MaxIterations -> 2000];];
check["09 Many-active n10", Round[First[FindMinimum[Evaluate[ma10], Evaluate[ma10s], Method -> "TNC", MaxIterations -> 2000]]]];

(* ---- 10 Booth, f* = 0 at (1,3) ---- *)
bench["10 Booth", FindMinimum[Evaluate[booth], {{x, 0}, {y, 0}}, Method -> "TNC"];];
check["10 Booth", Round[10^6 First[FindMinimum[Evaluate[booth], {{x, 0}, {y, 0}}, Method -> "TNC"]]]];

(* ---- 11 Wood/Colville 4D, f* = 0 ---- *)
bench["11 Wood", FindMinimum[Evaluate[wood], {{x, -3}, {y, -1}, {z, -3}, {w, -1}}, Method -> "TNC", MaxIterations -> 2000];];
check["11 Wood", Round[10^6 First[FindMinimum[Evaluate[wood], {{x, -3}, {y, -1}, {z, -3}, {w, -1}}, Method -> "TNC", MaxIterations -> 2000]]]];
