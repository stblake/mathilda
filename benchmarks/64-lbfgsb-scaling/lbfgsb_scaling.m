(* Experiment 64 -- L-BFGS-B scaling (FindMinimum, Method -> "LBFGSB").

   MEASURES src/numerical_calculus/findmin.c: the limited-memory BFGS local
   method (O(m*n) per iteration, m = 10) against scipy.optimize.minimize's
   L-BFGS-B, and against Mathilda's own full-memory QuasiNewton (O(n^2) per
   iteration) to show where limited memory earns its keep.

   The headline is the SCALING of time with dimension n on an ill-conditioned
   (condition 1e2) diagonal quadratic, which is unimodal so both systems and
   both methods reach the SAME optimum (0) -- the check is Round[objective], an
   integer, so a descent-path difference cannot desync the join. Bound-active
   cases use a separable box quadratic whose optimum is the corner (each x_i
   clamped to 1), value n. Small-n Rosenbrock rounds out the unconstrained set
   where the global optimum 0 is reachable (extended Rosenbrock is multimodal
   for n >= 4, so it is deliberately not used at large n).

   scipy is given the ANALYTIC Jacobian (jac=), matching Mathilda's exact
   compiled gradient, so the race is solver-vs-solver, not gradient-vs-FD. *)

Get["../harness.m"];

require[{"FindMinimum"}];

(* ---- problem builders (fresh scalar symbols z1..zn) ---------------------
   Objectives are built with Total[Table[...]], NOT Sum[...]: Sum evaluates its
   summand once with a symbolic iterator to attempt a closed form, and
   Symbol["z"<>ToString[i]] collapses to a single concrete symbol there, so Sum
   sees a geometric series and burns the whole budget on a (wrong) closed-form
   attempt. Table enumerates i = 1..n directly, materialising z1..zn. *)
cqObj[n_] := Total[Table[10^(2 (i - 1)/(n - 1)) Symbol["z" <> ToString[i]]^2, {i, 1, n}]];
cqSpec[n_] := Table[{Symbol["z" <> ToString[i]], 1.0}, {i, 1, n}];
bxObj[n_] := Total[Table[(Symbol["z" <> ToString[i]] - 2)^2, {i, 1, n}]];
bxSpec[n_] := Table[{Symbol["z" <> ToString[i]], 0, 0, 1}, {i, 1, n}];

(* Precompute the expanded objectives/specs OUTSIDE the timed region, so the
   benchmark measures the solve (which still pays Mathilda's symbolic-gradient
   + compile setup internally), not the Sum/Table expansion. *)
cqO50 = cqObj[50];   cqS50 = cqSpec[50];
cqO200 = cqObj[200]; cqS200 = cqSpec[200];
cqO1000 = cqObj[1000]; cqS1000 = cqSpec[1000];
bxO200 = bxObj[200]; bxS200 = bxSpec[200];
bxO1000 = bxObj[1000]; bxS1000 = bxSpec[1000];

(* ---- C-series: ill-conditioned quadratic, L-BFGS-B, scaling in n -------- *)
bench["C50 quadratic LBFGSB n=50",
  FindMinimum[Evaluate[cqO50], Evaluate[cqS50], Method -> "LBFGSB", MaxIterations -> 3000];];
check["C50 quadratic LBFGSB n=50",
  Round[First[FindMinimum[Evaluate[cqO50], Evaluate[cqS50], Method -> "LBFGSB", MaxIterations -> 3000]]]];

bench["C200 quadratic LBFGSB n=200",
  FindMinimum[Evaluate[cqO200], Evaluate[cqS200], Method -> "LBFGSB", MaxIterations -> 3000];];
check["C200 quadratic LBFGSB n=200",
  Round[First[FindMinimum[Evaluate[cqO200], Evaluate[cqS200], Method -> "LBFGSB", MaxIterations -> 3000]]]];

bench["C1000 quadratic LBFGSB n=1000",
  FindMinimum[Evaluate[cqO1000], Evaluate[cqS1000], Method -> "LBFGSB", MaxIterations -> 3000];];
check["C1000 quadratic LBFGSB n=1000",
  Round[First[FindMinimum[Evaluate[cqO1000], Evaluate[cqS1000], Method -> "LBFGSB", MaxIterations -> 3000]]]];

(* ---- Q-series: SAME quadratic, full-memory QuasiNewton (O(n^2)) --------- *)
bench["Q50 quadratic QuasiNewton n=50",
  FindMinimum[Evaluate[cqO50], Evaluate[cqS50], Method -> "QuasiNewton", MaxIterations -> 3000];];
check["Q50 quadratic QuasiNewton n=50",
  Round[First[FindMinimum[Evaluate[cqO50], Evaluate[cqS50], Method -> "QuasiNewton", MaxIterations -> 3000]]]];

bench["Q200 quadratic QuasiNewton n=200",
  FindMinimum[Evaluate[cqO200], Evaluate[cqS200], Method -> "QuasiNewton", MaxIterations -> 3000];];
check["Q200 quadratic QuasiNewton n=200",
  Round[First[FindMinimum[Evaluate[cqO200], Evaluate[cqS200], Method -> "QuasiNewton", MaxIterations -> 3000]]]];

bench["Q1000 quadratic QuasiNewton n=1000",
  FindMinimum[Evaluate[cqO1000], Evaluate[cqS1000], Method -> "QuasiNewton", MaxIterations -> 3000];];
check["Q1000 quadratic QuasiNewton n=1000",
  Round[First[FindMinimum[Evaluate[cqO1000], Evaluate[cqS1000], Method -> "QuasiNewton", MaxIterations -> 3000]]]];

(* ---- B-series: bound-constrained box (the "-B" path), optimum = n ------- *)
bench["B200 bound-box LBFGSB n=200",
  FindMinimum[Evaluate[bxO200], Evaluate[bxS200], Method -> "LBFGSB", MaxIterations -> 2000];];
check["B200 bound-box LBFGSB n=200",
  Round[First[FindMinimum[Evaluate[bxO200], Evaluate[bxS200], Method -> "LBFGSB", MaxIterations -> 2000]]]];

bench["B1000 bound-box LBFGSB n=1000",
  FindMinimum[Evaluate[bxO1000], Evaluate[bxS1000], Method -> "LBFGSB", MaxIterations -> 2000];];
check["B1000 bound-box LBFGSB n=1000",
  Round[First[FindMinimum[Evaluate[bxO1000], Evaluate[bxS1000], Method -> "LBFGSB", MaxIterations -> 2000]]]];

(* ---- R-series: small-n Rosenbrock, global optimum 0 -------------------- *)
bench["R2 rosenbrock LBFGSB n=2",
  FindMinimum[(1 - a)^2 + 100 (b - a^2)^2, {{a, -1.2}, {b, 1}}, Method -> "LBFGSB"];];
check["R2 rosenbrock LBFGSB n=2",
  Round[First[FindMinimum[(1 - a)^2 + 100 (b - a^2)^2, {{a, -1.2}, {b, 1}}, Method -> "LBFGSB"]]]];

rosObj5 = Total[Table[100 (Symbol["z" <> ToString[i + 1]] - Symbol["z" <> ToString[i]]^2)^2
              + (1 - Symbol["z" <> ToString[i]])^2, {i, 1, 4}]];
rosSpec5 = Table[{Symbol["z" <> ToString[i]], 0.5}, {i, 1, 5}];
bench["R5 rosenbrock LBFGSB n=5",
  FindMinimum[Evaluate[rosObj5], Evaluate[rosSpec5], Method -> "LBFGSB", MaxIterations -> 3000];];
check["R5 rosenbrock LBFGSB n=5",
  Round[First[FindMinimum[Evaluate[rosObj5], Evaluate[rosSpec5], Method -> "LBFGSB", MaxIterations -> 3000]]]];
