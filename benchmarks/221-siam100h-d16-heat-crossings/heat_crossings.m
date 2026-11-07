(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 221 -- SIAM 100-Digit Challenge HOLDOUT A, family hd16
   (generalises exp 191, App. D Problem 16, "Daubechies heat").

   Exp 191 solves u_t = u_xx on an interval from a rough initial profile and
   asks WHEN a functional of u(., t) (the graph length) reaches a level,
   stated as NDSolveValue (method of lines) + NIntegrate + FindRoot.  This
   family asks the same -- a heat-equation initial-boundary-value problem by
   NDSolveValue, a functional of the solution, FindRoot on the time -- for
   initial data whose solution has a series / image representation, so that
   the crossing time is known to 60+ digits.  Varied: smooth, kinked,
   discontinuous and boundary-incompatible data; Dirichlet, Neumann and
   inhomogeneous boundary values; point values and graph lengths; 1-D and 2-D.

   k  domain    u(x,0)                 BC                 t* where
   1  [0,1]     Sin[Pi x]              u = 0              graph length = 3/2
   2  [0,1]     1 - |2x - 1| (hat)     u = 0              u(1/2,t) = 1/2
   3  [0,1]     1 (incompatible)       u = 0              u(1/2,t) = 1/2
   4  [0,1]     x (1 - x)              u = 0              u(1/2,t) = 1/10
   5  [0,1]     1 on [1/4,3/4], else 0 u = 0              u(1/2,t) = 1/2
   6  [0,1]     x                      u_x = 0            u(0,t) = 1/3
   7  [0,1]     hat                    u = 0              graph length = 6/5
   8  [0,3]     1 on [1,2], else 0     u = 0              graph length = 4
   9  [0,1]     0                      u(0)=0, u(1)=1     u(1/4,t) = 1/10
   10 [0,1]^2   1 (incompatible)       u = 0              u(1/2,1/2,t) = 1/2
   (truth.tsv; graph length = Int Sqrt[1 + u_x^2] dx.)
   Rows hd16.k machine (NDSolveValue defaults, NIntegrate, secant FindRoot)
   and hd16.k wp50 (WorkingPrecision -> 60 throughout, as in exp 191).
   Python (scipy has no PDE solver), Maxima and FriCAS cannot state an
   initial-boundary-value problem with a built-in: all rows SKIP.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NDSolveValue", "NIntegrate", "FindRoot"}];

hat[x_] := 1 - Abs[2 x - 1];
ind[x_, a_, b_] := Piecewise[{{1, a <= x <= b}}, 0];

(* 1-D: {u0, L, bc, tmax, functional, level, t0, t1}; bc in {"D", "N", "D01"};
   functional "c"[x0] (point value) or "len" (graph length) *)
probs = {
  {Sin[Pi x], 1, "D", 1/5, "len", 3/2, 6/100, 7/100},
  {hat[x], 1, "D", 1/5, {"pt", 1/2}, 1/2, 45/1000, 55/1000},
  {1, 1, "D", 1/5, {"pt", 1/2}, 1/2, 9/100, 1/10},
  {x (1 - x), 1, "D", 1/5, {"pt", 1/2}, 1/10, 9/100, 1/10},
  {ind[x, 1/4, 3/4], 1, "D", 1/5, {"pt", 1/2}, 1/2, 55/1000, 65/1000},
  {x, 1, "N", 1/5, {"pt", 0}, 1/3, 85/1000, 95/1000},
  {hat[x], 1, "D", 1/5, "len", 6/5, 95/1000, 103/1000},
  {ind[x, 1, 2], 3, "D", 1/10, "len", 4, 17/1000, 19/1000},
  {0, 1, "D01", 1/5, {"pt", 1/4}, 1/10, 105/1000, 112/1000}};

bcs["D", L_] := {u[0, t] == 0, u[L, t] == 0};
bcs["N", L_] := {Derivative[1, 0][u][0, t] == 0, Derivative[1, 0][u][L, t] == 0};
bcs["D01", L_] := {u[0, t] == 0, u[L, t] == 1};

crossTime[{u0_, L_, bc_, tmax_, fun_, level_, t0_, t1_}, opts___] := Module[{sol, F},
  sol = NDSolveValue[Join[{D[u[x, t], t] == D[u[x, t], {x, 2}], u[x, 0] == u0}, bcs[bc, L]],
     u, {x, 0, L}, {t, 0, tmax}, opts];
  F[s_?NumericQ] := If[fun === "len",
     NIntegrate[Sqrt[1 + Derivative[1, 0][sol][x, s]^2], {x, 0, L}, opts],
     sol[fun[[2]], s]];
  t /. FindRoot[F[t] == level, {t, t0, t1}, opts]];

crossTime2D[opts___] := Module[{sol},
  sol = NDSolveValue[{D[u[x, y, t], t] == D[u[x, y, t], {x, 2}] + D[u[x, y, t], {y, 2}],
      u[x, y, 0] == 1, u[0, y, t] == 0, u[1, y, t] == 0, u[x, 0, t] == 0, u[x, 1, t] == 0},
     u, {x, 0, 1}, {y, 0, 1}, {t, 0, 1/5}, opts];
  t /. FindRoot[sol[1/2, 1/2, t] == 1/2, {t, 55/1000, 65/1000}, opts]];

Do[With[{p = probs[[k]], lab = "hd16." <> ToString[k]},
   vrun[lab <> " machine", crossTime[p], 10, 1, 60];
   vrun[lab <> " wp50", crossTime[p, WorkingPrecision -> 60], 48, 1, 90]],
  {k, Length[probs]}];
vrun["hd16.10 machine", crossTime2D[], 10, 1, 60];
vrun["hd16.10 wp50", crossTime2D[WorkingPrecision -> 60], 48, 1, 90];
