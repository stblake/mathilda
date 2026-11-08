(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 251 -- SIAM 100-Digit Challenge, SEALED HOLDOUT B, family D16
   (Appendix D Problem 16: heat equation from the Daubechies phi_2 data; the
   time at which the graph length of u(., t) falls to 5.4).

   Ten new problems with the capability exp 191 tests: the heat equation from
   ROUGH (discontinuous / kinked / incompatible) initial data, solved by the
   idiomatic PDE solver (NDSolveValue, method of lines), a functional of
   u(., t) (a point value, NIntegrate of u, u^2 or of the arc-length density
   Sqrt[1 + u_x^2]), and FindRoot (secant) for the time t* at which the
   functional takes a given value.  u_t = u_xx on [0,1] unless noted.
     bd16.1   u = 0 at both ends, u0 = 1:                        u(1/2, t) = 1/2     (easy)
     bd16.2   u = 0 at both ends, u0 = 1 - |2x - 1| (hat):        Int u dx = 1/4
     bd16.3   u = 0 at both ends, u0 = WaveletPsi[HaarWavelet[], x]
              (+1 on [0,1/2), -1 on [1/2,1)):                     arc length = 3
     bd16.4   u = 0 at both ends, u0 = x (jump at x = 1):          arc length = 2
     bd16.5   u_x = 0 at both ends (Neumann), u0 = UnitStep[x - 1/3]: u(0, t) = 1/3
     bd16.6   on [0, Pi], u = 0 at both ends, u0 = x (Pi - x):     Int u^2 dx = 1
     bd16.7   u_t = u_xx + 1 (source), u = 0 at both ends, u0 = 0: u(1/2, t) = 1/10
     bd16.8   u(0) = 0, u_x(1) + u(1) = 0 (Robin), u0 = 1:         u(1, t) = 1/4
     bd16.9   on the unit square, u = 0 on the boundary, u0 = 1:  u(1/2, 1/2, t) = 1/2 (2-D)
     bd16.10  u = 0 at both ends, u0 = WaveletPhi[HaarWavelet[], 3x - 1]
              (the box on [1/3, 2/3)):                             arc length = 28/10
              (harder: t* ~ 2e-4, the smoothing scale Sqrt[t*] ~ 0.014 must be resolved)
   Rows: "bd16.k machine" (default NDSolveValue / NIntegrate / FindRoot) and
   "bd16.k wp50" (WorkingPrecision -> 60 throughout, PrecisionGoal -> 50 in
   NDSolveValue, as exp 191).  Truths: truth.tsv (siam100_refs/bd16_ref.py:
   eigenfunction series == method of images / Laplace inversion).
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NDSolveValue", "NIntegrate", "FindRoot", "WaveletPsi", "WaveletPhi", "HaarWavelet"}];

pdeOpts[wp_] := If[wp === MachinePrecision, {}, {WorkingPrecision -> wp, PrecisionGoal -> wp - 10}];

(* t* with func[sol, s, wp] == target; sol = NDSolveValue of the 1-D problem on [0, len] *)
heatTime[src_, u0_, bcs_, len_, tmax_, func_, target_, {t0_, t1_}, wp_] := Module[{sol, F},
  sol = NDSolveValue[Join[{D[u[x, t], t] == D[u[x, t], {x, 2}] + src, u[x, 0] == u0}, bcs],
     u, {x, 0, len}, {t, 0, tmax}, Sequence @@ pdeOpts[wp]];
  F[s_?NumericQ] := func[sol, s, wp];
  t /. FindRoot[F[t] == target, {t, t0, t1}, WorkingPrecision -> wp]];

heatTime2D[target_, {t0_, t1_}, wp_] := Module[{sol, F},
  sol = NDSolveValue[{D[u[t, x, y], t] == D[u[t, x, y], x, x] + D[u[t, x, y], y, y],
      u[0, x, y] == 1, u[t, 0, y] == 0, u[t, 1, y] == 0, u[t, x, 0] == 0, u[t, x, 1] == 0},
     u, {t, 0, 1/5}, {x, 0, 1}, {y, 0, 1}, Sequence @@ pdeOpts[wp]];
  F[s_?NumericQ] := sol[s, 1/2, 1/2];
  t /. FindRoot[F[t] == target, {t, t0, t1}, WorkingPrecision -> wp]];

pointAt[x0_][sol_, s_, wp_] := sol[x0, s];
intU[len_][sol_, s_, wp_] := NIntegrate[sol[x, s], {x, 0, len}, WorkingPrecision -> wp];
intU2[len_][sol_, s_, wp_] := NIntegrate[sol[x, s]^2, {x, 0, len}, WorkingPrecision -> wp];
arcLen[len_][sol_, s_, wp_] := NIntegrate[Sqrt[1 + Derivative[1, 0][sol][x, s]^2], {x, 0, len},
   WorkingPrecision -> wp];

dir[len_] := {u[0, t] == 0, u[len, t] == 0};
neu = {Derivative[1, 0][u][0, t] == 0, Derivative[1, 0][u][1, t] == 0};
rob = {u[0, t] == 0, Derivative[1, 0][u][1, t] + u[1, t] == 0};

Do[With[{lab = p[[1]], f = p[[2]]},
    vrun[lab <> " machine", f[MachinePrecision], 10, 1, 90];
    vrun[lab <> " wp50", f[60], 48, 1, 60]],
  {p, {
    {"bd16.1", heatTime[0, 1, dir[1], 1, 1/2, pointAt[1/2], 1/2, {9/100, 1/10}, #] &},
    {"bd16.2", heatTime[0, 1 - Abs[2 x - 1], dir[1], 1, 1/2, intU[1], 1/4, {7/100, 8/100}, #] &},
    {"bd16.3", heatTime[0, WaveletPsi[HaarWavelet[], x], dir[1], 1, 1/10, arcLen[1], 3, {14/1000, 17/1000}, #] &},
    {"bd16.4", heatTime[0, x, dir[1], 1, 1/10, arcLen[1], 2, {27/10000, 32/10000}, #] &},
    {"bd16.5", heatTime[0, UnitStep[x - 1/3], neu, 1, 1/2, pointAt[0], 1/3, {55/1000, 65/1000}, #] &},
    {"bd16.6", heatTime[0, x (Pi - x), dir[Pi], Pi, 3, intU2[Pi], 1, {11/10, 12/10}, #] &},
    {"bd16.7", heatTime[1, 0, dir[1], 1, 1/2, pointAt[1/2], 1/10, {15/100, 18/100}, #] &},
    {"bd16.8", heatTime[0, 1, rob, 1, 1, pointAt[1], 1/4, {33/100, 37/100}, #] &},
    {"bd16.9", heatTime2D[1/2, {55/1000, 65/1000}, #] &},
    {"bd16.10", heatTime[0, WaveletPhi[HaarWavelet[], 3 x - 1], dir[1], 1, 1/50, arcLen[1], 28/10,
       {15/100000, 25/100000}, #] &}}}];
