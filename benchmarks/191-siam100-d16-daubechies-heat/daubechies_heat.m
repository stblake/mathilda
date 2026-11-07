(* bench-timeout: 400 *)
(* ==========================================================================
   Experiment 191 -- SIAM 100-Digit Challenge, Appendix D, Problem 16
   ("Daubechies heat").

     phi_2 is the Daubechies scaling function (support [0,3]).  Solve
     u_t = u_xx on [0,3], u(0,t) = u(3,t) = 0, u(x,0) = phi_2(x).  At what
     time does the length of the curve x -> u(x,t) in the x-y plane become 5.4?

         t* = 5.21278254742981720768...e-7          (truth.tsv)

   The graph of phi_2 has INFINITE length (L(t) ~ -0.08 log t as t -> 0), so
   5.4 is reached only at t ~ 5e-7, when the smoothing scale sqrt(t) ~ 7e-4
   is far below any coarse resolution of phi_2.

   Idiomatic statement (the PDE, as posed): phi_2 from the built-in
   WaveletPhi[DaubechiesWavelet[2], x]; the heat equation by NDSolveValue
   (method of lines); the curve length by NIntegrate of Sqrt[1 + u_x^2];
   the time by FindRoot (secant) on that length.
     d16.machine          length 5.4 (the problem)
     d16.machine-reduced  length 4.5 instead (t* ~ 5.39e-3), where the
                          smoothing scale ~0.07 makes the problem benign
     d16.wp50             length 5.4 with WorkingPrecision -> 60 throughout
   Note: WL's WaveletPhi is the piecewise-LINEAR interpolant of the cascade
   values at spacing 2^-8 (MaxRecursion default), whose graph length is only
   5.52; it cannot represent phi_2 at the sqrt(t_star) ~ 7e-4 scale, so the
   machine row is limited by the built-in itself, not by NDSolve.
   Python (scipy 1.17 has no Daubechies scaling function since
   scipy.signal.cascade/daub were removed; mpmath has none), Maxima and FriCAS
   have no phi_2 built in: all three SKIP every row.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"WaveletPhi", "DaubechiesWavelet", "NDSolveValue", "NIntegrate", "FindRoot"}];

(* t at which the curve length of u(., t) equals len; t0, t1 secant starts. *)
heatTime[len_, tmax_, t0_, t1_] := Module[{sol, arc},
  sol = NDSolveValue[{D[u[x, t], t] == D[u[x, t], {x, 2}],
      u[x, 0] == WaveletPhi[DaubechiesWavelet[2], x],
      u[0, t] == 0, u[3, t] == 0}, u, {x, 0, 3}, {t, 0, tmax}];
  arc[s_?NumericQ] := NIntegrate[Sqrt[1 + Derivative[1, 0][sol][x, s]^2], {x, 0, 3}];
  t /. FindRoot[arc[t] == len, {t, t0, t1}]];

heatTimeWP[len_, tmax_, t0_, t1_] := Module[{sol, arc},
  sol = NDSolveValue[{D[u[x, t], t] == D[u[x, t], {x, 2}],
      u[x, 0] == WaveletPhi[DaubechiesWavelet[2], x, WorkingPrecision -> 60],
      u[0, t] == 0, u[3, t] == 0}, u, {x, 0, 3}, {t, 0, tmax},
      WorkingPrecision -> 60, PrecisionGoal -> 50];
  arc[s_?NumericQ] := NIntegrate[Sqrt[1 + Derivative[1, 0][sol][x, s]^2], {x, 0, 3},
      WorkingPrecision -> 60, PrecisionGoal -> 50];
  t /. FindRoot[arc[t] == len, {t, t0, t1}, WorkingPrecision -> 60]];

vrun["d16.machine", heatTime[54/10, 10^-5, 4*^-7, 6*^-7], 10];
vrun["d16.machine-reduced", heatTime[45/10, 2/100, 4*^-3, 6*^-3], 10];
vrun["d16.wp50", heatTimeWP[54/10, 10^-5, 4*^-7, 6*^-7], 48];
