(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 203 -- SIAM 100-Digit Challenge, HOLDOUT A for Problem 8
   ("Hot Plate"): ten new heat-equation crossing times.  Each asks WHEN a
   temperature of a heat problem (u_t = Laplacian u, zero initial data and
   boundary value 1 -- or initial data 1 and zero boundary) reaches a level,
   stated as C8 is in experiment 173: the separated (eigenfunction) series at
   the point, summed with NSum, inside FindRoot.

   w(t)    = (4/Pi) Sum (-1)^k/(2k+1) Exp[-(2k+1)^2 Pi^2 t/4]   (centre of [-1,1],
             zero ends, initial 1; a box with all faces at 1 has 1 - prod w)
   W(x,t)  = the same with Cos[(2k+1) Pi x/2] (value at x)
     hc8.1   square [-1,1]^2, one side at 5: centre reaches 1/2 (5/4 (1 - w^2)).
     hc8.2   square, all four sides at 1: centre reaches 1/2   (1 - w^2).
     hc8.3   rod [-1,1], ends at 0, initial 1: midpoint falls to 1/2 (w).
     hc8.4   the same rod: MEAN temperature falls to 1/2
             (8/Pi^2) Sum Exp[-(2k+1)^2 Pi^2 t/4]/(2k+1)^2.
     hc8.5   cube [-1,1]^3, all faces at 1: centre reaches 1/2 (1 - w^3).
     hc8.6   rectangle [-2,2] x [-1,1], all sides at 1: centre reaches 1/2
             (1 - w(t/4) w(t)).
     hc8.7   unit disk, boundary at 1: centre reaches 1/2 -- the Bessel series
             1 - Sum_n 2/(j_n J1(j_n)) Exp[-j_n^2 t], j_n = BesselJZero[0, n].
     hc8.8   unit ball, surface at 1: centre reaches 1/2
             (1 + 2 Sum_{n>=1} (-1)^n Exp[-n^2 Pi^2 t]).
     hc8.9   square, all sides at 1: the point (1/2, 0) reaches 1/2
             (1 - W(1/2,t) w(t)).
     hc8.10  square, all sides at 1: centre reaches 1/100 -- an early time, a
             small target (the harder one: FindRoot on a residual of size 1e-2).
   Rows: "<id> machine" = NSum + FindRoot at machine precision; "<id> wp50" =
   both with WorkingPrecision -> 60.  True values: truth.tsv
   (siam100_refs/hc8_ref.py: eigen-series vs image series / inverse Laplace /
   Jacobi transformation).
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NSum", "FindRoot", "BesselJZero"}];

w[t_?NumericQ, o___] := 4/Pi NSum[(-1)^k/(2 k + 1) Exp[-(2 k + 1)^2 Pi^2 t/4], {k, 0, Infinity}, o];
ww[x_, t_?NumericQ, o___] := 4/Pi NSum[(-1)^k/(2 k + 1) Cos[(2 k + 1) Pi x/2] Exp[-(2 k + 1)^2 Pi^2 t/4],
                                       {k, 0, Infinity}, o];
mean[t_?NumericQ, o___] := 8/Pi^2 NSum[Exp[-(2 k + 1)^2 Pi^2 t/4]/(2 k + 1)^2, {k, 0, Infinity}, o];
disk[t_?NumericQ, o___] := 1 - NSum[2/(BesselJZero[0, n] BesselJ[1, BesselJZero[0, n]]) Exp[-BesselJZero[0, n]^2 t],
                                    {n, 1, Infinity}, o];
ball[t_?NumericQ, o___] := 1 + 2 NSum[(-1)^n Exp[-n^2 Pi^2 t], {n, 1, Infinity}, o];


vrun["hc8.1 machine", t /. FindRoot[5/4 (1 - w[t]^2) == 1/2, {t, 0.2}], 10, 3];
vrun["hc8.1 wp50", t /. FindRoot[5/4 (1 - w[t, WorkingPrecision -> 60]^2) == 1/2, {t, 2/10}, WorkingPrecision -> 60], 48];
vrun["hc8.2 machine", t /. FindRoot[1 - w[t]^2 == 1/2, {t, 0.3}], 10, 3];
vrun["hc8.2 wp50", t /. FindRoot[1 - w[t, WorkingPrecision -> 60]^2 == 1/2, {t, 3/10}, WorkingPrecision -> 60], 48];
vrun["hc8.3 machine", t /. FindRoot[w[t] == 1/2, {t, 0.3}], 10, 3];
vrun["hc8.3 wp50", t /. FindRoot[w[t, WorkingPrecision -> 60] == 1/2, {t, 3/10}, WorkingPrecision -> 60], 48];
vrun["hc8.4 machine", t /. FindRoot[mean[t] == 1/2, {t, 0.2}], 10, 3];
vrun["hc8.4 wp50", t /. FindRoot[mean[t, WorkingPrecision -> 60] == 1/2, {t, 2/10}, WorkingPrecision -> 60], 48];
vrun["hc8.5 machine", t /. FindRoot[1 - w[t]^3 == 1/2, {t, 0.2}], 10, 3];
vrun["hc8.5 wp50", t /. FindRoot[1 - w[t, WorkingPrecision -> 60]^3 == 1/2, {t, 2/10}, WorkingPrecision -> 60], 48];
vrun["hc8.6 machine", t /. FindRoot[1 - w[t/4] w[t] == 1/2, {t, 0.4}], 10, 3];
vrun["hc8.6 wp50", t /. FindRoot[1 - w[t/4, WorkingPrecision -> 60] w[t, WorkingPrecision -> 60] == 1/2, {t, 4/10}, WorkingPrecision -> 60], 48];
vrun["hc8.7 machine", t /. FindRoot[disk[t] == 1/2, {t, 0.15}], 10];
vrun["hc8.7 wp50", t /. FindRoot[disk[t, WorkingPrecision -> 60] == 1/2, {t, 15/100}, WorkingPrecision -> 60], 48];
vrun["hc8.8 machine", t /. FindRoot[ball[t] == 1/2, {t, 0.1}], 10, 3];
vrun["hc8.8 wp50", t /. FindRoot[ball[t, WorkingPrecision -> 60] == 1/2, {t, 1/10}, WorkingPrecision -> 60], 48];
vrun["hc8.9 machine", t /. FindRoot[1 - ww[1/2, t] w[t] == 1/2, {t, 0.2}], 10, 3];
vrun["hc8.9 wp50", t /. FindRoot[1 - ww[1/2, t, WorkingPrecision -> 60] w[t, WorkingPrecision -> 60] == 1/2, {t, 2/10}, WorkingPrecision -> 60], 48];
vrun["hc8.10 machine", t /. FindRoot[1 - w[t]^2 == 1/100, {t, 0.06}], 10, 3];
vrun["hc8.10 wp50", t /. FindRoot[1 - w[t, WorkingPrecision -> 60]^2 == 1/100, {t, 6/100}, WorkingPrecision -> 60], 48];
