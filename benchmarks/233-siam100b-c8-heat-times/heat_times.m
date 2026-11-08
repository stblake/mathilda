(* bench-timeout: 1700 *)
(* ==========================================================================
   Experiment 233 -- SIAM 100-Digit Challenge, HOLDOUT B for Problem 8
   ("Hot Plate").  Ten new problems that need the same capability: the time
   at which a heat-equation solution, written as its separated (Fourier /
   eigenfunction) series and summed with NSum, reaches a given level --
   FindRoot around NSum.  u_t = Laplacian u, zero initial data unless stated.

   Shorthand (centre value of the 1-D problem on an interval of length 2,
   zero ends, initial value 1):
       w(t)  = (4/Pi) Sum_k (-1)^k/(2k+1) Exp[-(2k+1)^2 Pi^2 t/4]
       w1(t) = (4/Pi) Sum_k (-1)^k/(2k+1) Exp[-(2k+1)^2 Pi^2 t]   (length 1)

     bc8.1  [-1,1]^2, side x = 1 held at 5, the others at 0 (the original
            plate): when does the centre reach u = 3/10?  (5/4)(1 - w^2) = 3/10
     bc8.2  rod [0,1], u(0) = 0, u(1) = 1: when does the midpoint reach 1/4?
            u = 1/2 - (2/Pi) Sum_k (-1)^k/(2k+1) Exp[-(2k+1)^2 Pi^2 t]
     bc8.3  cube [-1,1]^3, every face at 1: centre reaches 9/10.  1 - w^3
     bc8.4  unit disk, rim at 1: centre reaches 1/2.
            u = 1 - 2 Sum_k Exp[-j_k^2 t]/(j_k J1(j_k)), j_k = BesselJZero[0, k]
     bc8.5  unit ball, surface at 1: centre reaches 1/2.
            u = 1 + 2 Sum_{n>=1} (-1)^n Exp[-n^2 Pi^2 t]
     bc8.6  the original plate: when does the MEAN temperature reach 1?
            (5/4)(1 - wbar^2), wbar = (8/Pi^2) Sum_k Exp[-(2k+1)^2 Pi^2 t/4]/(2k+1)^2
     bc8.7  rectangle [0,2] x [0,1], every side at 1: centre reaches 1/2.
            1 - w(t) w1(t)
     bc8.8  rod [0,1], both ends at 1: centre reaches 1/10, the solution
            stated as its IMAGE series 2 Sum_n (-1)^n Erfc[(2n+1)/(4 Sqrt[t])]
     bc8.9  [0,1]^2, every side at 1: the point (1/4,1/4) reaches 1/2.
            1 - v^2, v = (4/Pi) Sum_k Sin[(2k+1) Pi/4]/(2k+1) Exp[-(2k+1)^2 Pi^2 t]
     bc8.10 the original plate: centre reaches 10^-6 (an early time: the
            level is a small difference 1 - w^2 of order 1).

   Tiers: "<id> machine" = FindRoot / NSum at machine precision;
   "<id> wp50" = both with WorkingPrecision -> 60.  True values and two
   independent methods each (image series, theta transform, numerical
   Laplace inversion): truth.tsv, siam100_refs/bc8_ref.py.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NSum", "FindRoot", "BesselJZero", "Erfc"}];

(* each series takes the working precision as its last argument *)
w[t_?NumericQ, p_] := 4/Pi NSum[(-1)^k/(2 k + 1) Exp[-(2 k + 1)^2 Pi^2 t/4], {k, 0, Infinity}, WorkingPrecision -> p];
w1[t_?NumericQ, p_] := 4/Pi NSum[(-1)^k/(2 k + 1) Exp[-(2 k + 1)^2 Pi^2 t], {k, 0, Infinity}, WorkingPrecision -> p];
u2[t_?NumericQ, p_] := 1/2 - 2/Pi NSum[(-1)^k/(2 k + 1) Exp[-(2 k + 1)^2 Pi^2 t], {k, 0, Infinity}, WorkingPrecision -> p];
u4[t_?NumericQ, p_] := 1 - 2 NSum[Exp[-BesselJZero[0, k]^2 t]/(BesselJZero[0, k] BesselJ[1, BesselJZero[0, k]]),
                                   {k, 1, Infinity}, WorkingPrecision -> p];
u5[t_?NumericQ, p_] := 1 + 2 NSum[(-1)^n Exp[-n^2 Pi^2 t], {n, 1, Infinity}, WorkingPrecision -> p];
wb[t_?NumericQ, p_] := 8/Pi^2 NSum[Exp[-(2 k + 1)^2 Pi^2 t/4]/(2 k + 1)^2, {k, 0, Infinity}, WorkingPrecision -> p];
u8[t_?NumericQ, p_] := 2 NSum[(-1)^n Erfc[(2 n + 1)/(4 Sqrt[t])], {n, 0, Infinity}, WorkingPrecision -> p];
v9[t_?NumericQ, p_] := 4/Pi NSum[Sin[(2 k + 1) Pi/4]/(2 k + 1) Exp[-(2 k + 1)^2 Pi^2 t], {k, 0, Infinity}, WorkingPrecision -> p];

vrun["bc8.1 machine", t /. FindRoot[5/4 (1 - w[t, MachinePrecision]^2) == 3/10, {t, N[15/100]}], 10, 1, 80];
vrun["bc8.1 wp50", t /. FindRoot[5/4 (1 - w[t, 60]^2) == 3/10, {t, 15/100}, WorkingPrecision -> 60], 48, 1, 80];
vrun["bc8.2 machine", t /. FindRoot[u2[t, MachinePrecision] == 1/4, {t, N[1/10]}], 10, 1, 80];
vrun["bc8.2 wp50", t /. FindRoot[u2[t, 60] == 1/4, {t, 1/10}, WorkingPrecision -> 60], 48, 1, 80];
vrun["bc8.3 machine", t /. FindRoot[1 - w[t, MachinePrecision]^3 == 9/10, {t, N[4/10]}], 10, 1, 80];
vrun["bc8.3 wp50", t /. FindRoot[1 - w[t, 60]^3 == 9/10, {t, 4/10}, WorkingPrecision -> 60], 48, 1, 80];
vrun["bc8.4 machine", t /. FindRoot[u4[t, MachinePrecision] == 1/2, {t, N[2/10]}], 10, 1, 80];
vrun["bc8.4 wp50", t /. FindRoot[u4[t, 60] == 1/2, {t, 2/10}, WorkingPrecision -> 60], 48, 1, 80];
vrun["bc8.5 machine", t /. FindRoot[u5[t, MachinePrecision] == 1/2, {t, N[14/100]}], 10, 1, 80];
vrun["bc8.5 wp50", t /. FindRoot[u5[t, 60] == 1/2, {t, 14/100}, WorkingPrecision -> 60], 48, 1, 80];
vrun["bc8.6 machine", t /. FindRoot[5/4 (1 - wb[t, MachinePrecision]^2) == 1, {t, N[24/100]}], 10, 1, 80];
vrun["bc8.6 wp50", t /. FindRoot[5/4 (1 - wb[t, 60]^2) == 1, {t, 24/100}, WorkingPrecision -> 60], 48, 1, 80];
vrun["bc8.7 machine", t /. FindRoot[1 - w[t, MachinePrecision] w1[t, MachinePrecision] == 1/2, {t, N[9/100]}], 10, 1, 80];
vrun["bc8.7 wp50", t /. FindRoot[1 - w[t, 60] w1[t, 60] == 1/2, {t, 9/100}, WorkingPrecision -> 60], 48, 1, 80];
vrun["bc8.8 machine", t /. FindRoot[u8[t, MachinePrecision] == 1/10, {t, N[3/100]}], 10, 1, 80];
vrun["bc8.8 wp50", t /. FindRoot[u8[t, 60] == 1/10, {t, 3/100}, WorkingPrecision -> 60], 48, 1, 80];
vrun["bc8.9 machine", t /. FindRoot[1 - v9[t, MachinePrecision]^2 == 1/2, {t, N[3/100]}], 10, 1, 80];
vrun["bc8.9 wp50", t /. FindRoot[1 - v9[t, 60]^2 == 1/2, {t, 3/100}, WorkingPrecision -> 60], 48, 1, 80];
vrun["bc8.10 machine", t /. FindRoot[5/4 (1 - w[t, MachinePrecision]^2) == 10^-6, {t, N[2/100]}], 10, 1, 80];
vrun["bc8.10 wp50", t /. FindRoot[5/4 (1 - w[t, 60]^2) == 10^-6, {t, 2/100}, WorkingPrecision -> 60], 48, 1, 80];
