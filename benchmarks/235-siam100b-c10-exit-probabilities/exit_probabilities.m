(* bench-timeout: 1600 *)
(* ==========================================================================
   Experiment 235 -- SIAM 100-Digit Challenge, HOLDOUT B for Problem 10
   ("Hitting the Ends").  Ten new problems that need the same capability:
   a Brownian exit probability (harmonic measure) = the value of a harmonic
   function at the start point, written as its separated (Fourier /
   Fourier-Bessel) series and summed with NSum.  Several values are tiny, so
   the digits asked for are RELATIVE digits of a small number.

     bc10.1  4 x 1 rectangle, start at the centre: P(exit through an end).
             (4/Pi) Sum_k (-1)^k/(2k+1) Sech[2 (2k+1) Pi]                  (easy)
     bc10.2  10 x 1 rectangle [-5,5] x [-1/2,1/2], start (3, 0): P(end).
             (4/Pi) Sum_k (-1)^k/(2k+1) Cosh[3 (2k+1) Pi]/Cosh[5 (2k+1) Pi]
     bc10.3  [0,2] x [0,1], start (1/2, 1/2): P(exit through the far end x = 2).
             (4/Pi) Sum_k (-1)^k/(2k+1) Sinh[(2k+1) Pi/2]/Sinh[2 (2k+1) Pi]
     bc10.4  3-D box [-1/2,1/2] x [-1,1] x [-3/2,3/2], centre: P(exit through
             one of the two 1 x 2 ends).  A DOUBLE series:
             (16/Pi^2) Sum_{m,n} (-1)^(m+n)/((2m+1)(2n+1))
                       Sech[(3/2) Pi Sqrt[(2m+1)^2 + ((2n+1)/2)^2]]
     bc10.5  half-strip {x > 0, 0 < y < 1}, start (1, 1/3): P(exit through the
             end x = 0).  (4/Pi) Sum_k Sin[(2k+1) Pi/3]/(2k+1) Exp[-(2k+1) Pi]
     bc10.6  10 x 1 rectangle, Brownian motion WITH DRIFT (generator
             Laplacian/2 + d/dx), start at the centre: P(end).
             (4/Pi) Sum_k (-1)^k/(2k+1) Cosh[5]/Cosh[5 Sqrt[1 + ((2k+1) Pi)^2]]
     bc10.7  circular sector of angle Pi/3 and radius 1, start r = 1/2 on the
             bisector: P(exit through the arc).
             (4/Pi) Sum_k (-1)^k/(2k+1) 8^-(2k+1)                          (easy)
     bc10.8  cylinder of radius 1/2 and length 10, start at the centre:
             P(exit through an end).  A Fourier-Bessel series:
             Sum_{k>=1} 2/(j_k J1(j_k)) Sech[10 j_k],  j_k = BesselJZero[0, k]
     bc10.9  unit square, start (1/2, 99/100), next to the top side: P(top).
             (4/Pi) Sum_k (-1)^k/(2k+1) Sinh[(2k+1) Pi 99/100]/Sinh[(2k+1) Pi]
             (terms decay only like Exp[-(2k+1) Pi/100]: thousands of terms
             or a convergence accelerator -- the hard one)
     bc10.10 3-D box 1 x 1 x 10, centre: P(end) (~7e-10; double series).
             (16/Pi^2) Sum_{m,n} (-1)^(m+n)/((2m+1)(2n+1))
                       Sech[5 Pi Sqrt[(2m+1)^2 + (2n+1)^2]]

   Tiers: "<id> machine" = NSum at machine precision; "<id> wp50" = NSum with
   WorkingPrecision -> 60.  True values (each by the stated series AND by the
   complementary expansion 1 - P(other exits) from the separation in the
   other coordinate, or a closed form): truth.tsv, siam100_refs/bc10_ref.py.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NSum", "BesselJZero"}];

s1[p_] := 4/Pi NSum[(-1)^k/(2 k + 1) Sech[2 (2 k + 1) Pi], {k, 0, Infinity}, WorkingPrecision -> p];
s2[p_] := 4/Pi NSum[(-1)^k/(2 k + 1) Cosh[3 (2 k + 1) Pi]/Cosh[5 (2 k + 1) Pi], {k, 0, Infinity}, WorkingPrecision -> p];
s3[p_] := 4/Pi NSum[(-1)^k/(2 k + 1) Sinh[(2 k + 1) Pi/2]/Sinh[2 (2 k + 1) Pi], {k, 0, Infinity}, WorkingPrecision -> p];
s4[p_] := 16/Pi^2 NSum[(-1)^(m + n)/((2 m + 1) (2 n + 1)) Sech[3/2 Pi Sqrt[(2 m + 1)^2 + ((2 n + 1)/2)^2]],
                       {m, 0, Infinity}, {n, 0, Infinity}, WorkingPrecision -> p];
s5[p_] := 4/Pi NSum[Sin[(2 k + 1) Pi/3]/(2 k + 1) Exp[-(2 k + 1) Pi], {k, 0, Infinity}, WorkingPrecision -> p];
s6[p_] := 4/Pi NSum[(-1)^k/(2 k + 1) Cosh[5]/Cosh[5 Sqrt[1 + ((2 k + 1) Pi)^2]], {k, 0, Infinity}, WorkingPrecision -> p];
s7[p_] := 4/Pi NSum[(-1)^k/(2 k + 1) 8^-(2 k + 1), {k, 0, Infinity}, WorkingPrecision -> p];
s8[p_] := NSum[2/(BesselJZero[0, k] BesselJ[1, BesselJZero[0, k]]) Sech[10 BesselJZero[0, k]],
               {k, 1, Infinity}, WorkingPrecision -> p];
s9[p_] := 4/Pi NSum[(-1)^k/(2 k + 1) Sinh[(2 k + 1) Pi 99/100]/Sinh[(2 k + 1) Pi], {k, 0, Infinity}, WorkingPrecision -> p];
s10[p_] := 16/Pi^2 NSum[(-1)^(m + n)/((2 m + 1) (2 n + 1)) Sech[5 Pi Sqrt[(2 m + 1)^2 + (2 n + 1)^2]],
                        {m, 0, Infinity}, {n, 0, Infinity}, WorkingPrecision -> p];

vrun["bc10.1 machine", s1[MachinePrecision], 10, 1, 75];
vrun["bc10.1 wp50", s1[60], 48, 1, 75];
vrun["bc10.2 machine", s2[MachinePrecision], 10, 1, 75];
vrun["bc10.2 wp50", s2[60], 48, 1, 75];
vrun["bc10.3 machine", s3[MachinePrecision], 10, 1, 75];
vrun["bc10.3 wp50", s3[60], 48, 1, 75];
vrun["bc10.4 machine", s4[MachinePrecision], 10, 1, 75];
vrun["bc10.4 wp50", s4[60], 48, 1, 75];
vrun["bc10.5 machine", s5[MachinePrecision], 10, 1, 75];
vrun["bc10.5 wp50", s5[60], 48, 1, 75];
vrun["bc10.6 machine", s6[MachinePrecision], 10, 1, 75];
vrun["bc10.6 wp50", s6[60], 48, 1, 75];
vrun["bc10.7 machine", s7[MachinePrecision], 10, 1, 75];
vrun["bc10.7 wp50", s7[60], 48, 1, 75];
vrun["bc10.8 machine", s8[MachinePrecision], 10, 1, 75];
vrun["bc10.8 wp50", s8[60], 48, 1, 75];
vrun["bc10.9 machine", s9[MachinePrecision], 10, 1, 75];
vrun["bc10.9 wp50", s9[60], 48, 1, 75];
vrun["bc10.10 machine", s10[MachinePrecision], 10, 1, 75];
vrun["bc10.10 wp50", s10[60], 48, 1, 75];
