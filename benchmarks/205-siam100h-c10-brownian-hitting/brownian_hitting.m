(* bench-timeout: 1200 *)
(* ==========================================================================
   Experiment 205 -- SIAM 100-Digit Challenge, HOLDOUT A for Problem 10
   ("Hitting the Ends"): ten new Brownian-motion problems whose answer is a
   value of a harmonic function (or of the expected exit time) on a box, each
   stated -- as C10 is in experiment 175 -- as its separated (Fourier) series
   summed with NSum.  Rectangles are [-a,a] x [-b,b] unless said otherwise.

     hc10.1  4x1 rectangle, start at the centre: P(hit an end x = +-2).
     hc10.2  10x1 rectangle, start (3, 0): P(hit an end).
     hc10.3  10x1 rectangle, start (4, 1/4): P(hit an end).
     hc10.4  2x1 rectangle, start at the centre: P(hit a LONG side y = +-1/2)
             (the series in x: slow geometric decay, sech((2k+1) pi/4)).
     hc10.5  unit square [0,1]^2, start (1/2, 1/4): P(hit the bottom side).
     hc10.6  semi-infinite strip [0,inf) x [0,1], start (1, 1/2): P(hit the
             end x = 0)  (= (4/pi) arctan e^-pi: an easy one).
     hc10.7  expected exit time E[tau] of standard BM (generator Delta/2) from
             the centre of the 2x1 rectangle: w = b^2 - y^2 - sum ..., Delta w = -2.
     hc10.8  3x1x1 box, start at the centre: P(hit an end face) -- a DOUBLE
             series (16/pi^2) sum_{j,k} (-1)^(j+k)/((2j+1)(2k+1)) sech(..).
     hc10.9  unit cube, start (1/2,1/2,1/4): P(hit the bottom face) -- a double
             series of sin * sinh-ratio terms.
     hc10.10 4x1 rectangle, start (1.9, 0), 0.1 from an end: the series in y
             decays only like exp(-0.1 pi (2k+1)) (the harder one).

   Off-centre cosh(c x)/cosh(c a) is written in the overflow-free form
   rc[c, x, a] = (e^{-c(a-x)} + e^{-c(a+x)})/(1 + e^{-2ca}) in every system.
   Rows: "<id> machine" = NSum at machine precision; "<id> wp50" = NSum with
   WorkingPrecision -> 60.  True values and provenance: truth.tsv
   (siam100_refs/hc10_ref.py: the series vs the complementary series in the
   other variable(s) / closed form).
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NSum", "Sech"}];

rc[c_, x_, a_] := (Exp[-c (a - x)] + Exp[-c (a + x)])/(1 + Exp[-2 c a]);

p1[o___]  := 4/Pi NSum[(-1)^k/(2 k + 1) Sech[2 (2 k + 1) Pi], {k, 0, Infinity}, o];
p2[o___]  := 4/Pi NSum[(-1)^k/(2 k + 1) rc[(2 k + 1) Pi, 3, 5], {k, 0, Infinity}, o];
p3[o___]  := 4/Pi NSum[(-1)^k/(2 k + 1) Cos[(2 k + 1) Pi/4] rc[(2 k + 1) Pi, 4, 5], {k, 0, Infinity}, o];
p4[o___]  := 4/Pi NSum[(-1)^k/(2 k + 1) Sech[(2 k + 1) Pi/4], {k, 0, Infinity}, o];
p5[o___]  := 4/Pi NSum[Sin[(2 j + 1) Pi/2]/(2 j + 1)
               (Exp[-(2 j + 1) Pi/4] - Exp[-7 (2 j + 1) Pi/4])/(1 - Exp[-2 (2 j + 1) Pi]), {j, 0, Infinity}, o];
p6[o___]  := 4/Pi NSum[(-1)^k Exp[-(2 k + 1) Pi]/(2 k + 1), {k, 0, Infinity}, o];
p7[o___]  := 1/4 - 8/Pi^3 NSum[(-1)^k/(2 k + 1)^3 Sech[(2 k + 1) Pi], {k, 0, Infinity}, o];
p8[o___]  := 16/Pi^2 NSum[(-1)^(j + k)/((2 j + 1) (2 k + 1))
               Sech[3/2 Pi Sqrt[(2 j + 1)^2 + (2 k + 1)^2]], {j, 0, Infinity}, {k, 0, Infinity}, o];
p9[o___]  := 16/Pi^2 NSum[(-1)^(j + k)/((2 j + 1) (2 k + 1))
               With[{r = Pi Sqrt[(2 j + 1)^2 + (2 k + 1)^2]}, (Exp[-r/4] - Exp[-7 r/4])/(1 - Exp[-2 r])],
               {j, 0, Infinity}, {k, 0, Infinity}, o];
p10[o___] := 4/Pi NSum[(-1)^k/(2 k + 1) rc[(2 k + 1) Pi, 19/10, 2], {k, 0, Infinity}, o];

vrun["hc10.1 machine", p1[], 10, 3];   vrun["hc10.1 wp50", p1[WorkingPrecision -> 60], 48];
vrun["hc10.2 machine", p2[], 10, 3];   vrun["hc10.2 wp50", p2[WorkingPrecision -> 60], 48];
vrun["hc10.3 machine", p3[], 10, 3];   vrun["hc10.3 wp50", p3[WorkingPrecision -> 60], 48];
vrun["hc10.4 machine", p4[], 10, 3];   vrun["hc10.4 wp50", p4[WorkingPrecision -> 60], 48];
vrun["hc10.5 machine", p5[], 10, 3];   vrun["hc10.5 wp50", p5[WorkingPrecision -> 60], 48];
vrun["hc10.6 machine", p6[], 10, 3];   vrun["hc10.6 wp50", p6[WorkingPrecision -> 60], 48];
vrun["hc10.7 machine", p7[], 10, 3];   vrun["hc10.7 wp50", p7[WorkingPrecision -> 60], 48];
vrun["hc10.8 machine", p8[], 10];      vrun["hc10.8 wp50", p8[WorkingPrecision -> 60], 48];
vrun["hc10.9 machine", p9[], 10];      vrun["hc10.9 wp50", p9[WorkingPrecision -> 60], 48];
vrun["hc10.10 machine", p10[], 10, 3]; vrun["hc10.10 wp50", p10[WorkingPrecision -> 60], 48];
