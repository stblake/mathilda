(* bench-timeout: 1200 *)
(* ==========================================================================
   Experiment 242 -- SIAM 100-Digit Challenge, SEALED HOLDOUT B for Appendix D
   Problem 7 (experiment 182, "area of a region").

   The capability: the area of a planar region {f > g}, stated as experiment
   182 states it -- in 1-D form: FindRoot for the points where the slice
   length vanishes or changes formula, then NIntegrate of the explicit slice
   length w, which has square-root (or cube-/fourth-root) behaviour at those
   points.  Ten new regions:

     bd7.1  {x^2 < y < cos x}                      (easy: smooth integrand)
              x0: cos x = x^2;  Int_{-x0}^{x0} (cos x - x^2) dx
     bd7.2  {e^x - 2 < y < 1 - x^2}                (easy: two roots, smooth)
              e^x + x^2 = 3 at a ~ -1.7, b ~ 0.8;  Int_a^b (3 - x^2 - e^x) dx
     bd7.3  {sin x sin y > 1/2} in [0,Pi]^2        (no root; sqrt ends)
              Int_{Pi/6}^{5Pi/6} (Pi - 2 ArcSin[1/(2 sin x)]) dx
     bd7.4  {exp(-(y - x^2)^2) > 1/2 + y^2/8}      (two kinks, as in 182)
              s = Sqrt[-Log[1/2 + y^2/8]];  yA, yB: s = |y| (-Log g = y^2);
              Int_{yA}^{yB} 2 Sqrt[y+s] dy + Int_{yB}^{2} 2(Sqrt[y+s] - Sqrt[y-s]) dy
     bd7.5  bounded oval of {y^2 < x^3 - 3x + 1}   (roots of a cubic)
              Int_{r1}^{r2} 2 Sqrt[x^3 - 3x + 1] dx, r1 ~ -1.9, r2 ~ 0.35
     bd7.6  {cos x > x^2 + y^4}                    (fourth-root ends)
              4 Int_0^{x0} (cos x - x^2)^(1/4) dx  (x0 as in bd7.1)
     bd7.7  {sin(x+y) > x^2 + y^2 - 1}             (rotated: u = (x+y)/Sqrt 2)
              h = 1 + sin(Sqrt[2] u) - u^2, roots u1 ~ -0.55, u2 ~ 1.4;
              Int_{u1}^{u2} 2 Sqrt[h] du
     bd7.8  {exp(-(2y + x^3)^2) > y^2/16 + exp(cos y - 1)}
              (cube-root slice, |y|^(1/3) point at y = 0, sqrt ends at
              g = 1, y ~ +-3.68):  w = Surd[s - 2y, 3] - Surd[-s - 2y, 3]
     bd7.9  {x^2 + y^2 < 9/4, cos x + cos y > 1}   (min of two widths; two
              crossings x1 ~ 0.75, x2 ~ 1.30 of Sqrt[9/4-x^2] = ArcCos[1-cos x])
     bd7.10 experiment 182 with y^2/100 in place of y^2/32 (HARDER: three
              islands, five island ends + the cube-root kink = six FindRoots)

   machine: default FindRoot / NIntegrate;  wp50: WorkingPrecision -> 60 in
   every FindRoot and NIntegrate, PrecisionGoal -> 50.  Starting values read
   off a plot (exact rationals in the wp50 rows).
   True values and provenance: truth.tsv (siam100_refs/bd7_ref.{py,gp}).
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NIntegrate", "FindRoot", "Surd"}];


(* 1 *)
a1M[] := Module[{x0 = x /. FindRoot[Cos[x] == x^2, {x, 0.8}]},
  NIntegrate[Cos[x] - x^2, {x, -x0, x0}]];
a150[] := Module[{x0 = x /. FindRoot[Cos[x] == x^2, {x, 8/10}, WorkingPrecision -> 60]},
  NIntegrate[Cos[x] - x^2, {x, -x0, x0}, WorkingPrecision -> 60, PrecisionGoal -> 50]];
(* 2 *)
a2M[] := Module[{a = x /. FindRoot[Exp[x] + x^2 == 3, {x, -1.7}],
                 b = x /. FindRoot[Exp[x] + x^2 == 3, {x, 0.8}]},
  NIntegrate[3 - x^2 - Exp[x], {x, a, b}]];
a250[] := Module[{a = x /. FindRoot[Exp[x] + x^2 == 3, {x, -17/10}, WorkingPrecision -> 60],
                  b = x /. FindRoot[Exp[x] + x^2 == 3, {x, 8/10}, WorkingPrecision -> 60]},
  NIntegrate[3 - x^2 - Exp[x], {x, a, b}, WorkingPrecision -> 60, PrecisionGoal -> 50]];
(* 3 *)
a3M[] := NIntegrate[Pi - 2 ArcSin[1/(2 Sin[x])], {x, Pi/6, 5 Pi/6}];
a350[] := NIntegrate[Pi - 2 ArcSin[1/(2 Sin[x])], {x, Pi/6, 5 Pi/6}, WorkingPrecision -> 60, PrecisionGoal -> 50];
(* 4 *)
s4[y_] := Sqrt[-Log[1/2 + y^2/8]];
a4M[] := Module[{yA = y /. FindRoot[-Log[1/2 + y^2/8] == y^2, {y, -0.75}],
                 yB = y /. FindRoot[-Log[1/2 + y^2/8] == y^2, {y, 0.75}]},
  NIntegrate[2 Sqrt[y + s4[y]], {y, yA, 0, yB}] +
  NIntegrate[2 (Sqrt[y + s4[y]] - Sqrt[y - s4[y]]), {y, yB, 2}]];
a450[] := Module[{yA = y /. FindRoot[-Log[1/2 + y^2/8] == y^2, {y, -3/4}, WorkingPrecision -> 60],
                  yB = y /. FindRoot[-Log[1/2 + y^2/8] == y^2, {y, 3/4}, WorkingPrecision -> 60]},
  NIntegrate[2 Sqrt[y + s4[y]], {y, yA, 0, yB}, WorkingPrecision -> 60, PrecisionGoal -> 50] +
  NIntegrate[2 (Sqrt[y + s4[y]] - Sqrt[y - s4[y]]), {y, yB, 2}, WorkingPrecision -> 60, PrecisionGoal -> 50]];
(* 5 *)
a5M[] := Module[{r1 = x /. FindRoot[x^3 - 3 x + 1 == 0, {x, -1.9}],
                 r2 = x /. FindRoot[x^3 - 3 x + 1 == 0, {x, 0.35}]},
  NIntegrate[2 Sqrt[x^3 - 3 x + 1], {x, r1, r2}]];
a550[] := Module[{r1 = x /. FindRoot[x^3 - 3 x + 1 == 0, {x, -19/10}, WorkingPrecision -> 60],
                  r2 = x /. FindRoot[x^3 - 3 x + 1 == 0, {x, 35/100}, WorkingPrecision -> 60]},
  NIntegrate[2 Sqrt[x^3 - 3 x + 1], {x, r1, r2}, WorkingPrecision -> 60, PrecisionGoal -> 50]];
(* 6 *)
a6M[] := Module[{x0 = x /. FindRoot[Cos[x] == x^2, {x, 0.8}]},
  4 NIntegrate[(Cos[x] - x^2)^(1/4), {x, 0, x0}]];
a650[] := Module[{x0 = x /. FindRoot[Cos[x] == x^2, {x, 8/10}, WorkingPrecision -> 60]},
  4 NIntegrate[(Cos[x] - x^2)^(1/4), {x, 0, x0}, WorkingPrecision -> 60, PrecisionGoal -> 50]];
(* 7 *)
h7[u_] := 1 + Sin[Sqrt[2] u] - u^2;
a7M[] := Module[{u1 = u /. FindRoot[h7[u] == 0, {u, -0.55}],
                 u2 = u /. FindRoot[h7[u] == 0, {u, 1.4}]},
  NIntegrate[2 Sqrt[h7[u]], {u, u1, u2}]];
a750[] := Module[{u1 = u /. FindRoot[h7[u] == 0, {u, -55/100}, WorkingPrecision -> 60],
                  u2 = u /. FindRoot[h7[u] == 0, {u, 14/10}, WorkingPrecision -> 60]},
  NIntegrate[2 Sqrt[h7[u]], {u, u1, u2}, WorkingPrecision -> 60, PrecisionGoal -> 50]];
(* 8 *)
g8[y_] := y^2/16 + Exp[Cos[y] - 1];
w8[y_] := Surd[Sqrt[-Log[g8[y]]] - 2 y, 3] - Surd[-Sqrt[-Log[g8[y]]] - 2 y, 3];
a8M[] := Module[{m = y /. FindRoot[g8[y] == 1, {y, -3.7}],
                 p = y /. FindRoot[g8[y] == 1, {y, 3.7}]},
  NIntegrate[w8[y], {y, m, 0, p}]];
a850[] := Module[{m = y /. FindRoot[g8[y] == 1, {y, -37/10}, WorkingPrecision -> 60],
                  p = y /. FindRoot[g8[y] == 1, {y, 37/10}, WorkingPrecision -> 60]},
  NIntegrate[w8[y], {y, m, 0, p}, WorkingPrecision -> 60, PrecisionGoal -> 50]];
(* 9 *)
a9M[] := Module[{x1 = x /. FindRoot[Sqrt[9/4 - x^2] == ArcCos[1 - Cos[x]], {x, 0.75}],
                 x2 = x /. FindRoot[Sqrt[9/4 - x^2] == ArcCos[1 - Cos[x]], {x, 1.3}]},
  4 (NIntegrate[Sqrt[9/4 - x^2], {x, 0, x1}] + NIntegrate[ArcCos[1 - Cos[x]], {x, x1, x2}] +
     NIntegrate[Sqrt[9/4 - x^2], {x, x2, 3/2}])];
a950[] := Module[{x1 = x /. FindRoot[Sqrt[9/4 - x^2] == ArcCos[1 - Cos[x]], {x, 3/4}, WorkingPrecision -> 60],
                  x2 = x /. FindRoot[Sqrt[9/4 - x^2] == ArcCos[1 - Cos[x]], {x, 13/10}, WorkingPrecision -> 60]},
  4 (NIntegrate[Sqrt[9/4 - x^2], {x, 0, x1}, WorkingPrecision -> 60, PrecisionGoal -> 50] + NIntegrate[ArcCos[1 - Cos[x]], {x, x1, x2}, WorkingPrecision -> 60, PrecisionGoal -> 50] +
     NIntegrate[Sqrt[9/4 - x^2], {x, x2, 3/2}, WorkingPrecision -> 60, PrecisionGoal -> 50])];
(* 10 *)
g10[y_] := y^2/100 + Exp[Sin[y]];
w10[y_] := Surd[-y + Sqrt[-Log[g10[y]]], 3] - Surd[-y - Sqrt[-Log[g10[y]]], 3];
a10M[] := Module[{y1, y2, y3, yc, y4, y5},
  {y1, y2, y3, y4, y5} = Table[y /. FindRoot[g10[y] == 1, {y, s}], {s, {-7.94, -7.04, -3.04, 3.25, 5.85}}];
  yc = y /. FindRoot[-Log[g10[y]] == y^2, {y, -0.86}];
  NIntegrate[w10[y], {y, y1, y2}] + NIntegrate[w10[y], {y, y3, yc, 0}] + NIntegrate[w10[y], {y, y4, y5}]];
a1050[] := Module[{y1, y2, y3, yc, y4, y5},
  {y1, y2, y3, y4, y5} = Table[y /. FindRoot[g10[y] == 1, {y, s}, WorkingPrecision -> 60], {s, {-794/100, -704/100, -304/100, 325/100, 585/100}}];
  yc = y /. FindRoot[-Log[g10[y]] == y^2, {y, -86/100}, WorkingPrecision -> 60];
  NIntegrate[w10[y], {y, y1, y2}, WorkingPrecision -> 60, PrecisionGoal -> 50] + NIntegrate[w10[y], {y, y3, yc, 0}, WorkingPrecision -> 60, PrecisionGoal -> 50] + NIntegrate[w10[y], {y, y4, y5}, WorkingPrecision -> 60, PrecisionGoal -> 50]];

vrun["bd7.1 machine", a1M[], 10, 3, 60];   vrun["bd7.1 wp50", a150[], 48, 1, 60];
vrun["bd7.2 machine", a2M[], 10, 3, 60];   vrun["bd7.2 wp50", a250[], 48, 1, 60];
vrun["bd7.3 machine", a3M[], 10, 3, 60];   vrun["bd7.3 wp50", a350[], 48, 1, 60];
vrun["bd7.4 machine", a4M[], 10, 3, 60];   vrun["bd7.4 wp50", a450[], 48, 1, 60];
vrun["bd7.5 machine", a5M[], 10, 3, 60];   vrun["bd7.5 wp50", a550[], 48, 1, 60];
vrun["bd7.6 machine", a6M[], 10, 3, 60];   vrun["bd7.6 wp50", a650[], 48, 1, 60];
vrun["bd7.7 machine", a7M[], 10, 3, 60];   vrun["bd7.7 wp50", a750[], 48, 1, 60];
vrun["bd7.8 machine", a8M[], 10, 3, 60];   vrun["bd7.8 wp50", a850[], 48, 1, 60];
vrun["bd7.9 machine", a9M[], 10, 3, 60];   vrun["bd7.9 wp50", a950[], 48, 1, 60];
vrun["bd7.10 machine", a10M[], 10, 3, 60]; vrun["bd7.10 wp50", a1050[], 48, 1, 60];
