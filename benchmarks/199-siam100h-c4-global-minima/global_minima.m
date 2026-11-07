(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 199 -- SIAM 100-Digit Challenge, holdout A for Problem 4
   (family hc4: the global minimum of a multimodal function on a box).
   The statement is C4's (experiment 169): the system's GLOBAL minimiser
   with the box as constraints, at machine precision and with
   WorkingPrecision -> 60.

     k   objective                                               box            min
     1   six-hump camel (4 - 2.1x^2 + x^4/3)x^2 + xy + (4y^2-4)y^2 [-3,3]x[-2,2]  -1.0316...  (easy)
     2   Branin                                                  [-5,10]x[0,15] 5/(4 pi)    (easy, 3 global minima)
     3   Styblinski-Tang, 3-D: Sum (u^4 - 16u^2 + 5u)/2          [-5,5]^3       -117.49...
     4   Shubert: Sum i cos((i+1)x+i) * Sum i cos((i+1)y+i)      [-10,10]^2     -186.73...  (18 global minima)
     5   Holder table -|sin x cos y e^|1 - r/pi||                [-10,10]^2     -19.208...  (non-smooth)
     6   cross-in-tray -(|sin x sin y e^|100 - r/pi|| + 1)^0.1/10^4 [-10,10]^2  -2.0626...
     7   C4 with frequencies 30, 40, 50, 60, 10                  [-1,1]^2       -3.2988...
     8   Eggholder                                               [-512,512]^2   -959.64...  (minimiser ON the edge x = 512)
     9   tilted Himmelblau (x^2+y-11)^2 + (x+y^2-7)^2 + x + y    [-5,5]^2       -7.0762...  (4 near-equal basins)
     10  C4 with frequencies 70, 80, 90, 100, 12                 [-1,1]^2       -3.3635...  (harder)

   True values: truth.tsv; references siam100_refs/hc4_ref.{py,wl}.
   Maxima and FriCAS have no global minimiser: SKIP.  mpmath has no
   minimiser; its wp50 row is the scipy global search followed by
   mpmath.findroot on the gradient (C4's wp50 statement).
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NMinimize"}];

f[1] = (4 - 21/10 x^2 + x^4/3) x^2 + x y + (-4 + 4 y^2) y^2;            b[1] = {{-3, 3}, {-2, 2}};
f[2] = (y - 51/10/(4 Pi^2) x^2 + 5/Pi x - 6)^2 + 10 (1 - 1/(8 Pi)) Cos[x] + 10; b[2] = {{-5, 10}, {0, 15}};
f[3] = Sum[(u^4 - 16 u^2 + 5 u)/2, {u, {x, y, z}}];                     b[3] = {{-5, 5}, {-5, 5}, {-5, 5}};
f[4] = Sum[i Cos[(i + 1) x + i], {i, 5}] Sum[i Cos[(i + 1) y + i], {i, 5}]; b[4] = {{-10, 10}, {-10, 10}};
f[5] = -Abs[Sin[x] Cos[y] Exp[Abs[1 - Sqrt[x^2 + y^2]/Pi]]];            b[5] = {{-10, 10}, {-10, 10}};
f[6] = -(Abs[Sin[x] Sin[y] Exp[Abs[100 - Sqrt[x^2 + y^2]/Pi]]] + 1)^(1/10)/10000; b[6] = {{-10, 10}, {-10, 10}};
f[7] = Exp[Sin[30 x]] + Sin[40 Exp[y]] + Sin[50 Sin[x]] + Sin[Sin[60 y]] - Sin[10 (x + y)] + (x^2 + y^2)/4; b[7] = {{-1, 1}, {-1, 1}};
f[8] = -(y + 47) Sin[Sqrt[Abs[x/2 + y + 47]]] - x Sin[Sqrt[Abs[x - (y + 47)]]]; b[8] = {{-512, 512}, {-512, 512}};
f[9] = (x^2 + y - 11)^2 + (x + y^2 - 7)^2 + x + y;                       b[9] = {{-5, 5}, {-5, 5}};
f[10] = Exp[Sin[70 x]] + Sin[80 Exp[y]] + Sin[90 Sin[x]] + Sin[Sin[100 y]] - Sin[12 (x + y)] + (x^2 + y^2)/4; b[10] = {{-1, 1}, {-1, 1}};

Do[
  With[{obj = f[k], vars = If[k == 3, {x, y, z}, {x, y}],
        cons = And @@ MapThread[#2[[1]] <= #1 <= #2[[2]] &, {If[k == 3, {x, y, z}, {x, y}], b[k]}]},
    vrun["hc4." <> ToString[k] <> " machine", First[NMinimize[{obj, cons}, vars]], 10, 1, 60];
    vrun["hc4." <> ToString[k] <> " wp50", First[NMinimize[{obj, cons}, vars,
                                                         WorkingPrecision -> 60]], 48, 1, 90]],
  {k, 10}];
