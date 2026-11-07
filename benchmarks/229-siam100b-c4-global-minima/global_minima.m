(* bench-timeout: 1800 *)
(* ==========================================================================
   Experiment 229 -- SIAM 100-Digit Challenge, sealed holdout B for Problem 4
   ("Think Globally, Act Locally"): ten more multimodal objectives whose
   GLOBAL minimum value over a box is asked.

     bc4.1   sin(10x) + cos(11y) + (x^2+y^2)/4                  [-1,1]^2
     bc4.2   (x^2+y-11)^2 + (x+y^2-7)^2 + 3x  (tilted Himmelblau) [-5,5]^2
     bc4.3   e^sin(30x) + sin(35 e^y) + sin(sin(40z)) - sin(10(x+y+z))
             + (x^2+y^2+z^2)/4                                  [-1,1]^3
     bc4.4   shifted Ackley: -20 exp(-sqrt(((x-3/10)^2+(y+2/5)^2)/2)/5)
             - exp((cos 2 pi x + cos 2 pi y)/2) + e + 20        [-5,5]^2
     bc4.5   Shubert g(x) g(y), g(t) = Sum_{i=1..5} i cos((i+1)t + i)
                                                                [-10,10]^2
     bc4.6   Eggholder -(y+47) sin sqrt|x/2+y+47| - x sin sqrt|x-(y+47)|
             (the minimum lies ON the boundary x = 512)        [-512,512]^2
     bc4.7   cos x + cos(sqrt2 x) + cos(sqrt3 x)  (1-D, quasi-periodic)
                                                                [0,100]
     bc4.8   (sin 25x + sin 25y)/2 + (x^2+y^2)/10
             - 2 exp(-300((x-3/5)^2 + (y+7/20)^2))  (narrow well) [-1,1]^2
     bc4.9   100(y-x^2)^2 + (1-x)^2 + (3/10) sin(15x) sin(15y)  [-2,2]^2
     bc4.10  4-D coupled Rastrigin: Sum (x_i^2 - 10 cos 2 pi x_i)
             + (x1 x2 + x2 x3 + x3 x4)/2 + (3/10) x1            [-5,5]^4

   Idiomatic statement: the system's GLOBAL minimiser with the box as
   constraints -- NMinimize at machine precision and with WorkingPrecision ->
   60.  Labels "bc4.k machine" / "bc4.k wp50".  True values: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NMinimize"}];

g5[t_] := Sum[i Cos[(i + 1) t + i], {i, 5}];
obj = {
  {Sin[10 x] + Cos[11 y] + (x^2 + y^2)/4, {x, y}, 1},
  {(x^2 + y - 11)^2 + (x + y^2 - 7)^2 + 3 x, {x, y}, 5},
  {Exp[Sin[30 x]] + Sin[35 Exp[y]] + Sin[Sin[40 z]] - Sin[10 (x + y + z)] +
     (x^2 + y^2 + z^2)/4, {x, y, z}, 1},
  {-20 Exp[-Sqrt[((x - 3/10)^2 + (y + 2/5)^2)/2]/5] -
     Exp[(Cos[2 Pi x] + Cos[2 Pi y])/2] + E + 20, {x, y}, 5},
  {g5[x] g5[y], {x, y}, 10},
  {-(y + 47) Sin[Sqrt[Abs[x/2 + y + 47]]] - x Sin[Sqrt[Abs[x - (y + 47)]]], {x, y}, 512},
  {Cos[x] + Cos[Sqrt[2] x] + Cos[Sqrt[3] x], {x}, {0, 100}},
  {(Sin[25 x] + Sin[25 y])/2 + (x^2 + y^2)/10 -
     2 Exp[-300 ((x - 3/5)^2 + (y + 7/20)^2)], {x, y}, 1},
  {100 (y - x^2)^2 + (1 - x)^2 + 3/10 Sin[15 x] Sin[15 y], {x, y}, 2},
  {Sum[v^2 - 10 Cos[2 Pi v], {v, {x, y, z, w}}] + (x y + y z + z w)/2 + 3/10 x,
     {x, y, z, w}, 5}};

(* the box: half-width h for [-h, h]^d, or an explicit {lo, hi} *)
box[vars_, h_] := Apply[And, If[ListQ[h], h[[1]] <= # <= h[[2]], -h <= # <= h] & /@ vars];

Do[With[{f = obj[[k, 1]], vars = obj[[k, 2]], cons = box[obj[[k, 2]], obj[[k, 3]]]},
    vrun["bc4." <> ToString[k] <> " machine", First[NMinimize[{f, cons}, vars]], 10, 3, 60];
    vrun["bc4." <> ToString[k] <> " wp50",
         First[NMinimize[{f, cons}, vars, WorkingPrecision -> 60]], 48, 1, 120]],
  {k, 10}];
