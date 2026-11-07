(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 224 -- SIAM 100-Digit Challenge HOLDOUT A, family hd19
   (generalises exp 194, App. D Problem 19, "optimal design").

   Exp 194 asks for the minimum of the compliance J(c) = Int u of a finite
   variable-thickness design problem (n x n cells, c constant per cell,
   Q1 finite elements, -div(c grad u) = 1, u = 0 on the boundary) over
   1 <= c <= 100 with a volume constraint, stated as FindMinimum with box and
   equality constraints on a LinearSolve black box.  This family asks the same
   on other instances, varying the mesh, the contrast, the volume fraction,
   the dimension, the load and the domain:

   k  dim  cells       c range    mean c    load   min J (truth.tsv)
   1  2    4 x 4       [1, 10]    11/2      1      4.7147300632e-3
   2  2    4 x 4       [1, 100]   103/4     1      9.9445766575e-4 (1/4 of volume high)
   3  2    6 x 6       [1, 100]   101/2     1      5.2900203547e-4
   4  2    5 x 5       [1, 100]   101/2     1      5.1177925859e-4 (odd n)
   5  1    8           [1, 100]   101/2     1      1/808 (-(c u')' = 1 on [0,1])
   6  1    16          [1, 4]     5/2       1      0.02574462890625 (dyadic)
   7  2    4 x 4       [1, 100]   101/2     x      1.2889450927e-4 (load breaks symmetry)
   8  2    8 x 4       [1, 100]   101/2     1      1.6452995387e-3 (rectangle 2 x 1, h = 1/4)
   9  2    4 x 4       [1, 2]     3/2       1      0.018589924291... (rational: bang-bang)
   10 2    6 x 6       [1, 1000]  1003/4    1      1.0609030346e-4 (high contrast)
   (The 1-D and the bang-bang optima are rational; SLSQP stops ~1e-4 short
   on 7 and 10 -- the references are KKT + a duality certificate.)

   Discretisation (the same in every system): square cells of side h;
   2-D element stiffness c_e (1/6)[[4,-1,-2,-1],[-1,4,-1,-2],[-2,-1,4,-1],
   [-1,-2,-1,4]] (corners counter-clockwise); 1-D element stiffness
   (c_e/h)[[1,-1],[-1,1]]; load F_i = h^dim f(node_i) (exact consistent
   load for f = 1 and f = x); J = F.U, K(c) U = F.  J is convex in c, so the
   minimum VALUE is unique.  Rows hd19.k machine and hd19.k wp50 (FindMinimum
   WorkingPrecision -> 60, PrecisionGoal -> 50).
   Python: scipy SLSQP (bounds + equality, numerical gradient), mpmath has no
   constrained optimiser -> wp50 SKIP.  Maxima: fmin_cobyla (machine only).
   FriCAS: no constrained optimiser -> SKIP.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"FindMinimum", "LinearSolve"}];

(* {dim, nx, ny, lo, hi, mean, load f[x]} *)
probs = {
  {2, 4, 4, 1, 10, 11/2, 1 &}, {2, 4, 4, 1, 100, 103/4, 1 &}, {2, 6, 6, 1, 100, 101/2, 1 &},
  {2, 5, 5, 1, 100, 101/2, 1 &}, {1, 8, 0, 1, 100, 101/2, 1 &}, {1, 16, 0, 1, 4, 5/2, 1 &},
  {2, 4, 4, 1, 100, 101/2, # &}, {2, 8, 4, 1, 100, 101/2, 1 &}, {2, 4, 4, 1, 2, 3/2, 1 &},
  {2, 6, 6, 1, 1000, 1003/4, 1 &}};

kq = {{4, -1, -2, -1}, {-1, 4, -1, -2}, {-2, -1, 4, -1}, {-1, -2, -1, 4}}/6;

(* mesh data {element node lists, local stiffness, load vector, #interior nodes} *)
mesh[{1, nx_, _, __}, f_] := With[{h = 1/nx},
  {Table[{If[a >= 1, a, 0], If[a + 1 <= nx - 1, a + 1, 0]}, {a, 0, nx - 1}],
   {{1, -1}, {-1, 1}}/h, Table[h f[h i], {i, nx - 1}], nx - 1}];
mesh[{2, nx_, ny_, __}, f_] := With[{h = 1/ny},
  Module[{idx},
   idx[i_, j_] := If[1 <= i <= nx - 1 && 1 <= j <= ny - 1, (i - 1) (ny - 1) + j, 0];
   {Flatten[Table[{idx[a, b], idx[a + 1, b], idx[a + 1, b + 1], idx[a, b + 1]},
       {a, 0, nx - 1}, {b, 0, ny - 1}], 1],
    kq, Flatten[Table[h^2 f[h i], {i, nx - 1}, {j, ny - 1}]], (nx - 1) (ny - 1)}]];

(* compliance F.U of the design cv (a numeric black box) *)
compliance[{els_, loc_, F_, m_}, cv_ /; VectorQ[cv, NumericQ]] := Module[{K, e, r, s, idx},
  K = ConstantArray[0 First[cv], {m, m}];
  Do[idx = els[[e]];
     Do[If[idx[[r]] > 0 && idx[[s]] > 0,
          K[[idx[[r]], idx[[s]]]] += cv[[e]] loc[[r, s]]], {r, Length[idx]}, {s, Length[idx]}],
     {e, Length[els]}];
  F . LinearSolve[K, F]];

designMin[p_, opts___] := Module[{ms = mesh[p, p[[7]]], ne, cv},
  ne = Length[ms[[1]]]; cv = Array[cc, ne];
  First[FindMinimum[{compliance[ms, cv], And @@ Thread[p[[4]] <= cv <= p[[5]]],
      Total[cv] == p[[6]] ne},
    Transpose[{cv, ConstantArray[p[[6]], ne]}], opts]]];

Do[With[{p = probs[[k]], lab = "hd19." <> ToString[k]},
   vrun[lab <> " machine", designMin[p], 10, 1, 90];
   vrun[lab <> " wp50", designMin[p, WorkingPrecision -> 60, PrecisionGoal -> 50], 48, 1, 120]],
  {k, Length[probs]}];
