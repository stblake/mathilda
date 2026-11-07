(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 244 -- SIAM 100-Digit Challenge, SEALED HOLDOUT B for Appendix D
   Problem 9 (experiment 184, "a triangular walk").

   The capability: the n-step distribution of a finite Markov chain as a
   matrix power -- MatrixPower of the machine-real transition matrix for the
   machine row, and the start vector pushed through n EXACT rational products
   v.P, then N[., 60], for the 50-digit row (exactly as experiment 184).
   Every walk steps to a neighbour of its current vertex chosen uniformly,
   except bd9.9.  Ten new walks (vertices / steps / question):

     bd9.1   path 0..20 (end vertices have 1 neighbour), from 0, 40 steps:
             P(end in 11..20)                                  (easy, 21 states)
     bd9.2   10x10 grid, 4-neighbour, from a corner, 30 steps: P(end in the far row)
     bd9.3   triangular lattice, 20 per side (210), from the top, 40 steps:
             P(end in the bottom row)                         (184 at a new size)
     bd9.4   6x6x6 cubic grid (216), from a corner, 24 steps: P(end on the far face)
     bd9.5   knight on 8x8 from a1, 20 moves: P(on h8)
     bd9.6   king on 8x8 from a1, 30 moves: P(on the 8th rank)
     bd9.7   hypercube Q_10 (1024), from 0, 40 steps: P(at the antipode)
     bd9.8   cycle C_17, 60 steps: P(back at the start)
     bd9.9   gambler's ruin on 0..30, up 2/5, down 3/5, 0 and 30 absorbing,
             from 15: P(absorbed at 30 within 200 steps)
     bd9.10  triangular lattice, 40 per side (820), 100 steps: P(bottom row) (HARDER)

   True values and provenance: truth.tsv (exact rationals).
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"MatrixPower", "SparseArray"}];

(* uniform walk on vertex list pts with neighbour function nb *)
mkP[pts_, nb_] := Module[{ix = AssociationThread[pts -> Range[Length[pts]]]},
  Normal[SparseArray[Flatten[Table[With[{l = nb[p]}, ({ix[p], ix[#]} -> 1/Length[l]) & /@ l],
    {p, pts}]], {Length[pts], Length[pts]}]]];
inGrid[q_, d_] := Min[q] >= 0 && And @@ MapThread[Less, {q, d}];
gridNb[d_, mv_, p_] := Select[(p + #) & /@ mv, inGrid[#, d] &];
gridPts[d_] := Tuples[Range[0, # - 1] & /@ d];
triPts[n_] := Flatten[Table[{r, c}, {r, 0, n - 1}, {c, 0, r}], 1];
triNb[n_, {r_, c_}] := Select[{{r, c - 1}, {r, c + 1}, {r - 1, c - 1}, {r - 1, c},
                               {r + 1, c}, {r + 1, c + 1}}, 0 <= #[[2]] <= #[[1]] < n &];
(* {P, start index, target indices, steps} *)
walk[pts_, nb_, s_, tg_, n_] := With[{ix = AssociationThread[pts -> Range[Length[pts]]]},
  {mkP[pts, nb], ix[s], ix /@ tg, n}];

king = DeleteCases[Tuples[{-1, 0, 1}, 2], {0, 0}];
knight = {{1, 2}, {2, 1}, {-1, 2}, {-2, 1}, {1, -2}, {2, -1}, {-1, -2}, {-2, -1}};
cube = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
W[1] := walk[gridPts[{21}], (gridNb[{21}, {{1}, {-1}}, #] &), {0}, Table[{i}, {i, 11, 20}], 40];
W[2] := walk[gridPts[{10, 10}], (gridNb[{10, 10}, {{1, 0}, {-1, 0}, {0, 1}, {0, -1}}, #] &), {0, 0},
             Table[{9, j}, {j, 0, 9}], 30];
W[3] := walk[triPts[20], (triNb[20, #] &), {0, 0}, Table[{19, c}, {c, 0, 19}], 40];
W[4] := walk[gridPts[{6, 6, 6}], (gridNb[{6, 6, 6}, cube, #] &), {0, 0, 0},
             Select[gridPts[{6, 6, 6}], First[#] == 5 &], 24];
W[5] := walk[gridPts[{8, 8}], (gridNb[{8, 8}, knight, #] &), {0, 0}, {{7, 7}}, 20];
W[6] := walk[gridPts[{8, 8}], (gridNb[{8, 8}, king, #] &), {0, 0}, Table[{i, 7}, {i, 0, 7}], 30];
W[7] := walk[Range[0, 1023], BitXor[#, 2^Range[0, 9]] &, 0, {1023}, 40];
W[8] := walk[Range[0, 16], {Mod[# + 1, 17], Mod[# - 1, 17]} &, 0, {0}, 60];
W[9] := {Normal[SparseArray[Join[{{1, 1} -> 1, {31, 31} -> 1},
            Table[{v + 1, v + 2} -> 2/5, {v, 1, 29}], Table[{v + 1, v} -> 3/5, {v, 1, 29}]], {31, 31}]],
         16, {31}, 200};
W[10] := walk[triPts[40], (triNb[40, #] &), {0, 0}, Table[{39, c}, {c, 0, 39}], 100];

Do[With[{w = W[k]},
  With[{P = w[[1]], s = w[[2]], tg = w[[3]], n = w[[4]], e = UnitVector[Length[w[[1]]], w[[2]]]},
    vrun["bd9." <> ToString[k] <> " machine", Total[MatrixPower[N[P], n][[s, tg]]], 10, 3, 60];
    vrun["bd9." <> ToString[k] <> " wp50", N[Total[Nest[#.P &, e, n][[tg]]], 60], 48, 1, 120]]],
  {k, 1, 10}];
