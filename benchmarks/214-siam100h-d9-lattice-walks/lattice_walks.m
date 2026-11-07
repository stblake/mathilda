(* bench-timeout: 1800 *)
(* ==========================================================================
   Experiment 214 -- SIAM 100-Digit Challenge HOLDOUT A, family D9
   (generalises experiment 184, "a triangular walk").

   The capability: the distribution of a finite random walk after a fixed
   number of steps, from the n x n row-stochastic transition matrix P --
   machine: the start row of MatrixPower[N[P], steps]; wp50: the EXACT start
   vector pushed through `steps` exact products v.P, then N[., 60].  Each
   answer is an exact rational.  A step goes to a uniformly chosen neighbour
   of the CURRENT point.  Ten instances, varied in graph, target and chain:

     hd9.1   triangle, 10 points/side (55),  20 steps from the top: P(bottom row)  (easy)
     hd9.2   triangle, 30/side (465),        60 steps: P(back at the top vertex)
     hd9.3   triangle, 20/side (210),        50 steps from (row 12, col 6):
                                             P(on the boundary)
     hd9.4   20 x 20 square grid (4-nbr),    60 steps from a corner: P(far row)
     hd9.5   12 x 12 king's-move grid (8-nbr), 30 steps from a corner:
                                             P(on the main diagonal)
     hd9.6   6 x 6 x 6 cube grid (6-nbr),    41 steps corner to opposite corner
     hd9.7   hypercube Q_8 (256),            30 steps: P(at the antipode)
     hd9.8   LAZY walk on the 30-triangle (stay w.p. 1/2), 60 steps: P(bottom row)
     hd9.9   30-triangle with the bottom row ABSORBING, 60 steps:
                                             P(hit the bottom row by step 60)
     hd9.10  triangle, 40/side (820),        80 steps: P(bottom row)  (harder)

   Lattice conventions as in 184: triangle point (r, c), 0 <= c <= r < n, has
   neighbours (r, c+-1), (r-1, c-1), (r-1, c), (r+1, c), (r+1, c+1) inside.
   The matrices are built (exact, dense) before the timed cases.  Caps: 30 s
   per machine run (warm-up + 3 timed), 90 s per wp50 run.
   True values and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"MatrixPower", "SparseArray", "AssociationThread"}];

(* transition matrix of the walk on pts with neighbour function nb;
   mode "plain" | "lazy" | an absorbing-point predicate *)
tmat[pts_, nb_, mode_] := Module[{ix = AssociationThread[pts -> Range[Length[pts]]]},
  Normal[SparseArray[Flatten[Table[With[{l = nb[p], i = ix[p]},
      Which[mode === "lazy", Prepend[({i, ix[#]} -> 1/(2 Length[l])) & /@ l, {i, i} -> 1/2],
            mode =!= "plain" && TrueQ[mode[p]], {{i, i} -> 1},
            True, ({i, ix[#]} -> 1/Length[l]) & /@ l]], {p, pts}], 1],
    {Length[pts], Length[pts]}]]];
pos[pts_, test_] := Flatten[Position[pts, _?test, {1}, Heads -> False]];

triPts[n_] := Flatten[Table[{r, c}, {r, 0, n - 1}, {c, 0, r}], 1];
triNb[n_][{r_, c_}] := Select[{{r, c - 1}, {r, c + 1}, {r - 1, c - 1}, {r - 1, c},
                               {r + 1, c}, {r + 1, c + 1}}, 0 <= #[[2]] <= #[[1]] < n &];
gridPts[n_] := Flatten[Table[{i, j}, {i, n}, {j, n}], 1];
gridNb[n_, ds_][p_] := Select[(p + #) & /@ ds, 1 <= #[[1]] <= n && 1 <= #[[2]] <= n &];
rook = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
king = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}, {1, 1}, {1, -1}, {-1, 1}, {-1, -1}};
cubePts[n_] := Flatten[Table[{i, j, k}, {i, n}, {j, n}, {k, n}], 2];
cubeNb[n_][p_] := Select[(p + #) & /@ {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}},
                         Min[#] >= 1 && Max[#] <= n &];
hypNb[d_][p_] := Table[BitXor[p, 2^i], {i, 0, d - 1}];

(* {P, start index, steps, target indices} per case *)
T10 = triPts[10]; T20 = triPts[20]; T30 = triPts[30]; T40 = triPts[40];
G20 = gridPts[20]; G12 = gridPts[12]; C6 = cubePts[6]; H8 = Range[0, 255];
prob[1] = {tmat[T10, triNb[10], "plain"], 1, 20, pos[T10, #[[1]] == 9 &]};
prob[2] = {tmat[T30, triNb[30], "plain"], 1, 60, {1}};
prob[3] = {tmat[T20, triNb[20], "plain"], First[pos[T20, # == {12, 6} &]], 50,
           pos[T20, #[[2]] == 0 || #[[2]] == #[[1]] || #[[1]] == 19 &]};
prob[4] = {tmat[G20, gridNb[20, rook], "plain"], 1, 60, pos[G20, #[[1]] == 20 &]};
prob[5] = {tmat[G12, gridNb[12, king], "plain"], 1, 30, pos[G12, #[[1]] == #[[2]] &]};
prob[6] = {tmat[C6, cubeNb[6], "plain"], 1, 41, {Length[C6]}};
prob[7] = {tmat[H8, hypNb[8], "plain"], 1, 30, {256}};
prob[8] = {tmat[T30, triNb[30], "lazy"], 1, 60, pos[T30, #[[1]] == 29 &]};
prob[9] = {tmat[T30, triNb[30], (#[[1]] == 29 &)], 1, 60, pos[T30, #[[1]] == 29 &]};
prob[10] = {tmat[T40, triNb[40], "plain"], 1, 80, pos[T40, #[[1]] == 39 &]};

walkM[{P_, s_, k_, tg_}] := Total[MatrixPower[N[P], k][[s, tg]]];
walk50[{P_, s_, k_, tg_}] := N[Total[Nest[#.P &, UnitVector[Length[P], s], k][[tg]]], 60];

vrun["hd9.1 machine", walkM[prob[1]], 10, 3, 30];
vrun["hd9.1 wp50", walk50[prob[1]], 48, 1, 90];
vrun["hd9.2 machine", walkM[prob[2]], 10, 3, 30];
vrun["hd9.2 wp50", walk50[prob[2]], 48, 1, 90];
vrun["hd9.3 machine", walkM[prob[3]], 10, 3, 30];
vrun["hd9.3 wp50", walk50[prob[3]], 48, 1, 90];
vrun["hd9.4 machine", walkM[prob[4]], 10, 3, 30];
vrun["hd9.4 wp50", walk50[prob[4]], 48, 1, 90];
vrun["hd9.5 machine", walkM[prob[5]], 10, 3, 30];
vrun["hd9.5 wp50", walk50[prob[5]], 48, 1, 90];
vrun["hd9.6 machine", walkM[prob[6]], 10, 3, 30];
vrun["hd9.6 wp50", walk50[prob[6]], 48, 1, 90];
vrun["hd9.7 machine", walkM[prob[7]], 10, 3, 30];
vrun["hd9.7 wp50", walk50[prob[7]], 48, 1, 90];
vrun["hd9.8 machine", walkM[prob[8]], 10, 3, 30];
vrun["hd9.8 wp50", walk50[prob[8]], 48, 1, 90];
vrun["hd9.9 machine", walkM[prob[9]], 10, 3, 30];
vrun["hd9.9 wp50", walk50[prob[9]], 48, 1, 90];
vrun["hd9.10 machine", walkM[prob[10]], 10, 3, 30];
vrun["hd9.10 wp50", walk50[prob[10]], 48, 1, 90];
