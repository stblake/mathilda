(* Holdout-A D9 family (experiment 214), method B (independent of hd9_ref.py's
   sparse dictionary propagation): the exact rational transition matrix as a
   SparseArray, MatrixPower[Transpose[P], steps, e_start] (exact), summed over
   the target set.  Prints the exact numerator/denominator. *)
tm[pts_, nb_, lazy_, absorb_] := Module[{ix = AssociationThread[pts -> Range[Length[pts]]], P},
  P = SparseArray[Flatten[Table[With[{l = nb[p]},
        Which[absorb[p], {{ix[p], ix[p]} -> 1},
              lazy, Prepend[({ix[p], ix[#]} -> 1/(2 Length[l])) & /@ l, {ix[p], ix[p]} -> 1/2],
              True, ({ix[p], ix[#]} -> 1/Length[l]) & /@ l]], {p, pts}], 1], {Length[pts], Length[pts]}];
  {ix, P}];
walk[pts_, nb_, start_, steps_, target_, lazy_: False, absorb_: (False &)] := Module[{ix, P, v},
  {ix, P} = tm[pts, nb, lazy, absorb];
  v = MatrixPower[Transpose[P], steps, UnitVector[Length[pts], ix[start]]];
  Total[v[[ix /@ Select[pts, target]]]]];
tri[n_] := {Flatten[Table[{r, c}, {r, 0, n - 1}, {c, 0, r}], 1],
   Function[p, Select[{{p[[1]], p[[2]] - 1}, {p[[1]], p[[2]] + 1}, {p[[1]] - 1, p[[2]] - 1}, {p[[1]] - 1, p[[2]]},
      {p[[1]] + 1, p[[2]]}, {p[[1]] + 1, p[[2]] + 1}}, 0 <= #[[2]] <= #[[1]] < n &]]};
grid[n_, king_] := {Flatten[Table[{i, j}, {i, n}, {j, n}], 1],
   Function[p, Select[p + # & /@ If[king, DeleteCases[Tuples[{-1, 0, 1}, 2], {0, 0}], {{1, 0}, {-1, 0}, {0, 1}, {0, -1}}],
      AllTrue[#, 1 <= # <= n &] &]]};
cube[n_] := {Tuples[Range[n], 3], Function[p, Select[p + # & /@ Join[IdentityMatrix[3], -IdentityMatrix[3]], AllTrue[#, 1 <= # <= n &] &]]};
hyper[d_] := {Range[0, 2^d - 1], Function[p, BitXor[p, 2^#] & /@ Range[0, d - 1]]};
r = <||>;
{pp, nn} = tri[10]; r[1] = walk[pp, nn, {0, 0}, 20, #[[1]] == 9 &];
{pp, nn} = tri[30]; r[2] = walk[pp, nn, {0, 0}, 60, # == {0, 0} &];
{pp, nn} = tri[20]; r[3] = walk[pp, nn, {12, 6}, 50, (#[[2]] == 0 || #[[2]] == #[[1]] || #[[1]] == 19) &];
{pp, nn} = grid[20, False]; r[4] = walk[pp, nn, {1, 1}, 60, #[[1]] == 20 &];
{pp, nn} = grid[12, True]; r[5] = walk[pp, nn, {1, 1}, 30, #[[1]] == #[[2]] &];
{pp, nn} = cube[6]; r[6] = walk[pp, nn, {1, 1, 1}, 41, # == {6, 6, 6} &];
{pp, nn} = hyper[8]; r[7] = walk[pp, nn, 0, 30, # == 255 &];
{pp, nn} = tri[30]; r[8] = walk[pp, nn, {0, 0}, 60, #[[1]] == 29 &, True];
{pp, nn} = tri[30]; r[9] = walk[pp, nn, {0, 0}, 60, #[[1]] == 29 &, False, (#[[1]] == 29 &)];
{pp, nn} = tri[40]; r[10] = walk[pp, nn, {0, 0}, 80, #[[1]] == 39 &];
Do[Print["hd9.", k, "  ", N[r[k], 70]]; Print["   exact ", InputForm[r[k]]], {k, 10}];
