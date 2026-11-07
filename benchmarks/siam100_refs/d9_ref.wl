(* Reference for d9: exact rational via MatrixPower of the exact transition matrix. *)
n = 30;
pts = Flatten[Table[{r, c}, {r, 0, n - 1}, {c, 0, r}], 1];
id = AssociationThread[pts -> Range[Length[pts]]];
nb[{r_, c_}] := Select[{{r, c - 1}, {r, c + 1}, {r - 1, c - 1}, {r - 1, c}, {r + 1, c}, {r + 1, c + 1}},
   0 <= #[[2]] <= #[[1]] < n &];
P = SparseArray[Flatten[Table[With[{l = nb[p]}, ({id[p], id[#]} -> 1/Length[l]) & /@ l], {p, pts}]], {Length[pts], Length[pts]}];
row = MatrixPower[Normal[P], 60][[1]];
ans = Total[row[[id /@ Select[pts, First[#] == n - 1 &]]]];
Print[ans];
Print[N[ans, 65]];
