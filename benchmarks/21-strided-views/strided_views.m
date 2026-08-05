(* Experiment 21 -- Strided views: Transpose and sliding windows.
   ROADMAP ITEM 2 in docs/experiments/README.md -- "a stride vector on NDArray
   that consumers honour, giving a free Transpose and a free sliding window",
   valued there at 1.36x on MLP training and ~3x on k-mer counting.

   This experiment exists so that item is RE-MEASURED every week rather than
   once.  NumPy has strided views, so its Transpose is O(1) and Mathilda's is a
   copy: the gap on the Transpose row is the value of the item, directly.

   NOTE the Python column deliberately materialises with .copy() where the
   Wolfram semantics materialise -- see the .py docstring.  Without that the
   comparison is a memory copy against a pointer change. *)

Get["../harness.m"];
Get["../data.m"];

require[{"Transpose", "Partition", "Take", "Part", "Reverse", "Flatten",
         "ArrayReshape", "RotateLeft"}];

m = rand01[{2000, 2000}];
v = rand01[{2000000}];

bench["Transpose 2000x2000", Transpose[m];];
check["Transpose 2000x2000", Total[Flatten[Transpose[{{1, 2, 3}, {4, 5, 6}}]]]];

bench["Transpose then Dot (fused?)", Transpose[m].m, 1];
check["Transpose then Dot (fused?)",
  Round[10^4 Total[Flatten[Transpose[N[spdIntMatrix[6]]].N[spdIntMatrix[6]]]]]];

(* A sliding window: with strided views this is free, without it is a copy of
   n*w elements. *)
bench["Partition window 8, offset 1", Partition[v, 8, 1];];
check["Partition window 8, offset 1",
  Length[Partition[Range[10], 8, 1]]];

benchIf["ArrayReshape 2x10^6 to 1000x2000", "ArrayReshape",
  ArrayReshape[v, {1000, 2000}];];
checkIf["ArrayReshape 2x10^6 to 1000x2000", "ArrayReshape",
  Total[Flatten[ArrayReshape[Range[6], {2, 3}]]]];

bench["Take rows 1;;1000 of 2000x2000", Take[m, 1000];];
check["Take rows 1;;1000 of 2000x2000", Length[Take[m, 1000]]];

bench["column slice m[[All, 1]]", m[[All, 1]];];
check["column slice m[[All, 1]]",
  Total[{{1, 2}, {3, 4}}[[All, 1]]]];
