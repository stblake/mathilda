(* Experiment 26 -- Structural operations: block moves and gathers.
   ROADMAP ITEM 9 -- "thread the block moves and the gather: Reverse, RotateLeft,
   Part", valued at 2.9-3.8x and noted as latency-bound, so threading is close to
   linear (experiments 6, 12, 19).

   These are pure memory motion -- no arithmetic at all -- so the comparison
   against NumPy is a comparison against memory bandwidth, which is the honest
   ceiling.  The gather row (x[[idx]]) is the one experiment 12 found costing
   PageRank 14.1 s: the INDEX arrived packed and the selector only accepted a
   plain List, so both operands were materialised. *)

Get["../harness.m"];
Get["../data.m"];

require[{"Reverse", "RotateLeft", "RotateRight", "Part", "Join", "Take",
         "Drop", "Riffle", "Flatten", "Sort", "Union"}];

n = 4000000;
v = rand01[{n}];
idx = RandomInteger[{1, n}, {n}];

bench["Reverse 4x10^6", Reverse[v];];
check["Reverse 4x10^6", Total[Reverse[Range[10]]]];

bench["RotateLeft 4x10^6 by 1000", RotateLeft[v, 1000];];
check["RotateLeft 4x10^6 by 1000", Total[RotateLeft[Range[10], 3]]];

(* THE GATHER.  A large operand meeting a large operand -- the case a
   per-value packing threshold cannot resolve. *)
bench["gather v[[idx]], 4x10^6", v[[idx]];];
check["gather v[[idx]], 4x10^6", Total[Range[10][[{2, 4, 6}]]]];

bench["Join two 2x10^6", Join[Take[v, 2000000], Take[v, -2000000]];];
check["Join two 2x10^6", Total[Join[Range[5], Range[5]]]];

bench["Take first half of 4x10^6", Take[v, 2000000];];
check["Take first half of 4x10^6", Total[Take[Range[10], 5]]];

bench["Sort 4x10^6", Sort[v], 1];
check["Sort 4x10^6", Total[Take[Sort[lcgInts[2000, 50000]], 10]]];

bench["Union of 4x10^6 integers", Union[idx], 1];
check["Union of 4x10^6 integers", Length[Union[{3, 1, 2, 3, 1}]]];
