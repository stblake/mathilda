(* Experiment 22 -- Small matrices: does LAPACK get reached at all?
   ROADMAP ITEM 3 -- "dispatch on 'is this a matrix of machine numbers', not 'is
   this already a buffer'", recorded as taking a 6x6 Inverse from 737 us to ~2 us
   and the Kalman row from 9.2 s to ~0.4 s.

   The mechanism: a small matrix sits BELOW the packing threshold, so it is a
   plain nested List, so the dispatch that looks for an existing buffer misses
   and the generic symbolic path runs.  Every row here is a small matrix
   operation repeated many times -- exactly the shape a Kalman filter or a
   ray tracer produces, and the shape a per-value threshold cannot see. *)

Get["../harness.m"];
Get["../data.m"];

require[{"Inverse", "Det", "LinearSolve", "Dot", "Eigenvalues", "Transpose"}];

m6 = N[spdIntMatrix[6]];
m3 = N[spdIntMatrix[3]];
v6 = N[Table[i, {i, 6}]];

bench["Inverse 6x6 x 2000", Do[Inverse[m6], {2000}];];
check["Inverse 6x6 x 2000", Round[10^6 Total[Flatten[Inverse[m6]]]]];

bench["Det 6x6 x 5000", Do[Det[m6], {5000}];];
check["Det 6x6 x 5000", Round[10^4 Det[m6]]];

bench["LinearSolve 6x6 x 2000", Do[LinearSolve[m6, v6], {2000}];];
check["LinearSolve 6x6 x 2000", Round[10^6 Total[LinearSolve[m6, v6]]]];

bench["Dot 6x6 x 6x6 x 10000", Do[m6.m6, {10000}];];
check["Dot 6x6 x 6x6 x 10000", Round[10^4 Total[Flatten[m6.m6]]]];

bench["Inverse 3x3 x 5000", Do[Inverse[m3], {5000}];];
check["Inverse 3x3 x 5000", Round[10^6 Total[Flatten[Inverse[m3]]]]];

bench["Eigenvalues 6x6 x 500", Do[Eigenvalues[m6], {500}];];
check["Eigenvalues 6x6 x 500", Round[10^4 Total[Sort[Re[Eigenvalues[m6]]]]]];
