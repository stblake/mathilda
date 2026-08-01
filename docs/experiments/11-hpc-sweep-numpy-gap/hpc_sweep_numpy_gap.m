(* ==========================================================================
   Experiment 11 -- Fourth sweep: closing the distance to NumPy
   ==========================================================================

   Runs unmodified in Mathilda and in Mathematica 14. See README.md.

       Mathilda      ./Mathilda -file hpc_sweep_numpy_gap.m
       Mathematica   wolframscript -file hpc_sweep_numpy_gap.m
       Python        python3 hpc_sweep_numpy_gap.py

   WHAT IT MEASURES.  This sweep did NOT start from applications.  The third
   sweep had left two application kernels far behind NumPy and neither was a
   linear-algebra row, so the hypothesis was that the remaining cost was in
   PRIMITIVES -- and specifically in two categories no earlier sweep had
   probed: STRUCTURAL operations (an O(1) or memcpy operation walking a boxed
   list) and SCANS.

   So this file is a probe table, not an application: one operation per line,
   each against its NumPy equivalent, at 10^6 float64.

   THE CONTROL IS THE POINT.  The last two rows -- a plain product and a Log
   -- were already at the memory floor before this sweep, which is what says
   that nothing above them is about arithmetic.  Read the table from the
   bottom up: everything expensive is structural or a scan.
   ========================================================================== *)

SetAttributes[bench, HoldRest];

(* bench[label, expr] -- one untimed warm-up, then the MINIMUM of three timed
   runs.  The minimum, not the mean: we are measuring the cost of the work, and
   every source of noise on a loaded machine can only add. *)
bench[label_String, expr_] := Module[{ts},
  expr;
  ts = Table[First[AbsoluteTiming[expr]], {3}];
  Print[StringPadRight[label, 52], ToString[Round[1000. Min[ts], 0.001]], " ms"]
];

check[label_String, value_] :=
  Print[StringPadRight[label, 52], "check = ", value];


n = 10^6;
v = RandomReal[{0, 1}, n];
k5 = RandomReal[{0, 1}, 5];
im = RandomReal[{0, 1}, {1024, 1024}];
k55 = RandomReal[{0, 1}, {5, 5}];

Print["Experiment 11 -- structural operations, scans and convolution"];
Print[""];
Print["-- the structural family: pure element MOVES --"];

bench["First[v]   (an O(1) element read)", First[v]];
bench["Last[v]",                           Last[v]];
bench["Most[v]",                           Most[v]];
bench["Rest[v]",                           Rest[v]];
bench["Drop[v, 250]",                      Drop[v, 250]];
bench["Part[v, 250 ;; -250]",              v[[250 ;; -250]]];
bench["Reverse[v]",                        Reverse[v]];
bench["RotateLeft[v, 3]",                  RotateLeft[v, 3]];

Print[""];
Print["-- scans --"];

bench["Accumulate[v]",                     Accumulate[v]];
bench["Differences[v]",                    Differences[v]];
bench["FoldList[Max, First[v], Rest[v]]",  FoldList[Max, First[v], Rest[v]]];
bench["FoldList[Max[#1,#2] &, ...]",       FoldList[Max[#1, #2] &, First[v], Rest[v]]];
bench["FoldList[Plus, 0., v]",             FoldList[Plus, 0., v]];
bench["FoldList[0.98 #1 + 0.02 #2 &, ...] (EMA)",
      FoldList[0.98 #1 + 0.02 #2 &, 0., v]];

Print[""];
Print["-- convolution --"];

bench["ListConvolve[k5, v]",               ListConvolve[k5, v]];
bench["ListCorrelate[k55, im]  (1024^2)",  ListCorrelate[k55, im]];

Print[""];
Print["-- clipping --"];
bench["Clip[v, {0.2, 0.8}]",               Clip[v, {0.2, 0.8}]];

Print[""];
Print["-- THE CONTROL: already at the memory floor before this sweep --"];
bench["v v",                                v v];
bench["Log[v]",                             Log[v]];
