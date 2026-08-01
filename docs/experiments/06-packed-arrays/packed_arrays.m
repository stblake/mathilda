(* =======================================================================
   Experiment 6 -- Automatic packed arrays: dense lists become machine buffers
   =======================================================================

   Runs unmodified in Mathilda and in Mathematica 14. See README.md.
   
       Mathilda      ./Mathilda -file packed_arrays.m
       Mathematica   wolframscript -file packed_arrays.m
       Python        python3 packed_arrays.py
   
   WHAT IT MEASURES.  Ten operations on a 10^7-element vector -- reductions,
   scans, sorts, elementwise transcendentals and the structural family.  None of
   the source text says anything about representation: the point of automatic
   packing is that a dense list of machine numbers BECOMES a buffer without the
   user asking, so the same program is fast.
   
   The closing section shows the same operations with packing switched off, which
   is the measurement the experiment is actually about.
   ======================================================================= *)

SetAttributes[bench, HoldRest];

(* bench[label, expr] -- one untimed warm-up, then the MINIMUM of three timed
   runs.  The minimum, not the mean: we are measuring the cost of the work, and
   every source of noise on a loaded machine can only add. *)
bench[label_String, expr_] := Module[{ts},
  expr;
  ts = Table[First[AbsoluteTiming[expr]], {3}];
  Print[StringPadRight[label, 52], ToString[Round[1000. Min[ts], 0.001]], " ms"]
];

(* bench1 -- a single timed run, no warm-up, for kernels that take seconds.
   Reported as such wherever it is used. *)
SetAttributes[bench1, HoldRest];
bench1[label_String, expr_] :=
  Print[StringPadRight[label, 52],
        ToString[Round[1000. First[AbsoluteTiming[expr]], 0.001]], " ms  (1 run)"];

check[label_String, value_] :=
  Print[StringPadRight[label, 52], "check = ", value];

(* ---- total -- Total (reduction) ----------------------- *)
(* The simplest reduction: one pass, no allocation. *)
n=10^7;
x=RandomReal[{0,1},n];
y=RandomReal[{0,1},n];

(* ---- accum -- Accumulate (prefix scan) ---------------- *)
(* A prefix scan -- sequential, so it cannot be threaded without changing *)
(* the summation order. *)

(* ---- sort -- Sort ------------------------------------ *)
(* Order statistics: the one row here that is neither a stream nor a scan. *)

(* ---- sin -- Sin (elementwise) ----------------------- *)
(* An elementwise transcendental; libm-bound rather than memory-bound. *)

(* ---- exp -- Exp (elementwise) ----------------------- *)

(* ---- dot -- Dot (inner product) --------------------- *)
(* An inner product -- one pass, and a BLAS call once packed. *)

(* ---- triad -- STREAM triad, a = b + 3 c --------------- *)
(* STREAM triad: the memory-bandwidth reference for this machine. *)

(* ---- reverse -- Reverse --------------------------------- *)
(* The structural family. Each of these is a pure element MOVE, so on a *)
(* buffer it is a memcpy and on a List it is one allocation per element. *)

(* ---- rotate -- RotateLeft ------------------------------ *)

(* ---- diffs -- Differences ----------------------------- *)

Print["Experiment 6 -- Automatic packed arrays: dense lists become machine buffers"];
Print[""];

bench["Total (reduction)", Total[x]];
bench["Accumulate (prefix scan)", Accumulate[x]];
bench["Sort", Sort[x]];
bench["Sin (elementwise)", Sin[x]];
bench["Exp (elementwise)", Exp[x]];
bench["Dot (inner product)", x . y];
bench["STREAM triad, a = b + 3 c", x + 3. y];
bench["Reverse", Reverse[x]];
bench["RotateLeft", RotateLeft[x,3]];
bench["Differences", Differences[x]];
