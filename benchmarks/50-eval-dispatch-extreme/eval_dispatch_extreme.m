(* ==========================================================================
   Experiment 50 -- Evaluator: dispatch / backtracking / rewrite at extreme size
   ==========================================================================

       cd benchmarks/50-eval-dispatch-extreme
       ../../Mathilda -file eval_dispatch_extreme.m
       wolframscript  -file eval_dispatch_extreme.m

   WHAT IT MEASURES.  The same evaluator machinery as experiments 40-42, pushed
   well past their sizes to confirm the per-attempt constant and the
   backtracking advantage hold under stress rather than only at the calibration
   point: naive Fibonacci at 31 (~4x the calls of fib[28]), Ackermann ack[3,6]
   (~20x ack[3,4]), the O(n^2) pairwise fold at Range[1200], and the equal-sum
   partition search at N=90.  All scalar recursion / symbolic rewriting -- no
   packing, no compilation.  The two heaviest single runs use benchOnce.
   `$RecursionLimit` is raised for ack[3,6], which descends ~500 frames.
   ========================================================================== *)

Get["../harness.m"];

$RecursionLimit = 100000;

(* ---- 1. Naive Fibonacci at 31 (single timed run) ----------------------- *)
Clear[fibn]; fibn[0] = 0; fibn[1] = 1; fibn[n_] := fibn[n - 1] + fibn[n - 2];
benchOnce["naive fib[31]", fibn[31];];
check["naive fib[31]", fibn[31]];

(* ---- 2. Ackermann ack[3,6] (harder than ack[3,4]) ---------------------- *)
Clear[ackn]; ackn[0, n_] := n + 1; ackn[m_, 0] := ackn[m - 1, 1];
ackn[m_, n_] := ackn[m - 1, ackn[m, n - 1]];
bench["ack[3,6] nested recursion", ackn[3, 6];];
check["ack[3,6] nested recursion", ackn[3, 6]];

(* ---- 3. Pairwise left-fold //. at Range[1200] (O(n^2) rewrite axis) ----- *)
bench["pairwise fold //. , Range[1200]",
  Range[1200] //. {a_, b_, r___} :> {a + b, r};];
check["pairwise fold //. , Range[1200]",
  First[Range[1200] //. {a_, b_, r___} :> {a + b, r}]];

(* ---- 4. Case M equal-sum blocks at N=90 (backtracking; single run) ------ *)
LM90 = Table[Mod[k, 4] + 1, {k, 90}];
benchOnce["Case M equal-sum blocks, N=90",
  ReplaceList[LM90, {___, b1__ /; Total[{b1}] == 5, b2__ /; Total[{b2}] == 5,
     b3__ /; Total[{b3}] == 5, ___} :> {b1, b2, b3}];];
check["Case M equal-sum blocks, N=90",
  Length[ReplaceList[LM90, {___, b1__ /; Total[{b1}] == 5,
     b2__ /; Total[{b2}] == 5, b3__ /; Total[{b3}] == 5, ___} :> {b1, b2, b3}]]];
