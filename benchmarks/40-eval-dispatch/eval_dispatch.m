(* ==========================================================================
   Experiment 40 -- Evaluator: recursive DownValue dispatch
   ==========================================================================

       cd benchmarks/40-eval-dispatch
       ../../Mathilda -file eval_dispatch.m
       wolframscript  -file eval_dispatch.m

   WHAT IT MEASURES.  The core evaluator's per-call dispatch cost through user
   DownValues -- pattern lookup, the DownValue dispatch index, MatchEnv
   allocation, and the per-step node reassembly in evaluate_step().  Every
   workload is scalar recursion or a scalar loop that NEVER builds a machine
   list, so nothing is diverted to the packed-array kernels (src/pack.c) or to
   Compile[].  This is the recursive expression-tree evaluator, nothing else.

   WHY IT IS HERE.  The pattern matcher already beats Mathematica 4x-135000x
   (MATHILDA_PATTERN_MATCHER_PERFORMANCE.md).  The one open evaluator gap is the
   per-attempt constant factor.  Naive fib / Ackermann are the cleanest probes
   of that constant: exponential/deep dispatch with a trivial body.
   ========================================================================== *)

Get["../harness.m"];

(* ---- 1. Naive Fibonacci: exponential DownValue dispatch ----------------- *)
Clear[fib]; fib[0] = 0; fib[1] = 1; fib[n_] := fib[n - 1] + fib[n - 2];
bench["fib[28] naive recursion", fib[28];];
check["fib[28] naive recursion", fib[28]];

(* ---- 2. Ackermann: deeply nested dispatch ------------------------------- *)
Clear[ack]; ack[0, n_] := n + 1; ack[m_, 0] := ack[m - 1, 1];
ack[m_, n_] := ack[m - 1, ack[m, n - 1]];
bench["ack[3,4] nested recursion", ack[3, 4];];
check["ack[3,4] nested recursion", ack[3, 4]];

(* ---- 3. Mutual recursion: linear dispatch depth ------------------------- *)
Clear[evenR, oddR]; evenR[0] = True; oddR[0] = False;
evenR[n_] := oddR[n - 1]; oddR[n_] := evenR[n - 1];
bench["mutual recursion depth 800", evenR[800];];
check["mutual recursion depth 800", evenR[800]];

(* ---- 4. Nest over a DownValue: tight rewrite loop ----------------------- *)
(* Nest is iterative (no C-stack recursion); each step is one DownValue match
   + one Plus.  Scalar integer accumulator, never a list. *)
Clear[incr]; incr[x_] := x + 1;
bench["Nest[incr, 0, 50000]", Nest[incr, 0, 50000];];
check["Nest[incr, 0, 50000]", Nest[incr, 0, 50000]];
