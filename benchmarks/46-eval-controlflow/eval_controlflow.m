(* ==========================================================================
   Experiment 46 -- Evaluator: evaluation-heavy control flow
   ==========================================================================

       cd benchmarks/46-eval-controlflow
       ../../Mathilda -file eval_controlflow.m
       wolframscript  -file eval_controlflow.m

   WHAT IT MEASURES.  Tight interpreted loops with small bodies -- Do / If /
   Set, guarded DownValues, and Rational FixedPoint iteration.  With a trivial
   body, the per-iteration interpreter constant IS the whole cost, so this is
   the broadest confirmation of the same per-attempt constant that Case B (exp
   43) isolates.  Every loop keeps a scalar accumulator; no machine list is
   built, so nothing packs or compiles.
   ========================================================================== *)

Get["../harness.m"];

(* ---- 1. Do loop with a symbolic Set accumulator (100k iterations) ------ *)
(* Result is 100000 a; hammers the Set OwnValue fast path + Plus each step. *)
bench["Do[s=s+a, 100000]", Module[{s = 0}, Do[s = s + a, {100000}]; s];];
check["Do[s=s+a, 100000]", (Module[{s = 0}, Do[s = s + a, {100000}]; s]) /. a -> 1];

(* ---- 2. Do + If + PrimeQ over 10000 integers --------------------------- *)
bench["Do[If[PrimeQ,...], 10000]",
  Module[{s = 0}, Do[If[PrimeQ[k], s = s + k], {k, 10000}]; s];];
check["Do[If[PrimeQ,...], 10000]",
  Module[{s = 0}, Do[If[PrimeQ[k], s = s + k], {k, 10000}]; s]];

(* ---- 3. Newton's method on exact Rationals (6 steps) ------------------- *)
bench["FixedPoint Newton sqrt2, Rational",
  FixedPoint[Together[# - (#^2 - 2)/(2 #)] &, 1, 6];];
check["FixedPoint Newton sqrt2, Rational",
  FixedPoint[Together[# - (#^2 - 2)/(2 #)] &, 1, 6]];

(* ---- 4. Collatz length via guarded recursive DownValues ---------------- *)
Clear[col];
col[1] := 0;
col[n_ /; EvenQ[n]] := 1 + col[n/2];
col[n_ /; OddQ[n]] := 1 + col[3 n + 1];
bench["Collatz col[27] guarded rules", col[27];];
check["Collatz col[27] guarded rules", col[27]];
