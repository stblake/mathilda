(* ==========================================================================
   Experiment 47 -- Evaluator: memoization and DownValue-table growth
   ==========================================================================

       cd benchmarks/47-eval-memoization
       ../../Mathilda -file eval_memoization.m
       wolframscript  -file eval_memoization.m

   WHAT IT MEASURES.  The cost of BUILDING and then dispatching through a
   DownValue table that grows to hundreds or thousands of entries: the
   `f[n_] := f[n] = ...` memoization idiom (a fresh DownValue per distinct
   argument, inserted into the dispatch index) and bulk `f[k] = ...` fact
   assertion.  Pure integer/bignum recursion and a scalar sum -- no machine
   list is ever built, so nothing packs (src/pack.c) or compiles.

   THE MEMO TRAP (read before editing).  Each timed expression CLEARS and
   REBUILDS its table from scratch, because a memo table persists across the
   harness's warm-up + 3 timed reps: timing a bare `f[1000]` would let reps 2-4
   measure an O(1) cached lookup, not the construction.  Clearing inside the
   timed body makes every rep pay the full build, which is the quantity of
   interest.  `$RecursionLimit` is raised because memoized fib[1000] descends
   1000 frames on the first (uncached) call; both engines honor the setting, so
   the comparison stays fair.
   ========================================================================== *)

Get["../harness.m"];

$RecursionLimit = 100000;

(* ---- 1. Memoized Fibonacci: 1000 DownValues, built then looked up ------- *)
bench["memo fib[1000]",
  (Clear[fibm]; fibm[0] = 0; fibm[1] = 1;
   fibm[n_] := fibm[n] = fibm[n - 1] + fibm[n - 2]; fibm[1000])];
check["memo fib[1000]",
  (Clear[fibm]; fibm[0] = 0; fibm[1] = 1;
   fibm[n_] := fibm[n] = fibm[n - 1] + fibm[n - 2]; fibm[1000])];

(* ---- 2. Memoized Pascal recurrence: 2-arg memo, ~2500 entries ----------- *)
bench["memo binom[100,50]",
  (Clear[bnm]; bnm[n_, 0] := 1; bnm[n_, n_] := 1;
   bnm[n_, k_] := bnm[n, k] = bnm[n - 1, k - 1] + bnm[n - 1, k]; bnm[100, 50])];
check["memo binom[100,50]",
  (Clear[bnm]; bnm[n_, 0] := 1; bnm[n_, n_] := 1;
   bnm[n_, k_] := bnm[n, k] = bnm[n - 1, k - 1] + bnm[n - 1, k]; bnm[100, 50])];

(* ---- 3. 5000 flat facts f[k]=k^2, then sum them (dispatch index at scale) *)
bench["build 5000 facts, sum",
  (Clear[fct]; Do[fct[k] = k^2, {k, 5000}]; Sum[fct[k], {k, 5000}])];
check["build 5000 facts, sum",
  (Clear[fct]; Do[fct[k] = k^2, {k, 5000}]; Sum[fct[k], {k, 5000}])];
