(* ==========================================================================
   Experiment 41 -- Evaluator: pattern matching + backtracking
   ==========================================================================

       cd benchmarks/41-eval-backtracking
       ../../Mathilda -file eval_backtracking.m
       wolframscript  -file eval_backtracking.m

   WHAT IT MEASURES.  The structural matcher's backtracking throughput on
   flat-list sequence patterns with guards -- src/match.c enumeration order,
   fail-fast guard evaluation, and the "last sequence blank consumes the
   remainder" pruning rule.  These are the cases where Mathilda's matcher is
   expected to win decisively: it prunes boundary searches Mathematica
   brute-forces (see MATHILDA_PATTERN_MATCHER_PERFORMANCE.md cases D and M).

   Inputs are DETERMINISTIC (Table[Mod[...]]), identical in both engines, so the
   two columns match the same list.  The integer lists pack once and are
   materialized on first matcher touch (MatchQ / ReplaceList are not
   packed-aware) -- a one-time O(n) cost, negligible against the ms-scale
   backtracking work, and confirmed by the MATHILDA_NO_PACK timing-invariance
   check in the writeup.

   NOTE ON SYNTAX.  Mathilda's parser does not accept the postfix `..`/`...`
   Repeated operators; `__`/`___`/`_` and the `Repeated[]` long form all work,
   and are what these patterns use.
   ========================================================================== *)

Get["../harness.m"];

(* ---- 1. Case D: palindromic sequence blanks (perf-doc D) ---------------- *)
L200 = Table[Mod[k, 7] + 1, {k, 200}];
bench["Case D palindrome, 200 ints",
  MatchQ[L200, {___, a___, b___, c___, b___, a___, ___} /; Length[{c}] > 10];];
check["Case D palindrome, 200 ints",
  MatchQ[L200, {___, a___, b___, c___, b___, a___, ___} /; Length[{c}] > 10]];

(* ---- 2. Case M: distinct equal-sum contiguous blocks, N=60 (perf-doc M) - *)
LM60 = Table[Mod[k, 4] + 1, {k, 60}];
bench["Case M equal-sum blocks, N=60",
  ReplaceList[LM60, {___, b1__ /; Total[{b1}] == 5, b2__ /; Total[{b2}] == 5,
     b3__ /; Total[{b3}] == 5, ___} :> {b1, b2, b3}];];
check["Case M equal-sum blocks, N=60",
  Length[ReplaceList[LM60, {___, b1__ /; Total[{b1}] == 5,
     b2__ /; Total[{b2}] == 5, b3__ /; Total[{b3}] == 5, ___} :> {b1, b2, b3}]]];

(* ---- 3. Case M at N=80: the scaling axis (benchOnce -- WL is seconds) ---- *)
LM80 = Table[Mod[k, 4] + 1, {k, 80}];
benchOnce["Case M equal-sum blocks, N=80",
  ReplaceList[LM80, {___, b1__ /; Total[{b1}] == 5, b2__ /; Total[{b2}] == 5,
     b3__ /; Total[{b3}] == 5, ___} :> {b1, b2, b3}];];
check["Case M equal-sum blocks, N=80",
  Length[ReplaceList[LM80, {___, b1__ /; Total[{b1}] == 5,
     b2__ /; Total[{b2}] == 5, b3__ /; Total[{b3}] == 5, ___} :> {b1, b2, b3}]]];
