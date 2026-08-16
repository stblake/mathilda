(* ==========================================================================
   Experiment 43 -- Evaluator: per-node arithmetic on symbolic/Rational data
   ==========================================================================

       cd benchmarks/43-eval-arithmetic
       ../../Mathilda -file eval_arithmetic.m
       wolframscript  -file eval_arithmetic.m

   WHAT IT MEASURES.  The interpreter's per-attempt constant factor on ordinary
   Length / Plus / Times / Equal, plus the Orderless canonicalization of large
   symbolic sums.  This is the flagship category: it isolates exactly the
   ~22us-vs-~14us-per-attempt gap that makes "Case B" 1.6x slower than
   Mathematica (MATHILDA_PATTERN_MATCHER_PERFORMANCE.md section B), and asks
   whether that constant shows up beyond the one synthetic case.

   Data is Range/symbolic/Rational: the Case-B lists are consumed by MatchQ
   (unaware -> materialized once), the sums are symbolic c[k] x (never pack), the
   Rational total is k/3 (Rationals never pack -- this is the same trick
   bench_eval.c uses for its calibration row).
   ========================================================================== *)

Get["../harness.m"];

(* ---- 1. Case B: typed sequence blanks + linear length guard ------------- *)
(* 580 guard evaluations; the matcher walks the identical search space as WL,
   so the difference is purely the per-attempt evaluator constant. *)
r500 = Range[500];
bench["Case B typed, |x|+|y|+|z| guard",
  MatchQ[r500, {xx : Repeated[_Integer], yy : Repeated[_Integer],
     zz : Repeated[_Integer]} /; 3 Length[{xx}] + 5 Length[{yy}] == Length[{zz}]];];
check["Case B typed, |x|+|y|+|z| guard",
  MatchQ[r500, {xx : Repeated[_Integer], yy : Repeated[_Integer],
     zz : Repeated[_Integer]} /; 3 Length[{xx}] + 5 Length[{yy}] == Length[{zz}]]];

(* ---- 2. Case B untyped: isolates the per-element _Integer re-check cost -- *)
bench["Case B untyped __ , same guard",
  MatchQ[r500, {xx__, yy__, zz__} /;
     3 Length[{xx}] + 5 Length[{yy}] == Length[{zz}]];];
check["Case B untyped __ , same guard",
  MatchQ[r500, {xx__, yy__, zz__} /;
     3 Length[{xx}] + 5 Length[{yy}] == Length[{zz}]]];

(* ---- 3. Orderless canonicalization of a 2000-term symbolic sum --------- *)
bench["Length[Plus @@ 2000 symbolic terms]",
  Length[Plus @@ Table[c[k] x, {k, 2000}]];];
check["Length[Plus @@ 2000 symbolic terms]",
  Length[Plus @@ Table[c[k] x, {k, 2000}]]];

(* ---- 4. Rational accumulation, 20000 terms (never packs) --------------- *)
bench["Total[Table[k/3, 20000]]", Total[Table[k/3, {k, 20000}]];];
check["Total[Table[k/3, 20000]]", Total[Table[k/3, {k, 20000}]]];

(* ---- 5. Nested symbolic Plus collapse, 4000 deep ----------------------- *)
bench["Nest[#+a&, x, 4000] then a->0", Nest[# + a &, x, 4000] /. a -> 0;];
check["Nest[#+a&, x, 4000] then a->0", (Nest[# + a &, x, 4000] /. a -> 0) === x];
