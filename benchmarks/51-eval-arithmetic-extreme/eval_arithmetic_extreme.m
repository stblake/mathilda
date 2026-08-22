(* ==========================================================================
   Experiment 51 -- Evaluator: arithmetic / control / bignum at extreme size
   ==========================================================================

       cd benchmarks/51-eval-arithmetic-extreme
       ../../Mathilda -file eval_arithmetic_extreme.m
       wolframscript  -file eval_arithmetic_extreme.m

   WHAT IT MEASURES.  Experiments 43 and 46, scaled up, plus a bignum axis:
   Orderless canonicalization of a 5000-term symbolic sum, a 300000-iteration
   Do loop with a symbolic Set accumulator, a 60000-term Rational total, and
   12000! evaluated as 12000 successive Times through the evaluation loop (each a
   GMP multiply on a growing operand).  Symbolic terms and Rationals never pack;
   the factorial keeps a scalar bignum accumulator, so all four stay on the
   interpreter rather than a machine kernel or Compile[].
   ========================================================================== *)

Get["../harness.m"];

(* ---- 1. Orderless canonicalization of a 5000-term symbolic sum --------- *)
bench["Length[Plus @@ 5000 symbolic terms]",
  Length[Plus @@ Table[cc[k] x, {k, 5000}]];];
check["Length[Plus @@ 5000 symbolic terms]",
  Length[Plus @@ Table[cc[k] x, {k, 5000}]]];

(* ---- 2. Do loop with a symbolic Set accumulator, 300000 iterations ----- *)
bench["Do[s=s+a, 300000]", Module[{s = 0}, Do[s = s + a, {300000}]; s];];
check["Do[s=s+a, 300000]", (Module[{s = 0}, Do[s = s + a, {300000}]; s]) /. a -> 1];

(* ---- 3. Rational accumulation, 60000 terms (never packs) --------------- *)
bench["Total[Table[k/3, 60000]]", Total[Table[k/3, {k, 60000}]];];
check["Total[Table[k/3, 60000]]", Total[Table[k/3, {k, 60000}]]];

(* ---- 4. 12000! via Fold[Times]: 12000 bignum multiplies through the loop - *)
bench["Fold[Times, 1, Range[12000]]", Fold[Times, 1, Range[12000]];];
check["Fold[Times, 1, Range[12000]]", Fold[Times, 1, Range[12000]] === 12000!];
