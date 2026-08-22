(* ==========================================================================
   Experiment 42 -- Evaluator: term rewriting to a fixed point (//.)
   ==========================================================================

       cd benchmarks/42-eval-rewrite
       ../../Mathilda -file eval_rewrite.m
       wolframscript  -file eval_rewrite.m

   WHAT IT MEASURES.  ReplaceRepeated (//.) drives the matcher + replace_bindings
   + the fixed-point convergence check once per pass, over many passes.  It is
   the purest test of "match an attempt, rewrite, re-scan" throughput.  The
   pairwise-fold row is the same one bench_eval.c tracks as its matcher-bound
   regression sentinel.

   All data is symbolic or a small integer list consumed by an unaware head
   (//. materializes any packed list), so nothing stays on a machine buffer.
   ========================================================================== *)

Get["../harness.m"];

(* ---- 1. Bubble sort by adjacent-swap rewrite --------------------------- *)
rev30 = Reverse[Range[30]];
bench["bubble sort //. , 30 elts",
  rev30 //. {p___, u_, v_, q___} /; u > v :> {p, v, u, q};];
check["bubble sort //. , 30 elts",
  (rev30 //. {p___, u_, v_, q___} /; u > v :> {p, v, u, q}) === Range[30]];

(* ---- 2. Pairwise left-fold to a scalar (bench_eval sentinel) ------------ *)
bench["pairwise fold //. , Range[400]",
  Range[400] //. {a_, b_, r___} :> {a + b, r};];
check["pairwise fold //. , Range[400]",
  First[Range[400] //. {a_, b_, r___} :> {a + b, r}]];

(* ---- 3. Same fold at 800: the O(n^2) rewrite axis ---------------------- *)
bench["pairwise fold //. , Range[800]",
  Range[800] //. {a_, b_, r___} :> {a + b, r};];
check["pairwise fold //. , Range[800]",
  First[Range[800] //. {a_, b_, r___} :> {a + b, r}]];

(* ---- 4. Peano addition by structural rewrite (purely symbolic) --------- *)
Clear[pl, sx];
peanoRules = {pl[0, n_] :> n, pl[sx[m_], n_] :> sx[pl[m, n]]};
n1 = Nest[sx, 0, 60]; n2 = Nest[sx, 0, 60];
bench["Peano add via //. , 60+60", pl[n1, n2] //. peanoRules;];
check["Peano add via //. , 60+60", (pl[n1, n2] //. peanoRules) === Nest[sx, 0, 120]];
