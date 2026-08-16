(* ==========================================================================
   Experiment 44 -- Evaluator: functional/structural ops over symbolic lists
   ==========================================================================

       cd benchmarks/44-eval-structural
       ../../Mathilda -file eval_structural.m
       wolframscript  -file eval_structural.m

   WHAT IT MEASURES.  The generic per-element dispatch loops behind Map / Fold /
   Count / Position over lists of SYMBOLIC function nodes f[k].  Because the
   elements have a symbolic head they can never pack (src/pack.c needs all
   machine int/real/bool), and Count / Position are not packed-aware anyway, so
   every element is dispatched through the evaluator + matcher -- exactly the
   generic loop this category is meant to time, not a machine kernel.
   ========================================================================== *)

Get["../harness.m"];

Clear[ff, sq, fp];
sq[ff[u_]] := u^2;              (* Map body: one DownValue per element *)
fp[acc_, ff[k_]] := acc + k;    (* Fold body: pattern-bind each element *)

S5 = Table[ff[k], {k, 5000}];   (* symbolic -> never packs *)

(* ---- 1. Map a DownValue over 5000 symbolic nodes ----------------------- *)
bench["Map[sq, 5000 f[k] nodes]", Map[sq, S5];];
check["Map[sq, 5000 f[k] nodes]", LeafCount[Map[sq, S5]]];

(* ---- 2. Fold a two-arg pattern over 5000 nodes ------------------------- *)
bench["Fold[fp, 0, 5000 f[k] nodes]", Fold[fp, 0, S5];];
check["Fold[fp, 0, 5000 f[k] nodes]", Fold[fp, 0, S5]];

(* ---- 3. Count with a PatternTest predicate ----------------------------- *)
bench["Count[5000 nodes, f[_?PrimeQ]]", Count[S5, ff[_?PrimeQ]];];
check["Count[5000 nodes, f[_?PrimeQ]]", Count[S5, ff[_?PrimeQ]]];

(* ---- 4. Position with a PatternTest predicate -------------------------- *)
S4 = Table[ff[k], {k, 4000}];
bench["Position[4000 nodes, f[_?EvenQ]]", Position[S4, ff[_?EvenQ]];];
check["Position[4000 nodes, f[_?EvenQ]]", Length[Position[S4, ff[_?EvenQ]]]];
