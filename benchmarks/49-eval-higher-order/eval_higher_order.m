(* ==========================================================================
   Experiment 49 -- Evaluator: higher-order / restructuring ops on symbolic data
   ==========================================================================

       cd benchmarks/49-eval-higher-order
       ../../Mathilda -file eval_higher_order.m
       wolframscript  -file eval_higher_order.m

   WHAT IT MEASURES.  Outer / Sort / FoldList / Apply -- the operators that
   build or reshape large expression trees element by element through the
   generic evaluator.  Every operand carries a symbolic head (aa[i], cf[..] q^..)
   so nothing can pack (src/pack.c needs all-machine data) and the work stays on
   the interpreter: Outer builds a 150x150 grid of nodes, Sort exercises the
   canonical-order comparator (expr_compare) over 3000 symbolic terms.

   NOTE ON THE SORT CHECK.  Mathilda and Mathematica do NOT agree on the
   canonical ORDER of arbitrary symbolic terms, so First/Last[Sort[...]] would
   CHECK-FAIL on a difference that is not a bug.  The check is instead the
   order-independent invariant Total[Sort[l]] === Total[l] -- true in both
   systems, and still falsified by a Sort that dropped or duplicated a term.
   ========================================================================== *)

Get["../harness.m"];

symL = Table[cf[Mod[k^2, 101]] q^Mod[k, 7], {k, 3000}];   (* symbolic -> never packs *)

(* ---- 1. Outer product: 150x150 = 22500 symbolic nodes ------------------ *)
bench["Outer 150x150 symbolic",
  Outer[gg, Table[aa[i], {i, 150}], Table[bb[j], {j, 150}]];];
check["Outer 150x150 symbolic",
  LeafCount[Outer[gg, Table[aa[i], {i, 150}], Table[bb[j], {j, 150}]]]];

(* ---- 2. Sort 3000 symbolic terms (canonical-order comparator) ---------- *)
bench["Sort 3000 symbolic", Sort[symL];];
check["Sort 3000 symbolic", Total[Sort[symL]] === Total[symL]];

(* ---- 3. FoldList: build 6000 nested symbolic accumulations ------------- *)
bench["FoldList 6000 symbolic", FoldList[hh, x0, Table[cc[k], {k, 6000}]];];
check["FoldList 6000 symbolic", Depth[FoldList[hh, x0, Table[cc[k], {k, 6000}]]]];

(* ---- 4. Apply at level 1: restructure 5000 triples --------------------- *)
bench["Apply Plus @ 5000 triples",
  Apply[Plus, Table[{aa[k], bb[k], cc[k]}, {k, 5000}], {1}];];
check["Apply Plus @ 5000 triples",
  LeafCount[Apply[Plus, Table[{aa[k], bb[k], cc[k]}, {k, 5000}], {1}]]];
