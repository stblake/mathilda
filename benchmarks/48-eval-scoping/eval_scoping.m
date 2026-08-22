(* ==========================================================================
   Experiment 48 -- Evaluator: lexical scoping and pure-function application
   ==========================================================================

       cd benchmarks/48-eval-scoping
       ../../Mathilda -file eval_scoping.m
       wolframscript  -file eval_scoping.m

   WHAT IT MEASURES.  The constructs that create and tear down bindings on every
   call: Module (fresh renamed locals per instantiation), With (constant
   substitution into a body), and pure-function / Slot application.  These run
   the IDENTICAL evaluator work on both engines -- there is no native kernel to
   divert to, unlike a builtin such as D[] -- so the comparison stays on the
   general evaluator, and every operand is symbolic or a scalar accumulator so
   nothing packs.  `$RecursionLimit` is raised for the depth-1500 Module
   recursion; both engines honor the setting.
   ========================================================================== *)

Get["../harness.m"];

$RecursionLimit = 100000;

(* ---- 1. Module-bodied recursion: fresh locals at every frame ----------- *)
Clear[gs]; gs[0] := 0; gs[n_] := Module[{p}, p = gs[n - 1]; p + n];
bench["Module recursion gs[1500]", gs[1500];];
check["Module recursion gs[1500]", gs[1500]];

(* ---- 2. With substitution inside a 30000-iteration loop ---------------- *)
bench["With loop, 30000",
  Module[{s = 0}, Do[With[{u = k}, s = s + u^2], {k, 30000}]; s];];
check["With loop, 30000",
  Module[{s = 0}, Do[With[{u = k}, s = s + u^2], {k, 30000}]; s]];

(* ---- 3. Module symbol generation: 30000 fresh lexical scopes ----------- *)
bench["Module gen, 30000 scopes",
  LeafCount[Table[Module[{a, b, c}, a + b + c], {30000}]];];
check["Module gen, 30000 scopes",
  LeafCount[Table[Module[{a, b, c}, a + b + c], {30000}]]];

(* ---- 4. Pure-function (Slot) application folded 6000 deep (symbolic) ---- *)
bench["Fold pure-function, 6000",
  Fold[(cs[#1, #2]) &, x0, Table[a[k], {k, 6000}]];];
check["Fold pure-function, 6000",
  Depth[Fold[(cs[#1, #2]) &, x0, Table[a[k], {k, 6000}]]]];
