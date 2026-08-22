(* ==========================================================================
   Experiment 45 -- Evaluator: expression construction + traversal
   ==========================================================================

       cd benchmarks/45-eval-construction
       ../../Mathilda -file eval_construction.m
       wolframscript  -file eval_construction.m

   WHAT IT MEASURES.  Node allocation and whole-tree traversal/rebuild on
   DISTINCT-node symbolic structures.  Distinct nodes matter: a self-similar
   Nest[{#,#}&, ...] collapses to a shared DAG under Mathilda's copy-on-write
   refcounting and would measure sharing rather than construction, so this uses
   a wide list of distinct node[k,k^2,k^3] terms (never packs -- symbolic head)
   and a deep single-arg nest.
   ========================================================================== *)

Get["../harness.m"];

Clear[node, gg];
wide = Table[node[k, k^2, k^3], {k, 15000}];   (* built once, outside timing *)

(* ---- 1. Construct 15000 distinct symbolic nodes ------------------------ *)
bench["build 15000 node[k,k^2,k^3]", Table[node[k, k^2, k^3], {k, 15000}];];
check["build 15000 node[k,k^2,k^3]", LeafCount[Table[node[k, k^2, k^3], {k, 15000}]]];

(* ---- 2. Full-tree leaf traversal -------------------------------------- *)
bench["LeafCount over 60001 leaves", LeafCount[wide];];
check["LeafCount over 60001 leaves", LeafCount[wide]];

(* ---- 3. Whole-tree rewrite/rebuild via ReplaceAll --------------------- *)
bench["ReplaceAll rebuild, 15000 nodes", wide /. node[a_, b_, c_] :> a + b + c;];
check["ReplaceAll rebuild, 15000 nodes", Total[wide /. node[a_, b_, c_] :> a + b + c]];

(* ---- 4. Deep single-arg nest: 900-level traversal --------------------- *)
bench["Depth[Nest[gg, x0, 900]]", Depth[Nest[gg, x0, 900]];];
check["Depth[Nest[gg, x0, 900]]", Depth[Nest[gg, x0, 900]]];
