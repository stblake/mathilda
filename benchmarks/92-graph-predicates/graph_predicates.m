(* Experiment 92 -- Graph structural predicates and ordering.
   MEASURES the v0.183 batch in src/graph/ (graphprops.c, membership.c,
   acyclic.c): UndirectedGraphQ, EmptyGraphQ, CompleteGraphQ, BipartiteGraphQ,
   VertexQ, EdgeQ, AcyclicGraphQ, TreeGraphQ, TopologicalSort.

   Every graph is deterministic, so all three systems build the same object:
     dag   -- 100000 vertices, i->i+1 plus forward chords i->i+2+Mod[7i,50]
     tree  -- 100000-vertex binary-heap tree, i<->Floor[i/2]
     bip   -- 100000-vertex path plus odd chords i<->i+3 (even cycles only)
     kn    -- CompleteGraph[400], 79800 edges

   Baseline is networkx (pure Python), as in experiment 29.  Checks run on
   small graphs whose answer is unique, so the systems cannot silently time
   different computations: TopologicalSort's check graph is a chain, which has
   exactly one topological order. *)

Get["../harness.m"];

require[{"UndirectedGraphQ", "EmptyGraphQ", "CompleteGraphQ", "BipartiteGraphQ",
         "VertexQ", "EdgeQ", "AcyclicGraphQ", "TreeGraphQ", "TopologicalSort"}];

nv = 100000;
dag  = Graph[Range[nv], Join[Table[i -> i + 1, {i, 1, nv - 1}],
                             Table[i -> i + 2 + Mod[7 i, 50], {i, 1, nv - 52}]]];
tree = Graph[Range[nv], Table[i <-> Floor[i/2], {i, 2, nv}]];
bip  = Graph[Range[nv], Join[Table[i <-> i + 1, {i, 1, nv - 1}],
                             Table[i <-> i + 3, {i, 1, nv - 3, 2}]]];
kn   = CompleteGraph[400];

benchIf["AcyclicGraphQ, 100000-vertex DAG", "AcyclicGraphQ", AcyclicGraphQ[dag]];
checkIf["AcyclicGraphQ, 100000-vertex DAG", "AcyclicGraphQ",
  {AcyclicGraphQ[Graph[{1 -> 2, 2 -> 3, 1 -> 3}]], AcyclicGraphQ[Graph[{1 -> 2, 2 -> 3, 3 -> 1}]]}];

benchIf["TopologicalSort, 100000-vertex DAG", "TopologicalSort", TopologicalSort[dag]];
checkIf["TopologicalSort, 100000-vertex DAG", "TopologicalSort",
  TopologicalSort[Graph[{3 -> 1, 1 -> 2, 2 -> 4}]]];

benchIf["TreeGraphQ, 100000-vertex tree", "TreeGraphQ", TreeGraphQ[tree]];
checkIf["TreeGraphQ, 100000-vertex tree", "TreeGraphQ",
  {TreeGraphQ[Graph[{1 <-> 2, 1 <-> 3}]], TreeGraphQ[Graph[{1 <-> 2, 2 <-> 3, 3 <-> 1}]]}];

benchIf["AcyclicGraphQ, 100000-vertex tree", "AcyclicGraphQ", AcyclicGraphQ[tree]];
checkIf["AcyclicGraphQ, 100000-vertex tree", "AcyclicGraphQ",
  {AcyclicGraphQ[Graph[{1 <-> 2, 3 <-> 4}]], AcyclicGraphQ[Graph[{1 <-> 2, 2 <-> 3, 3 <-> 1}]]}];

benchIf["BipartiteGraphQ, 100000 vertices", "BipartiteGraphQ", BipartiteGraphQ[bip]];
checkIf["BipartiteGraphQ, 100000 vertices", "BipartiteGraphQ",
  {BipartiteGraphQ[Graph[{1 <-> 2, 2 <-> 3, 3 <-> 4, 4 <-> 1}]],
   BipartiteGraphQ[Graph[{1 <-> 2, 2 <-> 3, 3 <-> 1}]]}];

benchIf["CompleteGraphQ, K_400", "CompleteGraphQ", CompleteGraphQ[kn]];
checkIf["CompleteGraphQ, K_400", "CompleteGraphQ",
  {CompleteGraphQ[Graph[{1 <-> 2, 2 <-> 3, 3 <-> 1}]], CompleteGraphQ[Graph[{1 <-> 2, 2 <-> 3}]]}];

benchIf["UndirectedGraphQ, 100000-vertex tree", "UndirectedGraphQ", UndirectedGraphQ[tree]];
checkIf["UndirectedGraphQ, 100000-vertex tree", "UndirectedGraphQ",
  {UndirectedGraphQ[Graph[{1 <-> 2}]], UndirectedGraphQ[Graph[{1 -> 2}]]}];

benchIf["EmptyGraphQ, 100000-vertex tree", "EmptyGraphQ", EmptyGraphQ[tree]];
checkIf["EmptyGraphQ, 100000-vertex tree", "EmptyGraphQ",
  {EmptyGraphQ[Graph[{1, 2}, {}]], EmptyGraphQ[Graph[{1 <-> 2}]]}];

benchIf["VertexQ, 100000-vertex DAG", "VertexQ", VertexQ[dag, nv]];
checkIf["VertexQ, 100000-vertex DAG", "VertexQ",
  {VertexQ[Graph[{1 -> 2}], 2], VertexQ[Graph[{1 -> 2}], 3]}];

benchIf["EdgeQ, 100000-vertex DAG", "EdgeQ", EdgeQ[dag, (nv - 1) -> nv]];
checkIf["EdgeQ, 100000-vertex DAG", "EdgeQ",
  {EdgeQ[Graph[{1 -> 2}], 1 -> 2], EdgeQ[Graph[{1 -> 2}], 2 -> 1]}];

(* ---- COLD cases ----------------------------------------------------------
   The cases above repeat one call on one graph object, so after the warm-up
   both CASes answer from a per-graph cache (Mathematica's atomic Graph caches
   its properties; Mathilda's validated-graph memo does the same).  These build
   a FRESH graph before every timed call -- construction is untimed -- so no
   cache can help and the column measures the algorithm itself. *)

dagE  = EdgeList[dag];  treeE = EdgeList[tree];  bipE = EdgeList[bip];
newDag[]  := Graph[Range[nv], dagE];
newTree[] := Graph[Range[nv], treeE];
newBip[]  := Graph[Range[nv], bipE];
newKn[]   := Graph[Range[400], EdgeList[kn]];

coldBench[label_String, mk_, f_] := Module[{ts, gg},
  ts = Table[gg = mk[]; First[AbsoluteTiming[f[gg]]], {5}];
  Print["BENCH\t", label, "\t", ToString[Round[1000. Min[ts], 0.001]]]];

coldBench["AcyclicGraphQ, 100000-vertex DAG (cold)", newDag, AcyclicGraphQ];
check["AcyclicGraphQ, 100000-vertex DAG (cold)", AcyclicGraphQ[newDag[]]];
coldBench["TopologicalSort, 100000-vertex DAG (cold)", newDag, TopologicalSort];
check["TopologicalSort, 100000-vertex DAG (cold)", Length[TopologicalSort[newDag[]]]];
coldBench["TreeGraphQ, 100000-vertex tree (cold)", newTree, TreeGraphQ];
check["TreeGraphQ, 100000-vertex tree (cold)", TreeGraphQ[newTree[]]];
coldBench["AcyclicGraphQ, 100000-vertex tree (cold)", newTree, AcyclicGraphQ];
check["AcyclicGraphQ, 100000-vertex tree (cold)", AcyclicGraphQ[newTree[]]];
coldBench["BipartiteGraphQ, 100000 vertices (cold)", newBip, BipartiteGraphQ];
check["BipartiteGraphQ, 100000 vertices (cold)", BipartiteGraphQ[newBip[]]];
coldBench["CompleteGraphQ, K_400 (cold)", newKn, CompleteGraphQ];
check["CompleteGraphQ, K_400 (cold)", CompleteGraphQ[newKn[]]];
coldBench["UndirectedGraphQ, 100000-vertex tree (cold)", newTree, UndirectedGraphQ];
check["UndirectedGraphQ, 100000-vertex tree (cold)", UndirectedGraphQ[newTree[]]];
coldBench["EmptyGraphQ, 100000-vertex tree (cold)", newTree, EmptyGraphQ];
check["EmptyGraphQ, 100000-vertex tree (cold)", EmptyGraphQ[newTree[]]];
coldBench["VertexQ, 100000-vertex DAG (cold)", newDag, VertexQ[#, nv] &];
check["VertexQ, 100000-vertex DAG (cold)", VertexQ[newDag[], nv]];
coldBench["EdgeQ, 100000-vertex DAG (cold)", newDag, EdgeQ[#, (nv - 1) -> nv] &];
check["EdgeQ, 100000-vertex DAG (cold)", EdgeQ[newDag[], (nv - 1) -> nv]];
