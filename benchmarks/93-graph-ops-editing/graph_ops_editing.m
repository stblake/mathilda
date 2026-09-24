(* Experiment 93 -- Graph editing, transforms, set operations, Euler tours.
   MEASURES the "ops" batch in src/graph/gops_*.c: VertexDelete, VertexAdd,
   EdgeAdd, EdgeDelete, Subgraph, IndexGraph, GraphUnion, GraphIntersection,
   GraphDifference, GraphDisjointUnion, GraphComplement, ReverseGraph,
   UndirectedGraph, DirectedGraph, LineGraph, FindEulerianCycle, FindPath,
   FindCycle.

   Every graph is deterministic, so all three systems build the same object:
     g   -- 100000 vertices: path i<->i+1 plus chords i<->i+2+Mod[7i,50]
            (199947 edges)
     g2  -- 100000 vertices: i<->i+3 (99997 edges)
     dg  -- g with every edge directed i->j (199947 arcs)
     eu  -- 50000-vertex circulant i<->i+1, i<->i+2 (mod n): 4-regular,
            connected, so Eulerian with 100000 edges
     c1k -- CycleGraph[1000]; its complement has 498500 edges
     h, h2 -- g and g2 on 2000 vertices, for GraphIntersection/GraphDifference:
            Mathematica 15's are roughly quadratic (about 60 s at 5000
            vertices, 230 s at 10000), so the 10^5 graphs would take hours

   Mathematica and Mathilda graphs are immutable, so each networkx case copies
   before mutating (g.copy() + remove_nodes_from, ...): the three columns time
   "produce the edited graph".  Checks are vertex/edge counts (and, for
   FindEulerianCycle, the tour length plus a validity bit computed by the
   system itself), which every system must agree on.

   COLD cases build a fresh graph (both graphs, for the set operations) before
   every timed call (construction is untimed), so no per-graph cache -- Mathematica's Graph object caches, and
   Mathilda's validated-graph memo -- can help. *)

Get["../harness.m"];

require[{"VertexDelete", "VertexAdd", "EdgeAdd", "EdgeDelete", "Subgraph",
         "IndexGraph", "GraphUnion", "GraphIntersection", "GraphDifference",
         "GraphDisjointUnion", "GraphComplement", "ReverseGraph",
         "UndirectedGraph", "DirectedGraph", "LineGraph", "FindEulerianCycle",
         "FindPath", "FindCycle"}];

nv = 100000;
gE  = Join[Table[i <-> i + 1, {i, 1, nv - 1}], Table[i <-> i + 2 + Mod[7 i, 50], {i, 1, nv - 52}]];
g2E = Table[i <-> i + 3, {i, 1, nv - 3}];
dgE = Join[Table[i -> i + 1, {i, 1, nv - 1}], Table[i -> i + 2 + Mod[7 i, 50], {i, 1, nv - 52}]];
ne = 50000;
euE = Join[Table[i <-> Mod[i, ne] + 1, {i, 1, ne}], Table[i <-> Mod[i + 1, ne] + 1, {i, 1, ne}]];
g   = Graph[Range[nv], gE];
g2  = Graph[Range[nv], g2E];
dg  = Graph[Range[nv], dgE];
eu  = Graph[Range[ne], euE];
c1k = CycleGraph[1000];
ns = 2000;       (* Mathematica's GraphIntersection/GraphDifference are ~quadratic *)
hE  = Join[Table[i <-> i + 1, {i, 1, ns - 1}], Table[i <-> i + 2 + Mod[7 i, 50], {i, 1, ns - 52}]];
h2E = Table[i <-> i + 3, {i, 1, ns - 3}];
h   = Graph[Range[ns], hE];
h2  = Graph[Range[ns], h2E];
delV = Range[1, nv, 100];                               (* 1000 vertices *)
addV = Range[nv + 1, nv + 10000];                       (* 10000 new vertices *)
addE = Table[i <-> i + 1000, {i, 1, 10000}];            (* 10000 new edges *)
delE = Table[i <-> i + 1, {i, 1, 20000, 2}];            (* 10000 existing edges *)
subV = Range[1, 50000];

cnt[h_] := {VertexCount[h], EdgeCount[h]};
tourOK[c_, gr_] := {Length[c], Sort[Sort /@ (List @@@ c)] === Sort[Sort /@ (List @@@ EdgeList[gr])]};

benchIf["VertexDelete, 1000 of 100000 vertices", "VertexDelete", VertexDelete[g, delV]];
checkIf["VertexDelete, 1000 of 100000 vertices", "VertexDelete", cnt[VertexDelete[g, delV]]];
benchIf["VertexAdd, 10000 vertices", "VertexAdd", VertexAdd[g, addV]];
checkIf["VertexAdd, 10000 vertices", "VertexAdd", cnt[VertexAdd[g, addV]]];
benchIf["EdgeAdd, 10000 edges", "EdgeAdd", EdgeAdd[g, addE]];
checkIf["EdgeAdd, 10000 edges", "EdgeAdd", cnt[EdgeAdd[g, addE]]];
benchIf["EdgeDelete, 10000 edges", "EdgeDelete", EdgeDelete[g, delE]];
checkIf["EdgeDelete, 10000 edges", "EdgeDelete", cnt[EdgeDelete[g, delE]]];
benchIf["Subgraph, 50000 of 100000 vertices", "Subgraph", Subgraph[g, subV]];
checkIf["Subgraph, 50000 of 100000 vertices", "Subgraph", cnt[Subgraph[g, subV]]];
benchIf["IndexGraph, 100000 vertices", "IndexGraph", IndexGraph[g, 0]];
checkIf["IndexGraph, 100000 vertices", "IndexGraph", cnt[IndexGraph[g, 0]]];
benchIf["GraphUnion, 2 x 100000 vertices", "GraphUnion", GraphUnion[g, g2]];
checkIf["GraphUnion, 2 x 100000 vertices", "GraphUnion", cnt[GraphUnion[g, g2]]];
benchIf["GraphIntersection, 2 x 2000 vertices", "GraphIntersection", GraphIntersection[h, h2]];
checkIf["GraphIntersection, 2 x 2000 vertices", "GraphIntersection", cnt[GraphIntersection[h, h2]]];
benchIf["GraphDifference, 2 x 2000 vertices", "GraphDifference", GraphDifference[h, h2]];
checkIf["GraphDifference, 2 x 2000 vertices", "GraphDifference", cnt[GraphDifference[h, h2]]];
benchIf["GraphDisjointUnion, 2 x 100000 vertices", "GraphDisjointUnion", GraphDisjointUnion[g, g2]];
checkIf["GraphDisjointUnion, 2 x 100000 vertices", "GraphDisjointUnion", cnt[GraphDisjointUnion[g, g2]]];
benchIf["GraphComplement, CycleGraph[1000]", "GraphComplement", GraphComplement[c1k]];
checkIf["GraphComplement, CycleGraph[1000]", "GraphComplement", cnt[GraphComplement[c1k]]];
benchIf["ReverseGraph, 199947 arcs", "ReverseGraph", ReverseGraph[dg]];
checkIf["ReverseGraph, 199947 arcs", "ReverseGraph", cnt[ReverseGraph[dg]]];
benchIf["UndirectedGraph, 199947 arcs", "UndirectedGraph", UndirectedGraph[dg]];
checkIf["UndirectedGraph, 199947 arcs", "UndirectedGraph", cnt[UndirectedGraph[dg]]];
benchIf["DirectedGraph, 199947 edges", "DirectedGraph", DirectedGraph[g]];
checkIf["DirectedGraph, 199947 edges", "DirectedGraph", cnt[DirectedGraph[g]]];
benchIf["LineGraph, 99997-edge graph", "LineGraph", LineGraph[g2]];
checkIf["LineGraph, 99997-edge graph", "LineGraph", cnt[LineGraph[g2]]];
benchIf["FindEulerianCycle, 100000 edges", "FindEulerianCycle", FindEulerianCycle[eu]];
checkIf["FindEulerianCycle, 100000 edges", "FindEulerianCycle", tourOK[First[FindEulerianCycle[eu]], eu]];
benchIf["FindPath, 100000 vertices", "FindPath", FindPath[g, 1, nv]];
checkIf["FindPath, 100000 vertices", "FindPath", {Length[FindPath[g, 1, nv]], Last[First[FindPath[g, 1, nv]]]}];
benchIf["FindCycle, 100000 vertices", "FindCycle", FindCycle[g]];
checkIf["FindCycle, 100000 vertices", "FindCycle", Length[FindCycle[g]]];

(* ---- COLD cases ---------------------------------------------------------- *)
newG[]  := Graph[Range[nv], gE];
newDG[] := Graph[Range[nv], dgE];
newEU[] := Graph[Range[ne], euE];
newC[]  := Graph[Range[1000], EdgeList[c1k]];
newPair[] := {Graph[Range[nv], gE], Graph[Range[nv], g2E]};
newSmallPair[] := {Graph[Range[ns], hE], Graph[Range[ns], h2E]};

coldBench[label_String, mk_, f_] := Module[{ts, gg},
  ts = Table[gg = mk[]; First[AbsoluteTiming[f[gg]]], {5}];
  Print["BENCH\t", label, "\t", ToString[Round[1000. Min[ts], 0.001]]]];

coldBench["VertexDelete, 1000 of 100000 vertices (cold)", newG, VertexDelete[#, delV] &];
check["VertexDelete, 1000 of 100000 vertices (cold)", cnt[VertexDelete[newG[], delV]]];
coldBench["EdgeAdd, 10000 edges (cold)", newG, EdgeAdd[#, addE] &];
check["EdgeAdd, 10000 edges (cold)", cnt[EdgeAdd[newG[], addE]]];
coldBench["EdgeDelete, 10000 edges (cold)", newG, EdgeDelete[#, delE] &];
check["EdgeDelete, 10000 edges (cold)", cnt[EdgeDelete[newG[], delE]]];
coldBench["Subgraph, 50000 of 100000 vertices (cold)", newG, Subgraph[#, subV] &];
check["Subgraph, 50000 of 100000 vertices (cold)", cnt[Subgraph[newG[], subV]]];
coldBench["GraphUnion, 2 x 100000 vertices (cold)", newPair, GraphUnion @@ # &];
check["GraphUnion, 2 x 100000 vertices (cold)", cnt[GraphUnion @@ newPair[]]];
coldBench["GraphDifference, 2 x 2000 vertices (cold)", newSmallPair, GraphDifference @@ # &];
check["GraphDifference, 2 x 2000 vertices (cold)", cnt[GraphDifference @@ newSmallPair[]]];
coldBench["GraphComplement, CycleGraph[1000] (cold)", newC, GraphComplement];
check["GraphComplement, CycleGraph[1000] (cold)", cnt[GraphComplement[newC[]]]];
coldBench["ReverseGraph, 199947 arcs (cold)", newDG, ReverseGraph];
check["ReverseGraph, 199947 arcs (cold)", cnt[ReverseGraph[newDG[]]]];
coldBench["FindEulerianCycle, 100000 edges (cold)", newEU, FindEulerianCycle];
check["FindEulerianCycle, 100000 edges (cold)", Length[First[FindEulerianCycle[newEU[]]]]];
coldBench["FindPath, 100000 vertices (cold)", newG, FindPath[#, 1, nv] &];
check["FindPath, 100000 vertices (cold)", Length[FindPath[newG[], 1, nv]]];
