(* Experiment 29 -- Graph operations.
   MEASURES src/graph/ (19 files): construction, degree sequences, connected
   components, shortest paths, spanning trees.

   NOT COVERED BY THE EXISTING CORPUS.  Experiment 12 measured graph ALGORITHMS
   written by hand out of gathers and segmented sums, deliberately avoiding a
   graph type, to isolate the gather.  This experiment measures the opposite: the
   src/graph/ builtins themselves, which nothing has ever timed.

   Baseline is networkx, which is pure Python -- so a Mathilda win here is a
   weaker claim than a win against a compiled library, and the .py docstring says
   so.  What the row set is really for is finding builtins that are absent or
   accidentally quadratic. *)

Get["../harness.m"];
Get["../data.m"];

require[{"Graph", "EdgeList", "VertexList", "VertexDegree", "AdjacencyMatrix",
         "ConnectedComponents", "FindShortestPath", "GraphDistance",
         "ConnectedGraphQ", "FindSpanningTree", "VertexCount", "EdgeCount"}];

(* A deterministic sparse graph: a cycle plus chords, 20000 vertices. *)
nv = 20000;
edges = Join[Table[i -> Mod[i, nv] + 1, {i, 1, nv}],
             Table[i -> Mod[3 i, nv] + 1, {i, 1, nv}]];
g = Graph[edges];

bench["Graph construction, 40000 edges", Graph[edges], 1];
check["Graph construction, 40000 edges", Length[edges]];

bench["VertexDegree, 20000 vertices", VertexDegree[g];];
check["VertexDegree, 20000 vertices",
  Total[VertexDegree[Graph[{1 -> 2, 2 -> 3, 3 -> 1}]]]];

bench["EdgeCount", EdgeCount[g];];
check["EdgeCount", EdgeCount[Graph[{1 -> 2, 2 -> 3, 3 -> 1}]]];

benchIf["ConnectedComponents", "ConnectedComponents",
  ConnectedComponents[g], 1];
(* Mathematica counts differently on a DIRECTED graph, so the check uses an
   undirected two-component graph where all three systems must agree. *)
checkIf["ConnectedComponents", "ConnectedComponents",
  Length[ConnectedComponents[Graph[{1 <-> 2, 3 <-> 4}]]]];

benchIf["GraphDistance from vertex 1", "GraphDistance",
  GraphDistance[g, 1], 1];
checkIf["GraphDistance from vertex 1", "GraphDistance",
  GraphDistance[Graph[{1 -> 2, 2 -> 3}], 1, 3]];

benchIf["FindShortestPath 1 to 10000", "FindShortestPath",
  FindShortestPath[g, 1, 10000], 1];
checkIf["FindShortestPath 1 to 10000", "FindShortestPath",
  Length[FindShortestPath[Graph[{1 -> 2, 2 -> 3}], 1, 3]]];

benchIf["AdjacencyMatrix, 20000 vertices", "AdjacencyMatrix",
  AdjacencyMatrix[g], 1];
checkIf["AdjacencyMatrix, 20000 vertices", "AdjacencyMatrix",
  Total[Flatten[AdjacencyMatrix[Graph[{1 -> 2, 2 -> 3, 3 -> 1}]]]]];
