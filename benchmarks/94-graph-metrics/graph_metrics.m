(* Experiment 94 -- Graph distances, centralities and clustering.
   MEASURES src/graph/gmet_*.c: GraphDistanceMatrix, BetweennessCentrality,
   ClosenessCentrality, PageRankCentrality, EigenvectorCentrality,
   GraphDiameter, GraphTriangleCount, LocalClusteringCoefficient,
   GlobalClusteringCoefficient, KatzCentrality.

   Every graph is deterministic (no RNG), so all three systems build the same
   object:
     sw[n]  -- ring i<->i+1 plus one long chord per vertex,
               i <-> i + 2 + Mod[7919 i^2 + 104729 i, n/2 - 2] (mod n);
               chord lengths < n/2 so there are no duplicate edges. Small world
               (diameter ~ log n), average degree 4.
     web    -- 100000 vertices, 5 out-arcs each, i -> i + k_t(i) (mod n) with
               the five offsets k_t in disjoint ranges (no duplicate arcs).
     band   -- 20000 vertices, i <-> i + d for d = 1..5: ~10^5 edges, rich in
               triangles.

   Checks are exact integers (rounded sums), so the systems cannot silently
   time different computations. networkx is the baseline (pure Python).

   COLD cases build a FRESH graph object before every timed call
   (construction untimed): Mathematica caches properties on its atomic Graph,
   and Mathilda caches results per graph node, so warm repeats compare caches. *)

Get["../harness.m"];

require[{"GraphDistanceMatrix", "BetweennessCentrality", "ClosenessCentrality",
         "PageRankCentrality", "EigenvectorCentrality", "GraphDiameter",
         "GraphTriangleCount", "LocalClusteringCoefficient",
         "GlobalClusteringCoefficient", "KatzCentrality"}];

swEdges[n_] := Join[Table[i <-> Mod[i, n] + 1, {i, n}],
  Table[i <-> Mod[i + 1 + Mod[7919 i^2 + 104729 i, Quotient[n, 2] - 2], n] + 1, {i, n}]];
nw = 100000;
webEdges = Flatten[Table[i -> Mod[i + t 19999 + Mod[7919 i^2 + 104729 i t + 31 t, 19998], nw] + 1,
                         {i, nw}, {t, 0, 4}]];
bandEdges = Flatten[Table[If[i + d <= 20000, i <-> i + d, Nothing], {i, 20000}, {d, 1, 5}]];

e10k = swEdges[10000];  e3k = swEdges[3000];
sw10k = Graph[Range[10000], e10k];
sw3k  = Graph[Range[3000], e3k];
web   = Graph[Range[nw], webEdges];
band  = Graph[Range[20000], bandEdges];

bench["BetweennessCentrality, 10^4-vertex small world", BetweennessCentrality[sw10k], 1];
check["BetweennessCentrality, 10^4-vertex small world", Round[Total[BetweennessCentrality[sw10k]]]];

bench["ClosenessCentrality, 10^4-vertex small world", ClosenessCentrality[sw10k]];
check["ClosenessCentrality, 10^4-vertex small world", Round[10^6 Total[ClosenessCentrality[sw10k]]]];

bench["GraphDiameter, 10^4-vertex small world", GraphDiameter[sw10k]];
check["GraphDiameter, 10^4-vertex small world", GraphDiameter[sw10k]];

bench["GraphDistanceMatrix, 3000-vertex small world", GraphDistanceMatrix[sw3k]];
check["GraphDistanceMatrix, 3000-vertex small world", Total[Flatten[GraphDistanceMatrix[sw3k]]]];

bench["PageRankCentrality, 10^5-vertex web graph", PageRankCentrality[web]];
check["PageRankCentrality, 10^5-vertex web graph", Round[10^8 Max[PageRankCentrality[web]]]];

bench["EigenvectorCentrality, 10^4-vertex small world", EigenvectorCentrality[sw10k]];
check["EigenvectorCentrality, 10^4-vertex small world", Round[10^8 Max[EigenvectorCentrality[sw10k]]]];

bench["KatzCentrality, 10^5-vertex web graph", KatzCentrality[web, 0.1]];
check["KatzCentrality, 10^5-vertex web graph", Round[10^9 Max[KatzCentrality[web, 0.1]]]];

bench["GraphTriangleCount, 10^5-edge band graph", GraphTriangleCount[band]];
check["GraphTriangleCount, 10^5-edge band graph", GraphTriangleCount[band]];

bench["GlobalClusteringCoefficient, 10^5-edge band graph", GlobalClusteringCoefficient[band]];
check["GlobalClusteringCoefficient, 10^5-edge band graph", GlobalClusteringCoefficient[band]];

bench["LocalClusteringCoefficient, 10^5-edge band graph", LocalClusteringCoefficient[band]];
check["LocalClusteringCoefficient, 10^5-edge band graph", Total[LocalClusteringCoefficient[band]]];

(* ---- COLD cases ---------------------------------------------------------- *)

newSw10k[] := Graph[Range[10000], e10k];
newSw3k[]  := Graph[Range[3000], e3k];
newWeb[]   := Graph[Range[nw], webEdges];
newBand[]  := Graph[Range[20000], bandEdges];

coldBench[label_String, mk_, f_] := Module[{ts, gg},
  ts = Table[gg = mk[]; First[AbsoluteTiming[f[gg]]], {3}];
  Print["BENCH\t", label, "\t", ToString[Round[1000. Min[ts], 0.001]]]];

coldBench["BetweennessCentrality, 10^4-vertex small world (cold)", newSw10k, BetweennessCentrality];
check["BetweennessCentrality, 10^4-vertex small world (cold)", Round[Total[BetweennessCentrality[newSw10k[]]]]];
coldBench["ClosenessCentrality, 10^4-vertex small world (cold)", newSw10k, ClosenessCentrality];
check["ClosenessCentrality, 10^4-vertex small world (cold)", Round[10^6 Total[ClosenessCentrality[newSw10k[]]]]];
coldBench["GraphDiameter, 10^4-vertex small world (cold)", newSw10k, GraphDiameter];
check["GraphDiameter, 10^4-vertex small world (cold)", GraphDiameter[newSw10k[]]];
coldBench["GraphDistanceMatrix, 3000-vertex small world (cold)", newSw3k, GraphDistanceMatrix];
check["GraphDistanceMatrix, 3000-vertex small world (cold)", Total[Flatten[GraphDistanceMatrix[newSw3k[]]]]];
coldBench["PageRankCentrality, 10^5-vertex web graph (cold)", newWeb, PageRankCentrality];
check["PageRankCentrality, 10^5-vertex web graph (cold)", Round[10^8 Max[PageRankCentrality[newWeb[]]]]];
coldBench["EigenvectorCentrality, 10^4-vertex small world (cold)", newSw10k, EigenvectorCentrality];
check["EigenvectorCentrality, 10^4-vertex small world (cold)", Round[10^8 Max[EigenvectorCentrality[newSw10k[]]]]];
coldBench["KatzCentrality, 10^5-vertex web graph (cold)", newWeb, KatzCentrality[#, 0.1] &];
check["KatzCentrality, 10^5-vertex web graph (cold)", Round[10^9 Max[KatzCentrality[newWeb[], 0.1]]]];
coldBench["GraphTriangleCount, 10^5-edge band graph (cold)", newBand, GraphTriangleCount];
check["GraphTriangleCount, 10^5-edge band graph (cold)", GraphTriangleCount[newBand[]]];
coldBench["GlobalClusteringCoefficient, 10^5-edge band graph (cold)", newBand, GlobalClusteringCoefficient];
check["GlobalClusteringCoefficient, 10^5-edge band graph (cold)", GlobalClusteringCoefficient[newBand[]]];
coldBench["LocalClusteringCoefficient, 10^5-edge band graph (cold)", newBand, LocalClusteringCoefficient];
check["LocalClusteringCoefficient, 10^5-edge band graph (cold)", Total[LocalClusteringCoefficient[newBand[]]]];
