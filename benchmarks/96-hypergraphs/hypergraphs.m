(* Experiment 96 -- the native hypergraph subsystem (src/graph/hyp_*.c).

   Mathematica 15 has no built-in hypergraphs, so its column is the Wolfram
   ecosystem's own hypergraph code -- the WolframInstitute/Hypergraph paclet,
   the Function Repository ResourceFunctions (HypergraphToGraph,
   ConnectedHypergraphQ, RandomHypergraph, TransversalHypergraph), and plain WL
   where neither has the function.  All of that is defined in mma_baseline.m,
   which only Mathematica loads; each definition there is tagged with its
   source.  The Python column is xgi 0.10 (pure Python over networkx/numpy/
   scipy); the minimum hitting set uses scipy's HiGHS MILP.

   Hypergraphs are deterministic so every system builds the same object:
     hB -- 20000 vertices, 100000 3-uniform hyperedges {a, b, c}:
           a = 7919 i mod n, b = a + 1 + (i mod 101), c = b + 1 + (3i mod 103)
     hL -- the same rule with n = 100000 (sparse: line graph / distances)
     hI -- the same rule, 1000 vertices x 10000 hyperedges (dense incidence)
     eT -- 16 vertices, 32 hyperedges (offsets mod 5 and 7): transversals
     eM -- 80 vertices, 200 LCG-random triples: minimum hitting set

   COLD cases build a fresh hypergraph object before each timed call
   (construction untimed): Mathilda's validated-hypergraph memo, like
   Mathematica's atomic Graph caches, would otherwise answer warm repeats from
   per-object state.  The paclet and xgi keep no such state, so for them warm
   and cold differ only by noise. *)

Get["../harness.m"];
If[$VersionNumber > 10, Get["mma_baseline.m"]];

require[{"Hypergraph", "HyperedgeSizes", "HypergraphDual", "HypergraphCliqueExpansion",
         "HypergraphStarExpansion", "HypergraphToGraph", "HypergraphLineGraph",
         "HypergraphConnectedComponents", "ConnectedHypergraphQ",
         "HyperedgeConnectedComponents", "HyperedgeDistance", "HypergraphDistance",
         "RandomHypergraph", "TransversalHypergraph", "FindMinimumTransversal"}];

mkE[n_, m_, p_, q_] := Table[With[{a = Mod[7919 i, n]},
    With[{b = Mod[a + 1 + Mod[i, p], n]}, {a, b, Mod[b + 1 + Mod[3 i, q], n]} + 1]], {i, m}];
lcg[n_, m_] := Partition[Mod[Quotient[Rest[NestList[Mod[1103515245 # + 12345, 2147483648] &,
    42, 3 m]], 65536], n] + 1, 3];

nB = 20000;  eB = mkE[nB, 100000, 101, 103];
nL = 100000; eL = mkE[nL, 100000, 101, 103];
nI = 1000;   eI = mkE[nI, 10000, 101, 103];
eT = mkE[16, 32, 5, 7];
eM = lcg[80, 200];

hB = Hypergraph[Range[nB], eB];
hL = Hypergraph[Range[nL], eL];
hI = Hypergraph[Range[nI], eI];

edgesOf[x_List] := x;
edgesOf[x_] := EdgeList[x];
finite[d_] := {Max[Cases[d, _Integer]], Count[d, _Integer]};

bench["Hypergraph construction, 10^5 hyperedges", Hypergraph[Range[nB], eB]];
check["Hypergraph construction, 10^5 hyperedges", Length[VertexList[Hypergraph[Range[nB], eB]]]];

bench["VertexDegree, 10^5 hyperedges", VertexDegree[hB]];
check["VertexDegree, 10^5 hyperedges", Total[VertexDegree[hB]]];

bench["HyperedgeSizes, 10^5 hyperedges", HyperedgeSizes[hB]];
check["HyperedgeSizes, 10^5 hyperedges", Total[HyperedgeSizes[hB]]];

bench["HypergraphDual, 10^5 hyperedges", HypergraphDual[hB]];
check["HypergraphDual, 10^5 hyperedges", Total[Length /@ edgesOf[HypergraphDual[hB]]]];

bench["HypergraphCliqueExpansion, 10^5 hyperedges", HypergraphCliqueExpansion[hB]];
check["HypergraphCliqueExpansion, 10^5 hyperedges", EdgeCount[HypergraphCliqueExpansion[hB]]];

bench["HypergraphStarExpansion, 10^5 hyperedges", HypergraphStarExpansion[hB]];
check["HypergraphStarExpansion, 10^5 hyperedges", EdgeCount[HypergraphStarExpansion[hB]]];

bench["HypergraphToGraph (FR), 10^5 hyperedges", HypergraphToGraph[eB]];
check["HypergraphToGraph (FR), 10^5 hyperedges", Length[Union[EdgeList[HypergraphToGraph[eB]]]]];

bench["HypergraphLineGraph, 10^5 hyperedges", HypergraphLineGraph[hL]];
check["HypergraphLineGraph, 10^5 hyperedges", EdgeCount[HypergraphLineGraph[hL]]];

bench["HypergraphConnectedComponents, 10^5 hyperedges", HypergraphConnectedComponents[hL]];
check["HypergraphConnectedComponents, 10^5 hyperedges",
  Sort[Length /@ HypergraphConnectedComponents[hL]]];

bench["ConnectedHypergraphQ (FR), 10^5 hyperedges", ConnectedHypergraphQ[eB]];
check["ConnectedHypergraphQ (FR), 10^5 hyperedges", ConnectedHypergraphQ[eB]];

bench["HyperedgeConnectedComponents s=2, 10^5 hyperedges", HyperedgeConnectedComponents[hB, 2]];
check["HyperedgeConnectedComponents s=2, 10^5 hyperedges", Length[HyperedgeConnectedComponents[hB, 2]]];

bench["HyperedgeDistance from e1, 10^5 hyperedges", HyperedgeDistance[hL, 1]];
check["HyperedgeDistance from e1, 10^5 hyperedges", finite[HyperedgeDistance[hL, 1]]];

bench["HypergraphDistance from v1, 10^5 hyperedges", HypergraphDistance[hL, 1]];
check["HypergraphDistance from v1, 10^5 hyperedges", finite[HypergraphDistance[hL, 1]]];

bench["IncidenceMatrix, 1000 x 10^4", IncidenceMatrix[hI]];
check["IncidenceMatrix, 1000 x 10^4", Total[IncidenceMatrix[hI], 2]];

bench["RandomHypergraph (FR form), 10^5 triples", RandomHypergraph[{nB, {100000, 3}}]];
check["RandomHypergraph (FR form), 10^5 triples", Length[edgesOf[RandomHypergraph[{nB, {100000, 3}}]]]];

benchOnce["TransversalHypergraph (FR), 16 vertices", TransversalHypergraph[eT]];
check["TransversalHypergraph (FR), 16 vertices", Sort[Sort /@ TransversalHypergraph[eT]] // Length];

bench["FindMinimumTransversal, 80 vertices x 200", FindMinimumTransversal[eM]];
check["FindMinimumTransversal, 80 vertices x 200", Length[FindMinimumTransversal[eM]]];

(* ---- COLD cases: a fresh object per timed call, construction untimed ---- *)

coldBench[label_String, mk_, f_] := Module[{ts, hh},
  ts = Table[hh = mk[]; First[AbsoluteTiming[f[hh]]], {3}];
  Print["BENCH\t", label, "\t", ToString[Round[1000. Min[ts], 0.001]]]];

newB[] := Hypergraph[Range[nB], eB];
newL[] := Hypergraph[Range[nL], eL];
newI[] := Hypergraph[Range[nI], eI];

coldBench["VertexDegree, 10^5 hyperedges (cold)", newB, VertexDegree];
check["VertexDegree, 10^5 hyperedges (cold)", Total[VertexDegree[newB[]]]];
coldBench["HyperedgeSizes, 10^5 hyperedges (cold)", newB, HyperedgeSizes];
check["HyperedgeSizes, 10^5 hyperedges (cold)", Total[HyperedgeSizes[newB[]]]];
coldBench["HypergraphDual, 10^5 hyperedges (cold)", newB, HypergraphDual];
check["HypergraphDual, 10^5 hyperedges (cold)", Total[Length /@ edgesOf[HypergraphDual[newB[]]]]];
coldBench["HypergraphCliqueExpansion, 10^5 hyperedges (cold)", newB, HypergraphCliqueExpansion];
check["HypergraphCliqueExpansion, 10^5 hyperedges (cold)", EdgeCount[HypergraphCliqueExpansion[newB[]]]];
coldBench["HypergraphStarExpansion, 10^5 hyperedges (cold)", newB, HypergraphStarExpansion];
check["HypergraphStarExpansion, 10^5 hyperedges (cold)", EdgeCount[HypergraphStarExpansion[newB[]]]];
coldBench["HypergraphLineGraph, 10^5 hyperedges (cold)", newL, HypergraphLineGraph];
check["HypergraphLineGraph, 10^5 hyperedges (cold)", EdgeCount[HypergraphLineGraph[newL[]]]];
coldBench["HypergraphConnectedComponents, 10^5 hyperedges (cold)", newL, HypergraphConnectedComponents];
check["HypergraphConnectedComponents, 10^5 hyperedges (cold)",
  Sort[Length /@ HypergraphConnectedComponents[newL[]]]];
coldBench["HyperedgeConnectedComponents s=2, 10^5 hyperedges (cold)", newB,
  HyperedgeConnectedComponents[#, 2] &];
check["HyperedgeConnectedComponents s=2, 10^5 hyperedges (cold)",
  Length[HyperedgeConnectedComponents[newB[], 2]]];
coldBench["HyperedgeDistance from e1, 10^5 hyperedges (cold)", newL, HyperedgeDistance[#, 1] &];
check["HyperedgeDistance from e1, 10^5 hyperedges (cold)", finite[HyperedgeDistance[newL[], 1]]];
coldBench["HypergraphDistance from v1, 10^5 hyperedges (cold)", newL, HypergraphDistance[#, 1] &];
check["HypergraphDistance from v1, 10^5 hyperedges (cold)", finite[HypergraphDistance[newL[], 1]]];
coldBench["IncidenceMatrix, 1000 x 10^4 (cold)", newI, IncidenceMatrix];
check["IncidenceMatrix, 1000 x 10^4 (cold)", Total[IncidenceMatrix[newI[]], 2]];
