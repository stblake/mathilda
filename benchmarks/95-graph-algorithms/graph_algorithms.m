(* Experiment 95 -- graph algorithms: flows/cuts, matching, cliques, covers,
   Hamiltonian cycles, isomorphism, planarity.
   MEASURES the algos stream (src/graph/galg_*.c).

   Every graph is deterministic: pseudo-random graphs come from the same
   linear congruential generator in all three systems (lcg below and in
   graph_algorithms.py), so the systems build identical objects:
     rnd[n, m, s]   first m distinct pairs of an LCG stream (simple, undirected)
     grid[w, h]     w x h grid; trigrid adds one diagonal per cell (planar,
                    non-bipartite triangulation)
     cubic[n, s]    random Hamiltonian cycle plus a random perfect matching
                    (duplicates dropped); the seed 28 at n = 10^4 has no
                    duplicate, so the graph is exactly 3-regular -- colour
                    refinement alone learns nothing, the hard case for
                    isomorphism -- and cubicP is its random relabelling
     m3[n, s]       union of three random perfect matchings (no planted
                    Hamiltonian cycle) -- the Hamiltonicity instance

   Warm cases repeat one call on one graph (both CASes may answer from a
   per-graph cache); COLD cases build a fresh graph before every timed call
   (construction untimed), so no cache helps.  Checks pin unique values
   (flow values, optimum sizes, truth values) so the columns time the same
   computation. *)

Get["../harness.m"];

require[{"FindMaximumFlow", "FindMinimumCut", "EdgeConnectivity", "FindIndependentEdgeSet",
         "FindEdgeCover", "FindClique", "FindIndependentVertexSet", "FindVertexCover",
         "FindHamiltonianCycle", "IsomorphicGraphQ", "FindGraphIsomorphism", "PlanarGraphQ"}];

lcg[len_, s_] := Rest[NestList[Mod[1103515245 # + 12345, 2147483648] &, s, len]];
rndPairs[n_, m_, s_] := Module[{r, a, b, p},
  r = lcg[2 m + m + 1000, s];
  a = Mod[Quotient[r[[1 ;; ;; 2]], 65536], n] + 1;
  b = Mod[Quotient[r[[2 ;; ;; 2]], 65536], n] + 1;
  p = DeleteDuplicates[Select[Transpose[{Min /@ Transpose[{a, b}], Max /@ Transpose[{a, b}]}],
                              #[[1]] < #[[2]] &]];
  Take[p, Min[m, Length[p]]]];
perm[n_, s_] := Ordering[lcg[n, s]];
gridPairs[w_, h_] := Join[
  Flatten[Table[{(i - 1) w + j, (i - 1) w + j + 1}, {i, h}, {j, w - 1}], 1],
  Flatten[Table[{(i - 1) w + j, i w + j}, {i, h - 1}, {j, w}], 1]];
triPairs[w_, h_] := Join[gridPairs[w, h],
  Flatten[Table[{(i - 1) w + j, i w + j + 1}, {i, h - 1}, {j, w - 1}], 1]];
cubicPairs[n_, s_] := Module[{p = perm[n, s], q = perm[n, s + 1], c, mt},
  c = Append[Transpose[{Most[p], Rest[p]}], {Last[p], First[p]}];
  mt = Partition[q, 2];
  DeleteDuplicates[Sort /@ Join[c, mt]]];
m3Pairs[n_, s_] := DeleteDuplicates[Sort /@ Join[Partition[perm[n, s], 2],
  Partition[perm[n, s + 1], 2], Partition[perm[n, s + 2], 2]]];
mk[n_, pairs_] := Graph[Range[n], UndirectedEdge @@@ pairs];

(* ---- the graphs ------------------------------------------------------- *)
nR = 20000; rP = rndPairs[nR, 100000, 7];            rnd1 = mk[nR, rP];
nC = 2000;  cP = rndPairs[nC, 20000, 11];             rnd2 = mk[nC, cP];
dP = rndPairs[200, 10000, 13];                        dense = mk[200, dP];
sP = rndPairs[100, 300, 17];                          sparse = mk[100, sP];
gw = 250; gh = 200; nG = gw gh; gP = gridPairs[gw, gh]; grid = mk[nG, gP];
fw = 500; fh = 400; nF = fw fh; fP = gridPairs[fw, fh]; gridF = mk[nF, fP];
capF = Table[Mod[7919 k, 100] + 1, {k, Length[fP]}];
nM = 200000; mP = rndPairs[nM, 1000000, 7];          rndM = mk[nM, mP];
capM = Table[Mod[104729 k, 1000] + 1, {k, Length[mP]}];
(* the capacity options are built once: the timed call is the flow, not the
   evaluation of a 10^6-element Rule *)
optF = EdgeCapacity -> capF; optM = EdgeCapacity -> capM;
tw = 316; nT = tw tw; tP = triPairs[tw, tw];          tri = mk[nT, tP];
nI = 10000; iP = cubicPairs[nI, 28]; rl = perm[nI, 29];
cub = mk[nI, iP];
cPP = SortBy[Map[rl[[#]] &, iP, {2}], First[lcg[1, #[[1]] + 7 #[[2]]]] &];
cubP = mk[nI, cPP];
hP = m3Pairs[400, 41];                                ham = mk[400, hP];

(* ---- warm cases ----------------------------------------------------------- *)
bench["FindMaximumFlow, 2*10^5-vertex grid, capacities", FindMaximumFlow[gridF, 1, nF, optF]];
check["FindMaximumFlow, 2*10^5-vertex grid, capacities", FindMaximumFlow[gridF, 1, nF, optF]];
bench["FindMaximumFlow, 10^6-edge random graph, capacities", FindMaximumFlow[rndM, 1, 2, optM]];
check["FindMaximumFlow, 10^6-edge random graph, capacities", FindMaximumFlow[rndM, 1, 2, optM]];
bench["FindMinimumCut, 2000 vertices 20000 edges", FindMinimumCut[rnd2]];
check["FindMinimumCut, 2000 vertices 20000 edges", First[FindMinimumCut[rnd2]]];
bench["EdgeConnectivity, 50000-vertex grid", EdgeConnectivity[grid]];
check["EdgeConnectivity, 50000-vertex grid", EdgeConnectivity[grid]];
bench["FindIndependentEdgeSet, 10^5-edge random graph", FindIndependentEdgeSet[rnd1]];
check["FindIndependentEdgeSet, 10^5-edge random graph", Length[FindIndependentEdgeSet[rnd1]]];
bench["FindIndependentEdgeSet, 10^5-vertex triangulated grid", FindIndependentEdgeSet[tri]];
check["FindIndependentEdgeSet, 10^5-vertex triangulated grid", Length[FindIndependentEdgeSet[tri]]];
bench["FindEdgeCover, 50000-vertex grid", FindEdgeCover[grid]];
check["FindEdgeCover, 50000-vertex grid", Length[FindEdgeCover[grid]]];
bench["FindClique, dense 200-vertex random graph", FindClique[dense]];
check["FindClique, dense 200-vertex random graph", Length[First[FindClique[dense]]]];
bench["FindIndependentVertexSet, 100 vertices 300 edges", FindIndependentVertexSet[sparse]];
check["FindIndependentVertexSet, 100 vertices 300 edges", Length[First[FindIndependentVertexSet[sparse]]]];
bench["FindVertexCover, 100 vertices 300 edges", FindVertexCover[sparse]];
check["FindVertexCover, 100 vertices 300 edges", Length[FindVertexCover[sparse]]];
bench["HamiltonianGraphQ, 400-vertex random cubic graph", HamiltonianGraphQ[ham]];
check["HamiltonianGraphQ, 400-vertex random cubic graph", HamiltonianGraphQ[ham]];
bench["IsomorphicGraphQ, 10^4-vertex 3-regular random graph, permuted", IsomorphicGraphQ[cub, cubP]];
check["IsomorphicGraphQ, 10^4-vertex 3-regular random graph, permuted", IsomorphicGraphQ[cub, cubP]];
bench["FindGraphIsomorphism, 10^4-vertex 3-regular random graph", FindGraphIsomorphism[cub, cubP]];
check["FindGraphIsomorphism, 10^4-vertex 3-regular random graph", Length[FindGraphIsomorphism[cub, cubP]]];
bench["PlanarGraphQ, 10^5-vertex triangulated grid", PlanarGraphQ[tri]];
check["PlanarGraphQ, 10^5-vertex triangulated grid", PlanarGraphQ[tri]];

(* ---- COLD cases: a fresh graph before every timed call ------------------- *)
coldBench[label_String, mk0_, f_] := Module[{ts, gg},
  ts = Table[gg = mk0[]; First[AbsoluteTiming[f[gg]]], {5}];
  Print["BENCH\t", label, "\t", ToString[Round[1000. Min[ts], 0.001]]]];

newGrid[] := mk[nG, gP]; newRnd1[] := mk[nR, rP]; newRnd2[] := mk[nC, cP];
newTri[] := mk[nT, tP]; newDense[] := mk[200, dP]; newSparse[] := mk[100, sP];
newCub[] := mk[nI, iP]; newCubPair[] := {mk[nI, iP], mk[nI, cPP]}; newGridF[] := mk[nF, fP]; newRndM[] := mk[nM, mP]; newHam[] := mk[400, hP];

coldBench["FindMaximumFlow, 2*10^5-vertex grid, capacities (cold)", newGridF,
  FindMaximumFlow[#, 1, nF, optF] &];
check["FindMaximumFlow, 2*10^5-vertex grid, capacities (cold)", FindMaximumFlow[newGridF[], 1, nF, optF]];
coldBench["FindMaximumFlow, 10^6-edge random graph, capacities (cold)", newRndM,
  FindMaximumFlow[#, 1, 2, optM] &];
check["FindMaximumFlow, 10^6-edge random graph, capacities (cold)", FindMaximumFlow[newRndM[], 1, 2, optM]];
coldBench["FindMinimumCut, 2000 vertices 20000 edges (cold)", newRnd2, FindMinimumCut];
check["FindMinimumCut, 2000 vertices 20000 edges (cold)", First[FindMinimumCut[newRnd2[]]]];
coldBench["FindIndependentEdgeSet, 10^5-edge random graph (cold)", newRnd1, FindIndependentEdgeSet];
check["FindIndependentEdgeSet, 10^5-edge random graph (cold)", Length[FindIndependentEdgeSet[newRnd1[]]]];
coldBench["FindIndependentEdgeSet, 10^5-vertex triangulated grid (cold)", newTri, FindIndependentEdgeSet];
check["FindIndependentEdgeSet, 10^5-vertex triangulated grid (cold)", Length[FindIndependentEdgeSet[newTri[]]]];
coldBench["FindClique, dense 200-vertex random graph (cold)", newDense, FindClique];
check["FindClique, dense 200-vertex random graph (cold)", Length[First[FindClique[newDense[]]]]];
coldBench["FindVertexCover, 100 vertices 300 edges (cold)", newSparse, FindVertexCover];
check["FindVertexCover, 100 vertices 300 edges (cold)", Length[FindVertexCover[newSparse[]]]];
coldBench["HamiltonianGraphQ, 400-vertex random cubic graph (cold)", newHam, HamiltonianGraphQ];
check["HamiltonianGraphQ, 400-vertex random cubic graph (cold)", HamiltonianGraphQ[newHam[]]];
(* both graphs fresh: neither system may reuse anything computed for either *)
coldBench["IsomorphicGraphQ, 10^4-vertex 3-regular random graph, permuted (cold)", newCubPair, IsomorphicGraphQ @@ # &];
check["IsomorphicGraphQ, 10^4-vertex 3-regular random graph, permuted (cold)", IsomorphicGraphQ @@ newCubPair[]];
coldBench["PlanarGraphQ, 10^5-vertex triangulated grid (cold)", newTri, PlanarGraphQ];
check["PlanarGraphQ, 10^5-vertex triangulated grid (cold)", PlanarGraphQ[newTri[]]];
