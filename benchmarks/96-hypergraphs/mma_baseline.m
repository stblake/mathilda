(* mma_baseline.m -- the Mathematica column of experiment 96, loaded ONLY by
   Mathematica (hypergraphs.m does `If[$VersionNumber > 10, Get[...]]`).

   Mathematica 15 has no built-in hypergraph symbols (Names["*Hypergraph*"] is
   empty).  The Wolfram-side baseline is therefore, in order of preference:

     1. the official WolframInstitute/Hypergraph paclet (PacletInstall[
        "WolframInstitute/Hypergraph"]) -- its Hypergraph object, VertexDegree,
        VertexList, EdgeList;
     2. Wolfram Function Repository ResourceFunctions with the same semantics
        -- HypergraphToGraph, ConnectedHypergraphQ, RandomHypergraph,
        TransversalHypergraph (fetched over the network on first use);
     3. where neither exists, the idiomatic plain-WL implementation a Wolfram
        user would write, using the fastest built-ins available (SparseArray
        products for line graphs, GroupBy for duals, the kernel's Graph
        algorithms, LinearOptimization for the minimum hitting set).

   Each definition below is marked with its source.  The paclet is loaded
   without putting its context on $ContextPath, so its HypergraphToGraph /
   ConnectedHypergraphQ / RandomHypergraph (different semantics from the FR
   functions of the same names) cannot shadow the definitions here. *)

Block[{$ContextPath = $ContextPath}, Needs["WolframInstitute`Hypergraph`"]];
HG = WolframInstitute`Hypergraph`Hypergraph;

(* [paclet] the object; VertexDegree/VertexList/EdgeList are paclet overloads *)
Hypergraph[v_List, e_List] := HG[v, e];

(* [plain WL] *)
HyperedgeSizes[h_] := Length /@ EdgeList[h];

(* [plain WL] dual as its edge lists (building a paclet object for the result
   would only add paclet-construction time to the Mathematica column) *)
HypergraphDual[h_] := Module[{e = EdgeList[h], inc},
  inc = GroupBy[Catenate[MapIndexed[Function[{s, j}, {#, First[j]} & /@ DeleteDuplicates[s]], e]],
                First -> Last];
  Lookup[inc, VertexList[h], {}]];

(* [plain WL] 2-section *)
HypergraphCliqueExpansion[h_] :=
  Graph[VertexList[h], UndirectedEdge @@@
    DeleteDuplicates[Sort /@ Catenate[Subsets[DeleteDuplicates[#], {2}] & /@ EdgeList[h]]]];

(* [plain WL] incidence bipartite graph *)
HypergraphStarExpansion[h_] := Module[{e = EdgeList[h]},
  Graph[Join[VertexList[h], Hyperedge /@ Range[Length[e]]],
        Catenate[MapIndexed[Function[{s, j}, UndirectedEdge[#, Hyperedge @@ j] & /@ DeleteDuplicates[s]], e]]]];

(* [FR] *)
HypergraphToGraph[e_List] := ResourceFunction["HypergraphToGraph"][e];
ConnectedHypergraphQ[e_List] := ResourceFunction["ConnectedHypergraphQ"][e];
RandomHypergraph[spec_List] := ResourceFunction["RandomHypergraph"][spec];
TransversalHypergraph[e_List] := ResourceFunction["TransversalHypergraph"][e];

(* [plain WL] sparse incidence matrix (vertices -> rows) *)
incSparse[h_] := Module[{v = VertexList[h], e = EdgeList[h], idx},
  idx = AssociationThread[v, Range[Length[v]]];
  SparseArray[Catenate[MapIndexed[Function[{s, j}, {idx[#], First[j]} -> 1 & /@ DeleteDuplicates[s]], e]],
              {Length[v], Length[e]}]];

(* [plain WL] s-line graph via the overlap matrix I^T.I *)
HypergraphLineGraph[h_, s_: 1] := Module[{im = incSparse[h], a, m},
  m = Length[EdgeList[h]];
  a = UpperTriangularize[Transpose[im].im, 1];
  Graph[Range[m], UndirectedEdge @@@ Keys[Select[Most[ArrayRules[a]], Last[#] >= s &]]]];

(* [plain WL] components through a path per hyperedge *)
HypergraphConnectedComponents[h_] :=
  ConnectedComponents[Graph[VertexList[h],
    UndirectedEdge @@@ Catenate[Partition[DeleteDuplicates[#], 2, 1] & /@ EdgeList[h]]]];

(* [plain WL] *)
HyperedgeConnectedComponents[h_, s_: 1] := ConnectedComponents[HypergraphLineGraph[h, s]];
(* Single-source distances by BreadthFirstScan: the kernel's GraphDistance[g, s]
   is quadratic here (8.5 s at 10^4 vertices, projecting to ~15 min at 10^5),
   so the Mathematica column uses the faster idiom (0.09 s at 10^4). *)
bfsDistances[g_, s_] := Module[{d = ConstantArray[Infinity, VertexCount[g]]},
  d[[VertexIndex[g, s]]] = 0;
  BreadthFirstScan[g, s, {"DiscoverVertex" -> ((d[[VertexIndex[g, #1]]] = #3) &)}];
  d];
HyperedgeDistance[h_, i_Integer] := bfsDistances[HypergraphLineGraph[h], i];
HypergraphDistance[h_, u_] := bfsDistances[HypergraphCliqueExpansion[h], u];

(* [plain WL] dense incidence matrix, the same value Mathilda returns *)
Unprotect[IncidenceMatrix];
IncidenceMatrix[h_HG] := Normal[incSparse[h]];
Protect[IncidenceMatrix];

(* [plain WL] minimum hitting set as a 0-1 integer program (the kernel's
   MILP solver).  NOTE when running by hand: redirect with a pipe or `>>`, not
   `>` -- the solver writes to the same stdout descriptor at its own offset and
   overwrites earlier lines of a truncating redirect. *)
FindMinimumTransversal[e_List] := Module[{v = Union @@ e, x, vars, sol},
  vars = x /@ v;
  sol = LinearOptimization[Total[vars],
          Join[(Total[x /@ #] >= 1) & /@ e, Thread[0 <= vars <= 1]],
          vars \[Element] Integers];
  Select[v, (x[#] /. sol) == 1 &]];
