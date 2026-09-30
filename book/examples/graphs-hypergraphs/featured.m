# Featured examples (text half; the pictures are in figures/graphs-hypergraphs).
moves = {{1, 2}, {2, 1}, {-1, 2}, {-2, 1}};
knight[n_] := Graph[Flatten[Table[If[1 <= i + m[[1]] <= n && 1 <= j + m[[2]] <= n, {i, j} <-> {i + m[[1]], j + m[[2]]}, Nothing], {i, n}, {j, n}, {m, moves}]]]
k8 = knight[8]
SortBy[Tally[VertexDegree[k8]], First]
Length[First[FindHamiltonianCycle[k8]]]
{FindHamiltonianCycle[knight[5]], Length[FindHamiltonianPath[knight[5]]]}
florence = Graph[{"Acciaiuoli" <-> "Medici", "Castellani" <-> "Peruzzi", "Castellani" <-> "Strozzi", "Castellani" <-> "Barbadori", "Medici" <-> "Barbadori", "Medici" <-> "Ridolfi", "Medici" <-> "Tornabuoni", "Medici" <-> "Albizzi", "Medici" <-> "Salviati", "Salviati" <-> "Pazzi", "Peruzzi" <-> "Strozzi", "Peruzzi" <-> "Bischeri", "Strozzi" <-> "Ridolfi", "Strozzi" <-> "Bischeri", "Ridolfi" <-> "Tornabuoni", "Tornabuoni" <-> "Guadagni", "Albizzi" <-> "Ginori", "Albizzi" <-> "Guadagni", "Bischeri" <-> "Guadagni", "Guadagni" <-> "Lamberteschi"}];
TakeLargestBy[Thread[VertexList[florence] -> BetweennessCentrality[florence]], Last, 3]
TakeLargestBy[Thread[VertexList[florence] -> VertexDegree[florence]], Last, 3]
FindVertexCut[florence]
FindShortestPath[florence, "Pazzi", "Lamberteschi"]
australia = Graph[{"WA", "NT", "SA", "Q", "NSW", "V", "T"}, {"WA" <-> "NT", "WA" <-> "SA", "NT" <-> "SA", "NT" <-> "Q", "SA" <-> "Q", "SA" <-> "NSW", "SA" <-> "V", "Q" <-> "NSW", "NSW" <-> "V"}];
Thread[VertexList[australia] -> FindVertexColoring[australia]]
papers = Hypergraph[{{"Erdos", "Renyi"}, {"Erdos", "Graham", "Chung"}, {"Graham", "Knuth", "Patashnik"}, {"Knuth", "Yao"}, {"Chung", "Yao", "Graham"}}];
HypergraphDistance[papers, "Renyi", "Patashnik"]
FindMinimumTransversal[papers]
IsomorphicGraphQ[HypergraphCliqueExpansion[Hypergraph[{{1, 2, 3}}]], HypergraphCliqueExpansion[Hypergraph[{{1, 2}, {2, 3}, {1, 3}}]]]
