# Matchings, covers, independent sets and cliques.
g = PetersenGraph[];
FindIndependentEdgeSet[g]
cover = FindVertexCover[g]
indep = First[FindIndependentVertexSet[g]]
Sort[Join[cover, indep]] == VertexList[g]
FindClique[g]
FindClique[CompleteGraph[5]]
Max[FindVertexColoring[g]]
