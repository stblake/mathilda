# Matrix views of a graph, feeding straight into linear algebra.
g = Graph[{1, 2, 3, 4}, {1 -> 2, 2 -> 3, 3 -> 4, 4 -> 1}];
AdjacencyMatrix[g]
Det[AdjacencyMatrix[g]]
Eigenvalues[AdjacencyMatrix[CycleGraph[4]]]
w = Graph[{1, 2, 3}, {1 -> 2, 2 -> 3}, EdgeWeight -> {5, 7}];
EdgeWeight[w]
WeightedAdjacencyMatrix[w]
EdgeWeight[Graph[{1, 2, 3}, {1 -> 2, 2 -> 3}]]
