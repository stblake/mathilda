# Hypergraph expansions, line graphs, components and distances.
h = Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}];
EdgeList[HypergraphCliqueExpansion[h]]
HypergraphStarExpansion[h]
EdgeList[HypergraphToGraph[Hypergraph[{{a, b, c}}]]]
HypergraphConnectedComponents[h]
{HypergraphDistance[h, 1, 6], HyperedgeDistance[h, 1, 3]}
f = Hypergraph[{{1, 2, 3}, {2, 3, 4}, {3, 4, 5}, {1, 5}}];
EdgeList[HypergraphLineGraph[f]]
EdgeList[HypergraphLineGraph[f, 2]]
HyperedgeConnectedComponents[f, 2]
