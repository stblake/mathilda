# Hypergraphs: construction and the shared accessors.
h = Hypergraph[{{1, 2, 3}, {3, 4}, {4, 5, 6}, {7}}]
InputForm[h]
VertexDegree[h]
HyperedgeSizes[h]
{HypergraphRank[h], HypergraphCorank[h], UniformHypergraphQ[h]}
IncidenceMatrix[h]
EdgeList[HypergraphDual[h]]
