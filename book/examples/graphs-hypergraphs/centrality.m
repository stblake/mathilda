# Centrality and clustering on a triangle with a pendant vertex.
k = Graph[{1 <-> 2, 2 <-> 3, 3 <-> 1, 3 <-> 4}]
DegreeCentrality[k]
ClosenessCentrality[k]
BetweennessCentrality[k]
PageRankCentrality[k]
EigenvectorCentrality[k]
LocalClusteringCoefficient[k]
{GlobalClusteringCoefficient[k], MeanClusteringCoefficient[k]}
