# Generators: the star, and reproducible random graphs.
StarGraph[5]
EdgeList[StarGraph[5]]
VertexDegree[StarGraph[5]]
SeedRandom[42];
r = RandomGraph[{6, 8}]
EdgeList[r]
Length[RandomGraph[{6, 5}, 3]]
