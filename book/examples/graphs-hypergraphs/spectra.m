# Laplacian and the matrix-tree theorem.
g = PetersenGraph[];
KirchhoffMatrix[CycleGraph[4]]
Eigenvalues[KirchhoffMatrix[g]]
Det[Drop[KirchhoffMatrix[g], {1}, {1}]]
Tally[Eigenvalues[AdjacencyMatrix[g]]]
