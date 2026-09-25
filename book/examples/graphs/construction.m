# Construction: build a directed 4-cycle and read it back.
g = Graph[{1, 2, 3, 4}, {1 -> 2, 2 -> 3, 3 -> 4, 4 -> 1}]
InputForm[g]
VertexList[g]
EdgeList[g]
EdgeCount[g]
Graph[{a -> b, b -> c, c -> a}]
EdgeList[Graph[{1, 2, 3}, {1 <-> 2, 2 <-> 3}]]
GraphQ[Graph[{1}, {1 -> 1}]]
