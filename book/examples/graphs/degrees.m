# In- and out-degrees of a directed graph.
d = Graph[{1, 2, 3, 4}, {1 -> 2, 2 -> 3, 3 -> 4, 4 -> 1, 1 -> 3}];
EdgeList[d]
VertexInDegree[d]
VertexOutDegree[d]
VertexInDegree[d, 3]
VertexOutDegree[d, 1]
