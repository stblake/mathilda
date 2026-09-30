# Directed graphs: dependency order, strong components, orientation.
dress = Graph[{socks -> shoes, trousers -> shoes, trousers -> belt, shirt -> belt, shirt -> tie, tie -> jacket, belt -> jacket}]
TopologicalSort[dress]
AcyclicGraphQ[dress]
AcyclicGraphQ[EdgeAdd[dress, jacket -> shirt]]
d = Graph[{1 -> 2, 2 -> 3, 3 -> 1, 3 -> 4, 4 -> 5, 5 -> 4, 6 -> 5}];
StronglyConnectedComponents[d]
EdgeList[ReverseGraph[Graph[{1 -> 2, 2 -> 3}]]]
b = IncidenceMatrix[Graph[{1 -> 2, 2 -> 3, 1 -> 3, 3 -> 4}]]
b . Transpose[b] == KirchhoffMatrix[UndirectedGraph[Graph[{1 -> 2, 2 -> 3, 1 -> 3, 3 -> 4}]]]
