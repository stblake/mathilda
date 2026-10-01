# Graphs are expressions: build them from data, rewrite them with rules.
follows = <|"alice" -> {"bob", "carol"}, "bob" -> {"carol"}, "dave" -> {"alice"}|>;
social = Graph[Flatten[KeyValueMap[Thread[#1 -> #2] &, follows]]]
EdgeList[social]
Select[VertexList[social], VertexInDegree[social, #] == 0 &]
EdgeRules[PathGraph[{1, 2, 3}]]
EdgeList[CycleGraph[4]] /. UndirectedEdge[u_, v_] :> UndirectedEdge[u, v + 10]
EdgeList[NeighborhoodGraph[PetersenGraph[], 1]]
Subgraph[PetersenGraph[], {1, 2, 3, 4, 5}]
