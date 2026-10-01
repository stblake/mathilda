# Editing, transforms and set operations.
VertexDelete[CycleGraph[5], 1]
EdgeList[EdgeAdd[PathGraph[{1, 2, 3}], 3 <-> 1]]
EdgeList[GraphComplement[CycleGraph[5]]]
IsomorphicGraphQ[CycleGraph[5], GraphComplement[CycleGraph[5]]]
EdgeList[LineGraph[StarGraph[4]]]
EdgeList[GraphUnion[PathGraph[{1, 2, 3}], PathGraph[{3, 4}]]]
