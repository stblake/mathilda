# Isomorphism, canonical forms, symmetry and planarity.
sq = Graph[{a <-> b, b <-> d, d <-> c, c <-> a}]
IsomorphicGraphQ[CycleGraph[4], sq]
FindGraphIsomorphism[CycleGraph[4], sq]
CanonicalGraph[CycleGraph[4]] === CanonicalGraph[sq]
GraphAutomorphismGroup[PetersenGraph[]]
PlanarGraphQ /@ {CompleteGraph[4], CompleteGraph[5], CompleteGraph[{3, 3}], GridGraph[{3, 3}]}
PlanarGraphQ[PetersenGraph[]]
