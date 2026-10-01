# Trees and Cayley's formula.
Table[Det[Drop[KirchhoffMatrix[CompleteGraph[n]], {1}, {1}]], {n, 2, 8}]
Table[n^(n - 2), {n, 2, 8}]
t = CompleteKaryTree[3, 2]
{TreeGraphQ[t], VertexCount[t], EdgeCount[t]}
FindSpanningTree[PetersenGraph[]]
TreeGraphQ[FindSpanningTree[PetersenGraph[]]]
