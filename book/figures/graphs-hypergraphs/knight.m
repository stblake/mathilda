# A closed knight's tour, drawn on the board.
moves = {{1, 2}, {2, 1}, {-1, 2}, {-2, 1}};
knight[n_] := Graph[Flatten[Table[If[1 <= i + m[[1]] <= n && 1 <= j + m[[2]] <= n, {i, j} <-> {i + m[[1]], j + m[[2]]}, Nothing], {i, n}, {j, n}, {m, moves}]]]
k8 = knight[8];
fig = GraphPlot[k8, VertexCoordinates -> Thread[VertexList[k8] -> N[VertexList[k8]]], EdgeStyle -> GrayLevel[0.85], GraphHighlight -> First[FindHamiltonianCycle[k8]], VertexSize -> 0.15]
