# Random graphs: the giant component and the degree distribution.
SeedRandom[1]; Table[{m, Max[Length /@ ConnectedComponents[RandomGraph[{1000, m}]]]}, {m, {250, 400, 500, 600, 750, 1000, 1500}}]
SeedRandom[2]; r = RandomGraph[{200, 400}];
Mean[N[VertexDegree[r]]]
SortBy[Tally[VertexDegree[r]], First]
Table[Round[200 Exp[-4.] 4^k/k!], {k, 0, 11}]
