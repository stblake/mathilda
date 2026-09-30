# The giant component of a random graph just above the threshold.
SeedRandom[5]; r = RandomGraph[{120, 80}];
giant = First[TakeLargestBy[ConnectedComponents[r], Length, 1]];
Length[giant]
fig = GraphPlot[r, GraphHighlight -> giant]
