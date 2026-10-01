# A perfect matching of the Petersen graph.
pent = Join[Table[0.5 {Cos[Pi/2 + 2 Pi k/5], Sin[Pi/2 + 2 Pi k/5]}, {k, 0, 4}], Table[{Cos[Pi/2 + 2 Pi k/5], Sin[Pi/2 + 2 Pi k/5]}, {k, 0, 4}]] // N;
g = PetersenGraph[];
fig = GraphPlot[g, VertexCoordinates -> pent, GraphHighlight -> FindIndependentEdgeSet[g]]
