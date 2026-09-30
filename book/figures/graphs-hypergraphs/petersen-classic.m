# The classic drawing: an outer pentagon and an inner pentagram.
pent = Join[Table[0.5 {Cos[Pi/2 + 2 Pi k/5], Sin[Pi/2 + 2 Pi k/5]}, {k, 0, 4}], Table[{Cos[Pi/2 + 2 Pi k/5], Sin[Pi/2 + 2 Pi k/5]}, {k, 0, 4}]] // N;
fig = GraphPlot[PetersenGraph[], VertexCoordinates -> pent]
