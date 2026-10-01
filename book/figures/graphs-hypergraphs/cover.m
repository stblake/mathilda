# A minimum vertex cover (red) and the maximum independent set it leaves.
pent = Join[Table[0.5 {Cos[Pi/2 + 2 Pi k/5], Sin[Pi/2 + 2 Pi k/5]}, {k, 0, 4}], Table[{Cos[Pi/2 + 2 Pi k/5], Sin[Pi/2 + 2 Pi k/5]}, {k, 0, 4}]] // N;
g = PetersenGraph[];
fig = GraphPlot[g, VertexCoordinates -> pent, GraphHighlight -> FindVertexCover[g]]
