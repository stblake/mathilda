# A minimum proper colouring: three colours.
pent = Join[Table[0.5 {Cos[Pi/2 + 2 Pi k/5], Sin[Pi/2 + 2 Pi k/5]}, {k, 0, 4}], Table[{Cos[Pi/2 + 2 Pi k/5], Sin[Pi/2 + 2 Pi k/5]}, {k, 0, 4}]] // N;
g = PetersenGraph[];
fig = GraphPlot[g, VertexCoordinates -> pent, VertexStyle -> Thread[VertexList[g] -> ({RGBColor[0.86, 0.3, 0.25], RGBColor[0.3, 0.6, 0.3], RGBColor[0.3, 0.45, 0.8]}[[#]] & /@ FindVertexColoring[g])], VertexSize -> 0.12]
