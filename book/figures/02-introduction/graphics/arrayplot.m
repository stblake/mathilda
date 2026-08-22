# format: png
# ArrayPlot: the GCD table -- shading reveals the number-theoretic structure.
fig = ArrayPlot[Table[GCD[i, j], {i, 1, 30}, {j, 1, 30}]]
