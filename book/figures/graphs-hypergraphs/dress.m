# The dressing order as a layered drawing.
dress = Graph[{socks -> shoes, trousers -> shoes, trousers -> belt, shirt -> belt, shirt -> tie, tie -> jacket, belt -> jacket}];
fig = GraphPlot[dress, VertexLabels -> "Name"]
