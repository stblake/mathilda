# Australia's regions, coloured so that neighbours differ.
australia = Graph[{"WA", "NT", "SA", "Q", "NSW", "V", "T"}, {"WA" <-> "NT", "WA" <-> "SA", "NT" <-> "SA", "NT" <-> "Q", "SA" <-> "Q", "SA" <-> "NSW", "SA" <-> "V", "Q" <-> "NSW", "NSW" <-> "V"}];
cols = {RGBColor[0.86, 0.3, 0.25], RGBColor[0.3, 0.6, 0.3], RGBColor[0.3, 0.45, 0.8]};
fig = GraphPlot[australia, VertexCoordinates -> {"WA" -> {0, 1}, "NT" -> {1.2, 1.8}, "SA" -> {1.3, 0.8}, "Q" -> {2.5, 1.7}, "NSW" -> {2.6, 0.7}, "V" -> {2.3, 0.1}, "T" -> {2.4, -0.7}}, VertexLabels -> "Name", VertexStyle -> Thread[VertexList[australia] -> cols[[FindVertexColoring[australia]]]], VertexSize -> 0.2]
