# Two elongated, tilted clouds: k-means splits them wrongly.
SeedRandom[2]; z = Partition[RandomVariate[NormalDistribution[], 600], 2];
pts = Join[Map[{{1.2, 0.}, {0.9, 0.35}} . # &, Take[z, 150]], Map[{{1.2, 0.}, {-0.9, 0.35}} . # + {1.5, 3.} &, Take[z, -150]]];
fig = ListPlot[FindClusters[pts, 2, Method -> "KMeans"], AspectRatio -> 1]
