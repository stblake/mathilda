# Contours of the fitted two-component mixture density.
SeedRandom[2]; z = Partition[RandomVariate[NormalDistribution[], 600], 2];
pts = Join[Map[{{1.2, 0.}, {0.9, 0.35}} . # &, Take[z, 150]], Map[{{1.2, 0.}, {-0.9, 0.35}} . # + {1.5, 3.} &, Take[z, -150]]];
d = LearnDistribution[pts, Method -> "GaussianMixture"];
fig = ContourPlot[PDF[d, {x, y}], {x, -4, 6}, {y, -2.5, 5.5}, AspectRatio -> 0.8]
