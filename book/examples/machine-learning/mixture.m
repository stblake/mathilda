# Elongated clusters: k-means against a Gaussian mixture fitted by EM.
SeedRandom[2]; z = Partition[RandomVariate[NormalDistribution[], 600], 2];
g1 = Map[{{1.2, 0.}, {0.9, 0.35}} . # &, Take[z, 150]]; g2 = Map[{{1.2, 0.}, {-0.9, 0.35}} . # + {1.5, 3.} &, Take[z, -150]]; pts = Join[g1, g2];
km = FindClusters[pts, 2, Method -> "KMeans"]; gm = FindClusters[pts, Method -> "GaussianMixture"];
{Map[Length, km], Map[Length, gm]}
{Length[Intersection[km[[1]], g1]], Length[Intersection[gm[[1]], g1]]}
d = LearnDistribution[pts, Method -> "GaussianMixture"]
d[[2, 1]]
d[[2, 2]]
PDF[d, {{0., 0.}, {1.5, 3.}, {10., 10.}}]
m = LearnDistribution[pts]
{PDF[m, {0.75, 1.5}], PDF[d, {0.75, 1.5}]}
ct = LearnDistribution[{"r", "r", "r", "b"}, Method -> "ContingencyTable"]
{PDF[ct, "r"], PDF[ct, "b"], PDF[ct, "g"]}
