# k-means with k = 3 on the unlabelled blobs; each cluster in its own colour.
SeedRandom[7]; blob[m_, s_, n_] := Transpose[{RandomVariate[NormalDistribution[m[[1]], s], n], RandomVariate[NormalDistribution[m[[2]], s], n]}];
pts = Join[blob[{0, 0}, 0.8, 30], blob[{3, 3}, 0.8, 30], blob[{0, 4}, 0.8, 30]];
clusters = FindClusters[pts, 3, Method -> "KMeans"];
fig = ListPlot[clusters, AspectRatio -> 1]
