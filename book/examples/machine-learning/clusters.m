# FindClusters: gaps on a line, spanning trees in the plane, k-means, and labels for new points.
FindClusters[{1, 2, 3, 10, 11, 12, 25}]
{FindClusters[{1, 2, 3, 10, 11, 12, 25}, 2], FindClusters[{1, 2, 10, 12, 3, 1, 13, 25}, UpTo[4]]}
SeedRandom[7]; blob[m_, s_, n_] := Transpose[{RandomVariate[NormalDistribution[m[[1]], s], n], RandomVariate[NormalDistribution[m[[2]], s], n]}];
pts = Join[blob[{0, 0}, 0.8, 30], blob[{3, 3}, 0.8, 30], blob[{0, 4}, 0.8, 30]];
Map[Length, FindClusters[pts]]
Map[Length, FindClusters[pts, 3]]
km = FindClusters[pts, 3, Method -> "KMeans"]; Map[Length, km]
Map[Mean, km]
FindClusters[pts, Method -> "KMeans"] // Head
labelled = Flatten[Table[Thread[km[[i]] -> i], {i, 3}]];
assign = Classify[labelled];
Map[assign, {{0, 0}, {3, 3}, {0, 4}}]
