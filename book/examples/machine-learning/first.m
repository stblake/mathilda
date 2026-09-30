# The running example: three labelled clouds and the default classifier.
SeedRandom[7]; blob[m_, s_, n_] := Transpose[{RandomVariate[NormalDistribution[m[[1]], s], n], RandomVariate[NormalDistribution[m[[2]], s], n]}];
{a, b, c} = {blob[{0, 0}, 0.8, 30], blob[{3, 3}, 0.8, 30], blob[{0, 4}, 0.8, 30]};
data = Join[Thread[a -> "A"], Thread[b -> "B"], Thread[c -> "C"]];
Length[data]
Take[data, 2]
nn = Classify[data]
nn[{0.5, 0.5}]
nn[{1.5, 2.}, "Probabilities"]
Map[nn, {{0, 0}, {3, 3}, {0, 4}}]
{nn["Classes"], nn["Method"], nn["FeatureCount"], nn["NeighborCount"]}
