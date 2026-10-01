# The running example: three Gaussian clouds of thirty points each.
SeedRandom[7]; blob[m_, s_, n_] := Transpose[{RandomVariate[NormalDistribution[m[[1]], s], n], RandomVariate[NormalDistribution[m[[2]], s], n]}];
{a, b, c} = {blob[{0, 0}, 0.8, 30], blob[{3, 3}, 0.8, 30], blob[{0, 4}, 0.8, 30]}; data = Join[Thread[a -> 1], Thread[b -> 2], Thread[c -> 3]];
fig = ListPlot[{a, b, c}, AspectRatio -> 1]
