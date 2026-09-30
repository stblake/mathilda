# Held-out accuracy of k-nearest-neighbours as k grows.
SeedRandom[7]; blob[m_, s_, n_] := Transpose[{RandomVariate[NormalDistribution[m[[1]], s], n], RandomVariate[NormalDistribution[m[[2]], s], n]}];
{a, b, c} = {blob[{0, 0}, 0.8, 30], blob[{3, 3}, 0.8, 30], blob[{0, 4}, 0.8, 30]}; data = Join[Thread[a -> "A"], Thread[b -> "B"], Thread[c -> "C"]];
SeedRandom[8]; test = Join[Thread[blob[{0, 0}, 0.8, 100] -> "A"], Thread[blob[{3, 3}, 0.8, 100] -> "B"], Thread[blob[{0, 4}, 0.8, 100] -> "C"]];
acc[f_, d_] := N[Mean[Boole[Map[f[First[#]] === Last[#] &, d]]]]
curve = Table[{k, acc[Classify[data, Method -> {"NearestNeighbors", "NeighborsNumber" -> k}], test]}, {k, 1, 60}];
fig = ListPlot[curve, Joined -> True]
