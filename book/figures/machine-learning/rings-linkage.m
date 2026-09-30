# Single linkage follows each ring around.
SeedRandom[5]; ring[r_, n_] := Table[With[{t = RandomReal[{0, 2 Pi}], s = r + RandomVariate[NormalDistribution[0, 0.1]]}, {s Cos[t], s Sin[t]}], {n}];
rings = Join[ring[1, 80], ring[3, 160]];
fig = ListPlot[FindClusters[rings, 2], AspectRatio -> 1]
