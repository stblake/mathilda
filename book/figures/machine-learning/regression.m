# Forty noisy points, the least-squares line, and the 3-nearest-neighbour predictor.
SeedRandom[11]; xs = RandomReal[{0, 10}, 40]; ys = 2.5 + 0.8 xs + RandomVariate[NormalDistribution[0, 1], 40];
{lin, knn} = {Predict[Thread[xs -> ys]], Predict[Thread[xs -> ys], Method -> "NearestNeighbors"]};
fig = Graphics[{PointSize[0.012], Point[Transpose[{xs, ys}]], Hue[0.6, 0.9, 0.7], Line[{{0, lin[0.]}, {10, lin[10.]}}], Hue[0.08, 0.9, 0.8], Line[Table[{x, knn[x]}, {x, 0., 10., 0.02}]]}, Frame -> True]
