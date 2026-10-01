# Least squares by Predict, checked two independent ways; then k-NN regression.
SeedRandom[11]; xs = RandomReal[{0, 10}, 40]; ys = 2.5 + 0.8 xs + RandomVariate[NormalDistribution[0, 1], 40];
lin = Predict[Thread[xs -> ys]]
{lin["Coefficients"], lin[5.]}
Fit[Transpose[{xs, ys}], {1, x}, x]
LeastSquares[Transpose[{ConstantArray[1., 40], xs}], ys]
Predict[Transpose[{xs, ys}]]["Coefficients"]
Predict[{{1., 1., 6.}, {2., 1., 8.}, {1., 2., 9.}, {3., 2., 13.}, {2., 3., 14.}}]["Coefficients"]
Predict[{{1., 2., 5.}, {2., 4., 7.}, {3., 6., 9.}, {4., 8., 11.}}]
knn = Predict[Thread[xs -> ys], Method -> "NearestNeighbors"]
{knn[5.], knn["NeighborCount"]}
Mean[Take[SortBy[Transpose[{xs, ys}], Abs[First[#] - 5.] &], 3][[All, 2]]]
Predict[Thread[xs -> ys], Method -> "RandomForest"] // Head
