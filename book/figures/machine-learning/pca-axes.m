# A correlated cloud with its two principal axes, each drawn two standard deviations long.
SeedRandom[4]; cloud = Map[{{2., 0.}, {1.2, 0.5}} . # + {5., 3.} &, Partition[RandomVariate[NormalDistribution[], 400], 2]];
{vals, vecs} = {Eigenvalues[Covariance[cloud]], Eigenvectors[Covariance[cloud]]}; m = Mean[cloud];
fig = Graphics[{PointSize[0.008], GrayLevel[0.5], Point[cloud], Thickness[0.006], Hue[0.0, 0.9, 0.8], Table[Line[{m - 2 Sqrt[vals[[i]]] vecs[[i]], m + 2 Sqrt[vals[[i]]] vecs[[i]]}], {i, 2}]}, Frame -> True, PlotRange -> {{0, 10}, {-0.5, 6.4}}]
