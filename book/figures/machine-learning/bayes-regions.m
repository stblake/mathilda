# Gaussian naive Bayes: the boundaries are conic sections.
SeedRandom[7]; blob[m_, s_, n_] := Transpose[{RandomVariate[NormalDistribution[m[[1]], s], n], RandomVariate[NormalDistribution[m[[2]], s], n]}];
{a, b, c} = {blob[{0, 0}, 0.8, 30], blob[{3, 3}, 0.8, 30], blob[{0, 4}, 0.8, 30]}; data = Join[Thread[a -> 1], Thread[b -> 2], Thread[c -> 3]];
picture[f_] := Graphics[{Table[{Hue[{0.6, 0.08, 0.33}[[f[{x, y}]]], 0.2, 1], Rectangle[{x, y}, {x + 0.13, y + 0.13}]}, {x, -2.5, 5.4, 0.125}, {y, -2.5, 6.4, 0.125}], PointSize[0.008], Table[{Hue[{0.6, 0.08, 0.33}[[k]], 0.9, 0.7], Point[{a, b, c}[[k]]]}, {k, 3}]}, Frame -> True]
fig = picture[Classify[data, Method -> "NaiveBayes"]]
