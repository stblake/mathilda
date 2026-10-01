# A noisy sine wave and the 5-nearest-neighbour predictor trained on it.
SeedRandom[9]; xs = Sort[RandomReal[{0, 2 Pi}, 120]]; ys = Sin[xs] + RandomVariate[NormalDistribution[0, 0.3], 120];
smoother = Predict[Thread[xs -> ys], Method -> {"NearestNeighbors", "NeighborsNumber" -> 5}];
fig = Graphics[{PointSize[0.008], GrayLevel[0.45], Point[Transpose[{xs, ys}]], Thickness[0.004], Hue[0.08, 0.9, 0.85], Line[Table[{x, smoother[x]}, {x, 0., N[2 Pi], 0.01}]], Hue[0.6, 0.9, 0.7], Line[Table[{x, Sin[x]}, {x, 0., N[2 Pi], 0.05}]]}, Frame -> True]
