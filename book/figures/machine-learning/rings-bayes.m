# Two concentric rings. Naive Bayes separates them; a linear model cannot.
SeedRandom[5]; ring[r_, n_] := Table[With[{t = RandomReal[{0, 2 Pi}], s = r + RandomVariate[NormalDistribution[0, 0.25]]}, {s Cos[t], s Sin[t]}], {n}];
{in, out} = {ring[1, 100], ring[2.5, 100]}; rd = Join[Thread[in -> 1], Thread[out -> 2]];
nb = Classify[rd, Method -> "NaiveBayes"];
fig = Graphics[{Table[{Hue[{0.6, 0.08}[[nb[{x, y}]]], 0.2, 1], Rectangle[{x, y}, {x + 0.11, y + 0.11}]}, {x, -3.5, 3.4, 0.1}, {y, -3.5, 3.4, 0.1}], PointSize[0.012], Hue[0.6, 0.9, 0.7], Point[in], Hue[0.08, 0.9, 0.7], Point[out]}, Frame -> True, AspectRatio -> 1]
