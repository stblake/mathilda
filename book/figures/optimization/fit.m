# Featured: fit a decaying exponential a e^(-k t) + c to noisy measurements by
# minimizing the sum of squared residuals.
SeedRandom[42]; data = Table[{t, 3 Exp[-0.7 t] + 0.5 + RandomReal[{-0.1, 0.1}]}, {t, 0, 6, 0.25}];
sse = Total[(a Exp[-k #[[1]]] + c - #[[2]])^2 & /@ data];
sol = FindMinimum[sse, {{a, 1}, {k, 1}, {c, 0}}]
model = a Exp[-k t] + c /. Last[sol];
fig = Graphics[{First[Plot[model, {t, 0, 6}]], Red, PointSize[0.012], Point[data]}, Frame -> True]
