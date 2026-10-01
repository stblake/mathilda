# Figure: the raw differential-evolution best (no local polish) on Rastrigin, as log10 of the
# objective against the generation budget MaxIterations.
rast = 20 + x^2 - 10 Cos[2 Pi x] + y^2 - 10 Cos[2 Pi y];
best[k_] := First[NMinimize[{rast, -5.12 <= x <= 5.12 && -5.12 <= y <= 5.12}, {x, y}, Method -> {"DifferentialEvolution", "PostProcess" -> False}, MaxIterations -> k]];
data = Table[{k, Log[10, best[k]]}, {k, 1, 100}];
fig = ListPlot[data, Joined -> True]
