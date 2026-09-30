# Figure: the quadratic penalty x^2 + mu max(0, 1 - x)^2 for minimize x^2 subject to x >= 1,
# with mu = 1, 10, 100. The minimizer mu/(1 + mu) creeps toward the constraint x = 1.
fig = Plot[Evaluate[Table[x^2 + mu Max[0, 1 - x]^2, {mu, {1, 10, 100}}]], {x, 0, 1.4}, PlotRange -> {{0, 1.4}, {0, 2.2}}]
