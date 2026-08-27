# Figure/cell: the heat equation u_t = u_xx solved by the method of lines. The initial
# sine profile decays; each curve is the temperature along the bar at a later instant.
sol = NDSolve[{D[u[t, x], t] == D[u[t, x], {x, 2}], u[0, x] == Sin[Pi x], u[t, 0] == 0, u[t, 1] == 0}, u, {t, 0, 0.1}, {x, 0, 1}]
fig = Plot[Evaluate[Table[u[t, x] /. First[sol], {t, 0, 0.1, 0.025}]], {x, 0, 1}]
