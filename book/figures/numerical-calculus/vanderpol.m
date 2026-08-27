# Figure/cell: the van der Pol relaxation oscillator, drawn as its phase-plane limit cycle.
# gen_figures.py shows In[1] the NDSolve and In[2] the plot command, then renders the plot
# as Out[2]. The final plot MUST be bound to `fig`; `fig = <expr>` is shown as just <expr>.
sol = NDSolve[{y''[x] - 3 (1 - y[x]^2) y'[x] + y[x] == 0, y[0] == 1, y'[0] == 0}, y, {x, 0, 20}]
fig = ParametricPlot[Evaluate[{y[x], y'[x]} /. sol], {x, 0, 20}]
