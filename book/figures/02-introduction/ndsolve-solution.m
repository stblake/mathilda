# Figure/cell: a damped harmonic oscillator solved by NDSolve, then plotted.
# gen_figures.py shows this as a REPL session -- In[1] the NDSolve, In[2] the
# Plot -- and renders the plot as Out[2] (the Export-ed PDF). The final plot
# MUST be bound to `fig`; `fig = <expr>` is shown to the reader as just <expr>.
sol = NDSolve[{y''[x] + 0.3 y'[x] + 4 y[x] == 0, y[0] == 1, y'[0] == 0}, y, {x, 0, 25}]
fig = Plot[Evaluate[y[x] /. sol], {x, 0, 25}]
