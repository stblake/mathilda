# Figure: the NDSolve solution of the van der Pol oscillator, plotted over [0, 20].
# gen_figures.py appends `Export["generated/figures/.../ndsolve-solution.pdf", fig]`
# and runs this through the real binary, so the figure is produced the same way
# every transcript is. Must bind the plot to `fig`.
sol = NDSolve[{y''[x] + 3 (y[x]^2 - 1) y'[x] + y[x] == 0, y[0] == 1, y'[0] == 0}, y, {x, 0, 20}];
fig = Plot[Evaluate[y[x] /. sol], {x, 0, 20}];
