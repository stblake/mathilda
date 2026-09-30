# format: png
# Figure: the Rastrigin function on [-5.12, 5.12]^2, an egg carton of local minima around
# a single global minimum at the origin.
fig = DensityPlot[20 + x^2 - 10 Cos[2 Pi x] + y^2 - 10 Cos[2 Pi y], {x, -5.12, 5.12}, {y, -5.12, 5.12}, PlotPoints -> 120]
