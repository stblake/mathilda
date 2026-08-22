# VectorPlot: the phase field of the undamped pendulum, y' = v, v' = -sin(y).
fig = VectorPlot[{y, -Sin[x]}, {x, -3, 3}, {y, -2, 2}]
