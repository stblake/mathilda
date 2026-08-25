# format: png
# ParametricPlot3D: a coil wound around a torus.
fig = ParametricPlot3D[{(2 + Cos[6 t]) Cos[t], (2 + Cos[6 t]) Sin[t], Sin[6 t]}, {t, 0, 2 Pi}]
