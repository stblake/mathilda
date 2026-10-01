# format: png
# Featured: the Ackley function, a nearly flat plateau riddled with small dimples around a
# single deep funnel at the origin.
ack = -20 Exp[-0.2 Sqrt[(x^2 + y^2)/2]] - Exp[(Cos[2 Pi x] + Cos[2 Pi y])/2] + 20 + E;
NMinimize[ack, {x, y}]
fig = Plot3D[Evaluate[ack], {x, -4, 4}, {y, -4, 4}, PlotPoints -> 80]
