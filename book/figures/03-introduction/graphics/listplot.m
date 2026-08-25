# ListPlot: a sampled, slowly beating signal as discrete points.
fig = ListPlot[Table[Sin[0.5 n] + Sin[0.6 n], {n, 0, 80}]]
