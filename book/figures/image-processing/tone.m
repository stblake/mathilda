# format: png
# Figure: a low-contrast image, its linear stretch, and its histogram equalisation.
dim = Image[Table[0.40 + 0.1 Boole[(x - 160)^2 + (y - 192)^2 < 11200] + 0.1 x/576, {y, 384}, {x, 576}]];
fig = ImageAssemble[{{dim}, {ImageAdjust[dim]}, {HistogramTransform[dim]}}]
