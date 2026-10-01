# format: png
# Figure: the exact Euclidean distance transform of two shapes, rescaled for display.
mask = Image[Table[Boole[(x - 160)^2 + (y - 192)^2 < 11200 || Abs[x - 400] < 80 && Abs[y - 192] < 120], {y, 384}, {x, 576}]];
fig = ImageAssemble[{mask, ImageAdjust[DistanceTransform[mask]]}]
