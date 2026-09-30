# format: png
# Figure: the Shi-Tomasi corner response of the shapes image, rescaled for display.
shapes = Image[Table[0.2 + 0.6 Boole[(x - 160)^2 + (y - 192)^2 < 11200] + 0.3 Boole[Abs[x - 400] < 80 && Abs[y - 192] < 120], {y, 384}, {x, 576}]];
fig = ImageAssemble[{shapes, ImageAdjust[CornerFilter[shapes, 8]]}]
