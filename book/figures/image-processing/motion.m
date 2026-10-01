# format: png
# Figure: a diagonal motion blur, a 61x61 kernel of rank 61 that cannot be factored.
shapes = Image[Table[0.2 + 0.6 Boole[(x - 160)^2 + (y - 192)^2 < 11200] + 0.3 Boole[Abs[x - 400] < 80 && Abs[y - 192] < 120], {y, 384}, {x, 576}]];
fig = ImageAssemble[{shapes, ImageConvolve[shapes, IdentityMatrix[61]/61.]}]
