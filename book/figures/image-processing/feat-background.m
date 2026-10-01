# format: png
# Featured: replace the background. A mask becomes the alpha channel, and the cut-out is
# composited onto a new sky.
fg = Image[Table[{0.9, 0.5 + 0.3 Sin[x/36.], 0.2}, {y, 384}, {x, 512}]];
mask = Image[Table[Boole[(x - 256)^2 + (y - 192)^2 < 22400 && (x - 256)^2 + (y - 120)^2 > 4000], {y, 384}, {x, 512}]];
sky = Image[Table[{0.3, 0.5, 0.6 + 0.4 y/384}, {y, 384}, {x, 512}]];
fig = ImageAssemble[{fg, ImageCompose[sky, SetAlphaChannel[fg, mask]]}]
