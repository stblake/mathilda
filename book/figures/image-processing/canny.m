# format: png
# Figure: Canny edges of a noisy image with no smoothing, and with radius 4.
shapes = Image[Table[0.2 + 0.6 Boole[(x - 160)^2 + (y - 192)^2 < 11200] + 0.3 Boole[Abs[x - 400] < 80 && Abs[y - 192] < 120], {y, 384}, {x, 576}]];
SeedRandom[2]; noisy = Image[ImageData[shapes] + RandomReal[{-0.15, 0.15}, {384, 576}]];
fig = ImageAssemble[{{noisy}, {EdgeDetect[noisy, 0]}, {EdgeDetect[noisy, 4]}}]
