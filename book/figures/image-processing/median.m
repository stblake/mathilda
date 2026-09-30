# format: png
# Figure: salt-and-pepper noise, then MeanFilter and MedianFilter at radius 1. Enlarged with
# "Nearest" so single noisy pixels stay visible.
shapes = Image[Table[0.2 + 0.6 Boole[(x - 40)^2 + (y - 48)^2 < 700] + 0.3 Boole[Abs[x - 100] < 20 && Abs[y - 48] < 30], {y, 96}, {x, 144}]];
SeedRandom[5]; salt = RandomChoice[{0.9, 0.05, 0.05} -> {0, 1, -1}, {96, 144}];
noisy = Image[Clip[ImageData[shapes] + salt, {0, 1}]];
fig = ImageResize[ImageAssemble[{{noisy}, {MeanFilter[noisy, 1]}, {MedianFilter[noisy, 1]}}], 576, Resampling -> "Nearest"]
