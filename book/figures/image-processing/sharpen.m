# format: png
# Figure: a blurred image, then sharpened by a Laplacian-type kernel. Enlarged with "Nearest"
# so the pixel-scale effect is visible.
shapes = Image[Table[0.2 + 0.6 Boole[(x - 40)^2 + (y - 48)^2 < 700] + 0.3 Boole[Abs[x - 100] < 20 && Abs[y - 48] < 30], {y, 96}, {x, 144}]];
soft = GaussianFilter[shapes, 3];
fig = ImageResize[ImageAssemble[{soft, ImageConvolve[soft, {{0, -1, 0}, {-1, 5, -1}, {0, -1, 0}}]}], 1152, Resampling -> "Nearest"]
