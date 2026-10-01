# format: png
# Figure: a checkerboard under a lighting ramp; Otsu's global threshold, then the local one.
lit = Image[Table[(0.1 + 0.8 x/768) (0.6 + 0.4 Mod[Quotient[x, 64] + Quotient[y, 64], 2]), {y, 256}, {x, 768}]];
fig = ImageAssemble[{{lit}, {Binarize[lit]}, {LocalAdaptiveBinarize[lit, 32]}}]
