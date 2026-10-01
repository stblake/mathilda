# format: png
# Featured: binarize a page photographed under uneven light. Glyph blocks on a lighting ramp;
# Otsu's single threshold, then a local threshold, which recovers the ink exactly.
SeedRandom[8]; ink = RandomInteger[1, {8, 24}];
page = Image[Table[(0.25 + 0.7 x/192) (1 - 0.6 Boole[Mod[y, 12] > 5 && Mod[x, 8] > 1 && ink[[Quotient[y, 12] + 1, Quotient[x, 8] + 1]] == 1]), {y, 96}, {x, 192}]];
fig = ImageResize[ImageAssemble[{{page}, {Binarize[page]}, {LocalAdaptiveBinarize[page, 8, {1, 0, -0.05}]}}], 768, Resampling -> "Nearest"]
