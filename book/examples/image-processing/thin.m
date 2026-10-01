# Thinning a 3-pixel-thick bar to its one-pixel skeleton, then pruning its ends.
bar = Image[Table[If[2 <= y <= 4 && 2 <= x <= 10, 1, 0], {y, 5}, {x, 11}]];
ImageData[Thinning[bar], "Bit"]
ImageData[Pruning[Thinning[bar], 2], "Bit"]
ImageData[Pruning[Image[{{0, 0, 0}, {0, 1, 0}, {0, 0, 0}}]], "Bit"]
ImageData[Thinning[Image[{{0.2, 0.7, 0.4}}]], "Bit"]
