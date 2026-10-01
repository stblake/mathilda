# format: png
# Featured: separate touching cells. Nine overlapping discs merge into six components; the cores
# of the distance transform (distance > 30) separate them again.
SeedRandom[21]; c = 3 RandomReal[{14, 146}, {9, 2}];
cells = Image[Table[Boole[Min[Norm[{x, y} - #] & /@ c] < 36], {y, 480}, {x, 480}]];
Max[MorphologicalComponents[cells]]
cores = Binarize[DistanceTransform[cells], 30];
Max[MorphologicalComponents[cores]]
fig = ImageAssemble[{cells, ImageAdjust[DistanceTransform[cells]], cores}]
