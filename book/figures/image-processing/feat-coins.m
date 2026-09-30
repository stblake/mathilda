# format: png
# Featured: count the coins on a tray. Twelve discs on a ribbed background.
SeedRandom[21]; centres = 3 Flatten[Table[{45 i, 45 j} - 10 + RandomReal[{-5, 5}, 2], {i, 4}, {j, 3}], 1];
tray = Image[Table[0.3 + 0.05 Sin[x/9.] + 0.5 Boole[Min[Norm[{x, y} - #] & /@ centres] < 39], {y, 450}, {x, 570}]];
count = Max[MorphologicalComponents[Binarize[tray]]]
fig = tray
