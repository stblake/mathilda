# One noisy sample of each digit, as the classifier sees it (0 = white, 1 = black).
glyphs = {{1,1,1, 1,0,1, 1,0,1, 1,0,1, 1,1,1}, {0,1,0, 1,1,0, 0,1,0, 0,1,0, 1,1,1}, {1,1,1, 0,0,1, 1,1,1, 1,0,0, 1,1,1}, {1,1,1, 0,0,1, 0,1,1, 0,0,1, 1,1,1}, {1,0,1, 1,0,1, 1,1,1, 0,0,1, 0,0,1}, {1,1,1, 1,0,0, 1,1,1, 0,0,1, 1,1,1}, {1,1,1, 1,0,0, 1,1,1, 1,0,1, 1,1,1}, {1,1,1, 0,0,1, 0,1,0, 0,1,0, 0,1,0}, {1,1,1, 1,0,1, 1,1,1, 1,0,1, 1,1,1}, {1,1,1, 1,0,1, 1,1,1, 0,0,1, 1,1,1}};
SeedRandom[3]; noisy[g_] := Clip[g + RandomVariate[NormalDistribution[0, 0.3], 15], {0, 1}];
tile[d_] := Map[Append[#, 0] &, Partition[noisy[glyphs[[d + 1]]], 3]];
fig = ArrayPlot[Join[Apply[Join, Transpose[Map[tile, Range[0, 4]]], {1}], {ConstantArray[0, 20]}, Apply[Join, Transpose[Map[tile, Range[5, 9]]], {1}]]]
