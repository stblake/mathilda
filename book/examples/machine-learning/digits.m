# Ten digits drawn on a 3 x 5 grid, flattened to 15 pixels; noisy copies are recognised by nearest neighbours.
glyphs = {{1,1,1, 1,0,1, 1,0,1, 1,0,1, 1,1,1}, {0,1,0, 1,1,0, 0,1,0, 0,1,0, 1,1,1}, {1,1,1, 0,0,1, 1,1,1, 1,0,0, 1,1,1}, {1,1,1, 0,0,1, 0,1,1, 0,0,1, 1,1,1}, {1,0,1, 1,0,1, 1,1,1, 0,0,1, 0,0,1}, {1,1,1, 1,0,0, 1,1,1, 0,0,1, 1,1,1}, {1,1,1, 1,0,0, 1,1,1, 1,0,1, 1,1,1}, {1,1,1, 0,0,1, 0,1,0, 0,1,0, 0,1,0}, {1,1,1, 1,0,1, 1,1,1, 1,0,1, 1,1,1}, {1,1,1, 1,0,1, 1,1,1, 0,0,1, 1,1,1}};
SeedRandom[3]; noisy[g_] := Clip[g + RandomVariate[NormalDistribution[0, 0.3], 15], {0, 1}];
reader = Classify[Flatten[Table[Thread[Table[noisy[glyphs[[d + 1]]], {20}] -> d], {d, 0, 9}]]];
test = Flatten[Table[Thread[Table[noisy[glyphs[[d + 1]]], {30}] -> d], {d, 0, 9}]];
N[Count[Map[reader[First[#]] === Last[#] &, test], True]/Length[test]]
reader[noisy[glyphs[[9]]], "Probabilities"]
mistakes = Select[Map[{Last[#], reader[First[#]]} &, test], First[#] =!= Last[#] &]; Take[Reverse[SortBy[Tally[mistakes], Last]], 4]
