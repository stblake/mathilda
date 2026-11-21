t = AbsoluteTiming[r = Integrate[2/((1 + 2 x) Log[1 + 2 x] (1 + Log[1 + 2 x])^(1/6)), x];][[1]];
Print["loglog ", t, " ", LeafCount[r]];
