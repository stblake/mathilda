t1 = AbsoluteTiming[r1 = Integrate[((1 + x^2)*(2 + x^4)^(1/6))/x, x, Method -> "Trager"];][[1]];
Print["trager ", t1, " ", LeafCount[r1]];
t2 = AbsoluteTiming[r2 = Integrate[2/((1 + 2 x) Log[1 + 2 x] (1 + Log[1 + 2 x])^(1/6)), x];][[1]];
Print["loglog ", t2, " ", LeafCount[r2]];
t3 = AbsoluteTiming[r3 = Integrate[(2*E^x^2*x)/(6*(1 + E^x^2)^(2/3) + 3*(1 + E^x^2)), x];][[1]];
Print["expcube ", t3, " ", LeafCount[r3]];
Print[$ModuleNumber];
