t = AbsoluteTiming[r = Integrate[(2*E^x^2*x)/(6*(1 + E^x^2)^(2/3) + 3*(1 + E^x^2)), x];][[1]];
Print["expcube ", t, " ", LeafCount[r], " ", InputForm[r]];
