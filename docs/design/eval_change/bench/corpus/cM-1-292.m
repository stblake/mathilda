t0 = TimeUsed[]; w0 = AbsoluteTime[];
r = Integrate[2/((1 + 2*x)*Log[1 + 2*x]*(1 + Log[1 + 2*x])^(1/6)), x, Method -> "Bronstein"];
Print["cM-1-292 cpu ", TimeUsed[] - t0, " wall ", AbsoluteTime[] - w0, " hash ", Hash[r], " leaves ", LeafCount[r]];
