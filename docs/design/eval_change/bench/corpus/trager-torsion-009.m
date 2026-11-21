t0 = TimeUsed[]; w0 = AbsoluteTime[];
r = Integrate[((1 + x^2)*(2 + x^4)^(1/6))/x, x, Method -> "Trager"];
Print["trager-torsion-009 cpu ", TimeUsed[] - t0, " wall ", AbsoluteTime[] - w0, " hash ", Hash[r], " leaves ", LeafCount[r]];
