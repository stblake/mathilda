t0 = TimeUsed[]; w0 = AbsoluteTime[];
r = Integrate[(1 + E^x)^(7/6)/(-2 + E^x), x];
Print["cM-3-21 cpu ", TimeUsed[] - t0, " wall ", AbsoluteTime[] - w0, " hash ", Hash[r], " leaves ", LeafCount[r]];
