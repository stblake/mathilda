t0 = TimeUsed[]; w0 = AbsoluteTime[];
r = Integrate[((2 + x)/(-3 + x))^(2/3)/x, x];
Print["trager-ratradicand-007 cpu ", TimeUsed[] - t0, " wall ", AbsoluteTime[] - w0, " hash ", Hash[r], " leaves ", LeafCount[r]];
