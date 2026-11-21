t0 = TimeUsed[]; w0 = AbsoluteTime[];
r = Integrate[(2*E^x^2*x)/(6*(1 + E^x^2)^(2/3) + 3*(1 + E^x^2)), x];
Print["expcube cpu ", TimeUsed[] - t0, " wall ", AbsoluteTime[] - w0, " hash ", Hash[r], " leaves ", LeafCount[r]];
