t0 = TimeUsed[]; w0 = AbsoluteTime[];
r = Integrate[((1 + E^x^2)*(1 + E^x^2 + E^(2*x^2))^(3/2)*x)/E^x^2, x, Method -> "Bronstein"];
Print["cM-3-152 cpu ", TimeUsed[] - t0, " wall ", AbsoluteTime[] - w0, " hash ", Hash[r], " leaves ", LeafCount[r]];
