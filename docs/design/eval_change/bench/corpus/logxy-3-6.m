t0 = TimeUsed[]; w0 = AbsoluteTime[];
r = Integrate[(-2 - 4*E^x^2*x - 3*(-2 + 2*E^x^2 + 2*x)^(2/3))/(6*(-2 + 2*E^x^2 + 2*x)^(2/3) + 3*x*(-2 + 2*E^x^2 + 2*x)^(2/3) + 3*(-2 + 2*E^x^2 + 2*x)), x, Method -> "Bronstein"];
Print["logxy-3-6 cpu ", TimeUsed[] - t0, " wall ", AbsoluteTime[] - w0, " hash ", Hash[r], " leaves ", LeafCount[r]];
