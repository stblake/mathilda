t0 = TimeUsed[]; w0 = AbsoluteTime[];
r = Integrate[(E^x^2 + 3*E^x^2*x^2 + 2*E^x^2*x^4 + E^x^2*x*(1 + E^x^2*x)^(1/3) + 2*E^x^2*x^3*(1 + E^x^2*x)^(1/3) + 3*(1 + E^x^2*x)^(2/3) + 3*x*(1 + E^x^2*x) + 3*(1 + E^x^2*x)^(4/3))/(3*x*(1 + E^x^2*x)^(2/3) + 3*(1 + E^x^2*x)), x, Method -> "Bronstein"];
Print["lccontent-2-3 cpu ", TimeUsed[] - t0, " wall ", AbsoluteTime[] - w0, " hash ", Hash[r], " leaves ", LeafCount[r]];
