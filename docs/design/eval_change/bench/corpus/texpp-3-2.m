t0 = TimeUsed[]; w0 = AbsoluteTime[];
r = Integrate[(E^x^2 + E^(2*x^2) + 2*E^x^2*x^2 + Log[Log[x]] + 2*E^x^2*Log[Log[x]] + E^(2*x^2)*Log[Log[x]])/(1 + E^x^2)^2, x];
Print["texpp-3-2 cpu ", TimeUsed[] - t0, " wall ", AbsoluteTime[] - w0, " hash ", Hash[r], " leaves ", LeafCount[r]];
