t0 = TimeUsed[]; w0 = AbsoluteTime[];
r = Integrate[(-3*(1 + E^x^2)^(2/3) + 6*(1 + E^x^2) + 3*(1 + E^x^2)^(4/3) - 2*E^x^2*x + 9*(1 + E^x^2)*x + 3*(1 + E^x^2)^(4/3)*x + 4*E^x^2*x^2 + 2*E^x^2*(1 + E^x^2)^(1/3)*x^2 + 3*(1 + E^x^2)*x^2 + 2*E^x^2*x^3)/(6*(1 + E^x^2)^(2/3) + 3*(1 + E^x^2) + 3*(1 + E^x^2)^(2/3)*x), x];
Print["mperturb-2-10 cpu ", TimeUsed[] - t0, " wall ", AbsoluteTime[] - w0, " hash ", Hash[r], " leaves ", LeafCount[r]];
