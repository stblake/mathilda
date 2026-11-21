t0 = TimeUsed[]; w0 = AbsoluteTime[];
f = (2*E^x^2*x)/(6*(1 + E^x^2)^(2/3) + 3*(1 + E^x^2)); r = Integrate[f, x];
Print["expcube-ish cpu ", TimeUsed[] - t0, " wall ", AbsoluteTime[] - w0, " hash ", Hash[r], " leaves ", LeafCount[r]];
chk = If[FreeQ[r, Integrate], Max[Abs[N[(D[r, x] - f) /. x -> #, 30] & /@ {3/7, 5/4, 7/3}]], "uneval"]; Print["check ", chk];
