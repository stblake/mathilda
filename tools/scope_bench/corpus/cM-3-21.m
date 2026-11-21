t0 = TimeUsed[]; w0 = AbsoluteTime[];
f = (1 + E^x)^(7/6)/(-2 + E^x); r = Integrate[f, x];
Print["cM-3-21 cpu ", TimeUsed[] - t0, " wall ", AbsoluteTime[] - w0, " hash ", Hash[r], " leaves ", LeafCount[r]];
chk = If[FreeQ[r, Integrate], Max[Abs[N[(D[r, x] - f) /. x -> #, 30] & /@ {3/7, 5/4, 7/3}]], "uneval"]; Print["check ", chk];
