t0 = TimeUsed[]; w0 = AbsoluteTime[];
f = 2/((1 + 2*x)*Log[1 + 2*x]*(1 + Log[1 + 2*x])^(1/6)); r = Integrate[f, x, Method -> "Bronstein"];
Print["cM-1-292 cpu ", TimeUsed[] - t0, " wall ", AbsoluteTime[] - w0, " hash ", Hash[r], " leaves ", LeafCount[r]];
chk = If[FreeQ[r, Integrate], Max[Abs[N[(D[r, x] - f) /. x -> #, 30] & /@ {3/7, 5/4, 7/3}]], "uneval"]; Print["check ", chk];
