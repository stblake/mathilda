t0 = TimeUsed[]; w0 = AbsoluteTime[];
f = ((2 + x)/(-3 + x))^(2/3)/x; r = Integrate[f, x];
Print["trager-ratradicand-007 cpu ", TimeUsed[] - t0, " wall ", AbsoluteTime[] - w0, " hash ", Hash[r], " leaves ", LeafCount[r]];
chk = If[FreeQ[r, Integrate], Max[Abs[N[(D[r, x] - f) /. x -> #, 30] & /@ {3/7, 5/4, 7/3}]], "uneval"]; Print["check ", chk];
