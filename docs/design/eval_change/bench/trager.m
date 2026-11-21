t0 = AbsoluteTime[];
t1 = AbsoluteTiming[r1 = Integrate[((1 + x^2)*(2 + x^4)^(1/6))/x, x, Method -> "Trager"];][[1]];
Print["trager ", t1, " wall-since-start ", AbsoluteTime[] - t0, " MN ", $ModuleNumber, " MemInUse ", MemoryInUse[]];
