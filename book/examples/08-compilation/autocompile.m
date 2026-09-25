# Section 8.5: auto-compilation. On by default; the switch changes speed, not answers.
$AutoCompilation
$AutoCompilation = False;
r1 = Table[Sin[x] Exp[-x/10.] + Sqrt[x + 1.], {x, 0., 5., 0.5}];
$AutoCompilation = True;
r2 = Table[Sin[x] Exp[-x/10.] + Sqrt[x + 1.], {x, 0., 5., 0.5}];
r1 == r2
$AutoCompilation = False;
tOff = First[AbsoluteTiming[Nest[3.5 # (1. - #) &, 0.31, 1000000]]];
$AutoCompilation = True;
tOn = First[AbsoluteTiming[Nest[3.5 # (1. - #) &, 0.32, 1000000]]];
Round[tOff / tOn]
