# Section 8.3: the same polynomial body, run two ways over a million points.
# poly is compiled; treePoly is an ordinary definition the evaluator rewrites
# node by node. Timing distinct expressions (poly vs treePoly) sidesteps any
# result caching. The printed ratio is machine-dependent; the order of it is not.
poly = Compile[{{x, _Real}}, (((x - 3.) x + 2.) x - 1.) x + 0.5];
treePoly[x_] := (((x - 3.) x + 2.) x - 1.) x + 0.5;
data = Range[1., 1000000.];
tCompiled = First[AbsoluteTiming[Map[poly, data]]];
tInterpreted = First[AbsoluteTiming[Map[treePoly, data]]];
Round[tInterpreted / tCompiled]
