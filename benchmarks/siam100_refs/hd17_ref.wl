(* Reference B for holdout family hd17 (exp 222): numeric confirmation of the
   closed-form nodal sets.  For each problem: NDEigensystem on a fine 2nd-order
   mesh (MaxCellMeasure 2e-4) -- the eigenvalue at the stated rank against its
   closed form, the max of |u_k| on the stated nodal set relative to max |u_k|,
   and the ContourPlot nodal length.  (The exact values: siam100_refs/hd18_ref.py.)
   Run: wolframscript -file hd17_ref.wl *)
j[n_, k_] := BesselJZero[n, k];
regions = {Rectangle[{0, 0}, {1, 1}], Rectangle[{0, 0}, {2, 1}], Polygon[{{0, 0}, {1, 0}, {1/2, Sqrt[3]/2}}],
  Disk[{0, 0}, 1], Disk[{0, 0}, 1, {0, Pi}], Disk[{0, 0}, 1, {0, Pi}], Disk[{0, 0}, 1], Rectangle[{0, 0}, {1, 1}],
  Disk[{0, 0}, 1], Polygon[{{0, 0}, {1, 0}, {0, 1}}]};
index = {4, 2, 4, 6, 4, 2, 2, 11, 15, 2};
lam = {8 Pi^2, 2 Pi^2, 64 Pi^2/3, j[0, 2]^2, j[1, 2]^2, j[2, 1]^2, j[1, 1]^2, 18 Pi^2, j[0, 3]^2, 10 Pi^2};
len = {2, 1, 3/2, 2 Pi j[0, 1]/j[0, 2], Pi j[1, 1]/j[1, 2], 1, 2, 4, 2 Pi (j[0, 1] + j[0, 2])/j[0, 3], Sqrt[2]/2};
(* sample points on the claimed nodal sets *)
pts = {Join[Table[{1/2, s}, {s, .05, .95, .1}], Table[{s, 1/2}, {s, .05, .95, .1}]],
  Table[{1, s}, {s, .05, .95, .1}],
  Join[Table[(1 - s) {1/4, Sqrt[3]/4} + s {3/4, Sqrt[3]/4}, {s, .05, .95, .1}],
    Table[(1 - s) {1/4, Sqrt[3]/4} + s {1/2, 0}, {s, .05, .95, .1}],
    Table[(1 - s) {1/2, 0} + s {3/4, Sqrt[3]/4}, {s, .05, .95, .1}]],
  Table[N[j[0, 1]/j[0, 2]] {Cos[a], Sin[a]}, {a, 0, 2 Pi, .5}],
  Table[N[j[1, 1]/j[1, 2]] {Cos[a], Sin[a]}, {a, .1, 3, .3}],
  Table[{0, s}, {s, .05, .95, .1}],
  {}, Join[Table[{1/3, s}, {s, .05, .95, .1}], Table[{s, 2/3}, {s, .05, .95, .1}]],
  Join[Table[N[j[0, 1]/j[0, 3]] {Cos[a], Sin[a]}, {a, 0, 6, .5}], Table[N[j[0, 2]/j[0, 3]] {Cos[a], Sin[a]}, {a, 0, 6, .5}]],
  Table[{s/2, s/2}, {s, .05, .95, .1}]};
out = OpenWrite[FileNameJoin[{DirectoryName[$InputFileName], "hd17_ref_B.out"}]];
Do[Module[{es, f, mx, cp, L},
   es = NDEigensystem[{-Laplacian[u[x, y], {x, y}], DirichletCondition[u[x, y] == 0, True]}, u[x, y],
     {x, y} \[Element] regions[[i]], index[[i]] + 1,
     Method -> {"PDEDiscretization" -> {"FiniteElement", "MeshOptions" -> {"MaxCellMeasure" -> 2*^-4}}}];
   f = es[[2, index[[i]]]];
   mx = NMaxValue[{Abs[f], {x, y} \[Element] regions[[i]]}, {x, y}];
   cp = ContourPlot[Evaluate[f == 0], {x, y} \[Element] regions[[i]], PlotPoints -> 60];
   L = Total[ArcLength /@ Cases[Normal[cp], _Line, Infinity]];
   WriteString[out, "hd17.", i, ": eigenvalues ", es[[1]], "  #", index[[i]], " rel.err ",
     ScientificForm[(es[[1, index[[i]]]] - lam[[i]])/lam[[i]], 2],
     "  max|u| on nodal set / max|u| = ", If[pts[[i]] === {}, "n/a (rotation-free diameter)",
       ScientificForm[Max[Abs[f /. Thread[{x, y} -> #]] & /@ pts[[i]]]/mx, 2]],
     "  contour length ", L, "  exact ", N[len[[i]], 20], "\n"]], {i, 10}];
Close[out];
