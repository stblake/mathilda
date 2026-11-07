(* Reference for SIAM100 App. D problem 17 (exp 192): length of the zero set
   of the third Dirichlet eigenfunction on the L-shaped domain of three unit
   squares,  Omega = (0,2)^2 \ [1,2]^2.

   Method 1 (exact): u = Sin[Pi x] Sin[Pi y] satisfies -Lap u = 2 Pi^2 u on
   Omega and vanishes on all of dOmega (every boundary edge lies on a line
   x in Z or y in Z), so 2 Pi^2 is an eigenvalue; checked symbolically below.
   Its interior zero set is {x = 1, 0 < y < 1} U {y = 1, 0 < x < 1}: length 2.
   Method 2 (numeric): fine-mesh FEM eigenvalues confirm 2 Pi^2 is the THIRD
   eigenvalue and SIMPLE (published L-shape values 9.6397238, 15.1972519,
   19.7392088 = 2 Pi^2, 29.5214811: Fox-Henrici-Moler 1967, Betcke-Trefethen
   2005), so the eigenfunction -- and its zero set -- is unique up to scale;
   and the FEM u_3 restricted to the two interior segments is ~0 relative to
   its maximum.  Run: wolframscript -file d17_ref.wl *)
u = Sin[Pi x] Sin[Pi y];
Print["-Lap u - 2Pi^2 u == 0: ", Simplify[-Laplacian[u, {x, y}] - 2 Pi^2 u] === 0];
Print["u on boundary edges: ", Simplify[{u /. y -> 0, u /. x -> 0, u /. x -> 2, u /. y -> 2,
   u /. x -> 1, u /. y -> 1}]];
Needs["NDSolve`FEM`"];
reg = Polygon[{{0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}}];
mesh = ToElementMesh[reg, MaxCellMeasure -> 0.0005, "MeshOrder" -> 2];
{vals, funs} = NDEigensystem[{-Laplacian[w[x, y], {x, y}], DirichletCondition[w[x, y] == 0, True]},
   w[x, y], {x, y} \[Element] mesh, 5];
Print["FEM eigenvalues: ", NumberForm[vals, 10], "   2Pi^2 = ", N[2 Pi^2, 12]];
f3 = funs[[3]];
mx = Max[Abs[Table[f3, {x, 0.05, 1.95, 0.05}, {y, 0.05, 0.95, 0.05}]]];
seg = Max[Abs[Join[Table[f3 /. {x -> 1, y -> s}, {s, 0.02, 0.98, 0.02}],
     Table[f3 /. {x -> s, y -> 1}, {s, 0.02, 0.98, 0.02}]]]];
Print["max |u3| on the two interior segments / max |u3| = ", seg/mx];
Print["nodal length (exact) = 2"];
