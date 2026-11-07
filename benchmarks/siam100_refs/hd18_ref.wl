(* Reference B for holdout family hd18 (exp 223): the closed forms in Wolfram
   arithmetic (Pi, BesselJZero, FindRoot of the annulus cross product) at 80
   digits, plus an FEM confirmation of each eigenvalue and its RANK
   (NDEigenvalues, 2nd-order mesh, MaxCellMeasure 2e-4).  Reference A:
   siam100_refs/hd18_ref.py (mpmath).  Run: wolframscript -file hd18_ref.wl *)
ann = k /. FindRoot[BesselJ[0, k/2] BesselY[0, k] - BesselJ[0, k] BesselY[0, k/2], {k, 62/10},
    WorkingPrecision -> 100];
exact = {64 Pi^2/3, 8 Pi^2, 5 Pi^2/4, 5 Pi^2, 112 Pi^2/9, BesselJZero[0, 1]^2, BesselJZero[1, 1]^2,
   BesselJZero[2, 1]^2, BesselJZero[3/2, 1]^2, ann^2};
regions = {Polygon[{{0, 0}, {1, 0}, {1/2, Sqrt[3]/2}}], Rectangle[{0, 0}, {1, 1}], Rectangle[{0, 0}, {2, 1}],
  Polygon[{{0, 0}, {1, 0}, {0, 1}}], Polygon[{{0, 0}, {1/2, 0}, {1/2, Sqrt[3]/2}}], Disk[{0, 0}, 1],
  Disk[{0, 0}, 1, {0, Pi}], Disk[{0, 0}, 1, {0, Pi/2}], Disk[{0, 0}, 1, {0, 2 Pi/3}], Annulus[{0, 0}, {1/2, 1}]};
index = {4, 4, 1, 1, 1, 1, 1, 1, 1, 1};
Do[Module[{ev = NDEigenvalues[{-Laplacian[u[x, y], {x, y}], DirichletCondition[u[x, y] == 0, True]},
      u[x, y], {x, y} \[Element] regions[[i]], index[[i]] + 1,
      Method -> {"PDEDiscretization" -> {"FiniteElement", "MeshOptions" -> {"MaxCellMeasure" -> 2*^-4}}}]},
   Print["hd18.", i, " exact ", N[exact[[i]], 72]];
   Print["       FEM first ", index[[i]] + 1, ": ", ev, "   rel.err of #", index[[i]], ": ",
     ScientificForm[(ev[[index[[i]]]] - exact[[i]])/exact[[i]], 2]]], {i, 10}];
