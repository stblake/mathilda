(* bench-timeout: 600 *)
(* ==========================================================================
   Experiment 222 -- SIAM 100-Digit Challenge HOLDOUT A, family hd17
   (generalises exp 192, App. D Problem 17, "L-shape nodal line").

   Exp 192 asks for the length of the zero set (in the open domain) of the
   eigenfunction of the 3rd Dirichlet eigenvalue of -Laplacian on the L-shape,
   stated as NDEigensystem + ContourPlot[f == 0] + ArcLength of the contour
   lines.  This family asks the same, on domains whose eigenfunctions are
   known in closed form, varying the shape, the index and the kind of nodal
   set (straight segments, a medial triangle, circles, a diameter, a
   semicircle, a diagonal):

   k  region                         index  nodal set              length
   1  unit square                    4      x = 1/2, y = 1/2        2
   2  rectangle 2 x 1                2      x = 1                   1
   3  equilateral triangle, side 1   4      medial triangle         3/2
   4  unit disk                      6      circle r = j01/j02      2 Pi j01/j02
   5  half disk (y > 0)              4      semicircle r = j11/j12  Pi j11/j12
   6  half disk (y > 0)              2      ray x = 0               1
   7  unit disk                      2      a diameter (lambda_2 is
                                            double; every eigenfunction's
                                            nodal set is a diameter)  2
   8  unit square                    11     x, y = 1/3, 2/3         4
   9  unit disk                      15     circles r = j01/j03,
                                            j02/j03         2 Pi (j01+j02)/j03
   10 right isosceles triangle        2     segment of y = x        Sqrt2/2
      (legs 1 on the axes)
      from (0,0) to (1/2,1/2)
   (Ranks and simplicity: siam100_refs/hd18_ref.py; the eigenvalues and
   nodal sets are confirmed against fine-mesh FEM in siam100_refs/hd17_ref.wl.
   Note: ContourPlot[u == 0] of an FEM eigenfunction also traces boundary
   edges where the interpolant is 0, which inflates the length on polygons --
   that is part of what the machine rows measure.)  Rows hd17.k machine and
   hd17.k wp50 (WorkingPrecision -> 60, not an option of NDEigensystem:
   expected N/A).  A missing NDEigensystem reports "unevaluated" (the guard
   keeps a contour of an unevaluated head from scoring as a length).  Python,
   Maxima and FriCAS have no FEM eigensolver: all rows SKIP (as in exp 192).
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NDEigensystem", "ContourPlot", "ArcLength"}];

regions = {
  Rectangle[{0, 0}, {1, 1}],
  Rectangle[{0, 0}, {2, 1}],
  Polygon[{{0, 0}, {1, 0}, {1/2, Sqrt[3]/2}}],
  Disk[{0, 0}, 1],
  Disk[{0, 0}, 1, {0, Pi}],
  Disk[{0, 0}, 1, {0, Pi}],
  Disk[{0, 0}, 1],
  Rectangle[{0, 0}, {1, 1}],
  Disk[{0, 0}, 1],
  Polygon[{{0, 0}, {1, 0}, {0, 1}}]};
index = {4, 2, 4, 6, 4, 2, 2, 11, 15, 2};

nodalLength[reg_, k_, opts___] := Module[{es, cp},
  es = NDEigensystem[{-Laplacian[u[x, y], {x, y}], DirichletCondition[u[x, y] == 0, True]},
     u[x, y], {x, y} \[Element] reg, k, opts];
  If[!MatchQ[es, {_List, _List}] || Length[es[[2]]] < k, "unevaluated",
    cp = ContourPlot[Evaluate[es[[2, k]] == 0], {x, y} \[Element] reg];
    Total[ArcLength /@ Cases[Normal[cp], _Line, Infinity]]]];

Do[With[{reg = regions[[i]], k = index[[i]], lab = "hd17." <> ToString[i]},
   vrun[lab <> " machine", nodalLength[reg, k], 10, 1, 60];
   vrun[lab <> " wp50", nodalLength[reg, k, WorkingPrecision -> 60], 48, 1, 60]],
  {i, Length[regions]}];
