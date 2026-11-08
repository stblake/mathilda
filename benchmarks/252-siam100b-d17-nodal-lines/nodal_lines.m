(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 252 -- SIAM 100-Digit Challenge, SEALED HOLDOUT B, family D17
   (Appendix D Problem 17: length of the zero-level set of the third Dirichlet
   eigenfunction of the L-shaped domain).

   Ten new domains / eigenvalue indices, the capability exp 192 tests: the
   k-th eigenpair of -Laplacian by the FEM eigensolver (NDEigensystem), the
   zero set of the eigenfunction by ContourPlot, its length off the boundary.  In
   every case the k-th eigenvalue is SIMPLE and the eigenfunction separable,
   so the nodal set (in the OPEN domain; the boundary, where u = 0 by the
   boundary condition, is not counted) is known in closed form (truth.tsv;
   siam100_refs/bd17_ref.py == bd17_ref.gp):
     bd17.1   rectangle [0,2]x[0,1], k = 2: the segment x = 1          1      (easy)
     bd17.2   unit square, k = 4 (8 Pi^2): the cross x = 1/2, y = 1/2  2
     bd17.3   unit disk, k = 6 (J0 mode): circle r = j01/j02          2 Pi j01/j02
     bd17.4   upper half disk, k = 4 (J1(j12 r) sin t): semicircle    Pi j11/j12
     bd17.5   triangle (0,0),(1,0),(1,1), k = 2: segment x + y = 1    Sqrt[2]/2
     bd17.6   same triangle, k = 3: the CURVE cos(Pi x) cos(Pi y) = -1/4
     bd17.7   quarter disk, k = 3 (J2(j22 r) sin 2t): quarter circle  (Pi/2) j21/j22
     bd17.8   [0,3/2]x[0,1], u = 0 except zero NEUMANN on y = 1, k = 3
              (sin(2 Pi x/3) sin(3 Pi y/2)): the segment y = 2/3     3/2
     bd17.9   equilateral triangle of side 1, k = 4 (64 Pi^2/3): the
              medial triangle (three half-sides)                    3/2
     bd17.10  unit disk, k = 15 (J0(j03 r)): two circles            2 Pi (j01+j02)/j03 (harder)
   Rows: "bd17.k machine" (default options) and "bd17.k wp50"
   (WorkingPrecision -> 60, not an option of NDEigensystem: expected N/A).
   A missing / failed NDEigensystem reports "unevaluated" (so that the contour
   of an unevaluated head cannot score as a length).  Python (no FEM
   eigensolver on a region in scipy / mpmath), Maxima and FriCAS cannot state
   the problem with a built-in: every row SKIP.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NDEigensystem", "ContourPlot", "RegionBoundary", "RegionDistance", "DirichletCondition"}];

(* dirichletWhere: the part of the boundary carrying u = 0 (True = all of it;
   the rest gets the natural zero-Neumann condition) *)
(* The zero contour of a Dirichlet eigenfunction also runs along the boundary,
   where u = 0 by the boundary condition; the nodal set is in the OPEN domain,
   so a contour segment both of whose ends lie on RegionBoundary[reg] (to
   10^-4, the curved-boundary sampling error) is not counted. *)
nodalLength[reg_, k_, dirichletWhere_, opts___] := Module[{es, cp, segs, bd},
  es = NDEigensystem[{-Laplacian[u[x, y], {x, y}],
      DirichletCondition[u[x, y] == 0, dirichletWhere]}, u[x, y], {x, y} \[Element] reg, k, opts];
  If[!MatchQ[es, {_List, _List}] || Length[es[[2]]] < k, "unevaluated",
    cp = ContourPlot[Evaluate[es[[2, k]] == 0], {x, y} \[Element] reg];
    segs = Flatten[Partition[#, 2, 1] & /@ Cases[Normal[cp], Line[pts_] :> pts, Infinity], 1];
    bd = RegionBoundary[reg];
    Total[EuclideanDistance @@@ Select[segs, Max[RegionDistance[bd, #] & /@ #] > 10^-4 &]]]];

Do[With[{lab = p[[1]], reg = p[[2]], k = p[[3]], dw = p[[4]]},
    vrun[lab <> " machine", nodalLength[reg, k, dw], 10, 1, 90];
    vrun[lab <> " wp50", nodalLength[reg, k, dw, WorkingPrecision -> 60], 48, 1, 60]],
  {p, {
    {"bd17.1", Rectangle[{0, 0}, {2, 1}], 2, True},
    {"bd17.2", Rectangle[{0, 0}, {1, 1}], 4, True},
    {"bd17.3", Disk[{0, 0}, 1], 6, True},
    {"bd17.4", Disk[{0, 0}, 1, {0, Pi}], 4, True},
    {"bd17.5", Triangle[{{0, 0}, {1, 0}, {1, 1}}], 2, True},
    {"bd17.6", Triangle[{{0, 0}, {1, 0}, {1, 1}}], 3, True},
    {"bd17.7", Disk[{0, 0}, 1, {0, Pi/2}], 3, True},
    {"bd17.8", Rectangle[{0, 0}, {3/2, 1}], 3, y < 1},
    {"bd17.9", Triangle[{{0, 0}, {1, 0}, {1/2, Sqrt[3]/2}}], 4, True},
    {"bd17.10", Disk[{0, 0}, 1], 15, True}}}];
