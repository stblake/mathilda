(* bench-timeout: 600 *)
(* ==========================================================================
   Experiment 223 -- SIAM 100-Digit Challenge HOLDOUT A, family hd18
   (generalises exp 193, App. D Problem 18, "Koch snowflake eigenvalue").

   Exp 193 asks for the smallest Dirichlet eigenvalue of -Laplacian on a
   polygonal (fractal-limit) domain, stated as NDEigenvalues (FEM) on the
   region.  This family asks the same capability -- a Dirichlet eigenvalue of
   -Laplacian on a planar region, by NDEigenvalues -- on regions whose
   spectrum is KNOWN (separable or Lame domains), varying the shape (polygon,
   curved, sector, non-simply-connected) and the eigenvalue index:

   k  region                                   index  lambda (truth.tsv)
   1  equilateral triangle, side 1             4      64 Pi^2/3
   2  unit square                              4      8 Pi^2
   3  rectangle 2 x 1                          1      5 Pi^2/4
   4  right isosceles triangle, legs 1         1      5 Pi^2
   5  30-60-90 triangle (half the side-1       1      112 Pi^2/9
      equilateral, cut along an altitude)
   6  unit disk                                1      j_{0,1}^2
   7  half disk                                1      j_{1,1}^2
   8  quarter disk                             1      j_{2,1}^2
   9  circular sector of angle 2Pi/3           1      j_{3/2,1}^2
   10 annulus 1/2 < r < 1                      1      root of a Bessel
                                                      cross product (squared)
   Rows hd18.k machine (default NDEigenvalues) and hd18.k wp50 (the same with
   WorkingPrecision -> 60, not an option of the FEM eigensolver: expected
   N/A).  A missing NDEigenvalues reports "unevaluated".  Python (no FEM
   eigensolver on a region in scipy/mpmath), Maxima and FriCAS have no
   built-in that states the problem: all rows SKIP (as in exp 193).
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NDEigenvalues", "Polygon", "Disk", "Annulus"}];

regions = {
  Polygon[{{0, 0}, {1, 0}, {1/2, Sqrt[3]/2}}],
  Rectangle[{0, 0}, {1, 1}],
  Rectangle[{0, 0}, {2, 1}],
  Polygon[{{0, 0}, {1, 0}, {0, 1}}],
  Polygon[{{0, 0}, {1/2, 0}, {1/2, Sqrt[3]/2}}],
  Disk[{0, 0}, 1],
  Disk[{0, 0}, 1, {0, Pi}],
  Disk[{0, 0}, 1, {0, Pi/2}],
  Disk[{0, 0}, 1, {0, 2 Pi/3}],
  Annulus[{0, 0}, {1/2, 1}]};
index = {4, 4, 1, 1, 1, 1, 1, 1, 1, 1};

lambdaK[reg_, k_, opts___] := Module[{ev},
  ev = NDEigenvalues[{-Laplacian[u[x, y], {x, y}], DirichletCondition[u[x, y] == 0, True]},
     u[x, y], {x, y} \[Element] reg, k, opts];
  If[ListQ[ev] && Length[ev] >= k, ev[[k]], "unevaluated"]];

Do[With[{reg = regions[[i]], k = index[[i]], lab = "hd18." <> ToString[i]},
   vrun[lab <> " machine", lambdaK[reg, k], 10, 1, 60];
   vrun[lab <> " wp50", lambdaK[reg, k, WorkingPrecision -> 60], 48, 1, 60]],
  {i, Length[regions]}];
