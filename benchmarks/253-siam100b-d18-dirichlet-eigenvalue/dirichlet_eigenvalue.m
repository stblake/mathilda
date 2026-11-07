(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 253 -- SIAM 100-Digit Challenge, SEALED HOLDOUT B, family D18
   (Appendix D Problem 18: the smallest Dirichlet eigenvalue of -Laplacian on
   the Koch snowflake).

   Ten new domains, the capability exp 193 tests -- the smallest eigenvalue
   of -Laplacian with u = 0 on the boundary, by the idiomatic FEM eigensolver
   NDEigenvalues on a geometric region -- chosen so that lambda_1 has a
   closed form (truth.tsv; siam100_refs/bd18_ref.py == bd18_ref.gp):
     bd18.1   unit square Rectangle[]                         2 Pi^2   (easy)
     bd18.2   right isosceles triangle, legs 1                5 Pi^2
     bd18.3   30-60-90 triangle (half of the equilateral
              triangle of side 1)                             112 Pi^2/9
     bd18.4   unit disk                                       BesselJZero[0,1]^2
     bd18.5   circular sector of angle 3Pi/2 (re-entrant
              corner: eigenfunction ~ r^(2/3))                BesselJZero[2/3,1]^2
     bd18.6   annulus 1/2 < r < 1                             (cross-product root)^2
     bd18.7   unit cube (3-D)                                 3 Pi^2
     bd18.8   unit ball (3-D, curved)                         Pi^2
     bd18.9   cylinder of radius 1 and height 1 (3-D)         BesselJZero[0,1]^2 + Pi^2
     bd18.10  circular sector of angle 7Pi/4 (sharper
              re-entrant corner, r^(4/7))                      BesselJZero[4/7,1]^2  (harder)
   Rows: "bd18.k machine" (default NDEigenvalues) and "bd18.k wp50"
   (WorkingPrecision -> 60, not an option of the FEM eigensolver: expected
   N/A, written so its behaviour is recorded).  A missing or failed
   NDEigenvalues reports "unevaluated".  Python (no FEM eigensolver on a
   region in scipy / mpmath), Maxima and FriCAS cannot state the problem
   with a built-in: every row SKIP.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NDEigenvalues", "DirichletCondition", "Laplacian"}];

lambda1[reg_, opts___] := Module[{vars, ev},
  vars = If[RegionEmbeddingDimension[reg] === 3, {x, y, z}, {x, y}];
  ev = NDEigenvalues[{-Laplacian[u @@ vars, vars], DirichletCondition[(u @@ vars) == 0, True]},
     u @@ vars, vars \[Element] reg, 1, opts];
  If[ListQ[ev] && Length[ev] >= 1, First[ev], "unevaluated"]];

Do[With[{lab = p[[1]], reg = p[[2]]},
    vrun[lab <> " machine", lambda1[reg], 10, 1, 90];
    vrun[lab <> " wp50", lambda1[reg, WorkingPrecision -> 60], 48, 1, 60]],
  {p, {
    {"bd18.1", Rectangle[{0, 0}, {1, 1}]},
    {"bd18.2", Triangle[{{0, 0}, {1, 0}, {0, 1}}]},
    {"bd18.3", Triangle[{{0, 0}, {1/2, 0}, {1/2, Sqrt[3]/2}}]},
    {"bd18.4", Disk[{0, 0}, 1]},
    {"bd18.5", Disk[{0, 0}, 1, {0, 3 Pi/2}]},
    {"bd18.6", Annulus[{0, 0}, {1/2, 1}]},
    {"bd18.7", Cuboid[{0, 0, 0}, {1, 1, 1}]},
    {"bd18.8", Ball[{0, 0, 0}, 1]},
    {"bd18.9", Cylinder[{{0, 0, 0}, {0, 0, 1}}, 1]},
    {"bd18.10", Disk[{0, 0}, 1, {0, 7 Pi/4}]}}}];
