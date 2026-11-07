(* bench-timeout: 300 *)
(* ==========================================================================
   Experiment 193 -- SIAM 100-Digit Challenge, Appendix D, Problem 18
   ("Koch snowflake eigenvalue").

     The smallest eigenvalue of -Laplacian with zero boundary conditions on
     the Koch snowflake built on an equilateral triangle of side 1.

   No closed form; the best published value is lambda_1 = 39.3485529 (Banjai
   2007: conformal maps of the level-j polygons + extrapolation in j; 13.1161843
   for the snowflake inscribed in the unit circle, times 3 for side 1, "doubt
   over the last digit"), inside the rigorous enclosure [39.34803, 39.34869]
   of Banjai & Boulton (2019).  Only 7 digits are established: truth.tsv keeps
   39.34855 and asks for 5 (known - 2).

   Idiomatic statement: NDEigenvalues (FEM) on the level-k polygonal
   approximation of the snowflake (3*4^k sides), built by the Koch rule:
     d18.machine          level k = 5 (3072 sides), scored against the
                          FRACTAL's lambda_1 (the question asked)
     d18.machine level=0  level 0 = the equilateral triangle of side 1, whose
                          lambda_1 = 16 Pi^2/3 exactly: a reduced instance
                          with a full-precision truth
     d18.wp50             level 5 with WorkingPrecision -> 60 (not an option of
                          NDEigenvalues: expected N/A; the 50-digit tier is
                          infeasible for every system)
   A missing NDEigenvalues reports "unevaluated".  Python (no FEM eigensolver
   on a region in scipy/mpmath), Maxima and FriCAS cannot state the problem
   with a built-in: all rows SKIP.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NDEigenvalues", "RotationMatrix", "Polygon"}];

(* one Koch refinement of a closed polygon (vertices in order, clockwise so
   that the bumps point outward) *)
kochStep[pts_] := Join @@ (Function[{a, b}, With[{d = (b - a)/3},
      {a, a + d, a + d + RotationMatrix[Pi/3].d, a + 2 d}]] @@@ Partition[pts, 2, 1, 1]);
koch[k_] := Polygon[Nest[kochStep, N[{{0, 0}, {1/2, Sqrt[3]/2}, {1, 0}}], k]];

lambda1[k_, opts___] := Module[{ev},
  ev = NDEigenvalues[{-Laplacian[u[x, y], {x, y}], DirichletCondition[u[x, y] == 0, True]},
     u[x, y], {x, y} \[Element] koch[k], 1, opts];
  If[ListQ[ev], First[ev], "unevaluated"]];

vrun["d18.machine", lambda1[5], 10];
vrun["d18.machine level=0", lambda1[0], 10, 3];
vrun["d18.wp50", lambda1[5, WorkingPrecision -> 60], 48];
