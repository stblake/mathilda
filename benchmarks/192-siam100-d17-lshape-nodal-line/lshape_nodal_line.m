(* bench-timeout: 300 *)
(* ==========================================================================
   Experiment 192 -- SIAM 100-Digit Challenge, Appendix D, Problem 17
   ("L-shape nodal line").

     u is an eigenfunction for the THIRD Dirichlet eigenvalue of the
     Laplacian on the L-shaped domain made of three unit squares,
     Omega = (0,2)^2 minus [1,2]^2.  What is the length of its zero-level set?

   Answer: EXACTLY 2.  lambda_3 = 2 Pi^2 is simple, with eigenfunction
   Sin[Pi x] Sin[Pi y] (it vanishes on every edge of Omega, all of which lie
   on integer lines); its zero set inside Omega is the two unit segments
   {x = 1, 0 < y < 1} and {y = 1, 0 < x < 1} that separate the squares.
   (Reading "zero-level set" as including the boundary, where u = 0 by the
   boundary condition, would add the perimeter 8: answer 10.  The benchmark
   uses the nodal set in the open domain, the usual meaning.)  Provenance:
   truth.tsv and siam100_refs/d17_ref.wl.

   Idiomatic statement: NDEigensystem (FEM) for the first three Dirichlet
   eigenpairs on the polygon; the zero set of the third eigenfunction by
   ContourPlot[f == 0, {x, y} in Omega]; its length by ArcLength of the
   contour lines.  Rows:
     d17.machine          the nodal length (the problem), machine precision
     d17.machine-reduced  lambda_3 itself (= 2 Pi^2) by NDEigenvalues: the
                          well-posed numeric step the nodal set depends on
     d17.wp50             the nodal length with WorkingPrecision -> 60 (not
                          an option of NDEigensystem: expected N/A)
   A missing NDEigensystem/NDEigenvalues is reported as "unevaluated" (the
   guard keeps a contour of an unevaluated head from scoring as a length 0).
   Python (no FEM eigensolver on a region in scipy/mpmath), Maxima and FriCAS
   have no built-in that states the problem: all rows SKIP.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NDEigensystem", "NDEigenvalues", "ContourPlot", "ArcLength"}];

lshape = Polygon[{{0, 0}, {2, 0}, {2, 1}, {1, 1}, {1, 2}, {0, 2}}];
dirichletL = {-Laplacian[u[x, y], {x, y}], DirichletCondition[u[x, y] == 0, True]};

nodalLength[opts___] := Module[{es, cp},
  es = NDEigensystem[dirichletL, u[x, y], {x, y} \[Element] lshape, 3, opts];
  If[!MatchQ[es, {_List, _List}], "unevaluated",
    cp = ContourPlot[Evaluate[es[[2, 3]] == 0], {x, y} \[Element] lshape];
    Total[ArcLength /@ Cases[Normal[cp], _Line, Infinity]]]];

lambda3[] := Module[{ev},
  ev = NDEigenvalues[dirichletL, u[x, y], {x, y} \[Element] lshape, 3];
  If[ListQ[ev], Last[ev], "unevaluated"]];

vrun["d17.machine", nodalLength[], 10];
vrun["d17.machine-reduced", lambda3[], 10, 3];
vrun["d17.wp50", nodalLength[WorkingPrecision -> 60], 48];
