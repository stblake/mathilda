(* bench-timeout: 300 *)
(* ==========================================================================
   Experiment 189 -- SIAM 100-Digit Challenge, Appendix D, Problem 14
   (Bornemann): a Poisson problem on a pentagon with mixed conditions.

       -Laplacian u = exp(alpha |x|^2) on the regular pentagon inscribed in the
       unit circle; u = 0 on four sides, du/dn = 0 on the fifth.  For which
       alpha does the integral of u along the Neumann side equal e^alpha ?
                                                    answer = -2.46283475604...

   Rotation invariance makes the choice of Neumann side immaterial; here the
   pentagon has a vertex at (0, 1) and the Neumann side is the bottom one,
   y = -Cos[Pi/5].  I(alpha) e^-alpha = Int e^{alpha(|x|^2-1)} w dx with w > 0
   the adjoint solution, which is strictly decreasing, so the root is unique.
   The two Dirichlet-Neumann corners carry an r^(5/6) singularity, which caps
   a quasi-uniform FEM at ~h^(5/3) (measured: MaxCellMeasure 1e-3 .. 5e-6
   gives 3.7e-4 .. 9.8e-6 error in I(-2.5)).

   Asked of each tier: the idiomatic statement is the finite-element pipeline
   NDSolve (pentagon region, DirichletCondition on four sides, the zero
   NeumannValue on the bottom -- the natural condition, written out) +
   NIntegrate along the bottom side + FindRoot (Brent on [-3, 0]).
     machine: mesh MaxCellMeasure -> 1e-4 (the default mesh gives ~1 digit).
     wp50:    the same with WorkingPrecision -> 60.  The FEM solver is
              machine-precision only, so this tier is N/A for every system
              (no system has an arbitrary-precision PDE solver for it).
   Pitfall met while writing this: the DirichletCondition predicate must be
   numeric (y > N[-Cos[Pi/5]] + 1e-9).  The exact predicate y > -Cos[Pi/5] +
   10^-6 silently marks the Neumann side Dirichlet as well (u = 0 there,
   I = 1e-17) in wolframscript, and the loose y > -0.8 leaves 1% of the two
   adjacent sides Neumann (error 1e-3 that does not shrink with the mesh).
   Likewise RegularPolygon[5] lists its vertices out of cyclic order; the
   pentagon is built from its vertices.
   True value and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NDSolve", "DirichletCondition", "NeumannValue", "NIntegrate", "FindRoot"}];

pent = Polygon[Table[{Cos[Pi/2 + 2 Pi k/5], Sin[Pi/2 + 2 Pi k/5]}, {k, 0, 4}]];

nint[a_?NumericQ, wp_] := Module[{s},
  s = u /. First[NDSolve[{-Laplacian[u[x, y], {x, y}] ==
        Exp[a (x^2 + y^2)] + NeumannValue[0, y < N[-Cos[Pi/5]] + 10.^-9],
      DirichletCondition[u[x, y] == 0, y > N[-Cos[Pi/5]] + 10.^-9]},
      u, {x, y} \[Element] pent,
      Method -> {"FiniteElement", "MeshOptions" -> {"MaxCellMeasure" -> 10.^-4}},
      WorkingPrecision -> wp]];
  NIntegrate[s[x, -Cos[Pi/5]], {x, -Sin[Pi/5], Sin[Pi/5]}, WorkingPrecision -> wp]];

vrun["d14.machine",
     ruleVal[FindRoot[nint[a, MachinePrecision] == Exp[a], {a, -3, 0}, Method -> "Brent"], a], 10];
vrun["d14.wp50",
     ruleVal[FindRoot[nint[a, 60] == Exp[a], {a, -3, 0}, Method -> "Brent",
                      WorkingPrecision -> 60], a], 10];
