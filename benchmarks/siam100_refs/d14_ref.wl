(* Reference check for SIAM100 App. D problem 14 (experiment 189), method C:
   WL finite elements at two mesh sizes + Richardson (rate h^(5/3), set by the
   r^(5/6) singularity at the two Dirichlet-Neumann corners), evaluated at the
   least-squares root a0 = -2.46283475604 (d14_ref.py).  Prints
   I(a0) - e^a0 for each mesh and the extrapolation, and the implied
   alpha error (the slope d/dalpha (I - e^alpha) is taken by a difference). *)
pent = Polygon[Table[{Cos[Pi/2 + 2 Pi k/5], Sin[Pi/2 + 2 Pi k/5]}, {k, 0, 4}]];
ii[a_?NumericQ, m_] := Module[{u, x, y, s},
  s = NDSolveValue[{-Laplacian[u[x, y], {x, y}] == Exp[a (x^2 + y^2)],
      DirichletCondition[u[x, y] == 0, y > N[-Cos[Pi/5]] + 10.^-9]}, u, {x, y} \[Element] pent,
      Method -> {"FiniteElement", "MeshOptions" -> {"MaxCellMeasure" -> m}}];
  NIntegrate[s[x, -Cos[Pi/5]], {x, -Sin[Pi/5], Sin[Pi/5]}, PrecisionGoal -> 12]];
as = -2.46283475604;
i1 = ii[as, 2*^-5]; i2 = ii[as, 5*^-6];
ir = i2 + (i2 - i1)/(2^(5/3) - 1);
slope = (ii[as + 0.01, 2*^-5] - Exp[as + 0.01] - (i1 - Exp[as]))/0.01;
Print[{i1 - Exp[as], i2 - Exp[as], ir - Exp[as]}];
Print["implied alpha (Richardson): ", InputForm[as - (ir - Exp[as])/slope]];
