(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 219 -- SIAM 100-Digit Challenge, holdout A, family d14
   (generalises App. D Problem 14, experiment 189: "for which alpha does the
   integral of u along the Neumann side equal e^alpha?").

   The capability: a PARAMETER SEARCH WRAPPED AROUND AN ELLIPTIC BOUNDARY-
   VALUE SOLVE -- find a parameter (in the source, the operator, or the
   domain itself) such that a functional of the solution (a point value, a
   boundary integral, the total integral) takes a given value.  The
   original's pentagon has only 12 known digits, so the instances are
   relatives with closed forms or rapidly convergent series (70 digits):

     1-D two-point problems (NDSolve):
     hd14.1  -u'' = e^{a x} on [0,1], u(0) = 0, u'(1) = 0;          u(1) = e^a/4
     hd14.2  -u'' + a u = 1 on [0,1], u(0) = 0, u'(1) = 0;          u(1) = 1/4
     hd14.3  -u'' = l e^u on [0,1], u(0) = u(1) = 0 (Bratu);       u(1/2) = 1/2
     2-D finite-element problems (NDSolve, "FiniteElement"):
     hd14.4  -Lap u = e^{a r^2} on the unit disk, u = 0;           u(0,0) = 1/2
     hd14.5  -Lap u = 1 on [0,s] x [0,1], u = 0;                   Int u = 1/20  (find the side s)
     hd14.6  -Lap u = 1 on [0,1] x [0,b], u = 0 on three sides, du/dn = 0 on y = b;
                                                                   Int_0^1 u(x,b) dx = 1/20
     hd14.7  -Lap u = 1 on the equilateral triangle of side s, u = 0;  Int u = 1/1000
     hd14.8  -Lap u = 1 on the annulus rho < r < 1, u = 0 on r = 1, du/dn = 0 on r = rho;
                                       Int_{r = rho} u ds = 1/5, the root with rho in (1/2, 9/10)
     hd14.9  -Lap u + a u = 1 on [0,1]^2, u = 0;                   Int u = 1/50
     hd14.10 -Lap u = e^{a x} on [0,1]^2, u = 0 on x = 0, y = 0, y = 1, du/dn = 0 on x = 1;
                                                                   Int_0^1 u(1,y) dy = e^a/10
             (the closest analogue of the original: a Neumann side, an
              exponential source, the side integral against e^a)

   Asked of each tier:
     machine: FindRoot (Brent, the bracket shown) on the parameter; 1-D rows
              NDSolve PrecisionGoal/AccuracyGoal 12; 2-D rows the FEM with
              MaxCellMeasure -> 10^-4 as in experiment 189 (the zero Neumann
              condition is the natural one; Dirichlet predicates are numeric,
              see the pitfall in experiment 189's header), functionals by
              NIntegrate of the interpolating solution.
     wp50:    the same with WorkingPrecision -> 60 (FindRoot too).  The FEM is
              machine-precision only in every system, so rows 4-10 of this tier
              are N/A (written so their cost/failure is recorded).
   True values and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NDSolve", "DirichletCondition", "NIntegrate", "FindRoot"}];

(* the solution function of a (one-function) NDSolve problem *)
solu[args___] := u /. First[NDSolve[args]];

o1[MachinePrecision] = Sequence[PrecisionGoal -> 12, AccuracyGoal -> 12];
o1[wp_] := Sequence[WorkingPrecision -> wp, PrecisionGoal -> 50, AccuracyGoal -> 50, MaxSteps -> Infinity];
fem = Method -> {"FiniteElement", "MeshOptions" -> {"MaxCellMeasure" -> 10.^-4}};
dc[cond_] := DirichletCondition[u[x, y] == 0, cond];
lap := -Laplacian[u[x, y], {x, y}];

g[1][a_?NumericQ, wp_] := solu[{-u''[x] == Exp[a x], u[0] == 0, u'[1] == 0}, u, {x, 0, 1}, o1[wp]][1] - Exp[a]/4;
g[2][a_?NumericQ, wp_] := solu[{-u''[x] + a u[x] == 1, u[0] == 0, u'[1] == 0}, u, {x, 0, 1}, o1[wp]][1] - 1/4;
g[3][a_?NumericQ, wp_] := solu[{-u''[x] == a Exp[u[x]], u[0] == 0, u[1] == 0}, u, {x, 0, 1}, o1[wp]][1/2] - 1/2;
g[4][a_?NumericQ, wp_] := solu[{lap == Exp[a (x^2 + y^2)], dc[True]}, u, {x, y} \[Element] Disk[],
    fem, WorkingPrecision -> wp][0, 0] - 1/2;
g[5][s_?NumericQ, wp_] := Module[{f = solu[{lap == 1, dc[True]}, u, {x, y} \[Element] Rectangle[{0, 0}, {s, 1}],
      fem, WorkingPrecision -> wp]}, NIntegrate[f[x, y], {x, 0, s}, {y, 0, 1}, WorkingPrecision -> wp]] - 1/20;
g[6][b_?NumericQ, wp_] := Module[{f = solu[{lap == 1, dc[y < b - 10.^-9]}, u,
      {x, y} \[Element] Rectangle[{0, 0}, {1, b}], fem, WorkingPrecision -> wp]},
    NIntegrate[f[x, b], {x, 0, 1}, WorkingPrecision -> wp]] - 1/20;
g[7][s_?NumericQ, wp_] := Module[{f = solu[{lap == 1, dc[True]}, u,
      {x, y} \[Element] Triangle[{{0, 0}, {s, 0}, {s/2, s Sqrt[3]/2}}], fem, WorkingPrecision -> wp]},
    NIntegrate[f[x, y], {x, y} \[Element] Triangle[{{0, 0}, {s, 0}, {s/2, s Sqrt[3]/2}}], WorkingPrecision -> wp]] - 1/1000;
g[8][r_?NumericQ, wp_] := Module[{f = solu[{lap == 1, dc[x^2 + y^2 > ((1 + r)/2)^2]}, u,
      {x, y} \[Element] Annulus[{0, 0}, {r, 1}], fem, WorkingPrecision -> wp]},
    NIntegrate[r f[r Cos[p], r Sin[p]], {p, 0, 2 Pi}, WorkingPrecision -> wp]] - 1/5;
g[9][a_?NumericQ, wp_] := Module[{f = solu[{lap + a u[x, y] == 1, dc[True]}, u,
      {x, y} \[Element] Rectangle[{0, 0}, {1, 1}], fem, WorkingPrecision -> wp]},
    NIntegrate[f[x, y], {x, 0, 1}, {y, 0, 1}, WorkingPrecision -> wp]] - 1/50;
g[10][a_?NumericQ, wp_] := Module[{f = solu[{lap == Exp[a x], dc[x < 1 - 10.^-9]}, u,
      {x, y} \[Element] Rectangle[{0, 0}, {1, 1}], fem, WorkingPrecision -> wp]},
    NIntegrate[f[1, y], {y, 0, 1}, WorkingPrecision -> wp] - Exp[a]/10];

brk = {{1/2, 5}, {1, 6}, {1, 33/10}, {1, 4}, {1, 2}, {1/10, 1}, {1/2, 1}, {1/2, 9/10}, {1, 30}, {-3, 0}};

Do[With[{k = k, lo = brk[[k, 1]], hi = brk[[k, 2]]},
  vrun["hd14." <> ToString[k] <> " machine",
       ruleVal[FindRoot[g[k][p, MachinePrecision] == 0, {p, N[lo], N[hi]}, Method -> "Brent"], p], 10, 1, 60];
  vrun["hd14." <> ToString[k] <> " wp50",
       ruleVal[FindRoot[g[k][p, 60] == 0, {p, lo, hi}, Method -> "Brent", WorkingPrecision -> 60,
                        AccuracyGoal -> 50, PrecisionGoal -> 50], p], 48, 1, 90]], {k, 10}];
