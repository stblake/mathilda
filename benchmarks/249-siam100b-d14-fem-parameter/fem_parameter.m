(* bench-timeout: 1800 *)
(* ==========================================================================
   Experiment 249 -- SIAM 100-Digit Challenge, sealed holdout B for Appendix D
   Problem 14 (pentagon Poisson problem with a Neumann side: find alpha with
   Int_{Neumann side} u = e^alpha).

   Capability probed: a finite-element boundary-value solve (NDSolveValue on a
   region, DirichletCondition / natural Neumann condition) inside a scalar
   root finder (FindRoot, Brent) on a functional of the solution (a point
   value, a boundary integral or a domain integral).  Every problem has a
   closed-form or rapidly convergent series solution, so its truth is known to
   >= 55 digits even though a FEM delivers far fewer:

     bd14.1  unit disk, -Lap u = exp(a r^2), u = 0: u(0) = 1/2.  a in [1, 4]
                                                                  a = 2.21977...
     bd14.2  rectangle [0,L]x[0,1], -Lap u = 1, u = 0 on x=0, x=L, y=0, du/dn = 0
             on y = 1: Int_0^L u(x,1) dx = 1/2.  L in [1.5, 4]    L = 2.22172...
     bd14.3  equilateral triangle of side s, -Lap u = 1, u = 0: Int u = 1/100.
             s in [0.5, 2]                                        s = 1.16586...
     bd14.4  annulus 1/2 < r < 1, -Lap u = exp(a r^2), u = 0 on r = 1, du/dn = 0
             on r = 1/2: Int_{r=1/2} u ds = e^a (the pentagon question on an
             annulus).  a in [-4, -1]                            a = -2.07529...
     bd14.5  unit disk Bratu -Lap u = lam e^u, u = 0 (NONLINEAR FEM): u(0) = 1
             on the lower branch.  lam in [1, 1.99]             lam = 1.90920...
     bd14.6  unit disk -Lap u + k^2 u = 1, u = 0: Int u = Pi/10.  k in [0.5, 3]
                                                                  k = 1.22936...
     bd14.7  rectangle [0,a]x[0,1], -Lap u = 1, u = 0: u(a/2, 1/2) = 1/10.
             a in [1, 3]                                          a = 1.47965...
     bd14.8  unit BALL (3-D FEM), -Lap u = exp(a r^2), u = 0: u(0) = 1/4.
             a in [0.5, 3]                                        a = 1.21355...
     bd14.9  unit square, -div((1 + a x) grad u) = 1, u = 0 on x = 0, 1, Neumann
             on y = 0, 1 (variable coefficient): Int u = 1/20.  a in [0.5, 5]
                                                                  a = 1.51708...
     bd14.10 unit BALL Bratu -Lap u = lam e^u, u = 0 (3-D NONLINEAR FEM):
             u(0) = 1.  lam in [1, 3.3] (fold at 3.32)          lam = 3.04599...

   machine : quadratic elements, MaxCellMeasure 10^-4 (2-D) / 10^-3 (3-D);
             Dirichlet predicates are numeric (see exp. 189 for the pitfall).
   wp50    : the same with WorkingPrecision -> 60 (NDSolveValue, NIntegrate,
             FindRoot).  No system has an arbitrary-precision finite-element
             solver, so this tier is expected N/A everywhere; it is written so
             that its cost / failure is recorded.
   True values and provenance: truth.tsv; references siam100_refs/bd14_ref.py.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NDSolveValue", "DirichletCondition", "NIntegrate", "FindRoot"}];

fem[2] = Method -> {"FiniteElement", "MeshOptions" -> {"MaxCellMeasure" -> 10.^-4}};
fem[3] = Method -> {"FiniteElement", "MeshOptions" -> {"MaxCellMeasure" -> 10.^-3}};
wpo[MachinePrecision] = {};
wpo[60] = {WorkingPrecision -> 60};
br[{a_, b_}, MachinePrecision] := N[{a, b}];
br[{a_, b_}, 60] := {a, b};
tri[s_] := Triangle[{{0, 0}, {s, 0}, {s/2, s Sqrt[3]/2}}];
sol2[pde_, dc_, reg_, wp_] := NDSolveValue[{pde, DirichletCondition[u[x, y] == 0, dc]}, u,
    {x, y} \[Element] reg, fem[2], Sequence @@ wpo[wp]];
sol3[pde_, reg_, wp_] := NDSolveValue[{pde, DirichletCondition[u[x, y, z] == 0, True]}, u,
    {x, y, z} \[Element] reg, fem[3], Sequence @@ wpo[wp]];
lap2 = Laplacian[u[x, y], {x, y}];
lap3 = Laplacian[u[x, y, z], {x, y, z}];

f[1][a_?NumericQ, wp_] := sol2[-lap2 == Exp[a (x^2 + y^2)], True, Disk[], wp][0, 0];
f[2][L_?NumericQ, wp_] := With[{s = sol2[-lap2 == 1, x < 10.^-9 || x > L - 10.^-9 || y < 10.^-9,
      Rectangle[{0, 0}, {L, 1}], wp]}, NIntegrate[s[x, 1], {x, 0, L}, Evaluate[Sequence @@ wpo[wp]]]];
f[3][sd_?NumericQ, wp_] := With[{s = sol2[-lap2 == 1, True, tri[sd], wp]},
    NIntegrate[s[x, y], {x, y} \[Element] tri[sd], Evaluate[Sequence @@ wpo[wp]]]];
f[4][a_?NumericQ, wp_] := With[{s = sol2[-lap2 == Exp[a (x^2 + y^2)], x^2 + y^2 > 0.5625,
      Annulus[{0, 0}, {1/2, 1}], wp]},
    NIntegrate[s[Cos[th]/2, Sin[th]/2]/2, {th, 0, 2 Pi}, Evaluate[Sequence @@ wpo[wp]]]];
f[5][lam_?NumericQ, wp_] := sol2[-lap2 == lam Exp[u[x, y]], True, Disk[], wp][0, 0];
f[6][k_?NumericQ, wp_] := With[{s = sol2[-lap2 + k^2 u[x, y] == 1, True, Disk[], wp]},
    NIntegrate[s[x, y], {x, y} \[Element] Disk[], Evaluate[Sequence @@ wpo[wp]]]];
f[7][a_?NumericQ, wp_] := sol2[-lap2 == 1, True, Rectangle[{0, 0}, {a, 1}], wp][a/2, 1/2];
f[8][a_?NumericQ, wp_] := sol3[-lap3 == Exp[a (x^2 + y^2 + z^2)], Ball[], wp][0, 0, 0];
f[9][a_?NumericQ, wp_] := With[{s = sol2[-Inactive[Div][(1 + a x) Inactive[Grad][u[x, y], {x, y}], {x, y}] == 1,
      x < 10.^-9 || x > 1 - 10.^-9, Rectangle[], wp]},
    NIntegrate[s[x, y], {x, y} \[Element] Rectangle[], Evaluate[Sequence @@ wpo[wp]]]];
f[10][lam_?NumericQ, wp_] := sol3[-lap3 == lam Exp[u[x, y, z]], Ball[], wp][0, 0, 0];

(* {target as a function of the parameter, bracket} *)
spec = {
  {1/2 &, {1, 4}}, {1/2 &, {3/2, 4}}, {1/100 &, {1/2, 2}}, {Exp[#] &, {-4, -1}},
  {1 &, {1, 199/100}}, {Pi/10 &, {1/2, 3}}, {1/10 &, {1, 3}}, {1/4 &, {1/2, 3}},
  {1/20 &, {1/2, 5}}, {1 &, {1, 33/10}}};

solve[k_, wp_] := Module[{p, b = br[spec[[k, 2]], wp]},
  ruleVal[FindRoot[f[k][p, wp] == spec[[k, 1]][p], {p, b[[1]], b[[2]]}, Method -> "Brent",
                   Evaluate[Sequence @@ wpo[wp]]], p]];

Do[vrun["bd14." <> ToString[k] <> " machine", solve[k, MachinePrecision], 10], {k, 1, 10}];
Do[vrun["bd14." <> ToString[k] <> " wp50", solve[k, 60], 48], {k, 1, 10}];
