(* bench-timeout: 900 *)
(* ==========================================================================
   Experiment 135 -- ISSAC '97 Systems Challenge, Problem 10
   ==========================================================================

       cd benchmarks/135-issac97-p10-ode-singularity
       ../../Mathilda -file ode_singularity.m

   PROBLEM.  y'' = x^3 + y^3 + y'^3, y(0) = y'(0) = 0.  y' blows up at a
   finite r ~ 1.6444 with a square-root branch point (y'' ~ y'^3 near r, so
   y ~ c + Sqrt[r - x]): y(r) is finite.  Find r and y(r) to 14 digits.

   METHODS (the notebook's two):
     naive  NDSolve straight at it (machine): it stops at the singularity,
            and the end of the solution's domain is r to ~6 digits;
     m1     solve for y(x) up to x = 6/5, then switch to the inverse function
            x(y), which is regular at the branch point:
                -x''(y) = (x^3 + y^3) x'(y)^3 + 1,
            and take r as the maximum of x(y) (FindMinimum of -x(y));
            machine and WorkingPrecision -> 40 (verbatim settings);
     m2     polar coordinates about (0, 1): x = rho Cos[phi],
            y = rho Sin[phi] + 1, so the solution curve stays O(1); the ODE
            for rho(phi) is derived symbolically (Together/Numerator/Factor,
            itself a case), then NDSolve and FindMinimum of -rho Cos[phi];
            machine and WorkingPrecision -> 30 (verbatim settings).
   CHECKS.  r and y(r) (leading digits): 13 at 30/40 digits, 6 at machine
   precision.  At 60 digits (the root of x'(y), a separate method)
   r = 1.644376690338840720193840833 and y(r) = 0.931938765110281053354;
   y(r) ends ...0281 at the 14th digit, so 14 digits sits on a rounding edge
   that a correct 40-digit solve can miss by one unit.
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"NDSolve", "FindMinimum", "InterpolatingFunction", "Together", "Factor"}];

ode = {y'[x] == z[x], z'[x] == x^3 + y[x]^3 + z[x]^3, y[0] == 0, z[0] == 0};

run["p10 naive NDSolve machine domain end",
    First[(y /. First[NDSolve[ode, {y, z}, {x, 0, 2}]])["Domain"]][[2]], sig[#, 5] &, 3];

(* ---- m1: y(x) to 6/5, then x(y) across the branch point ---------------- *)
m1[wp_] := Module[{opts, s1, y0, z0, s2, fm},
  opts = If[wp === MachinePrecision, {MaxSteps -> 3000},
            {WorkingPrecision -> wp, PrecisionGoal -> 30, AccuracyGoal -> 20, MaxSteps -> 3000}];
  s1 = NDSolve[ode, {y, z}, {x, 0, 6/5}, StartingStepSize -> 10^-3, Evaluate[Sequence @@ opts]];
  {y0, z0} = {y[6/5], z[6/5]} /. First[s1];
  s2 = NDSolve[{-xx''[t] == (xx[t]^3 + t^3) xx'[t]^3 + 1, xx[y0] == 6/5, xx'[y0] == 1/z0},
               xx, {t, y0, 1}, Evaluate[Sequence @@ opts]];
  fm = FindMinimum[Evaluate[-xx[t] /. First[s2]], {t, 9/10, 1},
                   Evaluate[Sequence @@ If[wp === MachinePrecision, {},
                                  {WorkingPrecision -> wp, PrecisionGoal -> 30, AccuracyGoal -> 30}]]];
  If[MatchQ[fm, {_?NumericQ, {_Rule}}], {-First[fm], t /. Last[fm]}, "failed"]];

pairChk[d_][v_] := If[MatchQ[v, {_?NumericQ, _?NumericQ}], {sig[v[[1]], d], sig[v[[2]], d]}, "failed"];

run["p10 m1 inverse function machine", m1[MachinePrecision], pairChk[6], 3];
run["p10 m1 inverse function wp40", m1[40], pairChk[13]];

(* ---- m2: polar coordinates about (0, 1) --------------------------------- *)
polarOde[] := Module[{X, Y, yp, ypp},
  X = rho[p] Cos[p]; Y = rho[p] Sin[p] + 1;
  yp = D[Y, p]/D[X, p]; ypp = D[yp, p]/D[X, p];
  Factor[Numerator[Together[ypp - (X^3 + Y^3 + yp^3)]]] == 0];
run["p10 m2 derive polar ODE (Together/Factor)", polarOde[],
    If[Head[#] === Equal, Length[Expand[First[#]]], "unevaluated"] &, 3];
newOde = polarOde[];

m2[wp_] := Module[{opts, s, fm, p0},
  opts = If[wp === MachinePrecision, {MaxSteps -> 5000},
            {WorkingPrecision -> wp, PrecisionGoal -> wp, AccuracyGoal -> 20, MaxSteps -> 5000}];
  s = NDSolve[{newOde, rho[-Pi/2] == 1, rho'[-Pi/2] == 0}, rho, {p, -Pi/2, Pi/4}, Evaluate[Sequence @@ opts]];
  fm = FindMinimum[Evaluate[-rho[p] Cos[p] /. First[s]], {p, -1/10, 1/10},
                   Evaluate[Sequence @@ If[wp === MachinePrecision, {}, {WorkingPrecision -> wp, AccuracyGoal -> 80}]]];
  If[MatchQ[fm, {_?NumericQ, {_Rule}}],
     p0 = p /. Last[fm]; {-First[fm], (rho[p0] Sin[p0] + 1) /. First[s]}, "failed"]];

run["p10 m2 polar machine", m2[MachinePrecision], pairChk[6], 3];
run["p10 m2 polar wp30", m2[30], pairChk[13]];
