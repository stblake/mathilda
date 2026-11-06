(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 130 -- ISSAC '97 Systems Challenge, Problem 5
   ==========================================================================

       cd benchmarks/130-issac97-p5-laguerre-root
       ../../Mathilda -file laguerre_root.m

   PROBLEM.  The largest zero of the Laguerre polynomial L_1000(x),
   3943.2473948452...  Near it |L_1000| ~ 1e867, so a double cannot even
   evaluate the function there: the problem is arbitrary precision, and the
   cost is dominated by evaluating a degree-1000 polynomial by its three-term
   recurrence at 100 digits.

   METHODS (the notebook's two, plus its verification):
     m1  FindRoot from the bound 4n + 2 = 4002: Newton with the analytic
         Jacobian -L^(1)_999 (verbatim), and the secant method from
         {4002, 4003} (verbatim), both at WorkingPrecision -> 100; Newton
         again at 200 digits for the cost-vs-digits curve;
     m1  the notebook's RootCount: the argument principle as a contour
         NIntegrate of -L^(1)_999/L_1000 round circles of radius 3940 and 3950
         (999 and 1000 roots), at 100 digits;
     m2  all 1000 roots in sequence, each FindRoot started by extrapolating
         the previous two (1000 FindRoots at 100 digits; 2400 s in 1997);
   MACHINE VARIANT.  The notebook's machine FindRoot on L_1000 (which a
   system must either escalate or fail), and the same Newton on L_100, whose
   largest zero 374.98... a double resolves (|L_100| <= 1e150 there).
   JACOBIAN FORM.  The 1997 code passed a scalar `Jacobian -> g[x]`; the
   portable spelling for one equation is the 1 x 1 matrix {{g[x]}} (a system
   that requires a matrix rejects the scalar and silently returns the start
   point), so every Newton row here uses the matrix form.
   NOT IN THE NOTEBOOK.  NSolve of L_1000 at 50 digits (all roots at once).
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"LaguerreL", "FindRoot", "NIntegrate", "NSolve", "SetPrecision"}];

f[(x_)?NumberQ] := LaguerreL[1000, SetPrecision[x, 100]];
g[(x_)?NumberQ] := -LaguerreL[999, 1, SetPrecision[x, 100]];
f2[(x_)?NumberQ] := LaguerreL[1000, SetPrecision[x, 200]];
g2[(x_)?NumberQ] := -LaguerreL[999, 1, SetPrecision[x, 200]];

run["p5 m1 Newton wp100",
    FindRoot[Evaluate[f[x]], {x, 4002}, WorkingPrecision -> 100, AccuracyGoal -> 1000,
             MaxIterations -> 50, Jacobian -> {{g[x]}}], sig[ruleVal[#, x], 40] &, 3];
run["p5 m1 secant wp100",
    FindRoot[Evaluate[f[x]], {x, 4002, 4003}, WorkingPrecision -> 100, AccuracyGoal -> 1000,
             MaxIterations -> 100], sig[ruleVal[#, x], 40] &, 3];
run["p5 m1 Newton wp200",
    FindRoot[Evaluate[f2[x]], {x, 4002}, WorkingPrecision -> 200, AccuracyGoal -> 2000,
             MaxIterations -> 50, Jacobian -> {{g2[x]}}], sig[ruleVal[#, x], 80] &, 3];

run["p5 machine FindRoot L1000", FindRoot[LaguerreL[1000, x], {x, 4002}],
    sig[ruleVal[#, x], 10] &];
run["p5 machine FindRoot L100",
    FindRoot[LaguerreL[100, x], {x, 402}, Jacobian -> {{-LaguerreL[99, 1, x]}}],
    sig[ruleVal[#, x], 10] &, 3];

rootCount[rho_] := Round[Im[NIntegrate[
    (-(LaguerreL[999, 1, rho E^(I phi)]/LaguerreL[1000, rho E^(I phi)])) I rho E^(I phi),
    {phi, 0, Pi, 2 Pi}, WorkingPrecision -> 100, PrecisionGoal -> 4]]/(2 Pi)];
run["p5 m1 RootCount contour wp100", {rootCount[3940], rootCount[3950]},
    If[MatchQ[#, {_Integer, _Integer}], #, "unevaluated"] &];

allRoots[] := Module[{r},
  {r[1], r[2]} = (x /. FindRoot[f[x], {x, #1}, WorkingPrecision -> 100, AccuracyGoal -> 1000,
                                MaxIterations -> 50, Jacobian -> {{g[x]}}] &) /@ {0, 1/100};
  Do[r[i] = x /. FindRoot[f[x], {x, r[i - 1] + (r[i - 1] - r[i - 2])/2,
                                   r[i - 1] + 3 (r[i - 1] - r[i - 2])/2},
                          WorkingPrecision -> 100, AccuracyGoal -> 1000, MaxIterations -> 50],
     {i, 3, 1000}];
  r[1000]];
run["p5 m2 all 1000 roots in sequence wp100", allRoots[], sig[#, 40] &, 1, 600];

run["p5 NSolve L1000 wp50",
    Max[x /. NSolve[LaguerreL[1000, x] == 0, x, WorkingPrecision -> 50]], sig[#, 30] &];
