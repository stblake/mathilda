(* bench-timeout: 480 *)
(* ==========================================================================
   Experiment 187 -- SIAM 100-Digit Challenge, Appendix D, Problem 12
   (Bornemann): a stiffly constrained particle.

       T = (x'^2 + y'^2)/2,  U = y + (eps^-2/2) (1 + alpha x^2) (x^2 + y^2 - 1)^2,
       (x, y)(0) = (0, 1),  (x', y')(0) = (1, 1).
       For which alpha does the particle first hit y = 0 at time 10, in the
       limit eps -> 0 ?                          answer = 10.6234113996...

   The eps -> 0 limit is NOT the naive constraint x^2 + y^2 = 1.  The initial
   velocity has a normal component 1, so the particle oscillates across the
   circle with frequency eps^-1 omega(x), omega = 2 Sqrt[1 + alpha x^2] (the
   normal Hessian of the constraint term at r = 1).  The normal action
   E_perp/omega is an adiabatic invariant (Rubin-Ungar 1957, Takens 1980;
   Bornemann, "Homogenization in time of singularly perturbed mechanical
   systems", LNM 1687), theta* = (1/2)/omega(0) = 1/4, and the limit motion is
   motion ON the circle under the homogenised potential
       U_hom = y + theta* omega(x) = y + Sqrt[1 + alpha x^2]/2.
   With x = Cos[th], y = Sin[th] (th(0) = Pi/2, th'(0) = -1):
       th'' = -Cos[th] + alpha Sin[th] Cos[th] / (2 Sqrt[1 + alpha Cos[th]^2]).
   (Checked against the full eps-system, eps = 1e-2 .. 2.5e-3: hitting times
   agree to O(eps^2); the naive constrained pendulum is off by 3-15%.)
   The correction builds a barrier at Sin[th] = 2 Sqrt[1 + alpha Cos[th]^2]/alpha;
   the particle (energy 2) crosses it iff alpha < alpha_c = 10.62347538..., and
   the hitting time grows like -log(alpha_c - alpha), so T(alpha) = 10 has
   exactly one root, just below alpha_c.

   Asked of each tier (the same limiting ODE everywhere):
     machine: FindRoot (Brent, bracket [10, 11]) on alpha of the first hitting
              time of Sin[th] = 0, found by NDSolveValue + WhenEvent
              (integration to t = 100; no hit -> domain end 100 > 10),
              NDSolve tolerances 1e-12 (as scipy rtol/atol 1e-12).
     wp50:    the same at WorkingPrecision 60, PrecisionGoal 50, but on the
              bracket [10.6234, 10.62347] that the event-less rows (mpmath,
              Maxima) need too.  At WorkingPrecision 60 wolframscript times
              out at 300 s with either bracket (as does mpmath: one 60-digit
              odefun solve to t = 10 costs ~90 s).
   True value and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NDSolve", "WhenEvent", "FindRoot"}];

hit[a_?NumericQ] := (th /. First[NDSolve[{
    th''[t] == -Cos[th[t]] + a Sin[th[t]] Cos[th[t]]/(2 Sqrt[1 + a Cos[th[t]]^2]),
    th[0] == Pi/2, th'[0] == -1,
    WhenEvent[Sin[th[t]] == 0, "StopIntegration"]}, th, {t, 0, 100},
    PrecisionGoal -> 12, AccuracyGoal -> 12]])["Domain"][[1, 2]];

hit60[a_?NumericQ] := (th /. First[NDSolve[{
    th''[t] == -Cos[th[t]] + a Sin[th[t]] Cos[th[t]]/(2 Sqrt[1 + a Cos[th[t]]^2]),
    th[0] == Pi/2, th'[0] == -1,
    WhenEvent[Sin[th[t]] == 0, "StopIntegration"]}, th, {t, 0, 100},
    WorkingPrecision -> 60, PrecisionGoal -> 50, AccuracyGoal -> 50,
    MaxSteps -> Infinity]])["Domain"][[1, 2]];

vrun["d12.machine",
     ruleVal[FindRoot[hit[a] == 10, {a, 10, 11}, Method -> "Brent"], a], 10];
vrun["d12.wp50",
     ruleVal[FindRoot[hit60[a] == 10, {a, 106234/10000, 1062347/100000}, Method -> "Brent",
                      WorkingPrecision -> 60, AccuracyGoal -> 50,
                      PrecisionGoal -> 50], a], 48, 1, 300];
