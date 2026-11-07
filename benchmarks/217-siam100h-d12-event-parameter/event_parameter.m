(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 217 -- SIAM 100-Digit Challenge, holdout A, family d12
   (generalises App. D Problem 12, experiment 187: "for which alpha does the
   particle first hit y = 0 at time 10").

   The capability: a PARAMETER SEARCH WRAPPED AROUND AN EVENT-LOCATING ODE
   SOLVE -- find a such that the FIRST time an event fires equals a given T.
   Ten instances that vary the kind of ODE (first/second order, autonomous /
   time-dependent coefficient, conservative / damped / self-excited, singular
   force), the kind of event (a position level, a zero, a turning point
   y' = 0) and the target time:

     hd12.1  y'' = -1 - a y'^2, y(0)=0, y'(0)=2; apex y' = 0 at t = 1       (quadratic drag)
     hd12.2  th'' = -a sin th, th(0)=0, th'(0)=1; th' = 0 at t = 2          (pendulum turning point)
     hd12.3  x'' = -x - a x^3, x(0)=1, x'(0)=0; x = 0 at t = 6/5           (Duffing)
     hd12.4  x'' + a x' + x = 0, x(0)=1, x'(0)=0; x = 0 at t = 2           (damped oscillator)
     hd12.5  y' = a + y^2, y(0)=0; y = 1 at t = 1/2                        (Riccati, first order)
     hd12.6  r'' = -a/r^2, r(0)=1, r'(0)=0; r = 1/2 at t = 1               (radial Kepler infall)
     hd12.7  y'' + a t y = 0, y(0)=1, y'(0)=0; y = 0 at t = 1              (Airy)
     hd12.8  y'' + (a + cos t) y = 0, y(0)=1, y'(0)=0; y = 0 at t = 1      (Hill / Mathieu type)
     hd12.9  th'' = -cos th + a sin th cos th/(2 Sqrt[1 + a cos^2 th]),
             th(0)=Pi/2, th'(0)=-1; sin th = 0 at t = 5                  (exp 187's limit ODE, T = 5)
     hd12.10 x'' - a (1 - x^2) x' + x = 0, x(0)=2, x'(0)=0; x = 0 at t = 2 (van der Pol)

   Every case has exactly one root in its bracket; closed forms exist for
   1-7 (see truth.tsv and siam100_refs/hd12_ref.py).

   Asked of each tier (the statement of exp 187):
     machine: FindRoot (Brent, the bracket shown) on a of the time the first
              event fires, found by NDSolve + WhenEvent["StopIntegration"]
              (t <= 100), NDSolve PrecisionGoal/AccuracyGoal 12;
     wp50:    the same at WorkingPrecision -> 60, PrecisionGoal/AccuracyGoal 50
              for both NDSolve and FindRoot.  Note: WL's DEFAULT NDSolve method
              fails at WorkingPrecision 60 on several of these ODEs
              (NDSolve::nderr / ndcf at t ~ 1e-9 .. 1.5, then FindRoot::bbrac),
              and Brent at 60 digits needs > 90 s on others; Method ->
              "Extrapolation" works and is what the reference uses
              (siam100_refs/hd12_ref.wl).  The case keeps the default method,
              as experiment 186 does.
   True values and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NDSolve", "WhenEvent", "FindRoot"}];

(* the time at which the first event fires (domain end; 100 if it never fires) *)
hitT[eqs_, ev_, f_, MachinePrecision] := (f /. First[NDSolve[Append[eqs, WhenEvent[ev, "StopIntegration"]],
    f, {t, 0, 100}, PrecisionGoal -> 12, AccuracyGoal -> 12]])["Domain"][[1, 2]];
hitT[eqs_, ev_, f_, wp_] := (f /. First[NDSolve[Append[eqs, WhenEvent[ev, "StopIntegration"]],
    f, {t, 0, 100}, WorkingPrecision -> wp, PrecisionGoal -> 50, AccuracyGoal -> 50,
    MaxSteps -> Infinity]])["Domain"][[1, 2]];

h[1][a_?NumericQ, wp_] := hitT[{y''[t] == -1 - a y'[t]^2, y[0] == 0, y'[0] == 2}, y'[t] == 0, y, wp];
h[2][a_?NumericQ, wp_] := hitT[{th''[t] == -a Sin[th[t]], th[0] == 0, th'[0] == 1}, th'[t] == 0, th, wp];
h[3][a_?NumericQ, wp_] := hitT[{x''[t] == -x[t] - a x[t]^3, x[0] == 1, x'[0] == 0}, x[t] == 0, x, wp];
h[4][a_?NumericQ, wp_] := hitT[{x''[t] + a x'[t] + x[t] == 0, x[0] == 1, x'[0] == 0}, x[t] == 0, x, wp];
h[5][a_?NumericQ, wp_] := hitT[{y'[t] == a + y[t]^2, y[0] == 0}, y[t] == 1, y, wp];
h[6][a_?NumericQ, wp_] := hitT[{r''[t] == -a/r[t]^2, r[0] == 1, r'[0] == 0}, r[t] == 1/2, r, wp];
h[7][a_?NumericQ, wp_] := hitT[{y''[t] + a t y[t] == 0, y[0] == 1, y'[0] == 0}, y[t] == 0, y, wp];
h[8][a_?NumericQ, wp_] := hitT[{y''[t] + (a + Cos[t]) y[t] == 0, y[0] == 1, y'[0] == 0}, y[t] == 0, y, wp];
h[9][a_?NumericQ, wp_] := hitT[{th''[t] == -Cos[th[t]] + a Sin[th[t]] Cos[th[t]]/(2 Sqrt[1 + a Cos[th[t]]^2]),
    th[0] == Pi/2, th'[0] == -1}, Sin[th[t]] == 0, th, wp];
h[10][a_?NumericQ, wp_] := hitT[{x''[t] - a (1 - x[t]^2) x'[t] + x[t] == 0, x[0] == 2, x'[0] == 0}, x[t] == 0, x, wp];

(* {target time, bracket} *)
spec = {{1, 1/2, 3}, {2, 2/5, 2}, {6/5, 0, 3}, {2, 1/10, 3/2}, {1/2, 1, 4},
        {1, 1/2, 1}, {1, 3, 15}, {1, 1/2, 3}, {5, 1057/100, 1058/100}, {2, 1/5, 2}};

Do[With[{k = k, T = spec[[k, 1]], lo = spec[[k, 2]], hi = spec[[k, 3]]},
  vrun["hd12." <> ToString[k] <> " machine",
       ruleVal[FindRoot[h[k][a, MachinePrecision] == T, {a, N[lo], N[hi]}, Method -> "Brent"], a], 10, 1, 60];
  vrun["hd12." <> ToString[k] <> " wp50",
       ruleVal[FindRoot[h[k][a, 60] == T, {a, lo, hi}, Method -> "Brent", WorkingPrecision -> 60,
                        AccuracyGoal -> 50, PrecisionGoal -> 50], a], 48, 1, 90]], {k, 10}];
