(* bench-timeout: 1800 *)
(* ==========================================================================
   Experiment 247 -- SIAM 100-Digit Challenge, sealed holdout B for Appendix D
   Problem 12 (find alpha such that the particle FIRST hits y = 0 at t = 10).

   Capability probed: a scalar root finder (FindRoot, Brent) on the time of
   the first ODE EVENT (NDSolve + WhenEvent["StopIntegration"], the event
   time read off as the end of the solution's domain; no event by t = 100 ->
   100).  Ten problems "find p such that the first event happens at t = T*":

     bd12.1  pendulum th'' = -sin th, th(0)=0, th'(0)=v; first turning point
             (th' = 0) at t = 2.  v in [1, 1.9]                v = 1.60481...
     bd12.2  vertical throw with quadratic drag y'' = -1 - y'|y'|/10, y(0)=0,
             y'(0)=v; lands (y < 0) at t = 3.  v in [0.5, 5]    v = 1.58827...
     bd12.3  van der Pol x'' = mu (1 - x^2) x' - x, x(0)=2, x'(0)=0; first
             x = 0 at t = 2.  mu in [0, 2]                     mu = 0.75203...
     bd12.4  Lotka-Volterra x' = x(1-y), y' = y(x-1), x(0)=a, y(0)=1; first
             x = 1 at t = 1.  a in [1.1, 6]                     a = 2.63280...
     bd12.5  Kepler r'' = -r/|r|^3 from (1,0) with velocity (0,v); first
             x = 0 at t = 13/10.  v in [0.2, 0.95]              v = 0.57711...
     bd12.6  damped oscillator x'' + 2 zeta x' + x = 0, x(0)=1, x'(0)=0;
             first x = 0 at t = 2.  zeta in [0, 0.9]         zeta = 0.31902...
     bd12.7  Riccati y' = alpha + y^2, y(0)=0; first y = 1 at t = 1.
             alpha in [0.1, 2]                              alpha = 0.74017...
     bd12.8  elastic pendulum (stiffness k, rest length 1, gravity 1) released
             at rest from (0.66, -0.88); first x = 0 at t = 2.  k in [2, 5]
                                                                k = 3.98319...
     bd12.9  SIR S' = -b S I, I' = b S I - I, S(0)=0.99, I(0)=0.01; epidemic
             peak (b S = 1) at t = 5.  b in [1.5, 3]            b = 1.81218...
     bd12.10 pendulum th(0)=0, th'(0)=v; first th = 3 at t = 4: the root sits
             just past the fold v_c = 2 sin(3/2) below which th never reaches
             3 (the event time jumps to 100), so the hitting time is steep and
             the bracket [1.9, 2.1] straddles a jump.   v = 1.99499781...
             Trap: just past the fold th rises above 3 only briefly (up- and
             down-crossing close together), so an event locator that only
             compares signs at step ends steps over both crossings (scipy's
             DOP853 does: 4 digits).

   machine : bracket above, NDSolve PrecisionGoal/AccuracyGoal 12.
   wp50    : the same at WorkingPrecision 60 (NDSolve PrecisionGoal 50,
             FindRoot Precision/AccuracyGoal 50), on a NARROW bracket around
             the root (the one the event-less rows -- mpmath, Maxima -- need
             for their fixed-time statement, as in exp. 187).
   True values and provenance: truth.tsv; references siam100_refs/bd12_ref.py.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NDSolve", "WhenEvent", "FindRoot"}];

opts[MachinePrecision] = {PrecisionGoal -> 12, AccuracyGoal -> 12};
opts[60] = {WorkingPrecision -> 60, PrecisionGoal -> 50, AccuracyGoal -> 50, MaxSteps -> Infinity};
fr[MachinePrecision] = {Method -> "Brent"};
fr[60] = {Method -> "Brent", WorkingPrecision -> 60, AccuracyGoal -> 50, PrecisionGoal -> 50};
br[{a_, b_}, MachinePrecision] := N[{a, b}];
br[{a_, b_}, 60] := {a, b};

(* end of the solution's domain = the first event time (or 100) *)
tEnd[eqs_, fns_List, wp_] :=
  (First[fns] /. First[NDSolve[eqs, fns, {t, 0, 100}, Sequence @@ opts[wp]]])["Domain"][[1, 2]];

ev[1][v_?NumericQ, wp_] := tEnd[{w''[t] == -Sin[w[t]], w[0] == 0, w'[0] == v,
    WhenEvent[w'[t] == 0, "StopIntegration"]}, {w}, wp];
ev[2][v_?NumericQ, wp_] := tEnd[{w''[t] == -1 - w'[t] Abs[w'[t]]/10, w[0] == 0, w'[0] == v,
    WhenEvent[w[t] < 0, "StopIntegration"]}, {w}, wp];
ev[3][mu_?NumericQ, wp_] := tEnd[{w''[t] == mu (1 - w[t]^2) w'[t] - w[t], w[0] == 2, w'[0] == 0,
    WhenEvent[w[t] == 0, "StopIntegration"]}, {w}, wp];
ev[4][a_?NumericQ, wp_] := tEnd[{x'[t] == x[t] (1 - y[t]), y'[t] == y[t] (x[t] - 1), x[0] == a, y[0] == 1,
    WhenEvent[x[t] == 1, "StopIntegration"]}, {x, y}, wp];
ev[5][v_?NumericQ, wp_] := tEnd[{x''[t] == -x[t]/(x[t]^2 + y[t]^2)^(3/2),
    y''[t] == -y[t]/(x[t]^2 + y[t]^2)^(3/2), x[0] == 1, y[0] == 0, x'[0] == 0, y'[0] == v,
    WhenEvent[x[t] == 0, "StopIntegration"]}, {x, y}, wp];
ev[6][z_?NumericQ, wp_] := tEnd[{w''[t] == -2 z w'[t] - w[t], w[0] == 1, w'[0] == 0,
    WhenEvent[w[t] == 0, "StopIntegration"]}, {w}, wp];
ev[7][a_?NumericQ, wp_] := tEnd[{w'[t] == a + w[t]^2, w[0] == 0,
    WhenEvent[w[t] == 1, "StopIntegration"]}, {w}, wp];
ev[8][k_?NumericQ, wp_] := tEnd[{
    x''[t] == -k (1 - 1/Sqrt[x[t]^2 + y[t]^2]) x[t],
    y''[t] == -k (1 - 1/Sqrt[x[t]^2 + y[t]^2]) y[t] - 1,
    x[0] == 66/100, y[0] == -88/100, x'[0] == 0, y'[0] == 0,
    WhenEvent[x[t] == 0, "StopIntegration"]}, {x, y}, wp];
ev[9][b_?NumericQ, wp_] := tEnd[{x'[t] == -b x[t] y[t], y'[t] == b x[t] y[t] - y[t],
    x[0] == 99/100, y[0] == 1/100, WhenEvent[b x[t] == 1, "StopIntegration"]}, {x, y}, wp];
ev[10][v_?NumericQ, wp_] := tEnd[{w''[t] == -Sin[w[t]], w[0] == 0, w'[0] == v,
    WhenEvent[w[t] == 3, "StopIntegration"]}, {w}, wp];

(* {target time, machine bracket, narrow wp50 bracket} *)
spec = {
  {2, {1, 19/10}, {160/100, 161/100}},
  {3, {1/2, 5}, {158/100, 159/100}},
  {2, {0, 2}, {75/100, 76/100}},
  {1, {11/10, 6}, {263/100, 264/100}},
  {13/10, {1/5, 19/20}, {57/100, 58/100}},
  {2, {0, 9/10}, {31/100, 32/100}},
  {1, {1/10, 2}, {74/100, 75/100}},
  {2, {2, 5}, {398/100, 399/100}},
  {5, {3/2, 3}, {181/100, 182/100}},
  {4, {19/10, 21/10}, {199499/100000, 1995/1000}}};

solve[k_, wp_] := Module[{p, T = spec[[k, 1]], b = br[spec[[k, If[wp === 60, 3, 2]]], wp]},
  ruleVal[FindRoot[ev[k][p, wp] == T, {p, b[[1]], b[[2]]}, Evaluate[Sequence @@ fr[wp]]], p]];

Do[vrun["bd12." <> ToString[k] <> " machine", solve[k, MachinePrecision], 10], {k, 1, 10}];
Do[vrun["bd12." <> ToString[k] <> " wp50", solve[k, 60], 48], {k, 1, 10}];
