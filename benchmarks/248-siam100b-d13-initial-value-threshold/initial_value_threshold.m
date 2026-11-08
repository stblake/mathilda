(* bench-timeout: 1800 *)
(* ==========================================================================
   Experiment 248 -- SIAM 100-Digit Challenge, sealed holdout B for Appendix D
   Problem 13 (find z0 with |u(50)| = 1 for a nonnormal nonlinear ODE).

   Capability probed: find an INITIAL value such that the solution meets a
   condition at a prescribed later time -- typically by passing close to an
   unstable equilibrium / stable manifold, so the answer is very sensitive to
   the integration.  Where the condition is "g first reaches 0 at t = T" the
   case is the event statement (NDSolve + WhenEvent, Brent on the event time,
   no event by t = 100 -> 100); where g(u(T)) is not monotone in time the case
   is the literal fixed-time statement (FindRoot, secant, on g(u(T)) = 0).

     bd13.1  y' = y^2 - y (unstable equilibrium y = 1); y first = 1/2 at t = 10.
             event, y0 in [0.99, 0.999999]          y0 = 1/(1+e^-10) = 0.99995...
     bd13.2  u' = A u, A = -I + 10 (superdiagonal) on R^3 (a nonnormal Jordan
             block), u(0) = (0, 0, z0); |u(10)| = 1.  fixed-time, secant from
             4, 5                                                z0 = 4.40441...
     bd13.3  Riccati y' = y^2 - t; y(3) = 0.  fixed-time, secant from 0.727,
             0.728 (a slightly larger y0 blows up before t = 3)  y0 = 0.72769...
     bd13.4  exp. 188's model (x' = -x + 10y + |u|(-0.7y - 0.03z), ...),
             u(0) = (0, 0, z0); |u| first = 1/2 at t = 40 (a hair above the
             saddle's stable manifold).  event, z0 in [0.0064, 0.0065]
                                                              z0 = 0.0064824426...
     bd13.5  damped pendulum th'' = -sin th - th'/10, th(0) = 0, th'(0) = v;
             th first = 4 (over the top) at t = 10.  event, v in [2.2, 2.3]
                                                                  v = 2.20403...
     bd13.6  radial Kepler r'' = -1/r^2, r(0) = 1, r'(0) = v; r first = 3 at
             t = 5 (just past the fold v = 2/Sqrt[3] where the apex is 3).
             event, v in [1.1, 1.3]                               v = 1.15503...
     bd13.7  Allee growth x' = x(1-x)(x - 3/10); x first = 9/10 at t = 10.
             event, x0 in [0.31, 0.5]                            x0 = 0.36834...
     bd13.8  Blasius f''' = -f f''/2, f(0) = f'(0) = 0, f''(0) = s; f'(10) = 1.
             fixed-time, secant from 0.3, 0.4                     s = 0.33205...
     bd13.9  seasonally harvested logistic x' = x(1-x) - (1 + sin 2 pi t)/10;
             x first = 1/2 at t = 10.  event, x0 in [0.12, 0.2]  x0 = 0.12885...
     bd13.10 Falkner-Skan (Hiemenz) f''' = -f f'' - 1 + f'^2, f(0) = f'(0) = 0,
             f''(0) = s; f'(10) = 1 (exponentially sensitive shooting).
             fixed-time, secant from 1.2, 1.25                    s = 1.23258...

   machine : as above, NDSolve PrecisionGoal/AccuracyGoal 12.
   wp50    : the same statements at WorkingPrecision 60 (NDSolve and FindRoot
             PrecisionGoal/AccuracyGoal 50) from a NARROW bracket / start pair
             around the root -- the one the event-less mpmath row needs too.
   True values and provenance: truth.tsv; references siam100_refs/bd13_ref.py.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NDSolve", "WhenEvent", "FindRoot"}];

opts[MachinePrecision] = {PrecisionGoal -> 12, AccuracyGoal -> 12};
opts[60] = {WorkingPrecision -> 60, PrecisionGoal -> 50, AccuracyGoal -> 50, MaxSteps -> Infinity};
frw[MachinePrecision] = {};
frw[60] = {WorkingPrecision -> 60, AccuracyGoal -> 50, PrecisionGoal -> 50};
nm[l_, MachinePrecision] := N[l];
nm[l_, 60] := l;

(* the event time, or 100 when the event never fires: the end of the domain is
   the event time only if the event condition holds there -- a solution that
   blows up first (bd13.9 below its threshold: w -> -Infinity) ends at its
   singularity, which is "no event" (100), as in the other halves *)
tEnd[eqs_, fns_List, wp_] := Module[{sol, te},
  sol = First[NDSolve[eqs, fns, {t, 0, 100}, Sequence @@ opts[wp]]];
  te = (First[fns] /. sol)["Domain"][[1, 2]];
  If[Abs[First[Cases[eqs, WhenEvent[a_ == b_, __] :> a - b]] /. t -> te /. sol] < 10^-6, te, 100]];
atT[eqs_, fns_List, T_, what_, wp_] := what /. First[NDSolve[eqs, fns, {t, 0, T}, Sequence @@ opts[wp]]];

c7 = {7/10, 3/100, 1/10};
model[z0_, c_] := {
    x'[t] == -x[t] + 10 y[t] + Sqrt[x[t]^2 + y[t]^2 + z[t]^2] (-c[[1]] y[t] - c[[2]] z[t]),
    y'[t] == -y[t] + 10 z[t] + Sqrt[x[t]^2 + y[t]^2 + z[t]^2] (c[[1]] y[t] - c[[3]] z[t]),
    z'[t] == -z[t] + Sqrt[x[t]^2 + y[t]^2 + z[t]^2] (c[[2]] x[t] + c[[3]] y[t]),
    x[0] == 0, y[0] == 0, z[0] == z0};

(* event problems: the event time *)
ev[1][p_?NumericQ, wp_] := tEnd[{w'[t] == w[t]^2 - w[t], w[0] == p,
    WhenEvent[w[t] == 1/2, "StopIntegration"]}, {w}, wp];
ev[4][p_?NumericQ, wp_] := tEnd[Join[model[p, If[wp === 60, c7, N[c7]]],
    {WhenEvent[x[t]^2 + y[t]^2 + z[t]^2 == 1/4, "StopIntegration"]}], {x, y, z}, wp];
ev[5][p_?NumericQ, wp_] := tEnd[{w''[t] == -Sin[w[t]] - w'[t]/10, w[0] == 0, w'[0] == p,
    WhenEvent[w[t] == 4, "StopIntegration"]}, {w}, wp];
ev[6][p_?NumericQ, wp_] := tEnd[{w''[t] == -1/w[t]^2, w[0] == 1, w'[0] == p,
    WhenEvent[w[t] == 3, "StopIntegration"]}, {w}, wp];
ev[7][p_?NumericQ, wp_] := tEnd[{w'[t] == w[t] (1 - w[t]) (w[t] - 3/10), w[0] == p,
    WhenEvent[w[t] == 9/10, "StopIntegration"]}, {w}, wp];
ev[9][p_?NumericQ, wp_] := tEnd[{w'[t] == w[t] (1 - w[t]) - (1 + Sin[2 Pi t])/10, w[0] == p,
    WhenEvent[w[t] == 1/2, "StopIntegration"]}, {w}, wp];
(* fixed-time problems: the residual at T *)
fx[2][p_?NumericQ, wp_] := Norm[atT[{x'[t] == -x[t] + 10 y[t], y'[t] == -y[t] + 10 z[t], z'[t] == -z[t],
    x[0] == 0, y[0] == 0, z[0] == p}, {x, y, z}, 10, {x[10], y[10], z[10]}, wp]] - 1;
fx[3][p_?NumericQ, wp_] := atT[{w'[t] == w[t]^2 - t, w[0] == p}, {w}, 3, w[3], wp];
fx[8][p_?NumericQ, wp_] := atT[{w'''[t] == -w[t] w''[t]/2, w[0] == 0, w'[0] == 0, w''[0] == p},
    {w}, 10, w'[10], wp] - 1;
fx[10][p_?NumericQ, wp_] := atT[{w'''[t] == -w[t] w''[t] - 1 + w'[t]^2, w[0] == 0, w'[0] == 0,
    w''[0] == p}, {w}, 10, w'[10], wp] - 1;

(* {kind, target time, machine bracket/starts, narrow wp50 bracket/starts} *)
spec = {
  {"ev", 10, {99/100, 999999/1000000}, {99995/100000, 99996/100000}},
  {"fx", 0, {4, 5}, {44/10, 45/10}},
  {"fx", 0, {727/1000, 728/1000}, {7276/10000, 7278/10000}},
  {"ev", 40, {64/10000, 65/10000}, {6482442672312/10^15, 6482442672313/10^15}},
  {"ev", 10, {22/10, 23/10}, {2204/1000, 22041/10000}},
  {"ev", 5, {11/10, 13/10}, {115503/100000, 115505/100000}},
  {"ev", 10, {31/100, 1/2}, {368/1000, 369/1000}},
  {"fx", 0, {3/10, 4/10}, {332/1000, 333/1000}},
  {"ev", 10, {12/100, 2/10}, {1288/10000, 1289/10000}},
  {"fx", 0, {12/10, 125/100}, {12325/10000, 12326/10000}}};

solve[k_, wp_] := Module[{p, s = spec[[k]], b},
  b = nm[s[[If[wp === 60, 4, 3]]], wp];
  If[s[[1]] === "ev",
    ruleVal[FindRoot[ev[k][p, wp] == s[[2]], {p, b[[1]], b[[2]]}, Method -> "Brent",
                     Evaluate[Sequence @@ frw[wp]]], p],
    ruleVal[FindRoot[fx[k][p, wp] == 0, {p, b[[1]], b[[2]]}, Evaluate[Sequence @@ frw[wp]]], p]]];

Do[vrun["bd13." <> ToString[k] <> " machine", solve[k, MachinePrecision], 10], {k, 1, 10}];
Do[vrun["bd13." <> ToString[k] <> " wp50", solve[k, 60], 48], {k, 1, 10}];
