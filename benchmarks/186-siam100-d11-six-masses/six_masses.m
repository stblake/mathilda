(* bench-timeout: 300 *)
(* ==========================================================================
   Experiment 186 -- SIAM 100-Digit Challenge, Appendix D Problem 11
   ("Six masses"; problem 13 on the book's web page).

     Six unit masses are fixed at (2,-1), (2,0), (2,1), (3,-1), (3,0), (3,1).
     A unit mass starts at (0,0) with speed 1 in direction theta and moves
     under the inverse-square attraction r^-2 of the six.  What is the
     shortest time in which it can reach (4,1)?
                                         Answer: 1.9160969183320049103...
                                         (theta = -0.27488125402627671...)

   Idiomatic statement: shooting.  NDSolve integrates
       r'' = Sum_i (m_i - r)/|m_i - r|^3,  r(0) = 0,  r'(0) = (cos th, sin th)
   to t = T, and FindRoot solves {x(T), y(T)} == {4, 1} for {th, T}.
   The problem has MANY solutions (a theta scan finds eight direct branches,
   e.g. th ~ 0.5513 arrives at T ~ 1.91763); the SHORTEST arrival is the
   th ~ -0.2749 branch.  Finding the branch is the scan (done once, by hand,
   in the reference: siam100_refs/d11_ref.wl and .py); every system's case starts the
   2-unknown root-finder at the same rounded guess th = -0.27, T = 1.92 and is
   scored on the T it converges to.
     machine : NDSolve PrecisionGoal/AccuracyGoal 13, FindRoot defaults;
     wp50    : NDSolve and FindRoot at WorkingPrecision -> 60
               (NDSolve PrecisionGoal/AccuracyGoal 55).  Note: WL's DEFAULT
               NDSolve method fails at WorkingPrecision 60-80 on this ODE
               (NDSolve::nderr / ndtol at t ~ 0.3-0.7), so the idiomatic wp50
               case times out in WL; Method -> "Extrapolation" works and is
               what the reference uses.  The case keeps the default method.
   True value and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NDSolve", "FindRoot"}];

masses = {{2, -1}, {2, 0}, {2, 1}, {3, -1}, {3, 0}, {3, 1}};
odes[th_] := {
   x''[t] == Sum[(m[[1]] - x[t])/((m[[1]] - x[t])^2 + (m[[2]] - y[t])^2)^(3/2), {m, masses}],
   y''[t] == Sum[(m[[2]] - y[t])/((m[[1]] - x[t])^2 + (m[[2]] - y[t])^2)^(3/2), {m, masses}],
   x[0] == 0, y[0] == 0, x'[0] == Cos[th], y'[0] == Sin[th]};

posM[th_?NumericQ, T_?NumericQ] := {x[T], y[T]} /. First[NDSolve[odes[th], {x, y}, {t, 0, T},
    PrecisionGoal -> 13, AccuracyGoal -> 13, MaxSteps -> Infinity]];
pos50[th_?NumericQ, T_?NumericQ] := {x[T], y[T]} /. First[NDSolve[odes[th], {x, y}, {t, 0, T},
    WorkingPrecision -> 60, PrecisionGoal -> 55, AccuracyGoal -> 55, MaxSteps -> Infinity]];

vrun["d11.machine", T /. FindRoot[posM[th, T] == {4, 1}, {{th, -0.27}, {T, 1.92}}], 10];
vrun["d11.wp50", T /. FindRoot[pos50[th, T] == {4, 1}, {{th, -27/100}, {T, 192/100}},
    WorkingPrecision -> 60], 48];
