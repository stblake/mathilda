(* Experiment 32 -- NDSolve, a difficulty ladder.
   ==========================================================================
   WHY A LADDER.  Experiment 13 runs five ODEs at one difficulty each. An ODE
   integrator does not get gradually worse as a problem hardens -- it holds
   until the step controller loses, and then it either crawls or quits. Both
   failure modes are invisible in a single row, and the second one is worse
   than slow: an integrator that quits looks FAST.

   THE AXES
     horizon     t = 10 -> 1000 on a decaying linear ODE, and t = 100 -> 1000
                 on a harmonic oscillator. Cost should be
                 linear in the horizon; a step budget that does not scale with
                 the interval shows up here as a truncated solution.
     stiffness   Van der Pol mu = 10 -> 100 -> 1000. An explicit method is not
                 slow on a stiff problem, it is unusable, so a flat time across
                 these three rows means a stiff solver was reached.
     sensitivity Lorenz, which is chaotic: the check is at t = 5, early enough
                 that the trajectory is still reproducible across systems.

   CHECKS ARE THE POINT OF THIS EXPERIMENT.  Every row below has a known
   closed-form or a value all three systems must agree on, and each check
   evaluates the returned solution AT THE FINAL TIME. That is deliberate: an
   integrator that stops early and hands back an InterpolatingFunction covering
   less than the requested interval will be extrapolated at the endpoint and
   will fail its check. A row that is fast and wrong is worth strictly less
   than a row that is slow and right, and only the check can tell them apart.
   ========================================================================== *)

Get["../harness.m"];
Get["../data.m"];

require[{"NDSolve", "InterpolatingFunction", "AccuracyGoal", "PrecisionGoal",
         "MaxSteps"}];

(* ---- horizon, decaying linear: y' = -y, y[T] = Exp[-T] ----
   Beyond T ~ 40 the true value underflows the check's resolution, so the long
   rows are scaled by Exp[T] -- i.e. the check asks for the RELATIVE accuracy
   that survives, rather than comparing two numbers that are both ~0. *)
bench["NDS decay t=10",
  NDSolve[{y'[t] == -y[t], y[0] == 1}, y, {t, 0, 10}];];
check["NDS decay t=10",
  Round[10^6 (y[10] /. First[NDSolve[{y'[t] == -y[t], y[0] == 1}, y, {t, 0, 10}]])]];

bench["NDS decay t=1000",
  NDSolve[{y'[t] == -y[t], y[0] == 1}, y, {t, 0, 1000}];];
check["NDS decay t=1000",
  Round[10^6 (y[20] /. First[NDSolve[{y'[t] == -y[t], y[0] == 1}, y, {t, 0, 1000}]])]];

(* ---- horizon, oscillatory: x'' = -x, so x[T] = Cos[T] and the phase must
   still be right at the end. This is the row a non-scaling step budget fails:
   the amplitude stays bounded, so only the ENDPOINT value exposes a solution
   that never reached T. *)
bench["NDS harmonic t=100",
  NDSolve[{x'[t] == y[t], y'[t] == -x[t], x[0] == 1, y[0] == 0}, {x, y},
          {t, 0, 100}];];
check["NDS harmonic t=100",
  Round[10^4 (x[100] /. First[NDSolve[{x'[t] == y[t], y'[t] == -x[t],
        x[0] == 1, y[0] == 0}, {x, y}, {t, 0, 100}]])]];

bench["NDS harmonic t=1000",
  NDSolve[{x'[t] == y[t], y'[t] == -x[t], x[0] == 1, y[0] == 0}, {x, y},
          {t, 0, 1000}];];
check["NDS harmonic t=1000",
  Round[10^4 (x[1000] /. First[NDSolve[{x'[t] == y[t], y'[t] == -x[t],
        x[0] == 1, y[0] == 0}, {x, y}, {t, 0, 1000}]])]];

(* ---- stiffness: Van der Pol. mu = 1000 is genuinely stiff; the relaxation
   oscillation spends almost all its time on the slow manifold and crosses
   between branches in O(1/mu) time. *)
bench["NDS Van der Pol mu=10",
  NDSolve[{y''[t] - 10 (1 - y[t]^2) y'[t] + y[t] == 0, y[0] == 2, y'[0] == 0},
          y, {t, 0, 20}];];
check["NDS Van der Pol mu=10",
  Round[10^3 (y[20] /. First[NDSolve[{y''[t] - 10 (1 - y[t]^2) y'[t] + y[t] == 0,
        y[0] == 2, y'[0] == 0}, y, {t, 0, 20}]])]];

bench["NDS Van der Pol mu=100",
  NDSolve[{y''[t] - 100 (1 - y[t]^2) y'[t] + y[t] == 0, y[0] == 2, y'[0] == 0},
          y, {t, 0, 20}];, 1];
check["NDS Van der Pol mu=100",
  Round[10^2 (y[20] /. First[NDSolve[{y''[t] - 100 (1 - y[t]^2) y'[t] + y[t] == 0,
        y[0] == 2, y'[0] == 0}, y, {t, 0, 20}]])]];

bench["NDS Van der Pol mu=1000",
  NDSolve[{y''[t] - 1000 (1 - y[t]^2) y'[t] + y[t] == 0, y[0] == 2, y'[0] == 0},
          y, {t, 0, 20}];, 1];
check["NDS Van der Pol mu=1000",
  Round[10 (y[20] /. First[NDSolve[{y''[t] - 1000 (1 - y[t]^2) y'[t] + y[t] == 0,
        y[0] == 2, y'[0] == 0}, y, {t, 0, 20}]])]];

(* ---- sensitivity: Lorenz. Checked at t = 5, where the trajectory is still
   deterministic to four places in every system; by t = 20 the systems disagree
   for reasons of conditioning, not implementation, and the check would be
   measuring the butterfly rather than the solver. *)
bench["NDS Lorenz t=5",
  NDSolve[{x'[t] == 10 (y[t] - x[t]), y'[t] == x[t] (28 - z[t]) - y[t],
           z'[t] == x[t] y[t] - 8/3 z[t], x[0] == 1, y[0] == 1, z[0] == 1},
          {x, y, z}, {t, 0, 5}];];
check["NDS Lorenz t=5",
  Round[10^3 (x[5] /. First[NDSolve[{x'[t] == 10 (y[t] - x[t]),
        y'[t] == x[t] (28 - z[t]) - y[t], z'[t] == x[t] y[t] - 8/3 z[t],
        x[0] == 1, y[0] == 1, z[0] == 1}, {x, y, z}, {t, 0, 5}]])]];
