(* bench-timeout: 600 *)
(* ==========================================================================
   Experiment 188 -- SIAM 100-Digit Challenge, Appendix D, Problem 13
   (Trefethen): a threshold in a nonnormal nonlinear ODE.

       u = (x, y, z),  u(0) = (0, 0, z0),  z0 >= 0,
       x' = -x + 10 y + |u| (-0.7 y - 0.03 z)
       y' = -y + 10 z + |u| ( 0.7 y - 0.1 z)
       z' = -z        + |u| ( 0.03 x + 0.1 y)
       If |u(50)| = 1, what is z0 ?                answer = 0.0064823961796...

   Landscape (siam100_refs/d13_landscape.py): u = 0 is stable, but the linear
   part is strongly nonnormal (transient growth ~100x).  There is a saddle
   equilibrium with |u*| = 0.21523 (one unstable eigenvalue 0.2995) whose stable
   manifold cuts the z-axis at z_c = 0.00648239617...  For z0 < z_c the orbit
   decays (max |u| < 0.2153, so |u(50)| = 1 is impossible); for z0 > z_c it
   leaves the saddle and blows up in finite time -- before t = 50 already for
   z0 > z_c (1 + ~1e-12).  So the question has exactly ONE answer (multiple
   solutions were a possibility the problem statement leaves open; there are
   none), a hair above z_c, and |u(50)| as a function of z0 is a needle: it is
   defined only on a relative window of width ~1e-12.  The decimal
   coefficients are exact (7/10, 3/100, 1/10).

   Asked of each tier: since |u| increases monotonically once the orbit leaves
   the saddle, |u(50)| = 1 is the same as "|u| first reaches 1 at t = 50",
   which is what an event-capable integrator states directly:
     machine: FindRoot (Brent, bracket [0.0064, 0.0065]) on z0 of the time at
              which |u|^2 = 1 (NDSolve + WhenEvent, t <= 100; never
              reached -> 100), NDSolve tolerances 1e-12.
     wp50:    the same at WorkingPrecision 60, PrecisionGoal 50, started
              (secant) from z0 = 0.00648239617962 / 0.00648239617963 -- the
              mpmath row has no events and must start inside the needle.
   True value and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NDSolve", "WhenEvent", "FindRoot"}];

reach[z0_?NumericQ] := (x /. First[NDSolve[{
    x'[t] == -x[t] + 10 y[t] + Sqrt[x[t]^2 + y[t]^2 + z[t]^2] (-0.7 y[t] - 0.03 z[t]),
    y'[t] == -y[t] + 10 z[t] + Sqrt[x[t]^2 + y[t]^2 + z[t]^2] (0.7 y[t] - 0.1 z[t]),
    z'[t] == -z[t] + Sqrt[x[t]^2 + y[t]^2 + z[t]^2] (0.03 x[t] + 0.1 y[t]),
    x[0] == 0, y[0] == 0, z[0] == z0,
    WhenEvent[x[t]^2 + y[t]^2 + z[t]^2 == 1, "StopIntegration"]}, x, {t, 0, 100},
    PrecisionGoal -> 12, AccuracyGoal -> 12]])["Domain"][[1, 2]];

reach60[z0_?NumericQ] := (x /. First[NDSolve[{
    x'[t] == -x[t] + 10 y[t] + Sqrt[x[t]^2 + y[t]^2 + z[t]^2] (-7/10 y[t] - 3/100 z[t]),
    y'[t] == -y[t] + 10 z[t] + Sqrt[x[t]^2 + y[t]^2 + z[t]^2] (7/10 y[t] - 1/10 z[t]),
    z'[t] == -z[t] + Sqrt[x[t]^2 + y[t]^2 + z[t]^2] (3/100 x[t] + 1/10 y[t]),
    x[0] == 0, y[0] == 0, z[0] == z0,
    WhenEvent[x[t]^2 + y[t]^2 + z[t]^2 == 1, "StopIntegration"]}, x, {t, 0, 100},
    WorkingPrecision -> 60, PrecisionGoal -> 50, AccuracyGoal -> 50,
    MaxSteps -> Infinity]])["Domain"][[1, 2]];

vrun["d13.machine",
     ruleVal[FindRoot[reach[z0] == 50, {z0, 0.0064, 0.0065}, Method -> "Brent"], z0], 10];
vrun["d13.wp50",
     ruleVal[FindRoot[reach60[z0] == 50, {z0, 648239617962/10^14, 648239617963/10^14},
                      WorkingPrecision -> 60, AccuracyGoal -> 50,
                      PrecisionGoal -> 50], z0], 48, 1, 400];
