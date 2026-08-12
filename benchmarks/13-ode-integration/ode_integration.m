(* ==========================================================================
   Experiment 13 -- Numerical ODE integration (NDSolve)
   ==========================================================================
   WHAT IT MEASURES.  src/numerical_calculus/ndsolve*.c -- explicit Runge-Kutta
   (ndsolve_rk.c), Adams (ndsolve_adams.c), and the implicit/stiff path
   (ndsolve_implicit.c), plus ndsolve_compile.c which lowers the right-hand side
   through Compile[].

   NON-STIFF AND STIFF ARE DIFFERENT MEASUREMENTS.  On a stiff system an
   explicit method is not slow, it is unusable: the step size collapses to keep
   stability. So the stiff row measures whether an implicit path EXISTS, and the
   check (the solution at the endpoint) is what proves the answer is real.

   Checks are the solution VALUE at the final time, rounded.
   ========================================================================== *)

Get["../harness.m"];
Get["../data.m"];

require[{"NDSolve", "Method", "MaxSteps", "WorkingPrecision", "InterpolatingFunction"}];

endv[sol_, f_, t1_] := Round[10^4 (f[t1] /. sol)];

(* 1. Exponential decay: the smoothest possible non-stiff problem. *)
bench["NDSolve y'=-y on [0,10]",
  NDSolve[{y'[t] == -y[t], y[0] == 1}, y, {t, 0, 10}];];
check["NDSolve y'=-y on [0,10]",
  endv[First[NDSolve[{y'[t] == -y[t], y[0] == 1}, y, {t, 0, 10}]], y, 10]];

(* 2. The harmonic oscillator as a 2-system: tests coupling, not stiffness. *)
bench["NDSolve harmonic oscillator",
  NDSolve[{x'[t] == y[t], y'[t] == -x[t], x[0] == 1, y[0] == 0},
          {x, y}, {t, 0, 20}];];
check["NDSolve harmonic oscillator",
  endv[First[NDSolve[{x'[t] == y[t], y'[t] == -x[t], x[0] == 1, y[0] == 0},
                     {x, y}, {t, 0, 20}]], x, 20]];

(* 3. Van der Pol at mu = 10: genuinely stiff. *)
bench["NDSolve Van der Pol mu=10",
  NDSolve[{x''[t] == 10 (1 - x[t]^2) x'[t] - x[t], x[0] == 2, x'[0] == 0},
          x, {t, 0, 20}], 1];
check["NDSolve Van der Pol mu=10",
  endv[First[NDSolve[{x''[t] == 10 (1 - x[t]^2) x'[t] - x[t],
                      x[0] == 2, x'[0] == 0}, x, {t, 0, 20}]], x, 20]];

(* 4. A nonlinear scalar equation with a known closed form: y' = y^2 t. *)
bench["NDSolve y'=y^2 t nonlinear",
  NDSolve[{y'[t] == y[t]^2 t, y[0] == 1/2}, y, {t, 0, 1}];];
check["NDSolve y'=y^2 t nonlinear",
  endv[First[NDSolve[{y'[t] == y[t]^2 t, y[0] == 1/2}, y, {t, 0, 1}]], y, 1]];

(* 5. A longer integration: 200 time units, so per-step cost dominates. *)
bench["NDSolve long horizon, t to 200",
  NDSolve[{y'[t] == -y[t]/10, y[0] == 1}, y, {t, 0, 200}];];
check["NDSolve long horizon, t to 200",
  endv[First[NDSolve[{y'[t] == -y[t]/10, y[0] == 1}, y, {t, 0, 200}]], y, 200]];
