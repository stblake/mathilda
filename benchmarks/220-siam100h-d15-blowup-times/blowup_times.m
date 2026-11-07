(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 220 -- SIAM 100-Digit Challenge, holdout A, family d15
   (generalises App. D Problem 15, experiment 190: "at what time does the
   solution of u_t = Laplacian u + e^u on the 3x3 square blow up?").

   The capability: integrate an evolution equation INTO a finite-time
   singularity and report where the solver stops ("step size effectively
   zero; singularity suspected"), i.e. the end of the solution's domain; plus
   one heat-equation crossing time (a PDE root problem).  The original's
   2-D blow-up has only ~8 known digits, so the instances are tractable
   relatives whose answers are known to 70+ digits:

     hd15.1  y' = 1 + y^2, y(0) = 0                       T = pi/2
     hd15.2  y' = t^2 + y^2, y(0) = 0                     T = Sqrt[2 j], J_{-1/4}(j) = 0
     hd15.3  y'' = y^2, y(0) = 1, y'(0) = 0               T = Sqrt[3/2] Beta[1/6,1/2]/3
     hd15.4  y' = e^y - y, y(0) = 0                       T = Int_0^Infinity dy/(e^y - y)
     hd15.5  x' = x^2 + y^2, y' = 2 x y, (x,y)(0) = (1, 1/2)   T = 2/3 (x+y solves u' = u^2)
     hd15.6  y' = y^2 + e^t, y(0) = 0                     T = first zero of a J0(2e^{t/2}) + b Y0(2e^{t/2})
     hd15.7  u_t = u_xx + u + u^2 on [0,1], u_x = 0 at x = 0, 1, u(x,0) = 1      T = Log[2]
     hd15.8  u_t = Laplacian u + e^u on [0,3]^2, du/dn = 0, u(x,y,0) = 0         T = 1
             (the original problem with Neumann instead of Dirichlet walls)
     hd15.9  u_t = u_xx, u(0,t) = u(1,t) = 0, u(x,0) = 4x(1-x):
             the time at which u(1/2, t) = 1/2   (NDSolve + FindRoot)         T = 0.0734153955...
     hd15.10 y'' = e^y, y(0) = y'(0) = 0                  T = Pi/Sqrt[2]

   Asked of each tier (the statement of exp 190):
     machine: NDSolve with default options integrated past the blow-up; the
              answer is the end of the solution's domain.  hd15.9: NDSolve
              of the PDE on t in [0, 1], then FindRoot (Brent on [1/100, 1/5]).
     wp50:    the same with WorkingPrecision -> 60.
   True values and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NDSolve", "FindRoot"}];

dom[sol_] := sol["Domain"][[1, 2]];
(* Same request in every system: the machine tier asks for the tolerances the
   .py half gives solve_ivp -- rtol = atol = 1e-12 for the ODE rows, 1e-10 for
   the method-of-lines PDE rows; the wp50 tier keeps NDSolve's defaults. *)
ode12[MachinePrecision] := Sequence[PrecisionGoal -> 12, AccuracyGoal -> 12];
ode12[_] := Sequence[];
pde10[MachinePrecision] := Sequence[PrecisionGoal -> 10, AccuracyGoal -> 10];
pde10[_] := Sequence[];

tinf[eqs_, f_, tmax_, wp_] := dom[f /. First[NDSolve[eqs, f, {t, 0, tmax}, WorkingPrecision -> wp, Evaluate[ode12[wp]]]]];

e[1] := {{y'[t] == 1 + y[t]^2, y[0] == 0}, y, 3};
e[2] := {{y'[t] == t^2 + y[t]^2, y[0] == 0}, y, 3};
e[3] := {{y''[t] == y[t]^2, y[0] == 1, y'[0] == 0}, y, 4};
e[4] := {{y'[t] == Exp[y[t]] - y[t], y[0] == 0}, y, 3};
e[5] := {{x'[t] == x[t]^2 + y[t]^2, y'[t] == 2 x[t] y[t], x[0] == 1, y[0] == 1/2}, x, 2};
e[6] := {{y'[t] == y[t]^2 + Exp[t], y[0] == 0}, y, 3};
e[10] := {{y''[t] == Exp[y[t]], y[0] == 0, y'[0] == 0}, y, 4};

pde7[wp_] := dom[u /. First[NDSolve[{D[u[t, x], t] == D[u[t, x], x, x] + u[t, x] + u[t, x]^2,
    u[0, x] == 1, Derivative[0, 1][u][t, 0] == 0, Derivative[0, 1][u][t, 1] == 0},
    u, {t, 0, 2}, {x, 0, 1}, WorkingPrecision -> wp, Evaluate[pde10[wp]]]]];
pde8[wp_] := dom[u /. First[NDSolve[{D[u[t, x, y], t] == D[u[t, x, y], x, x] + D[u[t, x, y], y, y] + Exp[u[t, x, y]],
    u[0, x, y] == 0, Derivative[0, 1, 0][u][t, 0, y] == 0, Derivative[0, 1, 0][u][t, 3, y] == 0,
    Derivative[0, 0, 1][u][t, x, 0] == 0, Derivative[0, 0, 1][u][t, x, 3] == 0},
    u, {t, 0, 2}, {x, 0, 3}, {y, 0, 3}, WorkingPrecision -> wp, Evaluate[pde10[wp]]]]];
cross9[wp_] := Module[{s},
  s = u /. First[NDSolve[{D[u[t, x], t] == D[u[t, x], x, x], u[0, x] == 4 x (1 - x), u[t, 0] == 0, u[t, 1] == 0},
    u, {t, 0, 1}, {x, 0, 1}, WorkingPrecision -> wp, Evaluate[pde10[wp]]]];
  ruleVal[FindRoot[s[tt, 1/2] == 1/2, {tt, 1/100, 1/5}, Method -> "Brent", WorkingPrecision -> wp], tt]];

Do[With[{k = k},
  Switch[k,
    7, vrun["hd15.7 machine", pde7[MachinePrecision], 10, 1, 60];
       vrun["hd15.7 wp50", pde7[60], 48, 1, 90],
    8, vrun["hd15.8 machine", pde8[MachinePrecision], 10, 1, 60];
       vrun["hd15.8 wp50", pde8[60], 48, 1, 90],
    9, vrun["hd15.9 machine", cross9[MachinePrecision], 10, 1, 60];
       vrun["hd15.9 wp50", cross9[60], 48, 1, 90],
    _, vrun["hd15." <> ToString[k] <> " machine", tinf[e[k][[1]], e[k][[2]], e[k][[3]], MachinePrecision], 10, 1, 60];
       vrun["hd15." <> ToString[k] <> " wp50", tinf[e[k][[1]], e[k][[2]], e[k][[3]], 60], 48, 1, 90]]], {k, 10}];
