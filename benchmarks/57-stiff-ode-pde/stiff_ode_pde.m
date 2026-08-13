(* ==========================================================================
   Experiment 57 -- Stiff ODEs and PDEs (group E)
   ==========================================================================
   WHAT IT MEASURES.  NDSolve past the non-stiff scalar IVP of experiment 13:
   the implicit/stiff integrators (src/numerical_calculus/ndsolve_implicit.c,
   ndsolve_adams.c) on a stiff scalar relaxation, the Robertson kinetics (the
   canonical stiff system), a linear ODE system, and the method-of-lines PDE
   path (ndsolve_mol.c, ndsolve_stencil.c) on the heat and wave equations.

   Baseline is scipy.integrate.solve_ivp with an implicit method (BDF/Radau/
   LSODA) for the stiff rows, and a hand-written method-of-lines discretisation
   for the PDE rows.

   Checks are chosen to be integrator-robust: steady-state and analytic values,
   never a point on a fast transient where two solvers legitimately disagree.
   The heat/wave rows check the exact separable solution exp(-pi^2 t) sin(pi x)
   and cos(pi t) sin(pi x) at (0.1, 0.5), which every accurate solver reproduces.
   ========================================================================== *)

Get["../harness.m"];
Get["../data.m"];

require[{"NDSolve"}];

(* ---- stiff / non-stiff ODEs -------------------------------------------- *)

bench["NDSolve stiff scalar relaxation",
      NDSolve[{y'[t] == -1000 y[t] + 3000 - 2000 Exp[-t], y[0] == 0}, y, {t, 0, 1}];];
check["NDSolve stiff scalar relaxation",
      Round[10^4 (y[1] /. First[
        NDSolve[{y'[t] == -1000 y[t] + 3000 - 2000 Exp[-t], y[0] == 0}, y, {t, 0, 1}]])]];

bench["NDSolve Robertson stiff kinetics",
      NDSolve[{a'[t] == -4/100 a[t] + 10^4 b[t] c[t],
               b'[t] == 4/100 a[t] - 10^4 b[t] c[t] - 3 10^7 b[t]^2,
               c'[t] == 3 10^7 b[t]^2, a[0] == 1, b[0] == 0, c[0] == 0},
              {a, b, c}, {t, 0, 1}];];
check["NDSolve Robertson stiff kinetics",
      Round[10^4 (c[1] /. First[
        NDSolve[{a'[t] == -4/100 a[t] + 10^4 b[t] c[t],
                 b'[t] == 4/100 a[t] - 10^4 b[t] c[t] - 3 10^7 b[t]^2,
                 c'[t] == 3 10^7 b[t]^2, a[0] == 1, b[0] == 0, c[0] == 0},
                {a, b, c}, {t, 0, 1}]])]];

bench["NDSolve linear system harmonic",
      NDSolve[{x'[t] == y[t], y'[t] == -x[t], x[0] == 1, y[0] == 0}, {x, y}, {t, 0, 1}];];
check["NDSolve linear system harmonic",
      Round[10^6 (x[1] /. First[
        NDSolve[{x'[t] == y[t], y'[t] == -x[t], x[0] == 1, y[0] == 0}, {x, y}, {t, 0, 1}]])]];

bench["NDSolve Van der Pol mu=10",
      NDSolve[{x''[t] == 10 (1 - x[t]^2) x'[t] - x[t], x[0] == 2, x'[0] == 0},
              x, {t, 0, 1}];];
check["NDSolve Van der Pol mu=10",
      Round[10^3 (x[1] /. First[
        NDSolve[{x''[t] == 10 (1 - x[t]^2) x'[t] - x[t], x[0] == 2, x'[0] == 0},
                x, {t, 0, 1}]])]];

(* ---- PDEs via method of lines ------------------------------------------ *)

bench["NDSolve heat PDE (method of lines)",
      NDSolve[{D[u[t, x], t] == D[u[t, x], {x, 2}], u[0, x] == Sin[Pi x],
               u[t, 0] == 0, u[t, 1] == 0}, u, {t, 0, 0.1}, {x, 0, 1}];];
check["NDSolve heat PDE (method of lines)",
      Round[10^3 (u[0.1, 0.5] /. First[
        NDSolve[{D[u[t, x], t] == D[u[t, x], {x, 2}], u[0, x] == Sin[Pi x],
                 u[t, 0] == 0, u[t, 1] == 0}, u, {t, 0, 0.1}, {x, 0, 1}]])]];

bench["NDSolve wave PDE (method of lines)",
      NDSolve[{D[u[t, x], {t, 2}] == D[u[t, x], {x, 2}], u[0, x] == Sin[Pi x],
               Derivative[1, 0][u][0, x] == 0, u[t, 0] == 0, u[t, 1] == 0},
              u, {t, 0, 0.1}, {x, 0, 1}];];
check["NDSolve wave PDE (method of lines)",
      Round[10^4 (u[0.1, 0.5] /. First[
        NDSolve[{D[u[t, x], {t, 2}] == D[u[t, x], {x, 2}], u[0, x] == Sin[Pi x],
                 Derivative[1, 0][u][0, x] == 0, u[t, 0] == 0, u[t, 1] == 0},
                u, {t, 0, 0.1}, {x, 0, 1}]])]];
