(* bench-timeout: 1800 *)
(* ==========================================================================
   Experiment 250 -- SIAM 100-Digit Challenge, sealed holdout B for Appendix D
   Problem 15 (blow-up time of u_t = Lap u + e^u on the 3x3 square).

   Capability probed: integrate a reaction-diffusion PDE (method of lines) or
   an ODE until the solution becomes singular, and report where NDSolve stops
   ("step size is effectively zero; singularity suspected") -- the end of the
   solution's domain.  Every problem is built to have an exact blow-up time
   (constant Neumann data -> an ODE; the substitution w = e^-u or w = 1/u ->
   the linear w_t = Lap w - 1, whose first zero is the blow-up; Bessel /
   Beta-function ODE blow-ups), so the truth has >= 55 digits:

     bd15.1  u_t = u_xx + u^2 + u on [0,1], u_x = 0 at both ends, u(0,x) = 1
                                                          t* = log 2 = 0.69314...
     bd15.2  u_t = u_xx - u_x^2 + e^u on [0,Pi], u = 0 on the boundary and at
             t = 0 (the original's structure in 1-D, plus a gradient term)
                                                                t* = 1.69527...
     bd15.3  u_t = Lap u - |grad u|^2 + e^u on [0,Pi]^2, Neumann, u(0) =
             -log(2 + cos x cos y) (blows up at the corners (0,Pi), (Pi,0))
                                               t* = 2 + W(-2e^-4)/2 = 1.98097...
     bd15.4  u_t = u_xx - 2 u_x^2/u + u^2 on [0,Pi], Neumann, u(0) = 1/(2+cos x)
             (blows up at x = Pi)                    t* = 2 + W(-e^-2) = 1.84140...
     bd15.5  QUENCHING u_t = u_xx + e^u/(1-u) on [0,1], Neumann, u(0) = 0
             (u reaches 1 with u_t -> oo)                     t* = 1/e = 0.36787...
     bd15.6  ODE y' = y^2 + t^2, y(0) = 0         t* = sqrt(2 j_{-1/4,1}) = 2.00314...
     bd15.7  ODE y'' = 6 y^2, y(0) = 1, y'(0) = 0 (double pole)
                                                      t* = B(1/6,1/2)/6 = 1.21432...
     bd15.8  SYSTEM u_t = u_xx + v^2, v_t = v_xx + u^2 on [0,1], Neumann,
             u(0) = 1, v(0) = 1/2                 t* = 2F1(1/3,2/3;4/3;7/8) = 1.31269...
     bd15.9  u_t = Lap u - |grad u|^2 + e^u on the DISK of radius 3, u = 0 on the
             boundary and at t = 0 (FEM in space)                  t* = 1.06374...
     bd15.10 u_t = Lap u - |grad u|^2 + e^u on [0,4]^2, u = 0 on the boundary
             and at t = 0 (the original's geometry, larger square)  t* = 1.61698...

   machine : NDSolve default options (no grid or tolerance tuning).
   wp50    : WorkingPrecision -> 60.  A fixed spatial grid limits every PDE
             row long before 50 digits (the ODE rows 6, 7 are attainable); the
             tier is written so that its cost / digits are recorded.
   True values and provenance: truth.tsv; references siam100_refs/bd15_ref.py.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NDSolve"}];

dom[sol_] := sol["Domain"][[1, 2]];
wpo[MachinePrecision] = {};
wpo[60] = {WorkingPrecision -> 60};
nx[f_] := Derivative[0, 1][f];
nxy[f_, i_] := (Derivative @@ ReplacePart[{0, 0, 0}, i -> 1])[f];

blow[1, wp_] := dom[u /. First[NDSolve[{D[u[t, x], t] == D[u[t, x], x, x] + u[t, x]^2 + u[t, x],
    u[0, x] == 1, nx[u][t, 0] == 0, nx[u][t, 1] == 0}, u, {t, 0, 2}, {x, 0, 1}, Sequence @@ wpo[wp]]]];
blow[2, wp_] := dom[u /. First[NDSolve[{D[u[t, x], t] == D[u[t, x], x, x] - D[u[t, x], x]^2 + Exp[u[t, x]],
    u[0, x] == 0, u[t, 0] == 0, u[t, Pi] == 0}, u, {t, 0, 3}, {x, 0, Pi}, Sequence @@ wpo[wp]]]];
blow[3, wp_] := dom[u /. First[NDSolve[{D[u[t, x, y], t] == D[u[t, x, y], x, x] + D[u[t, x, y], y, y] -
      D[u[t, x, y], x]^2 - D[u[t, x, y], y]^2 + Exp[u[t, x, y]],
    u[0, x, y] == -Log[2 + Cos[x] Cos[y]],
    nxy[u, 2][t, 0, y] == 0, nxy[u, 2][t, Pi, y] == 0, nxy[u, 3][t, x, 0] == 0, nxy[u, 3][t, x, Pi] == 0},
    u, {t, 0, 3}, {x, 0, Pi}, {y, 0, Pi}, Sequence @@ wpo[wp]]]];
blow[4, wp_] := dom[u /. First[NDSolve[{D[u[t, x], t] == D[u[t, x], x, x] - 2 D[u[t, x], x]^2/u[t, x] + u[t, x]^2,
    u[0, x] == 1/(2 + Cos[x]), nx[u][t, 0] == 0, nx[u][t, Pi] == 0}, u, {t, 0, 3}, {x, 0, Pi},
    Sequence @@ wpo[wp]]]];
blow[5, wp_] := dom[u /. First[NDSolve[{D[u[t, x], t] == D[u[t, x], x, x] + Exp[u[t, x]]/(1 - u[t, x]),
    u[0, x] == 0, nx[u][t, 0] == 0, nx[u][t, 1] == 0}, u, {t, 0, 1}, {x, 0, 1}, Sequence @@ wpo[wp]]]];
blow[6, wp_] := dom[y /. First[NDSolve[{y'[t] == y[t]^2 + t^2, y[0] == 0}, y, {t, 0, 3}, Sequence @@ wpo[wp]]]];
blow[7, wp_] := dom[y /. First[NDSolve[{y''[t] == 6 y[t]^2, y[0] == 1, y'[0] == 0}, y, {t, 0, 2},
    Sequence @@ wpo[wp]]]];
blow[8, wp_] := dom[u /. First[NDSolve[{D[u[t, x], t] == D[u[t, x], x, x] + v[t, x]^2,
    D[v[t, x], t] == D[v[t, x], x, x] + u[t, x]^2, u[0, x] == 1, v[0, x] == 1/2,
    nx[u][t, 0] == 0, nx[u][t, 1] == 0, nx[v][t, 0] == 0, nx[v][t, 1] == 0},
    {u, v}, {t, 0, 2}, {x, 0, 1}, Sequence @@ wpo[wp]]]];
blow[9, wp_] := dom[u /. First[NDSolve[{D[u[t, x, y], t] == Laplacian[u[t, x, y], {x, y}] -
      D[u[t, x, y], x]^2 - D[u[t, x, y], y]^2 + Exp[u[t, x, y]],
    u[0, x, y] == 0, DirichletCondition[u[t, x, y] == 0, True]},
    u, {t, 0, 2}, {x, y} \[Element] Disk[{0, 0}, 3], Sequence @@ wpo[wp]]]];
blow[10, wp_] := dom[u /. First[NDSolve[{D[u[t, x, y], t] == D[u[t, x, y], x, x] + D[u[t, x, y], y, y] -
      D[u[t, x, y], x]^2 - D[u[t, x, y], y]^2 + Exp[u[t, x, y]],
    u[0, x, y] == 0, u[t, 0, y] == 0, u[t, 4, y] == 0, u[t, x, 0] == 0, u[t, x, 4] == 0},
    u, {t, 0, 3}, {x, 0, 4}, {y, 0, 4}, Sequence @@ wpo[wp]]]];

Do[vrun["bd15." <> ToString[k] <> " machine", blow[k, MachinePrecision], 10], {k, 1, 10}];
Do[vrun["bd15." <> ToString[k] <> " wp50", blow[k, 60], 48], {k, 1, 10}];
