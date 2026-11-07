(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 218 -- SIAM 100-Digit Challenge, holdout A, family d13
   (generalises App. D Problem 13, experiment 188: "if |u(50)| = 1, what is
   z0?").

   The capability: SHOOTING ON INITIAL DATA -- find the initial value (or
   initial slope) for which the solution of an IVP meets a condition at a
   fixed final time.  Ten instances that vary the kind of ODE (scalar /
   system, first / second / third order, linear nonnormal, near-blow-up,
   conservative, a two-point BVP posed as shooting, a needle threshold):

     hd13.1  y' = y(1-y), y(0) = y0;  y(10) = 1/2                         (logistic)
     hd13.2  u' = A u, A = {{-1,10,0},{0,-1,10},{0,0,-1}}, u(0) = (0,0,z0);
             |u(5)| = 1                                                   (exp 188's linear part)
     hd13.3  y' = y^2 - y, y(0) = y0;  y(3) = 2                           (near blow-up)
     hd13.4  th'' = -sin th, th(0) = 0, th'(0) = v0;  th(2) = 1           (pendulum)
     hd13.5  y'' = -e^y, y(0) = 0, y'(0) = s;  y(1) = 0, smaller s        (Bratu by shooting)
     hd13.6  f''' + f f''/2 = 0, f(0) = f'(0) = 0, f''(0) = s;  f'(10) = 1  (Blasius)
     hd13.7  r'' = -r/|r|^3, r(0) = (1,0), r'(0) = (0,v0);  |r(2)| = 3/2   (Kepler)
     hd13.8  y' = t^2 + y^2, y(0) = y0;  y(1) = 2                         (Riccati)
     hd13.9  exp 188's ODE, u(0) = (0,0,z0);  |u(30)| = 1                (the original, T = 30)
     hd13.10 x'' + x + x^3 = 0, x(0) = x0, x'(0) = 0;  x(10) = 0, smallest x0 > 0
                                                                          (Duffing, 7th quarter period)
   Each bracket below holds exactly one root (checked); closed forms exist for
   1-5, 7, 8, 10 (truth.tsv, siam100_refs/hd13_ref.py).

   Asked of each tier:
     machine: FindRoot (Brent, the bracket shown) on the initial value of
              (state at T) - target, NDSolve PrecisionGoal/AccuracyGoal 12.
              hd13.9 is stated as exp 188 states it: the time at which |u|^2
              first reaches 1 (WhenEvent; never -> 100) equals 30 (|u| grows
              monotonically once the orbit leaves the saddle).
     wp50:    the same at WorkingPrecision -> 60, PrecisionGoal/AccuracyGoal 50
              for NDSolve and FindRoot.
   True values and provenance: truth.tsv.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NDSolve", "WhenEvent", "FindRoot"}];

(* the state at T from NDSolve: fin[eqs, {f1, ...}, T, wp] = {f1[T], ...} *)
fin[eqs_, fs_List, T_, wp_] := Through[(fs /. First[NDSolve[eqs, fs, {t, 0, T}, opts[wp]]])[T]];

opts[MachinePrecision] = Sequence[PrecisionGoal -> 12, AccuracyGoal -> 12];
opts[wp_] := Sequence[WorkingPrecision -> wp, PrecisionGoal -> 50, AccuracyGoal -> 50, MaxSteps -> Infinity];

g[1][p_?NumericQ, wp_] := First[fin[{y'[t] == y[t] (1 - y[t]), y[0] == p}, {y}, 10, wp]] - 1/2;
g[2][p_?NumericQ, wp_] := Norm[fin[{x'[t] == -x[t] + 10 y[t], y'[t] == -y[t] + 10 z[t], z'[t] == -z[t],
    x[0] == 0, y[0] == 0, z[0] == p}, {x, y, z}, 5, wp]] - 1;
g[3][p_?NumericQ, wp_] := First[fin[{y'[t] == y[t]^2 - y[t], y[0] == p}, {y}, 3, wp]] - 2;
g[4][p_?NumericQ, wp_] := First[fin[{th''[t] == -Sin[th[t]], th[0] == 0, th'[0] == p}, {th}, 2, wp]] - 1;
g[5][p_?NumericQ, wp_] := First[fin[{y''[t] == -Exp[y[t]], y[0] == 0, y'[0] == p}, {y}, 1, wp]];
g[6][p_?NumericQ, wp_] := (f' /. First[NDSolve[{f'''[t] + f[t] f''[t]/2 == 0, f[0] == 0, f'[0] == 0, f''[0] == p},
    f, {t, 0, 10}, opts[wp]]])[10] - 1;
g[7][p_?NumericQ, wp_] := Norm[fin[{x''[t] == -x[t]/(x[t]^2 + y[t]^2)^(3/2), y''[t] == -y[t]/(x[t]^2 + y[t]^2)^(3/2),
    x[0] == 1, y[0] == 0, x'[0] == 0, y'[0] == p}, {x, y}, 2, wp]] - 3/2;
g[8][p_?NumericQ, wp_] := First[fin[{y'[t] == t^2 + y[t]^2, y[0] == p}, {y}, 1, wp]] - 2;
g[9][p_?NumericQ, wp_] := (x /. First[NDSolve[{
    x'[t] == -x[t] + 10 y[t] + Sqrt[x[t]^2 + y[t]^2 + z[t]^2] (-7/10 y[t] - 3/100 z[t]),
    y'[t] == -y[t] + 10 z[t] + Sqrt[x[t]^2 + y[t]^2 + z[t]^2] (7/10 y[t] - 1/10 z[t]),
    z'[t] == -z[t] + Sqrt[x[t]^2 + y[t]^2 + z[t]^2] (3/100 x[t] + 1/10 y[t]),
    x[0] == 0, y[0] == 0, z[0] == p,
    WhenEvent[x[t]^2 + y[t]^2 + z[t]^2 == 1, "StopIntegration"]}, x, {t, 0, 100}, opts[wp]]])["Domain"][[1, 2]] - 30;
g[10][p_?NumericQ, wp_] := First[fin[{x''[t] + x[t] + x[t]^3 == 0, x[0] == p, x'[0] == 0}, {x}, 10, wp]];

brk = {{1/100000, 1/10000}, {1/20, 1/2}, {1, 104/100}, {1/2, 3/2}, {1/5, 1},
       {1/5, 1/2}, {11/10, 13/10}, {0, 9/10}, {64/10000, 65/10000}, {3/10, 6/10}};

Do[With[{k = k, lo = brk[[k, 1]], hi = brk[[k, 2]]},
  vrun["hd13." <> ToString[k] <> " machine",
       ruleVal[FindRoot[g[k][p, MachinePrecision] == 0, {p, N[lo], N[hi]}, Method -> "Brent"], p], 10, 1, 60];
  vrun["hd13." <> ToString[k] <> " wp50",
       ruleVal[FindRoot[g[k][p, 60] == 0, {p, lo, hi}, Method -> "Brent", WorkingPrecision -> 60,
                        AccuracyGoal -> 50, PrecisionGoal -> 50], p], 48, 1, 90]], {k, 10}];
