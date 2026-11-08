(* bench-timeout: 1800 *)
(* ==========================================================================
   Experiment 246 -- SIAM 100-Digit Challenge, sealed holdout B for Appendix D
   Problem 11 (six masses: FindRoot in two unknowns over an NDSolve shot).

   Capability probed: SHOOTING -- a vector root finder (FindRoot with two or
   three unknowns) wrapped around an ODE initial-value solve, with the final
   time itself one of the unknowns.  Ten problems, all started from the same
   rounded guesses in every system and scored on the reported unknown:

     bd11.1  projectile with linear drag  x'' = -x'/2, y'' = -1 - y'/2, from 0
             at speed 2: shortest time to reach (1, 1/2)          T = 0.74865...
     bd11.2  isotropic oscillator r'' = -r from (1,0) at speed 1: shortest
             time to reach (0, 1/2)                               T = 0.91173...
     bd11.3  Kepler r'' = -r/|r|^3 from (1,0) at speed 11/10: shortest time
             to reach (-1/2, 1)                                   T = 1.84358...
     bd11.4  Bratu y'' = -lam e^y, y(0)=0, y'(0)=s: (lam, s) with y(1) = 0 and
             y(1/2) = 1/2; report lam                           lam = 2.63490...
     bd11.5  clamped beam y'''' + 50(1+x) y = 1, y(0)=y'(0)=0: (y''(0), y'''(0))
             with y(1) = y'(1) = 0; report y''(0)                     0.073978...
     bd11.6  pendulum th'' = -sin th, th(0)=0, th'(0)=w: (w, T) with th(T) = 2
             and th'(T) = 0 (first turning point); report T     T = 2.08743...
     bd11.7  Henon-Heiles x'' = -x - 2xy, y'' = -y - x^2 + y^2, from 0 at
             speed 1/2: shortest time to reach (3/10, 1/5)        T = 0.83045...
     bd11.8  two unit masses at (-1,0), (1,0) (inverse-square), from (0,-2) at
             speed 4/5: shortest time to reach (1/2, 3/2)         T = 2.75823...
     bd11.9  3-D oscillator x''=-x, y''=-4y, z''=-9z from 0 at speed 1,
             direction (theta, phi): THREE unknowns (theta, phi, T), shortest
             time to reach (3/10, 1/5, 1/10)                      T = 0.40189...
     bd11.10 Arenstorf orbit (restricted three-body, mu = 0.012277471) from
             (0.994, 0) with velocity (0, v): (v, T) with y(T) = 0, x'(T) = 0
             (perpendicular crossing = half period); report T     T = 8.53260...
             (the harder one: a close lunar flyby at t = 0, long arc)

   machine : NDSolve PrecisionGoal/AccuracyGoal 13, FindRoot defaults.
   wp50    : NDSolve and FindRoot at WorkingPrecision -> 60 (NDSolve
             PrecisionGoal/AccuracyGoal 55), exact data and guesses.
   The guesses are rounded points on the shortest branch (found by a scan in
   the reference); the problems have other, longer solutions.
   True values and provenance: truth.tsv; references siam100_refs/bd11_ref.py.
   ========================================================================== *)
Get["../harness.m"]; Get["../siam100_common.m"];

require[{"NDSolve", "FindRoot"}];

opts[MachinePrecision] = {PrecisionGoal -> 13, AccuracyGoal -> 13, MaxSteps -> Infinity};
opts[60] = {WorkingPrecision -> 60, PrecisionGoal -> 55, AccuracyGoal -> 55, MaxSteps -> Infinity};
fr[MachinePrecision] = {};
fr[60] = {WorkingPrecision -> 60};
g[v_, MachinePrecision] := N[v];
g[v_, 60] := SetPrecision[v, 60];

(* planar point mass: x'' = ax, y'' = ay, start (x0,y0), speed v, heading th *)
posP[ax_, ay_, {x0_, y0_}, v_, th_?NumericQ, T_?NumericQ, wp_] :=
  {x[T], y[T]} /. First[NDSolve[{x''[t] == ax, y''[t] == ay, x[0] == x0, y[0] == y0,
     x'[0] == v Cos[th], y'[0] == v Sin[th]}, {x, y}, {t, 0, T}, Sequence @@ opts[wp]]];

planarT[ax_, ay_, x0_, v_, tgt_, {th0_, T0_}, wp_] :=
  T /. FindRoot[posP[ax, ay, x0, v, th, T, wp] == tgt, {{th, g[th0, wp]}, {T, g[T0, wp]}},
     Evaluate[Sequence @@ fr[wp]]];

twoMass = Sum[(p - {x[t], y[t]})/((p[[1]] - x[t])^2 + (p[[2]] - y[t])^2)^(3/2), {p, {{-1, 0}, {1, 0}}}];

bratu[lam_?NumericQ, s_?NumericQ, wp_] := {w[1], w[1/2]} /. First[NDSolve[
    {w''[t] == -lam Exp[w[t]], w[0] == 0, w'[0] == s}, w, {t, 0, 1}, Sequence @@ opts[wp]]];
beam[p_?NumericQ, q_?NumericQ, wp_] := {w[1], w'[1]} /. First[NDSolve[
    {w''''[t] == 1 - 50 (1 + t) w[t], w[0] == 0, w'[0] == 0, w''[0] == p, w'''[0] == q},
    w, {t, 0, 1}, Sequence @@ opts[wp]]];
pend[om_?NumericQ, T_?NumericQ, wp_] := {w[T], w'[T]} /. First[NDSolve[
    {w''[t] == -Sin[w[t]], w[0] == 0, w'[0] == om}, w, {t, 0, T}, Sequence @@ opts[wp]]];
osc3[th_?NumericQ, ph_?NumericQ, T_?NumericQ, wp_] := {x[T], y[T], z[T]} /. First[NDSolve[
    {x''[t] == -x[t], y''[t] == -4 y[t], z''[t] == -9 z[t], x[0] == 0, y[0] == 0, z[0] == 0,
     x'[0] == Cos[ph] Cos[th], y'[0] == Cos[ph] Sin[th], z'[0] == Sin[ph]},
    {x, y, z}, {t, 0, T}, Sequence @@ opts[wp]]];
mu = 12277471/10^9;
aren[v_?NumericQ, T_?NumericQ, wp_] := {y[T], x'[T]} /. First[NDSolve[
    {x''[t] == x[t] + 2 y'[t] - (1 - mu) (x[t] + mu)/((x[t] + mu)^2 + y[t]^2)^(3/2) -
                mu (x[t] - 1 + mu)/((x[t] - 1 + mu)^2 + y[t]^2)^(3/2),
     y''[t] == y[t] - 2 x'[t] - (1 - mu) y[t]/((x[t] + mu)^2 + y[t]^2)^(3/2) -
                mu y[t]/((x[t] - 1 + mu)^2 + y[t]^2)^(3/2),
     x[0] == 994/1000, y[0] == 0, x'[0] == 0, y'[0] == v}, {x, y}, {t, 0, T}, Sequence @@ opts[wp]]];

cases[wp_, tag_, d_] := (
  vrun["bd11.1 " <> tag, planarT[-x'[t]/2, -1 - y'[t]/2, {0, 0}, 2, {1, 1/2}, {0.45, 0.6}, wp], d];
  vrun["bd11.2 " <> tag, planarT[-x[t], -y[t], {1, 0}, 1, {0, 1/2}, {2.47, 0.91}, wp], d];
  vrun["bd11.3 " <> tag, planarT[-x[t]/(x[t]^2 + y[t]^2)^(3/2), -y[t]/(x[t]^2 + y[t]^2)^(3/2),
                                 {1, 0}, 11/10, {-1/2, 1}, {1.71, 1.84}, wp], d];
  vrun["bd11.4 " <> tag, lam /. FindRoot[bratu[lam, s, wp] == {0, 1/2},
                                 {{lam, g[2.7, wp]}, {s, g[1.9, wp]}}, Evaluate[Sequence @@ fr[wp]]], d];
  vrun["bd11.5 " <> tag, p /. FindRoot[beam[p, q, wp] == {0, 0},
                                 {{p, g[0.08, wp]}, {q, g[-0.5, wp]}}, Evaluate[Sequence @@ fr[wp]]], d];
  vrun["bd11.6 " <> tag, T /. FindRoot[pend[om, T, wp] == {2, 0},
                                 {{om, g[1.7, wp]}, {T, g[2.1, wp]}}, Evaluate[Sequence @@ fr[wp]]], d];
  vrun["bd11.7 " <> tag, planarT[-x[t] - 2 x[t] y[t], -y[t] - x[t]^2 + y[t]^2, {0, 0}, 1/2,
                                 {3/10, 1/5}, {0.5836, 0.8305}, wp], d];
  vrun["bd11.8 " <> tag, planarT[twoMass[[1]], twoMass[[2]], {0, -2}, 4/5, {1/2, 3/2},
                                 {1.50, 2.75}, wp], d];
  vrun["bd11.9 " <> tag, T /. FindRoot[osc3[th, ph, T, wp] == {3/10, 1/5, 1/10},
                                 {{th, g[0.6, wp]}, {ph, g[0.25, wp]}, {T, g[0.38, wp]}},
                                 Evaluate[Sequence @@ fr[wp]]], d];
  vrun["bd11.10 " <> tag, T /. FindRoot[aren[v, T, wp] == {0, 0},
                                 {{v, g[-2.0015851063790825, wp]}, {T, g[8.532608280078981, wp]}},
                                 Evaluate[Sequence @@ fr[wp]]], d]);

cases[MachinePrecision, "machine", 10];
cases[60, "wp50", 48];
