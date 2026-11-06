(* bench-timeout: 1800 *)
(* ==========================================================================
   Experiment 155 -- ISSAC '97 holdout B, family 10: ODE singularities and
   high-precision ODE values (the counterpart of Problem 10, the square-root
   blow-up of y'' = x^3 + y^3 + y'^3)
   ==========================================================================

       cd benchmarks/155-issac97b-ode-singularities
       ../../Mathilda -file ode_singularities.m

   The original has a movable square-root branch point.  These have a
   movable simple pole (a Riccati equation), a double pole, a logarithmic
   singularity, a regular singular starting point (Lane-Emden), a singular
   coefficient plus a boundary condition at infinity (Thomas-Fermi, Blasius),
   a near-separatrix period, chaos (Lorenz), an exponentially unstable
   forward integration (Airy), and repeated near-collisions (Kepler,
   e = 0.99).  True values (closed form where one exists), each established
   twice without either benchmarked system (mpmath odefun is a Taylor-series
   integrator, unrelated to either system's methods):

   b10.1  y' = x^2 + y^2, y(0) = 0: y = x J_{3/4}(x^2/2)/J_{-1/4}(x^2/2), a pole at
          x* = sqrt(2 j), j the first zero of J_{-1/4}: 2.00314735942688470800461097905
          (second: Taylor integration to 1.9, then w = 1/y, w' = -(x^2 w^2 + 1),
          to w = 0).
   b10.2  theta'' = -sin theta, theta(0) = 179 degrees: period
          T = 4 K(sin^2(theta0/2)) = 24.511115298106880440830513773 (second:
          4 x the first zero of theta(t) by Taylor integration, 30 digits agree).
   b10.3  Lorenz (10, 28, 8/3) from (1, 1, 1): x(5) = -6.512113699419599047000097,
          x(20) = 13.79319959512861883160497 (Taylor integration at 30 and 45
          digits; errors grow like e^(0.9 t)).
   b10.4  Lane-Emden n = 3, theta'' + (2/x) theta' + theta^3 = 0, theta(0) = 1:
          first zero xi_1 = 6.89684861937696037545452818712 (series start at
          x0 = 0.001 and 0.002, identical; published 6.896848619376960375,
          Horedt 1986).
   b10.5  Blasius f''' + f f''/2 = 0, f(0) = f'(0) = 0, f'(oo) = 1:
          f''(0) = 0.332057336215196298937180062 (Toepfer scaling: f''(0) =
          g'(oo)^(-3/2) for g''(0) = 1; published Boyd 1999 0.33205733621519630).
   b10.6  y'' = x y from Ai(0), Ai'(0): y(8) = Ai(8) = 4.69220761609923162564908e-8;
          machine instance y(3) = Ai(3) = 6.591139357460719e-3.  The recessive
          solution integrated forward: errors grow like Bi(x)/Ai(x), 2e3 at
          x = 3 and 2.5e13 at x = 8.
   b10.7  y'' = y^2, y(0) = 1, y'(0) = 0: a double pole at
          x* = sqrt(3/2) B(1/6, 1/2)/3 = 2.97447742540217556146948555029
          (energy: y'^2 = 2(y^3 - 1)/3; x* = sqrt(3/2) int_1^oo (y^3-1)^(-1/2) dy,
          y = s^(-1/3) gives the Beta function; second: tanh-sinh quadrature of
          sqrt6/sqrt(1 - v^6) on [0, 1]).
   b10.8  Thomas-Fermi y'' = y^(3/2)/sqrt(x), y(0) = 1, y(oo) = 0: the initial
          slope y'(0) = -1.588071022611375312718684508 (published, Boyd 2013;
          second: bisection shooting in IEEE double (DOP853, rtol 1e-13) on the
          sign of y' + 3y/x at X = 10 / 20 / 40: 1.58809300, 1.58807134,
          1.588071027, converging to the published value like X^-6.2, which
          also sets the method's digits: X = 40 is good to 8, X = 60 to 9).
   b10.9  y'' = e^y, y(0) = y'(0) = 0: y = -2 log cos(x/sqrt2), a logarithmic
          singularity at x* = pi/sqrt2 = 2.22144146907918312350794049503
          (second: u = e^(-y/2) satisfies u'' = -u/2, root of u by Taylor
          integration).
   b10.10 Kepler, mu = a = 1, from pericentre (x, y) = (1 - e, 0), velocity
          (0, sqrt((1+e)/(1-e))): position at t = 4 pi + 1/2 for e = 0.99
          (two pericentre passages at r = 0.01), from Kepler's equation
          E - e sin E = 1/2: (x, y) = (cos E - e, sqrt(1-e^2) sin E)
          = (-0.9057868130076083412421519, 0.1405662547299957358079998).
          Machine instance e = 9/10, t = 2 pi + 1/2:
          (-0.7146936458998492280827225, 0.4283406301936821537237888)
          (second: Taylor integration of the e = 0.9 orbit, 18 digits agree).

   METHODS.  NDSolve at machine precision and at WorkingPrecision 30-50, with
   the reformulations noted per problem (a switch to 1/y past the pole, an
   inverse function x(y) across the log singularity, a first-integral tail
   for the double pole, a series start at the regular singular point,
   Toepfer scaling, bisection shooting); FindRoot on the interpolating
   function for zeros; "naive" = the end of NDSolve's domain at the
   singularity (as in Problem 10).
   CHECKS.  Machine 3-6 digits (naive domain ends 3-4; NDSolve's goals are
   local, so the global error of a default machine solve is ~1e-7), arbitrary
   precision 8-15 (inside the propagated error: chaos and instability noted above).
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"NDSolve", "FindRoot", "InterpolatingFunction", "EllipticK", "AiryAi", "AiryAiPrime"}];

mc[d_][v_] := If[NumericQ[v], sig[v, d], "unevaluated"];
pairChk[d_][v_] := If[MatchQ[v, {_?NumericQ, _?NumericQ}], {sig[v[[1]], d], sig[v[[2]], d]}, "failed"];
(* default goals (half the working precision); NDSolve also gets MaxSteps -> 10^6 *)
wpOpts[wp_] := If[wp === MachinePrecision, {}, {WorkingPrecision -> wp}];
ndOpts[wp_] := Join[{MaxSteps -> 10^6}, wpOpts[wp]];
domEnd[s_, f_] := First[(f /. First[s])["Domain"]][[2]];

(* ---- b10.1 Riccati pole ----------------------------------------------- *)
ric = {y'[x] == x^2 + y[x]^2, y[0] == 0};
ricPole[wp_] := Module[{s1, y0, s2},
  s1 = NDSolve[ric, y, {x, 0, 19/10}, Evaluate[Sequence @@ ndOpts[wp]]];
  y0 = y[19/10] /. First[s1];
  s2 = NDSolve[{w'[x] == -(x^2 w[x]^2 + 1), w[19/10] == 1/y0}, w, {x, 19/10, 21/10}, Evaluate[Sequence @@ ndOpts[wp]]];
  ruleVal[FindRoot[(w /. First[s2])[x], {x, 2}, Evaluate[Sequence @@ wpOpts[wp]]], x]];
run["b10.1 naive domain end machine", domEnd[NDSolve[ric, y, {x, 0, 3}], y], mc[4], 3];
run["b10.1 reciprocal switch wp40", ricPole[40], mc[12], 3];
run["b10.1 reciprocal switch machine", ricPole[MachinePrecision], mc[6], 3];

(* ---- b10.2 pendulum near the separatrix ------------------------------- *)
pend[wp_] := Module[{th0 = 179 Pi/180, s},
  s = NDSolve[{th''[t] == -Sin[th[t]], th[0] == th0, th'[0] == 0}, th, {t, 0, 7}, Evaluate[Sequence @@ ndOpts[wp]]];
  4 ruleVal[FindRoot[(th /. First[s])[t], {t, 6}, Evaluate[Sequence @@ wpOpts[wp]]], t]];
run["b10.2 EllipticK exact", N[4 EllipticK[Sin[179 Pi/360]^2], 30], mc[15], 3];
run["b10.2 NDSolve quarter period wp40", pend[40], mc[15], 3];
run["b10.2 NDSolve quarter period machine", pend[MachinePrecision], mc[6], 3];

(* ---- b10.3 Lorenz ----------------------------------------------------- *)
lorenzX[tEnd_, wp_] := Module[{s},
  s = NDSolve[{a'[t] == 10 (b[t] - a[t]), b'[t] == a[t] (28 - c[t]) - b[t], c'[t] == a[t] b[t] - 8/3 c[t],
               a[0] == 1, b[0] == 1, c[0] == 1}, {a, b, c}, {t, 0, tEnd},
              Evaluate[Sequence @@ ndOpts[wp]]];
  a[tEnd] /. First[s]];
run["b10.3 NDSolve x(20) wp40", lorenzX[20, 40], mc[8]];
run["b10.3 NDSolve x(5) machine", lorenzX[5, MachinePrecision], mc[4], 3];

(* ---- b10.4 Lane-Emden n = 3 (regular singular start) ------------------ *)
laneEmden[wp_] := Module[{x0 = 1/1000, th0, dth0, s},
  th0 = 1 - x0^2/6 + 3 x0^4/120 - 3 19 x0^6/15120 + 3 (122 9 - 183 3 + 70) x0^8/3265920;
  dth0 = D[1 - x^2/6 + 3 x^4/120 - 3 19 x^6/15120 + 3 (122 9 - 183 3 + 70) x^8/3265920, x] /. x -> x0;
  s = NDSolve[{u''[x] == -u[x]^3 - 2 u'[x]/x, u[x0] == th0, u'[x0] == dth0}, u, {x, x0, 7}, Evaluate[Sequence @@ ndOpts[wp]]];
  ruleVal[FindRoot[(u /. First[s])[x], {x, 69/10}, Evaluate[Sequence @@ wpOpts[wp]]], x]];
run["b10.4 series start + FindRoot wp30", laneEmden[30], mc[13], 3];
run["b10.4 series start + FindRoot machine", laneEmden[MachinePrecision], mc[6], 3];

(* ---- b10.5 Blasius by Toepfer scaling --------------------------------- *)
blasius[wp_] := Module[{s},
  s = NDSolve[{g'''[x] == -g[x] g''[x]/2, g[0] == 0, g'[0] == 0, g''[0] == 1}, g, {x, 0, 12}, Evaluate[Sequence @@ ndOpts[wp]]];
  ((g' /. First[s])[12])^(-3/2)];
run["b10.5 Toepfer wp30", blasius[30], mc[13], 3];
run["b10.5 Toepfer machine", blasius[MachinePrecision], mc[6], 3];

(* ---- b10.6 Airy, the recessive solution forward ----------------------- *)
airyFwd[xEnd_, wp_] := Module[{s},
  s = NDSolve[{v''[x] == x v[x], v[0] == AiryAi[0], v'[0] == AiryAiPrime[0]}, v, {x, 0, xEnd}, Evaluate[Sequence @@ ndOpts[wp]]];
  v[xEnd] /. First[s]];
run["b10.6 NDSolve y(8) wp50", airyFwd[8, 50], mc[8], 3];
run["b10.6 NDSolve y(3) machine", airyFwd[3, MachinePrecision], mc[3], 3];

(* ---- b10.7 double pole of y'' = y^2 ----------------------------------- *)
(* x* = X + sqrt(3/2) int_Y^oo (s^3 - 1)^(-1/2) ds, Y = y(X) (energy integral) *)
dpole[wp_] := Module[{s, yy, xx = 297/100},
  s = NDSolve[{q''[x] == q[x]^2, q[0] == 1, q'[0] == 0}, q, {x, 0, xx}, Evaluate[Sequence @@ ndOpts[wp]]];
  yy = q[xx] /. First[s];
  xx + Sqrt[3/2] (2 yy^(-1/2) + yy^(-7/2)/7 + 3 yy^(-13/2)/52)];
run["b10.7 naive domain end machine", domEnd[NDSolve[{q''[x] == q[x]^2, q[0] == 1, q'[0] == 0}, q, {x, 0, 4}], q], mc[3], 3];
run["b10.7 energy tail wp30", dpole[30], mc[13], 3];
run["b10.7 energy tail machine", dpole[MachinePrecision], mc[6], 3];

(* ---- b10.8 Thomas-Fermi by bisection shooting ------------------------- *)
tfShoot[xEnd_, its_, wp_] := Module[{lo = 158/100, hi = 160/100, mid, s, x0 = 1/10^4, sEnd, xe},
  Do[mid = (lo + hi)/2;
     s = Quiet[NDSolve[{p''[x] == Max[p[x], 0]^(3/2)/Sqrt[x], p[x0] == 1 - mid x0 + 4/3 x0^(3/2),
                        p'[x0] == -mid + 2 Sqrt[x0]}, p, {x, x0, xEnd},
                       Evaluate[Sequence @@ ndOpts[wp]]]];
     (* the growing perturbation ~ x^4.77 has y' + 3y/x > 0 (a blow-up: too shallow),
        the 144/x^3 solution ~ 0, a zero crossing < 0; read at the end of the domain *)
     xe = Min[domEnd[s, p], xEnd];
     sEnd = (p'[xe] + 3 p[xe]/xe) /. First[s]; If[sEnd < 0, hi = mid, lo = mid],
     {its}];
  -(lo + hi)/2];
run["b10.8 bisection shooting wp30", tfShoot[60, 34, 30], mc[9], 1, 600];
run["b10.8 bisection shooting machine", tfShoot[40, 30, MachinePrecision], mc[6], 1, 300];

(* ---- b10.9 logarithmic singularity of y'' = e^y ----------------------- *)
(* y(x) to x = 6/5, then x(y): dx/dy = 1/v, dv/dy = e^y/v, to y = Y; tail sqrt2 e^(-Y/2) *)
bratu[yEnd_, wp_] := Module[{s1, y0, v0, s2},
  s1 = NDSolve[{r''[x] == Exp[r[x]], r[0] == 0, r'[0] == 0}, r, {x, 0, 6/5}, Evaluate[Sequence @@ ndOpts[wp]]];
  {y0, v0} = {r[6/5], r'[6/5]} /. First[s1];
  s2 = NDSolve[{xi'[yv] == 1/vv[yv], vv'[yv] == Exp[yv]/vv[yv], xi[y0] == 6/5, vv[y0] == v0},
               {xi, vv}, {yv, y0, yEnd}, Evaluate[Sequence @@ ndOpts[wp]]];
  (xi[yEnd] /. First[s2]) + Sqrt[2] Exp[-yEnd/2]];
run["b10.9 naive domain end machine", domEnd[NDSolve[{r''[x] == Exp[r[x]], r[0] == 0, r'[0] == 0}, r, {x, 0, 3}], r], mc[4], 3];
run["b10.9 inverse function wp30", bratu[80, 30], mc[13], 3];
run["b10.9 inverse function machine", bratu[40, MachinePrecision], mc[6], 3];

(* ---- b10.10 Kepler orbit through pericentre --------------------------- *)
kepler[e_, tEnd_, wp_] := Module[{s},
  s = NDSolve[{kx''[t] == -kx[t]/(kx[t]^2 + ky[t]^2)^(3/2), ky''[t] == -ky[t]/(kx[t]^2 + ky[t]^2)^(3/2),
               kx[0] == 1 - e, ky[0] == 0, kx'[0] == 0, ky'[0] == Sqrt[(1 + e)/(1 - e)]},
              {kx, ky}, {t, 0, tEnd}, Evaluate[Sequence @@ ndOpts[wp]]];
  {kx[tEnd], ky[tEnd]} /. First[s]];
run["b10.10 NDSolve e=0.99 wp30", kepler[99/100, 4 Pi + 1/2, 30], pairChk[9]];
run["b10.10 NDSolve e=0.9 machine", kepler[9/10, 2 Pi + 1/2, MachinePrecision], pairChk[4], 3];
