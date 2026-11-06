(* bench-timeout: 900 *)
(* ==========================================================================
   Experiment 145 -- ISSAC '97 holdout, family 10: ODE singularities and targets
   ==========================================================================

       cd benchmarks/145-issac97h-ode-singularities
       ../../Mathilda -file ode_singularities.m

   HOLDOUT.  Ten problems probing what experiment 135 (the blow-up point of
   y'' = x^3 + y^3 + y'^3 and the value of y there) probes -- locating a
   movable singularity, starting at or reaching a singular point, and
   hitting a boundary-value target, with NDSolve plus FindRoot on the
   solution -- on different equations.  See benchmarks/ISSAC97_HOLDOUT.md.

   Each problem: a machine-precision case (default NDSolve, checked to 6
   digits; 5 for the naive domain end) and a WorkingPrecision -> 30 case
   (PrecisionGoal and AccuracyGoal 20, checked to 15 digits).

   PROBLEMS, METHODS AND TRUE VALUES (the true value by a closed form and/or
   a 35-40 digit Taylor-series ODE integration in mpmath, odefun; where both
   exist they agree to all digits shown):
   h10.1  y' = x^2 + y^2, y(0) = 0: blow-up point
          = Sqrt[2 j], j the first zero of BesselJ[-1/4, t] (y = x J_{3/4}(x^2/2)/J_{-1/4}(x^2/2))
          = 2.003147359426884708004610979                    (closed form, odefun)
          method: integrate y to x = 1, continue with u = 1/y (u' = -(x^2 u^2 + 1),
          regular through the pole), FindRoot u = 0.
   h10.2  y' = 1 + y^3, y(0) = 0: blow-up point = Int_0^oo dy/(1 + y^3)
          = 2 Pi/(3 Sqrt[3]) = 1.209199576156145233729385505 (closed form, quadrature)
          method: x as a function of s = y/(1 + y) in [0, 1),
          x'(s) = (1 - s)/((1 - s)^3 + s^3), x(1) by NDSolve.
   h10.3  Lane-Emden n = 3, y'' + 2 y'/x + y^3 = 0, y(0) = 1, y'(0) = 0
          (regular singular point at 0): first zero xi_1
          = 6.896848619376960375454528187  (published, Horedt; odefun from the
          80-term power series at x = 1/2)
          method: start at x = 10^-6 from the series 1 - x^2/6 + n x^4/120,
          NDSolve to 7, FindRoot y = 0.
   h10.4  Lane-Emden n = 3/2 (y^(3/2) written Abs[y]^(3/2)): xi_1
          = 3.653753736219122424609094280  (published, Horedt; odefun)
   h10.5  Bratu y'' + E^y = 0, y(0) = y(1) = 0, lower branch: y'(0)
          = theta Tanh[theta/4], theta = Sqrt[2] Cosh[theta/4] (lower root)
          = 0.5493527287752708190186831942 (closed form; odefun residual 1e-42)
          method: shooting, FindRoot on y(1; s) = 0 by the secant method.
   h10.6  Blasius f''' + f f''/2 = 0, f(0) = f'(0) = 0, f'(oo) = 1: f''(0)
          = 0.3320573362151962989371800620 (published, Boyd 1999; odefun with
          the scaling g''(0) = 1, kappa = g'(L)^(-3/2), L = 10 and 14 agree to 19 digits)
          machine: shooting on f'(10) = 1 (secant); wp30: the scaling method at L = 12.
   h10.7  Gauss x(1-x) y'' + (c - (a+b+1) x) y' - a b y = 0, a = 1/3, b = 1/2,
          c = 2, y(0) = 1: the value at the singular point x = 1
          = Gamma[c] Gamma[c-a-b]/(Gamma[c-a] Gamma[c-b])
          = 1.159595266963928365769992052 (Gauss's theorem; Hypergeometric2F1[a, b, c, 1])
          method: start at x = 10^-6 from the series 1 + (a b/c) x, NDSolve to 1.
   h10.8  y' = x + y^2, y(0) = 0: blow-up point = first zero of
          Bi'(0) Ai(-x) - Ai'(0) Bi(-x) = 1.986352707430472813471818390
          (closed form, odefun) -- machine: the naive NDSolve's domain end;
          wp30: the reciprocal continuation of h10.1.
   h10.9  y'' = x + y^2 + y'^3, y(0) = y'(0) = 0: y' blows up at r with y(r)
          finite (y'' ~ y'^3, a square-root branch point).  method: y(x) to
          x = 3/2, then the inverse function X(t), t = y, which is regular at
          the branch point: X'' = -(X + t^2) X'^3 - 1, and r = X(ts) where
          X'(ts) = 0.  True values (odefun, switching at x = 3/2 and at 7/5,
          agreeing to all digits): r = 1.617137534968669344001011396,
          y(r) = 1.149567134607272215450543590.
   h10.10 Bessel x y'' + y' + x y = 0, y(0) = 1, y'(0) = 0 (regular singular
          point at 0): y(10) = BesselJ[0, 10] = -0.2459357644513483351977608625
          (closed form) -- start at x = 10^-6 from 1 - x^2/4.
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"NDSolve", "FindRoot", "InterpolatingFunction"}];

mp = MachinePrecision;
opts[wp_] := If[wp === mp, {MaxSteps -> 100000},
  {WorkingPrecision -> wp, PrecisionGoal -> 20, AccuracyGoal -> 20, MaxSteps -> 100000}];
fropts[wp_] := If[wp === mp, {}, {WorkingPrecision -> wp}];
num[v_] := If[NumericQ[v], v, "failed"];

(* ---- h10.1 / h10.8: Riccati blow-up through the reciprocal ------------- *)
recip[rhs_, urhs_, wp_] := Module[{s1, y1, s2},
  s1 = NDSolve[{y'[x] == rhs[x, y[x]], y[0] == 0}, y, {x, 0, 1}, Evaluate[Sequence @@ opts[wp]]];
  y1 = y[1] /. First[s1];
  s2 = NDSolve[{u'[x] == urhs[x, u[x]], u[1] == 1/y1}, u, {x, 1, 23/10}, Evaluate[Sequence @@ opts[wp]]];
  num[ruleVal[FindRoot[Evaluate[u[x] /. First[s2]], {x, 2}, Evaluate[Sequence @@ fropts[wp]]], x]]];
run["h10.1 reciprocal machine", recip[#1^2 + #2^2 &, -(#1^2 #2^2 + 1) &, mp], sig[#, 6] &, 3];
run["h10.1 reciprocal wp30",    recip[#1^2 + #2^2 &, -(#1^2 #2^2 + 1) &, 30], sig[#, 15] &];

run["h10.8 naive NDSolve domain end machine",
    First[(y /. First[NDSolve[{y'[x] == x + y[x]^2, y[0] == 0}, y, {x, 0, 3}]])["Domain"]][[2]],
    sig[#, 5] &, 3];
run["h10.8 reciprocal wp30", recip[#1 + #2^2 &, -(#1 #2^2 + 1) &, 30], sig[#, 15] &];

(* ---- h10.2: blow-up as x(s) at s = 1 ------------------------------------ *)
bu2[wp_] := Module[{s = NDSolve[{xs'[t] == (1 - t)/((1 - t)^3 + t^3), xs[0] == 0}, xs, {t, 0, 1},
     Evaluate[Sequence @@ opts[wp]]]}, num[xs[1] /. First[s]]];
run["h10.2 substitution machine", bu2[mp], sig[#, 6] &, 3];
run["h10.2 substitution wp30",    bu2[30], sig[#, 15] &];

(* ---- h10.3 / h10.4: Lane-Emden first zero from a series start ---------- *)
laneEmden[n_, xend_, guess_, wp_] := Module[{x0 = 10^-6, s},
  s = NDSolve[{y''[x] + 2 y'[x]/x + Abs[y[x]]^n == 0,
               y[x0] == 1 - x0^2/6 + n x0^4/120, y'[x0] == -x0/3 + n x0^3/30},
              y, {x, x0, xend}, Evaluate[Sequence @@ opts[wp]]];
  num[ruleVal[FindRoot[Evaluate[y[x] /. First[s]], {x, guess}, Evaluate[Sequence @@ fropts[wp]]], x]]];
run["h10.3 Lane-Emden n=3 machine", laneEmden[3, 7, 69/10, mp], sig[#, 6] &, 3];
run["h10.3 Lane-Emden n=3 wp30",    laneEmden[3, 7, 69/10, 30], sig[#, 15] &];
run["h10.4 Lane-Emden n=3/2 machine", laneEmden[3/2, 37/10, 365/100, mp], sig[#, 6] &, 3];
run["h10.4 Lane-Emden n=3/2 wp30",    laneEmden[3/2, 37/10, 365/100, 30], sig[#, 15] &];

(* ---- h10.5: Bratu shooting ---------------------------------------------- *)
bratuEnd[s_?NumericQ, wp_] := y[1] /. First[NDSolve[{y''[x] + Exp[y[x]] == 0, y[0] == 0, y'[0] == s},
    y, {x, 0, 1}, Evaluate[Sequence @@ opts[wp]]]];
bratu[wp_] := num[ruleVal[FindRoot[bratuEnd[s, wp], {s, 1/2, 6/10}, Evaluate[Sequence @@ fropts[wp]]], s]];
run["h10.5 Bratu shooting machine", bratu[mp], sig[#, 6] &, 3];
run["h10.5 Bratu shooting wp30",    bratu[30], sig[#, 15] &];

(* ---- h10.6: Blasius ------------------------------------------------------ *)
blasiusEnd[s_?NumericQ] := f'[10] /. First[NDSolve[{f'''[x] + f[x] f''[x]/2 == 0, f[0] == 0, f'[0] == 0,
    f''[0] == s}, f, {x, 0, 10}, MaxSteps -> 100000]];
run["h10.6 Blasius shooting machine", num[ruleVal[FindRoot[blasiusEnd[s] == 1, {s, 3/10, 35/100}], s]],
    sig[#, 6] &, 3];
blasiusScaled[wp_] := Module[{sol = NDSolve[{g'''[x] + g[x] g''[x]/2 == 0, g[0] == 0, g'[0] == 0, g''[0] == 1},
     g, {x, 0, 12}, Evaluate[Sequence @@ opts[wp]]]}, num[(g'[12] /. First[sol])^(-3/2)]];
run["h10.6 Blasius scaling wp30", blasiusScaled[30], sig[#, 15] &];

(* ---- h10.7: value at the singular point of the Gauss equation ---------- *)
gauss[wp_] := Module[{a = 1/3, b = 1/2, c = 2, x0 = 10^-6, s},
  s = NDSolve[{x (1 - x) y''[x] + (c - (a + b + 1) x) y'[x] - a b y[x] == 0,
               y[x0] == 1 + a b x0/c, y'[x0] == a b/c + a (a + 1) b (b + 1) x0/(2 c (c + 1))},
              y, {x, x0, 1}, Evaluate[Sequence @@ opts[wp]]];
  num[y[1] /. First[s]]];
run["h10.7 Gauss value at x=1 machine", gauss[mp], sig[#, 6] &, 3];
run["h10.7 Gauss value at x=1 wp30",    gauss[30], sig[#, 15] &];

(* ---- h10.9: square-root branch point through the inverse function ------ *)
branch[wp_] := Module[{s1, y0, z0, s2, ts},
  s1 = NDSolve[{y'[x] == z[x], z'[x] == x + y[x]^2 + z[x]^3, y[0] == 0, z[0] == 0}, {y, z}, {x, 0, 3/2},
               Evaluate[Sequence @@ opts[wp]]];
  {y0, z0} = {y[3/2], z[3/2]} /. First[s1];
  s2 = NDSolve[{xx''[t] == -(xx[t] + t^2) xx'[t]^3 - 1, xx[y0] == 3/2, xx'[y0] == 1/z0},
               xx, {t, y0, y0 + 1}, Evaluate[Sequence @@ opts[wp]]];
  ts = ruleVal[FindRoot[Evaluate[xx'[t] /. First[s2]], {t, y0 + 3/10}, Evaluate[Sequence @@ fropts[wp]]], t];
  If[NumericQ[ts], {xx[ts] /. First[s2], ts}, "failed"]];
pairChk[d_][v_] := If[MatchQ[v, {_?NumericQ, _?NumericQ}], {sig[v[[1]], d], sig[v[[2]], d]}, "failed"];
run["h10.9 inverse function machine", branch[mp], pairChk[6], 3];
run["h10.9 inverse function wp30",    branch[30], pairChk[15]];

(* ---- h10.10: Bessel from the regular singular point --------------------- *)
bessel0[wp_] := Module[{x0 = 10^-6, s},
  s = NDSolve[{x y''[x] + y'[x] + x y[x] == 0, y[x0] == 1 - x0^2/4, y'[x0] == -x0/2 + x0^3/16},
              y, {x, x0, 10}, Evaluate[Sequence @@ opts[wp]]];
  num[y[10] /. First[s]]];
run["h10.10 Bessel y(10) machine", bessel0[mp], sig[#, 6] &, 3];
run["h10.10 Bessel y(10) wp30",    bessel0[30], sig[#, 15] &];
