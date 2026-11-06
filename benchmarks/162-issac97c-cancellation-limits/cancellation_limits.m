(* bench-timeout: 1500 *)
(* ==========================================================================
   Experiment 162 -- ISSAC '97 holdout C (SEALED), family c7:
   limits with deep cancellation (counterpart of Problem 7)
   ==========================================================================

       cd benchmarks/162-issac97c-cancellation-limits
       ../../Mathilda -file cancellation_limits.m

   Ten limits in which the leading terms cancel and the answer lives in a
   second- or third-order term: a Gamma ratio past its 1/(8x) term, the
   Laurent expansion of Zeta at 1 past its pole and constant (the Stieltjes
   constant gamma_1), an iterated logarithm (whose correction decays only like
   1/Log[x]), a nested radical, Erf o Sin against Sin o Erf, the harmonic
   number past the 1/(2n) term (at non-integer n), a power tower at infinity
   and at 0+, Gamma at its pole past the -EulerGamma term, and the
   dilogarithm at its branch point.

   THREE ROWS PER LIMIT.
     exact    Limit[f, x -> x0], checked to 20 digits of the exact answer.
     machine  N[f /. x -> P1] with P1 chosen so that machine arithmetic still
              carries the digits checked: the leading terms cancel to at most
              ~1e7, so a double keeps >= 8 digits.  This is a REDUCED instance
              (at the wp point a double returns noise: -1e40, 0., 5.8e19).
     wp30     N[f /. x -> P2, 30] at a point P2 where the cancellation is 20
              to 60 digits deep (adaptive-precision evaluation of an exact
              expression).  Checked to 20 digits of the value AT P2, which is
              not the limit: at P2 the value still differs from it, in the 2nd
              digit (c7.3, whose correction is 1/(2 Log[x])) down to the 27th
              (c7.7, c7.9).

   TRUE ANSWERS.  Limits: A = the asymptotic expansion derived by hand (below)
   or sympy series (c7.5); B = mpmath at 300 digits at a point 30-60 orders
   closer to x0.  Point values: A = mpmath at 150 digits, B = mpmath at 250
   digits (agree to >= 1e-80).  Digits: machine 5-7 (5 for c7.2, whose Zeta[1.001] loses 7 digits),
   wp30 20, exact 20.

   c7.1  x^2 (Gamma[x+1/2]/(Sqrt[x] Gamma[x]) - 1 + 1/(8x)), x -> oo.
         Gamma[x+1/2]/Gamma[x] = Sqrt[x] (1 - 1/(8x) + 1/(128 x^2) + ...)
         -> 1/128.  P1 = 100+1/3: 0.00786110073873594954270274898838;
         P2 = 10^20+1/3: 0.007812500000000000000048828125.
   c7.2  (Zeta[s] - 1/(s-1) - EulerGamma)/(s-1), s -> 1: Zeta[s] = 1/(s-1)
         + gamma - gamma_1 (s-1) + ... -> -StieltjesGamma[1]
         = 0.072815845483676724860586375874901.  P1 = 1 + 10^-3:
         0.0728109999598714357921509401701; P2 = 1 + 10^-20:
         0.0728158454836767248605379240589.
   c7.3  x^2 Log[x] (Log[Log[x+1]] - Log[Log[x]] - 1/(x Log[x])), x -> oo:
         = -1/2 - 1/(2 Log[x]) + O(1/x) -> -1/2.  P1 = 100:
         -0.604034620028079957818435419433; P2 = 10^30:
         -0.507238241365054197127518815315.
   c7.4  Sqrt[x] (Sqrt[x + Sqrt[x + Sqrt[x]]] - Sqrt[x] - 1/2), x -> oo:
         Sqrt[x + Sqrt[x + Sqrt[x]]] = Sqrt[x] + 1/2 + 1/(8 Sqrt[x]) + ...
         -> 1/8.  P1 = 10^4: 0.123758547060251570021092889914; P2 = 10^40:
         0.12499999999999999999875.
   c7.5  (Erf[Sin[x]] - Sin[Erf[x]])/x^3, x -> 0: the series of the numerator
         starts (4 - Pi)/(3 Pi^(3/2)) x^3 -> 0.051386301650969986602749494317.
         P1 = 10^-3: 0.051386244424042204412835357969; P2 = 10^-12:
         0.0513863016509699866027494370897.
   c7.6  n^2 (HarmonicNumber[n] - Log[n] - EulerGamma - 1/(2n)), n -> oo:
         H_n = Log[n] + gamma + 1/(2n) - 1/(12 n^2) + ... -> -1/12.  At
         non-integer n (H_n = PolyGamma[0, n+1] + gamma): P1 = 100+1/3:
         -0.0833325055670544577553412563859; P2 = 10^20+1/3:
         -0.0833333333333333333333333333333 (the next term is 1/(120 n^4)).
   c7.7  (x^(x^(1/x)) - x)/Log[x]^2, x -> oo: x^(1/x) = 1 + Log[x]/x + ...,
         so x^(x^(1/x)) = x + Log[x]^2 + O(Log[x]^4/x) -> 1.  P1 = 10^8:
         1.00000178871303624146689042945; P2 = 10^30:
         1.00000000000000000000000000242.
   c7.8  (Gamma[x] - 1/x + EulerGamma)/x, x -> 0: Gamma[x] = 1/x - gamma
         + (gamma^2 + Pi^2/6) x/2 + ... -> (EulerGamma^2 + Pi^2/6)/2
         = 0.98905599532797255539539565150063.  P1 = 10^-3:
         0.988149496998975588713389726403; P2 = 10^-20:
         0.98905599532797255538632086074.
   c7.9  (x^(x^x) - x)/(x^2 Log[x]^2), x -> 0+: x^x = 1 + x Log[x] + ...,
         x^(x^x) = x (1 + x Log[x]^2 + ...) -> 1.  P1 = 10^-8:
         1.0000016045056037476186910802; P2 = 10^-30:
         1.00000000000000000000000000235.
   c7.10 (PolyLog[2, x] - Pi^2/6)/(1 - x) - Log[1 - x], x -> 1-: with
         u = 1 - x, Li2(x) = Pi^2/6 - Log[x] Log[u] - Li2(u) gives
         -1 + (u/2) Log[u] - u/4 + ... -> -1.  P1 = 1 - 10^-6:
         -1.00000715775999526695053169905; P2 = 1 - 10^-20:
         -1.0000000000000000002327585093.
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"Limit", "Gamma", "Zeta", "Erf", "HarmonicNumber", "PolyLog", "StieltjesGamma"}];

(* The ten expressions in x (plain expressions, no SubValues). *)
ex = {x^2 (Gamma[x + 1/2]/(Sqrt[x] Gamma[x]) - 1 + 1/(8 x)),
      (Zeta[x] - 1/(x - 1) - EulerGamma)/(x - 1),
      x^2 Log[x] (Log[Log[x + 1]] - Log[Log[x]] - 1/(x Log[x])),
      Sqrt[x] (Sqrt[x + Sqrt[x + Sqrt[x]]] - Sqrt[x] - 1/2),
      (Erf[Sin[x]] - Sin[Erf[x]])/x^3,
      x^2 (HarmonicNumber[x] - Log[x] - EulerGamma - 1/(2 x)),
      (x^(x^(1/x)) - x)/Log[x]^2,
      (Gamma[x] - 1/x + EulerGamma)/x,
      (x^(x^x) - x)/(x^2 Log[x]^2),
      (PolyLog[2, x] - Pi^2/6)/(1 - x) - Log[1 - x]};

x0 = {Infinity, 1, Infinity, Infinity, 0, Infinity, Infinity, 0, 0, 1};
dir = {Automatic, Automatic, Automatic, Automatic, Automatic, Automatic, Automatic, Automatic,
       -1, 1};
p1 = {100 + 1/3, 1 + 10^-3, 100, 10^4, 10^-3, 100 + 1/3, 10^8, 10^-3, 10^-8, 1 - 10^-6};
p2 = {10^20 + 1/3, 1 + 10^-20, 10^30, 10^40, 10^-12, 10^20 + 1/3, 10^30, 10^-20, 10^-30, 1 - 10^-20};
dm = {6, 5, 6, 6, 6, 6, 7, 6, 7, 7};

Do[
  With[{e = ex[[k]], a = x0[[k]], dr = dir[[k]], q1 = p1[[k]], q2 = p2[[k]], d1 = dm[[k]],
        lab = "c7." <> ToString[k]},
    run[lab <> " Limit exact",
        If[dr === Automatic, Limit[e, x -> a], Limit[e, x -> a, Direction -> dr]], sig[#, 20] &];
    run[lab <> " N machine", N[e /. x -> q1], sig[#, d1] &];
    run[lab <> " N wp30", N[e /. x -> q2, 30], sig[#, 20] &]],
  {k, 10}];
