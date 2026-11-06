(* bench-timeout: 2400 *)
(* ==========================================================================
   Experiment 165 -- ISSAC '97 holdout C (SEALED), family c10:
   ODE singularities and high-precision ODE values (counterpart of Problem 10)
   ==========================================================================

       cd benchmarks/165-issac97c-ode-singularities
       ../../Mathilda -file ode_singularities.m

   Ten initial/boundary-value problems: the first real pole of Painleve I
   (a double pole) and of Painleve II (a simple pole), the blow-up of the
   Emden-Fowler equation y'' = x y^2, the first zero crossing of a stiff Van
   der Pol oscillator (mu = 50: a relaxation jump after t ~ 41), a singularly
   perturbed boundary-value problem (boundary layer of width eps at x = 0), a
   chaotic Rossler orbit and a chaotic forced Duffing (Ueda) orbit carried to
   high precision, Kummer's equation started next to its regular singular
   point, a Painleve IV value, and the Hastings-McLeod solution of Painleve II
   (integrated from the Airy asymptote at L back to 0).

   POLES (c10.1-3).  NDSolve runs to the end of its domain (it stops at the
   pole xe); then at x1 = xe - delta the Laurent expansion gives
   p = x1 + k y(x1)/y'(x1), k = 2 for a double pole, 1 for a simple pole: for
   Painleve I y = (x-p)^-2 + O((x-p)^2), so the error is O(delta^5); for
   Painleve II y = +-1/(x-p) + O(x-p), error O(delta^3); for Emden-Fowler
   y = (6/p)(x-p)^-2 - (3.6/p^2)(x-p)^-1 + ..., error O(delta^2).  delta is
   chosen per row so the correction error is below the digits checked.

   TRUE ANSWERS (A = mpmath odefun, an arbitrary-precision Taylor-series
   integrator, at 40-80 digits with the run repeated at two precisions or two
   deltas; B = scipy DOP853 / Radau / LSODA at rtol 1e-12..1e-13 (10-14
   digits), a closed form, or the published value).  Digits: machine 6
   (5 for the two chaotic orbits, whose default-tolerance error is amplified
   ~e^3.5 by the Lyapunov growth, and for the stiff crossing time c10.4, whose
   global error accumulates over the slow phase), wp 15 (12 for the two
   WorkingPrecision -> 32 rows, c10.4 and c10.5).

   c10.1  Painleve I y'' = 6 y^2 + x, y(0) = y'(0) = 0: first pole
          p = 2.61557120988237384514980859223 (A deltas 1e-6/1e-7 agree to
          1e-30; B scipy 2.615571209882378).
   c10.2  Painleve II y'' = 2 y^3 + x y, y(0) = 1, y'(0) = 0: first pole
          p = 1.26277886619324886642763675001 (A deltas 1e-10/1e-11; B scipy
          1.262778866193255).
   c10.3  Emden-Fowler y'' = x y^2, y(0) = 1, y'(0) = 0: blow-up at
          p = 2.84089482764194948767064961441 (A deltas 1e-12/1e-15 from an
          iterated correction, errors falling as delta^2: 1.1e-13, 1.1e-19,
          1.1e-25; B scipy 2.8408948276422, 12 digits).
   c10.4  y'' - 50 (1 - y^2) y' + y = 0, y(0) = 2, y'(0) = 0: first t with
          y(t) = 0: t = 40.9452029284009547381334925019 (A odefun at 30 and 40
          digits; B Radau 40.94520292840087, LSODA 40.9452029275).
   c10.5  eps y'' + (1 + x) y' + y = 0, y(0) = 0, y(1) = 1: integrating once,
          eps y' + (1 + x) y = C, so y = C/eps E^-phi(x) int_0^x E^phi,
          phi = (x + x^2/2)/eps, an Erfi closed form.  y(1/2): machine
          eps = 1/100: 1.33597477915866113864683884112; wp32 eps = 1/1000:
          1.33359306918159219890529130807 (A Erfi form; B mpmath quadrature).
   c10.6  Rossler x' = -y - z, y' = x + y/5, z' = 1/5 + z (x - 57/10), all
          started at 1: x(50) = -7.06153226337260248572175221451 (machine; B
          scipy -7.06153226337858); x(100) = 9.65823343087484059218820188529 (wp40; odefun at
          60 and 80 digits agree, scipy 9.6582334278).
   c10.7  Ueda x'' + x'/20 + x^3 = (15/2) Cos[t], x(0) = 3, x'(0) = 4:
          x(20) = -0.725120577945468811698443108149 (machine; B scipy
          -0.725120577944944); x(40) = 0.595491880642409614465877717234 (wp40; odefun at 60
          and 80 digits agree, scipy 0.59549188064145).
   c10.8  Kummer x y'' + (1/2 - x) y' - y/3 = 0, regular at 0, started at
          x0 = 1/1000 from 11 terms of its series: y(5) = M(1/3, 1/2, 5)
          = 77.4157788604427077180591019652 (A mpmath hyp1f1; B PARI hypergeom).
   c10.9  Painleve IV (alpha = beta = 0) y'' = y'^2/(2y) + 3y^3/2 + 4x y^2
          + 2x^2 y, y(0) = 1, y'(0) = 0: y(1/2) = 1.33910862429232758493332543647
          (A odefun at 40 and 60 digits; B scipy 1.339108624292324).
   c10.10 Hastings-McLeod: y'' = x y + 2 y^3 with y ~ AiryAi[x] as x -> oo,
          integrated from L back to 0 (machine L = 8, wp40 L = 20):
          y(0) = 0.367061551548078427747792113176 (A odefun from L = 20 and 24;
          B the published 0.36706155154807843).
   ========================================================================== *)

Get["../harness.m"]; Get["../issac97_common.m"];
require[{"NDSolve", "FindRoot", "InterpolatingFunction", "AiryAi", "AiryAiPrime"}];

wpo[wp_] := If[wp === MachinePrecision, {MaxSteps -> 200000},
               {WorkingPrecision -> wp, PrecisionGoal -> wp/2, AccuracyGoal -> wp/2, MaxSteps -> 200000}];

(* ---- c10.1-3 poles: domain end, then the Laurent correction ------------- *)
pole[rhs_, y0_, k_, delta_, wp_] := Module[{s, f, xe, x1},
  s = Quiet[NDSolve[{y''[x] == rhs, y[0] == y0[[1]], y'[0] == y0[[2]]}, y, {x, 0, 4},
                    Evaluate[Sequence @@ wpo[wp]]]];
  f = y /. First[s];
  xe = First[f["Domain"]][[2]];
  x1 = xe - delta;
  x1 + k f[x1]/f'[x1]];
run["c10.1 NDSolve pole machine", pole[6 y[x]^2 + x, {0, 0}, 2, 10^-3, MachinePrecision], sig[#, 6] &];
run["c10.1 NDSolve pole wp40", pole[6 y[x]^2 + x, {0, 0}, 2, 10^-6, 40], sig[#, 15] &];
run["c10.2 NDSolve pole machine", pole[2 y[x]^3 + x y[x], {1, 0}, 1, 10^-6, MachinePrecision], sig[#, 6] &];
run["c10.2 NDSolve pole wp40", pole[2 y[x]^3 + x y[x], {1, 0}, 1, 10^-11, 40], sig[#, 15] &];
run["c10.3 NDSolve pole machine", pole[x y[x]^2, {1, 0}, 2, 10^-8, MachinePrecision], sig[#, 6] &];
run["c10.3 NDSolve pole wp40", pole[x y[x]^2, {1, 0}, 2, 10^-16, 40], sig[#, 15] &];

(* ---- c10.4 stiff Van der Pol: first zero ------------------------------- *)
vdp[wp_] := Module[{s, f},
  s = NDSolve[{y''[t] - 50 (1 - y[t]^2) y'[t] + y[t] == 0, y[0] == 2, y'[0] == 0}, y, {t, 0, 45},
              Evaluate[Sequence @@ wpo[wp]]];
  f = y /. First[s];
  ruleVal[FindRoot[f[t], {t, 409/10, 41}, Evaluate[Sequence @@ If[wp === MachinePrecision, {},
                                                {WorkingPrecision -> wp}]]], t]];
run["c10.4 NDSolve stiff machine", vdp[MachinePrecision], sig[#, 5] &];
run["c10.4 NDSolve stiff wp32", vdp[32], sig[#, 12] &];

(* ---- c10.5 boundary layer ---------------------------------------------- *)
bl[eps_, wp_] := (y[1/2] /. First[NDSolve[{eps y''[x] + (1 + x) y'[x] + y[x] == 0, y[0] == 0, y[1] == 1},
                                         y, {x, 0, 1}, Evaluate[Sequence @@ wpo[wp]]]]);
run["c10.5 NDSolve BVP machine eps=1/100", bl[1/100, MachinePrecision], sig[#, 6] &];
run["c10.5 NDSolve BVP wp32 eps=1/1000", bl[1/1000, 32], sig[#, 12] &];

(* ---- c10.6 Rossler ------------------------------------------------------ *)
ross[T_, wp_] := (a[T] /. First[NDSolve[{a'[t] == -b[t] - c[t], b'[t] == a[t] + b[t]/5,
                                         c'[t] == 1/5 + c[t] (a[t] - 57/10), a[0] == 1, b[0] == 1, c[0] == 1},
                                        {a, b, c}, {t, 0, T}, Evaluate[Sequence @@ wpo[wp]]]]);
run["c10.6 NDSolve chaotic machine T=50", ross[50, MachinePrecision], sig[#, 5] &];
run["c10.6 NDSolve chaotic wp40 T=100", ross[100, 40], sig[#, 15] &];

(* ---- c10.7 Ueda --------------------------------------------------------- *)
ueda[T_, wp_] := (y[T] /. First[NDSolve[{y''[t] + y'[t]/20 + y[t]^3 == 15/2 Cos[t], y[0] == 3, y'[0] == 4},
                                        y, {t, 0, T}, Evaluate[Sequence @@ wpo[wp]]]]);
run["c10.7 NDSolve chaotic machine T=20", ueda[20, MachinePrecision], sig[#, 5] &];
run["c10.7 NDSolve chaotic wp40 T=40", ueda[40, 40], sig[#, 15] &];

(* ---- c10.8 Kummer next to its regular singular point -------------------- *)
kser = Sum[Pochhammer[1/3, k]/(Pochhammer[1/2, k] k!) x^k, {k, 0, 10}];
kum[wp_] := (
  y[5] /. First[NDSolve[{x y''[x] + (1/2 - x) y'[x] - y[x]/3 == 0,
                         y[1/1000] == (kser /. x -> 1/1000), y'[1/1000] == (D[kser, x] /. x -> 1/1000)},
                        y, {x, 1/1000, 5}, Evaluate[Sequence @@ wpo[wp]]]]);
run["c10.8 NDSolve singular start machine", kum[MachinePrecision], sig[#, 6] &];
run["c10.8 NDSolve singular start wp40", kum[40], sig[#, 15] &];

(* ---- c10.9 Painleve IV -------------------------------------------------- *)
p4[wp_] := (y[1/2] /. First[NDSolve[{y''[x] == y'[x]^2/(2 y[x]) + 3 y[x]^3/2 + 4 x y[x]^2 + 2 x^2 y[x],
                                     y[0] == 1, y'[0] == 0}, y, {x, 0, 1/2}, Evaluate[Sequence @@ wpo[wp]]]]);
run["c10.9 NDSolve machine", p4[MachinePrecision], sig[#, 6] &];
run["c10.9 NDSolve wp40", p4[40], sig[#, 15] &];

(* ---- c10.10 Hastings-McLeod from the Airy asymptote --------------------- *)
(* y(L) = AiryAi[L] is 1e-7 (L = 8) or 2e-27 (L = 20): an absolute error
   tolerance would swamp it, so the error control is purely relative. *)
hm[L_, wp_] := (y[0] /. First[NDSolve[{y''[x] == x y[x] + 2 y[x]^3, y[L] == AiryAi[L], y'[L] == AiryAiPrime[L]},
                                      y, {x, L, 0}, AccuracyGoal -> Infinity, MaxSteps -> 200000,
                                      Evaluate[Sequence @@ If[wp === MachinePrecision, {},
                                                {WorkingPrecision -> wp, PrecisionGoal -> wp/2}]]]]);
run["c10.10 NDSolve backward machine L=8", hm[8, MachinePrecision], sig[#, 6] &];
run["c10.10 NDSolve backward wp40 L=20", hm[20, 40], sig[#, 15] &];
