(* KummerLogGamma.m -- definite integrals over [0, 1] of Log Gamma times a
   trigonometric polynomial of period 1, by Kummer's Fourier series.

   Loaded lazily by the definite-integration stage "KummerLogGamma"
   (src/calculus/integrate.c), which calls KummerLogGamma`Run[f, x] only for
   the interval {x, 0, 1}.

   THEOREM (Kummer 1847).  On 0 < x < 1,
     Log Gamma(x) = Log(2 Pi)/2
                    + Sum_{n >= 1} ( Cos(2 Pi n x)/(2 n)
                                     + (EulerGamma + Log(2 Pi n)) Sin(2 Pi n x)/(Pi n) ),
   a Fourier series that converges in L^2 (0, 1).  Its coefficients are the
   moments of Log Gamma against the characters E^(2 Pi I k x), k an integer:
     m(0) = INT_0^1 Log Gamma(x) dx = Log(2 Pi)/2                     (Raabe)
     m(k) = INT_0^1 Log Gamma(x) E^(2 Pi I k x) dx
          = 1/(4 |k|) + I Sign(k) (EulerGamma + Log(2 Pi |k|))/(2 Pi |k|).
   So for any trigonometric polynomial T(x) = Sum_k c_k E^(2 Pi I k x) of
   period 1 (Sin and Cos of integer multiples of 2 Pi x, their powers and
   products, constants), INT_0^1 Log Gamma(x) T(x) dx = Sum_k c_k m(k):
   exact, no quadrature.  Log Gamma(1 - x) is reduced to Log Gamma(x) by
   x -> 1 - x, which maps [0, 1] onto itself and T onto a trigonometric
   polynomial of the same kind.  A term of the integrand free of Log Gamma
   is handed back to Integrate.

   The integrand must be of degree one in Log Gamma of the integration
   variable; anything else -- another argument, Log Gamma squared, a
   non-trigonometric x dependence, a frequency that is not an integer
   multiple of 2 Pi -- is $Failed (a decline). *)

KummerLogGamma`savedContextPath = $ContextPath;
$ContextPath = {"System`"};
Begin["KummerLogGamma`Private`"];

(* the moment of Log Gamma against E^(2 Pi I k x) over [0, 1] *)
moment[0] = Log[2 Pi]/2;
moment[k_Integer] := 1/(4 Abs[k]) + I Sign[k] (EulerGamma + Log[2 Pi Abs[k]])/(2 Pi Abs[k]);

hasLG[f_, a_] := !FreeQ[f, LogGamma[a]] || !FreeQ[f, Log[Gamma[a]]];

(* the integral of LogGamma[x] t over [0, 1], or $Failed *)
(* a term c z^e of the Laurent polynomial in z = E^(2 Pi I x) -> {c, e},
   or $Failed when it is not one (a non-integer frequency, a Log[z], ...) *)
monomial[term_, z_] := Module[{e = Exponent[term, z], c},
  If[!IntegerQ[e], $Failed,
    c = Expand[term z^(-e)];
    If[FreeQ[c, z], {c, e}, $Failed]]];

weighted[t_, x_] := Module[{z, tz, terms, mons},
  If[!FreeQ[t, LogGamma] || !FreeQ[t, Gamma], Return[$Failed]];
  If[FreeQ[t, x], Return[t moment[0]]];
  (* E^(p0 + p1 x) -> E^p0 z^(p1/(2 Pi I)): a linear exponent becomes a power
     of z = E^(2 Pi I x); any other x dependence is left in place and fails *)
  tz = Expand[Expand[TrigToExp[t]] /. Power[E, p_] :>
         With[{q = Expand[p]},
           If[PolynomialQ[q, x] && Exponent[q, x] <= 1,
             E^(q /. x -> 0) z^Together[Coefficient[q, x, 1]/(2 Pi I)], E^q]]];
  If[!FreeQ[tz, x], Return[$Failed]];
  terms = If[Head[tz] === Plus, List @@ tz, {tz}];
  mons = monomial[#, z] & /@ terms;
  If[!FreeQ[mons, $Failed], Return[$Failed]];
  Expand[Total[(#[[1]] moment[#[[2]]]) & /@ mons]]];

KummerLogGamma`Run[f0_, x_Symbol] := Module[{f = f0, L, g, w, rest, v, r},
  If[!hasLG[f, x] && hasLG[f, 1 - x], f = f /. x -> 1 - x];
  If[!hasLG[f, x], Return[$Failed]];
  g = Expand[f /. {LogGamma[x] -> L, Log[Gamma[x]] -> L}];
  If[!PolynomialQ[g, L] || Exponent[g, L] =!= 1, Return[$Failed]];
  w = Coefficient[g, L, 1];
  rest = Coefficient[g, L, 0];
  If[!FreeQ[{w, rest}, L], Return[$Failed]];
  v = weighted[w, x];
  If[v === $Failed, Return[$Failed]];
  If[rest =!= 0,
    r = Integrate[rest, {x, 0, 1}];
    If[!FreeQ[r, Integrate], Return[$Failed]];
    v = v + r];
  v];

End[];
$ContextPath = KummerLogGamma`savedContextPath;
