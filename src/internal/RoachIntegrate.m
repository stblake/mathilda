(* RoachIntegrate.m -- definite integration by Meijer G functions:
   Integrate[f, {x, a, b}, Method -> "Roach"].

   Kelly Roach, "Meijer G Function Representations" (ISSAC '97) sec. 5 and
   "Solving Integrals with the Quantum Computer Algebra System" (2008). The
   integration theorems and their validity conditions are those of
   Prudnikov, Brychkov, Marichev, Integrals and Series vol. 3, sec. 2.24.1
   (the 35-case condition set, as tabulated at functions.wolfram.com
   07.34.21.0011 and implemented in SymPy's meijerint, BSD licence). The
   G -> named functions step is Roach's Formula algorithm
   (src/special_functions/hrep_core.c), reached through FunctionExpand.

   Loaded lazily by the definite stage METHOD_ROACH (src/calculus/integrate.c).
   Every symbol is qualified in Roach`.

   PIPELINE (Roach`Run):
     1. interval -> (0, Infinity): x -> x + a; (-Infinity, b) by x -> -x;
        the whole line split at 0; a finite [a, b] by x -> a + (b - a) x with
        the indicator Roach`Ind[x] = UnitStep[1 - x];
     2. linearity over the expanded integrand;
     3. each term C x^s F1(c1 x^q1) F2(c2 x^q2) with every factor a G function
        (Roach`GForm, the lookup table of PBM sec. 8.4): the Mellin transform
        of one G, or the convolution of two G (with argument inflation for a
        rational ratio of powers), under the validity conditions;
     4. fallbacks: t = E^-x for exponential integrands, t = Tan[x] for
        trigonometric ones on [0, Pi/2), integration by parts for an inverse
        function times a power;
     5. the resulting G -> FunctionExpand -> Refine under the assumptions.

   CONDITIONS. Each is decided under the call's assumptions: True drops it,
   False makes the term decline, anything else is kept and the answer is a
   ConditionalExpression. A result of a substitution route is checked by
   NIntegrate at a FindInstance point of the assumptions. *)

(* Every unqualified symbol below (pattern variables, locals) lives in
   Roach`Private`. In Global` a pattern variable a_ would be the user's own
   symbol a, and the integrand's parameters would be captured. *)
Roach`savedContextPath = $ContextPath;
$ContextPath = {"System`"};
Begin["Roach`Private`"];

(* ---------------------------------------------------------------- *)
(* Small helpers                                                    *)
(* ---------------------------------------------------------------- *)

(* {c, q} with e == c x^q, c free of x, q a nonzero rational; $Failed otherwise *)
Roach`CQ[e_, x_] := Module[{q, c},
  If[FreeQ[e, x], Return[$Failed]];
  q = Together[x D[e, x]/e];
  If[!FreeQ[q, x] || q === 0, Return[$Failed]];
  c = Together[e/x^q];
  If[!FreeQ[c, x], c = Simplify[PowerExpand[c]]];
  If[!FreeQ[c, x], $Failed, {c, q}]]

(* G term: Roach`g[C, {an, ap}, {bm, bq}, c, q] = C G^{m,n}_{p,q}(c x^q | a; b) *)
Roach`mk[C_, an_, ap_, bm_, bq_, arg_, x_] := With[{cq = Roach`CQ[arg, x]},
  If[cq === $Failed, $Failed, {Roach`g[C, {an, ap}, {bm, bq}, cq[[1]], cq[[2]]]}]]

(* Roach`Ind0[t] is inert and stands for UnitStep[1 - t] *)

(* True when t == x^q with q > 0: then UnitStep[1 - t] == UnitStep[1 - x] on
   x > 0, so a G-form carrying theta(1 - t) may stand for Roach`Ind0[x]. With
   a coefficient c != 1 it would be theta(1 - c x^q), a different interval. *)
Roach`unitQ[t_, x_] := With[{cq = Roach`CQ[t, x]},
  cq =!= $Failed && cq[[1]] === 1 && TrueQ[cq[[2]] > 0]]

(* ---------------------------------------------------------------- *)
(* Function -> G table (PBM vol. 3 sec. 8.4)                         *)
(* Returns a list of G terms whose sum is the factor, or $Failed.    *)
(* ---------------------------------------------------------------- *)

Roach`GForm[f_, x_] := Roach`GForm1[f, x]

Roach`GForm1[Power[E, t_], x_] := Roach`mk[1, {}, {}, {0}, {}, -t, x]
Roach`GForm1[Power[b_ + t_, ma_], x_] /; FreeQ[b, x] && FreeQ[ma, x] && !FreeQ[t, x] &&
    Roach`CQ[t, x] =!= $Failed && !MatchQ[ma, _Integer?Positive] :=
  With[{a = -ma}, Roach`mk[b^(-a)/Gamma[a], {1 - a}, {}, {0}, {}, t/b, x]]
(* (1 - t)^(a-1) on 0 < t < 1, with the interval indicator *)
Roach`GForm1[Power[1 - t_, am1_] Roach`Ind0[x_], x_] /; FreeQ[am1, x] && Roach`unitQ[t, x] :=
  With[{a = am1 + 1}, Roach`mk[Gamma[a], {}, {a}, {0}, {}, t, x]]
Roach`GForm1[Roach`Ind0[x_], x_] := Roach`mk[1, {}, {1}, {0}, {}, x, x]
Roach`GForm1[Power[Log[t_], n_Integer?Positive] Roach`Ind0[x_], x_] /; Roach`unitQ[t, x] :=
  Roach`mk[(-1)^n n!, {}, ConstantArray[1, n + 1], ConstantArray[0, n + 1], {}, t, x]
Roach`GForm1[Log[t_] Roach`Ind0[x_], x_] /; Roach`unitQ[t, x] :=
  Roach`mk[-1, {}, {1, 1}, {0, 0}, {}, t, x]
Roach`GForm1[Log[1 + t_], x_] /; Roach`CQ[t, x] =!= $Failed :=
  Roach`mk[1, {1, 1}, {}, {1}, {0}, t, x]
Roach`GForm1[Sin[t_], x_] := Roach`mk[Sqrt[Pi], {}, {}, {1/2}, {0}, t^2/4, x]
Roach`GForm1[Cos[t_], x_] := Roach`mk[Sqrt[Pi], {}, {}, {0}, {1/2}, t^2/4, x]
Roach`GForm1[Sinh[t_], x_] := Roach`mk[Pi^(3/2), {}, {1}, {1/2}, {1, 0}, t^2/4, x]
Roach`GForm1[Cosh[t_], x_] := Roach`mk[Pi^(3/2), {}, {1/2}, {0}, {1/2, 1/2}, t^2/4, x]
Roach`GForm1[ArcTan[t_], x_] := Roach`mk[1/2, {1/2, 1}, {}, {1/2}, {0}, t^2, x]
(* ArcCot(t) = ArcTan(1/t) for t > 0: the ArcTan G at 1/t^2, flipped *)
Roach`GForm1[ArcCot[t_], x_] := Roach`mk[1/2, {1/2}, {1}, {1/2, 0}, {}, t^2, x]
Roach`GForm1[ExpIntegralEi[mt_], x_] /; Roach`CQ[-mt, x] =!= $Failed &&
    TrueQ[Simplify[Roach`CQ[-mt, x][[1]] > 0]] :=
  Roach`mk[-1, {}, {1}, {0, 0}, {}, -mt, x]       (* Ei(-t) = -E1(t), t > 0 *)
Roach`GForm1[ExpIntegralE[a_, t_], x_] /; FreeQ[a, x] := Roach`mk[1, {}, {a}, {a - 1, 0}, {}, t, x]
Roach`GForm1[SinIntegral[t_], x_] := Roach`mk[Sqrt[Pi]/2, {1}, {}, {1/2}, {0, 0}, t^2/4, x]
Roach`GForm1[CosIntegral[t_], x_] := Roach`mk[-Sqrt[Pi]/2, {}, {1}, {0, 0}, {1/2}, t^2/4, x]
Roach`GForm1[Erf[t_], x_] := Roach`mk[1/Sqrt[Pi], {1}, {}, {1/2}, {0}, t^2, x]
Roach`GForm1[Erfc[t_], x_] := Roach`mk[1/Sqrt[Pi], {}, {1}, {0, 1/2}, {}, t^2, x]
Roach`GForm1[Gamma[a_, t_], x_] /; FreeQ[a, x] := Roach`mk[1, {}, {1}, {a, 0}, {}, t, x]
Roach`GForm1[FresnelS[t_], x_] := Roach`mk[1/2, {1}, {}, {3/4}, {0, 1/4}, Pi^2 t^4/16, x]
Roach`GForm1[FresnelC[t_], x_] := Roach`mk[1/2, {1}, {}, {1/4}, {0, 3/4}, Pi^2 t^4/16, x]
Roach`GForm1[BesselJ[a_, t_], x_] /; FreeQ[a, x] := Roach`mk[1, {}, {}, {a/2}, {-a/2}, t^2/4, x]
Roach`GForm1[BesselY[a_, t_], x_] /; FreeQ[a, x] :=
  Roach`mk[1, {}, {-(a + 1)/2}, {a/2, -a/2}, {-(a + 1)/2}, t^2/4, x]
Roach`GForm1[BesselI[a_, t_], x_] /; FreeQ[a, x] :=
  Roach`mk[Pi, {}, {(1 + a)/2}, {a/2}, {-a/2, (1 + a)/2}, -t^2/4, x]
Roach`GForm1[BesselK[a_, t_], x_] /; FreeQ[a, x] := Roach`mk[1/2, {}, {}, {a/2, -a/2}, {}, t^2/4, x]
Roach`GForm1[StruveH[a_, t_], x_] /; FreeQ[a, x] :=
  Roach`mk[1, {(a + 1)/2}, {}, {(a + 1)/2}, {-a/2, a/2}, t^2/4, x]
Roach`GForm1[EllipticK[t_], x_] := Roach`mk[1/2, {1/2, 1/2}, {}, {0}, {0}, -t, x]
Roach`GForm1[EllipticE[t_], x_] := Roach`mk[-1/4, {1/2, 3/2}, {}, {0}, {0}, -t, x]
(* products of Bessel functions of one argument *)
Roach`GForm1[Power[BesselJ[a_, t_], 2], x_] /; FreeQ[a, x] :=
  Roach`mk[1/Sqrt[Pi], {1/2}, {}, {a}, {-a, 0}, t^2, x]
Roach`GForm1[f_, x_] := $Failed

(* ---------------------------------------------------------------- *)
(* G-function bookkeeping                                           *)
(* ---------------------------------------------------------------- *)

Roach`mnpq[Roach`g[_, {an_, ap_}, {bm_, bq_}, _, _]] :=
  {Length[bm], Length[an], Length[an] + Length[ap], Length[bm] + Length[bq]}
Roach`nu[Roach`g[_, {an_, ap_}, {bm_, bq_}, _, _]] := Total[Join[bm, bq]] - Total[Join[an, ap]]
Roach`delta[g_Roach`g] := With[{k = Roach`mnpq[g]}, k[[1]] + k[[2]] - (k[[3]] + k[[4]])/2]
Roach`shift[Roach`g[C_, {an_, ap_}, {bm_, bq_}, c_, q_], e_] :=
  Roach`g[C, {an + e, ap + e}, {bm + e, bq + e}, c, q]
(* G(1/z) -> G'(z) *)
Roach`flip[Roach`g[C_, {an_, ap_}, {bm_, bq_}, c_, q_]] :=
  Roach`g[C, {1 - bm, 1 - bq}, {1 - an, 1 - ap}, 1/c, -q]
(* G(z) = C' G'(z^n): Gauss multiplication (Luke p. 150 eq. 5) *)
Roach`inflate[g : Roach`g[C_, {an_, ap_}, {bm_, bq_}, c_, q_], n_Integer] := Module[{v, Cn, inf},
  v = (Length[an] + Length[ap]) - (Length[bm] + Length[bq]);
  Cn = n^(1 + Roach`nu[g] + v/2)/(2 Pi)^((n - 1) Roach`delta[g]);
  inf[l_] := Flatten[Table[(a + i)/n, {a, l}, {i, 0, n - 1}]];
  Roach`g[C Cn, {inf[an], inf[ap]}, {inf[bm], inf[bq]}, c^n n^(n v), q n]]

Roach`MG[Roach`g[C_, {an_, ap_}, {bm_, bq_}, c_, q_], z_] := C MeijerG[{an, ap}, {bm, bq}, z]

(* Mellin transform of G at w: Prod Gamma(b_j + w) Gamma(1 - a_j - w) / ... *)
Roach`mellin[Roach`g[_, {an_, ap_}, {bm_, bq_}, _, _], w_] :=
  Times @@ (Gamma[# + w] & /@ bm) Times @@ (Gamma[1 - # - w] & /@ an) /
   (Times @@ (Gamma[1 - # - w] & /@ bq) Times @@ (Gamma[# + w] & /@ ap))

(* ---------------------------------------------------------------- *)
(* Validity conditions                                              *)
(* ---------------------------------------------------------------- *)

Roach`argOf[c_] := Arg[c]

(* int_0^oo G(eta x) dx exists (PBM 2.24.2, Luke 5.6.1); helper -> only at oo *)
Roach`cond1[g0 : Roach`g[_, {an_, ap_}, {bm_, bq_}, eta_, 1], helper_ : False] := Module[
  {m, n, p, q, delta, nu, c3, c3s, c4, conds = {}, ae = Roach`argOf[eta], s, g = g0},
  {m, n, p, q} = Roach`mnpq[g];
  If[p > q, Return[Roach`cond1[Roach`g[1, {1 - bm, 1 - bq}, {1 - an, 1 - ap}, 1/eta, 1], helper]]];
  delta = Roach`delta[g]; nu = Roach`nu[g];
  c3 = And @@ Join[(-Re[#] < 1) & /@ bm, (1 < 1 - Re[#]) & /@ an];
  c3s = And[c3, And @@ ((-Re[#] < 1) & /@ bq), And @@ ((1 < 1 - Re[#]) & /@ ap)];
  c4 = (-Re[nu] + (q + 1 - p)/2 > q - p);
  With[{tmp = {delta > 0, Abs[ae] < delta Pi}, extra = If[helper, {}, {eta != 0, c3}]},
    AppendTo[conds, And @@ Join[{1 <= n, p < q, 1 <= m}, tmp, extra]];
    AppendTo[conds, And @@ Join[{1 <= p, 1 <= m, q == p + 1, !(n == 0 && m == p + 1)}, tmp, extra]];
    AppendTo[conds, And @@ Join[{1 <= p, q == p},
      Table[Abs[ae] != (delta - 2 k) Pi, {k, 0, Ceiling[delta/2]}], tmp, extra]]];
  AppendTo[conds, And @@ Join[{n == 0, p + 1 <= m, m <= q, Abs[ae] < delta Pi}, If[helper, {}, {c3}]]];
  AppendTo[conds, And @@ Join[{p < q, 1 <= m, delta > 0, Abs[ae] == delta Pi}, If[helper, {}, {c3, c4}]]];
  AppendTo[conds, And @@ Join[{p <= q - 2, delta == 0, Abs[ae] == 0}, If[helper, {}, {c3, c4}]]];
  s = Total[MapThread[#2 - #1 &, {Take[Join[an, ap], Min[p, q]], Take[Join[bm, bq], Min[p, q]]}]];
  AppendTo[conds, And @@ Join[{p == q, delta == 0, ae == 0, eta != 0}, If[helper, {}, {c3}], {Re[s] < 0}]];
  AppendTo[conds, And @@ Join[{delta > 0, Abs[ae] < delta Pi}, If[helper, {}, {c3}]]];
  Or @@ conds]

(* int_0^oo G1(sigma x) G2(omega x) dx: PBM 2.24.1.1, conditions 1-35 (k = l = 1) *)
Roach`cond2[g1 : Roach`g[_, {c1n_, c1p_}, {d1m_, d1q_}, sigma_, 1],
            g2 : Roach`g[_, {a2n_, a2p_}, {b2m_, b2q_}, omega_, 1]] := Module[
  {s, t, u, v, m, n, p, q, bstar, cstar, rho, mu, phi, eta, psi, theta,
   C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, C12, C13, C14, C15, conds = {},
   asig = Roach`argOf[sigma], aom = Roach`argOf[omega], zos, zso, z0, lc, r, mt1, mt2},
  {s, t, u, v} = Roach`mnpq[g1]; {m, n, p, q} = Roach`mnpq[g2];
  bstar = s + t - (u + v)/2; cstar = m + n - (p + q)/2;
  rho = Roach`nu[g1] + (u - v)/2 + 1; mu = Roach`nu[g2] + (p - q)/2 + 1;
  phi = q - p - (v - u); eta = 1 - (v - u) - mu - rho;
  psi = If[q == p, Indeterminate, (Pi (q - m - n) + Abs[aom])/(q - p)];
  theta = If[v == u, Indeterminate, (Pi (v - s - t) + Abs[asig])/(v - u)];
  C1 = And @@ Flatten[{Table[!(IntegerQ[i - j] && i - j > 0), {i, c1n}, {j, d1m}],
                       Table[!(IntegerQ[i - j] && i - j > 0), {i, a2n}, {j, b2m}]}];
  C2 = And @@ Flatten[Table[Re[1 + i + j] > 0, {i, d1m}, {j, b2m}]];
  C3 = And @@ Flatten[Table[Re[1 + i + j] < 2, {i, c1n}, {j, a2n}]];
  C4 = And @@ ((p - q) Re[#] - Re[mu] > -3/2 & /@ c1n);
  C5 = And @@ ((p - q) Re[1 + #] - Re[mu] > -3/2 & /@ d1m);
  C6 = And @@ ((u - v) Re[#] - Re[rho] > -3/2 & /@ a2n);
  C7 = And @@ ((u - v) Re[1 + #] - Re[rho] > -3/2 & /@ b2m);
  C8 = Abs[phi] + 2 Re[(rho - 1) (q - p) + (v - u) (q - p) + (mu - 1) (v - u)] > 0;
  C9 = Abs[phi] - 2 Re[(rho - 1) (q - p) + (v - u) (q - p) + (mu - 1) (v - u)] > 0;
  C10 = Abs[asig] < bstar Pi; C11 = Abs[asig] == bstar Pi;
  C12 = Abs[aom] < cstar Pi;  C13 = Abs[aom] == cstar Pi;
  z0 = Exp[-(bstar + cstar) Pi I]; zos = z0 omega/sigma; zso = z0 sigma/omega;
  C14 = Or[And[phi == 0, bstar - 1 + cstar <= 0,
               Or[And[zos != 1, Abs[Arg[1 - zos]] < Pi], And[Re[mu + rho + v - u] < 1, zos == 1]]],
           And[phi == 0, cstar - 1 + bstar <= 0,
               Or[And[zso != 1, Abs[Arg[1 - zso]] < Pi], And[Re[mu + rho + q - p] < 1, zso == 1]]]];
  lc = If[q == p || v == u, Indeterminate,
          (q - p) Abs[omega]^(1/(q - p)) Cos[psi] + (v - u) Abs[sigma]^(1/(v - u)) Cos[theta]];
  C15 = If[lc === Indeterminate, False, lc > 0];
  Roach`pos[e_] := TrueQ[e > 0];
  conds = {
   And[m n s t != 0, Roach`pos[bstar], Roach`pos[cstar], C1, C2, C3, C10, C12],
   And[u == v, bstar == 0, Roach`pos[cstar], sigma > 0, Re[rho] < 1, C1, C2, C3, C12],
   And[p == q, cstar == 0, Roach`pos[bstar], omega > 0, Re[mu] < 1, C1, C2, C3, C10],
   And[p == q, u == v, bstar == 0, cstar == 0, sigma > 0, omega > 0, Re[mu] < 1, Re[rho] < 1,
       sigma != omega, C1, C2, C3],
   And[p == q, u == v, bstar == 0, cstar == 0, sigma > 0, omega > 0, Re[mu + rho] < 1,
       omega != sigma, C1, C2, C3],
   And[p > q, Roach`pos[s], Roach`pos[bstar], cstar >= 0, C1, C2, C3, C5, C10, C13],
   And[p < q, Roach`pos[t], Roach`pos[bstar], cstar >= 0, C1, C2, C3, C4, C10, C13],
   And[u > v, Roach`pos[m], Roach`pos[cstar], bstar >= 0, C1, C2, C3, C7, C11, C12],
   And[u < v, Roach`pos[n], Roach`pos[cstar], bstar >= 0, C1, C2, C3, C6, C11, C12],
   And[p > q, u == v, bstar == 0, cstar >= 0, sigma > 0, Re[rho] < 1, C1, C2, C3, C5, C13],
   And[p < q, u == v, bstar == 0, cstar >= 0, sigma > 0, Re[rho] < 1, C1, C2, C3, C4, C13],
   And[p == q, u > v, bstar >= 0, cstar == 0, omega > 0, Re[mu] < 1, C1, C2, C3, C7, C11],
   And[p == q, u < v, bstar >= 0, cstar == 0, omega > 0, Re[mu] < 1, C1, C2, C3, C6, C11],
   And[p < q, u > v, bstar >= 0, cstar >= 0, C1, C2, C3, C4, C7, C11, C13],
   And[p > q, u < v, bstar >= 0, cstar >= 0, C1, C2, C3, C5, C6, C11, C13],
   And[p > q, u > v, bstar >= 0, cstar >= 0, C1, C2, C3, C5, C7, C8, C11, C13, C14],
   And[p < q, u < v, bstar >= 0, cstar >= 0, C1, C2, C3, C4, C6, C9, C11, C13, C14],
   And[t == 0, Roach`pos[s], Roach`pos[bstar], Roach`pos[phi], C1, C2, C10],
   And[s == 0, Roach`pos[t], Roach`pos[bstar], TrueQ[phi < 0], C1, C3, C10],
   And[n == 0, Roach`pos[m], Roach`pos[cstar], TrueQ[phi < 0], C1, C2, C12],
   And[m == 0, Roach`pos[n], Roach`pos[cstar], Roach`pos[phi], C1, C3, C12],
   And[s t == 0, Roach`pos[bstar], Roach`pos[cstar], C1, C2, C3, C10, C12],
   And[m n == 0, Roach`pos[bstar], Roach`pos[cstar], C1, C2, C3, C10, C12]};
  (* Luke 1969: G1 decays exponentially, the Mellin transform of G2 exists *)
  mt1 = Roach`cond1[Roach`g[1, {c1n, c1p}, {d1m, d1q}, sigma, 1], True];
  mt2 = Roach`cond1[Roach`g[1, {a2n, a2p}, {b2m, b2q}, omega, 1], True];
  conds = Join[conds, {
   And[mt2, t == 0, u < s, Roach`pos[bstar], C10, C1, C2, C3],
   And[mt2, s == 0, v < t, Roach`pos[bstar], C10, C1, C2, C3],
   And[mt1, n == 0, p < m, Roach`pos[cstar], C12, C1, C2, C3],
   And[mt1, m == 0, q < n, Roach`pos[cstar], C12, C1, C2, C3]}];
  r = Or @@ conds;
  If[r =!= False, Return[r]];
  conds = Join[conds, {
   And[m + n > p, t == 0, phi == 0, Roach`pos[s], Roach`pos[bstar], TrueQ[cstar < 0],
       Abs[aom] < (m + n - p + 1) Pi, C1, C2, C10, C14, C15],
   And[m + n > q, s == 0, phi == 0, Roach`pos[t], Roach`pos[bstar], TrueQ[cstar < 0],
       Abs[aom] < (m + n - q + 1) Pi, C1, C3, C10, C14, C15],
   And[p == q - 1, t == 0, phi == 0, Roach`pos[s], Roach`pos[bstar], cstar >= 0,
       cstar Pi < Abs[aom], C1, C2, C10, C14, C15],
   And[p == q + 1, s == 0, phi == 0, Roach`pos[t], Roach`pos[bstar], cstar >= 0,
       cstar Pi < Abs[aom], C1, C3, C10, C14, C15],
   And[p < q - 1, t == 0, phi == 0, Roach`pos[s], Roach`pos[bstar], cstar >= 0,
       cstar Pi < Abs[aom], Abs[aom] < (m + n - p + 1) Pi, C1, C2, C10, C14, C15],
   And[p > q + 1, s == 0, phi == 0, Roach`pos[t], Roach`pos[bstar], cstar >= 0,
       cstar Pi < Abs[aom], Abs[aom] < (m + n - q + 1) Pi, C1, C3, C10, C14, C15],
   And[n == 0, phi == 0, s + t > 0, Roach`pos[m], Roach`pos[cstar], TrueQ[bstar < 0],
       Abs[asig] < (s + t - u + 1) Pi, C1, C2, C12, C14, C15],
   And[m == 0, phi == 0, s + t > v, Roach`pos[n], Roach`pos[cstar], TrueQ[bstar < 0],
       Abs[asig] < (s + t - v + 1) Pi, C1, C3, C12, C14, C15],
   And[n == 0, phi == 0, u == v - 1, Roach`pos[m], Roach`pos[cstar], bstar >= 0,
       bstar Pi < Abs[asig], Abs[asig] < (bstar + 1) Pi, C1, C2, C12, C14, C15],
   And[m == 0, phi == 0, u == v + 1, Roach`pos[n], Roach`pos[cstar], bstar >= 0,
       bstar Pi < Abs[asig], Abs[asig] < (bstar + 1) Pi, C1, C3, C12, C14, C15],
   And[n == 0, phi == 0, u < v - 1, Roach`pos[m], Roach`pos[cstar], bstar >= 0,
       bstar Pi < Abs[asig], Abs[asig] < (s + t - u + 1) Pi, C1, C2, C12, C14, C15],
   And[m == 0, phi == 0, u > v + 1, Roach`pos[n], Roach`pos[cstar], bstar >= 0,
       bstar Pi < Abs[asig], Abs[asig] < (s + t - v + 1) Pi, C1, C3, C12, C14, C15]}];
  Or @@ conds]

(* decide a condition under the assumptions: True | False | residual.
   Each relational atom is refined on its own, which is fast and catches what
   the assumptions imply (b > 0 && b < c from c > b > 0). The And/Or/Not
   structure then folds by evaluation. *)
Roach`atom[r_, assum_] := Module[{v = TimeConstrained[Refine[r, assum], 1, r]},
  If[v === True || v === False, v,
    With[{w = TimeConstrained[Simplify[r, assum], 1, v]}, If[w === True || w === False, w, v]]]]
Roach`decide[c_, assum_] := Module[{r},
  r = c /. rel : (_Less | _Greater | _LessEqual | _GreaterEqual | _Equal | _Unequal) :>
            Roach`atom[rel, assum];
  If[r === True || r === False, r, TimeConstrained[Simplify[r, assum], 2, r]]]

(* ---------------------------------------------------------------- *)
(* The two theorems                                                 *)
(* ---------------------------------------------------------------- *)

(* the sign of an exponent under the call's assumptions: 1, -1 or 0 (unknown) *)
Roach`$assum = True;
Roach`sgn[q_] := Which[TrueQ[Refine[q > 0, Roach`$assum]], 1,
                       TrueQ[Refine[q < 0, Roach`$assum]], -1, True, 0]

(* int_0^oo K x^s G(c x^q) dx  ->  {value, condition} *)
Roach`int1[K_, s_, g0 : Roach`g[C_, ap_, bp_, c_, q_], x_] := Module[{e, g, gs, sq = Roach`sgn[q]},
  (* t = x^q: int x^s G(c x^q) dx = (1/|q|) int t^((s+1)/q - 1) G(c t) dt *)
  If[sq == 0, Return[$Failed]];
  g = If[sq < 0, Roach`flip[g0], g0];
  e = (s + 1)/(sq q) - 1;
  gs = Roach`shift[Roach`g[1, g[[2]], g[[3]], g[[4]], 1], e];
  (* int_0^oo t^e G(c t) dt = c^(-e-1) M[G](e + 1) *)
  {K C/(sq q) g[[4]]^(-e - 1) Roach`mellin[Roach`g[1, g[[2]], g[[3]], 1, 1], e + 1],
   Roach`cond1[Roach`g[1, gs[[2]], gs[[3]], g[[4]], 1]]}]

(* int_0^oo K x^s G1(c1 x^q1) G2(c2 x^q2) dx -> {value, condition} (Saxena) *)
Roach`int2[K_, s_, g10_, g20_, x_] := Module[
  {g1 = g10, g2 = g20, q1, q2, tau, r1, r2, b, e, a1, a2, h1, h2, val, cond,
   sig, om, s1, t1, u1, v1},
  If[g1[[5]] < 0, g1 = Roach`flip[g1]];
  If[g2[[5]] < 0, g2 = Roach`flip[g2]];
  q1 = g1[[5]]; q2 = g2[[5]];
  If[!MatchQ[q1, _Integer | _Rational] || !MatchQ[q2, _Integer | _Rational], Return[$Failed]];
  (* common power tau = lcm of numerators over gcd of denominators *)
  tau = LCM[Numerator[q1] Denominator[q2], Numerator[q2] Denominator[q1]];
  r1 = tau/(Numerator[q1] Denominator[q2]); r2 = tau/(Numerator[q2] Denominator[q1]);
  If[r1 != 1, g1 = Roach`inflate[g1, r1]];
  If[r2 != 1, g2 = Roach`inflate[g2, r2]];
  b = g1[[5]];                       (* == g2[[5]] *)
  (* t = x^b, absorb t^e into g1 *)
  e = (s + 1)/b - 1;
  a1 = g1[[4]]; a2 = g2[[4]];
  h1 = Roach`g[1, g1[[2]] + e, g1[[3]] + e, a1, 1];
  h2 = Roach`g[1, g2[[2]], g2[[3]], a2, 1];
  (* int_0^oo G1(a1 t) G2(a2 t) dt = (1/a1) G(a2/a1 | ...) (Luke 5.6.2 eq. 1) *)
  val = K g1[[1]] g2[[1]]/(Abs[b] a1^e) / a1 *
    MeijerG[{Join[-h1[[3, 1]], h2[[2, 1]]], Join[h2[[2, 2]], -h1[[3, 2]]]},
            {Join[-h1[[2, 1]], h2[[3, 1]]], Join[h2[[3, 2]], -h1[[2, 2]]]}, a2/a1];
  cond = Roach`cond2[h1, h2];
  {val, cond}]

(* ---------------------------------------------------------------- *)
(* One term: C x^s F1 F2                                            *)
(* ---------------------------------------------------------------- *)

(* (k0 + k1 x^q)^e -> k0^e (1 + (k1/k0) x^q)^e when k0 > 0 under the
   assumptions, so that the table's (1 + t)^-a and (1 - t)^(a-1) forms match. *)
Roach`normBinom[Power[B_Plus, e_], x_] /; FreeQ[e, x] && Length[B] == 2 && !FreeQ[B, x] :=
  Module[{k0 = Select[B, FreeQ[#, x] &], k1, cq},
    k1 = B - k0;
    Which[
      k0 === 0 || TrueQ[k0 == 1], Power[B, e],
      TrueQ[Refine[k0 > 0, Roach`$assum]], Roach`pw[k0, e] Power[Expand[1 + k1/k0], e],
      (* (k1 x^q + k0)^e = k1^e x^(q e) (1 + (k0/k1) x^-q)^e for k1 > 0, x > 0 *)
      (cq = Roach`CQ[k1, x]) =!= $Failed && TrueQ[Refine[cq[[1]] > 0, Roach`$assum]],
        Roach`pw[cq[[1]], e] Roach`pw[x, cq[[2]] e] Power[Expand[1 + k0/k1], e],
      True, Power[B, e]]]
Roach`normBinom[e_, x_] := e

(* split a product into constant, power of x, and the other factors *)
Roach`split[f0_, x_] := Module[{f, fs, C = 1, s = 0, rest = {}},
  (* an integer power of a sum with negative powers of x: one fraction first,
     exact for an integer exponent (1/(1/t - t) = t/(1 - t^2)) *)
  f = f0 /. Power[B_Plus, e_Integer] /; !FreeQ[B, Power[x, _?Negative]] :> Together[B]^e;
  f = f /. pw : Power[_Plus, _] :> Roach`normBinom[pw, x];
  f = Expand[f /. Roach`pw -> Power, x];
  fs = If[Head[f] === Times, List @@ f, {f}];
  Do[Which[
      FreeQ[fi, x], C = C fi,
      fi === x, s = s + 1,
      MatchQ[fi, Power[x, _]] && FreeQ[fi[[2]], x], s = s + fi[[2]],
      MatchQ[fi, Power[E, _]], rest = Append[rest, fi],
      True, rest = Append[rest, fi]], {fi, fs}];
  (* E^(a + b) -> E^a E^b; several exponentials -> one *)
  {C, s, rest}]

(* the indicator joins a (1 - t)^a or Log factor, or stands alone *)
Roach`indQ[f_, x_] := MatchQ[f, Power[1 - _, _] | Power[Log[_], _] | Log[_] | Power[1 + _, _]]
Roach`attachInd[rest_, x_] := Module[{r = rest, i, all},
  If[FreeQ[r, Roach`Ind0], Return[{r}]];
  r = DeleteCases[r, Roach`Ind0[x]];
  (* theta^2 = theta: every factor that has an indicator form may carry it *)
  all = If[Count[r, _?(Roach`indQ[#, x] &)] >= 2,
          {If[Roach`indQ[#, x], # Roach`Ind0[x], #] & /@ r}, {}];
  Join[all, Table[If[i == 0, Append[r, Roach`Ind0[x]],
     ReplacePart[r, i -> r[[i]] Roach`Ind0[x]]], {i, 0, Length[r]}]]]

(* Log[x]^k x^s = d^k/ds^k x^s: a power of the log of the variable is
   removed from the term, the term is solved with x^(s + eps), and the result
   is differentiated k times in eps at 0. The validity conditions are open in
   s, so they hold at eps = 0 when they hold in a neighbourhood. This keeps a
   second indicator G off the integrand. A second indicator jumping at the
   same point (sigma = omega) is a case the PBM conditions exclude. *)
Roach`logPower[f_, x_] := Module[{fs = If[Head[f] === Times, List @@ f, {f}], k = 0, rest = {}},
  Do[Which[
      fi === Log[x], k += 1,
      MatchQ[fi, Power[Log[x], _Integer?Positive]], k += fi[[2]],
      True, AppendTo[rest, fi]], {fi, fs}];
  {k, Times @@ rest}]

(* On (0, 1): t^a (-Log t)^k / (1 - s t^m), s = +1 or -1, is a Lerch sum,
     int_0^1 t^a (-Log t)^k t^(m j) dt = k!/(a + 1 + m j)^(k+1),
   so the term is k! m^(-k-1) Phi(s, k+1, (a+1)/m): Zeta[k+1, (a+1)/m] for
   s = 1, and the alternating sum
   2^(-k-1) (Zeta[k+1, (a+1)/(2m)] - Zeta[k+1, (a+1+m)/(2m)]) for s = -1.
   Here 1/(1 - t^m) alone is not integrable at 1, so this is not two G
   functions. *)
Roach`lerchTerm[f_, x_, assum_] := Module[{C, s, rest, k = 0, m, sg, r = {}, a, h},
  If[FreeQ[f, Roach`Ind0[x]], Return[$Failed]];
  {C, s, rest} = Roach`split[f /. Roach`Ind0[x] -> 1, x];
  Do[Which[
      fi === Log[x], k += 1,
      MatchQ[fi, Power[Log[x], _Integer?Positive]], k += fi[[2]],
      MatchQ[fi, Power[1 + c_. Power[x, mm_.], -1] /; MatchQ[c, 1 | -1] && IntegerQ[mm] && mm > 0],
        AppendTo[r, fi],
      True, Return[$Failed, Module]], {fi, rest}];
  If[k < 1 || Length[r] != 1, Return[$Failed]];
  h = r[[1, 1]];                                   (* 1 + c x^m *)
  m = Exponent[h, x]; sg = -Coefficient[h, x, m];
  a = s;
  (* (-Log t)^k = (-1)^k Log[t]^k *)
  {C (-1)^k k! m^(-k - 1) If[sg === 1, Zeta[k + 1, (a + 1)/m],
     2^(-k - 1) (Zeta[k + 1, (a + 1)/(2 m)] - Zeta[k + 1, (a + 1 + m)/(2 m)])],
   Re[a] > -1}]

Roach`term[f_, x_, assum_] := Module[{k, g, r, v},
  r = Roach`term0[f, x, assum];
  If[r =!= $Failed, Return[r]];
  r = Roach`lerchTerm[f, x, assum];
  If[r =!= $Failed, Return[{r[[1]], Roach`decide[r[[2]], assum]}]];
  {k, g} = Roach`logPower[f, x];
  If[k == 0, Return[$Failed]];
  r = Roach`term0[g x^Roach`eps, x, assum];
  If[r === $Failed, Return[$Failed]];
  (* a derivative in a parameter of a pFq / MeijerG is not a closed form *)
  If[!FreeQ[r[[1]], HypergeometricPFQ | MeijerG], Return[$Failed]];
  v = Quiet[D[r[[1]], {Roach`eps, k}] /. Roach`eps -> 0];
  If[!FreeQ[v, Derivative], Return[$Failed]];
  {v, r[[2]] /. Roach`eps -> 0}]

Roach`term0[f_, x_, assum_] := Module[{C, s, rest, groups, gs, res = $Failed, r, val, cond, out},
  Roach`$assum = assum;
  {C, s, rest} = Roach`split[f, x];
  (* merge exponentials *)
  With[{ex = Cases[rest, Power[E, _]]},
    If[Length[ex] > 1, rest = Append[DeleteCases[rest, Power[E, _]], E^Total[ex[[All, 2]]]]]];
  (* an exponential of a sum of x-powers splits into factors *)
  rest = Flatten[rest /. Power[E, p_Plus] /; !FreeQ[p, x] :> (E^# & /@ List @@ p)];
  rest = Flatten[{rest}];
  Catch[
   Do[
     gs = Roach`GForm[#, x] & /@ grp;
     If[MemberQ[gs, $Failed], Continue[]];
     If[!And @@ (Length[#] == 1 & /@ gs), Continue[]];
     gs = First /@ gs;
     Switch[Length[gs],
      0, Continue[],
      1, r = Roach`int1[C, s, gs[[1]], x],
      2, r = Roach`int2[C, s, gs[[1]], gs[[2]], x],
      _, Continue[]];
     If[r === $Failed, Continue[]];
     {val, cond} = r;
     cond = Roach`decide[cond, assum];
     If[cond === False, Continue[]];
     Throw[{val, cond}],
     {grp, Roach`attachInd[rest, x]}];
   $Failed]]

(* ---------------------------------------------------------------- *)
(* (0, oo) driver with linearity and fallbacks                      *)
(* ---------------------------------------------------------------- *)

Roach`expandForms[f_, x_] := DeleteDuplicates[{f, Expand[f],
   Expand[TrigExpand[f]], Expand[f /. {Csch[u_] :> 2/(E^u - E^-u), Sech[u_] :> 2/(E^u + E^-u)}]}]

Roach`int0oo[f_, x_, assum_] := Module[{terms, res, val, cond = True, out = $Failed},
  Catch[
   Do[
     terms = If[Head[g] === Plus, List @@ g, {g}];
     res = Roach`term[#, x, assum] & /@ terms;
     If[MemberQ[res, $Failed], Continue[]];
     val = Total[res[[All, 1]]];
     cond = Roach`decide[And @@ res[[All, 2]], assum];
     If[cond === False, Continue[]];
     Throw[{val, cond}],
     {g, Roach`expandForms[f, x]}];
   $Failed]]

(* t = E^-x on (0, oo): int f(x) dx = int_0^1 f(-Log t) dt/t *)
Roach`expSubst[f_, x_, assum_] := Module[{t, g, r},
  If[FreeQ[f, Power[E, _?(!FreeQ[#, x] &)]] && FreeQ[f, Csch | Sech | Sinh | Cosh | Tanh | Coth], Return[$Failed]];
  g = f /. {Csch[u_] :> 2/(E^u - E^-u), Sech[u_] :> 2/(E^u + E^-u),
            Sinh[u_] :> (E^u - E^-u)/2, Cosh[u_] :> (E^u + E^-u)/2};
  (* E^(c x) becomes t^-c on evaluation. No PowerExpand: it is not valid
     on every branch, and Roach`normBinom does the rewriting that is valid
     for 0 < t < 1. *)
  g = (g /. x -> -Log[t])/t;
  If[!FreeQ[g /. Log[t] -> Roach`L, Log] || !FreeQ[g, Power[E, _]], Return[$Failed]];
  r = Roach`int0oo[Expand[g] Roach`Ind0[t], t, And[assum, 0 < t < 1]];
  r]

(* t = Tan[x] on [0, b], b <= Pi/2 *)
Roach`tanSubst[f_, x_, b_, assum_] := Module[{t, g},
  If[FreeQ[f, Sin | Cos | Tan | Sec | Csc | Cot], Return[$Failed]];
  If[!TrueQ[0 < b <= Pi/2], Return[$Failed]];
  g = f /. {Cos[2 x] -> (1 - t^2)/(1 + t^2), Sin[2 x] -> 2 t/(1 + t^2),
            Cos[x] -> 1/Sqrt[1 + t^2], Sin[x] -> t/Sqrt[1 + t^2], Tan[x] -> t,
            Sec[x] -> Sqrt[1 + t^2], Csc[x] -> Sqrt[1 + t^2]/t, Cot[x] -> 1/t};
  If[!FreeQ[g, x], Return[$Failed]];
  (* all bases are positive for 0 < t < Tan[b] <= Infinity, so PowerExpand is exact here *)
  g = PowerExpand[g/(1 + t^2)];
  Roach`definite[g, t, 0, If[b === Pi/2, Infinity, Tan[b]], And[assum, 0 < t]]]

(* integration by parts: u(x) x^s with u an inverse function of x *)
Roach`ibp[f_, x_, a_, b_, assum_] := Module[{u, v, du, bt, r},
  u = Cases[f, (ArcSin | ArcCos | ArcTan | ArcCot)[_?(!FreeQ[#, x] &)], {0, Infinity}];
  If[Length[u] != 1, Return[$Failed]];
  u = u[[1]];
  v = Together[f/u];
  If[!FreeQ[v, u[[0]]], Return[$Failed]];
  v = Integrate[v, x];
  If[!FreeQ[v, Integrate], Return[$Failed]];
  bt = Limit[u v, x -> b, Direction -> "FromBelow"] - Limit[u v, x -> a, Direction -> "FromAbove"];
  If[!FreeQ[bt, Limit | Indeterminate | DirectedInfinity | ComplexInfinity], Return[$Failed]];
  du = Together[D[u, x] v];
  r = Roach`definite[du, x, a, b, assum];
  If[r === $Failed, $Failed, {bt - r[[1]], r[[2]]}]]

(* Log[c y^q] -> Log[c] + q Log[y] for y > 0 and c provably positive (the
   rescale x -> b y of a finite interval produces Log[b y]) *)
Roach`logSplit[e_, y_, assum_] := e /. Log[u_] :> With[{cq = Roach`CQ[u, y]},
    Log[cq[[1]]] + cq[[2]] Log[y] /; cq =!= $Failed && cq[[1]] =!= 1 &&
      TrueQ[Quiet[Simplify[cq[[1]] > 0, assum]]]]

(* any interval -> {value, condition} or $Failed *)
Roach`definite[f_, x_, a_, b_, assum_] := Module[{y, r, r1, r2},
  Which[
   a === b, {0, True},
   a === 0 && b === Infinity,
     r = Roach`int0oo[f, x, And[assum, x > 0]];
     If[r === $Failed, r = Roach`expSubst[f, x, assum]];
     If[r === $Failed, r = Roach`ibp[f, x, a, b, assum]];
     r,
   a === -Infinity && b === Infinity,
     r1 = Roach`definite[f, x, 0, Infinity, assum];
     If[r1 === $Failed, Return[$Failed]];
     r2 = Roach`definite[f /. x -> -x, x, 0, Infinity, assum];
     If[r2 === $Failed, Return[$Failed]];
     {r1[[1]] + r2[[1]], Roach`decide[r1[[2]] && r2[[2]], assum]},
   a === -Infinity,
     Roach`definite[f /. x -> -x, x, -b, Infinity, assum],
   b === Infinity,
     Roach`definite[f /. x -> x + a, x, 0, Infinity, assum],
   True,
     (* finite: x -> a + (b - a) y on 0 < y < 1 *)
     r = Roach`int0oo[Expand[Roach`logSplit[(b - a) (f /. x -> a + (b - a) y), y,
                                             And[assum, 0 < y]]] Roach`Ind0[y], y,
                      And[assum, 0 < y < 1]];
     If[r === $Failed && a === 0, r = Roach`tanSubst[f, x, b, assum]];
     If[r === $Failed, r = Roach`ibp[f, x, a, b, assum]];
     r]]

(* simplify the G result *)
Roach`finish[v_, assum_] := Module[{r = FunctionExpand[v]},
  r = TimeConstrained[Refine[r, assum], 3, r];
  If[FreeQ[r, MeijerG | HypergeometricPFQ],
    r = TimeConstrained[Simplify[r, assum], 5, r]];
  r]

(* Check an answer numerically at a point of the assumptions. A quadrature
   can itself be wrong on an oscillatory integrand, so the answer passes when
   it agrees with EITHER of two independent evaluations: the default
   NIntegrate, or Gauss-Kronrod over a split of the interval. Each is capped
   in time. The point is generic (Roach`genericPt), never a special value
   such as 1. An answer that is not numeric at the point is rejected. When
   no point is found, or neither quadrature gives a number, the check is
   inconclusive and does not reject. *)
Roach`splitPts[a_, b_] := Which[
  a === -Infinity || b === -Infinity, {a, b},
  b === Infinity, Join[{a}, a + {1, 4, 16, 64, 256}, {Infinity}],
  True, Subdivide[a, b, 4]]
(* a point satisfying the assumptions, away from the special values 0, 1, 2
   at which a lost factor such as Log[z] or z^s would go unseen: generic
   rationals first, FindInstance only when none of them satisfies assum *)
Roach`genericPt[vars_, assum_] := Module[{cands, pt, isInt, gen},
  (* a variable assumed integer gets a generic integer (3, 4, ...), not 1 *)
  isInt = (!FreeQ[assum, Element[#, Integers]]) & /@ vars;
  (* real variables: i + k + 1 + 2/7 + 1/(11 i + 3) (fractional part in (0.28, 0.36]), never an integer or a
     half-integer, distinct per variable; then the negatives; then values in (0, 1) *)
  gen[i_, k_] := i + k + 1 + 2/7 + 1/(11 i + 3);
  cands = Table[Thread[vars -> MapThread[If[#2, 2 + #1 + k, gen[#1, k]] &,
                                         {Range[Length[vars]], isInt}]], {k, 0, 3}];
  cands = Join[cands,
    {Thread[vars -> MapThread[If[#2, -2 - #1, -gen[#1, 0]] &, {Range[Length[vars]], isInt}]],
     Thread[vars -> MapThread[If[#2, 2 + #1, (3 + 2 #1)/(11 + 3 #1)] &, {Range[Length[vars]], isInt}]]}];
  Do[If[TrueQ[Quiet[Simplify[assum /. c]]], Return[c, Module]], {c, cands}];
  pt = Quiet[FindInstance[assum, vars, Reals]];
  If[MatchQ[pt, {{__Rule}, ___}], First[pt], $Failed]]
Roach`agrees[vv_, ni_] := NumericQ[ni] && Abs[vv - ni] < 10^-8 (1 + Abs[ni])
Roach`verify[f_, x_, a_, b_, val_, assum_] := Module[{vars, pt, vv, n1, n2, fa, aa, bb},
  vars = Complement[Cases[{f, val, a, b}, _Symbol?(!ValueQ[#] && Context[#] === "Global`" &), Infinity], {x}];
  pt = If[vars === {}, {}, Roach`genericPt[vars, assum]];
  If[pt === $Failed, Return[True]];
  vv = Quiet[N[val /. pt]];
  (* an answer with no value at a point the assumptions allow (a removable
     Gamma pole at integer n, a pFq N cannot continue) is not usable there *)
  If[!NumericQ[vv], Return[False]];
  fa = f /. pt; aa = a /. pt; bb = b /. pt;
  (* Block[{x}, ...]: the variable is restored even when a time cap aborts
     the quadrature mid-sample *)
  n1 = Block[{x}, Quiet[TimeConstrained[NIntegrate[fa, {x, aa, bb}], 5, $Failed]]];
  If[Roach`agrees[vv, n1], Return[True]];
  n2 = Block[{x}, Quiet[TimeConstrained[
         With[{p = Roach`splitPts[aa, bb]},
           Total[Table[NIntegrate[fa, {x, p[[i]], p[[i + 1]]}], {i, Length[p] - 1}]]], 5, $Failed]]];
  If[Roach`agrees[vv, n2], Return[True]];
  (* both quadratures inconclusive -> do not reject; any numeric disagreement -> reject *)
  !(NumericQ[n1] || NumericQ[n2])]

Roach`Run[f_, x_, a_, b_, assum_] := Module[{r, v, c},
  r = Roach`definite[f, x, a, b, assum];
  If[r === $Failed, Return[$Failed]];
  {v, c} = r;
  v = Roach`finish[v, assum];
  If[!Roach`verify[f, x, a, b, v, assum], Return[$Failed]];
  If[c === True, v, ConditionalExpression[v, c]]]

End[];
$ContextPath = Roach`savedContextPath;
