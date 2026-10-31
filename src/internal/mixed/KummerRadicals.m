(* ::Package:: *)

(* ========================================================================== *)
(* KummerRadicals.m -- products of flat radicals, split by Kummer class.      *)
(*                                                                            *)
(* Shared by ParallelMixedTower (src/internal/mixed/ParallelMixed.m) and      *)
(* Trager (src/internal/trager/Pipeline.m): an integrand rational in x and    *)
(* several flat radicals b^e (b rational in x) is written f = Sum g M, M a    *)
(* radical monomial, and each piece g M is integrated through ONE radical P   *)
(* with M = kappa P, kappa a locally constant root of unity.  Both integrators *)
(* handle one radical; neither handles a compositum.                          *)
(* ========================================================================== *)

BeginPackage["KummerRadicals`", {"Mathilda`"}];

KummerSplit::usage = "KummerSplit[f, x] writes f, rational in x and flat radicals b^e (b rational in x, free of radicals), as {{M, g, single}, ...} with f = Sum g M, g rational, M a radical monomial (1 for the rational part); single is True when M is the only monomial of its Kummer class. None when f has fewer than two radicands, a nested radical, or a radical sum in the denominator.";

PowerFreeRoot::usage = "PowerFreeRoot[R, L] is {root, sf, c}: R^(1/L) written as root = c sf^(1/L), sf L-th-power free, equal to R^(1/L) up to a locally constant L-th root of unity.";

KummerIntegrate::usage = "KummerIntegrate[pieces, f, x, integrators] integrates f = Sum g M (pieces from KummerSplit) piece by piece. Each integrator is a function h -> antiderivative | {\"not elementary\", ...} | {\"failed\", ...} of an integrand h in x with one radical; they are tried in order. Returns the antiderivative (each piece verified by PossibleZeroQ over complex x), a {\"not elementary\", ...} certificate, or {\"failed\", reason}.";

Begin["`Private`"];

(* KummerSplit[f, x]
   Each radicand b becomes one symbol s_b = b^(1/q_b), q_b the lcm of its
   exponents' denominators (b^(k/q) = (b^(1/q))^k exactly for the principal
   power), so f is rational in x and the s_b.  The class of a monomial
   prod s_b^n_b is its exponent vector over the irreducible factors p of the
   radicands, mod 1: two monomials are in one class iff their quotient is in
   C(x) up to a locally constant root of unity (squarefree coprime p are
   multiplicatively independent modulo m-th powers in C(x)).  Two different
   monomials in one class (x/(Sqrt[x - 1] Sqrt[x + 1]) and Sqrt[x^2 - 1]: their
   quotient is a sign that changes across the plane) stay separate pieces,
   marked not single. *)
KummerSplit[f_, x_Symbol] := Module[{pows, bases, qs, syms, rules, facs, vecs, g, num, den, dmono, terms, n, mono, key,
    res = <||>, cls = <||>},
  pows = DeleteDuplicates[Cases[{f}, Power[b_, e_Rational] /; ! FreeQ[b, x], {0, Infinity}]];
  bases = DeleteDuplicates[pows[[All, 1]]];
  If[Length[bases] < 2, Return[None]];
  If[AnyTrue[bases, ! FreeQ[#, Power[_?(! FreeQ[#, x] &), _Rational]] &], Return[None]];
  qs = (LCM @@ Denominator[Cases[pows, Power[#, e_] :> e]]) & /@ bases;
  syms = Table[Unique["s"], {Length[bases]}];
  rules = Flatten[MapThread[{Power[#1, e_Rational] :> #3^(e #2)} &, {bases, qs, syms}]];
  g = Together[f /. rules];
  If[! FreeQ[g, Power[_?(! FreeQ[#, x] &), _Rational]], Return[None]];
  facs = DeleteDuplicates[Select[Flatten[(First /@ Rest[FactorList[#]]) & /@
            Flatten[{Numerator[Together[#]], Denominator[Together[#]]} & /@ bases]], ! FreeQ[#, x] &]];
  vecs = MapThread[Function[{b, q}, With[{t = Together[b]},
      (kummerMult[Numerator[t], #, x] - kummerMult[Denominator[t], #, x])/q & /@ facs]], {bases, qs}];
  num = Expand[Numerator[g]]; den = Denominator[g];
  dmono = Times @@ (#^Exponent[den, #] & /@ syms);
  den = Cancel[den/dmono];
  If[! FreeQ[den, Alternatives @@ syms], Return[None]];
  terms = If[Head[num] === Plus, List @@ num, {num}];
  Do[n = Exponent[t, #] - Exponent[dmono, #] & /@ syms;
    mono = Times @@ (syms^n);
    key = Mod[n . vecs, 1];
    cls[key] = Union[Lookup[cls, Key[key], {}], {mono}];
    res[mono] = Lookup[res, Key[mono], 0] + t/Times @@ (syms^(n + (Exponent[dmono, #] & /@ syms))),
    {t, terms}];
  KeyValueMap[Function[{mono, coef},
      {mono /. Thread[syms -> MapThread[#1^(1/#2) &, {bases, qs}]], Together[coef/den],
       Length[SelectFirst[Values[cls], MemberQ[#, mono] &]] == 1}], res]];

(* the multiplicity of the irreducible factor p in the polynomial r *)
kummerMult[r_, p_, x_] := Module[{k = 0, t = r},
  While[Together[PolynomialRemainder[t, p, x]] === 0, t = PolynomialQuotient[t, p, x]; k++]; k];

(* PowerFreeRoot[R, L]: num den^(L-1) = sq^L sf, sf L-th-power-free, so
   R^(1/L) = kappa sq sf^(1/L)/den -- the form whose structure (x^(k-1) g(x^k),
   a conic) an integrator sees. *)
PowerFreeRoot[R_, L_] := Module[{num = Numerator[Together[R]], den = Denominator[Together[R]], sq = 1, sf = 1},
  Do[sq *= fc[[1]]^Quotient[fc[[2]], L]; sf *= fc[[1]]^Mod[fc[[2]], L],
    {fc, FactorSquareFreeList[Expand[num den^(L - 1)]]}];
  {sq sf^(1/L)/den, sf, sq/den}];

(* KummerIntegrate[pieces, f, x, integrators]
   For a piece, L is the least integer with M^L = R rational, and g P with
   P = PowerFreeRoot[R, L] (one radical) is integrated; M = kappa P with kappa
   a locally constant L-th root of unity, so an antiderivative F_P of g P
   becomes one of g M written through M (sf^(1/L) -> M/c, exact since
   (M/c)^L = sf), else (M/P) F_P -- each candidate verified by PossibleZeroQ
   over complex x.  A piece proven non-elementary makes f non-elementary when
   M is alone in its class: the projection onto that Kummer character commutes
   with D and keeps Liouville's form.  In a class with several monomials the
   pieces can cancel, so there it is only a decline. *)
KummerIntegrate[pieces_, f_, x_Symbol, integrators_List] := Module[{total = 0, M, g, single, L, R, pf, P, r, F},
  Do[{M, g, single} = piece;
    L = LCM @@ Append[Denominator /@ Cases[{M}, Power[_, e_Rational] :> e, {0, Infinity}], 1];
    R = Together[M^L];
    If[! FreeQ[R, Power[_?(! FreeQ[#, x] &), _Rational]], Return[{"failed", "Kummer piece: monomial power not rational"}, Module]];
    pf = If[L == 1, {R, 1, R}, PowerFreeRoot[R, L]];
    P = pf[[1]];
    r = {"failed", "no integrator"};
    Catch[Scan[Function[integ,
        r = integ[g P];
        If[! MatchQ[r, {"failed", ___}], Throw[True, "kint"]]], integrators], "kint"];
    If[MatchQ[r, {"not elementary", ___}],
      Return[If[single, r, {"failed", "a Kummer piece is not elementary, but its class has several monomials"}], Module]];
    If[ListQ[r], Return[r, Module]];
    F = SelectFirst[
      {r, r /. Power[s_, e_Rational] /; Together[s - pf[[2]]] === 0 :> (M/pf[[3]])^(e L), (M/P) r},
      TrueQ[Quiet[PossibleZeroQ[D[#, x] - g M, Assumptions -> Element[x, Complexes]]]] &, None];
    If[F === None, Return[{"failed", "verification failed"}, Module]];
    total += F,
    {piece, pieces}];
  total];

End[];
EndPackage[];
