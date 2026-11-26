(* easy.m -- is an integrand in one of the "easy" classes?  (benchmark 262)

   easyClass[f] (f an integrand in x) returns the first class that applies:

     "rational"  f is a rational function of x.
     "linrad"    f is in Q(x, y), y^m = (a x + b)/(c x + d): every radical of f
                 is a power of ONE linear-fractional base (up to a constant
                 factor or a reciprocal).  The substitution u^m = (a x + b)/(c x + d)
                 makes f rational in u.
     "euler"     f is in Q(x, Sqrt[q]), q of degree 2 in x: every radical of f is
                 a half-integer power of ONE quadratic (up to a constant factor or
                 a reciprocal).  An Euler substitution makes f rational.
     "dd"        Integrate[f, x, Method -> "DerivativeDivides"] closes f within
                 ddTimeout seconds with an answer that PossibleZeroQ confirms.
                 That is the classical derivative-divides method: f = c h(u) u'
                 for a kernel u(x), reduced to Integrate[h[u], u].
     "none"      none of these.

   The first three are decided syntactically on f itself.  "linrad" and "euler"
   are reported as "linrad-sub" / "euler-sub" when f also has a Log/Exp (or other
   transcendental) kernel but every radical of f is still a radical of one
   linear-fractional or quadratic base in x: the same substitution then leaves an
   integrand with no radical at all (Log[2 x + 1] Sqrt[x + 3] becomes, with
   u^2 = x + 3, a rational function of u times Log[2 u^2 - 5]), so the radical is
   not what makes such an integral hard.

   ddClosesQ is exported for a caller that has already ruled out the
   syntactic classes. *)
(* 120 s, not 30: the classification runs several Mathilda processes at once, and a
   timeout under that load is a "none".  261's cM-113-232 was classed "none" at v0.6732
   and closes in 0.3 s alone at v0.6750 (a timeout or the older method: not known).  A
   class is only as current as the build that made it: easy.tsv records the version,
   and every corpus is re-classified when DerivativeDivides changes. *)
ddTimeout = 120;

easyTranscendentalHeads = Log | Exp | Sin | Cos | Tan | Cot | Sec | Csc | Sinh | Cosh | Tanh | Coth | Sech | Csch |
  ArcTan | ArcTanh | ArcSin | ArcCos | ArcSinh | ArcCosh | ArcCot | ArcCoth | ArcSec | ArcCsc;

easyAlgebraicQ[f_] := FreeQ[f, easyTranscendentalHeads] && FreeQ[f, Power[b_, e_] /; !FreeQ[e, x]] &&
  FreeQ[f, Power[E, _]];

(* the radical bases of f, each Together-ed: {base, root index} *)
easyRadicals[f_] := DeleteDuplicates[Cases[{f}, Power[b_, q_Rational] /; !FreeQ[b, x] :> {Together[b], Denominator[q]},
  {0, Infinity}]];

(* b2 is b1 up to a constant factor or a reciprocal *)
easySameBase[b1_, b2_] := TrueQ[Quiet[PossibleZeroQ[Together[D[b2/b1, x]]]]] ||
  TrueQ[Quiet[PossibleZeroQ[Together[D[b2 b1, x]]]]];

easyOneBase[rs_] := AllTrue[Rest[rs], easySameBase[rs[[1, 1]], #[[1]]] &];

easyLinFracQ[b_] := Module[{n = Numerator[b], d = Denominator[b]},
  PolynomialQ[n, x] && PolynomialQ[d, x] && Exponent[n, x] <= 1 && Exponent[d, x] <= 1];

(* quadratic up to a square factor of the constant: Sqrt[c q] with q quadratic *)
easyQuadQ[b_] := Module[{n = Numerator[b], d = Denominator[b]},
  PolynomialQ[n, x] && PolynomialQ[d, x] && ((Exponent[n, x] <= 2 && Exponent[d, x] === 0) ||
    (Exponent[n, x] === 0 && Exponent[d, x] <= 2))];

ddClosesQ[f_] := Module[{F},
  F = Quiet[TimeConstrained[Integrate[f, x, Method -> "DerivativeDivides"], ddTimeout, $Aborted]];
  F =!= $Aborted && FreeQ[F, Integrate] &&
    TrueQ[Quiet[TimeConstrained[PossibleZeroQ[D[F, x] - f], ddTimeout, False]]]];

(* ---------------- rationalisation by substitution ----------------
   A radical is removed by a substitution when its radicand is linear-fractional
   in x (any root index m: x = (d u^m - b)/(a - c u^m)) or quadratic with index 2
   (Euler: Sqrt[A x^2 + B x + C] = s x + u, s = Sqrt[A], x = (u^2 - C)/(B - 2 s u)).
   easyRationalise applies such substitutions to an innermost radical (one whose
   radicand has no radical in it) and repeats, so a nested or paired radical that
   a chain of these substitutions removes -- Sqrt[x + Sqrt[x^2 + c]],
   Sqrt[x + 1]/(Sqrt[x + 1] + Sqrt[x]) -- is caught as well.  The integrand is
   carried as g(u) x'(u), renamed back to x after each step.  Every radical with
   the radicand (up to a constant factor or a reciprocal) is rewritten at once,
   by Factor and PowerExpand of the substituted radicand: PowerExpand's branch
   choice does not matter here, only whether a radical is left. *)
easyInnermost[f_] := Select[easyRadicals[f], FreeQ[#[[1]], Power[b_, _Rational] /; !FreeQ[b, x]] &];

easyParam[b_, m_] := Module[{n = Numerator[b], d = Denominator[b], t, cn, cd, s},
  If[easyLinFracQ[b],
    {cn, cd} = {Coefficient[n, x, 1], Coefficient[d, x, 1]};
    (* (cn x + n0)/(cd x + d0) = t^m *)
    Return[Together[(t^m (d /. x -> 0) - (n /. x -> 0))/(cn - cd t^m)] /. t -> x]];
  If[m === 2 && easyQuadQ[b] && Exponent[d, x] === 0,
    {cn, cd, s} = {Coefficient[n, x, 1], n /. x -> 0, Sqrt[Coefficient[n, x, 2]]};
    If[s === 0, Return[$Failed]];
    Return[Together[(t^2 - cd)/(cn - 2 s t)] /. t -> x]];
  $Failed];

(* bottom-up, so an inner radical is collapsed before its enclosing radicand is
   factored (a single ReplaceAll would stop at the outer one) *)
easyCollapse[t_][Power[b_, q_Rational]] /; !FreeQ[b, t] := PowerExpand[Power[Factor[Together[b]], q]];
easyCollapse[t_][e_] := e;
easySubstitute[f_, phi_] := Module[{t, g},
  g = (f /. x -> (phi /. x -> t)) D[phi /. x -> t, t];
  g = MapAll[easyCollapse[t], g];
  g /. t -> x];

(* True when a chain of at most `depth` substitutions leaves no radical in x *)
easyRationalisableQ[f_, depth_: 4] := Module[{rs = easyInnermost[f], cands},
  If[easyRadicals[f] === {}, Return[True]];
  If[depth <= 0 || !AllTrue[rs, easyAlgebraicQ[#[[1]]] &], Return[False]];
  cands = DeleteCases[Function[r, Module[{grp = Select[easyRadicals[f], easySameBase[r[[1]], #[[1]]] &]},
      easyParam[r[[1]], LCM @@ (Last /@ grp)]]] /@ rs, $Failed];
  AnyTrue[cands, Function[phi, Module[{g = Quiet[TimeConstrained[easySubstitute[f, phi], 10, $Failed]]},
      g =!= $Failed && easyRationalisableQ[g, depth - 1]]]]];

(* a radicand with a repeated factor (Sqrt[(1 + x^2)^2]): a degenerate draw *)
easyDegenerateQ[f_] := AnyTrue[easyRadicals[f], Function[r, Module[{n = Numerator[r[[1]]], d = Denominator[r[[1]]]},
  (PolynomialQ[n, x] && Exponent[PolynomialGCD[n, D[n, x]], x] > 0) ||
  (PolynomialQ[d, x] && Exponent[PolynomialGCD[d, D[d, x]], x] > 0)]]];

(* the syntactic classes only: "rational", "degenerate", "logradical", "linrad", "euler" (a
   single radicand), "nested" (a chain of those substitutions), each with "-sub"
   when f also has a transcendental kernel; or "none" *)
easySyntacticClass[f_] := Module[{rs = easyRadicals[f], alg = easyAlgebraicQ[f], sfx},
  If[alg && rs === {}, Return["rational"]];
  If[easyDegenerateQ[f], Return["degenerate"]];
  (* every radical sits in a Log of a product of powers: Log[g^p h^q] is
     p Log[g] + q Log[h] up to a constant, so the radical is not one *)
  If[rs =!= {} && easyRadicals[f /. Log[a_] :> PowerExpand[Log[Factor[a]]]] === {}, Return["logradical"]];
  sfx = If[alg, "", "-sub"];
  If[rs =!= {} && AllTrue[rs, easyAlgebraicQ[#[[1]]] &] && easyOneBase[rs],
    If[easyLinFracQ[rs[[1, 1]]], Return["linrad" <> sfx]];
    If[AllTrue[rs, #[[2]] === 2 &] && easyQuadQ[rs[[1, 1]]], Return["euler" <> sfx]]];
  If[rs =!= {} && TrueQ[Quiet[easyRationalisableQ[f]]], Return["nested" <> sfx]];
  "none"];

easyClass[f_] := Module[{c = easySyntacticClass[f]},
  If[c =!= "none", c, If[ddClosesQ[f], "dd", "none"]]];
