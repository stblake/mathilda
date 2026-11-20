(* gen_campaign.m -- benchmark 259, the transcendental decision corpus: the
   families of the C completeness campaign's hunts (tools/bronstein_hunts), as a
   shootout corpus, so that FriCAS and the other systems run on exactly the
   integrands the campaign was built against.

   Run:  COUNT=40 ./Mathilda -file tools/integration_shootout/gen_campaign.m \
           > benchmarks/259-transcendental-decision/corpus.tsv

   Every family has two halves (COUNT each):
   - "<fam>-elem": f = D[F] with F random in the family's kernels, so f has an
     elementary integral (a "not elementary" from any system is a false one) --
     except in "prim", whose kernels are Erf / Ei / li: there F is Liouvillian and
     the half is expect "unknown";
   - "<fam>-rand": random integrands from the same kernels, expect "unknown"
     (decided by the runs themselves).
   The parameter a of the Tanh and tangent hunts is 3/2 here (the points must
   evaluate).  The tan-over-Q(i) and Gaussian families are complex-valued: their
   points need every log argument and radicand positive, not f real.  Columns and
   points as in gen.m (benchmark 256); class T throughout. *)
x =.;
count = ToExpression[Environment["COUNT"] /. $Failed -> "40"];
r[a_, b_] := RandomInteger[{a, b}];
nz[a_, b_] := RandomChoice[{-1, 1}] RandomInteger[{a, b}];
Get["tools/integration_shootout/common.m"];

(* complex families: f may be complex at the points; logs/radicands positive *)
okPointC[f_, p_] := Module[{bases, fv},
  bases = Join[Cases[{f}, Power[b_, _Rational] :> b, {0, Infinity}], Cases[{f}, Log[u_] :> u, {0, Infinity}]];
  fv = Quiet[N[f /. x -> p, 30]];
  NumberQ[fv] && Abs[fv] < 10^8 &&
    AllTrue[Quiet[N[bases /. x -> p, 30]], NumberQ[#] && (Im[#] != 0 || # > 1/100) &]];
emitC[id_, fam_, expect_, f_] := Module[{pts, fv},
  If[!FreeQ[f, x] && LeafCount[f] <= 600 && !KeyExistsQ[emitted, f],
    pts = TimeConstrained[Take[Select[cands, okPointC[f, #] &], UpTo[3]], 30, {}];
    If[Length[pts] === 3,
      emitted[f] = True;
      fv = N[f /. x -> #, 40] & /@ pts;
      Print["CASE\t", id, "\tT\t", fam, "\t", expect, "\t", ToString[InputForm[f]], "\tNone\t",
        ToString[InputForm[pts]], "\t", ToString[InputForm[fv]], "\t", ToString[FullForm[f]]]]]];

a = 3/2;
(* {name, complex?, F generator, f generator}; each a function of no arguments *)
tanhT = {Tanh[x], Tanh[x^2], Tanh[2 x], Tanh[x + 1], Tanh[x]/x};
tanhK[T_] := {E^T, E^(x T), Log[1 - T^2], Log[1 + T], Log[x + T], Log[x], E^x, 1, T, T^2, 1/(1 + T), 1/(1 - T^2), 1/T};
tanT = {Tan[x], Tan[x^2], Tan[x]/x, Tan[2 x], Tan[x + 1]};
tanK[T_] := {E^T, E^(T^2), E^(x T), E^(T/x), Log[1 + T^2], Log[T], Log[x + T], Log[x], E^x, 1, T, T^2, 1/(1 + T), 1/(1 + T^2), 1/T};
taniT = {Tan[x], Tan[x^2], Tan[2 x], Tan[x + 1]};
taniK[T_] := {Log[x + I], Log[T - I], Log[T + I], 1/(T - I), 1/(T + I), E^(I x), E^T, Log[x], 1, T, T^2, 1/(1 + T^2)};
gE = {E^(I x), E^(I (1 + x)), E^(-I (1 + x)), E^(I (2 + x)), E^(2 I x), E^(I + 2 I x), E^x, E^(1 + x)};
gT = {Tan[1 + x], Tan[2 + x], Tan[x], Tanh[1 + x]};
gC = {1, E^I, E^(-I), E, 2 I, 1 + I};
gK := {RandomChoice[gE], RandomChoice[gT], Log[x + I], Log[RandomChoice[gT] + I], x, 1, 1/(1 + RandomChoice[gE]), 1/(RandomChoice[gT] - I)};
limL = {Log[x], Log[1 + x], Log[x^2 + 1], Log[1 + E^x], Log[x + E^x]};
limG = {E^x, E^(2 x), E^(x^2), Tanh[x], Tanh[2 x], Tan[x], Coth[x], 1/(1 + E^x), E^x/(1 + E^(2 x))};
primK = {Erf[x], Erfi[x], ExpIntegralEi[x], LogIntegral[x], Erf[x^2], ExpIntegralEi[2 x], E^(-x^2), E^x, Log[x], x, 1};
depK = {{Log[x], Log[1 + x], Log[x^2 - 1], Log[x - 1], Log[x + 1], Log[2 + 2 x]}, {Log[x], E^x, Log[x E^x], Log[1 + x], Log[(1 + x) E^(2 x)]}};
f12K = {E^x, E^(x^2), E^(x + x^2), E^(2 x + x^2), E^(x - x^2), E^(3 x/2), E^(3 x/2 + x^2), Log[x], Log[1 + E^x], E^(x^2 + x^3), E^(x^3)};

fams = {
  {"tanh", False,
    Function[{}, With[{T = RandomChoice[tanhT]}, With[{ks = tanhK[T]}, Sum[RandomChoice[{-1, 1, 2, 1/2, a}] x^r[0, 1] RandomChoice[ks] RandomChoice[ks]/RandomChoice[{1, x, 1 + x, 1 + T, 2 + T^2, T}], {r[1, 2]}]]]],
    Function[{}, With[{T = RandomChoice[tanhT]}, RandomChoice[{1, x, T, 1 - T^2, 1/x, x T}] RandomChoice[tanhK[T]] RandomChoice[{1, T, x, 1 - T^2}]]]},
  {"tan", False,
    Function[{}, With[{T = RandomChoice[tanT]}, With[{ks = tanK[T]}, Sum[RandomChoice[{-1, 1, 2, 1/2, a}] x^r[0, 1] RandomChoice[ks] RandomChoice[ks]/RandomChoice[{1, x, 1 + x, 1 + T, 1 + T^2, T}], {r[1, 2]}]]]],
    Function[{}, With[{T = RandomChoice[tanT]}, RandomChoice[{1, x, T, 1 + T^2, 1/x, x T}] RandomChoice[tanK[T]] RandomChoice[{1, T, x, 1 + T^2}]]]},
  {"tani", True,
    Function[{}, With[{T = RandomChoice[taniT]}, With[{ks = taniK[T]}, Sum[RandomChoice[{-1, 1, I, -I, 1 + I, 1/2}] x^r[0, 1] RandomChoice[ks] RandomChoice[ks]/RandomChoice[{1, x, x + I, 1 + T, T - I}], {r[1, 2]}]]]],
    Function[{}, With[{T = RandomChoice[taniT]}, RandomChoice[{1, x, I, T - I, 1/(x + I)}] RandomChoice[taniK[T]] RandomChoice[{1, T, I x, T + I}]]]},
  {"gauss", True,
    Function[{}, Sum[RandomChoice[gC] x^r[0, 1] RandomChoice[gK] RandomChoice[gK]/RandomChoice[{1, x, x + I, 1 + RandomChoice[gE]}], {r[1, 2]}]],
    Function[{}, RandomChoice[gC] RandomChoice[{1, x, 1/x, 1/(x + I), RandomChoice[gT]}] RandomChoice[gE] RandomChoice[{1, RandomChoice[gT], RandomChoice[gE], x}]]},
  {"lim", False,
    Function[{}, With[{G = RandomChoice[limG], P = RandomChoice[limL]}, Sum[RandomChoice[{-1, 1, 2, 1/2}] x^r[0, 1] RandomChoice[{P, P^2, 1}] RandomChoice[{G, G^2, Log[1 + G^2], Log[1 + G], 1/(1 + G)}]/RandomChoice[{1, x, 1 + x, 1 + G^2}], {r[1, 2]}]]],
    Function[{}, With[{G = RandomChoice[limG], P = RandomChoice[limL]}, RandomChoice[{1, x, 1/x, 1/(1 + x)}] RandomChoice[{P, P^2}] RandomChoice[{G, G^2, 1/(1 + G), D[G, x], 1 - G^2}]]]},
  {"prim", False,
    Function[{}, Sum[RandomChoice[{-1, 1, 2, 1/2}] x^r[0, 2] RandomChoice[primK] RandomChoice[primK]/RandomChoice[{1, x, 1 + x}], {r[1, 2]}]],
    Function[{}, RandomChoice[{1, x, 1/x, 1/(1 + x)}] RandomChoice[primK] RandomChoice[{E^(-x^2), E^x, Log[x], 1/x, x}]]},
  {"dep", False,
    Function[{}, With[{ks = RandomChoice[depK]}, Sum[RandomChoice[{-1, 1, 2}] x^r[0, 2] RandomChoice[ks]^r[1, 2] RandomChoice[ks], {r[1, 3]}]]],
    Function[{}, With[{ks = RandomChoice[depK]}, RandomChoice[{1, x, 1/x, 1/(1 + x)}] RandomChoice[ks] RandomChoice[Append[ks, 1/x]]]]},
  {"f12", False,
    Function[{}, Sum[RandomChoice[{-1, 1, 2, 1/2}] x^r[0, 2] RandomChoice[f12K]^r[1, 2] RandomChoice[f12K]/RandomChoice[{1, 1 + x, x, 1 + E^x, 1 + E^(x^2)}], {r[1, 3]}]],
    Function[{}, RandomChoice[{1, x, 1/x, 1/(1 + x)}] RandomChoice[f12K] RandomChoice[f12K]/RandomChoice[{1, 1 + E^x, 1 + x}]]}
};

Do[
  Module[{nm = fam[[1]], cx = fam[[2]], k, tries, f, id},
    SeedRandom[25900 + StringLength[nm] 7 + First[ToCharacterCode[nm]]];
    {k, tries} = {0, 0};
    While[k < count && tries < 6 count, tries++;
      f = Quiet[D[fam[[3]][], x]];
      If[f =!= 0 && Head[f] =!= D, id = nm <> "-elem-" <> ToString[k + 1];
        before = Length[emitted];
        With[{ex = If[nm === "prim", "unknown", "elem"]},
          If[cx, emitC[id, nm <> "-elem", ex, f], emit[id, "T", nm <> "-elem", ex, f, None]]];
        If[Length[emitted] > before, k++]]];
    {k, tries} = {0, 0};
    While[k < count && tries < 6 count, tries++;
      f = Quiet[fam[[4]][]];
      If[!FreeQ[f, x], id = nm <> "-rand-" <> ToString[k + 1];
        before = Length[emitted];
        If[cx, emitC[id, nm <> "-rand", "unknown", f], emit[id, "T", nm <> "-rand", "unknown", f, None]];
        If[Length[emitted] > before, k++]]]],
  {fam, fams}];
