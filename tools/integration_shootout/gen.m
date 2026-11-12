(* gen.m -- the integration-shootout corpus (benchmarks/256-integration-shootout).

   Run:  SEED=1 COUNT=40 ./Mathilda -file tools/integration_shootout/gen.m > corpus.tsv

   One line per case, tab separated:
     CASE id class family expect f F points fvals fullform
   - class: R rational | T transcendental | A algebraic | M mixed (radical on top
     of a Log/Exp tower) | B radical at the bottom (kernels over Q(x, y)).
   - expect: "elem" (f = D[F] for the elementary F of column 7, or a theorem),
     "nonelem" (a theorem: Liouville / Chebyshev / a holomorphic differential),
     "unknown" (perturbed and random cases: decided by the runs themselves).
   - f, F: InputForm; F is None when not known.
   - points: three rationals where f is finite and real and every radicand and
     log argument is positive (so every system's principal branches agree there);
     fvals: f at the points to 40 digits.  Every system verifies ITS OWN answer by
     comparing its derivative at the points with fvals.
   - fullform: FullForm[f], parsed by tools/integration_shootout/translate.py into
     each system's syntax.
   COUNT is the number of cases per generated family; the algebraic class also
   takes every case of tests/data/trager_stress_corpus.m. *)
seed = ToExpression[Environment["SEED"]]; count = ToExpression[Environment["COUNT"]];
If[!IntegerQ[seed], seed = 1]; If[!IntegerQ[count], count = 40];
SeedRandom[seed];
here = DirectoryName[$InputFileName];
Get[FileNameJoin[{here, "families.m"}]];

pick[l_List] := RandomChoice[l];
poly[d_] := Sum[r[-3, 3] x^k, {k, 0, d - 1}] + nz[1, 3] x^d;
sqfreePoly[d_] := Module[{p}, While[True, p = poly[d]; If[Exponent[p, x] === d && PolynomialGCD[p, D[p, x]] === 1, Return[p]]]];

(* ---------------- R: rational ---------------- *)
ratFamilies = {
  "rherm" -> Function[{}, Module[{q = sqfreePoly[pick[{2, 3}]]}, poly[pick[{1, 2, 3}]]/q^pick[{1, 2, 3}] + r[0, 2] Log[q]]],
  "rlog" -> Function[{}, nz[-3, 3] Log[sqfreePoly[3]] + r[-2, 2] ArcTan[(nz[1, 3] x + r[-2, 2])/nz[1, 3]] + r[0, 1] x^2]
};

(* ---------------- T: transcendental ---------------- *)
logArg[] := pick[{x, x + 1, x^2 + 1, 2 x + 3}];
transFamilies = {
  "tlog" -> Function[{}, With[{L = Log[logArg[]]},
    pick[{r[1, 3] x L^2 + r[-2, 2] L, L/(L + nz[1, 3]), Log[L + nz[1, 3]] + x L, (x + 1) L^3, L^2/x}]]],
  "texp" -> Function[{}, With[{t = pick[{E^x, E^(2 x), E^(-x), E^(x^2), E^(x/2)}]},
    pick[{x t/(1 + t), Log[1 + nz[1, 3] t], t/(x + nz[1, 3]), (x^2 + r[-2, 2]) t, Log[t + x]}]]],
  "tmixed" -> Function[{}, pick[{Log[x] E^x, E^x/Log[x + 2], x Log[1 + E^x], E^(Log[x]^2), Log[x + E^x] x, E^x Log[x]^2/(1 + x)}]],
  "tnested" -> Function[{}, pick[{Log[Log[x]] Log[x], E^(E^x), Log[1 + Log[x]^2], x Log[Log[x + 1]], E^(E^x) x, Log[Log[x]]^2}]],
  "ttrig" -> Function[{}, pick[{Sin[x]/(1 + Cos[x]), x Tan[x], Log[1 + Sin[x]^2], E^x Sin[x], x ArcTan[x], Tan[x]^2 x, Sin[x]^3 Cos[x]^2}]],
  "tconstk" -> Function[{}, With[{L = Log[nz[2, 3] x + nz[2, 3]], t = E^(nz[1, 2] x + nz[1, 2])}, pick[{x L, L^2, t/(1 + t), Log[1 + t] + L}]]]
};
perturbT[] := pick[{E^x/x, 1/Log[x], E^(x^2), Log[x]/(x + 1), Sin[x]/x, E^x/(x + 1), Log[Log[x]]}];

(* ---------------- A: algebraic, generated ---------------- *)
algFamilies = {
  "agen" -> Function[{}, Module[{n = pick[{2, 2, 3, 4}], y},
    y = poly[pick[{1, 2, 3}]]^(1/n);
    pick[{poly[pick[{0, 1, 2}]] y, nz[-2, 2] Log[poly[1] + r[1, 2] y] + r[0, 1] x y, y/(x + nz[1, 3]) + r[0, 1] y^(n - 1)}]]]
};

(* ---------------- B: radical at the bottom ---------------- *)
botFamilies = {
  "bbottom" -> Function[{}, With[{a = nz[1, 3]}, pick[{Sqrt[x^2 + a] Log[x], E^x Sqrt[1 + x^2], Log[x] Log[x + Sqrt[x^2 + a]],
    x^(1/3) Log[x], Log[1 + Sqrt[x + a]] x, E^x/Sqrt[x + a] + 2 E^x Sqrt[x + a]}]]]
};
perturbB[] := pick[{E^x/Sqrt[1 + x^2], Log[x]/Sqrt[x^2 + 1], Sqrt[x] E^x, Log[x] Sqrt[1 + x^3]}];

(* theorem-backed non-elementary classics *)
classicNonelem = {"T" -> E^(x^2), "T" -> E^x/x, "T" -> 1/Log[x], "T" -> Sin[x]/x, "T" -> Log[x]/(x + 1),
  "T" -> E^(E^x) / x, "T" -> x Tan[x], "T" -> Log[Log[x]], "M" -> Sqrt[Log[x]], "M" -> 1/Sqrt[Log[x]],
  "M" -> Log[1 + E^x]^(1/3)/(1 + Log[1 + E^x]), "A" -> 1/Sqrt[x^3 + 1], "A" -> Sqrt[1 - x^4], "B" -> E^x/Sqrt[x]};

(* ---------------- points and output ---------------- *)
cands = {3/7, 5/4, 7/3, 13/5, 2/9, 11/3, 17/4, 6/5, 29/7, 1/11, 23/3, 31/9, 9/7, 8/3};
okPoint[f_, p_] := Module[{bases, fv},
  bases = Join[Cases[{f}, Power[b_, _Rational] :> b, {0, Infinity}], Cases[{f}, Log[u_] :> u, {0, Infinity}]];
  fv = Quiet[N[f /. x -> p, 30]];
  NumberQ[fv] && Im[fv] == 0 && Abs[fv] < 10^8 &&
    AllTrue[Quiet[N[bases /. x -> p, 30]], NumberQ[#] && Im[#] == 0 && # > 1/100 &]];
points[f_] := Take[Select[cands, okPoint[f, #] &], UpTo[3]];

emitted = <||>;
emit[id_, cls_, fam_, expect_, f0_, F_] := Module[{f = f0, pts, fv},
  If[!FreeQ[f, x] && LeafCount[f] <= 600 && !KeyExistsQ[emitted, f],
    pts = TimeConstrained[points[f], 30, {}];
    If[Length[pts] === 3,
      emitted[f] = True;
      fv = N[f /. x -> #, 40] & /@ pts;
      Print["CASE\t", id, "\t", cls, "\t", fam, "\t", expect, "\t", ToString[InputForm[f]], "\t",
        If[F === None, "None", ToString[InputForm[F]]], "\t", ToString[InputForm[pts]], "\t",
        ToString[InputForm[fv]], "\t", ToString[FullForm[f]]]]]];

genFamily[cls_, fam_, fun_, expect_, pert_] := Do[
  Module[{F = Quiet[fun[]], f},
    f = Quiet[Together[D[F, x]]];
    If[pert =!= None, f = Together[f + pert[]]; emit[StringJoin[fam, "p-", ToString[seed], "-", ToString[i]], cls, fam <> "-perturbed", "unknown", f, None],
      emit[StringJoin[fam, "-", ToString[seed], "-", ToString[i]], cls, fam, expect, f, F]]],
  {i, count}];

Do[genFamily["R", fm[[1]], fm[[2]], "elem", None], {fm, ratFamilies}];
Do[genFamily["T", fm[[1]], fm[[2]], "elem", None], {fm, transFamilies}];
genFamily["T", "tlog", transFamilies[[1, 2]], "unknown", perturbT];
genFamily["T", "texp", transFamilies[[2, 2]], "unknown", perturbT];
Do[genFamily["A", fm[[1]], fm[[2]], "elem", None], {fm, algFamilies}];
Do[genFamily["B", fm[[1]], fm[[2]], "elem", None], {fm, botFamilies}];
genFamily["B", "bbottom", botFamilies[[1, 2]], "unknown", Function[{}, perturbB[]]];
(* M: the Bronstein Shape A families, then random / perturbed radical integrands *)
Do[genFamily["M", fm[[1]], Function[{}, fm[[2]][theta[], RandomChoice[{2, 2, 2, 3, 3, 4}]]], "elem", None], {fm, mixedFamilies}];
Do[Module[{th = theta[], n = RandomChoice[{2, 2, 3}], Y, f},
   Y = lin[th]^(1/n);
   f = (r[-2, 2] + r[-1, 1] x + nz[-2, 2] th)/(r[1, 2] x + RandomChoice[{0, 1}] th^r[1, 2] + nz[1, 3]) Y^RandomChoice[{1, -1}];
   emit[StringJoin["mrandom-", ToString[seed], "-", ToString[i]], "M", "mrandom", "unknown", f, None]],
  {i, count}];
Do[Module[{th = theta[], n = RandomChoice[{2, 2, 3}], fam = RandomChoice[Take[mixedFamilies, 6]], F, Y},
   F = fam[[2]][th, n];
   Y = First[Cases[F, Power[b_, q_Rational] :> b^(1/Denominator[q]), {0, Infinity}], lin[th]^(1/n)];
   emit[StringJoin["mperturb-", ToString[seed], "-", ToString[i]], "M", "mperturb", "unknown",
     Together[D[F, x] + RandomChoice[{x Y, Y/x, th Y, Y/(x + 1)}]], None]],
  {i, count}];
(* B, Shape B on curves of genus >= 1 and cube roots (Bronstein's component and
   curve-monomial methods, v0.6527-28).  Generated AFTER every older family, so
   adding them left the earlier families' random draws, and their cases, as
   they were. *)
shapeBFamilies = {
  "bgenus" -> Function[{}, With[{a = nz[1, 3]}, pick[{E^x Sqrt[1 + x^3], x Log[x] (x^3 + a)^(1/3),
    Log[x]^2 Sqrt[x^4 + a], Log[Log[x] + Sqrt[x^3 + a]], E^(2 x) (x^2 + a)^(1/3)/x,
    Sqrt[x^3 - x + a] Log[1 + x], E^x Log[x] Sqrt[1 + x^3]}]]],
  "bcurve" -> Function[{}, With[{a = nz[1, 3]}, pick[{E^(x Sqrt[x^3 + a])/x, x Log[x + Sqrt[x^3 + a]],
    Log[x + Sqrt[x^3 + a]]^2, Sqrt[x^3 + a] E^Sqrt[x^3 + a], (x + Sqrt[x^4 + a]) Log[x^2 + Sqrt[x^4 + a]],
    (2 + 2 Sqrt[x^3 + 1]) E^(x Sqrt[x^3 + 1])/x}]]]
};
perturbBG[] := pick[{E^x/Sqrt[1 + x^3], Log[x] Sqrt[1 + x^3], Log[x]/Sqrt[x^3 + 2], E^Sqrt[1 + x^3],
  Log[x + Sqrt[1 + x^3]], E^x (x^2 + 1)^(1/3)}];
Do[genFamily["B", fm[[1]], fm[[2]], "elem", None], {fm, shapeBFamilies}];
genFamily["B", "bgenus", shapeBFamilies[[1, 2]], "unknown", Function[{}, perturbBG[]]];
genFamily["B", "bcurve", shapeBFamilies[[2, 2]], "unknown", Function[{}, perturbBG[]]];
(* theorem-labelled classics, and the Trager corpus (first seed only) *)
If[seed === 1,
  Do[emit["classic-" <> ToString[i], classicNonelem[[i, 1]], "classic-nonelem", "nonelem", classicNonelem[[i, 2]], None],
    {i, Length[classicNonelem]}];
  Get[FileNameJoin[{here, "..", "..", "tests", "data", "trager_stress_corpus.m"}]];
  Do[With[{c = TragerStressCorpus[[i]]},
     emit["trager-" <> c[[1]], "A", "trager-" <> c[[2]], Switch[c[[4]], "elem", "elem", "nonelem", "nonelem", _, "unknown"],
       ToExpression[c[[3]]], None]],
    {i, Length[TragerStressCorpus]}]];
