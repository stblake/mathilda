(* common.m -- the corpus-generation helpers shared by the integration-shootout
   generators (gen.m, benchmark 256; gen_compact.m, benchmark 257): random
   draws, the evaluation points and the one-line-per-case TSV output.  Read
   after families.m (r, nz); defines nothing that draws randomness at load. *)
pick[l_List] := RandomChoice[l];
poly[d_] := Sum[r[-3, 3] x^k, {k, 0, d - 1}] + nz[1, 3] x^d;
sqfreePoly[d_] := Module[{p}, While[True, p = poly[d]; If[Exponent[p, x] === d && PolynomialGCD[p, D[p, x]] === 1, Return[p]]]];

(* ---------------- points and output ---------------- *)
cands = {3/7, 5/4, 7/3, 13/5, 2/9, 11/3, 17/4, 6/5, 29/7, 1/11, 23/3, 31/9, 9/7, 8/3};
okPoint[f_, p_] := Module[{bases, fv},
  bases = Join[Cases[{f}, Power[b_, _Rational] :> b, {0, Infinity}], Cases[{f}, Log[u_] :> u, {0, Infinity}]];
  fv = Quiet[N[f /. x -> p, 30]];
  NumberQ[fv] && Im[fv] == 0 && Abs[fv] < 10^8 &&
    AllTrue[Quiet[N[bases /. x -> p, 30]], NumberQ[#] && Im[#] == 0 && # > 1/100 &]];
points[f_] := Take[Select[cands, okPoint[f, #] &], UpTo[3]];

emitted = <||>;
(* the integrand is emitted without its constant factors (stripConstantsK, below), and the
   reference scaled to match *)
emit[id_, cls_, fam_, expect_, f0_, F0_] := Module[{kf = stripConstantsK[f0], f, F, pts, fv},
  f = kf[[2]]; F = If[F0 === None, None, scDivide[F0, kf[[1]]]];
  If[!FreeQ[f, x] && LeafCount[f] <= 600 && !KeyExistsQ[emitted, f],
    pts = TimeConstrained[points[f], 30, {}];
    If[Length[pts] === 3,
      emitted[f] = True;
      fv = N[f /. x -> #, 40] & /@ pts;
      Print["CASE\t", id, "\t", cls, "\t", fam, "\t", expect, "\t", ToString[InputForm[f]], "\t",
        If[F === None, "None", ToString[InputForm[F]]], "\t", ToString[InputForm[pts]], "\t",
        ToString[InputForm[fv]], "\t", ToString[FullForm[f]]]]]];


(* ---------------- constant factors ----------------
   stripConstants[f]: f with every constant factor removed.  A constant factor of an
   integrand is no part of the integration problem (Integrate[c f, x] is c Integrate[f, x])
   and only clutters the corpus.  Two kinds are removed:
     - the numeric factors of the top-level product (-1/2 in -1/2 (x^2 - 3 x - 2)/Sqrt[g]);
     - the numeric content of a sum that is a factor, or the base of a power that is one
       (Sqrt[4 + 2 E^(-2 x) + 8 x + 4 x^2] = Sqrt[2] Sqrt[2 + E^(-2 x) + 4 x + 2 x^2]), which
       is then a constant factor too.  The content is taken positive (a negative one
       negates the sum term by term), so c^r is the principal power of c > 0 and the
       rewrite holds for every complex x.
   Sums deeper inside (the argument of a Log or an Exp) are not touched: a constant
   there is part of the function. *)
scPull[b_Plus] := Module[{l = FactorTermsList[b], c, h},
  If[Length[l] =!= 2 || !NumericQ[l[[1]]] || l[[1]] === 1, Return[{1, b}]];
  {c, h} = l;
  If[!(Head[h] === Plus || !FreeQ[h, x]), Return[{1, b}]];
  If[TrueQ[c < 0], c = -c; h = If[Head[h] === Plus, -# & /@ h, -h]];
  {c, h}];
(* scFactor[e] -> {constant, rest}: e = constant rest *)
scFactor[b_Plus] := scPull[b];
scFactor[Power[b_Plus, r_?NumericQ]] := With[{l = scPull[b]}, {l[[1]]^r, Power[l[[2]], r]}];
scFactor[e_?NumericQ] := {e, 1};
scFactor[e_] := {1, e};
(* scStep[{k, f}] -> {k', f'} with k f = k' f' *)
scStep[{k_, f_}] := Module[{l = If[Head[f] === Times, scFactor /@ List @@ f, {scFactor[f]}]},
  {k Times @@ l[[All, 1]], Times @@ l[[All, 2]]}];
(* stripConstantsK[f] -> {k, g}: f = k g, g free of constant factors *)
stripConstantsK[f_] := FixedPoint[scStep, {1, f}, 4];
stripConstants[f_] := stripConstantsK[f][[2]];
(* scDivide[F, k]: the reference F/k, divided term by term when F is a sum *)
scDivide[F_, k_] := If[Head[F] === Plus, #/k & /@ F, F/k];
