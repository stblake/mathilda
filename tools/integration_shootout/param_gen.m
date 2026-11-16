(* param_gen.m -- the parametric corpus of the Phase 4.4 step-6 A/B run
   (benchmark 258, tools/integration_shootout/param_ab.py).

   Read after a data file that sets
     pgCases = {{id, expect, Hold[f], Hold[F] | None}, ...};
   For each case one integer coefficient k (the middle one, never an exponent of
   a power) becomes k + a, a symbolic parameter:
     - an `elem` case with a recorded antiderivative F: F_a = F with k -> k + a
       and f_a = Together[D[F_a, x]], still elementary over Q(a) (expect elem);
     - otherwise f itself is parametrised (expect unknown).
   Three evaluation points are chosen as for the shootout corpus (common.m), at
   the parameter value a = 1/7 that the derivative check uses.  One line per
   case:  PCASE<TAB>id<TAB>expect<TAB>f_a (InputForm)<TAB>points
   Nothing here draws randomness: the corpus is a function of the source. *)
pgA0 = 1/7;
pgCands = {3/7, 5/4, 7/3, 13/5, 2/9, 11/3, 17/4, 6/5, 29/7, 1/11, 23/3, 31/9, 9/7, 8/3};
pgOkPoint[f_, p_] := Module[{bases, fv},
  bases = Join[Cases[{f}, Power[b_, _Rational] :> b, {0, Infinity}], Cases[{f}, Log[u_] :> u, {0, Infinity}]];
  fv = Quiet[N[f /. x -> p, 30]];
  NumberQ[fv] && Im[fv] == 0 && Abs[fv] < 10^8 &&
    AllTrue[Quiet[N[bases /. x -> p, 30]], NumberQ[#] && Im[#] == 0 && # > 1/100 &]];

(* an integer coefficient that is not the exponent of a power *)
pgOkPos[e_, p_List] := Length[p] >= 1 && !(Last[p] === 2 && MatchQ[Extract[e, Most[p]], _Power]);
pgParamize[e_] := Module[{pos},
  pos = Select[Position[e, _Integer, {1, Infinity}, Heads -> False], pgOkPos[e, #] &];
  If[pos === {}, None,
    With[{q = pos[[Ceiling[Length[pos]/2]]]}, ReplacePart[e, q -> Extract[e, q] + a]]]];

pgOne[{id_, expect_, Hold[f_], F0_}] := Module[{Fa, fa, ex = "unknown", pts},
  fa = None;
  If[expect === "elem" && MatchQ[F0, Hold[_]],
    Fa = pgParamize[First[F0]];
    If[Fa =!= None && !FreeQ[Fa, a],
      fa = Quiet[TimeConstrained[Together[D[Fa, x]], 20, None]];
      If[fa =!= None && !FreeQ[fa, a] && LeafCount[fa] <= 800, ex = "elem", fa = None]]];
  If[fa === None, fa = pgParamize[f]];
  If[fa === None || FreeQ[fa, a] || FreeQ[fa, x], Return[Null]];
  pts = Quiet[TimeConstrained[Take[Select[pgCands, pgOkPoint[fa /. a -> pgA0, #] &], UpTo[3]], 30, {}]];
  If[Length[pts] =!= 3, Return[Null]];
  Print["PCASE\t", id, "\t", ex, "\t", ToString[InputForm[fa]], "\t", ToString[InputForm[pts]]]];

Scan[pgOne, pgCases];
