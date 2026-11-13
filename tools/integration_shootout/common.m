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
emit[id_, cls_, fam_, expect_, f0_, F_] := Module[{f = f0, pts, fv},
  If[!FreeQ[f, x] && LeafCount[f] <= 600 && !KeyExistsQ[emitted, f],
    pts = TimeConstrained[points[f], 30, {}];
    If[Length[pts] === 3,
      emitted[f] = True;
      fv = N[f /. x -> #, 40] & /@ pts;
      Print["CASE\t", id, "\t", cls, "\t", fam, "\t", expect, "\t", ToString[InputForm[f]], "\t",
        If[F === None, "None", ToString[InputForm[F]]], "\t", ToString[InputForm[pts]], "\t",
        ToString[InputForm[fv]], "\t", ToString[FullForm[f]]]]]];

