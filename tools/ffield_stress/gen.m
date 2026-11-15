(* tools/ffield_stress/gen.m -- random cases for the ffield / polynomial-matrix
   kernels (Phase 4.5), printed for the independent checkers (check.py runs
   PARI/GP and the reference system on them).

   ./Mathilda -file tools/ffield_stress/gen.m        (env FFS_SEED, FFS_N)

   Fields: Q, GF(p), Q(theta), Q(a), Q(a, b), Q(theta)(a), Q(theta)(a, b).
   Operations: QuotientRemainder, ExtendedGCD (with planted common factors),
   HermiteForm (square and tall, with planted unimodular transformations),
   DetAdjugate.  One line per case:
     FFS <id> <field> <op> <mu> <modulus> <parameters> <inputs> <output>
   tab separated, every field in InputForm; output "DECLINED" when the kernel
   left the call unevaluated. *)
ffsSeed = With[{s = Environment["FFS_SEED"]}, If[StringQ[s], ToExpression[s], 20261115]];
ffsN = With[{s = Environment["FFS_N"]}, If[StringQ[s], ToExpression[s], 30]];
SeedRandom[ffsSeed];
ffsFields = {
  <|"name" -> "Q", "opts" -> {}, "p" -> {}, "mu" -> None, "mod" -> 0|>,
  <|"name" -> "GF101", "opts" -> {Modulus -> 101}, "p" -> {}, "mu" -> None, "mod" -> 101|>,
  <|"name" -> "GF10007", "opts" -> {Modulus -> 10007}, "p" -> {}, "mu" -> None, "mod" -> 10007|>,
  <|"name" -> "NF2", "opts" -> {Extension -> {t, t^2 - 2}}, "p" -> {}, "mu" -> t^2 - 2, "mod" -> 0|>,
  <|"name" -> "NF3", "opts" -> {Extension -> {t, t^3 - t - 1}}, "p" -> {}, "mu" -> t^3 - t - 1, "mod" -> 0|>,
  <|"name" -> "Fn1", "opts" -> {"Parameters" -> {a}}, "p" -> {a}, "mu" -> None, "mod" -> 0|>,
  <|"name" -> "Fn2", "opts" -> {"Parameters" -> {a, b}}, "p" -> {a, b}, "mu" -> None, "mod" -> 0|>,
  <|"name" -> "FnX1", "opts" -> {"Parameters" -> {a}, Extension -> {t, t^2 + 1}}, "p" -> {a}, "mu" -> t^2 + 1, "mod" -> 0|>,
  <|"name" -> "FnX3", "opts" -> {"Parameters" -> {a}, Extension -> {t, t^3 - 2}}, "p" -> {a}, "mu" -> t^3 - 2, "mod" -> 0|>,
  <|"name" -> "FnX2p", "opts" -> {"Parameters" -> {a, b}, Extension -> {t, t^2 - 3}}, "p" -> {a, b}, "mu" -> t^2 - 3, "mod" -> 0|>};

ffsRat[] := RandomInteger[{-5, 5}]/RandomChoice[{1, 1, 1, 2, 3}];
(* a random element of the field *)
ffsElem[F_] := Module[{n = If[F["mu"] === None, 1, Exponent[F["mu"], t]], c},
  c = Sum[ffsScalar[F] t^j, {j, 0, n - 1}];
  If[F["mod"] > 0, Mod[c, F["mod"]], c]];
ffsScalar[F_] := If[F["p"] === {}, ffsRat[],
  Module[{num, den},
    num = Sum[ffsRat[] Product[v^RandomInteger[{0, 1}], {v, F["p"]}], {2}];
    den = If[RandomReal[] < 0.6, 1, 1 + Sum[RandomInteger[{-2, 2}] v, {v, F["p"]}]];
    If[den === 0, num, Together[num/den]]]];
ffsPoly[F_, d_] := Sum[ffsElem[F] x^k, {k, 0, d}];
ffsNZ[F_, d_] := Module[{p = ffsPoly[F, d]}, While[Expand[p] === 0 || Exponent[p, x] < d, p = ffsPoly[F, d] + x^d]; p];
ffsRed[F_, e_] := Module[{r = Expand[e]},
  If[F["mu"] =!= None, r = Expand[PolynomialRemainder[r, F["mu"], t]]];
  If[F["mod"] > 0, r = PolynomialMod[r, F["mod"]]];
  r];
ffsMat[F_, r_, c_, d_] := Table[ffsRed[F, ffsPoly[F, RandomInteger[{0, d}]]], {r}, {c}];
(* a unimodular matrix: a product of elementary operations *)
ffsUnimod[F_, n_] := Module[{U = IdentityMatrix[n], i, j},
  Do[{i, j} = RandomSample[Range[n], 2];
     U[[i]] = ffsRed[F, #] & /@ (U[[i]] + ffsPoly[F, RandomInteger[{0, 1}]] U[[j]]), {2 n}];
  U];

SetAttributes[ffsOut, HoldAll];
ffsOut[id_, F_, op_, spec_, ins_, out0_] := With[{out = TimeConstrained[out0, 20, "TIMEOUT"]}, Print["FFS\t", id, "\t", F["name"], "\t", op, "\t",
  ToString[F["mu"], InputForm], "\t", F["mod"], "\t", ToString[F["p"], InputForm], "\t", ToString[ins, InputForm], "\t",
  Which[out === "TIMEOUT", "TIMEOUT",
    MatchQ[Head[out], _Symbol] && StringStartsQ[ToString[Head[out]], "FLINT`"], "DECLINED",
    True, ToString[out, InputForm]]]];

id = 0;
Do[
  Module[{opts = F["opts"], A, B, C, M, U, H, sq},
    (* quotient and remainder *)
    Do[A = ffsRed[F, ffsPoly[F, RandomInteger[{1, 5}]]]; B = ffsRed[F, ffsNZ[F, RandomInteger[{1, 3}]]];
      ffsOut[++id, F, "divrem", F, {A, B}, FLINT`PolynomialQuotientRemainder[A, B, x, Sequence @@ opts]], {ffsN}];
    (* extended gcd, half with a planted common factor *)
    Do[C = If[OddQ[k], ffsRed[F, ffsNZ[F, RandomInteger[{1, 2}]]], 1];
      A = ffsRed[F, C ffsNZ[F, RandomInteger[{1, 3}]]]; B = ffsRed[F, C ffsNZ[F, RandomInteger[{1, 3}]]];
      ffsOut[++id, F, "xgcd", F, {A, B}, FLINT`PolynomialExtendedGCD[A, B, x, Sequence @@ opts]], {k, ffsN}];
    (* Hermite forms: square, tall, and planted U.H *)
    Do[sq = If[F["p"] === {}, RandomInteger[{2, 3}], 2];
      M = Which[
        Mod[k, 3] == 0, ffsMat[F, sq, sq, If[F["p"] === {}, 2, 1]],
        Mod[k, 3] == 1, ffsMat[F, sq + 1, sq, If[F["p"] === {}, 2, 1]],
        True, Module[{T = UpperTriangularize[ffsMat[F, sq, sq, 1]]},
          Do[T[[i, i]] = ffsRed[F, ffsNZ[F, RandomInteger[{0, 2}]]], {i, sq}];
          ffsRed[F, #] & /@ (ffsUnimod[F, sq] . T) // Map[ffsRed[F, #] &, #, {2}] &]];
      ffsOut[++id, F, "hnf", F, {M}, FLINT`PolynomialHermiteForm[M, x, Sequence @@ opts]], {k, ffsN}];
    (* determinant and adjugate *)
    Do[sq = If[F["p"] === {}, RandomInteger[{2, 3}], 2]; M = ffsMat[F, sq, sq, If[F["p"] === {}, 2, 1]];
      ffsOut[++id, F, "detadj", F, {M}, FLINT`PolynomialDetAdjugate[M, x, Sequence @@ opts]], {ffsN}]],
  {F, ffsFields}];
Print["FFSDONE\t", id];
