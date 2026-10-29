(* ::Package:: *)

(* ========================================================================== *)
(* Hermite.m -- algebraic Hermite reduction, [T] Ch 4 Sec 2 pp.44-48;          *)
(*              [B] Sec 2.1-2.2, eq. (12); [T] Appendix A cases 1-2.           *)
(*                                                                            *)
(* f = Sum_i A_i w_i / D.  While D has a squarefree factor V of multiplicity  *)
(* j >= 2 (D = U V^j, gcd(U, V) = 1), find B_i with deg B_i < deg V such that *)
(*     f - (Sum_i B_i w_i / V^(j-1))'                                         *)
(* has V to multiplicity at most j-1 in its denominator.  Because the basis   *)
(* is diagonal (w_i' = s_i w_i) the congruence of [T] eq. (11) decouples:     *)
(*     A_i == B_i U (V s_i - (j-1) V')   (mod V),        one per coordinate.  *)
(* The factor U (V s_i - (j-1) V') is a rational function whose denominator   *)
(* is prime to V, and it is invertible modulo V: at a root of V not dividing  *)
(* G this is (j-1) V' != 0; at a root of a factor p^e of G it is             *)
(* V/p p' (e/m - (j-1)) != 0 since 0 < e < m ([T] p.77, case 2).              *)
(* Termination: every step lowers the largest multiplicity.                   *)
(* ========================================================================== *)

hermiteReduce::usage = "hermiteReduce[B, f] performs Hermite reduction over Q of the element f (coordinates on the basis B, rational functions of x). It returns <|\"Algebraic\" -> v, \"Remainder\" -> r|> with f = v' + r and r having a squarefree common denominator (both as coordinate lists).";
hermiteReduce[B_Association, f0_List] := Module[
  {F = fieldQ[], x = B["x"], m = B["m"], f = f0, v, D0, sq, V, j, U, k, Vp, Acoef, Bcoef, H, it = 0, maxit},
  v = afZero[B];
  (* every step lowers the multiplicity of one squarefree factor of the     *)
  (* initial denominator by one: Sum (e - 1) deg steps at most (an           *)
  (* invariant derived from the input, not a limit)                         *)
  maxit = 1 + Total[(#[[2]] - 1) Exponent[#[[1]], x] & /@
    Select[Rest[FactorSquareFreeList[commonDenominator[f0, x]]], !FreeQ[#[[1]], x] &]];
  While[True,
    If[++it > maxit, throwInternal["HermiteIterationBound"]];
    D0 = commonDenominator[f, x];
    sq = Select[Rest[FactorSquareFreeList[D0]], !FreeQ[#[[1]], x] && #[[2]] >= 2 &];
    If[sq === {}, Break[]];
    {V, j} = Last[SortBy[sq, Last]];
    U = Together[D0/V^j];
    k = j - 1;
    Vp = D[V, x];
    Acoef = Expand[Together[# D0]] & /@ f;
    Bcoef = Table[
      hermiteSolve[Acoef[[i + 1]], Together[U (V B["s"][[i + 1]] - k Vp)], V, x],
      {i, 0, m - 1}];
    H = Together[#/V^k] & /@ Bcoef;
    v = afPlus[F, v, H];
    f = afMinus[F, f, afD[F, B, H]]];
  <|"Algebraic" -> v, "Remainder" -> f|>];

commonDenominator::usage = "commonDenominator[f, x] is the least common multiple of the denominators of the coordinates of f.";
commonDenominator[f_List, x_Symbol] := Module[{dens = Denominator[Together[#]] & /@ f},
  Expand[PolynomialLCM @@ Join[{1}, dens]]];

hermiteSolve::usage = "hermiteSolve[A, C, V, x] solves B C == A (mod V) for the polynomial B of degree < deg V, where C is a rational function whose denominator is prime to V and whose numerator is invertible modulo V. A non-invertible C is an internal error (it contradicts [T] Ch 4 Sec 2).";
hermiteSolve[A_, C_, V_, x_Symbol] := Module[{num = Numerator[C], den = Denominator[C], eg},
  If[PolynomialGCD[den, V] =!= 1 && Exponent[PolynomialGCD[den, V], x] > 0,
    throwInternal["HermiteDenominatorNotCoprime", <|"C" -> C, "V" -> V|>]];
  eg = PolynomialExtendedGCD[num, V, x];
  If[Exponent[eg[[1]], x] > 0,
    throwInternal["HermiteNonInvertibleCongruence", <|"C" -> C, "V" -> V|>]];
  Expand[PolynomialRemainder[Expand[A den eg[[2, 1]]/eg[[1]]], V, x]]];

squarefreeDenominatorQ::usage = "squarefreeDenominatorQ[f, x] is True when the common denominator of the coordinates of f is squarefree.";
squarefreeDenominatorQ[f_List, x_Symbol] := Module[{D0 = commonDenominator[f, x]},
  Exponent[PolynomialGCD[D0, D[D0, x]], x] <= 0];
