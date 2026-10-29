(* ::Package:: *)

(* ========================================================================== *)
(* Arithmetic.m -- arithmetic in K(x)[Y]/(Y^m - G) on the integral basis.     *)
(*                                                                            *)
(* An element is the list of its coordinates {c_0, ..., c_(m-1)} on the      *)
(* basis w_i = Y^i/d_i of Basis.m, each a rational function over the          *)
(* coefficient field F (Field.m), kept canonical with rc.  An element with    *)
(* polynomial coordinates is INTEGRAL over K[x] (it lies in the integral      *)
(* closure), which is how ideals (Ideal.m) are represented.                   *)
(* ========================================================================== *)

afZero::usage = "afZero[B] is the zero element on the basis B.";
afZero[B_Association] := Table[0, {B["m"]}];

afOne::usage = "afOne[B] is the element 1 = w_0.";
afOne[B_Association] := UnitVector[B["m"], 1];

afUnit::usage = "afUnit[B, i] is the basis element w_i (i = 0, ..., m-1) as a coordinate list.";
afUnit[B_Association, i_Integer] := UnitVector[B["m"], i + 1];

afCanon::usage = "afCanon[F, a] puts every coordinate of a in canonical form over F.";
afCanon[F_Association, a_List] := rc[F, #] & /@ a;

afZeroQ::usage = "afZeroQ[F, a] is True when the element a is zero over F (exact).";
afZeroQ[F_Association, a_List] := AllTrue[a, rcZeroQ[F, #] &];

afPlus::usage = "afPlus[F, a, b] is a + b.";
afPlus[F_Association, a_List, b_List] := MapThread[rc[F, #1 + #2] &, {a, b}];

afMinus::usage = "afMinus[F, a, b] is a - b.";
afMinus[F_Association, a_List, b_List] := MapThread[rc[F, #1 - #2] &, {a, b}];

afScale::usage = "afScale[F, r, a] is r a for r in F(x).";
afScale[F_Association, r_, a_List] := rc[F, r #] & /@ a;

afTimes::usage = "afTimes[F, B, a, b] is the product a b, computed in the basis with w_i w_j = factor w_k (the table B[\"mult\"]).";
afTimes[F_Association, B_Association, a_List, b_List] := Module[{m = B["m"], out, kf},
  out = Table[0, {m}];
  Do[
    If[a[[i + 1]] =!= 0 && b[[j + 1]] =!= 0,
      kf = B["mult"][[i + 1, j + 1]];
      out[[kf[[1]] + 1]] += a[[i + 1]] b[[j + 1]] kf[[2]]],
    {i, 0, m - 1}, {j, 0, m - 1}];
  afCanon[F, out]];

afTimesPoly::usage = "afTimesPoly[F, B, a, b] is afTimes for elements with POLYNOMIAL coordinates (integral elements, e.g. rows of ideals): the multiplication factors being polynomials, the product is formed by expansion and one coefficient reduction per coordinate, without the rational-function normal form.";
afTimesPoly[F_Association, B_Association, a_List, b_List] := Module[{m = B["m"], out, kf},
  out = Table[0, {m}];
  Do[
    If[a[[i + 1]] =!= 0 && b[[j + 1]] =!= 0,
      kf = B["mult"][[i + 1, j + 1]];
      out[[kf[[1]] + 1]] += a[[i + 1]] b[[j + 1]] kf[[2]]],
    {i, 0, m - 1}, {j, 0, m - 1}];
  kc[F, #] & /@ out];

afPower::usage = "afPower[F, B, a, k] is a^k for an integer k >= 0 (binary powering).";
afPower[F_Association, B_Association, a_List, k_Integer?NonNegative] := Module[{r = afOne[B], base = a, e = k},
  While[e > 0,
    If[OddQ[e], r = afTimes[F, B, r, base]];
    e = Quotient[e, 2];
    If[e > 0, base = afTimes[F, B, base, base]]];
  r];

afD::usage = "afD[F, B, a] is the derivative d/dx of a: (c_i' + c_i s_i) on the basis, the basis being diagonal for the derivation.";
afD[F_Association, B_Association, a_List] := Module[{x = B["x"]},
  Table[rc[F, D[a[[i]], x] + a[[i]] B["s"][[i]]], {i, B["m"]}]];

afFromY::usage = "afFromY[F, B, e, Y] is the element given by the expression e, a polynomial in the symbol Y of degree < m with coefficients in F(x): c_i = d_i Coefficient[e, Y, i].";
afFromY[F_Association, B_Association, e_, Y_Symbol] := Module[{ex = Expand[Numerator[Together[e]]], den = Denominator[Together[e]]},
  If[!FreeQ[den, Y], throwInternal["afFromYDenominatorInY"]];
  Table[rc[F, B["d"][[i + 1]] Coefficient[ex, Y, i]/den], {i, 0, B["m"] - 1}]];

afToY::usage = "afToY[B, a, Y] is the expression Sum_i c_i Y^i / d_i of the element a in the symbol Y.";
afToY[B_Association, a_List, Y_Symbol] := Sum[a[[i + 1]] Y^i/B["d"][[i + 1]], {i, 0, B["m"] - 1}];

afMultMatrix::usage = "afMultMatrix[F, B, a] is the m x m matrix whose row k holds the coordinates of w_k a.";
afMultMatrix[F_Association, B_Association, a_List] :=
  Table[afTimes[F, B, afUnit[B, k], a], {k, 0, B["m"] - 1}];

afNorm::usage = "afNorm[F, B, a] is the norm of a from K(x)(Y) to K(x): the determinant of its multiplication matrix.";
afNorm[F_Association, B_Association, a_List] := rc[F, Det[afMultMatrix[F, B, a]]];

afInverse::usage = "afInverse[B, a] is 1/a over Q: the solution c of c.M = e_0 with M the multiplication matrix of a (the ring is a field since Y^m - G is irreducible).";
afInverse[B_Association, a_List] := Module[{M = afMultMatrix[fieldQ[], B, a], L, P, det, adj1},
  (* c . M = e_0 means c is the first row of M^(-1).  Clear denominators   *)
  (* (M = P/L, P polynomial) and use the cofactors of P: the first row of  *)
  (* P^(-1) is (C_11, ..., C_m1)/det P.  This avoids LinearSolve over Q(x), *)
  (* which stalled for m >= 4 on rational coordinates.                     *)
  L = PolynomialLCM @@ Join[{1}, Denominator[Together[#]] & /@ Flatten[M]];
  P = Map[Expand[Together[# L]] &, M, {2}];
  det = Expand[Det[P]];
  If[det === 0, throwInternal["afInverseZeroDivisor"]];
  adj1 = If[B["m"] === 1, {1}, Table[(-1)^(1 + j) Det[Drop[P, {j}, {1}]], {j, B["m"]}]];
  Together[# L/det] & /@ adj1];

afPolynomialQ::usage = "afPolynomialQ[F, a] is True when every coordinate of a is a polynomial in x (a is integral over K[x]).";
afPolynomialQ[F_Association, a_List] := AllTrue[a, FreeQ[Denominator[rc[F, #]], _Symbol] &];

basisMod::usage = "basisMod[B, p] reduces the basis descriptor B (over Q) modulo the prime p: G, d_i, s_i and the multiplication factors are reduced over GF(p). Under good reduction ([T] p.68) the result is the integral basis of the reduced curve.";
basisMod[B_Association, p_Integer] := Module[{F = fieldGF[p], x = B["x"]},
  Join[B, <|
    "G" -> gfMod[B["G"], p],
    "d" -> (gfMod[#, p] & /@ B["d"]),
    "s" -> (rc[F, #] & /@ B["s"]),
    "mult" -> Map[{#[[1]], gfMod[#[[2]], p]} &, B["mult"], {2}],
    "factors" -> Select[FactorList[gfMod[B["G"], p], Modulus -> p], !FreeQ[#[[1]], x] &],
    "p" -> p|>]];
