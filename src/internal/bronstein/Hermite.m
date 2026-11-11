(* ::Package:: *)

(* ========================================================================== *)
(* Hermite.m -- Hermite reduction at the NORMAL places, in the one live       *)
(*              coordinate of a component ([B] pp.132-134, system (S),        *)
(*              Lemma 2.1, Prop 2.2; [T] App A Sec 2.1-2.2 p.76).             *)
(*                                                                            *)
(* The component is c Y with Y^m = G over k(theta), D Y = s Y, s = DG/(mG),   *)
(* c = A/D0.  For a normal squarefree factor V of D0 of multiplicity j >= 2   *)
(* (D0 = U V^j, gcd(U, V) = 1) find B in k[theta], deg B < deg V, with        *)
(*     c Y - D(B Y / V^(j-1))                                                 *)
(* having V to multiplicity at most j-1 in its denominator.  With             *)
(* D(B Y/V^(j-1)) = [V DB - (j-1) B DV + s V B]/V^j Y, system (S) reduces to  *)
(* the single congruence                                                      *)
(*     A == B U (V s - (j-1) DV)   (mod V).                                   *)
(* Invertible (Prop 2.2): at a root of V not dividing G, (j-1) DV != 0 mod V  *)
(* since V is normal; at a root of a factor p^e of G, V s - (j-1) DV ==       *)
(* (V/p) Dp (e/m - (j-1)) != 0 since 0 < e < m.  The special part of D0 (the  *)
(* power of theta, exponential case) is left to the algebraic-part Risch DE.  *)
(* ========================================================================== *)

brHermiteSolve::usage = "brHermiteSolve[T, A, C, V] solves B C == A (mod V) in k[theta] for B of degree < deg V, C a rational function whose denominator is prime to V and whose numerator is invertible modulo V. A non-invertible C is an internal error (it contradicts [B] Prop 2.2).";
brHermiteSolve[T_Association, A_, C_, V_] := Module[{th = T["theta"], tc = Together[C], num, den, eg},
  num = Numerator[tc]; den = Denominator[tc];
  If[Exponent[PolynomialGCD[den, V], th] > 0,
    throwInternal["brHermiteDenominatorNotCoprime", <|"C" -> C, "V" -> V|>]];
  eg = PolynomialExtendedGCD[num, V, th];
  If[!FreeQ[eg[[1]], th] || eg[[1]] === 0,
    throwInternal["brHermiteNonInvertibleCongruence", <|"C" -> C, "V" -> V|>]];
  brCollectTheta[T, PolynomialRemainder[Together[A den eg[[2, 1]]/eg[[1]]], V, th]]];

brCollectTheta::usage = "brCollectTheta[T, p] is the polynomial p in theta with each coefficient (an element of k) put through Together.";
brCollectTheta[T_Association, p_] := Module[{th = T["theta"], t = Together[p], cl},
  If[t === 0, Return[0]];
  If[!FreeQ[Denominator[t], th], throwInternal["brCollectThetaNotPolynomial", <|"p" -> p|>]];
  cl = CoefficientList[Expand[Numerator[t]], th];
  Sum[Together[cl[[i + 1]]/Denominator[t]] th^i, {i, 0, Length[cl] - 1}]];

brNormalMultipleFactors::usage = "brNormalMultipleFactors[T, D0] is the list {{V, j}, ...} of the normal squarefree factors V (involving theta) of the polynomial D0 with multiplicity j >= 2 (from Risch`SplitSquarefreeFactor).";
brNormalMultipleFactors[T_Association, D0_] := Module[{ns = First[brSplitSquarefree[T, D0]]},
  Select[Table[{ns[[j]], j}, {j, Length[ns]}], !FreeQ[#[[1]], T["theta"]] && #[[2]] >= 2 &]];

brHermiteReduce::usage = "brHermiteReduce[T, cv, c] Hermite-reduces the component c Y (cv the curve record brCurve: Y^m = G, D Y = s Y) at the normal places. It returns <|\"Algebraic\" -> b, \"Remainder\" -> r|> with c Y = D(b Y) + r Y, i.e. c = D b + s b + r, b in k(theta), and the normal part of the denominator of r squarefree.";
brHermiteReduce[T_Association, cv_Association, c0_] := Module[
  {th = T["theta"], s = cv["s"], c = Together[c0], b = 0, D0, mf, V, j, U, A, C, Bv, H, it = 0, maxit},
  maxit = 1 + Total[(#[[2]] - 1) Exponent[#[[1]], th] & /@ brNormalMultipleFactors[T, Denominator[c]]];
  While[True,
    If[++it > maxit, throwInternal["brHermiteIterationBound"]];
    D0 = Denominator[c];
    mf = brNormalMultipleFactors[T, D0];
    If[mf === {}, Break[]];
    {V, j} = Last[SortBy[mf, Last]];
    U = Together[D0/V^j];
    A = Expand[Numerator[c]];
    C = Together[U (V s - (j - 1) brD[T, V])];
    Bv = brHermiteSolve[T, A, C, V];
    H = Together[Bv/V^(j - 1)];
    b = Together[b + H];
    c = Together[c - (brD[T, H] + s H)]];
  <|"Algebraic" -> b, "Remainder" -> c|>];

brNormalDenominatorSquarefreeQ::usage = "brNormalDenominatorSquarefreeQ[T, c] is True when the normal part of the denominator of c (in theta) is squarefree.";
brNormalDenominatorSquarefreeQ[T_Association, c_] := brNormalMultipleFactors[T, Denominator[Together[c]]] === {};
