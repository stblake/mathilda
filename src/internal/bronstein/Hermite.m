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
brHermiteK::usage = "brHermiteK[T] is the coefficient field k = Q(x, lower)(K0) of the tower T over the constant field K0 = Q(constants) in use, as fieldFn[{x, lower...}, {g, mu, gamma}] with the constants' representations (constNF): the exact arithmetic for Hermite reduction in theta with radical constants (Phase 4.3b C6).";
brHermiteK[T_Association] := Module[{cf = constNF[], K},
  K = fieldFn[Join[{T["x"]}, T["lower"]], {cf["t"], cf["mu"], cf["theta"]}];
  Join[K, KeyTake[cf, {"consts", "I"}]]];

brHermiteBack::usage = "brHermiteBack[K, e] maps an element of brHermiteK[T] back to the radicals: the generator replaced by gamma, combined over Q(constants).";
brHermiteBack[K_Association, e_] := tragerTogether[e /. K["t"] -> K["theta"]];

brCancelK::usage = "brCancelK[T, e] is the rational function e in theta over k with common theta-factors of numerator and denominator cancelled EXACTLY over k(K0) (brHermiteK, kGCD) when radical constants are in use; tragerTogether[e] otherwise. The combination over Q[g] alone does not cancel a factor that divides only modulo mu, so Hermite's multiplicities never fell (brHermiteIterationBound, Phase 4.3b C6).";
brCancelK[T_Association, e_] := Module[{th = T["theta"], t = tragerTogether[e], K, num, den, g},
  If[!brConstantCoefficientsQ[t] || FreeQ[Denominator[t], th], Return[t]];
  K = brHermiteK[T];
  num = kc[K, Numerator[t]]; den = kc[K, Denominator[t]];
  g = kGCD[K, num, den, th];
  If[FreeQ[g, th], Return[t]];
  brHermiteBack[K, kQuo[K, num, g, th]]/brHermiteBack[K, kQuo[K, den, g, th]]];

brHermiteSolve[T_Association, A_, C_, V_] := Module[{th = T["theta"], tc = tragerTogether[C], num, den, eg, K, sol},
  num = Numerator[tc]; den = Denominator[tc];
  (* with radical constants: the extended gcd and the remainder over
     k = Q(x, lower)(K0) exactly (brHermiteK) *)
  If[brConstantCoefficientsQ[{A, C, V}],
    K = brHermiteK[T];
    If[Exponent[kGCD[K, kc[K, den], kc[K, V], th], th] > 0,
      throwInternal["brHermiteDenominatorNotCoprime", <|"C" -> C, "V" -> V|>]];
    eg = kExtGCD[K, kc[K, num], kc[K, V], th];
    If[!FreeQ[eg[[1]], th] || eg[[1]] === 0,
      throwInternal["brHermiteNonInvertibleCongruence", <|"C" -> C, "V" -> V|>]];
    sol = kRem[K, kc[K, A den eg[[2]]], kc[K, V], th];
    Return[brCollectTheta[T, brHermiteBack[K, sol]], Module]];
  If[Exponent[tragerGCD[den, V], th] > 0,
    throwInternal["brHermiteDenominatorNotCoprime", <|"C" -> C, "V" -> V|>]];
  eg = PolynomialExtendedGCD[num, V, th];
  If[!FreeQ[eg[[1]], th] || eg[[1]] === 0,
    throwInternal["brHermiteNonInvertibleCongruence", <|"C" -> C, "V" -> V|>]];
  brCollectTheta[T, PolynomialRemainder[tragerTogether[A den eg[[2, 1]]/eg[[1]]], V, th]]];

brCollectTheta::usage = "brCollectTheta[T, p] is the polynomial p in theta with each coefficient (an element of k) put through Together.";
brCollectTheta[T_Association, p_] := Module[{th = T["theta"], t = tragerTogether[p], cl},
  If[t === 0, Return[0]];
  If[!FreeQ[Denominator[t], th], throwInternal["brCollectThetaNotPolynomial", <|"p" -> p|>]];
  cl = CoefficientList[Expand[Numerator[t]], th];
  Sum[tragerTogether[cl[[i + 1]]/Denominator[t]] th^i, {i, 0, Length[cl] - 1}]];

brNormalMultipleFactors::usage = "brNormalMultipleFactors[T, D0] is the list {{V, j}, ...} of the normal squarefree factors V (involving theta) of the polynomial D0 with multiplicity j >= 2 (from Risch`SplitSquarefreeFactor).";
brNormalMultipleFactors[T_Association, D0_] := Module[{ns = First[brSplitSquarefree[T, D0]]},
  Select[Table[{ns[[j]], j}, {j, Length[ns]}], !FreeQ[#[[1]], T["theta"]] && #[[2]] >= 2 &]];

brHermiteReduce::usage = "brHermiteReduce[T, cv, c] Hermite-reduces the component c Y (cv the curve record brCurve: Y^m = G, D Y = s Y) at the normal places. It returns <|\"Algebraic\" -> b, \"Remainder\" -> r|> with c Y = D(b Y) + r Y, i.e. c = D b + s b + r, b in k(theta), and the normal part of the denominator of r squarefree.";
brHermiteReduce[T_Association, cv_Association, c0_] := Module[
  {th = T["theta"], s = cv["s"], c = brCancelK[T, c0], b = 0, D0, mf, V, j, U, A, C, Bv, H, it = 0, maxit},
  maxit = 1 + Total[(#[[2]] - 1) Exponent[#[[1]], th] & /@ brNormalMultipleFactors[T, Denominator[c]]];
  While[True,
    If[++it > maxit, throwInternal["brHermiteIterationBound"]];
    D0 = Denominator[c];
    mf = brNormalMultipleFactors[T, D0];
    If[mf === {}, Break[]];
    {V, j} = Last[SortBy[mf, Last]];
    (* the exact quotient over k(K0) when constants are in use: Together over
       Q[g] cannot cancel V = theta - 2^(2/3), which divides D0 only modulo mu *)
    U = If[brConstantCoefficientsQ[{D0, V}],
      With[{K = brHermiteK[T]}, brHermiteBack[K, kQuo[K, kc[K, D0], kc[K, V^j], th]]],
      tragerTogether[D0/V^j]];
    A = Expand[Numerator[c]];
    C = tragerTogether[U (V s - (j - 1) brD[T, V])];
    Bv = brHermiteSolve[T, A, C, V];
    H = tragerTogether[Bv/V^(j - 1)];
    b = brCancelK[T, b + H];
    c = brCancelK[T, c - (brD[T, H] + s H)]];
  <|"Algebraic" -> b, "Remainder" -> c|>];

brNormalDenominatorSquarefreeQ::usage = "brNormalDenominatorSquarefreeQ[T, c] is True when the normal part of the denominator of c (in theta) is squarefree.";
brNormalDenominatorSquarefreeQ[T_Association, c_] := brNormalMultipleFactors[T, Denominator[tragerTogether[c]]] === {};
