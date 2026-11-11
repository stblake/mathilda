(* ::Package:: *)

(* ========================================================================== *)
(* Puiseux.m -- explicit expansions at the special places theta = 0 and        *)
(*              theta = infinity of a component c Y, Y^m = G(theta).          *)
(*                                                                            *)
(* For a binomial curve no Newton polygon is needed ([B] "Global approach"):  *)
(*   at infinity, G = lc theta^g Ghat with Ghat = 1 + O(1/theta), so          *)
(*       Y = zeta gamma theta^(g/m) Ghat^(1/m),   gamma^m = lc,               *)
(*   gcd(m, g) places of ramification r = m/gcd(m, g), and Ghat^(1/m) is the  *)
(*   binomial series in 1/theta with coefficients in k;                       *)
(*   at 0, G = theta^e G0 with G0(0) != 0, so Y = zeta gamma0 theta^(e/m)      *)
(*   G0n^(1/m) with gamma0^m = G0(0), G0n = G0/G0(0) = 1 + O(theta).          *)
(* So the expansion of c Y at a place is zeta gamma theta^(g/m) (resp.         *)
(* theta^(e/m)) times a Laurent series in 1/theta (resp. theta) over k, which *)
(* brLaurentInf / brLaurentZero compute to a requested order.                 *)
(* ========================================================================== *)

brLaurentInf::usage = "brLaurentInf[T, e, n] is the list {{d, a_d}, {d-1, a_(d-1)}, ...} of the coefficients of the expansion of the element e of k(theta) at theta = infinity, e = Sum_{j <= d} a_j theta^j, from the leading degree d down to degree n (inclusive); {} for e = 0. Coefficients are in k.";
brLaurentInf[T_Association, e_, n_Integer] := Module[{th = T["theta"], t = Together[e], num, den, d, u, ser, cl, lead},
  If[t === 0, Return[{}]];
  num = Expand[Numerator[t]]; den = Expand[Denominator[t]];
  d = Exponent[num, th] - Exponent[den, th];
  If[n > d, Return[{}]];
  (* e = theta^d * (num/theta^dn)/(den/theta^dd), a power series in u = 1/theta *)
  u = Unique["Bronstein`Private`brU"];
  ser = Together[(Expand[num/th^Exponent[num, th]] /. th -> 1/u)/(Expand[den/th^Exponent[den, th]] /. th -> 1/u)];
  cl = brSeriesCoefficients[ser, u, d - n];
  Table[{d - i, cl[[i + 1]]}, {i, 0, d - n}]];

brLaurentZero::usage = "brLaurentZero[T, e, n] is the list {{v, a_v}, {v+1, a_(v+1)}, ...} of the coefficients of the Laurent expansion of e (in k(theta)) at theta = 0, from its order v up to degree n (inclusive); {} for e = 0.";
brLaurentZero[T_Association, e_, n_Integer] := Module[{th = T["theta"], t = Together[e], num, den, v, ser, cl},
  If[t === 0, Return[{}]];
  num = Expand[Numerator[t]]; den = Expand[Denominator[t]];
  v = brThetaOrder[num, th] - brThetaOrder[den, th];
  If[n < v, Return[{}]];
  ser = Together[Expand[num/th^brThetaOrder[num, th]]/Expand[den/th^brThetaOrder[den, th]]];
  cl = brSeriesCoefficients[ser, th, n - v];
  Table[{v + i, cl[[i + 1]]}, {i, 0, n - v}]];

brThetaOrder::usage = "brThetaOrder[p, th] is the order of the polynomial p at th = 0 (the lowest power of th with a nonzero coefficient); 0 for a nonzero constant.";
brThetaOrder[p_, th_Symbol] := Module[{cl = CoefficientList[Expand[p], th], i = 0},
  (* a plain scan: FirstPosition's default level spec also matches the whole list *)
  If[Together[Expand[p]] === 0, throwInternal["brThetaOrderZero"]];
  While[i < Length[cl] && Together[cl[[i + 1]]] === 0, i++];
  i];

brSeriesCoefficients::usage = "brSeriesCoefficients[r, u, n] is the list of the first n+1 Taylor coefficients at u = 0 of the rational function r of u (regular at 0; coefficients are rational in the other symbols), computed by the power-series division recurrence (no Series call).";
brSeriesCoefficients[r_, u_Symbol, n_Integer] := Module[{t = Together[r], P, Q, p, q, out = {}, s},
  P = Expand[Numerator[t]]; Q = Expand[Denominator[t]];
  p = PadRight[CoefficientList[P, u], n + 1]; q = PadRight[CoefficientList[Q, u], n + 1];
  If[Together[q[[1]]] === 0, throwInternal["brSeriesCoefficientsPole", <|"r" -> r|>]];
  Do[
    s = Together[(p[[i + 1]] - Sum[q[[j + 1]] out[[i - j + 1]], {j, 1, i}])/q[[1]]];
    AppendTo[out, s],
    {i, 0, n}];
  out];

brBinomialSeries::usage = "brBinomialSeries[a, alpha, u, n] is the list of the first n+1 coefficients of (1 + a(u))^alpha at u = 0, where a is a polynomial (or a coefficient list {a_1, a_2, ...}) in u with a(0) = 0 and alpha a rational number.";
brBinomialSeries[a_, alpha_, u_Symbol, n_Integer] := Module[{ac, out = {1}, s},
  (* (1+a)^alpha: c_0 = 1, i c_i = Sum_{j=1..i} (alpha j - (i - j)) a_j c_(i-j) (J.C.P. Miller recurrence) *)
  ac = PadRight[If[ListQ[a], Prepend[a, 0], CoefficientList[Expand[a], u]], n + 1];
  If[Together[ac[[1]]] =!= 0, throwInternal["brBinomialSeriesConstantTerm"]];
  Do[
    s = Together[Sum[(alpha j - (i - j)) ac[[j + 1]] out[[i - j + 1]], {j, 1, i}]/i];
    AppendTo[out, s],
    {i, 1, n}];
  out];

brPlaceDataInf::usage = "brPlaceDataInf[T, cv] describes the places of Y^m = G above theta = infinity: <|\"g\" -> deg G, \"lc\" -> lc(G), \"count\" -> gcd(m, g), \"r\" -> m/gcd(m, g), \"unramified\" -> (m divides g)|>.";
brPlaceDataInf[T_Association, cv_Association] := Module[{th = T["theta"], G = cv["G"], m = cv["m"], g},
  g = Exponent[G, th];
  <|"g" -> g, "lc" -> Coefficient[G, th, g], "count" -> GCD[m, g], "r" -> m/GCD[m, g], "unramified" -> Divisible[g, m]|>];

brPlaceDataZero::usage = "brPlaceDataZero[T, cv] describes the places of Y^m = G above theta = 0: <|\"e\" -> multiplicity of theta in G, \"G0\" -> the coefficient G0(0) of the lowest power, \"count\" -> gcd(m, e), \"r\" -> m/gcd(m, e), \"unramified\" -> (e == 0)|>.";
brPlaceDataZero[T_Association, cv_Association] := Module[{th = T["theta"], G = cv["G"], m = cv["m"], e},
  e = brThetaOrder[G, th];
  <|"e" -> e, "G0" -> Coefficient[G, th, e], "count" -> GCD[m, e], "r" -> m/GCD[m, e], "unramified" -> (e === 0)|>];

brRootFactorInfAt::usage = "brRootFactorInfAt[T, cv, n] is the list of the first n+1 coefficients (in powers of 1/theta) of Ghat^(1/m), Ghat = G/(lc theta^g) = 1 + O(1/theta): Y = zeta gamma theta^(g/m) Sum_i coeffs[[i+1]] theta^(-i).";
brRootFactorInfAt[T_Association, cv_Association, n_Integer] := Module[{th = T["theta"], G = cv["G"], g, lc, cl, u},
  g = Exponent[G, th]; lc = Coefficient[G, th, g];
  cl = Reverse[CoefficientList[Expand[G], th]];   (* lc, G_(g-1), ..., G_0 *)
  u = Unique["Bronstein`Private`brU"];
  brBinomialSeries[Together[Rest[cl]/lc], 1/cv["m"], u, n]];
