(* ::Package:: *)

(* ========================================================================== *)
(* Poly.m -- polynomials in theta over k = Q(x)(t_1..t_r): content, the       *)
(*           normal/special split ([B] Def 1.3, Lemma 1.5-1.8; Bronstein,     *)
(*           Symbolic Integration I, Sec 3.5) and the tower derivation.       *)
(*                                                                            *)
(* Elements of k[theta] are expressions in x, the lower monomials and theta;  *)
(* "content" is the gcd over Q[x, lower] of the theta-coefficients of the     *)
(* cleared numerator, so p = content * primitive with primitive in            *)
(* Q[x, lower][theta] of trivial content.  The splits go through the C        *)
(* builtins Risch`SplitFactor / Risch`SplitSquarefreeFactor with the tower    *)
(* derivation T["rules"] (theta a plain symbol, checked in Phase 1d).          *)
(* ========================================================================== *)

brD::usage = "brD[T, e] is the tower derivation of e (Risch`Derivation with T[\"rules\"]); the same as brTowerD.";
brD[T_Association, e_] := brTowerD[T, e];

brThetaContent::usage = "brThetaContent[T, p] is {c, pp} with p = c pp for the element p of k[theta]: c in k, pp in Q[x, lower][theta] with no non-unit content over Q[x, lower] and positive leading numerical coefficient. p = 0 gives {0, 0}.";
brThetaContent[T_Association, p_] := Module[{th = T["theta"], t = Together[p], num, den, cs, g, pp, lc},
  If[t === 0, Return[{0, 0}]];
  num = Expand[Numerator[t]]; den = Denominator[t];
  If[!PolynomialQ[num, th] || !FreeQ[den, th], throwInternal["brThetaContentNotPolynomial", <|"p" -> p|>]];
  cs = Select[CoefficientList[num, th], # =!= 0 &];
  g = If[Length[cs] === 1, First[cs], PolynomialGCD @@ cs];
  pp = Expand[Together[num/g]];
  lc = brNumericalLC[pp];
  {Together[g lc/den], Expand[pp/lc]}];

brNumericalLC::usage = "brNumericalLC[p] is the numerical coefficient of the leading term of the polynomial p in its variables (sign and rational content normalisation helper).";
brNumericalLC[p_] := Module[{rules = CoefficientRules[Expand[p], Variables[p]]},
  If[rules === {}, 1, Sign[Last[First[rules]]]]];

brSplit::usage = "brSplit[T, p] is {pn, ps}: p = pn ps with ps special and every squarefree factor of pn normal with respect to the tower derivation in theta (Risch`SplitFactor). For a logarithmic theta every polynomial is normal ([B] Lemma 1.7); for an exponential theta the special part is the power of theta.";
brSplit[T_Association, p_] := Module[{r = Risch`SplitFactor[p, T["theta"], T["rules"]]},
  If[!ListQ[r] || Length[r] =!= 2, throwInternal["brSplitFailed", <|"p" -> p|>]];
  r];

brSplitSquarefree::usage = "brSplitSquarefree[T, p] is {{N1, N2, ...}, {S1, S2, ...}}: squarefree factorizations of the normal and special parts of p, p = (N1 N2^2 ...)(S1 S2^2 ...) (Risch`SplitSquarefreeFactor).";
brSplitSquarefree[T_Association, p_] := Module[{r},
  (* over Q(constants) (Phase 4.3b C6) the C builtin returns {{}, {}} for a
     polynomial with radical coefficients -- the whole denominator lost, Hermite
     reduced nothing, and two elementary integrands were certified NOT elementary *)
  If[brConstantCoefficientsQ[p], Return[brSplitSquarefreeNF[T, p]]];
  r = Risch`SplitSquarefreeFactor[p, T["theta"], T["rules"]];
  If[!ListQ[r] || Length[r] =!= 2, throwInternal["brSplitSquarefreeFailed", <|"p" -> p|>]];
  r];

brConstantCoefficientsQ::usage = "brConstantCoefficientsQ[p] is True when radical constants are in use (tragerAlgConstants) and p carries one: the arithmetic must then be over Q(constants) (Phase 4.3b C6).";
brConstantCoefficientsQ[p_] := tragerAlgConstants =!= {} && !FreeQ[p, Power[b_, _Rational] /; rationalNumberQ[b]];

brSplitSquarefreeNF::usage = "brSplitSquarefreeNF[T, p] is brSplitSquarefree over the constant field Q(constants) ([B] SplitSquarefreeFactor): the squarefree factorisation of p over the extension (tragerSquareFreeList), theta-free content dropped, each part q_e of multiplicity e split as q_e = N_e S_e with S_e = gcd_theta(q_e, D q_e) (its special part) and N_e = q_e/S_e. {{N_1, ..., N_E}, {S_1, ..., S_E}} as the C builtin returns.";
brSplitSquarefreeNF[T_Association, p_] := Module[{th = T["theta"], sfl, E0, Ns, Ss},
  sfl = Select[tragerSquareFreeList[Expand[p]], !FreeQ[#[[1]], th] &];
  If[sfl === {}, Return[{{}, {}}]];
  E0 = Max[Last /@ sfl];
  Ns = Table[1, {E0}]; Ss = Table[1, {E0}];
  Do[With[{q = f[[1]], e = f[[2]]},
      With[{g = brSpecialPartNF[T, q]},
        Ss[[e]] = Expand[Ss[[e]] g];
        Ns[[e]] = Expand[Ns[[e]] tragerTogether[q/g]]]],
    {f, sfl}];
  {Ns, Ss}];

brSpecialPartNF::usage = "brSpecialPartNF[T, q] is the special part gcd(q, D q) of the squarefree q in k[theta] over Q(constants), normalised to the theta-dependent part (1 when it is free of theta).";
brSpecialPartNF[T_Association, q_] := Module[{th = T["theta"], dq = Numerator[tragerTogether[brD[T, q]]], g},
  g = tragerGCD[q, dq];
  If[FreeQ[g, th], 1, g]];

brNormalQ::usage = "brNormalQ[T, p] is True when p is normal in theta: gcd(p, D p) has degree 0 in theta.";
brNormalQ[T_Association, p_] := If[brConstantCoefficientsQ[p],
  FreeQ[tragerGCD[p, Numerator[tragerTogether[brD[T, p]]]], T["theta"]],
  TrueQ[Risch`NormalQ[p, T["theta"], T["rules"]]]];

brSpecialQ::usage = "brSpecialQ[T, p] is True when p is special in theta: p divides D p in k[theta]. A polynomial free of theta (a unit of k[theta]) is special (Risch`SpecialQ answers False for it).";
brSpecialQ[T_Association, p_] := FreeQ[Together[p], T["theta"]] ||
  If[brConstantCoefficientsQ[p],
    FreeQ[Denominator[tragerTogether[brD[T, p]/p]], T["theta"]],
    TrueQ[Risch`SpecialQ[p, T["theta"], T["rules"]]]];

brSpecialPlaces::usage = "brSpecialPlaces[T] is the list of the special places of k(theta) ([B] Def 1.4): {\"Infinity\"} for a logarithmic theta, {0, \"Infinity\"} for an exponential theta.";
brSpecialPlaces[T_Association] := If[T["kind"] === "Exp", {0, "Infinity"}, {"Infinity"}];
