(* ::Package:: *)

(* ========================================================================== *)
(* Residues.m -- the residues tau_p of a normal (resp. reduced) component c Y. *)
(*                                                                            *)
(* Normal places ([B] Def 2.3, Prop 2.4-2.6): with c = N/(Dn S), Dn the       *)
(* normal and S the special part of the denominator,                          *)
(*     R(z) = res_theta(prim_z(N^m G - z^m (D Dn)^m), Bp),                    *)
(* Bp = Dn/gcd(Dn, G) the branch-free part (a pole of the component at a      *)
(* branch place has residue 0, [T] Ch 5 Sec 4), D the tower derivation.       *)
(* The nonzero roots of R are the residues at the (unramified) places where   *)
(* c Y has a pole.  Resultant criterion (Prop 2.9): if c Y has an elementary   *)
(* integral then R = alpha P with P monic of CONSTANT coefficients;           *)
(* otherwise the integral is not elementary (a certificate).                  *)
(*                                                                            *)
(* Special places.                                                            *)
(*   theta = log(eta), theta = infinity (Def 2.10, Prop 2.11): at a place p of *)
(*     ramification r, a = coefficient of theta^-1 of the Puiseux expansion;  *)
(*     a = v' + c eta'/eta (v in k-bar, c constant) gives tau_p = -r c, and no *)
(*     such (v, c) proves non-elementarity.  For a component, c Y ~ theta^(g/m)*)
(*     (Laurent series): the theta^-1 term exists only when m | g (r = 1, m   *)
(*     places p_zeta, a = zeta gamma a1, gamma^m = lc(G)).                    *)
(*   theta = exp(eta), theta = 0 and infinity (Def 2.21, Prop 2.22): for a    *)
(*     REDUCED f, with a fixed special place pbar, v_p(f) - v_pbar(f) =       *)
(*     c eta' + Sum c_i d_i'/d_i gives tau_p = eps r c (eps = 1 above 0, -1    *)
(*     above infinity).  The value of c Y at a place is a theta^0 coefficient: *)
(*     nonzero only at unramified places (theta not dividing G for 0, m | g   *)
(*     for infinity), where it is zeta gamma w, w in k.                       *)
(* The algebraic constants gamma are supported when gamma^m = const * u^m, u  *)
(* in k (so gamma = const^(1/m) u); otherwise the stage declines.             *)
(* ========================================================================== *)

brResultantR::usage = "brResultantR[T, cv, c, z] is R(z) of [B] Prop 2.6 for the normal component c Y (cv the curve record): a polynomial in z with coefficients in Q[x, lower], or 1 when c Y has no pole at an unramified normal place.";
brResultantR[T_Association, cv_Association, c_, z_Symbol] := Module[
  {th = T["theta"], m = cv["m"], G = cv["G"], t = Together[c], den, pn, Dn, Bp, N, E, cl, cont, prim, R},
  den = Denominator[t];
  pn = Numerator[Together[First[brSplit[T, den]]]];
  Dn = If[FreeQ[pn, th], 1, pn];
  If[Dn === 1, Return[1]];
  Bp = Together[Dn/PolynomialGCD[Dn, Expand[G]]];
  Bp = Numerator[Bp];
  If[FreeQ[Bp, th], Return[1]];
  N = Together[t Dn];
  E = Expand[Numerator[Together[N^m G - z^m brD[T, Dn]^m]]];
  cl = DeleteCases[CoefficientList[E, z], 0];
  cont = If[Length[cl] === 1, First[cl], PolynomialGCD @@ cl];
  prim = Expand[Together[E/cont]];
  R = Expand[Resultant[prim, Bp, th]];
  If[R === 0, throwInternal["brResultantRZero", <|"c" -> c|>]];
  R];

brResultantCriterion::usage = "brResultantCriterion[T, R, z] applies [B] Prop 2.9: it returns {\"Constant\", P} with P(z) the monic polynomial of R (whose coefficients are then constants: free of x and the tower symbols), or {\"NonConstant\", Rmonic} when some coefficient of the monic R is not constant -- a certificate that the integral is not elementary.";
brResultantCriterion[T_Association, R_, z_Symbol] := Module[{vars = Prepend[T["syms"], T["x"]], d, P},
  If[FreeQ[R, z], Return[{"Constant", 1}]];
  d = Exponent[R, z];
  P = Expand[Together[R/Coefficient[R, z, d]]];
  P = Sum[Together[Coefficient[P, z, i]] z^i, {i, 0, d}];
  If[AllTrue[CoefficientList[P, z], !brDependsQ[#, vars] &], {"Constant", P}, {"NonConstant", P}]];

brNonzeroResiduePolynomial::usage = "brNonzeroResiduePolynomial[P, z] is the monic squarefree polynomial over Q (or Q-bar constants) whose roots are the distinct NONZERO roots of P(z); 1 when there are none.";
brNonzeroResiduePolynomial[P_, z_Symbol] := Module[{q = Expand[P], sq},
  If[FreeQ[q, z], Return[1]];
  While[Expand[q /. z -> 0] === 0 && !FreeQ[q, z], q = Expand[Together[q/z]]];
  If[FreeQ[q, z], Return[1]];
  sq = Together[q/PolynomialGCD[q, D[q, z]]];
  Expand[Together[sq/Coefficient[sq, z, Exponent[sq, z]]]]];

brNormalResidues::usage = "brNormalResidues[T, cv, c] computes the residues of the normal component c Y at the normal places. It returns <|\"R\" -> R(z), \"P\" -> the nonzero-residue polynomial (monic, squarefree, rational coefficients), \"z\" -> z|>, or throws NonElementary when R fails the resultant criterion ([B] Prop 2.9) and Decline when the residues are not in Q-bar of a supported form.";
brNormalResidues[T_Association, cv_Association, c_] := Module[{z = Unique["Bronstein`Private`brZ"], R, crit, P},
  R = brResultantR[T, cv, c, z];
  crit = brResultantCriterion[T, R, z];
  If[First[crit] === "NonConstant",
    throwNonElementary["resultant criterion ([B] Prop 2.9): R(z) has non-constant coefficients", <|"R" -> R|>]];
  P = brNonzeroResiduePolynomial[crit[[2]], z];
  If[!FreeQ[P, _Complex | _Real] || (!FreeQ[P, z] && !AllTrue[CoefficientList[P, z], rationalNumberQ]),
    throwDecline["residue polynomial with non-rational coefficients"]];
  <|"R" -> R, "P" -> P, "z" -> z|>];

brGammaForm::usage = "brGammaForm[lc, m] writes gamma = lc^(1/m) as kappa^(1/m) u with kappa a rational constant and u in k, returning {kappa, u}, when lc = kappa u^m; None otherwise.";
brGammaForm[lc_, m_Integer] := Module[{u = brMthPowerCofactor[lc, m], kappa},
  If[u === None, Return[None]];
  kappa = Together[lc/u^m];
  If[!FreeQ[kappa, _Symbol], Return[None]];
  {kappa, u}];

brCoefficientAtInf::usage = "brCoefficientAtInf[T, cv, c, j] is the coefficient w in k of theta^j in c theta^(g/m) Ghat^(1/m) (m | g): the expansion of c Y at a place p_zeta above infinity is zeta gamma Sum_j w_j theta^j.";
brCoefficientAtInf[T_Association, cv_Association, c_, j_Integer] := Module[{inf = brPlaceDataInf[T, cv], gm, lau, rho, d},
  gm = inf["g"]/cv["m"];
  lau = brLaurentInf[T, c, j - gm];
  If[lau === {}, Return[0]];
  d = lau[[1, 1]];
  rho = brRootFactorInfAt[T, cv, Max[0, d + gm - j]];
  Together[Total[(#[[2]] If[#[[1]] + gm - j >= 0, rho[[#[[1]] + gm - j + 1]], 0]) & /@ lau]]];

brCoefficientAtZero::usage = "brCoefficientAtZero[T, cv, c, j] is the coefficient w in k of theta^j in c G0n^(1/m) (theta not dividing G, G0n = G/G(0)): the expansion of c Y at a place p_zeta above 0 is zeta gamma0 Sum_j w_j theta^j, gamma0^m = G(0).";
brCoefficientAtZero[T_Association, cv_Association, c_, j_Integer] := Module[{th = T["theta"], G = cv["G"], G0, lau, v, rho, u},
  G0 = Coefficient[G, th, 0];
  lau = brLaurentZero[T, c, j];
  If[lau === {}, Return[0]];
  v = lau[[1, 1]];
  u = Unique["Bronstein`Private`brU"];
  rho = brBinomialSeries[Together[Rest[CoefficientList[Expand[G], th]]/G0], 1/cv["m"], u, Max[0, j - v]];
  Together[Total[(#[[2]] If[j - #[[1]] >= 0, rho[[j - #[[1]] + 1]], 0]) & /@ lau]]];

brLogInfinityResidues::usage = "brLogInfinityResidues[T, cv, c] is the list of the residues tau_p at the places above theta = infinity of the normal component c Y, theta a logarithm ([B] Def 2.10). It returns <|\"ramified\" -> True, \"tau\" -> 0|> when m does not divide deg G (no theta^-1 term: tau = 0 at every such place), or <|\"ramified\" -> False, \"kappa\" -> kappa, \"u\" -> u, \"lambda\" -> lambda|> meaning tau_(p_zeta) = -zeta kappa^(1/m) lambda at the place where Y/theta^(g/m) -> zeta kappa^(1/m) u. Throws NonElementary when a = v' + c eta'/eta has no solution, Decline when undecided or when lc(G) is not a constant times an m-th power in k.";
brLogInfinityResidues[T_Association, cv_Association, c_] := Module[{inf = brPlaceDataInf[T, cv], gf, a1, li},
  If[!inf["unramified"], Return[<|"ramified" -> True, "tau" -> 0|>]];
  a1 = brCoefficientAtInf[T, cv, c, -1];
  If[a1 === 0, Return[<|"ramified" -> False, "kappa" -> 1, "u" -> 1, "lambda" -> 0|>]];
  gf = brGammaForm[inf["lc"], cv["m"]];
  If[gf === None, throwDecline["residue at infinity over a radical extension of k (lc(G) not an m-th power)"]];
  li = brLimitedIntegrate[brLowerTower[T], Together[gf[[2]] a1], brEta[T]];
  Switch[First[li],
    "NoSolution", throwNonElementary["residue at infinity undefined ([B] Def 2.10, Prop 2.11)"],
    "Unknown", throwDecline["residue at infinity: " <> li[[2]]]];
  <|"ramified" -> False, "kappa" -> gf[[1]], "u" -> gf[[2]], "lambda" -> li[[3]]|>];

brExpValueIntegral::usage = "brExpValueIntegral[T, w] decides whether the element w of k satisfies w = c eta' + Sum c_i d_i'/d_i (c, c_i constants, d_i in k-bar): it returns {\"Solution\", c}, {\"NoSolution\"} (proven) or {\"Unknown\", reason}. The k-part of Int w must be c eta + constant.";
brExpValueIntegral[T_Association, w_] := Module[{TL = brLowerTower[T], r, sp, eta = Last[T["args"]], cc},
  If[Together[w] === 0, Return[{"Solution", 0}]];
  r = brIntegrateInK[TL, w];
  Switch[First[r],
    "NonElementary", Return[{"NoSolution"}],
    "Unknown", Return[r]];
  sp = brSplitLogs[TL, r[[2]]];
  If[sp === None, Return[{"Unknown", "value integral is not k-part plus logarithms"}]];
  cc = Together[brD[TL, sp[[1]]]/brD[TL, eta]];
  If[brDependsQ[cc, Prepend[TL["syms"], T["x"]]], {"NoSolution"}, {"Solution", cc}]];

brExpSpecialResidues::usage = "brExpSpecialResidues[T, cv, c] computes the residues at the special places of the REDUCED component c Y, theta an exponential ([B] Def 2.21). It returns <|\"zero\" -> data0, \"inf\" -> dataInf|>, each either <|\"ramified\" -> True|> (every place there has value 0) or <|\"ramified\" -> False, \"kappa\", \"u\", \"c\" -> c_w|> where the value at the place zeta is zeta kappa^(1/m) u w and Int w = c_w eta + logs; with pbar the place zeta = 1 above 0 (or above infinity when 0 is ramified), tau_p = eps r (c_p - c_pbar). Throws NonElementary when a value difference has no integral of the required form, Decline when undecided or an algebraic constant is unsupported.";
brExpSpecialResidues[T_Association, cv_Association, c_] := Module[{z = brPlaceDataZero[T, cv], inf = brPlaceDataInf[T, cv], one},
  one[unram_, lc_, wfun_] := Module[{w, gf, r},
    If[!unram, Return[<|"ramified" -> True|>, Module]];
    w = wfun[];
    If[w === 0, Return[<|"ramified" -> False, "kappa" -> 1, "u" -> 1, "c" -> 0, "w" -> 0|>, Module]];
    gf = brGammaForm[lc, cv["m"]];
    If[gf === None, throwDecline["special value over a radical extension of k"]];
    r = brExpValueIntegral[T, Together[gf[[2]] w]];
    Switch[First[r],
      "NoSolution", throwNonElementary["residue at a special place undefined ([B] Def 2.21, Prop 2.22)"],
      "Unknown", throwDecline["special residue: " <> r[[2]]]];
    <|"ramified" -> False, "kappa" -> gf[[1]], "u" -> gf[[2]], "c" -> r[[2]], "w" -> w|>];
  <|"zero" -> one[z["unramified"], z["G0"], brCoefficientAtZero[T, cv, c, 0] &],
    "inf" -> one[inf["unramified"], inf["lc"], brCoefficientAtInf[T, cv, c, 0] &]|>];

brSpecialResiduesZeroQ::usage = "brSpecialResiduesZeroQ[T, data] is True when every residue at a special place is zero: data from brLogInfinityResidues (log) or brExpSpecialResidues (exp).";
brSpecialResiduesZeroQ[T_Association, data_Association] := If[T["kind"] === "Log",
  TrueQ[data["ramified"]] || data["lambda"] === 0,
  AllTrue[{data["zero"], data["inf"]}, TrueQ[#["ramified"]] || #["c"] === 0 &]];
