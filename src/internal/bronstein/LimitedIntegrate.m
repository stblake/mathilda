(* ::Package:: *)

(* ========================================================================== *)
(* LimitedIntegrate.m -- integration problems in the coefficient field k:      *)
(*                                                                            *)
(*   brIntegrateInK     Int a dx for a in k, as k-part + new logarithms;       *)
(*   brParamIntegrate   find constants kappa_i and w in k with                *)
(*                      w' = rho + Sum kappa_i sigma_i (the parametric form   *)
(*                      of [T] App A eqs (3)-(5) and [B] (E_0));              *)
(*   brLimitedIntegrate find lambda constant and v in k with a = v' + lambda w *)
(*                      ([B] Def 2.10: a = v' + c eta'/eta).                  *)
(*                                                                            *)
(* Method.  Each piece is integrated by the transcendental Risch algorithm    *)
(* (Integrate`RischTranscendental over the tower k), which is a decision      *)
(* procedure: its answer is an antiderivative, unique up to a constant.  An   *)
(* integral lies in k up to a constant iff its new logarithms cancel.  The    *)
(* logarithms of all pieces are reduced to a basis L_1..L_p independent over  *)
(* k by the structure theorem (Risch`LogReducible: a False answer proves a    *)
(* new monomial), so "w' = rho + Sum kappa_i sigma_i has a solution in k" is  *)
(* the linear system  coeff_rho(L_l) + Sum kappa_i coeff_i(L_l) = 0 over the  *)
(* constants.  A non-elementary rho with elementary sigma_i is a certified    *)
(* "no" (Risch`ElementaryIntegralQ False); a non-elementary sigma_i, an       *)
(* ArcTan/RootSum/complex log, or an undecided integral is "Unknown".         *)
(* ========================================================================== *)

brNonElementaryHeads::usage = "brNonElementaryHeads[] is the pattern of heads whose presence in an answer of the integrator means a non-elementary (or not decided) antiderivative.";
brNonElementaryHeads[] := Integrate | Integrate`RischTranscendental | Erf | Erfi | Erfc | ExpIntegralEi |
  LogIntegral | PolyLog | Gamma | SinIntegral | CosIntegral | FresnelS | FresnelC | Hypergeometric2F1;

brIntegrateInK::usage = "brIntegrateInK[T, a] integrates a, an element of the field of the tower T (written in its symbols), by the transcendental Risch algorithm. It returns {\"Elementary\", G} with G the antiderivative in the symbols of T (new logarithms remain as Log[u], u in the symbols), {\"NonElementary\"} (proven), or {\"Unknown\", reason}.";
brIntegrateInK[T_Association, a_] := Module[{aS, G, q},
  If[Together[a] === 0, Return[{"Elementary", 0}]];
  aS = brToSurface[T, a];
  G = Quiet[Integrate`RischTranscendental[aS, T["x"]]];
  If[Head[G] === Integrate`RischTranscendental || !FreeQ[G, brNonElementaryHeads[]],
    q = Quiet[Risch`ElementaryIntegralQ[aS, T["x"]]];
    Return[If[q === False, {"NonElementary"}, {"Unknown", "integral in k not decided"}]]];
  If[!FreeQ[G, ArcTan | RootSum | Root | _Complex | ArcTanh],
    Return[{"Unknown", "integral in k with a non-rational logarithm"}]];
  {"Elementary", brFromSurface[T, G]}];

brSplitLogs::usage = "brSplitLogs[T, G] writes G = g0 + Sum c_j Log[u_j] with g0 in the field of T, c_j rational constants and u_j depending on x or the symbols. It returns {g0, {{c_1, u_1}, ...}} or None when G is not of that form (a non-linear occurrence of a logarithm, a non-constant coefficient).";
brSplitLogs[T_Association, G_] := Module[{vars = Prepend[T["syms"], T["x"]], logs, L, e, cs, g0},
  logs = DeleteDuplicates[Cases[G, Log[u_] /; brDependsQ[u, vars], {0, Infinity}]];
  If[logs === {}, Return[If[brRationalInQ[G, vars], {Together[G], {}}, None]]];
  L = Table[Unique["Bronstein`Private`brL"], {Length[logs]}];
  e = Expand[G /. Thread[logs -> L]];
  If[!PolynomialQ[e, L] || Max[Exponent[e, #] & /@ L] > 1, Return[None]];
  cs = Coefficient[e, #, 1] & /@ L;
  g0 = Together[e /. Thread[L -> 0]];
  If[!AllTrue[cs, FreeQ[#, Alternatives @@ Join[vars, L]] && FreeQ[#, _Complex] &], Return[None]];
  If[!brRationalInQ[g0, vars], Return[None]];
  {g0, Transpose[{cs, First /@ logs}]}];

brLogBasisReduce::usage = "brLogBasisReduce[T, us] reduces the logarithms Log[u] (u in us, written in the symbols of T) over k: it returns {basis, coords, kpart} where basis is a list of arguments b_l whose logarithms are independent over k (each a new monomial by Risch`LogReducible), coords[[i]] the rational vector with D Log[u_i] = Sum coords[[i, l]] D Log[b_l] + D kpart[[i]], and kpart[[i]] in k (a combination of the tower's monomials). Log identities are used at the level of derivatives only, so no branch is chosen.";
brLogBasisReduce[T_Association, us_List] := Module[{x = T["x"], monos = T["monos"], basis = {}, bspecs = {},
    coords = {}, kpart = {}, spec, r, nm, nb, row, kp, Lsym},
  Do[
    spec = Join[brMonoSpec[monos], bspecs];
    r = If[spec === {}, False, Risch`LogReducible[u, x, spec]];
    If[r === False,
      (* a new independent logarithm *)
      Lsym = Unique["Bronstein`Private`brLB"];
      AppendTo[basis, u];
      AppendTo[bspecs, {Lsym, "Log", Together[Risch`Derivation[u, brRules[x, monos]]/u]}];
      AppendTo[coords, UnitVector[Length[basis], Length[basis]]];
      AppendTo[kpart, 0],
      If[!ListQ[r] || Length[r] =!= Length[spec],
        throwUnimplemented["LimitedIntegrate:logrel", "logarithm relation not decided", <|"u" -> u|>]];
      nm = Length[monos]; nb = Length[basis];
      kp = Sum[r[[i]] If[monos[[i]]["kind"] === "Exp", monos[[i]]["arg"], monos[[i]]["sym"]], {i, nm}];
      AppendTo[coords, Take[r, {nm + 1, nm + nb}]];
      AppendTo[kpart, Together[kp]]],
    {u, us}];
  nb = Length[basis];
  coords = PadRight[#, nb] & /@ coords;
  {basis, coords, kpart}];

brParamIntegrate::usage = "brParamIntegrate[T, rho, sigmas] decides for which constants kappa_i the element rho + Sum kappa_i sigma_i of the field of T has an antiderivative in that field. It returns {\"Solution\", w0, ws, eqs, kappas}: kappas the symbols kappa_i, eqs the linear equations they must satisfy (True-free list; {} when unconstrained), and for every solution the antiderivative is w0 + Sum kappa_i ws[[i]] (up to a constant); or {\"NoSolution\"} (proven) or {\"Unknown\", reason}.";
brParamIntegrate[T_Association, rho_, sigmas_List] := Module[{pieces, ints, splits, us, basis, coords, kpart, red,
    kappas, logc, eqs, w0, ws, nbs, sol},
  pieces = Prepend[sigmas, rho];
  ints = brIntegrateInK[T, #] & /@ pieces;
  Which[
    AnyTrue[Rest[ints], First[#] =!= "Elementary" &], Return[{"Unknown", "a parameter's integral is not elementary or not decided"}],
    First[ints[[1]]] === "NonElementary", Return[{"NoSolution"}],
    First[ints[[1]]] =!= "Elementary", Return[ints[[1]]]];
  splits = brSplitLogs[T, #[[2]]] & /@ ints;
  If[MemberQ[splits, None], Return[{"Unknown", "integral in k is not k-part plus logarithms"}]];
  us = DeleteDuplicates[Flatten[(Last /@ #[[2]]) & /@ splits]];
  {basis, coords, kpart} = brLogBasisReduce[T, us];
  red = Association[Thread[us -> Transpose[{coords, kpart}]]];
  nbs = Length[basis];
  (* per piece: the vector of coefficients on the new-log basis and the k-part *)
  logc = Table[
    Total[Prepend[(#[[1]] red[#[[2]]][[1]]) & /@ splits[[j, 2]], ConstantArray[0, nbs]]], {j, Length[pieces]}];
  w0 = Together[splits[[1, 1]] + Total[(#[[1]] red[#[[2]]][[2]]) & /@ splits[[1, 2]]]];
  ws = Table[Together[splits[[j, 1]] + Total[(#[[1]] red[#[[2]]][[2]]) & /@ splits[[j, 2]]]], {j, 2, Length[pieces]}];
  kappas = Table[Symbol["Bronstein`Private`brKappa" <> ToString[i]], {i, Length[sigmas]}];
  eqs = If[nbs === 0, {}, Table[logc[[1, l]] + Sum[kappas[[i]] logc[[i + 1, l]], {i, Length[sigmas]}] == 0, {l, nbs}]];
  eqs = DeleteCases[eqs, True];
  If[!brLinearConstantSystemQ[eqs, kappas], Return[{"Unknown", "constant system not well formed"}]];
  If[MemberQ[eqs, False], Return[{"NoSolution"}]];
  If[eqs =!= {} && kappas =!= {} && (sol = Solve[eqs, kappas]) === {}, Return[{"NoSolution"}]];
  If[eqs =!= {} && kappas === {}, Return[{"NoSolution"}]];
  {"Solution", w0, ws, eqs, kappas}];

brLimitedIntegrate::usage = "brLimitedIntegrate[T, a, w] finds a constant lambda and v in the field of T with a = v' + lambda w (w in that field, not a derivative there). It returns {\"Solution\", v, lambda}, {\"NoSolution\"} (proven) or {\"Unknown\", reason}.";
brLimitedIntegrate[T_Association, a_, w_] := Module[{r, sol, lam},
  r = brParamIntegrate[T, a, {-w}];
  If[First[r] =!= "Solution", Return[r]];
  sol = If[r[[4]] === {}, {}, Solve[r[[4]], r[[5]]]];
  Which[
    r[[4]] === {}, {"Solution", r[[2]], 0},
    sol === {}, {"NoSolution"},
    True,
      lam = r[[5, 1]] /. First[sol];
      If[!FreeQ[lam, r[[5, 1]]], {"Unknown", "limited integration: lambda undetermined"},
        {"Solution", Together[r[[2]] + lam r[[3, 1]]], lam}]]];
