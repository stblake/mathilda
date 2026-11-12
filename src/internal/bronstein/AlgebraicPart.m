(* ::Package:: *)

(* ========================================================================== *)
(* AlgebraicPart.m -- the algebraic part of a component: D b + s b = c over   *)
(*                    k(theta) ([B] Prop 2.14-2.15, system (S_log2) for an    *)
(*                    n-th root; [T] App A eqs (3)-(6)).                      *)
(*                                                                            *)
(* For the component c Y (Y^m = G, D Y = s Y) the algebraic part is b Y with  *)
(* D(b Y) = (D b + s b) Y, so it is the Risch DE  D b + s b = c  over k(theta),*)
(* whose solution is unique (the homogeneous solution G^(-1/m) is not in      *)
(* k(theta) for a reduced curve).                                             *)
(*                                                                            *)
(* 1. Risch`RischDECertified over the whole tower (brRischDE): a solution or  *)
(*    a proof of none.                                                        *)
(* 2. When it is undecided and theta is a LOGARITHM, the decision procedure   *)
(*    of [B] Prop 2.15 / [T] App A eqs (3)-(5): after Hermite reduction and   *)
(*    removal of the logarithmic part, b is a POLYNOMIAL in theta, of degree  *)
(*    at most N = floor(max(deg c + 1, 1 - deg G/m)) (Prop 2.14: deg(b Y) <=  *)
(*    max(deg f + 1, 1) at the infinite places).  The coefficients b_N..b_0   *)
(*    are found top-down: the theta-degree j + deg G + deg U equation of      *)
(*        m G U D b + U (D G) b = m G A,     c = A/U,                         *)
(*    is the Risch DE  b_j' + beta b_j = rho_j in k, beta = lcG'/(m lcG),      *)
(*    rho_j depending on b_N..b_(j+1).  When the homogeneous equation has a   *)
(*    solution h in k (lcG = constant * u^m, h = 1/u), b_j = h (w_j + kappa_j) *)
(*    with a free constant kappa_j, and w_j' = rho_j/h is a parametric         *)
(*    integration problem in k (brParamIntegrate); the lower equations fix    *)
(*    the constants (a linear system over the constant field).                *)
(* ========================================================================== *)

brAlgebraicPart::usage = "brAlgebraicPart[T, cv, c] solves D b + s b = c for b in k(theta) (cv the curve record of Y^m = G, s = cv[\"s\"]). It returns {\"Solution\", b} (b verified exactly), {\"NoSolution\"} (proven: the component c Y is not D of an element b Y) or {\"Unknown\", reason}. c must be the normal remainder of Hermite reduction minus the log part (the precondition of [B] Prop 2.15).";
brAlgebraicPart[T_Association, cv_Association, c_] := Module[{s = cv["s"], r, b, own},
  If[Together[c] === 0, Return[{"Solution", 0}]];
  r = brRischDE[brToSurface[T, s], brToSurface[T, c], T["x"], T["kernels"]];
  Which[
    First[r] === "Solution",
      b = Together[brFromSurface[T, r[[2]]]];
      If[brAlgebraicResidualZeroQ[T, s, c, b], {"Solution", b},
        throwInternal["brAlgebraicPartConversion", <|"b" -> r[[2]]|>]],
    First[r] === "NoSolution", r,
    T["kind"] === "Log",
      own = brPolyRDELog[T, cv, c];
      If[First[own] === "Unknown", {"Unknown", r[[2]] <> "; " <> own[[2]]}, own],
    True, r]];

brAlgebraicResidualZeroQ::usage = "brAlgebraicResidualZeroQ[T, s, c, b] is True when D b + s b - c is exactly zero in k(theta).";
brAlgebraicResidualZeroQ[T_Association, s_, c_, b_] := Together[brD[T, b] + s b - c] === 0;

brMthPowerCofactor::usage = "brMthPowerCofactor[e, m] is u when the rational function e (over Q, in any variables) is a rational constant times u^m, u a rational function; None otherwise.";
brMthPowerCofactor[e_, m_Integer] := Module[{t = Together[e], fl},
  fl = Join[FactorList[Numerator[t]], {#[[1]], -#[[2]]} & /@ FactorList[Denominator[t]]];
  fl = Select[fl, Variables[#[[1]]] =!= {} &];
  If[AllTrue[fl, Divisible[#[[2]], m] &], Times @@ (#[[1]]^(#[[2]]/m) & /@ fl), None]];

brKSolveRDE::usage = "brKSolveRDE[TL, beta, rho] solves y' + beta y = rho for y in the field k of the lower tower TL (brRischDE over its kernels). Returns brRischDE's answer with the solution written in the symbols of TL.";
brKSolveRDE[TL_Association, beta_, rho_] := Module[{r},
  r = brRischDE[brToSurface[TL, beta], brToSurface[TL, rho], TL["x"], TL["kernels"]];
  If[First[r] === "Solution", {"Solution", Together[brFromSurface[TL, r[[2]]]]}, r]];

brPolyRDELog::usage = "brPolyRDELog[T, cv, c] decides D b + s b = c for b in k[theta] when theta is a logarithm, by the top-down method of [T] App A eqs (3)-(5) (see the header). Returns {\"Solution\", b}, {\"NoSolution\"} (proven: no polynomial solution, hence none at all for a normal c with the log part removed, [B] Prop 2.15) or {\"Unknown\", reason}.";
brPolyRDELog[T_Association, cv_Association, c_] := Module[
  {th = T["theta"], x = T["x"], m = cv["m"], G = cv["G"], s = cv["s"], TL, vars, tc, A, U, g, du, dc, N,
   lcG, lcU, beta, u, h, Dth, bs, dbs, P, cl, known = {}, kappas = {}, cons = {}, nk = 0, e, eq, rho, rho0,
   sig, alpha, sol, val, pr, w, rest, eqs, sols, B},
  TL = brLowerTower[T];
  vars = Prepend[TL["syms"], x];
  tc = Together[c]; A = Expand[Numerator[tc]]; U = Expand[Denominator[tc]];
  If[!PolynomialQ[A, th] || !PolynomialQ[U, th], Return[{"Unknown", "c is not a rational function of theta"}]];
  g = Exponent[G, th]; du = Exponent[U, th]; dc = Exponent[A, th] - du;
  N = Floor[Max[dc + 1, 1 - g/m]];
  If[N < 0, Return[If[tc === 0, {"Solution", 0}, {"NoSolution"}]]];
  lcG = Coefficient[G, th, g]; lcU = Coefficient[U, th, du];
  beta = Together[brD[TL, lcG]/(m lcG)];
  u = brMthPowerCofactor[lcG, m];
  h = If[u === None, None, Together[1/u]];
  Dth = brEta[T];
  bs = Table[Symbol["Bronstein`Private`brB" <> ToString[j]], {j, 0, N}];
  dbs = Table[Symbol["Bronstein`Private`brDB" <> ToString[j]], {j, 0, N}];
  P = Expand[Together[m G U Sum[dbs[[j + 1]] th^j + j bs[[j + 1]] Dth th^(j - 1), {j, 0, N}] +
          U brD[T, G] Sum[bs[[j + 1]] th^j, {j, 0, N}] - m G A] /. th^(-1) -> 0];
  cl = CoefficientList[P, th];
  cl = PadRight[cl, N + g + du + 1];
  Do[
    e = j + g + du;
    eq = Together[cl[[e + 1]] /. known];
    alpha = Coefficient[Expand[Numerator[eq]], dbs[[j + 1]]]/Denominator[eq];
    If[Together[alpha] === 0, throwInternal["brPolyRDELogLeadingZero", <|"j" -> j|>]];
    rho = Together[-(eq /. {dbs[[j + 1]] -> 0, bs[[j + 1]] -> 0})/alpha];
    (* the b_j coefficient divided by alpha must be beta *)
    If[Together[Coefficient[Expand[Numerator[eq]], bs[[j + 1]]]/Denominator[eq]/alpha - beta] =!= 0,
      throwInternal["brPolyRDELogBeta", <|"j" -> j|>]];
    If[h === None,
      If[kappas =!= {} && !FreeQ[rho, Alternatives @@ kappas], throwInternal["brPolyRDELogStrayKappa"]];
      sol = brKSolveRDE[TL, beta, rho];
      If[First[sol] =!= "Solution", Return[If[First[sol] === "NoSolution", {"NoSolution"}, {"Unknown", "coefficient RDE: " <> sol[[2]]}]]];
      val = sol[[2]],
      (* b_j = h w, w' = rho/h = rho0 + Sum kappa_i sig_i *)
      rho0 = Together[(rho /. Thread[kappas -> 0])/h];
      sig = Together[Coefficient[Expand[Numerator[Together[rho]]], #]/Denominator[Together[rho]]/h] & /@ kappas;
      pr = brParamIntegrate[TL, rho0, sig];
      If[First[pr] =!= "Solution", Return[If[First[pr] === "NoSolution", {"NoSolution"}, {"Unknown", "coefficient integration: " <> pr[[2]]}]]];
      (* the integration's own parameter symbols are our kappas, in order *)
      cons = Join[cons, pr[[4]] /. Thread[pr[[5]] -> kappas]];
      nk++;
      AppendTo[kappas, Symbol["Bronstein`Private`brK" <> ToString[nk]]];
      w = pr[[2]] + Total[MapThread[Times, {Most[kappas], pr[[3]]}]] + Last[kappas];
      val = Together[h w]];
    known = Join[known, {bs[[j + 1]] -> val, dbs[[j + 1]] -> Together[brD[TL, val]]}],
    {j, N, 0, -1}];
  (* the remaining equations (theta-degrees below g + du) fix the constants *)
  rest = Together[(cl[[# + 1]] /. known)] & /@ Range[0, g + du - 1];
  eqs = Join[cons, Flatten[brIdenticallyZero[#, vars] & /@ rest]];
  eqs = DeleteCases[eqs, True];
  (* a "no" only from a well-formed linear system over the constants *)
  If[!brLinearConstantSystemQ[eqs, kappas], Return[{"Unknown", "constant system not well formed"}]];
  If[MemberQ[eqs, False], Return[{"NoSolution"}]];
  B = Sum[(bs[[j + 1]] /. known) th^j, {j, 0, N}];
  If[eqs =!= {},
    If[kappas === {}, Return[{"NoSolution"}]];
    sols = Solve[eqs, kappas];
    If[sols === {}, Return[{"NoSolution"}]];
    B = B /. First[sols]];
  B = Together[B /. Thread[kappas -> 0]];
  If[!brAlgebraicResidualZeroQ[T, s, c, B], throwInternal["brPolyRDELogResidual", <|"B" -> B|>]];
  {"Solution", B}];

brIdenticallyZero::usage = "brIdenticallyZero[e, vars] is the list of equations (over the constants: symbols not in vars) that make the rational function e of vars vanish identically: every coefficient of its numerator, as a polynomial in vars, equal to 0.";
brIdenticallyZero[e_, vars_List] := Module[{num = Expand[Numerator[Together[e]]]},
  If[num === 0, {}, (# == 0) & /@ (Last /@ CoefficientRules[num, vars])]];

brLinearConstantSystemQ::usage = "brLinearConstantSystemQ[eqs, kappas] is True when every element of eqs is True, False or lhs == rhs with lhs - rhs a polynomial of degree <= 1 in the symbols kappas whose coefficients are explicit numbers (algebraic constants allowed): the only systems from which an inconsistency is read as a proof.";
brLinearConstantSystemQ[eqs_List, kappas_List] := AllTrue[eqs,
  Function[e, e === True || e === False ||
    (Head[e] === Equal && Length[e] === 2 &&
      With[{d = Expand[e[[1]] - e[[2]]]},
        PolynomialQ[d, kappas] && (kappas === {} || Max[Exponent[d, #] & /@ kappas] <= 1) &&
        Complement[Variables[d], kappas] === {} && FreeQ[d, Integrate | Risch`Derivation]])]];

