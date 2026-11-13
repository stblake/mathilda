(* ::Package:: *)

(* ========================================================================== *)
(* Pipeline.m -- the per-component driver and the public entry point          *)
(*               BronsteinIntegrate (BRONSTEIN_PLAN.md "Per-component         *)
(*               pipeline").                                                  *)
(*                                                                            *)
(* Shape A: f = Sum_i a_i yh^i over k(theta) (Curve.m).  Coordinate 0 is an    *)
(* element of k(theta), integrated by the transcendental Risch algorithm.     *)
(* Each component a_i yh^i is c Y on its curve Y^m = G ([T] App A Sec 1: f is *)
(* elementary iff every component is), and                                    *)
(*     Int c Y = b Y + L1,                                                    *)
(* b = Hermite part + special-place part + algebraic part, L1 the weight-1     *)
(* log part.  Every certificate on the way is a theorem of [B]/[T]: the        *)
(* resultant criterion, an undefined special residue, a residue divisor of    *)
(* infinite order, a missing reduction or algebraic part (proven by the       *)
(* certified Risch DE), a non-elementary coordinate 0.  Every antiderivative   *)
(* is checked by PossibleZeroQ[D[F, x] - f] before it is returned; a failure  *)
(* there is an internal error (a decline), never an answer.                   *)
(* ========================================================================== *)

BronsteinIntegrate::usage = "BronsteinIntegrate[f, x] integrates f, an elementary function of x built from rational operations, logarithms, exponentials and ONE radical (n-th root) whose radicand involves the top transcendental monomial (Bronstein 1990, the radical on top of a Log/Exp tower; a purely algebraic integrand goes to Trager's algorithm and a purely transcendental one to the transcendental Risch algorithm). It returns an antiderivative (verified), {\"not elementary\", reason} (a certificate), or {\"failed\", reason} (out of scope or undecided).";

brIntegrateTranscendental::usage = "brIntegrateTranscendental[f, x] is the integral of the transcendental (radical-free) f by Integrate`RischTranscendental, verified; it throws NonElementary when Risch`ElementaryIntegralQ proves there is none and Decline when undecided.";
brIntegrateTranscendental[f_, x_Symbol] := Module[{F, q},
  If[Together[f] === 0, Return[0]];
  F = Quiet[Integrate`RischTranscendental[f, x]];
  If[Head[F] === Integrate`RischTranscendental || !FreeQ[F, brNonElementaryHeads[]],
    q = Quiet[Risch`ElementaryIntegralQ[f, x]];
    If[q === False, throwNonElementary["transcendental part has no elementary integral (Risch)"],
      throwUnimplemented["Pipeline:transcendental", "transcendental part not decided (Integrate`RischTranscendental)"]]];
  F];

brTermExpression::usage = "brTermExpression[F, e] is the expression e over the field F with the residue-field generator replaced by its algebraic number (Trager's fieldValue for a number field; the theta of a fieldFn extension).";
brTermExpression[F_, e_] := If[F === None || !fieldExtQ[F], e, rc[F, e] /. F["t"] -> F["theta"]];

brSplitLogs::usage = "brSplitLogs[F, logs, Y, vars] rewrites the log terms {{c, u}, ...} of a component (u over the residue field F, its generator t a symbol) in fewer, smaller logarithms: an argument free of the curve generator Y (a norm term, an element of k(theta)(t)) is split into its irreducible factors over Q[vars, t], c Log[u] -> Sum e_i c Log[p_i], dropping the factors free of vars (x and the tower symbols: they change F only by a constant); then terms with the same argument are merged (coefficients added in F, zero sums dropped). Every step changes Sum c Log[u] by a constant at most, so its derivative -- what brComponentIdentityQ checks -- is unchanged.";
brSplitLogs[F_, logs_List, Y_Symbol, vars_List] := Module[{out = {}, fl, t, merged},
  Do[
    If[!FreeQ[lg[[2]], Y], AppendTo[out, lg]; Continue[]];
    t = Together[lg[[2]]];
    fl = Join[FactorList[Numerator[t]], {#[[1]], -#[[2]]} & /@ FactorList[Denominator[t]]];
    Do[If[brDependsQ[q[[1]], vars], AppendTo[out, {lg[[1]] q[[2]], q[[1]]}]], {q, fl}],
    {lg, logs}];
  merged = {#[[1, 2]], #} & /@ GatherBy[out, #[[2]] &];
  merged = {If[F === None, Together[Total[#[[2, All, 1]]]], kc[F, Together[Total[#[[2, All, 1]]]]]], #[[1]]} & /@ merged;
  Select[merged, Together[#[[1]]] =!= 0 &]];

brCollectLogTerms::usage = "brCollectLogTerms[e, x] merges the logarithms of an antiderivative written as a sum: the terms c Log[a] with c free of x are grouped by their argument a and their coefficients added (a zero sum drops the term). An identity.";
brCollectLogTerms[e_, x_Symbol] := Module[{terms = If[Head[e] === Plus, List @@ e, {e}], logs, rest, g},
  logs = Select[terms, MatchQ[#, (c_. Log[a_]) /; FreeQ[c, x] && !FreeQ[a, x]] &];
  rest = Complement[terms, logs];
  logs = Replace[#, c_. Log[a_] :> {c, a}] & /@ logs;
  g = GatherBy[logs, #[[2]] &];
  Total[rest] + Total[With[{c = Together[Total[#[[All, 1]]]]}, If[c === 0, 0, c Log[#[[1, 2]]]]] & /@ g]];

brIntegrateComponent::usage = "brIntegrateComponent[T, cv, c] integrates the component c Y on the curve cv (Y^m = G). It returns <|\"Algebraic\" -> b, \"Logs\" -> {{coef, arg}, ...}|> with Int c Y = b Y + Sum coef Log[arg], arg expressions in the tower symbols and the symbol Ysym of cv (\"Y\" -> Ysym); it throws the failures of its stages.";
brIntegrateComponent[T_Association, cv_Association, c_] := Module[{hr, sr, lp, ap, b, Ys = Unique["Bronstein`Private`brYc"], logs = {}, F, B, orbits},
  If[!brReducedCurveQ[T, cv["m"], cv["G"]], throwInternal["brCurveNotReduced", <|"m" -> cv["m"], "G" -> cv["G"]|>]];
  hr = brHermiteReduce[T, cv, c];
  sr = brSpecialReduce[T, cv, hr["Remainder"]];
  lp = brComponentLogPart[T, cv, sr["Remainder"]];
  ap = brAlgebraicPart[T, cv, lp["remainder"]];
  Switch[First[ap],
    "NoSolution", throwNonElementary["no algebraic part: the Risch DE D b + s b = c has no solution ([B] Prop 2.15/2.18)"],
    "Unknown", throwUnimplemented["AlgebraicPart:undecided", "algebraic part undecided: " <> ap[[2]]]];
  b = Together[hr["Algebraic"] + sr["Algebraic"] + ap[[2]]];
  F = lp["F"]; B = lp["B"];
  Do[
    logs = Join[logs, {
      {t[[1]], afToY[B, t[[2]], Ys]},
      {-t[[1]]/cv["m"], t[[3]]}}],
    {t, lp["terms"]}];
  logs = brSplitLogs[F, logs, Ys, Prepend[T["syms"], T["x"]]];
  (* per-orbit log terms (brOrbitLogPart): RootSums over the orbit
     polynomials, each over its own field; the exact check takes their
     traced derivative "OrbitD" instead of differentiating them *)
  orbits = brOrbitLogExpression[B, #, Ys] & /@ Lookup[lp, "orbits", {}];
  (* raw: coefficients and logands over F (the residue-field generator t still a
     symbol), as the exact verification needs them; "Logs" the presented ones *)
  <|"Algebraic" -> b, "RawLogs" -> logs, "Field" -> F,
    "OrbitD" -> Lookup[lp, "orbitD", None], "OrbitLogs" -> orbits,
    "Logs" -> Join[({brTermExpression[F, #[[1]]], brTermExpression[F, #[[2]]]} & /@ logs), orbits], "Y" -> Ys|>];

brAssembleComponent::usage = "brAssembleComponent[T, res, Yexpr, x] writes the component integral res (brIntegrateComponent) in the surface: the curve generator replaced by Yexpr (an expression in the tower symbols), the tower symbols by their kernels; logarithms of arguments free of x (constants of integration) are dropped.";
brAssembleComponent[T_Association, res_Association, Yexpr_, x_Symbol] := Module[{Y = res["Y"], alg, logs},
  alg = brToSurface[T, Together[res["Algebraic"] Y /. Y -> Yexpr]];
  logs = Total[If[FreeQ[brToSurface[T, #[[2]] /. Y -> Yexpr], x], 0,
      If[Length[#] === 3, #[[1]] brToSurface[T, #[[2]] /. Y -> Yexpr],
        #[[1]] Log[brToSurface[T, Together[#[[2]] /. Y -> Yexpr]]]]] & /@ res["Logs"]];
  alg + logs];

brComponentIdentityQ::usage = "brComponentIdentityQ[T, cv, c, res] verifies EXACTLY that res integrates the component c Y on its curve Y^m = G (cv): with the tower derivation extended by D Y = s Y, s = DG/(mG), the numerator of D(b Y + Sum c_j Log[u_j]) - c Y reduces to 0 modulo Y^m - G and the residue field's minimal polynomial mu(t). This is an algebraic identity in k(theta)(Y); writing Y as yh^i/s (or G^(1/m)) and the symbols as surface kernels preserves derivatives (every kernel's derivative is its rule), so it is the statement PossibleZeroQ[D[F, x] - f] tests, decided exactly.";
brComponentIdentityQ[T_Association, cv_Association, c_, res_Association] := Module[
  {Y = res["Y"], m = cv["m"], G = cv["G"], rules, Fi, R, num, F = res["Field"]},
  rules = Append[T["rules"], Y -> cv["s"] Y];
  Fi = res["Algebraic"] Y + Total[(#[[1]] Log[#[[2]]]) & /@ res["RawLogs"]];
  R = Risch`Derivation[Fi, rules] - c Y;
  (* per-orbit RootSum terms: their derivative, traced exactly from each
     orbit field (brOrbitLogPart), on the basis w_i = Y^i/d_i *)
  If[ListQ[res["OrbitD"]], R = R + afToY[cv["B"], res["OrbitD"], Y]];
  num = Expand[Numerator[Together[R]]];
  num = Expand[PolynomialRemainder[num, Y^m - G, Y]];
  If[F =!= None && fieldExtQ[F], num = Expand[PolynomialRemainder[num, F["mu"], F["t"]]]];
  Together[num] === 0];

brShapeA::usage = "brShapeA[T, f, x] integrates the Shape A integrand f (tower record T): coordinate 0 by the transcendental Risch algorithm, every component by brIntegrateComponent; it returns the (unverified) antiderivative.";
brShapeA[T_Association, f_, x_Symbol] := Module[{cf, n, g, yh, total, cc, cv, res, Yexpr, rho, a0, exact = True},
  cf = brCanonicalForm[T];
  n = cf["n"]; g = cf["g"]; yh = cf["y"];
  a0 = brToSurface[T, cf["coeffs"][[1]]];
  total = brIntegrateTranscendental[a0, x];
  (* coordinate 0 is small: checked on its own *)
  If[!TrueQ[antiderivativeQ[a0, total, x]], throwInternal["brCoordinate0Check", <|"a0" -> a0|>]];
  Do[
    If[cf["coeffs"][[i + 1]] =!= 0,
      cc = brComponentCurve[T, n, g, i];
      If[cc["m"] === 1,
        (* the component is rational in theta times the locally constant rho *)
        rho = brToSurface[T, (cf["yexpr"]^i/cc["s"])/cc["G"]^(1/cc["m"])];
        exact = False;
        total = total + rho brIntegrateTranscendental[brToSurface[T, Together[cf["coeffs"][[i + 1]] cc["s"] cc["G"]]], x],
        cv = brCurve[T, cc["m"], cc["G"]];
        res = brIntegrateComponent[T, cv, Together[cf["coeffs"][[i + 1]] cc["s"]]];
        If[!brComponentIdentityQ[T, cv, Together[cf["coeffs"][[i + 1]] cc["s"]], res],
          throwInternal["brComponentIdentity", <|"i" -> i|>]];
        (* RootSum terms: the identity used their traced derivative, so the
           assembled answer is checked by PossibleZeroQ as well *)
        If[res["OrbitLogs"] =!= {}, exact = False];
        If[cc["exact"],
          rho = 1; Yexpr = Together[cf["yexpr"]^i/cc["s"]],
          Yexpr = cc["G"]^(1/cc["m"]);
          rho = brToSurface[T, (cf["yexpr"]^i/cc["s"])/Yexpr]];
        total = total + rho brAssembleComponent[T, res, Yexpr, x]]],
    {i, 1, n - 1}];
  <|"F" -> brCollectLogTerms[total, x], "verified" -> exact|>];

brIntegrateBody::usage = "brIntegrateBody[f, x] is the body of BronsteinIntegrate: it throws instead of returning failure lists.";
brIntegrateBody[f0_, x_Symbol] := Module[{f = f0 /. Exp[a_] :> E^a, T, F, r, fl, g0, u = Unique["Bronstein`Private`brUf"], verified = False},
  (* Shape B with a radical inside a kernel, on a genus-0 curve: a rational
     change of variable removes the radical, kernels and all, and the
     transcendental algorithm in u decides.  With radical-free kernels the
     component method (brShapeB) keeps the user's kernels instead. *)
  brCheckHeads[f, x];
  g0 = If[FreeQ[f, (Log[a_] | Power[E, a_]) /; brRadicalQ[a, {x}]], None, brGenus0Flatten[f, x, u]];
  If[g0 =!= None,
    (* a certificate of the parametrised integral is one for f (Q(x)(y) = Q(u));
       a decline falls through to the curve-monomial method below *)
    r = brCatch[brCurveCanon[brIntegrateTranscendental[g0["g"], u] /. g0["back"], x, g0["b"], g0["n"]]];
    If[brFailQ[r] && !tragerRetryableQ[r], Throw[r, "TragerFail"]];
    If[!brFailQ[r],
      If[!antiderivativeQ[f, r, x], throwInternal["AntiderivativeCheckFailed", <|"f" -> f, "F" -> r|>]];
      Return[brPresentable[f, r, x], Module]]];
  T = brTower[f, x];
  F = Switch[T["shape"],
    "Rational" | "Transcendental", brIntegrateTranscendental[f, x],
    "Algebraic",
      r = Trager`TragerIntegrate[f, x, "Trager"];
      Which[
        MatchQ[r, {"not elementary", _}], throwNonElementary["Trager: " <> r[[2]]],
        MatchQ[r, {"implementation incomplete", _}], throwIncomplete["Trager: " <> r[[2]]],
        MatchQ[r, {"failed", _}], throwUnimplemented["Trager", "Trager: " <> r[[2]]],
        True, r],
    "B", With[{sb = brShapeB[T, f, x]}, sb["F"]],
    "BCurve", With[{sb = brShapeBMono[T, f, x]}, sb["F"]],
    "A", With[{sa = brShapeA[T, f, x]}, If[sa["verified"], verified = True]; sa["F"]],
    _, throwInternal["brUnknownShape"]];
  If[!verified && !antiderivativeQ[f, F, x], throwInternal["AntiderivativeCheckFailed", <|"f" -> f, "F" -> F|>]];
  brPresentable[f, F, x]];

brPresentable::usage = "brPresentable[f, F, x] is Trager's real-form presentation of the verified antiderivative F (presentable, which returns the real form only when antiderivativeQ confirms it: the presentation may rewrite logarithms, and Log[2 + 2x] -> Log[1 + x] inside a radical changed the answer), else F itself. It is optional, so it is bounded by brPresentationTime[]: a presentation that overran an enclosing budget (the Automatic cascade's $BronsteinCascadeBudget) aborted the whole call and lost the verified answer -- lccontent-2-3 in the 2026-11-12 shootout, solved in 2.1 s, spent 5 s proving its ArcTan form and missed the 10 s budget under load.";
brPresentable[f_, F_, x_Symbol] := Quiet[TimeConstrained[presentable[f, F, x], brPresentationTime[], F]];

brPresentationTime::usage = "brPresentationTime[] is the time the real-form presentation may take: brOption[\"PresentationTime\"], and at most a third of what is left of an enclosing TimeConstrained (Internal`TimeConstrainedRemaining[]).";
brPresentationTime[] := With[{left = Internal`TimeConstrainedRemaining[]},
  If[NumericQ[left], Min[brOption["PresentationTime"], Max[left/3, 0.01]], brOption["PresentationTime"]]];

BronsteinScope[f_, x_Symbol] := Module[{r},
  If[FreeQ[f, Power[b_, _Rational] /; !FreeQ[b, x]] && FreeQ[f, Root | RootSum], Return["NotApplicable"]];
  r = brCatch[brTower[f, x]];
  Which[
    !brFailQ[r], "InScope",
    MemberQ[{"Incomplete", "Unimplemented"}, First[brFailTag[r]]], failToList[r],
    True, "NotApplicable"]];

BronsteinIntegrate[f_, x_Symbol] := Module[{r},
  r = brCatch[brIntegrateBody[f, x]];
  If[brFailQ[r], failToList[r], r]];
