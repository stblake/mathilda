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
      throwDecline["transcendental part not decided"]]];
  F];

brTermExpression::usage = "brTermExpression[F, e] is the expression e over the field F with the residue-field generator replaced by its algebraic number (Trager's fieldValue for a number field; the theta of a fieldFn extension).";
brTermExpression[F_, e_] := If[F === None || !fieldExtQ[F], e, rc[F, e] /. F["t"] -> F["theta"]];

brIntegrateComponent::usage = "brIntegrateComponent[T, cv, c] integrates the component c Y on the curve cv (Y^m = G). It returns <|\"Algebraic\" -> b, \"Logs\" -> {{coef, arg}, ...}|> with Int c Y = b Y + Sum coef Log[arg], arg expressions in the tower symbols and the symbol Ysym of cv (\"Y\" -> Ysym); it throws the failures of its stages.";
brIntegrateComponent[T_Association, cv_Association, c_] := Module[{hr, sr, lp, ap, b, Ys = Unique["Bronstein`Private`brYc"], logs = {}, F, B},
  If[!brReducedCurveQ[T, cv["m"], cv["G"]], throwInternal["brCurveNotReduced", <|"m" -> cv["m"], "G" -> cv["G"]|>]];
  hr = brHermiteReduce[T, cv, c];
  sr = brSpecialReduce[T, cv, hr["Remainder"]];
  lp = brComponentLogPart[T, cv, sr["Remainder"]];
  ap = brAlgebraicPart[T, cv, lp["remainder"]];
  Switch[First[ap],
    "NoSolution", throwNonElementary["no algebraic part: the Risch DE D b + s b = c has no solution ([B] Prop 2.15/2.18)"],
    "Unknown", throwDecline["algebraic part undecided: " <> ap[[2]]]];
  b = Together[hr["Algebraic"] + sr["Algebraic"] + ap[[2]]];
  F = lp["F"]; B = lp["B"];
  Do[
    logs = Join[logs, {
      {brTermExpression[F, t[[1]]], brTermExpression[F, afToY[B, t[[2]], Ys]]},
      {brTermExpression[F, -t[[1]]/cv["m"]], brTermExpression[F, t[[3]]]}}],
    {t, lp["terms"]}];
  <|"Algebraic" -> b, "Logs" -> logs, "Y" -> Ys|>];

brAssembleComponent::usage = "brAssembleComponent[T, res, Yexpr, x] writes the component integral res (brIntegrateComponent) in the surface: the curve generator replaced by Yexpr (an expression in the tower symbols), the tower symbols by their kernels; logarithms of arguments free of x (constants of integration) are dropped.";
brAssembleComponent[T_Association, res_Association, Yexpr_, x_Symbol] := Module[{Y = res["Y"], alg, logs},
  alg = brToSurface[T, Together[res["Algebraic"] Y /. Y -> Yexpr]];
  logs = Total[If[FreeQ[brToSurface[T, #[[2]] /. Y -> Yexpr], x], 0,
      #[[1]] Log[brToSurface[T, Together[#[[2]] /. Y -> Yexpr]]]] & /@ res["Logs"]];
  alg + logs];

brShapeA::usage = "brShapeA[T, f, x] integrates the Shape A integrand f (tower record T): coordinate 0 by the transcendental Risch algorithm, every component by brIntegrateComponent; it returns the (unverified) antiderivative.";
brShapeA[T_Association, f_, x_Symbol] := Module[{cf, n, g, yh, total, cc, cv, res, Yexpr, rho},
  cf = brCanonicalForm[T];
  n = cf["n"]; g = cf["g"]; yh = cf["y"];
  total = brIntegrateTranscendental[brToSurface[T, cf["coeffs"][[1]]], x];
  Do[
    If[cf["coeffs"][[i + 1]] =!= 0,
      cc = brComponentCurve[T, n, g, i];
      If[cc["m"] === 1,
        (* the component is rational in theta times the locally constant rho *)
        rho = brToSurface[T, (cf["yexpr"]^i/cc["s"])/cc["G"]^(1/cc["m"])];
        total = total + rho brIntegrateTranscendental[brToSurface[T, Together[cf["coeffs"][[i + 1]] cc["s"] cc["G"]]], x],
        cv = brCurve[T, cc["m"], cc["G"]];
        res = brIntegrateComponent[T, cv, Together[cf["coeffs"][[i + 1]] cc["s"]]];
        If[cc["exact"],
          rho = 1; Yexpr = Together[cf["yexpr"]^i/cc["s"]],
          Yexpr = cc["G"]^(1/cc["m"]);
          rho = brToSurface[T, (cf["yexpr"]^i/cc["s"])/Yexpr]];
        total = total + rho brAssembleComponent[T, res, Yexpr, x]]],
    {i, 1, n - 1}];
  total];

brIntegrateBody::usage = "brIntegrateBody[f, x] is the body of BronsteinIntegrate: it throws instead of returning failure lists.";
brIntegrateBody[f0_, x_Symbol] := Module[{f = f0 /. Exp[a_] :> E^a, T, F, r, fl, u = Unique["Bronstein`Private`brUf"]},
  T = brTower[f, x];
  F = Switch[T["shape"],
    "Rational" | "Transcendental", brIntegrateTranscendental[f, x],
    "Algebraic",
      r = Trager`TragerIntegrate[f, x, "Trager"];
      Which[
        MatchQ[r, {"not elementary", _}], throwNonElementary["Trager: " <> r[[2]]],
        MatchQ[r, {"failed", _}], throwDecline["Trager: " <> r[[2]]],
        True, r],
    "B",
      fl = brFlattenLinear[f, x, T["radical"]["b"] /. T["back"], T["radical"]["n"], u];
      If[fl === None, throwDecline["radical over Q(x) below the transcendental tower (Shape B, a later phase)"]];
      r = brIntegrateTranscendental[fl["g"], u];
      r /. fl["back"],
    "A", brShapeA[T, f, x],
    _, throwInternal["brUnknownShape"]];
  If[!antiderivativeQ[f, F, x], throwInternal["AntiderivativeCheckFailed", <|"f" -> f, "F" -> F|>]];
  presentable[f, F, x]];

BronsteinIntegrate[f_, x_Symbol] := Module[{r},
  r = brCatch[brIntegrateBody[f, x]];
  If[brFailQ[r], failToList[r], r]];
