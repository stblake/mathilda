(* ========================================================================== *)
(* GeneralisedGoursat -- symmetry reduction of Integral F(t)/R(t)^p dt,         *)
(* p in {1/2, 1/3, 2/3, 1/4, 3/4}: Integrate's Method -> "GeneralisedGoursat".   *)
(*                                                                            *)
(* A faithful port of the sympy reference implementation goursat.py that       *)
(* accompanies "A Generalisation of Goursat's Algorithm for Integration in     *)
(* Finite Terms" (cube-root paper, 2026-11-26).  Function names follow         *)
(* goursat.py, camelCased, one file per section:                               *)
(*                                                                            *)
(*   goursat/Common.m      K, polyrem, descend, coords, normalize, tidy        *)
(*   goursat/RatInt.m      rat_int, _rat_int_alg                               *)
(*   goursat/Hermite.m     hermite, classify_reduced                           *)
(*   goursat/Sqrt.m        p = 1/2: involutions, F0 trace, conics, order 3,     *)
(*                         second level                                        *)
(*   goursat/CubicField.m  p = 1/2 over an irreducible resolvent cubic          *)
(*   goursat/Cbrt.m        p = 1/3, 2/3: covariants, eigen-decomposition,       *)
(*                         second level                                        *)
(*   goursat/Qrt.m         p = 1/4, 3/4: binomial quartic radicands             *)
(*   goursat/Verify.m      exact verification by formal differentiation        *)
(*   (this file)           driver goursat_integrate, front end                  *)
(*                         integrate_algebraic, verification, Integrate entry   *)
(*                                                                            *)
(* Design (from the reference):                                                *)
(*   * the roots of R are never computed -- only rational covariants of R; for  *)
(*     p = 1/2 the roots of the resolvent cubic, for p = 1/3, 2/3 the two roots *)
(*     of the Hessian;                                                         *)
(*   * eigen-projections are replaced by descent: linear algebra over Q(x), x a  *)
(*     rational generator of the invariant field, with polynomial remainders;  *)
(*   * the answer is written in t and the single radical y = R^(1/n), so it is  *)
(*     verified by formal differentiation (y' = R'/(n y^(n-1))) with no branch   *)
(*     ambiguity;                                                              *)
(*   * an honest three-valued verdict: the residual is Hermite-reduced and       *)
(*     classified -- zero: elementary; non-zero of the second kind: certified    *)
(*     non-elementary; third kind: undecided.  Failure of the symmetry          *)
(*     criterion alone is never reported as non-elementarity.                   *)
(* ========================================================================== *)

BeginPackage["GeneralisedGoursat`", {"Mathilda`"}];

GoursatIntegrate::usage = "GoursatIntegrate[f, x] is Integrate's Method -> \"GeneralisedGoursat\" stage. f must lie in Q(x, R^(1/n)) for a single squarefree polynomial radicand R with n in {2, 3, 4} (genus 1: deg R in {3, 4} for n = 2, {2, 3} for n = 3, a binomial quartic a(x - e)^4 + b for n = 4). It returns an antiderivative (verified by formal differentiation before it is returned), {\"not elementary\", partial, residual} when the remaining integral of residual is certified non-elementary (f = D[partial, x] + residual), or $Failed when it declines or cannot decide.";

GoursatReduce::usage = "GoursatReduce[F, R, t, p] reduces Integral F(t)/R(t)^p dt, p in {1/2, 1/3, 2/3, 1/4, 3/4}, F rational over Q(algebraic numbers), R squarefree. It returns an Association with keys \"Status\" (\"elementary\", \"non-elementary\" or \"undecided\"), \"Antiderivative\" (in t and the symbol \"y\" -> y, y^n = R), \"Residual\" (the reduced F~ with Integral F dt/y^k = Antiderivative + Integral F~ dt/y^k), \"Kind\", \"n\", \"k\", \"R\", \"y\" and \"Log\".";

GoursatIntegrateAlgebraic::usage = "GoursatIntegrateAlgebraic[f, t] integrates f in Q(t, R^(1/n)), n in {2, 3, 4}, by splitting it into eigencomponents G0 + Sum_j Gj y^j under y -> zeta y and reducing each. Returns <|\"Status\", \"Antiderivative\", \"Residual\", \"Parts\"|> with Antiderivative and Residual in t and R^(1/n) (f = D[Antiderivative, t] + Residual).";

GoursatVerify::usage = "GoursatVerify[res, F, t] checks a GoursatReduce result exactly: D[Antiderivative] + Residual/y^k - F/y^k reduces to 0 in Q(t)[y]/(y^n - R), differentiating with y' = R'/(n y^(n-1)). Returns True or False.";

GoursatVerifyExpr::usage = "GoursatVerifyExpr[res, f, t] checks a GoursatIntegrateAlgebraic result: D[Antiderivative, t] + Residual - f is zero (formally in y, else by PossibleZeroQ). Returns True or False.";

GoursatHermite::usage = "GoursatHermite[F, R, t, n, k] is Hermite reduction on y^n = R: {A, P, N, Ds} with F = A' R + (1 - k/n) A R' + P + N/Ds, deg P <= deg R - 2, Ds squarefree and coprime to R, deg N < deg Ds.";
GoursatClassify::usage = "GoursatClassify[P, N, R, t, n, k] is the kind (\"zero\", \"second-kind\", \"third-kind\") of the reduced differential (P + N/Ds) dt/y^k.";
GoursatF0::usage = "GoursatF0[F, R, t] is the Klein-four-group average (F + Sum F o S_theta)/4 for y^2 = R, computed as a trace over the resolvent cubic without its roots.";
GoursatInvolutions::usage = "GoursatInvolutions[R, t] lists the involutions {L, N, M}, S(t) = (N t - M)/(L t - N), pairing the branch points of y^2 = R (resolvent cubic with rational or quadratic roots).";
GoursatCubicCovariants::usage = "GoursatCubicCovariants[R, t] is {H, G}, the Hessian and Jacobian covariants of R as a binary cubic.";
GoursatQuarticCovariants::usage = "GoursatQuarticCovariants[R, t] is {H, T}, the Hessian and sextic covariant T = R' H - R H' of R as a binary quartic.";
GoursatDescend::usage = "GoursatDescend[g, X, t, x] is G with G(X(t)) = g(t), for g known to be a rational function of X.";
GoursatCoords::usage = "GoursatCoords[F, basis, X, t, x] are the coordinates of F in a basis of Q(t) over Q(X(t)).";
GoursatNormalize::usage = "GoursatNormalize[E, R, t, y, n] writes a rational expression in (t, y), y^n = R, as Sum_j c_j(t) y^j.";
$GoursatDebug::usage = "$GoursatDebug (default False) set to True makes every GeneralisedGoursat stage print its name and wall time (to stderr), indented by recursion depth.";

GoursatRationalIntegrate::usage = "GoursatRationalIntegrate[rho, m] is the verified antiderivative of a rational function over Q(algebraic numbers) used by every reduction.";

$GoursatDebug = False;

Begin["`Private`"];

LoadModule["goursat/Common.m"];
LoadModule["goursat/RatInt.m"];
LoadModule["goursat/Hermite.m"];
LoadModule["goursat/Sqrt.m"];
LoadModule["goursat/CubicField.m"];
LoadModule["goursat/Cbrt.m"];
LoadModule["goursat/Qrt.m"];
LoadModule["goursat/Verify.m"];

maxDepth[] := 2;            (* recursion bound of the p = 1/2 second level *)
secondLevelQ[] := True;     (* False stops after the first-level reductions *)

nk[p_] := Switch[p, 1/2, {2, 1}, 1/3, {3, 1}, 2/3, {3, 2}, 1/4, {4, 1}, 3/4, {4, 3},
                 _, ggDecline["unsupported exponent"]];

mkResult[status_, anti_, res_, kind_, n_, k_, R_, y_, log_] :=
  <|"Status" -> status, "Antiderivative" -> anti, "Residual" -> res, "Kind" -> kind,
    "n" -> n, "k" -> k, "R" -> R, "y" -> y, "Log" -> log|>;

(* ---- driver: goursat_integrate ------------------------------------------- *)
goursatIntegrate[F0_, R0_, t_, p_, y_, depth_] := Module[
  {n, k, F, R, d, A, P, Nn, Ds, anti, Fred, log, j0, pre, part0, Fin, lg0, part, Fres,
   lg, o3, A2, P2, N2, Ds2, kind, res, sub, back, status},
  {n, k} = nk[p];
  F = gK[F0]; R = Expand[R0];
  d = Exponent[R, t];
  If[!PolynomialQ[R, t] || Exponent[PolynomialGCD[R, D[R, t]], t] != 0,
    ggDecline["R must be a squarefree polynomial"]];
  If[!((n == 2 && MemberQ[{3, 4}, d]) || (n == 3 && MemberQ[{2, 3}, d]) ||
       (n == 4 && quarticBinomial[R, t] =!= None)),
    ggDecline["unsupported radicand"]];

  {A, P, Nn, Ds} = stage["hermite", hermite[F, R, t, n, k]];
  anti = A y^(n - k);
  Fred = gK[P + Nn/Ds];
  log = {<|"hermiteAlgebraicPart" -> A, "hermiteReduced" -> Fred|>};
  If[Fred === 0, Return[mkResult["elementary", anti, 0, "zero", n, k, R, y, log], Module]];

  Which[
    n == 2,
      (* j = 0: Goursat's order-3 substitutions.  The eigencomponent is peeled off *)
      (* before the Klein-group step (simplest output) and again from the          *)
      (* V4-invariant residual (where one order-3 element stands for all).         *)
      j0 = isZero[quarticI[R, t]];
      pre = If[j0, order3Reduce[Fred, R, t, y], None];
      {part0, Fin, lg0} = If[pre === None, {0, Fred, {}}, pre];
      {part, Fres, lg} = If[Fin =!= 0, stage["sqrtReduce", sqrtReduce[Fin, R, t, y]], {0, Fin, {}}];
      part = part0 + part; lg = Join[lg0, lg];
      Which[
        j0 && pre === None,
          AppendTo[lg, <|"note" -> "j = 0 but no rational branch point: order-3 reduction not attempted"|>],
        j0 && Fres =!= 0,
          o3 = order3Reduce[Fres, R, t, y];
          part = part + o3[[1]]; Fres = o3[[2]]; lg = Join[lg, o3[[3]]]],
    n == 3,
      {part, Fres, lg} = stage["cbrtReduce", cbrtReduce[Fred, R, t, y, k]],
    True,
      {part, Fres, lg} = stage["qrtReduce", qrtReduce[Fred, R, t, y, k]]];
  anti = anti + part;
  log = Join[log, lg];

  {A2, P2, N2, Ds2} = stage["hermite (residual)", hermite[Fres, R, t, n, k]];
  kind = classifyReduced[P2, N2, R, t, n, k];
  res = gK[P2 + N2/Ds2];

  (* second level: the residual is again a Goursat (square-root) integral *)
  If[secondLevelQ[] && kind === "third-kind" &&
     ((n == 2 && depth < maxDepth[]) || n == 3) && !specialisedUndecidedQ[Fres, R, t, y, n, k, depth],
    {sub, back} = If[n == 2, sqrtSecondLevel[Fres, R, t, y, depth],
                             cbrtSecondLevel[Fres, R, t, y, k]];
    AppendTo[log, <|"secondLevel" -> sub["Status"], "PhiResidual" -> sub["Residual"]|>];
    If[sub["Status"] =!= "undecided",
      anti = tidy[anti, R, t, y, n] + tidyArgs[sub["Antiderivative"] /. back, R, t, y, n];
      If[sub["Status"] === "elementary",
        Return[mkResult["elementary", anti, 0, "zero", n, k, R, y, log], Module]];
      Return[mkResult["non-elementary", anti, gK[Fres], sub["Kind"], n, k, R, y, log], Module]]];

  anti = stage["tidy", tidy[anti + A2 y^(n - k), R, t, y, n]];
  status = Switch[kind, "zero", "elementary", "second-kind", "non-elementary", _, "undecided"];
  mkResult[status, anti, res, kind, n, k, R, y, log]];

(* ---- parameters at the second level ----------------------------------------- *)
(* Over Q(params) the second level is a chain of multiplication-by-2 descents   *)
(* whose integrands grow at every step: (2 - (k + 1) t)/((1 - (k + 1) t)        *)
(* (t (1 - t)(1 - k t))^(2/3)) reached a 7000-leaf residual at depth 2 and ran   *)
(* the stage's 30 s budget out, to end undecided.  The verdict of the generic   *)
(* problem is that of almost every specialisation: if the residue divisor is    *)
(* torsion over Q(params), it stays torsion (of an order dividing) at every     *)
(* parameter value where the curve and the poles stay distinct.  So the second  *)
(* level is first run at two parameter values; when both end undecided the     *)
(* parametric one is not attempted and the reduction ends undecided -- a        *)
(* decline, never an answer or a certificate.  A specialisation that degenerates *)
(* (R no longer squarefree, a pole on the branch locus) or fails is no evidence, *)
(* and the parametric second level runs as before.                             *)
specialisationValues[] := {{7/11, 13/5, -5/3, 11/7, 3/13}, {-9/7, 17/19, 5/2, -2/11, 23/9}};

specialisedUndecidedQ[Fres_, R_, t_, y_, n_, k_, depth_] := Module[{ps, verdicts},
  ps = Complement[Cases[Variables[{Fres, R}], _Symbol], {t, y}];
  If[ps === {} || Length[ps] > 5, Return[False, Module]];
  verdicts = Table[
    With[{rules = Thread[ps -> Take[vals, Length[ps]]]},
      Module[{Fs = gK[Fres /. rules], Rs = Expand[R /. rules], r},
        r = If[!PolynomialQ[Rs, t] || Exponent[Rs, t] != Exponent[R, t] ||
               Exponent[PolynomialGCD[Rs, D[Rs, t]], t] != 0 || Fs === 0,
          "degenerate",
          Quiet[TimeConstrained[Catch[
            (If[n == 2, sqrtSecondLevel[Fs, Rs, t, y, depth], cbrtSecondLevel[Fs, Rs, t, y, k]][[1]])["Status"],
            "GoursatFail"], 10, "timeout"]]];
        r]],
    {vals, specialisationValues[]}];
  AllTrue[verdicts, # === "undecided" &]];

(* ---- front end: integrate_algebraic ---------------------------------------- *)
radicalsOf[e_, t_] := Union[Cases[e, Power[b_, q_Rational] /; !FreeQ[b, t], {0, Infinity}]];

integrateAlgebraic[expr_, t_] := Module[
  {rads, n, bases, R, y, y2, E1, cs, anti, resid = 0, parts = {}, stat = {}, r, rad, status},
  rads = radicalsOf[expr, t];
  If[rads === {},
    Return[<|"Status" -> "elementary", "Antiderivative" -> ratInt[gK[expr], t],
             "Residual" -> 0, "Parts" -> {}|>, Module]];
  n = LCM @@ (Denominator[#[[2]]] & /@ rads);
  If[!MemberQ[{2, 3, 4}, n], ggDecline["only square, cube and fourth roots"]];
  bases = Union[Expand[#[[1]]] & /@ rads];
  If[Length[bases] != 1 || !PolynomialQ[First[bases], t],
    ggDecline["need a single polynomial radicand"]];
  R = First[bases];
  E1 = expr /. Power[b_, q_Rational] /; !FreeQ[b, t] :> y^(q n);
  If[radicalsOf[E1, t] =!= {},
    ggDecline["radical left after substitution"]];
  cs = CoefficientList[stage["normalize", normalizeAlgebraic[E1, R, t, y, n]], y];
  cs = Join[cs, ConstantArray[0, n - Length[cs]]];
  anti = ratInt[gK[cs[[1]]], t];
  Do[If[cs[[j + 1]] =!= 0,
       If[n == 4 && j == 2,
         (* y^2 = sqrt(R): a square-root integral *)
         r = goursatIntegrate[gK[cs[[j + 1]] R], R, t, 1/2, y2, 0];
         anti += r["Antiderivative"] /. y2 -> y^2;
         resid += r["Residual"]/y^2,
         r = goursatIntegrate[gK[cs[[j + 1]] R], R, t, (n - j)/n, y, 0];
         anti += r["Antiderivative"];
         resid += r["Residual"]/y^(n - j)];
       AppendTo[parts, r]; AppendTo[stat, r["Status"]]],
     {j, 1, n - 1}];
  status = Which[MemberQ[stat, "non-elementary"], "non-elementary",
                 MemberQ[stat, "undecided"], "undecided",
                 True, "elementary"];
  rad = R^(1/n);
  <|"Status" -> status, "Antiderivative" -> (anti /. y -> rad),
    "Residual" -> (resid /. y -> rad), "Parts" -> parts, "R" -> R, "n" -> n, "y" -> y, "t" -> t,
    "AntiderivativeY" -> anti, "ResidualY" -> resid|>];

(* ---- public API ------------------------------------------------------------ *)
catchG[e_] := Catch[e, "GoursatFail"];
SetAttributes[catchG, HoldAll];

GoursatReduce[F_, R_, t_Symbol, p_] := Module[{y, r},
  r = catchG[goursatIntegrate[F, R, t, p, Symbol["GeneralisedGoursat`y"], 0]];
  If[MatchQ[r, _ggFail], $Failed, r]];

GoursatIntegrateAlgebraic[f_, t_Symbol] := Module[{r},
  r = catchG[integrateAlgebraic[f, t]];
  If[MatchQ[r, _ggFail], $Failed, r]];

GoursatVerify[res_Association, F_, t_Symbol] := verifyResult[res, F, t];
GoursatVerifyExpr[res_Association, f_, t_Symbol] := verifyExprResult[res, f, t];

GoursatHermite[F_, R_, t_Symbol, n_Integer, k_Integer] := catchG[hermite[F, R, t, n, k]];
GoursatClassify[P_, Nn_, R_, t_Symbol, n_Integer, k_Integer] := classifyReduced[P, Nn, R, t, n, k];
GoursatF0[F_, R_, t_Symbol] := catchG[goursatF0[F, R, t]];
GoursatInvolutions[R_, t_Symbol] := catchG[involutions[R, t]];
GoursatCubicCovariants[R_, t_Symbol] := cubicCovariants[R, t];
GoursatQuarticCovariants[R_, t_Symbol] := quarticCovariants[R, t];
GoursatDescend[g_, X_, t_Symbol, x_Symbol] := catchG[descend[g, X, t, x]];
GoursatCoords[F_, basis_List, X_, t_Symbol, x_Symbol] := catchG[coords[F, basis, X, t, x]];
GoursatNormalize[E_, R_, t_Symbol, y_Symbol, n_Integer] := catchG[normalizeAlgebraic[E, R, t, y, n]];
GoursatRationalIntegrate[rho_, m_Symbol] := catchG[ratInt[rho, m]];

(* The Integrate stage.  A verified antiderivative, a certified partial answer, *)
(* or $Failed (decline: out of shape, undecided, unverified).                   *)
GoursatIntegrate[f_, x_Symbol] := Block[{$ggDefer = True, $ggDeferred = {}}, Module[{r, vd = False},
  r = stage["integrateAlgebraic", catchG[integrateAlgebraic[f, x]]];
  If[!AssociationQ[r], Return[$Failed, Module]];
  (* The verdict is in.  An undecided reduction ends here, without integrating *)
  (* its pieces (RatInt.m, deferred integration).  Otherwise the answer is      *)
  (* verified before its pieces are integrated, by the chain rule (Verify.m,     *)
  (* verifyDeferred: an algebraic identity, where the closed form would put      *)
  (* logarithms and RootSums of the pieces to PossibleZeroQ), and the deferred   *)
  (* integrals are computed; the closed form is verified only if that failed.   *)
  If[MemberQ[{"elementary", "non-elementary"}, r["Status"]],
    vd = stage["verify (deferred)", verifyDeferred[r, f, x]];
    r = stage["integrate pieces", catchG[resolveAnswer[r]]];
    If[!AssociationQ[r], Return[$Failed, Module]]];
  Switch[r["Status"],
    "elementary",
      If[vd || stage["verify", verifyExprResult[r, f, x]], r["Antiderivative"], $Failed],
    "non-elementary",
      If[vd || stage["verify", verifyExprResult[r, f, x]], {"not elementary", r["Antiderivative"], r["Residual"]}, $Failed],
    _, $Failed]]];

(* Is there anything for tidy to do: algebraic terms to combine (the rational   *)
(* pieces of a resolved answer can cancel to 0), or a logarithm argument with y *)
(* in a denominator (a second-level back-substitution)?  A first-level answer   *)
(* has neither, and its tidy only re-evaluated it (0.1 s over Q(k)(I Sqrt[3])). *)
needsTidyQ[ay_, y_] := Module[{terms = If[Head[ay] === Plus, List @@ ay, {ay}]},
  Count[terms, tm_ /; !transcendentalQ[tm]] > 1 ||
  !FreeQ[Cases[ay, (Log | ArcTan | ArcTanh | ArcCot | ArcCoth)[a_] :> a, {0, Infinity}],
         Power[b_, e_?Negative] /; !FreeQ[b, y]]];

(* An integrateAlgebraic result with its deferred integrals computed, the     *)
(* logarithm arguments brought to normal form (the tidy step the reduction     *)
(* applied to the pieces before they were deferred).                          *)
resolveAnswer[r_] := Module[{ay, n = r["n"], R = r["R"], y = r["y"]},
  If[!KeyExistsQ[r, "AntiderivativeY"],
    Return[Append[r, "Antiderivative" -> resolveDeferred[r["Antiderivative"]]], Module]];
  ay = resolveDeferred[r["AntiderivativeY"]];
  If[ay =!= r["AntiderivativeY"] && needsTidyQ[ay, y],
    ay = stage["tidy (resolved)", tidy[ay, R, r["t"], y, n]]];
  Join[r, <|"AntiderivativeY" -> ay, "Antiderivative" -> (ay /. y -> R^(1/n))|>]];

End[];
EndPackage[];
