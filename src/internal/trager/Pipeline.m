(* ::Package:: *)

(* ========================================================================== *)
(* Pipeline.m -- the per-component driver and the public entry points.        *)
(*                                                                            *)
(* TragerIntegrate[f, x, method]                                              *)
(*   1. canonicalForm: f = Sum_i a_i(x) yh^i, yh^n = g in Q[x] (Validate.m).  *)
(*   2. a_0 is a rational function: integrated by Integrate (exact).          *)
(*   3. each a_i yh^i, i >= 1, lives on its own curve Y^m = G (Reduce.m) and   *)
(*      is integrated there by integrateComponent:                            *)
(*        a. move infinity to a finite point if the differential has a pole   *)
(*           there (Shift.m);                                                 *)
(*        b. Hermite reduction (Hermite.m);                                   *)
(*        c. a pole at infinity left in the remainder proves the integral is  *)
(*           not elementary ([T] Ch 4 Sec 3 Theorem 1; [B] p.17 step 5);      *)
(*        d. residues (Residues.m); none while the remainder is nonzero       *)
(*           proves it is not elementary (a nonzero holomorphic differential);*)
(*        e. logarithmic part (TragerLogTerms.m or KauersLogTerms.m);         *)
(*        f. EXACT verification r = Sum c du/u in K(x)(Y); for Trager's        *)
(*           construction a nonzero difference is a holomorphic differential  *)
(*           and proves non-elementarity ([T] p.71), for Kauers' heuristic it  *)
(*           is a decline;                                                    *)
(*        g. back to x (inverse shift) and to the user's radical.             *)
(*   4. a final PossibleZeroQ[D[result] - f] guards the reassembly;  *)
(*      a failure there is an internal error (decline), never an answer.      *)
(* ========================================================================== *)

integrateComponent::usage = "integrateComponent[c, x, m, G, method] integrates c(x) Y dx on the curve Y^m = G(x) (reduced and absolutely irreducible; see reducedCurveQ). method is \"Trager\" or \"Kauers\". It returns <|\"Algebraic\" -> v, \"Logs\" -> {{coef, arg}, ..., {coef, rootsum, \"RootSum\"}, ...}, \"Y\" -> Y|>, the integral being v + Sum coef Log[arg] + Sum coef rootsum as expressions in x and the symbol Y (algebraic numbers explicit). Throws NonElementary / Decline / Internal failures.";
integrateComponent[c_, x_Symbol, m_Integer, G_, method_String] := Module[
  {F = fieldQ[], B0, f, sh = None, z, w, B, her, v, r, cr, num, den, genus, res, ctx, Y, vExpr, logExprs},
  If[!reducedCurveQ[m, G, x], throwInternal["CurveNotReduced", <|"m" -> m, "G" -> G|>]];
  Y = Unique["Y"];
  genus = curveGenus[m, G, x];
  B0 = curveBasis[m, G, x];
  f = PadRight[{0, Together[c]}, m];
  If[regularAtInfinityQ[F, B0, f],
    w = x; B = B0,
    z = Unique["z"];
    sh = shiftCurve[c, m, G, x, z];
    w = z; B = curveBasis[m, sh["G"], z];
    f = PadRight[{0, sh["c"]}, m]];
  her = hermiteReduce[B, f];
  v = her["Algebraic"]; r = her["Remainder"];
  If[!AllTrue[Delete[r, 2], # === 0 &] || !AllTrue[Delete[v, 2], # === 0 &],
    throwInternal["ComponentLeftCoordinateOne"]];
  If[!regularAtInfinityQ[F, B, r],
    throwNonElementary["the Hermite remainder has a pole at infinity"]];
  cr = Together[r[[2]]];
  logExprs = {};
  If[cr =!= 0,
    num = Numerator[cr]; den = Denominator[cr];
    res = componentResiduePolynomial[num, den, m, B["G"], w];
    If[res["R"] === 1,
      throwNonElementary["the Hermite remainder is a nonzero holomorphic differential (no residues)"]];
    ctx = Join[res, <|"B" -> B, "genus" -> genus, "N" -> num, "E" -> den,
                      "Ep" -> Expand[D[den, w]], "Ysym" -> Unique["Yk"]|>];
    logExprs = Switch[method,
      "Trager", tragerLogPart[ctx, r, Y],
      "Kauers", kauersLogPart[ctx, r, Y],
      _, throwInternal["UnknownMethod"]]];
  vExpr = afToY[B, v, Y];
  If[sh =!= None,
    vExpr = unshift[vExpr, sh, x, z, Y];
    logExprs = ReplacePart[#, 2 -> unshift[#[[2]], sh, x, z, Y]] & /@ logExprs];
  <|"Algebraic" -> Together[vExpr], "Logs" -> logExprs, "Y" -> Y|>];

tragerLogPart::usage = "tragerLogPart[ctx, r, Y] is the logarithmic part of the remainder r by Trager's construction, verified exactly, as log-term expressions in x and Y (Pipeline format). When the field of all residues is small (product of factorials of the factor degrees of R at most 8) the Q-basis construction runs first; otherwise, or when it declines, the per-orbit construction (TragerLogTerms.m); the Q-basis construction over the splitting field is the complete fallback when an orbit divisor is not torsion. A nonzero holomorphic difference proves non-elementarity ([T] p.71).";
tragerLogPart[ctx_Association, r_List, Y_Symbol] := Module[{Z = ctx["Z"], degs, crude, basisFirst, attempt},
  degs = Exponent[#[[1]], Z] & /@ Select[FactorList[ctx["R"]], !FreeQ[#[[1]], Z] &];
  crude = Times @@ (Factorial /@ degs);
  (* Genus 0: every orbit divisor is principal, so one orbit at a time over  *)
  (* Q[Z]/(q) always succeeds -- use the Q-basis construction first only for *)
  (* a small multiquadratic residue field (simplest logands).  Positive      *)
  (* genus: a single orbit divisor need not be torsion when the residues     *)
  (* satisfy Q-linear relations (gamma, omega gamma, omega^2 gamma), so the  *)
  (* complete Q-basis construction comes first unless the splitting field is *)
  (* hopeless; the orbit construction is then the fallback.                  *)
  basisFirst = If[ctx["genus"] === 0, Max[degs] <= 2 && crude <= 8, crude <= 5040];
  (* A certificate first when the field of all residues may be larger than   *)
  (* that of one residue: trace functionals live over Q(gamma) (it throws a  *)
  (* NonElementary certificate, or returns without a verdict).               *)
  If[ctx["genus"] =!= 0 && crude > Max[degs], traceCertificate[ctx]];
  attempt = If[basisFirst,
    catchTrager[tragerBasisLogPart[ctx, r, Y]],
    tragerFail["Decline", "skipped", <||>]];
  If[!tragerFailQ[attempt], Return[attempt]];
  If[attempt[[1]] =!= "Decline", Throw[attempt, "TragerFail"]];
  If[Max[degs] > 1,
    attempt = tragerOrbitLogPart[ctx, r, Y];
    If[attempt =!= $Failed, Return[attempt]]];
  If[basisFirst, Throw[tragerFail["Decline", "logarithmic part: no construction applies", <||>], "TragerFail"]];
  tragerBasisLogPart[ctx, r, Y]];

traceCertificate::usage = "traceCertificate[ctx] looks for a proof of non-elementarity that needs only the field Q(gamma) of ONE residue gamma, not the field of all of them. If the integral is elementary, r = Sum c_j du_j/u_j, and for every Q-linear lambda on the span of the residues Sum_i lambda(alpha_i) Delta(alpha_i) = Sum_j lambda(c_j) div(u_j) is a torsion divisor. lambda_gamma(alpha) = Tr(alpha gamma) is Q-linear and constant on the roots of each irreducible factor f_k of R over Q(gamma), so its divisor Sum_k lambda_k Delta(f_k) is defined over Q(gamma), of degree deg q instead of that of the splitting field. A residue gamma is taken from each irreducible factor q of R in turn; when its divisor is shown of infinite order (orders modulo two primes disagree, or the only possible order N fails torsionScreen) a NonElementary certificate is thrown. Otherwise (all lambda_gamma torsion, or zero) it returns Null without a verdict and the full construction follows.";
traceCertificate[ctx_Association] := Module[{Z = ctx["Z"], qs},
  qs = SortBy[Expand[#/Coefficient[#, Z, Exponent[#, Z]]] & /@ (First /@ Select[FactorList[ctx["R"]], !FreeQ[#[[1]], Z] &]),
    Exponent[#, Z] &];
  Do[traceCertificateAt[ctx, q, qs], {q, qs}]];

traceCertificateAt::usage = "traceCertificateAt[ctx, q, qs] is one step of traceCertificate: gamma a root of q (rational, or the generator of Q[t]/(q)), and the factors f_k of every q' in qs over Q(gamma). For each g = gamma^j (j = 0 .. deg q - 1) the functional lambda_g(alpha) = Tr(alpha g) is constant on the roots of f_k, lambda_k = Tr_{Q(gamma)/Q}(g s_k)/(deg q deg f_k) (s_k the root sum of f_k: the mean of the conjugates of alpha g, Tr_{K/Q} up to a common factor); g ranges over a basis of Q(gamma) because one power can vanish identically (power sums of R that are zero). Each nonzero lambda, up to proportionality, is scaled to integers and its divisor Sum_k lambda_k Delta(f_k) tested for torsion over Q(gamma).";
traceCertificateAt[ctx_Association, q_, qs_List] := Module[
  {Z = ctx["Z"], t = Unique["tg"], K, gamma, n, facs, lams, octx},
  n = Exponent[q, Z];
  If[n === 1,
    K = fieldQ[]; gamma = -Coefficient[q, Z, 0]; facs = qs,
    K = fieldNF[t, q /. Z -> t, None]; gamma = t;
    facs = Join @@ (factorOverField[K, #, Z] & /@ qs)];
  lams = Table[
    Table[With[{d = Exponent[f, Z]}, kTrace[K, gamma^j (-Coefficient[f, Z, d - 1])]/(n d)], {f, facs}],
    {j, 0, n - 1}];
  (* integral primitive representatives, zero and proportional ones dropped *)
  lams = DeleteDuplicates[Select[
    (With[{v = # LCM @@ (Denominator /@ #)}, v/(GCD @@ v) Sign[First[Select[v, # =!= 0 &]]]] &) /@
      Select[lams, !AllTrue[#, # === 0 &] &], True &]];
  octx = Join[ctx, <|"F" -> K, "residues" -> {}, "refc" -> None, "facs" -> facs|>];
  Do[traceDivisorTest[octx, <|"facs" -> facs, "Z" -> Z, "b" -> b, "ref" -> None|>, n], {b, lams}]];

traceDivisorTest::usage = "traceDivisorTest[octx, spec, n] tests the divisor spec of a trace functional for torsion over the field octx[\"F\"] (degree n): a NonElementary certificate when its orders modulo two good primes disagree (torsionOrderSpec) or its only possible order N is refuted at a further good prime (torsionScreen); Null without a verdict otherwise. An internal error is passed on.";
traceDivisorTest[octx_Association, spec_Association, n_Integer] := Module[{r},
  r = catchTrager[Module[{Nn = torsionOrderSpec[octx, spec]}, {Nn, torsionScreen[octx, spec, Nn]}]];
  Which[
    tragerFailQ[r] && r[[1]] === "NonElementary",
      throwNonElementary["a trace functional of the residues gives a divisor of infinite order over Q(gamma), deg " <>
        ToString[n] <> " (" <> r[[2]] <> ")", r[[3]]],
    tragerFailQ[r] && r[[1]] === "Internal", Throw[r, "TragerFail"],
    tragerFailQ[r], Null,                 (* a decline: no verdict, carry on *)
    r[[2]] === False,
      throwNonElementary["a trace functional of the residues gives a divisor of infinite order over Q(gamma), deg " <>
        ToString[n] <> ": its only possible order " <> ToString[r[[1]]] <> " is refuted at a further good prime"],
    True, Null]];

tragerBasisLogPart::usage = "tragerBasisLogPart[ctx, r, Y] is the Q-basis construction (TragerLogTerms.m, tragerLogTerms) over the field of all residues, verified exactly; a nonzero holomorphic difference proves non-elementarity.";
tragerBasisLogPart[ctx_Association, r_List, Y_Symbol] := Module[{B = ctx["B"], full, logs},
  full = Join[ctx, componentResidues[ctx["N"], ctx["E"], B["m"], B["G"], B["x"]]];
  logs = tragerLogTerms[full];
  If[!afZeroQ[full["F"], logSumResidual[full["F"], B, r, logs]],
    throwNonElementary["the remainder minus the logarithmic part is a nonzero holomorphic differential"]];
  tragerLogExpressions[full["F"], B, logs, Y]];

tragerOrbitLogPart::usage = "tragerOrbitLogPart[ctx, r, Y] is the per-orbit construction (tragerOrbitLogTerms), verified exactly, or $Failed when an orbit divisor is not torsion (not a certificate).";
tragerOrbitLogPart[ctx_Association, r_List, Y_Symbol] := Module[{orb = tragerOrbitLogTerms[ctx]},
  If[orb === $Failed, Return[$Failed]];
  If[!afZeroQ[fieldQ[], kauersResidual[ctx, r, orb]],
    throwNonElementary["the remainder minus the logarithmic part is a nonzero holomorphic differential"]];
  kauersLogExpressions[ctx, orb, Y]];

kauersLogPart::usage = "kauersLogPart[ctx, r, Y] is the logarithmic part of the remainder r by Kauers' heuristic, verified exactly, as log-term expressions in x and Y; a Decline when the heuristic does not close.";
kauersLogPart[ctx_Association, r_List, Y_Symbol] := Module[{logs = kauersLogTerms[ctx]},
  If[!afZeroQ[fieldQ[], kauersResidual[ctx, r, logs]],
    throwDecline["Kauers' heuristic: logarithmic part incomplete"]];
  kauersLogExpressions[ctx, logs, Y]];

fieldValue::usage = "fieldValue[F, a] is the explicit algebraic number of the element a of F (the generator t replaced by the embedding theta).";
fieldValue[F_Association, a_] := If[F["Type"] === "NF", kc[F, a] /. F["t"] -> F["theta"], a];

tragerLogExpressions::usage = "tragerLogExpressions[F, B, logs, Y] converts the log terms {{c_j, u_j}} (c_j in F, u_j coordinate lists) to {{coef, arg}} with explicit algebraic numbers and arg an expression in x and Y.";
tragerLogExpressions[F_Association, B_Association, logs_List, Y_Symbol] :=
  {fieldValue[F, #[[1]]], Together[fieldValue[F, afToY[B, #[[2]], Y]]]} & /@ logs;

integrateCurve::usage = "integrateCurve[c, x, m, G, method] is integrateComponent extended to the degenerate curve m = 1 (G = 1, a rational component, integrated directly). Returns the same Association.";
integrateCurve[c_, x_Symbol, 1, G_, method_String] :=
  <|"Algebraic" -> rationalIntegral[c G, x], "Logs" -> {}, "Y" -> Unique["Y"]|>;
integrateCurve[c_, x_Symbol, m_Integer, G_, method_String] := integrateComponent[c, x, m, G, method];

rationalIntegral::usage = "rationalIntegral[r, x] is the integral of the rational function r over Q (Integrate, which is exact on rational functions).";
rationalIntegral[r_, x_Symbol] := Module[{res},
  If[Together[r] === 0, Return[0]];
  res = Integrate[r, x];
  If[!FreeQ[res, Integrate], throwInternal["RationalIntegralUnevaluated", <|"r" -> r|>]];
  res];

componentExpression::usage = "componentExpression[res, sub] assembles the integral of one component from its Association: the curve generator res[\"Y\"] is replaced by sub (an expression in the radical symbol), and every logand and the algebraic part are passed through tidy.";
componentExpression[res_Association, sub_, tidy_] := Module[{Y = res["Y"]},
  tidy[res["Algebraic"] /. Y -> sub] +
    Total[If[Length[#] === 3,
        #[[1]] tidyRootSum[#[[2]] /. Y -> sub, tidy],
        #[[1]] Log[tidy[#[[2]] /. Y -> sub]]] & /@ res["Logs"]]];

tidyRootSum::usage = "tidyRootSum[RootSum[q, Function[body]], tidy] applies tidy to every logand inside the (held) body of the RootSum: the slot is replaced by a fresh symbol, the logands are tidied, and the slot is put back.";
tidyRootSum[RootSum[q_, Function[body_]], tidy_] := Module[{s = Unique["s"], nb},
  nb = (body /. Slot[1] -> s) /. Log[a_] :> Log[tidy[a]];
  With[{b2 = nb /. s -> Slot[1]}, RootSum[q, Function[b2]]]];
tidyRootSum[e_, tidy_] := e;

reduceRadicalPowers::usage = "reduceRadicalPowers[e, yh, n, g] rewrites the rational function e of x and the symbol yh, where yh^n = g, with numerator and denominator reduced modulo yh^n - g (every power of yh below n).";
reduceRadicalPowers[e_, yh_Symbol, n_Integer, g_] := Module[{t = Together[e]},
  Together[reduceYDegree[Numerator[t], yh, n, g]/reduceYDegree[Denominator[t], yh, n, g]]];

TragerIntegrateCurve[a_, x_Symbol, m_Integer, G_, method_String : "Trager"] := Module[{r},
  r = catchTrager[integrateComponent[a, x, m, G, method]];
  If[tragerFailQ[r], failToList[r], r]];

TragerIntegrate[f_, x_Symbol, method_String : "Trager"] := Module[{r},
  r = catchTrager[tragerIntegrateBody[f, x, method]];
  If[tragerFailQ[r], failToList[r], r]];

tragerIntegrateBody::usage = "tragerIntegrateBody[f, x, method] is the body of TragerIntegrate: it throws instead of returning failure lists.";
tragerIntegrateBody[f_, x_Symbol, method_String] := Module[
  {cf, n, g, yh, yexpr, coeffs, total, cc, res, Yexpr, rho, part},
  If[!MemberQ[{"Trager", "Kauers"}, method], throwDecline["unknown method " <> method]];
  cf = canonicalForm[f, x];
  n = cf["n"]; coeffs = cf["coeffs"];
  If[n === 1, throwDecline["no radical in the integrand (a rational function)"]];
  g = cf["g"]; yh = cf["y"]; yexpr = radicalExpression[cf];
  total = rationalIntegral[coeffs[[1]], x];
  Do[
    If[coeffs[[i + 1]] =!= 0,
      cc = componentCurve[n, g, i, x];
      res = integrateCurve[Together[coeffs[[i + 1]] cc["s"]], x, cc["m"], cc["G"], method];
      (* The generator in terms of the user's radical.  On an exact curve Y = *)
      (* yh^i/s is expressed in the symbol yh and every power of yh reduced    *)
      (* modulo yh^n = g before yh is replaced by the radical, which keeps the *)
      (* logands in the lowest powers of the user's radical.                 *)
      If[cc["exact"],
        rho = 1;
        part = componentExpression[res, yh^i/cc["s"], reduceRadicalPowers[#, yh, n, g] &] /. yh -> yexpr,
        Yexpr = If[cc["m"] === 1, 1, cc["G"]^(1/cc["m"])];
        rho = (yexpr^i/cc["s"])/Yexpr;
        part = componentExpression[res, Yexpr, Together]];
      total = total + rho part],
    {i, 1, n - 1}];
  total = evaluateRootSumBodies[total];
  spotCheck[f, total, x];
  presentable[f, total, x]];

evaluateRootSumBodies::usage = "evaluateRootSumBodies[e] evaluates the (held) body of every RootSum in e, so that substitutions made into it -- the user's radical for the curve generator -- are simplified like the rest of the answer.";
evaluateRootSumBodies[e_] := e /. RootSum[q_, Function[b_]] :> Module[{s = Unique["s"], b1},
  b1 = b /. Slot[1] -> s;
  With[{b2 = b1 /. s -> Slot[1]}, RootSum[q, Function[b2]]]];

presentable::usage = "presentable[f, F, x] is the real form of the verified antiderivative F (LogToArcTan.m) when that passes antiderivativeQ, else F itself. The real form is not attempted on answers above tragerOption[\"RealFormMaxLeaves\"] leaves and is given tragerOption[\"RealFormTime\"] seconds.";
presentable[f_, Fx_, x_Symbol] := Module[{real},
  (* The real form is presentation only (Fx is already verified): skip it on *)
  (* large answers and bound its time, rather than spend minutes in Together *)
  (* over nested radicals.                                                   *)
  If[FreeQ[Fx, Complex] && FreeQ[Fx, RootSum], Return[Quiet[realForm[Fx, x]]]];
  If[LeafCount[Fx] > tragerOption["RealFormMaxLeaves"], Return[Fx]];
  real = TimeConstrained[Quiet[realForm[Fx, x]], tragerOption["RealFormTime"], Fx];
  If[real =!= Fx && antiderivativeQ[f, real, x], real, Fx]];

spotCheck::usage = "spotCheck[f, F, x] confirms that D[F, x] = f by PossibleZeroQ[D[F, x] - f] alone -- no numeric fallback: a correct answer PossibleZeroQ cannot confirm is a PossibleZeroQ bug, fixed there (v0.6091). The per-component results were already verified exactly in K(x)(Y), so this guards the reassembly into the user's radical; a failure is an internal error (decline), never an answer.";
spotCheck[f_, Fx_, x_Symbol] := If[antiderivativeQ[f, Fx, x], True,
  throwInternal["AntiderivativeCheckFailed", <|"f" -> f, "F" -> Fx|>]];

antiderivativeQ::usage = "antiderivativeQ[f, F, x] is True when PossibleZeroQ[D[F, x] - f, Assumptions -> Element[x, Complexes]] is True: F must be an antiderivative for all complex x, so x is sampled at complex points too.";
antiderivativeQ[f_, Fx_, x_Symbol] := TrueQ[Quiet[PossibleZeroQ[D[Fx, x] - f, Assumptions -> Element[x, Complexes]]]];

