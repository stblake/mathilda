(* ::Package:: *)

(* ========================================================================== *)
(* KauersLogTerms.m -- the logarithmic part by Kauers' Groebner-basis          *)
(*                     heuristic, [K] (ISSAC 2008).                            *)
(*                                                                            *)
(* Remainder r = (N/E) Y dx on Y^m = G.  For an irreducible factor q(Z) of the *)
(* residue polynomial, the ideal                                              *)
(*     J_q = < Y^m - G(x), Bp(x), N(x) Y - Z E'(x), q(Z) >   in Q[x, Y, Z]     *)
(* describes, for each root gamma of q, the affine places with residue gamma  *)
(* (the analogue of gcd(D, A - gamma D') in the rational case).  The          *)
(* heuristic looks for the least k such that J_q^k is principal modulo the    *)
(* curve and q:                                                               *)
(*     J_q^k + <Y^m - G, q> = < p, Y^m - G, q >                                *)
(* for an element p of the Groebner basis of J_q^k.  Then p(x, Y, gamma) has   *)
(* divisor k Delta(gamma) on the affine curve, and the candidate is           *)
(*     Sum_q (1/k_q) Sum_{q(gamma) = 0} gamma Log[p_q(x, Y, gamma)]            *)
(*         = Sum_q (1/k_q) RootSum[q, # Log[p_q(x, Y, #)] &].                  *)
(* Nothing guarantees the poles at infinity cancel, or that a power is        *)
(* principal at all, so the candidate is verified exactly (by traces over     *)
(* Q[Z]/(q)) and the method declines when it does not close.  Groebner bases  *)
(* are lexicographic in x > Y > Z, an elimination order for x and Y (the      *)
(* block order of [K] is not needed).                                         *)
(* ========================================================================== *)

kauersLogTerms::usage = "kauersLogTerms[ctx] (ctx[\"Ysym\"] the symbol for the curve generator) runs Kauers' heuristic on every irreducible factor q of the residue polynomial of the component ctx. It returns {{q, k, p}, ...}: q(Z) the factor, k the ideal power, p(x, Y, Z) the logand. Throws a Decline when some factor has no principal power up to tragerOption[\"KauersMaxPower\"].";
kauersLogTerms[ctx_Association] := Module[{B = ctx["B"], x, Y = ctx["Ysym"], Z = ctx["Z"], m, curve, facs, out = {}},
  x = B["x"]; m = B["m"];
  curve = Y^m - B["G"];
  facs = First /@ Select[tragerFactorList[ctx["R"]], !FreeQ[#[[1]], Z] &];
  Do[
    Module[{base, J, Jk, k = 1, p = None},
      base = {curve, q};
      J = groebnerLex[{curve, ctx["Bp"], ctx["N"] Y - Z ctx["Ep"], q}, {x, Y, Z}];
      Jk = J;
      While[p === None && k <= tragerOption["KauersMaxPower"],
        p = kauersPrincipal[Jk, base, {x, Y, Z}];
        If[p === None,
          k++;
          Jk = groebnerLex[Join[Flatten[Outer[Times, Jk, J]], base], {x, Y, Z}]]];
      If[p === None, throwDecline["Kauers' heuristic: no principal power of the residue ideal"]];
      AppendTo[out, {q, k, p}]],
    {q, facs}];
  out];

groebnerLex::usage = "groebnerLex[polys, vars] is the reduced lexicographic Groebner basis of polys (vars ordered from highest), each element made monic, sorted canonically so that bases can be compared with ===.";
groebnerLex[polys_List, vars_List] := Module[{gb = GroebnerBasis[polys, vars]},
  Sort[Expand[#/leadingCoefficientLex[#, vars]] & /@ gb]];

leadingCoefficientLex::usage = "leadingCoefficientLex[p, vars] is the coefficient of the lexicographically leading monomial of p.";
leadingCoefficientLex[p_, vars_List] := Module[{rules = CoefficientRules[Expand[p], vars]},
  If[rules === {}, 1, Last[Last[SortBy[rules, First]]]]];

kauersPrincipal::usage = "kauersPrincipal[Jk, base, vars] looks for p in the Groebner basis Jk with GroebnerBasis[base U {p}] equal to Jk (the ideal is principal modulo base); returns p or None.";
kauersPrincipal[Jk_List, base_List, vars_List] := Module[{cands},
  cands = SortBy[Complement[Jk, groebnerLex[base, vars]], LeafCount];
  SelectFirst[cands, groebnerLex[Append[base, #], vars] === Jk &, None]];

kauersResidual::usage = "kauersResidual[ctx, r, logs] is r - Sum_q (1/k_q) Tr_{Q[Z]/(q)}(Z p_q'/p_q) as an element over Q: zero exactly when the Kauers candidate integrates the remainder r.";
kauersResidual[ctx_Association, r_List, logs_List] := Module[{B = ctx["B"], acc = r},
  Do[
    acc = afMinus[baseField[], acc, afScale[baseField[], 1/lg[[2]], kauersTrace[ctx, lg[[1]], lg[[3]]]]],
    {lg, logs}];
  acc];

kauersTrace::usage = "kauersTrace[ctx, q, p] is Sum_{q(gamma) = 0} gamma p'(gamma)/p(gamma) as an element of Q(x)(Y): the trace from K = Q[Z]/(q) of Z p'/p, computed exactly (inverse of p by its adjugate, denominator rationalised by its norm, trace of Z^j by the companion matrix).";
kauersTrace[ctx_Association, q_, p_] := If[Exponent[q, ctx["Z"]] === 1,
  kauersTraceRational[ctx, q, p],
  kauersTraceNF[ctx, q, p]];

kauersTraceRational::usage = "kauersTraceRational[ctx, q, p] is kauersTrace for a linear factor q: the value gamma p'/p at its rational root gamma.";
kauersTraceRational[ctx_Association, q_, p_] := Module[{B = ctx["B"], Z = ctx["Z"], g0, pp},
  g0 = Z /. Solve[q == 0, Z][[1]];
  pp = afFromY[baseField[], B, p /. Z -> g0, ctx["Ysym"]];
  afScale[baseField[], g0, afTimes[baseField[], B, afD[baseField[], B, pp], afInverse[B, pp]]]];

kauersTraceNF::usage = "kauersTraceNF[ctx, q, p] is kauersTrace for an irreducible factor q of degree >= 2, computed over K = Q[Z]/(q). p may be a rational function of x, Y and Z: Log of numerator minus Log of denominator.";
kauersTraceNF[ctx_Association, q_, p_] := Module[
  {B = ctx["B"], Z = ctx["Z"], K, tp = Together[p], L, comp, d, tr, rt},
  (* over a number field of constants (q with radical or Root coefficients):
     the trace through the absolute field (kauersTraceRel), not the relative
     presentation Q[Z]/(q) with literal algebraic constants, whose norms and
     exact quotients ran as symbolic radical arithmetic (G-perf-4) *)
  rt = kauersTraceRel[ctx, q, p];
  If[rt =!= None, Return[rt, Module]];
  K = extField[Z, tragerMonic[q, Z]];
  L = afScale[K, Z, afMinus[K, logDerivativeNF[K, B, Numerator[tp], ctx["Ysym"]],
                               logDerivativeNF[K, B, Denominator[tp], ctx["Ysym"]]]];
  d = Exponent[K["mu"], Z];
  comp = companionMatrix[K["mu"], Z];
  tr = Table[Tr[MatrixPower[comp, j]], {j, 0, d - 1}];
  Table[
    With[{nd = rcNumDen[K, L[[i]]]},
      Together[(PadRight[CoefficientList[nd[[1]], Z], d] . tr)/nd[[2]]]],
    {i, B["m"]}]];

logDerivativeNF::usage = "logDerivativeNF[K, B, P, Y] is P'/P as an element over K (coordinates with prime-field denominators) for P a POLYNOMIAL in x, Y and the field generator: the inverse of P is the first row of the adjugate of its multiplication matrix over its determinant, and that determinant is rationalised by its norm.";
logDerivativeNF[K_Association, B_Association, P_, Y_Symbol] := Module[
  {x = B["x"], pa, M, det, adj1, nrm, cof, pinv, dN, dD},
  If[FreeQ[P, x] && FreeQ[P, Y], Return[afZero[B]]];
  pa = afFromY[K, B, P, Y];
  M = afMultMatrix[K, B, pa];
  (* det = dN/dD with dN, dD polynomials over K: 1/det = dD (Norm dN/dN)/Norm dN,
     the exact quotient taken on the POLYNOMIAL numerator.  The integral basis
     can give det a denominator ((1 + x)^2 for 1/((x^2 + I) Sqrt[x + 1])), and
     schoolbook division of the rational Norm det over a Gaussian residue field
     never terminated (G-adv-t2) *)
  {dN, dD} = rcNumDen[K, Det[M]];
  adj1 = If[B["m"] === 1, {1}, Table[kc[K, (-1)^(1 + j) Det[Drop[M, {j}, {1}]]], {j, B["m"]}]];
  nrm = kNorm[K, dN];
  cof = kExactQuo[K, nrm, dN, x];
  pinv = rc[K, # dD cof/nrm] & /@ adj1;
  afTimes[K, B, afD[K, B, pa], pinv]];

kauersTraceRel::usage = "kauersTraceRel[ctx, q, p] is kauersTraceNF for a factor q whose coefficients carry algebraic constants (radicals, Root objects) while tragerAlgConstants is in use: Tr_{K/K0}(gamma p'/p) computed in the ABSOLUTE field K = K0(gamma) of fieldOfRoot (number-field kernels, a rational norm) and taken down to K0 = Q(constants) by relTrace. The log argument p, written relative to K0 (a polynomial in the root symbol Z, relToZ), is mapped to K by Z -> gamma and the constants' representations (fieldSubI). None when q has no algebraic constant, or the absolute field has no relative structure (the caller keeps the relative presentation).";
kauersTraceRel[ctx_Association, q_, p_] := Module[
  {B = ctx["B"], Z = ctx["Z"], qm, fr, F, gamma, tp, toF, L},
  If[tragerAlgConstants === {} ||
     FreeQ[q, (Power[b_, _Rational] /; rationalNumberQ[b]) | _Root], Return[None, Module]];
  qm = tragerMonic[q, Z];
  fr = kauersFieldOfRoot[qm, Z];
  F = fr["F"]; gamma = fr["gamma"];
  If[!KeyExistsQ[F, "rel"] || F["Type"] =!= "NF", Return[None, Module]];
  tp = Together[p];
  (* every algebraic constant in the generator: radicals by fieldSubI, and
     Root objects (Root[-16 + #^3 &, 1] = 2 2^(1/3), as relToZ writes theta0's
     powers) by fieldConstRep -- one left literal keeps the norm out of FLINT
     (symbolic Root arithmetic again); one that F cannot represent: None *)
  toF[e_] := Module[{r = Expand[e /. Z -> gamma], roots, reps},
    roots = DeleteDuplicates[Cases[r, _Root, {0, Infinity}]];
    reps = Catch[(# -> fieldConstRep[F, #]) & /@ roots, "TragerFail"];
    If[!ListQ[reps] || !FreeQ[reps[[All, 2]], Root], Throw[None, "kauersTraceRel"]];
    Expand[fieldSubI[F, r /. reps]]];
  Catch[
    L = afScale[F, gamma, afMinus[F, logDerivativeNF[F, B, toF[Numerator[tp]], ctx["Ysym"]],
                                     logDerivativeNF[F, B, toF[Denominator[tp]], ctx["Ysym"]]]];
    (* K0's generator comes back from relTrace as the Root it is in the field
       (16^(1/3) as Root[-16 + #^3 &, 1]); written as the radical, so that the
       residual against the integrand (spelled in radicals) is in one
       spelling *)
    tragerTogether[relTrace[F, #] /. r_Root :> With[{v = binomialRootRadical[r]}, If[v === None, r, v]]] & /@ L,
    "kauersTraceRel"]];

kauersFieldOfRoot::usage = "kauersFieldOfRoot[q, Z] is fieldOfRoot[q, Z, Z], memoised for the current constants (the orbit construction built the same field moments before).";
kauersFieldOfRoot[q_, Z_Symbol] := Module[{key = {q, Z, tragerAlgConstants}},
  If[KeyExistsQ[kauersFieldMemo, key], kauersFieldMemo[key],
    If[Length[kauersFieldMemo] > 32, kauersFieldMemo = <||>];
    kauersFieldMemo[key] = fieldOfRoot[q, Z, Z]]];
kauersFieldMemo = <||>;

companionMatrix::usage = "companionMatrix[mu, Z] is the companion matrix of the monic polynomial mu(Z); Tr of its j-th power is the j-th power sum of the roots.";
companionMatrix[mu_, Z_Symbol] := Module[{d = Exponent[mu, Z], c = CoefficientList[mu, Z]},
  Table[Which[j === d, -c[[i]], i === j + 1, 1, True, 0], {i, d}, {j, d}]];

kauersLogExpressions::usage = "kauersLogExpressions[ctx, logs, Y] converts the Kauers log terms {{q, k, p}} into expressions in x and Y: (gamma/k) Log[p] for a linear q, (1/k) RootSum[q, # Log[p] &] otherwise.";
kauersLogExpressions[ctx_Association, logs_List, Y_Symbol] := Module[{Z = ctx["Z"]},
  Map[
    Function[lg,
      With[{q = lg[[1]], k = lg[[2]], p = lg[[3]] /. ctx["Ysym"] -> Y},
        If[Exponent[q, Z] === 1,
          With[{g0 = Z /. Solve[q == 0, Z][[1]]}, {g0/k, Together[p /. Z -> g0]}],
          With[{qf = q /. Z -> Slot[1], body = Slot[1] Log[p /. Z -> Slot[1]]},
            {1/k, RootSum[Function[qf], Function[body]], "RootSum"}]]]],
    logs]];
