(* ::Package:: *)

(* ========================================================================== *)
(* Common.m -- the failure protocol and the package options.                  *)
(*                                                                            *)
(* Bronstein uses Trager's protocol unchanged (trager/Common.m): a stage that  *)
(* cannot continue throws tragerFail[class, reason, info] with tag            *)
(* "TragerFail" via throwNonElementary / throwDecline / throwInternal, and the *)
(* entry point catches it with catchTrager and maps it with failToList.  The  *)
(* rule that keeps certificates sound: throwNonElementary only where a        *)
(* theorem of [B] or [T] proves non-elementarity AND every sub-decision it     *)
(* rests on was itself certified (a Risch DE "NoSolution" from                *)
(* Risch`RischDECertified, never a $Failed).                                  *)
(* ========================================================================== *)

brOption::usage = "brOption[name] is the default value of an internal setting of the Bronstein package. \"SpecializationTries\": number of candidate specializations per round of the growing-height search; \"SpecializationHeight\": bound on the numerators of the first round of rationals substituted for x and the lower monomials (it doubles every round); \"SpecializationRounds\": rounds of the growing-height search before an Internal error; \"PresentationTime\": seconds for the real-form presentation and its re-verification (past it the exactly verified form is returned); \"ArcTanFormTime\": seconds ParallelMixedTower may take on the real ArcTan form of a transcendental part (brIntegrateArcTanForm).";
brOption["SpecializationTries"] = 12;
brOption["SpecializationHeight"] = 97;
brOption["SpecializationRounds"] = 8;
brOption["PresentationTime"] = 5;
brOption["SupportSplitLeaves"] = 20000;
brOption["SupportMergeTime"] = 4;
brOption["KernelSubstitutionTime"] = 3;   (* brKernelSubstitution: the algebraic integral in u = theta, theta a logarithm *)
brOption["KernelSubstitutionExpTime"] = 0.6;   (* ... theta an exponential: every 257 cM success closes in <= 0.31 s; cM-4-306, which Trager in u needs 2.5 s for, falls to the 1.8 s tower sooner *)
brOption["ArcTanFormTime"] = 20;

brCatch::usage = "brCatch[expr] evaluates expr and returns its value, or the tragerFail[...] object thrown by any stage inside it (Trager's catchTrager).";
SetAttributes[brCatch, HoldFirst];
brCatch[expr_] := catchTrager[expr];

brFailQ::usage = "brFailQ[e] is True when e is a tragerFail[...] object.";
brFailQ[e_] := tragerFailQ[e];

brFailTag::usage = "brFailTag[e] is {class, reason} of the tragerFail object e (class \"NonElementary\", \"Decline\" or \"Internal\"), or None when e is not one. Used by the unit suites to assert the exact failure.";
brFailTag[tragerFail[class_, reason_, _]] := {class, reason};
brFailTag[_] := None;

brWithParams::usage = "brWithParams[vars] is vars followed by the symbolic parameters tragerParams not already in it: the transcendental symbols of a coefficient field over the constant field Q(params) (Phase 4.4).";
brWithParams[vars_List] := Join[vars, Complement[tragerParams, vars]];

brFn::usage = "brFn[vars] is fieldFn over the symbols vars and the symbolic parameters (brWithParams). brFn[vars, t, mu] adjoins a root t of the monic mu(t): its embedding theta is the first Root of mu when mu is free of the parameters, else None (an element is then rendered over the orbit, never by its numeric value; Phase 4.4).";
brFn[vars_List] := fieldFn[brWithParams[vars]];
brFn[vars_List, t_Symbol, mu_] := fieldFn[brWithParams[vars],
  If[tragerParams === {} || FreeQ[mu, Alternatives @@ tragerParams],
    {t, mu, Root[Function[Evaluate[mu /. t -> #]], 1]},
    With[{mm = Collect[Together[mu/Coefficient[mu, t, Exponent[mu, t]]], t, Together]}, {t, mm, brParamRoot[mm, t]}]]];

brParamRoot::usage = "brParamRoot[mu, t] is an explicit root of the monic irreducible mu(t) over Q(params), in radicals, for the presentation of an element of Q(params)[t]/(mu) (any root will do: the identities of the field hold at every root, and the answer is verified): the root of a linear mu, c^(1/n) for a binomial t^n - c, the root (-b + Sqrt[b^2 - 4 c])/2 of a quadratic t^2 + b t + c; None otherwise (Phase 4.4).";
brParamRoot[mu_, t_Symbol] := Module[{n = Exponent[mu, t], cl = CoefficientList[mu, t]},
  Which[
    n === 1, Together[-cl[[1]]],
    n >= 2 && AllTrue[cl[[2 ;; n]], Together[#] === 0 &], Together[-cl[[1]]]^(1/n),
    n === 2, (-cl[[2]] + Sqrt[Together[cl[[2]]^2 - 4 cl[[1]]]])/2,
    True, None]];
