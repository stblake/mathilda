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

brOption::usage = "brOption[name] is the default value of an internal setting of the Bronstein package. \"SpecializationTries\": number of random specializations tried before a good one is declared unavailable (a decline, never a certificate); \"SpecializationHeight\": bound on the numerators/denominators of the random rationals substituted for x and the lower monomials; \"PresentationTime\": seconds for the real-form presentation and its re-verification (past it the exactly verified form is returned).";
brOption["SpecializationTries"] = 12;
brOption["SpecializationHeight"] = 97;
brOption["PresentationTime"] = 5;

brCatch::usage = "brCatch[expr] evaluates expr and returns its value, or the tragerFail[...] object thrown by any stage inside it (Trager's catchTrager).";
SetAttributes[brCatch, HoldFirst];
brCatch[expr_] := catchTrager[expr];

brFailQ::usage = "brFailQ[e] is True when e is a tragerFail[...] object.";
brFailQ[e_] := tragerFailQ[e];

brFailTag::usage = "brFailTag[e] is {class, reason} of the tragerFail object e (class \"NonElementary\", \"Decline\" or \"Internal\"), or None when e is not one. Used by the unit suites to assert the exact failure.";
brFailTag[tragerFail[class_, reason_, _]] := {class, reason};
brFailTag[_] := None;
