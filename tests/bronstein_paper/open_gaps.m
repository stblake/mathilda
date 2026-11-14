(* ::Package:: *)
(* ========================================================================== *)
(* tests/bronstein_paper/open_gaps.m -- the register of KNOWN OPEN GAPS of the *)
(* Bronstein completion campaign (and the lower-layer defects it found), each  *)
(* with an example and its CURRENT behaviour, so follow-up work can find them  *)
(* and a fix is noticed: tests/bronstein_paper/open_gaps_run.m reports every  *)
(* entry STILL-OPEN (behaviour as recorded) or CHANGED (closed or moved: move *)
(* it to examples.m with its fix, or update the record).                      *)
(*                                                                            *)
(*   "id"       stable label                                                  *)
(*   "area"     Bronstein | Trager | Risch (C) | Core | Cascade | Presentation *)
(*   "gap"      what is missing or wrong                                      *)
(*   "kind"     "Integrate" (BronsteinIntegrate[f, x]) or "Expr" (f is any    *)
(*              held expression; "current" is its InputForm today)            *)
(*   "f"        held input                                                    *)
(*   "current"  Integrate: "incomplete" | "timeout" | "internal" | "slow"     *)
(*              (answers, but over the budget) | "presentation" (answers,     *)
(*              verified, but the form needs work); Expr: the InputForm string *)
(*   "budget"   seconds for the check (a "timeout"/"slow" entry is over it)   *)
(*   "followup" the planned work (plan phase or fix)                          *)
(*   "source"   where it was found                                            *)
(* Gaps without a single example (the ratchet sites, whole plan phases) are   *)
(* listed in openGapsStructural below.                                        *)
(* ========================================================================== *)

bronsteinOpenGaps = {

(* ---- performance --------------------------------------------------------- *)
<|"id" -> "G-perf-4", "area" -> "Trager", "gap" -> "a residue field over Q(constants) of large degree: x/((2^(1/3) x^3 - 1) Sqrt[x^2 + 2]) (genus 0, elementary) times out, where x/((2 x^3 - 1) Sqrt[x^2 + 2]) takes 0.1 s and the linear 2^(1/3) x - 1 takes 2.6 s",
  "kind" -> "Integrate", "f" -> Hold[x/((2^(1/3) x^3 - 1) Sqrt[x^2 + 2])], "current" -> "timeout", "budget" -> 60,
  "followup" -> "the residue field of an orbit over K0 is relative of small degree; build it as a tower over K0 (or per orbit) instead of the absolute splitting field joined with the constants", "source" -> "Phase 4.3b C4 battery"|>,
<|"id" -> "G-perf-1", "area" -> "Bronstein", "gap" -> "a three-radical compositum (Kummer group Z/30): Trager on the degree-30 curve is too slow",
  "kind" -> "Integrate", "f" -> Hold[Sqrt[x] (x + 1)^(1/5) (x + 2)^(1/3)], "current" -> "timeout", "budget" -> 60,
  "followup" -> "profile Trager on high-degree Kummer curves; smaller presentations of the compositum", "source" -> "gap log (E6b), 2026-11-14"|>,
<|"id" -> "G-perf-2", "area" -> "Bronstein", "gap" -> "two log terms with special residues on a genus-1 cubic: 225 s",
  "kind" -> "Integrate", "f" -> Hold[Together[D[Log[Log[x] + (1 + Log[x] + Log[x]^3)^(1/3)] + 2 Log[Log[x]^2 + Log[x] (1 + Log[x] + Log[x]^3)^(1/3) + (1 + Log[x] + Log[x]^3)^(2/3)], x]]],
  "current" -> "timeout", "budget" -> 60, "followup" -> "profile the E3 torsion / principal-generator steps", "source" -> "E3 tests, 2026-11-14"|>,
<|"id" -> "G-perf-3", "area" -> "Bronstein", "gap" -> "the fourth root of the tangent times out (13 of 120 random m = 3, 4 components also time out at 30 s)",
  "kind" -> "Integrate", "f" -> Hold[Tan[x]^(1/4)], "current" -> "timeout", "budget" -> 60,
  "followup" -> "profile m = 4 over Q(i); the random-component sweep of E6e", "source" -> "review of E4.2; E6e sweep"|>,

(* ---- Bronstein scope ----------------------------------------------------- *)
<|"id" -> "G-sign-1", "area" -> "Bronstein", "gap" -> "radicals equal up to a sign seen only on a branch cut: undecided (E6d)",
  "kind" -> "Integrate", "f" -> Hold[(1 + E^(x^2))/2 + (1 - E^(x^2))/2 Sqrt[x + 1] Sqrt[1/(x + 1)]], "current" -> "incomplete", "budget" -> 60,
  "followup" -> "prove eps = 1 on every open set (no open region with eps = -1) and integrate f_+ alone", "source" -> "review of E6d"|>,
<|"id" -> "G-trig-1", "area" -> "Bronstein", "gap" -> "TrigToExp writes ArcCosh with Sqrt[x - 1] Sqrt[x + 1]: two radicals (E6d over a constant-sign split)",
  "kind" -> "Integrate", "f" -> Hold[ArcCosh[x]/Sqrt[x^2 - 1]], "current" -> "incomplete", "budget" -> 60,
  "followup" -> "a front-end rewrite Sqrt[x - 1] Sqrt[x + 1] = eps Sqrt[x^2 - 1] (E6d)", "source" -> "review of E4.2"|>,
<|"id" -> "G-qi-1", "area" -> "Bronstein", "gap" -> "residue at infinity of a logarithm: the integral in k has a logarithm with non-rational coefficient (Q(i))",
  "kind" -> "Integrate", "f" -> Hold[Sqrt[I + Log[x]^2]], "current" -> "incomplete", "budget" -> 60,
  "followup" -> "brLimitedIntegrate over Q(i) with algebraic log coefficients (Phase 3 LimitedIntegrate)", "source" -> "Phase 4.1b battery"|>,
<|"id" -> "G-nf-2", "area" -> "Bronstein", "gap" -> "radicals equal over a larger constant field (Sqrt[2 x] = Sqrt[2] Sqrt[x])",
  "kind" -> "Integrate", "f" -> Hold[Sqrt[2 x] + Sqrt[x]], "current" -> "incomplete", "budget" -> 60,
  "followup" -> "Phase 4.3 (E6b/E6d over a number field)", "source" -> "E6b"|>,
<|"id" -> "G-red-1", "area" -> "Bronstein", "gap" -> "two reducible logarithms: kappa values may be algebraically dependent (E6a needs one)",
  "kind" -> "Integrate", "f" -> Hold[(Log[x^2] - 2 Log[x]) (Log[x^3] - 3 Log[x]) Sqrt[Log[x]]], "current" -> "incomplete", "budget" -> 60,
  "followup" -> "E6a for several kappa (their values are logs of algebraic numbers: Baker's theorem on linear forms)", "source" -> "E6a"|>,
<|"id" -> "G-mid-1", "area" -> "Bronstein", "gap" -> "a radical in the middle of the tower",
  "kind" -> "Integrate", "f" -> Hold[Log[1 + Sqrt[Log[x]]]/x], "current" -> "incomplete", "budget" -> 60,
  "followup" -> "Phase 5 (E1, E2)", "source" -> "plan"|>,
<|"id" -> "G-param-1", "area" -> "Bronstein", "gap" -> "symbolic parameters",
  "kind" -> "Integrate", "f" -> Hold[Sqrt[a + Log[x]]/x], "current" -> "incomplete", "budget" -> 60,
  "followup" -> "Phase 4.4 (E4: generic verdicts)", "source" -> "plan"|>,

(* ---- presentation -------------------------------------------------------- *)
<|"id" -> "G-pres-1", "area" -> "Presentation", "gap" -> "answers over Q(i) carry Root objects of the residue field (no radicals, no real form)",
  "kind" -> "Integrate", "f" -> Hold[Sqrt[Tan[x]]], "current" -> "presentation", "budget" -> 60,
  "followup" -> "ToRadicals on residue-field generators; a verified real form for Q(i) answers of real integrands", "source" -> "Phase 4.2"|>,

(* ---- lower layers (found by the campaign's reviews) --------------------- *)
<|"id" -> "G-core-8", "area" -> "Core", "gap" -> "FactorSquareFreeList without Extension does not split a polynomial with algebraic coefficients at all (z^2 - 2^(1/3) z^3 is reported squarefree); Trager passes Extension -> Automatic (tragerSquareFreeList)",
  "kind" -> "Expr", "f" -> Hold[Length[FactorSquareFreeList[z^2 - 2^(1/3) z^3]]], "current" -> "2", "budget" -> 60,
  "followup" -> "the default path should treat algebraic numbers as independent generators (split z^2); lifting them as indeterminates in poly.c (as for Root) rewrote opaque nested radicals, so it was withdrawn", "source" -> "Phase 4.3b"|>,
<|"id" -> "G-core-9", "area" -> "Core", "gap" -> "Together without Extension leaves a sum of fractions uncombined when their common denominator needs Sqrt[2]^2 = 2 (Trager now passes Extension -> Automatic)",
  "kind" -> "Expr", "f" -> Hold[Head[Together[3/(2 z^3 + 2 Sqrt[2] z^6) + (-3 + 3 z + z^2 - 3 Sqrt[2] z^3 - 2 Sqrt[2] z^5)/(2 z^3 - 2 z^4 + 4 Sqrt[2] z^6 - 2 Sqrt[2] z^7 + 4 z^9)]]], "current" -> "Plus", "budget" -> 60,
  "followup" -> "combine over Q(constants) when the divide-back over Q[constants] fails", "source" -> "Phase 4.3b"|>,
<|"id" -> "G-core-4", "area" -> "Risch (C)", "gap" -> "Log[x^2] is not its own monomial in the C transcendental integrator (the retry fixed in v0.6547 no longer answers wrongly)",
  "kind" -> "Expr", "f" -> Hold[Head[Integrate`RischTranscendental[Log[x^2]/x, x]]], "current" -> "Integrate`RischTranscendental", "budget" -> 60,
  "followup" -> "tower builder: a logarithm of a power as a monomial, or the kappa split of E6a in C", "source" -> "Phase 3 part 1"|>,
(* ---- found by the adversarial reviews of Bronstein, Trager and PMT (2026-11-14) *)
<|"id" -> "G-adv-b5", "area" -> "Bronstein", "gap" -> "an inverse trig times a logarithm times out (elementary: D[ArcSin[x] Log[x]])",
  "kind" -> "Integrate", "f" -> Hold[Together[D[ArcSin[x] Log[x], x]]], "current" -> "timeout", "budget" -> 60,
  "followup" -> "profile; the kappa split and the trig front end on reducible logs", "source" -> "adversarial review of Bronstein"|>,
<|"id" -> "G-adv-c1", "area" -> "Risch (C)", "gap" -> "a transcendental constant (Log[2]) beside the tower: the C decider's no is now demoted to undecided, so this correct certificate is lost",
  "kind" -> "Expr", "f" -> Hold[Head[Quiet[Risch`ElementaryIntegralQ[Log[2]/(1 + u^2) + Log[u]/(1 + u^2), u]]]], "current" -> "Risch`ElementaryIntegralQ", "budget" -> 60,
  "followup" -> "treat transcendental constants as parameters of the constant field (Phase 4.4)", "source" -> "adversarial review of Bronstein"|>,
<|"id" -> "G-adv-p1", "area" -> "ParallelMixedTower", "gap" -> "the tower builder does not prove its generators independent (Gaussian exponent combinations, an exponential beside a flattened radical, relations lost on rebuild); certificates now need independent confirmation, which drops PMT's unique ones",
  "kind" -> "Expr", "f" -> Hold[Quiet[Integrate`ParallelMixedTower[(I E^(I x) + 2 E^(2 x) + I E^((2 + I) x))/(1 + E^(2 x)), x]][[1]]], "current" -> "\"failed\"", "budget" -> 120,
  "followup" -> "structure-theorem independence in BuildTower over Q(i) and with flattened radicals; then certify without confirmation", "source" -> "adversarial review of PMT"|>,
<|"id" -> "G-adv-t2", "area" -> "Trager", "gap" -> "a hang over Q(i) with residue fields of degree >= 2: kDivRemNF is handed a rational (not polynomial) dividend in normalizeAtInfinity",
  "kind" -> "Integrate", "f" -> Hold[1/((x^2 + I) Sqrt[x + 1])], "current" -> "timeout", "budget" -> 60,
  "followup" -> "make the dividend polynomial (clear the denominator) before kDivRem over Q(i)", "source" -> "adversarial review of Trager"|>,
<|"id" -> "G-adv-t3", "area" -> "Trager", "gap" -> "a Gaussian binomial times out (elementary by Chebyshev)",
  "kind" -> "Integrate", "f" -> Hold[(1 + (1 + I) x^4)^(3/4)], "current" -> "timeout", "budget" -> 60,
  "followup" -> "after G-adv-t2", "source" -> "adversarial review of Trager"|>,
<|"id" -> "G-adv-t4", "area" -> "Trager", "gap" -> "a rational integrand over Q(i) is declined as not algebraic (the rational integrator is not reached)",
  "kind" -> "Expr", "f" -> Hold[Trager`TragerIntegrate[1/(x^2 + I), x][[1]]], "current" -> "\"implementation incomplete\"", "budget" -> 60,
  "followup" -> "route a radical-free Q(i) integrand to the rational integrator", "source" -> "adversarial review of Trager"|>,
<|"id" -> "G-nf-3", "area" -> "Bronstein", "gap" -> "certificates over a number field are demoted (NF:certificate): Trager's base field is fieldQ with Expand, exact for Gaussian rationals but not for radicals (1/(1 + Sqrt[2]) is not rationalised)",
  "kind" -> "Integrate", "f" -> Hold[Sqrt[Sqrt[2] + Log[x]]], "current" -> "incomplete", "budget" -> 60,
  "followup" -> "Phase 4.3b: the constant field as fieldNF[t, mu, gamma] through Trager's integrateComponent and Bronstein's k arithmetic", "source" -> "Phase 4.3"|>,
<|"id" -> "G-nf-6", "area" -> "Trager", "gap" -> "torsion over a number field: decided at degree-1 primes of the constant field since v0.6559 (orders mod 17 and 71 disagree), but the certificate is still demoted (NF:certificate)",
  "kind" -> "Integrate", "f" -> Hold[1/((x - Sqrt[2]) Sqrt[x^3 + 1])], "current" -> "incomplete", "budget" -> 60,
  "followup" -> "Phase 4.3b C3/C6/C8: exact residue/log fields over K0, then remove the demotion after an adversarial review", "source" -> "Phase 4.3"|>,

};

openGapsStructural = {
  {"ratchet", "41 throwUnimplemented sites and 12 uncertain C Risch DE paths (tools/check_bronstein_completeness.py, BASELINE_M / BASELINE_C)", "Phases 2-5"},
  {"Phase 2", "ParamRischDE and LimitedIntegrate with new logarithms (Ch. 7); rewrite dependent Log/Exp towers; rt_field_integrate completeness (m > 1 new logs, mixed resultants)", "Phase 2 remainder"},
  {"E6e", "LogPart:weights at an exponential's special places (Bronstein's Def 2.21 residues): the weight theorem not yet written for them", "Phase 3"},
  {"E5", "several ramified places above one special point (Divisors:multi-ramified)", "Phase 3"},
  {"4.1", "Trager:residue-field-QI: the tower method of residueField over Q(i)", "Phase 4.1"},
  {"policy", "PossibleZeroQ samples the real line by design (PossibleZeroQ[Log[E^x] - x] is True); verifiers of complex antiderivatives use Assumptions -> Element[x, Complexes] (was G-core-1)", "none: a documented convention"},
  {"certificates", "certificates without an external cross-check (FriCAS times out or reports incomplete): Sqrt[Log[x] + I], 1/((x - I) Sqrt[x^3 + 1]), (1 + Log[x])/(x (1 + Log[x] + Log[x]^3)^(1/3))", "an independent checker (another CAS or a proof certificate)"},
  {"Trager", "Method -> \"TragerKauers\" floods about 200 Thread::tdlen messages on Gaussian binomials (kauersPrincipal differentiates an unevaluated GroebnerBasis)", "KauersLogTerms"},
  {"tests", "integrate_cascade_trager_bronstein_tests:89 (Check of 5/((x - 3) Sqrt[x^3 - 2])) intermittently sees a firing, also when the suite runs alone, with $MessageList empty and no anonymous mth_msg_note_fired (so not the Root memo): a timing-dependent Quiet/Check frame leak, suspected around a TimeConstrained unwind; repl_pty_test and image_tests hit their alarms under load", "test infrastructure"},
  {"design", "goodPrimes searches primes without a bound: data that no prime is good for (inconsistent input) would loop; characteristic-0 prechecks guard the known cases", "a proof-backed bound or an Internal error"}
};
