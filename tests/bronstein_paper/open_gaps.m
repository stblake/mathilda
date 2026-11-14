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
<|"id" -> "G-sign-2", "area" -> "Bronstein", "gap" -> "Sqrt[Tan] + Sqrt[Cot]: the sign split's own check fails (an Internal error)",
  "kind" -> "Integrate", "f" -> Hold[Sqrt[Tan[x]] + Sqrt[Cot[x]]], "current" -> "internal", "budget" -> 60,
  "followup" -> "debug SignSplitCheckFailed over Q(i)", "source" -> "review of E4.2"|>,
<|"id" -> "G-trig-1", "area" -> "Bronstein", "gap" -> "TrigToExp writes ArcCosh with Sqrt[x - 1] Sqrt[x + 1]: two radicals (E6d over a constant-sign split)",
  "kind" -> "Integrate", "f" -> Hold[ArcCosh[x]/Sqrt[x^2 - 1]], "current" -> "incomplete", "budget" -> 60,
  "followup" -> "a front-end rewrite Sqrt[x - 1] Sqrt[x + 1] = eps Sqrt[x^2 - 1] (E6d)", "source" -> "review of E4.2"|>,
<|"id" -> "G-qi-1", "area" -> "Bronstein", "gap" -> "residue at infinity of a logarithm: the integral in k has a logarithm with non-rational coefficient (Q(i))",
  "kind" -> "Integrate", "f" -> Hold[Sqrt[I + Log[x]^2]], "current" -> "incomplete", "budget" -> 60,
  "followup" -> "brLimitedIntegrate over Q(i) with algebraic log coefficients (Phase 3 LimitedIntegrate)", "source" -> "Phase 4.1b battery"|>,
<|"id" -> "G-nf-1", "area" -> "Bronstein", "gap" -> "number-field constants (Sqrt[2], ...): Tower:algebraic-constant",
  "kind" -> "Integrate", "f" -> Hold[Sqrt[2] Sqrt[Log[x]] Log[x + 1]], "current" -> "incomplete", "budget" -> 60,
  "followup" -> "Phase 4.3", "source" -> "plan"|>,
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
<|"id" -> "G-core-1", "area" -> "Core", "gap" -> "PossibleZeroQ false zero: Log[E^x] - x is -2 Pi I for Im x in (Pi, 3 Pi]",
  "kind" -> "Expr", "f" -> Hold[PossibleZeroQ[Log[E^x] - x]], "current" -> "True", "budget" -> 60,
  "followup" -> "PossibleZeroQ: sample off the strip |Im x| < Pi for Log of exponentials", "source" -> "review of E6a"|>,
<|"id" -> "G-core-2", "area" -> "Cascade", "gap" -> "the Automatic cascade returns Indeterminate (Bronstein certifies it)",
  "kind" -> "Expr", "f" -> Hold[Quiet[Integrate[(Log[x^2] - 2 Log[x]) Exp[x^2], x]]], "current" -> "Indeterminate", "budget" -> 120,
  "followup" -> "find the stage that returns Indeterminate", "source" -> "review of E6a"|>,
<|"id" -> "G-core-3", "area" -> "Cascade", "gap" -> "an uncaught Throw[..., \"PIM\"] from the Automatic cascade",
  "kind" -> "Expr", "f" -> Hold[Head[Quiet[Integrate[Exp[x + Log[x^2]/4], x]]]], "current" -> "Throw", "budget" -> 120,
  "followup" -> "catch the PIM throw at its stage boundary", "source" -> "review of E6c"|>,
<|"id" -> "G-core-4", "area" -> "Risch (C)", "gap" -> "Log[x^2] is not its own monomial in the C transcendental integrator (the retry fixed in v0.6547 no longer answers wrongly)",
  "kind" -> "Expr", "f" -> Hold[Head[Integrate`RischTranscendental[Log[x^2]/x, x]]], "current" -> "Integrate`RischTranscendental", "budget" -> 60,
  "followup" -> "tower builder: a logarithm of a power as a monomial, or the kappa split of E6a in C", "source" -> "Phase 3 part 1"|>,
<|"id" -> "G-core-5", "area" -> "Core", "gap" -> "machine ArcSech at negative reals disagrees with its exact form (latent: no integrand reaches it)",
  "kind" -> "Expr", "f" -> Hold[Sign[Im[ArcSech[-3.7]]] == Sign[Im[N[ArcSech[-37/10]]]]], "current" -> "False", "budget" -> 60,
  "followup" -> "fix the machine ArcSech branch for x < 0", "source" -> "review of E4.2"|>,
(* ---- found by the adversarial reviews of Bronstein, Trager and PMT (2026-11-14) *)
<|"id" -> "G-adv-b2", "area" -> "Bronstein", "gap" -> "Internal brLogPartRemainderNotRational: an explicit I with residues in Q(i) over a Log tower",
  "kind" -> "Integrate", "f" -> Hold[(I/2)/(x Sqrt[Log[x]]) + 1/(2 x Sqrt[Log[x]] (1 + Log[x]))], "current" -> "internal", "budget" -> 60,
  "followup" -> "the log-part remainder over a residue field containing i (fieldSubI before the t-freeness check, LogPart.m)", "source" -> "adversarial review of Bronstein"|>,
<|"id" -> "G-adv-b3", "area" -> "Bronstein", "gap" -> "Internal AntiderivativeCheckFailed: a reducible logarithm Log[x^2] beside Log[1 + x]",
  "kind" -> "Integrate", "f" -> Hold[Log[x^2]/(1 + x) + 2 Log[1 + x]/x], "current" -> "internal", "budget" -> 60,
  "followup" -> "debug brIntegrateKappaSplit on this case", "source" -> "adversarial review of Bronstein"|>,
<|"id" -> "G-adv-b5", "area" -> "Bronstein", "gap" -> "an inverse trig times a logarithm times out (elementary: D[ArcSin[x] Log[x]])",
  "kind" -> "Integrate", "f" -> Hold[Together[D[ArcSin[x] Log[x], x]]], "current" -> "timeout", "budget" -> 60,
  "followup" -> "profile; the kappa split and the trig front end on reducible logs", "source" -> "adversarial review of Bronstein"|>,
<|"id" -> "G-adv-c1", "area" -> "Risch (C)", "gap" -> "a transcendental constant (Log[2]) beside the tower: the C decider's no is now demoted to undecided, so this correct certificate is lost",
  "kind" -> "Expr", "f" -> Hold[Head[Quiet[Risch`ElementaryIntegralQ[Log[2]/(1 + u^2) + Log[u]/(1 + u^2), u]]]], "current" -> "Risch`ElementaryIntegralQ", "budget" -> 60,
  "followup" -> "treat transcendental constants as parameters of the constant field (Phase 4.4)", "source" -> "adversarial review of Bronstein"|>,
<|"id" -> "G-adv-p1", "area" -> "ParallelMixedTower", "gap" -> "the tower builder does not prove its generators independent (Gaussian exponent combinations, an exponential beside a flattened radical, relations lost on rebuild); certificates now need independent confirmation, which drops PMT's unique ones",
  "kind" -> "Expr", "f" -> Hold[Quiet[Integrate`ParallelMixedTower[(I E^(I x) + 2 E^(2 x) + I E^((2 + I) x))/(1 + E^(2 x)), x]][[1]]], "current" -> "\"failed\"", "budget" -> 120,
  "followup" -> "structure-theorem independence in BuildTower over Q(i) and with flattened radicals; then certify without confirmation", "source" -> "adversarial review of PMT"|>,
<|"id" -> "G-adv-p4", "area" -> "ParallelMixedTower", "gap" -> "an undocumented status (needs torsion realisation, milestone iii) for elementary inputs",
  "kind" -> "Expr", "f" -> Hold[Quiet[Integrate`ParallelMixedTower[Together[D[ArcTan[Log[x]] + Sqrt[x^3 - x] (2 x - 3), x]], x]][[1]]], "current" -> "\"needs torsion realisation (milestone iii)\"", "budget" -> 120,
  "followup" -> "map it to failed (or realise the torsion) and document the protocol", "source" -> "adversarial review of PMT"|>,
<|"id" -> "G-core-6", "area" -> "Core", "gap" -> "Together does not combine fractions that carry Root objects (worked around in PMT by naming them)",
  "kind" -> "Expr", "f" -> Hold[With[{r = Root[1 + #1 + #1^4 &, 1]}, Head[Together[D[(1 + t + 2 t r + 2 t r^2 + 4 r^3 + 4 t r^3 + 2 t r^5)/(1 + t + 4 r^3 + 4 t r^3), t]]]]], "current" -> "Plus", "budget" -> 60,
  "followup" -> "Together over Q(Root): treat Root objects as algebraic generators", "source" -> "adversarial review of PMT"|>,
<|"id" -> "G-adv-t2", "area" -> "Trager", "gap" -> "a hang over Q(i) with residue fields of degree >= 2: kDivRemNF is handed a rational (not polynomial) dividend in normalizeAtInfinity",
  "kind" -> "Integrate", "f" -> Hold[1/((x^2 + I) Sqrt[x + 1])], "current" -> "timeout", "budget" -> 60,
  "followup" -> "make the dividend polynomial (clear the denominator) before kDivRem over Q(i)", "source" -> "adversarial review of Trager"|>,
<|"id" -> "G-adv-t3", "area" -> "Trager", "gap" -> "a Gaussian binomial times out (elementary by Chebyshev)",
  "kind" -> "Integrate", "f" -> Hold[(1 + (1 + I) x^4)^(3/4)], "current" -> "timeout", "budget" -> 60,
  "followup" -> "after G-adv-t2", "source" -> "adversarial review of Trager"|>,
<|"id" -> "G-adv-t4", "area" -> "Trager", "gap" -> "a rational integrand over Q(i) is declined as not algebraic (the rational integrator is not reached)",
  "kind" -> "Expr", "f" -> Hold[Trager`TragerIntegrate[1/(x^2 + I), x][[1]]], "current" -> "\"implementation incomplete\"", "budget" -> 60,
  "followup" -> "route a radical-free Q(i) integrand to the rational integrator", "source" -> "adversarial review of Trager"|>,
<|"id" -> "G-core-7", "area" -> "Core", "gap" -> "Gather with a test function is not implemented (the two-argument form stays unevaluated; PMT's generator count worked around it)",
  "kind" -> "Expr", "f" -> Hold[Head[Gather[{1, 2, 3}, Mod[#1, 2] === Mod[#2, 2] &]]], "current" -> "Gather", "budget" -> 60,
  "followup" -> "implement Gather[list, test]", "source" -> "adversarial review of PMT (fix)"|>
};

openGapsStructural = {
  {"ratchet", "41 throwUnimplemented sites and 12 uncertain C Risch DE paths (tools/check_bronstein_completeness.py, BASELINE_M / BASELINE_C)", "Phases 2-5"},
  {"Phase 2", "ParamRischDE and LimitedIntegrate with new logarithms (Ch. 7); rewrite dependent Log/Exp towers; rt_field_integrate completeness (m > 1 new logs, mixed resultants)", "Phase 2 remainder"},
  {"E6e", "LogPart:weights at an exponential's special places (Bronstein's Def 2.21 residues): the weight theorem not yet written for them", "Phase 3"},
  {"E5", "several ramified places above one special point (Divisors:multi-ramified)", "Phase 3"},
  {"4.1", "Trager:residue-field-QI: the tower method of residueField over Q(i)", "Phase 4.1"},
  {"certificates", "certificates without an external cross-check (FriCAS times out or reports incomplete): Sqrt[Log[x] + I], 1/((x - I) Sqrt[x^3 + 1]), (1 + Log[x])/(x (1 + Log[x] + Log[x]^3)^(1/3))", "an independent checker (another CAS or a proof certificate)"},
  {"Trager", "Method -> \"TragerKauers\" floods about 200 Thread::tdlen messages on Gaussian binomials (kauersPrincipal differentiates an unevaluated GroebnerBasis)", "KauersLogTerms"},
  {"tests", "integrate_cascade_trager_bronstein_tests:89 emits a message only under full parallel ctest load; repl_pty_test and image_tests hit their alarms under load", "test infrastructure"},
  {"design", "goodPrimes searches primes without a bound: data that no prime is good for (inconsistent input) would loop; characteristic-0 prechecks guard the known cases", "a proof-backed bound or an Internal error"}
};
