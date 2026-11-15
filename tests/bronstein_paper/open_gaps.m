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
<|"id" -> "G-nf-10", "area" -> "Bronstein", "gap" -> "random elementary integrands over Q(constants) (seed 20261114, 40 cases): 17 verified, 15 time out, 6 internal (brResultantRZero, brComponentIdentity, brHermiteIterationBound before brCancelK), 2 undecided -- no wrong answer and no false certificate",
  "kind" -> "Expr", "f" -> Hold[Length[{}]], "current" -> "0", "budget" -> 10,
  "followup" -> "rerun the battery (scratch c6r.m pattern) after each 4.3b step; promote it to a nightly stress ctest", "source" -> "Phase 4.3b C6"|>,

(* ---- Bronstein scope ----------------------------------------------------- *)
<|"id" -> "G-sign-1", "area" -> "Bronstein", "gap" -> "radicals equal up to a sign seen only on a branch cut: undecided (E6d)",
  "kind" -> "Integrate", "f" -> Hold[(1 + E^(x^2))/2 + (1 - E^(x^2))/2 Sqrt[x + 1] Sqrt[1/(x + 1)]], "current" -> "incomplete", "budget" -> 60,
  "followup" -> "prove eps = 1 on every open set (no open region with eps = -1: the connected components of the complement of the cuts, e.g. by a cylindrical decomposition in (Re x, Im x)) and integrate f_+ alone; brSignMerge (v0.6567) already reduces it to the single sign Sqrt[x + 1] Sqrt[1/(x + 1)]", "source" -> "review of E6d"|>,
<|"id" -> "G-mid-1", "area" -> "Bronstein", "gap" -> "a radical in the middle of the tower",
  "kind" -> "Integrate", "f" -> Hold[Log[1 + Sqrt[Log[x]]]/x], "current" -> "incomplete", "budget" -> 60,
  "followup" -> "Phase 5 (E1, E2)", "source" -> "plan"|>,

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
(* ---- found by the adversarial reviews of Bronstein, Trager and PMT (2026-11-14) *)
<|"id" -> "G-adv-c1", "area" -> "Risch (C)", "gap" -> "a transcendental constant (Log[2]) beside the tower: the C decider's no is now demoted to undecided, so this correct certificate is lost",
  "kind" -> "Expr", "f" -> Hold[Head[Quiet[Risch`ElementaryIntegralQ[Log[2]/(1 + u^2) + Log[u]/(1 + u^2), u]]]], "current" -> "Risch`ElementaryIntegralQ", "budget" -> 60,
  "followup" -> "treat transcendental constants as parameters of the constant field (Phase 4.4)", "source" -> "adversarial review of Bronstein"|>,
<|"id" -> "G-adv-p1", "area" -> "ParallelMixedTower", "gap" -> "the tower builder does not prove its generators independent (Gaussian exponent combinations, an exponential beside a flattened radical, relations lost on rebuild); certificates now need independent confirmation, which drops PMT's unique ones",
  "kind" -> "Expr", "f" -> Hold[Quiet[Integrate`ParallelMixedTower[(I E^(I x) + 2 E^(2 x) + I E^((2 + I) x))/(1 + E^(2 x)), x]][[1]]], "current" -> "\"failed\"", "budget" -> 120,
  "followup" -> "structure-theorem independence in BuildTower over Q(i) and with flattened radicals; then certify without confirmation", "source" -> "adversarial review of PMT"|>
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
  {"design", "goodPrimes searches primes without a bound: data that no prime is good for (inconsistent input) would loop; characteristic-0 prechecks guard the known cases", "a proof-backed bound or an Internal error"},
  {"Core", "Together over Gaussian coefficients: fixed in v0.6571 (two stages: combine with I free, then one reduction over Q(I); 29.6 s -> 0.66 s, random sums 3x); what remains is the reduction over Q(I) itself, 10x the symbol route of brComponentIdentityQ (0.07 s), and a unit (-1, +-I) left where the input put it, so Together over Q(I) is canonical only up to a unit", "a unit normalisation for Gaussian Together (a convention: leading coefficient of the denominator positive real?)"},
  {"Phase 4.5", "the multi-modular kernels (ff_modular.c) decline when an entry involves a parameter (the direct ff_fnx elimination runs): evaluation and interpolation in the parameters (G-perf-4 closed v0.6573 through the absolute field)", "Phase 4.5 remainder"}
};
