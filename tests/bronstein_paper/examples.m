(* ::Package:: *)
(* ========================================================================== *)
(* tests/bronstein_paper/examples.m -- the gap log of the Bronstein completion *)
(* campaign, as checked examples (the corpus behind the JSC paper's tables).  *)
(*                                                                            *)
(* One record per example:                                                    *)
(*   "id"      stable label, <section>-<n>                                    *)
(*   "gap"     the ratchet site id (or campaign item) the example closed      *)
(*   "theory"  the section of docs/notes/bronstein_extensions.tex             *)
(*   "version" the release that closed it                                     *)
(*   "kind"    "Integrate"  BronsteinIntegrate[f, x]                          *)
(*             "EIQ"        Risch`ElementaryIntegralQ[f, x]                   *)
(*             "RDE"        Risch`RischDECertified[f, g, x] (f = {f, g})      *)
(*             "Scope"      BronsteinScope[f, x]                              *)
(*   "f"       the input, held                                                *)
(*   "expect"  "elementary" (a verified antiderivative), "nonelementary" (a  *)
(*             certificate), "incomplete", True/False (EIQ), a solution or     *)
(*             "NoSolution" (RDE), "InScope"                                  *)
(*   "before"  what the system did before the version (for the paper)         *)
(*   "note"    why the example is non-trivial                                 *)
(* Every example is checked by tests/bronstein_paper/run.m (ctest            *)
(* bronstein_paper_examples); tools/bronstein_paper_table.py renders them     *)
(* into docs/notes/bronstein_examples_table.tex for the theory note.          *)
(* Add the examples of a gap in the same commit as its fix.                   *)
(* ========================================================================== *)

bronsteinPaperExamples = {

(* ---- Phase 0: soundness (v0.6545) --------------------------------------- *)
<|"id" -> "P0-1", "gap" -> "Phase 0: unrun Risch DE read as a proof", "theory" -> "--", "version" -> "v0.6545",
  "kind" -> "EIQ", "f" -> Hold[Cos[E^x]], "expect" -> False,
  "before" -> "False, by accident: the base solver could not form the equation",
  "note" -> "the lower monomial E^x reached the base-field solver in surface form; the empty result was read as a theorem"|>,
<|"id" -> "P0-2", "gap" -> "Phase 0: unrun Risch DE read as a proof", "theory" -> "--", "version" -> "v0.6545",
  "kind" -> "EIQ", "f" -> Hold[E^x Cos[E^x]], "expect" -> True,
  "before" -> "True", "note" -> "the same tower, elementary: Sin[E^x]"|>,

(* ---- Phase 1: implementation incomplete (v0.6545) ------------------------ *)
<|"id" -> "P1-1", "gap" -> "Phase 1: Integrate::incomplete", "theory" -> "--", "version" -> "v0.6545",
  "kind" -> "Scope", "f" -> Hold[Sqrt[Log[x] + Sqrt[Log[x]]]], "expect" -> "incomplete",
  "before" -> "unevaluated, silent", "note" -> "nested radicals: outside one simple radical extension, reported as such"|>,
<|"id" -> "P1-2", "gap" -> "Phase 1: Kummer cyclicity", "theory" -> "--", "version" -> "v0.6545",
  "kind" -> "Scope", "f" -> Hold[Sqrt[x] + Sqrt[x + 1]], "expect" -> "incomplete",
  "before" -> "unevaluated, silent", "note" -> "Z/2 x Z/2: a genuine compositum"|>,

(* ---- Phase 2: the tower Risch DE (v0.6546) ------------------------------- *)
<|"id" -> "P2-1", "gap" -> "RdeBoundDegreePrim (deg b = deg a)", "theory" -> "A2", "version" -> "v0.6546",
  "kind" -> "RDE", "f" -> Hold[{1/x, 2 Log[x]/x^2}], "expect" -> Hold[Log[x]^2/x],
  "before" -> "undecided", "note" -> "the solution has degree dc - da + 1 in Log[x]: the bound an off-by-one first got wrong"|>,

(* ---- E6b: a cyclic compositum as one radical (v0.6547) ------------------- *)
<|"id" -> "E6b-1", "gap" -> "Tower:cyclic-compositum", "theory" -> "E6b", "version" -> "v0.6547",
  "kind" -> "Integrate", "f" -> Hold[(Sqrt[x] + (x + 1)^(1/3))/(x (x + 1))], "expect" -> "elementary",
  "before" -> "implementation incomplete", "note" -> "y^6 = x^3 (x+1)^2, y = Sqrt[x] (x+1)^(1/3); answer in the user's radicals"|>,
<|"id" -> "E6b-2", "gap" -> "Tower:cyclic-compositum", "theory" -> "E6b", "version" -> "v0.6547",
  "kind" -> "Integrate", "f" -> Hold[Sqrt[x] (x + 1)^(1/3)], "expect" -> "nonelementary",
  "before" -> "implementation incomplete", "note" -> "Chebyshev's binomial case, certified by the Risch-Trager-Bronstein machinery on y^6 = x^3 (x+1)^2"|>,
<|"id" -> "E6b-3", "gap" -> "Tower:cyclic-compositum", "theory" -> "E6b", "version" -> "v0.6547",
  "kind" -> "Integrate", "f" -> Hold[D[Sqrt[x - 1] (x + 1)^(1/3) E^x, x]], "expect" -> "elementary",
  "before" -> "implementation incomplete", "note" -> "the compositum at the bottom of a transcendental tower (Shape B, E^x on top)"|>,
<|"id" -> "E6b-5", "gap" -> "Tower:cyclic-compositum", "theory" -> "E6b", "version" -> "v0.6547",
  "kind" -> "Integrate", "f" -> Hold[Sqrt[x] (x + 1)^(1/3) Log[x]], "expect" -> "nonelementary",
  "before" -> "implementation incomplete", "note" -> "Log[x] over the compositum: certified by the coefficient Risch DE over the base"|>,
<|"id" -> "E6b-6", "gap" -> "Tower:cyclic-compositum", "theory" -> "E6b", "version" -> "v0.6547",
  "kind" -> "Integrate", "f" -> Hold[1/(x Sqrt[Log[x]] (Log[x] + 1)^(1/3))], "expect" -> "nonelementary",
  "before" -> "implementation incomplete", "note" -> "Chebyshev's case in u = Log[x] (u^(-1/2) (1 + u)^(-1/3)), decided in the tower"|>,
<|"id" -> "E6b-4", "gap" -> "Tower:cyclic-compositum", "theory" -> "E6b", "version" -> "v0.6547",
  "kind" -> "Integrate", "f" -> Hold[Sqrt[Log[x]] (Log[x] + 1)^(1/3)/x], "expect" -> "nonelementary",
  "before" -> "implementation incomplete", "note" -> "the compositum over the transcendental tower Q(x)(Log[x])"|>,

(* ---- E6d: one radical up to a locally constant sign (v0.6549) ----------- *)
<|"id" -> "E6d-1", "gap" -> "Tower:cyclic-compositum (sign)", "theory" -> "E6d", "version" -> "v0.6549",
  "kind" -> "Integrate", "f" -> Hold[Sqrt[x^2 (x + 1)]/x^2 + 1/Sqrt[x + 1]], "expect" -> "elementary",
  "before" -> "implementation incomplete", "note" -> "Sqrt[x^2 (x+1)] = eps x Sqrt[x+1], eps = +-1 locally constant; the answer is a sum over sign regions"|>,
<|"id" -> "E6d-2", "gap" -> "Tower:cyclic-compositum (sign)", "theory" -> "E6d", "version" -> "v0.6549",
  "kind" -> "Integrate", "f" -> Hold[Sqrt[x^2 (x^3 + 1)] - x Sqrt[x^3 + 1]], "expect" -> "nonelementary",
  "before" -> "implementation incomplete", "note" -> "identically 0 where eps = 1, the elliptic -2 x Sqrt[x^3+1] where eps = -1"|>,
<|"id" -> "E6d-3", "gap" -> "Tower:cyclic-compositum (sign)", "theory" -> "E6d", "version" -> "v0.6549",
  "kind" -> "Integrate", "f" -> Hold[Sqrt[x^2 (x^3 + 1)] Sqrt[x^3 + 1]/x], "expect" -> "elementary",
  "before" -> "implementation incomplete", "note" -> "the product is eps (x^3 + 1): rational on each region"|>,
<|"id" -> "E6d-4", "gap" -> "Tower:cyclic-compositum (sign)", "theory" -> "E6d", "version" -> "v0.6549",
  "kind" -> "Integrate", "f" -> Hold[Sqrt[Log[x]^2 (Log[x] + 1)] + Sqrt[Log[x] + 1]/x], "expect" -> "nonelementary",
  "before" -> "implementation incomplete", "note" -> "over the tower: Log[x] Sqrt[Log[x] + 1] on the region eps = 1"|>,

<|"id" -> "E6d-5", "gap" -> "Tower:cyclic-compositum (sign)", "theory" -> "E6d", "version" -> "v0.6550",
  "kind" -> "Integrate", "f" -> Hold[(1 + E^(x^2))/2 + (1 - E^(x^2))/2 Sqrt[x + 1] Sqrt[1/(x + 1)]], "expect" -> "incomplete",
  "before" -> "not elementary (v0.6549: false, a sign seen only on a branch cut)",
  "note" -> "equals 1 on every open set; the sign -1 occurs only on the cut x < -1 and must never certify"|>,

(* ---- E6c: log-linear exponents (v0.6547) --------------------------------- *)
<|"id" -> "E6c-1", "gap" -> "Tower:exp-of-log-radical", "theory" -> "E6c", "version" -> "v0.6547",
  "kind" -> "Integrate", "f" -> Hold[E^(x + Log[x]/2)/x], "expect" -> "nonelementary",
  "before" -> "undecided (a new hyperexponential monomial whose Risch DE needs a radical)", "note" -> "E^x/Sqrt[x] in disguise"|>,
<|"id" -> "E6c-4", "gap" -> "Tower:exp-of-log-radical", "theory" -> "E6c", "version" -> "v0.6547",
  "kind" -> "Integrate", "f" -> Hold[E^(x^2 + Log[Log[x]]/2)], "expect" -> "nonelementary",
  "before" -> "undecided", "note" -> "E^(x^2) Sqrt[Log[x]]: Shape A with the exponential below the radicand's monomial"|>,
<|"id" -> "E6c-2", "gap" -> "Tower:exp-of-log-radical", "theory" -> "E6c", "version" -> "v0.6547",
  "kind" -> "Integrate", "f" -> Hold[D[E^(x^2 + Log[x + 1]/3), x]], "expect" -> "elementary",
  "before" -> "undecided", "note" -> "E^(x^2) (x+1)^(1/3): a radical at the bottom under an exponential"|>,
<|"id" -> "E6c-3", "gap" -> "Tower:exp-of-log-radical", "theory" -> "E6c", "version" -> "v0.6547",
  "kind" -> "Integrate", "f" -> Hold[D[E^(Log[Log[x]]/2 + x), x]], "expect" -> "elementary",
  "before" -> "implementation incomplete (radical in the middle of the tower)", "note" -> "Log[Log[x]] is rewritten away and must be pruned"|>,

(* ---- E6a: reducible logarithms (v0.6547) --------------------------------- *)
<|"id" -> "E6a-1", "gap" -> "Tower:reducible-log", "theory" -> "E6a", "version" -> "v0.6547",
  "kind" -> "Integrate", "f" -> Hold[Log[x^2] Sqrt[Log[x]]/x], "expect" -> "elementary",
  "before" -> "implementation incomplete", "note" -> "answer carries the branch difference Log[x^2] - 2 Log[x] explicitly"|>,
<|"id" -> "E6a-2", "gap" -> "Tower:reducible-log", "theory" -> "E6a", "version" -> "v0.6547",
  "kind" -> "Integrate", "f" -> Hold[(Log[x^2] - 2 Log[x]) Sqrt[Log[x]]], "expect" -> "nonelementary",
  "before" -> "implementation incomplete", "note" -> "zero on Re x > 0; certified by Hermite-Lindemann + specialisation (Thm E6a(3))"|>,
<|"id" -> "E6a-3", "gap" -> "Tower:reducible-log", "theory" -> "E6a", "version" -> "v0.6547",
  "kind" -> "Integrate", "f" -> Hold[Log[x (x + 1)] Sqrt[Log[x] + Log[x + 1]] (1/x + 1/(x + 1))], "expect" -> "elementary",
  "before" -> "implementation incomplete", "note" -> "a product logarithm over two log monomials"|>,
<|"id" -> "E6a-4", "gap" -> "Tower:reducible-log", "theory" -> "E6a", "version" -> "v0.6547",
  "kind" -> "Integrate", "f" -> Hold[(Log[E^x] - x) Sqrt[Log[x]] + 1/x], "expect" -> "nonelementary",
  "before" -> "Log[x] from the Automatic cascade (wrong: PossibleZeroQ[Log[E^x] - x] is True)",
  "note" -> "the branch difference vanishes only on |Im x| < Pi; certified at a sample with large imaginary part"|>,

(* ---- E3: special places, positive genus (v0.6548) ------------------------ *)
<|"id" -> "E3-1", "gap" -> "LogPart:special-genus", "theory" -> "E3", "version" -> "v0.6548",
  "kind" -> "Integrate", "f" -> Hold[1/(x (1 + Log[x]^3)^(1/3))], "expect" -> "elementary",
  "before" -> "implementation incomplete", "note" -> "genus 1 (Fermat cubic over Q(x)(Log[x])); residues at theta = infinity"|>,
<|"id" -> "E3-2", "gap" -> "LogPart:special-genus", "theory" -> "E3", "version" -> "v0.6548",
  "kind" -> "Integrate", "f" -> Hold[E^x/(1 + E^(3 x))^(1/3)], "expect" -> "elementary",
  "before" -> "implementation incomplete", "note" -> "the exponential analogue: special places over theta = 0 and infinity"|>,
<|"id" -> "E3-3", "gap" -> "LogPart:special-genus", "theory" -> "E3", "version" -> "v0.6548",
  "kind" -> "Integrate", "f" -> Hold[1/(x (1 + Log[x]^4)^(1/4))], "expect" -> "elementary",
  "before" -> "implementation incomplete", "note" -> "m = 4, genus 3"|>,
<|"id" -> "E3-4", "gap" -> "LogPart:special-genus", "theory" -> "E3", "version" -> "v0.6548",
  "kind" -> "Integrate", "f" -> Hold[D[(1 + Log[x] + Log[x]^3)^(1/3) + Log[Log[x] + (1 + Log[x] + Log[x]^3)^(1/3)], x]],
  "expect" -> "elementary", "before" -> "implementation incomplete",
  "note" -> "algebraic part plus a log part with finite and special residues on y^3 = L^3 + L + 1"|>,
<|"id" -> "E3-5", "gap" -> "LogPart:special-genus", "theory" -> "E3", "version" -> "v0.6548",
  "kind" -> "Integrate", "f" -> Hold[(1 + Log[x])/(x (1 + Log[x] + Log[x]^3)^(1/3))], "expect" -> "nonelementary",
  "before" -> "implementation incomplete", "note" -> "the residue divisor has infinite order (orders mod 13 and 37 disagree)"|>,
<|"id" -> "E3-6", "gap" -> "LogPart:special-genus", "theory" -> "E3", "version" -> "v0.6548",
  "kind" -> "Integrate", "f" -> Hold[D[x + Log[-E^x + (1 + E^(3 x))^(1/3)], x]], "expect" -> "elementary",
  "before" -> "implementation incomplete", "note" -> "exponential tower, genus 1"|>,
(* ---- E4.1: the constant field Q(i) (v0.6551) ------------------------------ *)
<|"id" -> "E41-1", "gap" -> "Tower:complex, Trager:constants", "theory" -> "E4.1", "version" -> "v0.6551",
  "kind" -> "Integrate", "f" -> Hold[Sqrt[x + I]/x], "expect" -> "elementary",
  "before" -> "implementation incomplete", "note" -> "a curve over Q(i)(x); FriCAS agrees"|>,
<|"id" -> "E41-2", "gap" -> "Tower:complex", "theory" -> "E4.1", "version" -> "v0.6551",
  "kind" -> "Integrate", "f" -> Hold[D[Log[E^(I x) + Sqrt[1 + E^(2 I x)]], x]], "expect" -> "elementary",
  "before" -> "implementation incomplete", "note" -> "Shape A over the exponential E^(I x): an inverse-hyperbolic of E^(I x)"|>,
<|"id" -> "E41-3", "gap" -> "Tower:complex", "theory" -> "E4.1", "version" -> "v0.6551",
  "kind" -> "Integrate", "f" -> Hold[1/Sqrt[E^(2 I x) + 1]], "expect" -> "elementary",
  "before" -> "implementation incomplete", "note" -> "residues at the special places of E^(I x) over Q(i)"|>,
<|"id" -> "E41-4", "gap" -> "Tower:complex", "theory" -> "E4.1", "version" -> "v0.6551",
  "kind" -> "Integrate", "f" -> Hold[(1 + I) Sqrt[x]/(x^2 + I)], "expect" -> "elementary",
  "before" -> "implementation incomplete", "note" -> "residue field Q(i)(zeta_8): i adjoined and written in the primitive element"|>,
<|"id" -> "E41-5", "gap" -> "Tower:complex", "theory" -> "E4.1", "version" -> "v0.6551",
  "kind" -> "Integrate", "f" -> Hold[Sqrt[Log[x] + I]], "expect" -> "nonelementary",
  "before" -> "implementation incomplete", "note" -> "FriCAS: implementation incomplete (constant residues); Risch DE b' + b (1 + 1/(2 (u + I))) = 1 has no rational solution"|>
};
