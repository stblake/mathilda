(* ====================================================================== *)
(* Refute.m -- counterexample search, refutations and Refute (loaded by    *)
(* Lean.m after the backends).                                              *)
(*                                                                          *)
(* A refutation proves the NEGATION of a claim's certified statement S     *)
(* (the obligation's "Statement", what a backend would prove) from a        *)
(* concrete witness w: Lean's kernel checks                                 *)
(*   theorem : not (S) := fun h => absurd (h w...) (by decide +kernel)      *)
(* Faithfulness (plan section 3.6): w must lie in the quantifier domain of  *)
(* S -- integers for Int, rationals for Rat and for the generic char-0 /    *)
(* ordered-field statements (Q is such a field, so a rational witness       *)
(* refutes "for every field ...", and Q is inside C and R, so the user's    *)
(* claim is refuted too) -- every hypothesis of S must hold at w, its       *)
(* conclusion must fail there, and every denominator of the claim must be  *)
(* nonzero at w (a point where the rational claim is undefined is not a     *)
(* counterexample, even when it refutes the cleared statement).  All of it  *)
(* is checked by Mathilda's exact arithmetic before a refutation is built;  *)
(* FindInstance's domain argument is never trusted (it can return x -> -2 I *)
(* for Rationals).  An algebraic witness (-2 I, -Sqrt[2]) has no kernel     *)
(* proof in core Lean: such a refutation is Exact-only.                     *)
(*                                                                          *)
(* Certify: when a kernel backend REJECTS an item, refuteAfterRejection     *)
(* searches a counterexample (bounded, "Counterexample" -> Automatic = 2 s), *)
(* runs every refutation of the call in one Lean file and attaches          *)
(* "Counterexample", "Refuted", "RefutationTheorem" and "RefutationSource". *)
(* Refute[claim] states the claim, searches ("SearchTime", 5 s) and gives a *)
(* ProofObject of Kind "Refutation".                                        *)
(*                                                                          *)
(* The messages of this file are defined here (Lean.m's line budget); the   *)
(* file opens the package itself (every symbol is in Lean`Private`; Refute *)
(* is exported by Lean.m).                                                  *)
(* ====================================================================== *)

BeginPackage["Lean`"];
Begin["`Private`"];

Certify::refuted = "`1` is false: `2` is a counterexample (see obj[\"RefutationTheorem\"])";
Refute::none = "no counterexample to `1` was found within the search bounds; this is not a proof of the claim";
Refute::nostmt = "`1` is not a claim Refute can state and refute: an identity a == b, an implication of equations, a ForAll or inequality, Sum[t, {k, a, n}] == C or Integrate[f, x] == F (primality, Tier 1, optimality and computation certificates are not refutable)";
Refute::optx = "unknown option `1`; the options of Refute are `2`";
Refute::optvalue = "the option name `1` is a symbol with the value `2`, which would replace the name; Clear[`1`] or write the name as the string \"`3`\"";
Refute::optdefault = "the default Assumptions -> `1` was set with SetOptions[Refute, ...], which evaluated it before Refute could check it; pass Assumptions in the call instead, or reset it with SetOptions[Refute, Assumptions -> {}]";
Refute::engine = "unknown Engine `1`; Refute uses \"Lean\" or \"Exact\"";

(* Test-only hook: when not None, the search returns this rule list
   unchecked (witnessOK must then refuse a non-witness). *)
$witnessOverride = None;

(* ---------------------------------------------------------------------- *)
(* The refutation problem of an IR                                          *)
(* ---------------------------------------------------------------------- *)

(* <|"Form" -> "Forall" | "Closed" | "Induction" | "Antiderivative",
     "Vars", "Body" -> Hold[proposition in Vars], "Dens" -> {Hold[d]...},
     "DomainForm", "Antecedents", "Statement" -> held S, "Shift"|>, or
   $Failed when the obligation is not refutable.  Body is false at w exactly
   when w is a counterexample to S (a ForAll[vs, h, s] body is Implies[h, s];
   an implication chain is false when every antecedent holds and the goal
   fails). *)
refProblem[ir_] := refProblemType[ir["Obligation"]["Type"], ir];

$refutableTypes = {"Identity", "PolyIdentities", "IdealImplication", "WeightedSOS",
  "RootsComplete", "ClosedFact"};

refProblemType[ty_, ir_] /; MemberQ[$refutableTypes, ty] &&
    !(ty === "PolyIdentities" && ir["Kind"] === "Antiderivative") :=
  Module[{ob = ir["Obligation"], vars, body, k, dens},
  {vars, body, k} = stmtParts[ob["Statement"]];
  dens = Join[Cases[ir["Hypotheses"], HoldPattern[HoldForm[Unequal[d_, 0]]] :> Hold[d]],
    Hold /@ denominatorBases[Lookup[ob, "Pairs", {}]]];
  <|"Form" -> If[vars === {}, "Closed", "Forall"], "Vars" -> vars, "Body" -> body,
    "Dens" -> DeleteDuplicates[dens], "DomainForm" -> ir["DomainForm"], "Antecedents" -> k,
    "Statement" -> ob["Statement"], "ClosedFact" -> (ty === "ClosedFact")|>];

(* A closed-form sum: the claim at the integers n >= a (S over m : Nat is
   the claim at n = m + a). *)
refProblemType["Induction", ir_] := Module[{ob = ir["Obligation"], vars, body},
  {vars, body} = Most[stmtParts[ob["Statement"]]];
  <|"Form" -> "Induction", "Vars" -> vars, "Body" -> body, "Dens" -> {},
    "DomainForm" -> "Int", "Antecedents" -> 0, "Statement" -> ob["Statement"],
    "Shift" -> ob["LowerLimit"]|>];

(* F is an antiderivative of f (Tier 0 claim or computation): D[F, x] - f at
   rational points where both sides are rational numbers; Exact-only. *)
refProblemType["PolyIdentities", ir_] /; ir["Kind"] === "Antiderivative" &&
    ir["Certificate"]["Form"] === "ClearedDerivative" :=
  With[{c = ir["Certificate"]}, antiderivProblem[ir["Result"], c["Integrand"][[1]]/c["Integrand"][[2]],
    c["Variable"]]];
refProblemType["Antiderivative", ir_] := With[{ob = ir["Obligation"]},
  antiderivProblem[ob["F"], ob["f"], ob["Variable"]]];
refProblemType[_, _] := $Failed;

antiderivProblem[F_, f_, x_] := With[{dF = Quiet[D[F, x]], ff = f, xx = x},
  <|"Form" -> "Antiderivative", "Vars" -> {x}, "Body" -> Hold[antiderivSame[dF, ff]],
    "Dens" -> {}, "DomainForm" -> "Rat", "Antecedents" -> 0,
    "Statement" -> HoldForm[ForAll[xx, D[F, xx] == ff]]|>];
(* Both sides at a point: rational numbers compared exactly, else undecided
   (a pole, a surd or a transcendental value is never used). *)
antiderivSame[a_, b_] := If[ratNumberQ[a] && ratNumberQ[b], a === b, refUndecided];

(* {vars, Hold[body], antecedents} of a held statement. *)
stmtParts[HoldForm[ForAll[vs_, h_, s_]]] := {Flatten[{vs}], Hold[Implies[h, s]], 1};
stmtParts[HoldForm[ForAll[vs_, s_]]] := {Flatten[{vs}], Hold[s], implDepth[Hold[s]]};
stmtParts[HoldForm[s_]] := {{}, Hold[s], 0};
implDepth[Hold[Implies[_, r_]]] := 1 + implDepth[Hold[r]];
implDepth[_] := 0;

(* ---------------------------------------------------------------------- *)
(* Exact truth at a witness                                                 *)
(* ---------------------------------------------------------------------- *)

(* Relations decided by exact arithmetic: a difference that normalises to an
   exact number (rational or complex rational) is zero or not; an algebraic
   difference is zero when RootReduce gives 0 and nonzero only when
   RootReduce does not give 0 AND its 50-digit value is away from 0; an
   order relation needs a rational difference.  Anything else is
   refUndecided, which keeps the proposition from evaluating to False. *)
exactNumberQ[r_] := ratNumberQ[r] || MatchQ[r, Complex[_?ratNumberQ, _?ratNumberQ]];
refNorm[d_] := With[{r = Quiet[Together[Expand[d]]]}, If[exactNumberQ[r], r, Quiet[RootReduce[r]]]];
refZero[d_] := With[{r = refNorm[d]}, Which[
  r === 0, True,
  exactNumberQ[r], False,
  algebraicQ[r] && With[{v = Quiet[N[r, 50]]}, NumberQ[v] && Abs[v] > 10^-30], False,
  True, refUndecided]];
refSign[d_] := With[{r = refNorm[d]}, If[ratNumberQ[r], r, refUndecided]];
refEqual[a_, b_] := refZero[a - b];
refUnequal[a_, b_] := With[{z = refZero[a - b]}, If[BooleanQ[z], !z, refUndecided]];
refOrder[rel_, a_, b_] := With[{s = refSign[a - b]}, If[ratNumberQ[s], rel[s, 0], refUndecided]];
$refRelations = {Equal -> refEqual, Unequal -> refUnequal,
  GreaterEqual -> (refOrder[GreaterEqual, #1, #2] &), LessEqual -> (refOrder[LessEqual, #1, #2] &),
  Greater -> (refOrder[Greater, #1, #2] &), Less -> (refOrder[Less, #1, #2] &)};

(* An algebraic number built from rationals, I, + * and rational powers. *)
algebraicQ[v_?ratNumberQ] := True;
algebraicQ[Complex[a_, b_]] := ratNumberQ[a] && ratNumberQ[b];
algebraicQ[Power[b_, e_?ratNumberQ]] := algebraicQ[b];
algebraicQ[(Plus | Times)[a___]] := AllTrue[{a}, algebraicQ];
algebraicQ[_] := False;

(* The witness's values lie in the quantifier domain of S (kernel: integers
   for Int, rationals otherwise; Exact-only: algebraic numbers for the char-0
   form, real algebraic numbers for the ordered form). *)
valueOK["Int", v_] := IntegerQ[v];
valueOK["Rat", v_] := ratNumberQ[v];
valueOK["CharZeroField", v_] := algebraicQ[v];
valueOK["OrderedField", v_] := algebraicQ[v] && FreeQ[v, Complex];
valueOK[_, _] := False;
kernelWitnessQ[prob_, w_] := AllTrue[w[[All, 2]], If[prob["DomainForm"] === "Int", IntegerQ,
  ratNumberQ]] && prob["Form"] =!= "Antiderivative";

(* The full check of a witness (every refutation is built only after it):
   the variables of S in order, domain membership, every denominator of the
   claim nonzero, the proposition False -- all by exact arithmetic. *)
witnessOK[prob_, w_] := ListQ[w] && MatchQ[w, {___Rule}] &&
  If[prob["Form"] === "Closed", w === {} && assumeOK[prob, w] && closedTruth[prob] === False,
    w[[All, 1]] === prob["Vars"] &&
    AllTrue[w[[All, 2]], valueOK[prob["DomainForm"], #] &] &&
    AllTrue[prob["Dens"], refZero[ReleaseHold[# /. w]] === False &] &&
    assumeOK[prob, w] &&
    Quiet[ReleaseHold[prob["Body"] /. w /. $refRelations]] === False];

(* The Domain and Assumptions of the call, for a claim whose statement S
   omits them (only identity claims state them; S is then stronger than the
   claim, so a proof of S proves the claim, but a counterexample to S
   refutes the claim only when it satisfies them).  The problem gains
   "Assume" -> Hold[condition]; every witness must satisfy it exactly. *)
withAssume[prob_, ir_, o_] := Module[{c = {}},
  If[MemberQ[{"Identity", "Factorization"}, ir["Kind"]], Return[prob, Module]];
  If[MemberQ[{Integers, Rationals, Reals, Complexes}, o["Domain"]],
    c = Join[c, Element[#, o["Domain"]] & /@ prob["Vars"]]];
  If[!MatchQ[o["Assumptions"], {} | True | Automatic],
    c = Join[c, Flatten[{o["Assumptions"]} /. And -> List]]];
  If[c === {}, prob, Join[prob, <|"Assume" -> Hold @@ {And @@ c},
    "AssumeText" -> ToString[And @@ c, InputForm]|>]]];
assumeOK[prob_, w_] := !KeyExistsQ[prob, "Assume"] ||
  TrueQ[Quiet[ReleaseHold[prob["Assume"] /. w /. $refRelations]]];

closedTruth[prob_] := If[TrueQ[prob["ClosedFact"]], Quiet[exClosed[prob["Statement"]]],
  Quiet[ReleaseHold[prob["Body"] /. $refRelations]]];

(* The quick test of a rational grid point: denominators, then the
   proposition by Mathilda's (exact rational) evaluation. *)
fastWitnessQ[prob_, w_] := Quiet[If[AllTrue[prob["Dens"], ReleaseHold[# /. w] =!= 0 &] &&
    assumeOK[prob, w], ReleaseHold[prob["Body"] /. w] === False, False]];

(* ---------------------------------------------------------------------- *)
(* The search (deterministic, no randomness)                                *)
(* ---------------------------------------------------------------------- *)

(* Grid values by increasing height. *)
$intGrid = {0, -1, 1, -2, 2, -3, 3, -4, 4, -5, 5};
$ratGrid = {0, -1, 1, -2, 2, -3, 3, -1/2, 1/2, -3/2, 3/2, -1/3, 1/3, -2/3, 2/3};
$gridMax = 4000;

(* The first grid point (tuples by increasing largest height) that is a
   witness, or None. *)
gridSearch[prob_, vals_List] := Module[{n = Length[prob["Vars"]], count = 0, w},
  Do[
    Do[If[MemberQ[t, m],
        count++;
        If[count > $gridMax, Return[None, Module]];
        w = Thread[prob["Vars"] -> vals[[t + 1]]];
        If[fastWitnessQ[prob, w], Return[w, Module]]],
      {t, Tuples[Range[0, m], n]}],
    {m, 0, Length[vals] - 1}];
  None];

(* The negated body for FindInstance: Implies[a, b] -> a && not b, and the
   negation pushed into the relations. *)
SetAttributes[negHeld, HoldAll];
negHeld[Implies[a_, b_]] := a && negHeld[b];
negHeld[Equal[a_, b_]] := a != b;
negHeld[Unequal[a_, b_]] := a == b;
negHeld[GreaterEqual[a_, b_]] := a < b;
negHeld[LessEqual[a_, b_]] := a > b;
negHeld[Greater[a_, b_]] := a <= b;
negHeld[Less[a_, b_]] := a >= b;
negHeld[e_] := !e;

(* FindInstance of the negated body (denominators nonzero) in dom, as a rule
   list in the variables' order; None when it gives nothing. *)
findWitness[prob_, dom_, t_] := Module[{cond, r},
  If[t <= 0, Return[None, Module]];
  cond = (negHeld @@ prob["Body"]) && (And @@ ((ReleaseHold[#] != 0) & /@ prob["Dens"])) &&
    ReleaseHold[Lookup[prob, "Assume", Hold[True]]];
  r = TimeConstrained[Quiet[If[dom === All, FindInstance[cond, prob["Vars"]],
    FindInstance[cond, prob["Vars"], dom]]], t, None];
  If[MatchQ[r, {{__Rule}, ___}], Thread[prob["Vars"] -> (prob["Vars"] /. First[r])], None]];

(* A witness rule list (kernel-eligible first), or None, within budget
   seconds: sums at n = a .. a + 30; the exact grid; FindInstance over
   Integers / Rationals (domain re-checked); then, Exact-only, FindInstance
   over C (char-0 form) or the reals (ordered form). *)
searchCounterexample[prob_, budget_] := Module[{t0 = AbsoluteTime[], left, w, df},
  If[$witnessOverride =!= None, Return[$witnessOverride, Module]];
  left := budget - (AbsoluteTime[] - t0);  (* the time still available *)
  df = prob["DomainForm"];
  w = TimeConstrained[Switch[prob["Form"],
      "Closed", If[closedTruth[prob] === False, {}, None],
      "Induction", Module[{a = prob["Shift"], hit = None},
        Do[If[fastWitnessQ[prob, {First[prob["Vars"]] -> n0}],
          hit = {First[prob["Vars"]] -> n0}; Break[]], {n0, a, a + 30}];
        hit],
      _, gridSearch[prob, If[df === "Int", $intGrid, $ratGrid]]], budget, None];
  If[w =!= None || prob["Form"] =!= "Forall", Return[w, Module]];
  w = findWitness[prob, If[df === "Int", Integers, Rationals], left];
  If[w =!= None && witnessOK[prob, w] && kernelWitnessQ[prob, w], Return[w, Module]];
  w = Switch[df, "CharZeroField", findWitness[prob, All, left],
    "OrderedField", findWitness[prob, Reals, left], _, None];
  If[w =!= None && witnessOK[prob, w], w, None]];

(* ---------------------------------------------------------------------- *)
(* Refutation IR and trusted steps                                          *)
(* ---------------------------------------------------------------------- *)

notHeld[s_HoldForm] := s /. HoldForm[e_] :> HoldForm[Not[e]];

$refDensStep = "the witness lies in the claim's domain of definition (every denominator is \
nonzero there; checked by Mathilda's exact arithmetic)";
$refStmtStep = "the refuted statement is the claim as Certify states it (\"RefutedStatement\"; \
a rational claim with its denominators cleared): it is false at the witness, so the claim is \
false there (Mathilda's statement of the claim, not checked by the kernel)";
$refSumStep = "Lean's recursive definitions of the term and the partial sum (in the source) \
are the claim's term and Sum from the lower limit (Mathilda's translation)";
$refExactStep = "Mathilda's exact arithmetic evaluates the statement at the witness: an exact \
check by the system that found it, not a kernel proof";
$refAntiStep = "D[F, x] and f are rational numbers at the witness and differ (Mathilda's D \
and exact arithmetic), so F is not an antiderivative of f; the claim Integrate[f, x] == F is \
certified as D[F, x] == f";
$refAssumeStep = "the witness satisfies `1` (the call's Domain and Assumptions, which the \
refuted statement omits; checked by Mathilda's exact arithmetic), so it is a counterexample to \
the claim as stated";
$refAlgStep = "the nonzero tests at an algebraic witness: RootReduce and a 50-digit numerical \
value (no kernel proof in core Lean for this witness)";

refutationIR[ir_, prob_, w_] := <|"Kind" -> "Refutation", "Input" -> ir["Input"],
  "Refuted" -> ir, "Problem" -> prob, "Witness" -> w, "Kernel" -> kernelWitnessQ[prob, w],
  "Statement" -> prob["Statement"], "Theorem" -> notHeld[prob["Statement"]],
  "TrustedSteps" -> Join[{$refStmtStep}, If[prob["Dens"] === {}, {}, {$refDensStep}],
    If[KeyExistsQ[prob, "Assume"], {StringReplace[$refAssumeStep, "`1`" -> prob["AssumeText"]]}, {}],
    If[prob["Form"] === "Induction", {$refSumStep}, {}]]|>;

(* The budget of the automatic search. *)
searchBudget[t_?NumericQ] /; t > 0 := t;
searchBudget[_] := 2;

(* ---------------------------------------------------------------------- *)
(* Certify: the counterexample after a kernel rejection                     *)
(* ---------------------------------------------------------------------- *)

(* out: the ProofObjects (aligned with the call), items/vs/engines as in
   certifyItems.  For every item a kernel engine rejected (verdict False):
   search, check the witness, refute (one Lean run for all refutations of
   the call), attach the keys and issue Certify::refuted. *)
refuteAfterRejection[out_, items_List, vs_List, engines_List, o_] := Module[{res = out,
    kern, ker, rirs = <||>, w, prob, ir, j, run, rv, keys, a, enc},
  kern = Select[Range[Length[engines]], $backends[engines[[#]]]["Trust"] === "Kernel" &];
  ker = SelectFirst[engines[[kern]], KeyExistsQ[$backends[#], "Refute"] &, None];
  Do[
    If[AnyTrue[kern, vs[[#, j]]["Verified"] === False &],
      ir = items[[j]]["G"]; prob = refProblem[ir];
      If[AssociationQ[prob], prob = withAssume[prob, ir, o]];
      If[AssociationQ[prob],
        w = searchCounterexample[prob, searchBudget[o["Counterexample"]]];
        If[w =!= None && witnessOK[prob, w],
          rirs[j] = Join[refutationIR[ir, prob, w], <|"LeanStatement" -> ir["Statement"],
            "Preamble" -> Lookup[ir, "Preamble", ""]|>]]]],
    {j, Length[items]}];
  If[rirs === <||>, Return[res, Module]];
  (* The kernel refutations of the call, in one run. *)
  run = Select[Keys[rirs], rirs[#]["Kernel"] && ker =!= None &];
  rv = <||>;
  If[run =!= {},
    enc = <|"Index" -> items[[#]]["Index"], "G" -> $backends[ker]["Refute"][rirs[#], o]|> & /@ run;
    rv = AssociationThread[run -> Quiet[$backends[ker]["Run"][enc, o], {Certify::rejected}]]];
  Do[
    With[{r = rirs[j], v = Lookup[rv, j, None]},
      keys = <|"Counterexample" -> r["Witness"],
        "Refuted" -> If[AssociationQ[v] && v["Verified"] === True, True, Missing["Exact", True]],
        "RefutationTheorem" -> r["Theorem"]|>;
      If[AssociationQ[v], keys = Append[keys, "RefutationSource" -> v["Source"]]];
      Message[Certify::refuted, items[[j]]["G"]["Input"], r["Witness"]];
      a = Join[First[res[[items[[j]]["Index"]]]], keys];
      res[[items[[j]]["Index"]]] = issueObject[a]],
    {j, Keys[rirs]}];
  res];

issueObject[a_] := With[{b = a}, $issued[Hash[b]] = True; ProofObject[b]];

(* ---------------------------------------------------------------------- *)
(* Refute                                                                   *)
(* ---------------------------------------------------------------------- *)

Options[Refute] = {"Engine" -> "Lean", "Domain" -> Automatic, Assumptions -> {},
  TimeConstraint -> 60, "Verify" -> True, "SearchTime" -> 5, "LeanExecutable" -> Automatic,
  "LeanProject" -> Automatic};

SetAttributes[Refute, HoldAll];

Refute[c_, opts:OptionsPattern[]] := If[heldOptionProblemFor[Refute, HoldComplete[opts]] === $Failed,
  $Failed, refuteClaim[HoldComplete[c], {opts}]];

(* Refute's options as the association the generators and backends read. *)
refOptions[opts_List] := Module[{keys = {"Engine", "Domain", "Assumptions", "TimeConstraint",
    "Verify", "SearchTime", "LeanExecutable", "LeanProject"}, v},
  If[optionNameProblem[Refute, opts] =!= None, Return[$Failed, Module]];
  v = AssociationThread[keys, OptionValue[Refute, opts, #] & /@ keys];
  Join[v, <|"Tier" -> Automatic, "Method" -> Automatic, "Counterexample" -> False,
    "Engines" -> {v["Engine"]}|>]];

$refuteShapes = HoldComplete[_Equal | _Implies | _ForAll | _GreaterEqual | _LessEqual];

(* The IR of a claim, statement only (no certificate search), or $Failed
   after a message. *)
refuteIR[HoldComplete[c:Equal[Integrate[f_, x_Symbol], F_]], o_] :=
  <|"Kind" -> "Antiderivative", "Input" -> HoldForm[c], "Hypotheses" -> {}, "Domain" -> Reals,
    "DomainForm" -> "Rat", "Obligation" -> <|"Type" -> "Antiderivative", "F" -> F, "f" -> f,
      "Variable" -> x|>|>;
refuteIR[hc_HoldComplete, o_] := If[!MatchQ[hc, $refuteShapes],
  Message[Refute::nostmt, Replace[hc, HoldComplete[e_] :> HoldForm[e]]]; $Failed,
  Block[{$statementOnly = True}, generate[hc, o]]];

refuteClaim[hc_HoldComplete, opts_List] := Module[{o, rw, ir, prob, w, t0, ts, rir, input},
  t0 = AbsoluteTime[];
  o = refOptions[opts];
  If[o === $Failed, Return[$Failed, Module]];
  If[!MemberQ[{"Lean", "Exact"}, o["Engine"]],
    Message[Refute::engine, o["Engine"]]; Return[$Failed, Module]];
  If[!AnyTrue[$domains, # === o["Domain"] &],
    Message[Certify::domain, ToString[o["Domain"], InputForm]]; Return[$Failed, Module]];
  input = Replace[hc, HoldComplete[e_] :> HoldForm[e]];
  rw = rewriteInput[hc];
  If[rw === $Failed, Return[$Failed, Module]];
  ir = refuteIR[First[rw], o];
  If[!AssociationQ[ir], Return[$Failed, Module]];
  prob = refProblem[ir];
  If[!AssociationQ[prob], Message[Refute::nostmt, input]; Return[$Failed, Module]];
  (* A Domain / Assumptions the statement omits restricts the witness. *)
  prob = withAssume[prob, ir, o];
  w = searchCounterexample[prob, searchBudget[o["SearchTime"]]];
  If[w === None || !witnessOK[prob, w], Message[Refute::none, input]; Return[$Failed, Module]];
  rir = Append[refutationIR[ir, prob, w], "Input" -> input];
  ts = AbsoluteTime[] - t0;
  If[o["Engine"] === "Lean" && rir["Kernel"], refuteLean[rir, ir, o, ts],
    refuteExact[rir, o, ts]]];

(* The Lean refutation: the statement's Lean text from the Lean encoder,
   then the backend's "Refute" text and one run. *)
refuteLean[rir_, ir_, o_, ts_] := Module[{art, g, item, v, named},
  art = leanEncode[ir, o];
  If[!AssociationQ[art], Return[$Failed, Module]];
  g = $backends["Lean"]["Refute"][Join[rir, <|"LeanStatement" -> art["Statement"],
    "Preamble" -> Lookup[art, "Preamble", ""]|>], o];
  item = <|"Index" -> 1, "G" -> g|>;
  v = First[Quiet[$backends["Lean"]["Run"][{item}, o], {Certify::rejected}]];
  named = leanNamed[item];
  issueObject[Join[First[mkProofObject[Join[refutationG[rir, ir], <|"Proof" -> named["Proof"]|>],
      Append[v, "Timing" -> v["Timing"] + ts], "Lean", 0, v["LeanVersion"], v["Source"],
      named["Statement"]]],
    <|"Counterexample" -> rir["Witness"], "RefutedStatement" -> rir["Statement"]|>]]];

(* The Exact refutation: the backend re-evaluates the witness. *)
refuteExact[rir_, o_, ts_] := With[{v = $backends["Exact"]["Refute"][rir, o], g = refutationG[rir, rir["Refuted"]]},
  issueObject[<|"Kind" -> "Refutation", "Input" -> g["Input"], "Result" -> False,
    "Theorem" -> g["Theorem"], "Certificate" -> g["Certificate"], "Checks" -> v["Checks"],
    "Verified" -> v["Verified"], "Engine" -> "Exact", "Tier" -> "Exact",
    "Domain" -> g["Domain"], "Hypotheses" -> {},
    "TrustedSteps" -> DeleteDuplicates[Join[g["TrustedSteps"], {$refExactStep},
      If[AllTrue[rir["Witness"][[All, 2]], ratNumberQ], {}, {$refAlgStep}],
      If[rir["Problem"]["Form"] === "Antiderivative", {$refAntiStep}, {}]]],
    "Timing" -> v["Timing"] + ts, "Cached" -> False,
    "Counterexample" -> rir["Witness"], "RefutedStatement" -> rir["Statement"]|>]];

(* The generator-like association of a refutation (mkProofObject's input). *)
refutationG[rir_, ir_] := <|"Kind" -> "Refutation", "Input" -> rir["Input"], "Result" -> False,
  "Theorem" -> rir["Theorem"],
  "Certificate" -> <|"Witness" -> rir["Witness"], "Form" -> rir["Problem"]["Form"],
    "RefutedKind" -> ir["Kind"], "DomainForm" -> rir["Problem"]["DomainForm"],
    "Antecedents" -> rir["Problem"]["Antecedents"]|>,
  "Domain" -> Lookup[ir, "Domain", Rationals], "Hypotheses" -> {},
  "TrustedSteps" -> rir["TrustedSteps"]|>;

End[];

EndPackage[];
