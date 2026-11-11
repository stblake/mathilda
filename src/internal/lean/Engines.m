(* ====================================================================== *)
(* Engines.m -- engine selection, the heuristic engine and multi-engine    *)
(* consensus for Certify (loaded by Lean.m).                                *)
(*                                                                          *)
(* "Engine" is one registered backend name or a non-empty duplicate-free    *)
(* list of them.  With a list every engine runs (list order) and the        *)
(* ProofObject carries the combined verdict, "Engines" (name -> verdict)    *)
(* and "EngineDetails".  Only a backend with "Trust" -> "Kernel" may make   *)
(* the combined verdict True.  All messages are in Lean.m.  The file opens  *)
(* the package itself (every symbol is in Lean`Private`).                   *)
(* ====================================================================== *)

BeginPackage["Lean`"];
Begin["`Private`"];

(* ---------------------------------------------------------------------- *)
(* Engine selection                                                        *)
(* ---------------------------------------------------------------------- *)

(* The engines of an "Engine" value, or $Failed (Certify::engine). *)
engineList[e_String] := If[KeyExistsQ[$backends, e], {e}, $Failed];
engineList[l_List] := If[l =!= {} && Length[DeleteDuplicates[l]] === Length[l] &&
  AllTrue[l, StringQ[#] && KeyExistsQ[$backends, #] &], l, $Failed];
engineList[_] := $Failed;

(* Certify::engineopt for an explicit Method or "Tier" that no listed engine
   reads (its "Options"); LeanExecutable, LeanProject and TimeConstraint are
   ignored silently. *)
engineOptMessages[o_, engines_List] := Do[
  If[o[r[[1]]] =!= Automatic &&
      !AnyTrue[engines, MemberQ[Lookup[$backends[#], "Options", {}], r[[1]]] &],
    Message[Certify::engineopt, r[[2]] -> o[r[[1]]],
      If[Length[engines] === 1, First[engines], engines]]],
  {r, {{"Method", Method}, {"Tier", "Tier"}}}];

(* A True from a backend that is not a kernel is shown as Missing[name, True]. *)
guardVerdicts[e_, vs_List] := If[$backends[e]["Trust"] === "Kernel", vs,
  Map[If[#["Verified"] === True, Append[#, "Verified" -> Missing[e, True]], #] &, vs]];

backendSteps[e_, g_] := Lookup[$backends[e], "TrustedSteps", ({} &)][g];

(* ---------------------------------------------------------------------- *)
(* Engine -> "Mathilda": the heuristic                                     *)
(* ---------------------------------------------------------------------- *)

(* PossibleZeroQ / Simplify on lhs - rhs of the claim (lists compare
   entrywise); a kind with no claim is not checked. *)
heuristicZeroQ[d_] := TrueQ[Quiet[PossibleZeroQ[d]]] || TrueQ[Quiet[Simplify[d]] === 0];
(* "Verified" is never True without Lean: the heuristic's own answer is
   carried as Missing["Heuristic", True | False]. *)
heuristicVerdict[g_] := Module[{v = Missing["NotChecked"], t = 0.},
  If[ListQ[g["Claim"]],
    t = First[AbsoluteTiming[v = Missing["Heuristic", AllTrue[
      Flatten[{g["Claim"][[1]] - g["Claim"][[2]]}], heuristicZeroQ]]]]];
  <|"Verified" -> v, "Axioms" -> {}, "Logs" -> {}, "Timing" -> t|>];

$heuristicStep = "Mathilda's PossibleZeroQ / Simplify (a heuristic check, not a proof)";

(* The heuristic shows the Lean statement and batch source (it reads the
   Lean encoding, hence "Tier" and Method). *)
heuristicRun[items_List, o_] := MapThread[Append[heuristicVerdict[#1["G"]], "Source" -> #2] &,
  {items, leanBatchSources[items, o]}];
heuristicPresent[item_, v_, o_] := With[{g = leanNamed[item]},
  mkProofObject[Append[g, "TrustedSteps" -> {$heuristicStep}],
    v, "Mathilda", "Heuristic", Missing["NotChecked"], v["Source"], g["Statement"]]];

registerBackend["Mathilda", <|"Name" -> "Mathilda", "Trust" -> "Heuristic",
  "Encode" -> leanEncode, "Run" -> heuristicRun, "Present" -> heuristicPresent,
  "Label" -> ("Mathilda Heuristic" &), "TrustedSteps" -> ({$heuristicStep} &),
  "Options" -> {"Tier", "Method"}|>];

(* ---------------------------------------------------------------------- *)
(* Multi-engine consensus                                                  *)
(* ---------------------------------------------------------------------- *)

(* The combined verdict of the verdict values vs (aligned with engines).
   Heuristic engines only contribute their entry.  Over the others, with K
   the kernel engines:
     1. a False or Missing[_, False]: False; Certify::disagree when a kernel
        says True, or a kernel says False while another engine says
        Missing[_, True];
     2. K non-empty and every kernel True: True;
     3. K empty: Missing["Consensus", True] when every verdict is
        Missing[_, True], else Missing["NotChecked"];
     4. else the first kernel verdict that is not True (Missing[...]).
   Only heuristic engines: the first one's verdict. *)
consensus[g_, vs_List, engines_List] := Module[{trust, idx, kern, falseQ},
  trust = $backends[#]["Trust"] & /@ engines;
  idx = Select[Range[Length[engines]], trust[[#]] =!= "Heuristic" &];
  If[idx === {}, Return[First[vs], Module]];
  kern = Select[idx, trust[[#]] === "Kernel" &];
  falseQ[v_] := v === False || MatchQ[v, Missing[_, False]];
  Which[
    AnyTrue[vs[[idx]], falseQ],
      If[MemberQ[vs[[kern]], True] || (MemberQ[vs[[kern]], False] &&
          AnyTrue[vs[[Complement[idx, kern]]], MatchQ[#, Missing[_, True]] &]),
        Message[Certify::disagree, g["Input"], AssociationThread[engines -> vs]]];
      False,
    kern =!= {} && AllTrue[vs[[kern]], # === True &], True,
    kern === {}, If[AllTrue[vs[[idx]], MatchQ[#, Missing[_, True]] &],
      Missing["Consensus", True], Missing["NotChecked"]],
    True, SelectFirst[vs[[kern]], # =!= True &]]];

(* The multi-engine ProofObject: the presentation of Lean (if listed), else
   of Exact (if listed), else of the first engine, with the combined
   "Verified", "Engine" -> the list, "Engines", "EngineDetails", the summed
   "Timing" and "TrustedSteps" = the generator's plus the backend-level
   steps of the kernel engines (of every engine when none is a kernel, and
   then without the Lean-specific Lucas step). *)
mkMultiObject[item_, vs_List, engines_List, o_] := Module[{g = item["G"], base, a, kern},
  base = Which[MemberQ[engines, "Lean"], "Lean", MemberQ[engines, "Exact"], "Exact",
    True, First[engines]];
  a = First[$backends[base]["Present"][item, vs[[First[Flatten[Position[engines, base]]]]], o]];
  kern = Select[engines, $backends[#]["Trust"] === "Kernel" &];
  a = Join[a, <|"Verified" -> consensus[g, #["Verified"] & /@ vs, engines],
    "Engine" -> engines,
    "TrustedSteps" -> DeleteDuplicates[Join[
      If[kern === {}, DeleteCases[g["TrustedSteps"], $lucasStep], g["TrustedSteps"]],
      Flatten[backendSteps[#, g] & /@ If[kern === {}, engines, kern]]]],
    "Timing" -> Total[Select[#["Timing"] & /@ vs, NumericQ]],
    "Engines" -> AssociationThread[engines -> (#["Verified"] & /@ vs)],
    "EngineDetails" -> AssociationThread[engines -> MapThread[Join[
      <|"Verified" -> #2["Verified"], "Timing" -> #2["Timing"],
        "TrustedSteps" -> backendSteps[#1, g], "Logs" -> Lookup[#2, "Logs", {}]|>,
      KeyDrop[#2, {"Verified", "Timing", "Logs", "Source"}]] &, {engines, vs}]]|>];
  With[{b = a}, $issued[Hash[b]] = True; ProofObject[b]]];

End[];

EndPackage[];
