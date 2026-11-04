(* ::Package:: *)
(* ------------------------------------------------------------------ *)
(* run_df_cases.m -- driver for tools/check_differentialforms.py.      *)
(*                                                                    *)
(* Loads the DifferentialForms package (Zizza/Jentschura, issue #125), *)
(* then the case corpus, and evaluates every case under               *)
(* TimeConstrained with its messages collected. Uses only portable     *)
(* language, so the same driver can record expected answers on another *)
(* kernel binary; all parameters come from the environment:            *)
(*                                                                    *)
(*   DF_PKG       path of DifferentialForms.m                          *)
(*   DF_CORPUS    path of tests/differentialforms/corpus.m             *)
(*   DF_TIMEOUT   per-case seconds (default 60)                        *)
(*   DF_ONLY      optional comma-separated case ids                    *)
(*   DF_EXPECTED  optional .m file setting $DFExpected (compare mode)  *)
(*                                                                    *)
(* Dump mode (no DF_EXPECTED) prints, per case,                       *)
(*   DFR <tab> id <tab> kind <tab> FullForm <tab> messages             *)
(* Compare mode prints                                                *)
(*   DF <tab> verdict <tab> id <tab> kind <tab> got <tab> messages     *)
(* where verdict is PASS, FAIL, MSG (right answer, unexpected message) *)
(* or TIMEOUT.                                                         *)
(*                                                                    *)
(* Case kinds (see corpus.m):                                         *)
(*   "recorded"  compared with the recorded expected answer            *)
(*   "identity"  must evaluate to True                                 *)
(*   "expect"    compared with the corpus' own HoldComplete expected   *)
(* ------------------------------------------------------------------ *)

$DFPkg = Environment["DF_PKG"];
$DFCorpus = Environment["DF_CORPUS"];
$DFTimeout = With[{t = Environment["DF_TIMEOUT"]}, If[StringQ[t], ToExpression[t], 60]];
$DFOnly = With[{o = Environment["DF_ONLY"]},
  If[StringQ[o] && o =!= "", StringSplit[o, ","], All]];
$DFExpectedFile = Environment["DF_EXPECTED"];

If[!StringQ[$DFPkg] || !StringQ[$DFCorpus],
  WriteString["stdout", "DFERROR\tDF_PKG and DF_CORPUS must be set\n"]; Exit[2]];

(* Loading may warn (a kernel whose TensorProduct is Protected rejects the
   package's rules for it); the load is not under test. *)
Quiet[Get[$DFPkg]];
If[!MemberQ[$ContextPath, "DifferentialForms`"],
  WriteString["stdout", "DFERROR\tpackage did not load\n"]; Exit[3]];

Get[$DFCorpus];
If[!ListQ[$DFCases], WriteString["stdout", "DFERROR\tcorpus did not define $DFCases\n"]; Exit[3]];

compareMode = StringQ[$DFExpectedFile] && $DFExpectedFile =!= "";
If[compareMode, Get[$DFExpectedFile]];

oneLine[s_String] := StringReplace[s, {"\t" -> " ", "\n" -> " ", "\r" -> " "}];
full[e_] := oneLine[ToString[FullForm[e]]];
(* Message names, as "Head::tag" strings, independent of context. *)
msgNames[msgs_List] := Sort[DeleteDuplicates[Replace[msgs,
  {HoldForm[MessageName[s_, tag_, ___]] :> SymbolName[Unevaluated[s]] <> "::" <> tag,
   other_ :> ToString[other, InputForm]}, {1}]]];

(* Evaluate a held case with its messages and a time limit. *)
run[h_HoldComplete] := Module[{r, msgs},
  r = Block[{$MessageList = {}},
    With[{v = TimeConstrained[ReleaseHold[h], $DFTimeout, $DFTimedOut]},
      msgs = $MessageList; v]];
  {r, msgs}];

(* Semantic equality of two package results. Exact structure first; then the
   difference must vanish after the coefficients of every altform/form/vec
   (the package's linear objects, which merge under Plus) are simplified. *)
dfAlt = DifferentialForms`Private`altform;
dfForm = DifferentialForms`Private`form;
dfVec = DifferentialForms`Private`vec;
dfZeroQ[e_] := Module[{z},
  z = Expand[e] /. (h : (dfAlt | dfForm | dfVec))[c_, v_] :> h[Simplify[c], v];
  z = Simplify[z];
  AllTrue[Flatten[{z}], # === 0 &]];
dfSameQ[a_, b_] := a === b ||
  (Head[a] === List && Head[b] === List && Length[a] === Length[b] &&
     And @@ MapThread[dfSameQ, {a, b}]) ||
  (Head[a] =!= List && Head[b] =!= List && TrueQ[Quiet[dfZeroQ[a - b]]]);

emitDump[id_, kind_, r_, msgs_] :=
  WriteString["stdout", "DFR\t", id, "\t", kind, "\t", full[r], "\t",
    oneLine[StringRiffle[msgNames[msgs], ";"]], "\n"];

emit[verdict_, id_, kind_, r_, msgs_] :=
  WriteString["stdout", "DF\t", verdict, "\t", id, "\t", kind, "\t",
    StringTake[#, Min[StringLength[#], 2000]] &[full[r]], "\t",
    oneLine[StringRiffle[msgNames[msgs], ";"]], "\n"];

(* The loop variable is $-prefixed: Do binds it dynamically, and a plain name
   such as c would capture the corpus' own free symbol c. *)
Do[
  Module[{id = $dfCase[[1]], h = $dfCase[[2]], kind = $dfCase[[3]], r, msgs, ok, exp, refMsgs, extra},
    If[$DFOnly === All || MemberQ[$DFOnly, id],
      {r, msgs} = run[h];
      If[!compareMode,
        emitDump[id, kind, r, msgs],
        (* compare mode *)
        refMsgs = {};
        ok = Which[
          r === $DFTimedOut, "TIMEOUT",
          kind === "identity", If[r === True, "PASS", "FAIL"],
          kind === "expect", If[dfSameQ[r, ReleaseHold[$dfCase[[4]]]], "PASS", "FAIL"],
          kind === "recorded",
            exp = Lookup[$DFExpected, id, Missing[]];
            If[MissingQ[exp], "NOREF",
              refMsgs = exp[[2]];
              If[dfSameQ[r, ToExpression[exp[[1]]]], "PASS", "FAIL"]],
          True, "FAIL"];
        (* An answer that is right but came with a message the recorded run
           (or, for identity/expect cases, the corpus) did not expect is a MSG. *)
        extra = Complement[msgNames[msgs],
          Join[refMsgs, If[Length[$dfCase] >= 5, $dfCase[[5]], {}]]];
        If[ok === "PASS" && extra =!= {}, ok = "MSG"];
        emit[ok, id, kind, r, msgs]]]],
  {$dfCase, $DFCases}];

WriteString["stdout", "DFEND\n"];
