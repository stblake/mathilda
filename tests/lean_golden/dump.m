(* ====================================================================== *)
(* tests/lean_golden/dump.m -- goldenDump[] records the golden corpus.     *)
(*                                                                          *)
(* goldenDump[] evaluates every item of goldenCorpus (corpus.m) in order   *)
(* and returns one String, per item:                                        *)
(*   === <i> <the item, InputForm>                                          *)
(*   Messages: <the messages the item issued, InputForm>                    *)
(*   Head: ProofObject | List[<n>] | <head>                                 *)
(* then, per ProofObject (each element of a list, after "--- <j>"), every   *)
(* key of obj["Properties"] except "Timing" as "<Key>: <value, InputForm>"  *)
(* and its summary box with the timing replaced by "T s"; any other value   *)
(* as "Value: <InputForm>".  The inlined prelude of "LeanSource" is         *)
(* replaced by <<PRELUDE>> after checking that the source contains it       *)
(* (PRELUDE-MISMATCH otherwise).  Items run without Quiet, so               *)
(* $MessageList sees their messages; $LeanProject and $LeanExecutable are   *)
(* fixed, so the record does not depend on the environment.                 *)
(* ====================================================================== *)

goldenText[e_] := ToString[e, InputForm];

goldenSource[src_String] := With[{p = Lean`Private`prelude[]},
  If[StringQ[p] && p =!= "" && StringContainsQ[src, p],
    StringReplace[src, p -> "<<PRELUDE>>"], "PRELUDE-MISMATCH " <> src]];
goldenSource[other_] := other;

goldenSummary[obj_] := StringReplace[ToString[obj],
  RegularExpression["\\| [^|]* s\\]$"] -> "| T s]"];

goldenValue[obj_ProofObject] := Flatten[{
  Function[k, k <> ": " <> goldenText[If[k === "LeanSource", goldenSource[obj[k]], obj[k]]]] /@
    DeleteCases[obj["Properties"], "Timing"],
  "Summary: " <> goldenSummary[obj]}];
goldenValue[other_] := {"Value: " <> goldenText[other]};

goldenRecord[i_Integer, h_Hold] := Module[{r, msgs, body},
  msgs = Block[{$MessageList = {}}, r = ReleaseHold[h]; $MessageList];
  body = Which[
    Head[r] === ProofObject, Prepend[goldenValue[r], "Head: ProofObject"],
    ListQ[r], Flatten[{"Head: List[" <> ToString[Length[r]] <> "]",
      Table[{"--- " <> ToString[j], goldenValue[r[[j]]]}, {j, Length[r]}]}],
    True, Prepend[goldenValue[r], "Head: " <> goldenText[Head[r]]]];
  StringRiffle[Flatten[{"=== " <> ToString[i] <> " " <> goldenText[HoldForm @@ h],
    "Messages: " <> goldenText[msgs], body}], "\n"]];

goldenDump[] := Block[{$LeanProject = None, $LeanExecutable = Automatic},
  SetOptions[Certify, Assumptions -> {}];
  StringJoin[Table[goldenRecord[i, Extract[goldenCorpus, i, Hold]] <> "\n",
    {i, Length[goldenCorpus]}]]];
