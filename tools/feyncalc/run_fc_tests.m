(* ::Package:: *)
(* ================================================================== *)
(* run_fc_tests.m -- run one FeynCalc unit-test driver (.mt) under     *)
(* Mathilda, unmodified, with a minimal MUnit`Test.                    *)
(*                                                                    *)
(*   FC_ROOT=<abs path of packages/feyncalc> FC_MT=<abs path of .mt>  *)
(*   [FC_ONLYSUB="re1|re2"] [FC_TIMEOUT=60] Mathilda -file run_fc_tests.m *)
(*                                                                    *)
(* Normally driven by tools/check_feyncalc.py.  Prints one line per   *)
(* test on stdout:                                                    *)
(*   FC <TAB> verdict <TAB> id <TAB> micros                           *)
(* verdict: PASS, FAIL (wrong value), MSG (right value, unexpected    *)
(* messages), ABORT, TIMEOUT, SKIP (typesetting: out of scope).       *)
(* A non-pass is followed by a detail line                            *)
(*   FCD <TAB> id <TAB> got <TAB> expected <TAB> messages             *)
(* and the run ends with "FCEND".                                     *)
(*                                                                    *)
(* Semantics follow MUnit's Test: HoldAllComplete; the input is       *)
(* evaluated, then the expected value, and they are compared with     *)
(* EquivalenceFunction (default SameQ); the messages the input emits  *)
(* are compared with the expected-messages argument (default {}) or   *)
(* passed to MessagesEquivalenceFunction.  One deliberate difference: *)
(* the drivers build most tests as test[in, out, ...] /. test -> Test, *)
(* which evaluates every input eagerly, outside Test.  `test` is made *)
(* HoldAllComplete before the driver is read, so each input runs      *)
(* inside its own Test -- under CheckAbort, TimeConstrained and the   *)
(* message capture -- and one bad test cannot abort its whole group.  *)
(* ================================================================== *)

$FCHarnessRoot = Environment["FC_ROOT"];
$FCHarnessMT = Environment["FC_MT"];
$FCHarnessTimeout = With[{t = Environment["FC_TIMEOUT"]},
  If[StringQ[t], ToExpression[t], 60]];

If[!StringQ[$FCHarnessRoot] || !StringQ[$FCHarnessMT],
  WriteString["stdout", "FCERROR\tFC_ROOT and FC_MT must be set\n"]; Exit[2]];

BeginPackage["MUnit`"];
Test::usage = "Test[input, expected, messages, opts] checks that input evaluates to expected.";
TestID::usage = "TestID is an option to Test naming the test.";
EquivalenceFunction::usage = "EquivalenceFunction is an option to Test (default SameQ).";
MessagesEquivalenceFunction::usage = "MessagesEquivalenceFunction is an option to Test.";
TestFailureMessage::usage = "TestFailureMessage is an option to Test.";
Begin["`Private`"];

SetAttributes[Test, HoldAllComplete];

count = 0;

(* One line with no tabs or newlines, at most n characters. *)
flat[e_, n_: 400] := Module[{s = ToString[e, InputForm]},
  s = StringReplace[s, {"\t" -> " ", "\n" -> " ", "\r" -> " "}];
  If[StringLength[s] > n, StringTake[s, n] <> "...", s]];

msgString[l_List] := StringRiffle[ToString /@ l, ","];
msgString[_] := "";

(* The time goes out as integer microseconds: the same digits in both kernels,
   whatever their Real printing; check_feyncalc.py converts it to seconds. *)
emit[verdict_, id_, secs_] :=
  WriteString["stdout", "FC\t", verdict, "\t", id, "\t", ToString[Round[10^6 secs]], "\n"];
detail[id_, got_, exp_, msgs_] :=
  WriteString["stdout", "FCD\t", id, "\t", got, "\t", exp, "\t", msgString[msgs], "\n"];

(* Evaluate held code under the per-test budget; an abort or a timeout
   becomes a sentinel instead of unwinding the driver. *)
guarded[h_Hold] := CheckAbort[
  TimeConstrained[ReleaseHold[h], $FCHarnessTimeoutValue, $timedOut], $aborted];

Test[input_, expected_, rest___] := Module[{args = Hold[rest], opts, id},
  count++;
  opts = Cases[args, (Rule | RuleDelayed)[_, _], {1}];
  id = Replace[TestID /. opts, TestID -> "unnamed-" <> ToString[count]];
  If[!StringQ[id], id = ToString[id]];
  If[StringMatchQ[id, "fcstSharedObjectsTypesetting*"],
    emit["SKIP", id, 0],
    runOne[Hold[input], Hold[expected], args, opts, id]];
  Null];

runOne[input_Hold, expected_Hold, args_Hold, opts_, id_] := Module[
  {pos, eqf, mef, emsgs, t0, actual, msgs = {}, exp = $none, valueOK, msgsOK, verdict},
  $FCHarnessTimeoutValue = Global`$FCHarnessTimeout;
  (* The positional argument after the expected value: expected messages. *)
  pos = Cases[args, Except[(Rule | RuleDelayed)[_, _]], {1}, Heads -> False];
  eqf = Replace[EquivalenceFunction /. opts, EquivalenceFunction -> SameQ];
  mef = Replace[MessagesEquivalenceFunction /. opts, MessagesEquivalenceFunction -> Automatic];
  emsgs = If[pos === {}, {}, guarded[Hold @@ {First[pos]}]];
  t0 = AbsoluteTime[];
  actual = Block[{$MessageList = {}},
    With[{r = guarded[input]}, msgs = $MessageList; r]];
  Which[
    actual === $aborted, verdict = "ABORT",
    actual === $timedOut, verdict = "TIMEOUT",
    True,
      exp = guarded[expected];
      valueOK = TrueQ[guarded[Hold[eqf[actual, exp]]]];
      msgsOK = If[mef === Automatic,
        (ToString /@ msgs) === (ToString /@ Flatten[{emsgs}]),
        TrueQ[guarded[Hold[mef[msgs, emsgs]]]]];
      verdict = Which[!valueOK, "FAIL", !msgsOK, "MSG", True, "PASS"]];
  emit[verdict, id, AbsoluteTime[] - t0];
  If[verdict =!= "PASS",
    detail[id, If[verdict === "FAIL" || verdict === "MSG", flat[actual], verdict],
      If[exp === $none, "", flat[exp]], msgs]];
  Null];

End[];
EndPackage[];

(* The drivers' eager `test[...]` placeholder: held, so it reaches Test
   unevaluated (see the header). *)
SetAttributes[test, HoldAllComplete];

$FeynCalcStartupMessages = False;
$FeynCalcPackageRoot = $FCHarnessRoot;
Get[FileNameJoin[{$FCHarnessRoot, "load.m"}]];
If[!MemberQ[$Packages, "FeynCalc`"],
  WriteString["stdout", "FCERROR\tFeynCalc did not load\n"]; Exit[3]];

(* The upstream TestSuite.m hook, built as it builds it from onlySubTest:
   "a|b" -> {RegularExpression[".*a.*"], RegularExpression[".*b.*"]}, or "". *)
$OnlySubTest = With[{o = Environment["FC_ONLYSUB"]},
  If[StringQ[o] && o =!= "",
    RegularExpression[(".*" <> # <> ".*")] & /@ StringSplit[o, "|"], ""]];

Get[$FCHarnessMT];
WriteString["stdout", "FCEND\n"];
Exit[0];
