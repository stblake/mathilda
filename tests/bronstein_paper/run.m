(* tests/bronstein_paper/run.m -- check every example of examples.m.
   ./Mathilda -file tests/bronstein_paper/run.m  prints one line per example
   (PASS|FAIL <id> <seconds> <detail>) and "SUMMARY pass fail". *)
Integrate[Sqrt[x], x, Method -> "Bronstein"];   (* loads the package *)
Get["tests/bronstein_paper/examples.m"];
bpCheck[r_Association] := Module[{f = ReleaseHold[r["f"]], res, ok, det = ""},
  res = TimeConstrained[Switch[r["kind"],
      "Integrate", Bronstein`BronsteinIntegrate[f, x],
      "EIQ", Risch`ElementaryIntegralQ[f, x],
      "RDE", Risch`RischDECertified[f[[1]], f[[2]], x],
      "Scope", Bronstein`BronsteinScope[f, x]], 300, "TIMEOUT"];
  ok = Switch[r["kind"],
    "Integrate", Switch[r["expect"],
       "elementary", !ListQ[res] && res =!= "TIMEOUT" && TrueQ[PossibleZeroQ[D[res, x] - f]],
       "nonelementary", MatchQ[res, {"not elementary", _}],
       "incomplete", MatchQ[res, {"implementation incomplete", _}]],
    "EIQ", res === r["expect"],
    "RDE", If[r["expect"] === "NoSolution", res === "NoSolution", TrueQ[PossibleZeroQ[res - ReleaseHold[r["expect"]]]]],
    "Scope", If[r["expect"] === "incomplete", MatchQ[res, {"implementation incomplete", _}], res === r["expect"]]];
  If[!ok, det = StringTake[ToString[res, InputForm], UpTo[160]]];
  {ok, det}];
bpPass = 0; bpFail = 0;
Do[With[{t = AbsoluteTiming[bpCheck[r]]},
    If[t[[2, 1]], bpPass++, bpFail++];
    Print[If[t[[2, 1]], "PASS ", "FAIL "], r["id"], " ", NumberForm[t[[1]], {6, 2}], " ", t[[2, 2]]];
    (* one tab-separated row for tools/bronstein_paper_table.py *)
    Print["ROW\t", StringRiffle[{r["id"], r["gap"], r["theory"], r["version"], r["kind"], ToString[r["expect"], InputForm],
      StringReplace[ToString[r["f"], InputForm], {"Hold[" ~~ s__ ~~ "]" :> s}], ToString[NumberForm[t[[1]], {6, 2}]],
      If[t[[2, 1]], "PASS", "FAIL"], r["before"], r["note"]}, "\t"]]],
  {r, bronsteinPaperExamples}];
Print["SUMMARY ", bpPass, " ", bpFail];
