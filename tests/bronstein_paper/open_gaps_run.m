(* tests/bronstein_paper/open_gaps_run.m -- check the open-gap register:
   STILL-OPEN when an entry behaves as recorded, CHANGED otherwise (a gap
   closed or moved: move it to examples.m with its fix, or update the record).
   ./Mathilda -file tests/bronstein_paper/open_gaps_run.m *)
Integrate[Sqrt[x], x, Method -> "Bronstein"];
Get["tests/bronstein_paper/open_gaps.m"];
ogClass[r_Association] := Module[{f, res, t},
  If[r["kind"] === "Expr",
    res = TimeConstrained[Catch[ReleaseHold[r["f"]], _, "Throw" &], r["budget"], "TIMEOUT"];
    Return[If[res === "Throw", "Throw", ToString[res, InputForm]]]];
  f = ReleaseHold[r["f"]];
  {t, res} = AbsoluteTiming[TimeConstrained[Bronstein`BronsteinIntegrate[f, x], r["budget"], "TIMEOUT"]];
  Which[
    res === "TIMEOUT", "timeout",
    MatchQ[res, {"implementation incomplete", _}], "incomplete",
    MatchQ[res, {"failed", _String}] && StringStartsQ[res[[2]], "internal"], "internal",
    MatchQ[res, {"not elementary", _}], "nonelementary",
    ListQ[res], "failed",
    !FreeQ[res, Root | RootSum], "presentation",
    True, "answer"]];
ogOpen = 0; ogChanged = 0;
Do[With[{c = ogClass[r]},
    If[c === r["current"], ogOpen++, ogChanged++];
    Print[If[c === r["current"], "STILL-OPEN ", "CHANGED    "], r["id"], "  recorded: ", r["current"], "  now: ", c];
    Print["GAPROW\t", StringRiffle[{r["id"], r["area"], r["gap"], StringReplace[ToString[r["f"], InputForm], {"Hold[" ~~ s__ ~~ "]" :> s}],
      r["current"], r["followup"], r["source"], If[c === r["current"], "open", "changed"]}, "\t"]]],
  {r, bronsteinOpenGaps}];
Do[Print["GAPSTRUCT\t", StringRiffle[g, "\t"]], {g, openGapsStructural}];
Print["OPENGAPS ", ogOpen, " ", ogChanged];
