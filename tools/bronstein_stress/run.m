(* run.m -- integrate the cases of the file CASES (lines CASE<TAB>id<TAB>fam<TAB>f)
   with BronsteinIntegrate (METHOD=Bronstein) or plain Integrate (METHOD=Automatic),
   and print one RESULT line per case:
     RESULT<TAB>id<TAB>fam<TAB>class<TAB>seconds<TAB>detail
   class: OK | WRONG | FALSE_CERT | DECLINE | TIMEOUT | UNEVAL.            *)
LoadModule["bronstein/Bronstein.m"];
method = Environment["METHOD"]; If[!StringQ[method], method = "Bronstein"];
budget = ToExpression[Environment["BUDGET"]]; If[!NumberQ[budget], budget = 30];
lines = Select[ReadList[Environment["CASES"], String], StringLength[#] > 5 && StringTake[#, 5] === "CASE\t" &];
Do[Module[{parts = StringSplit[ln, "\t"], id, fam, f, t0, res, cls, det = ""},
   {id, fam} = parts[[2 ;; 3]]; f = ToExpression[parts[[4]]];
   t0 = AbsoluteTime[];
   res = TimeConstrained[Quiet[If[method === "Bronstein", BronsteinIntegrate[f, x], Integrate[f, x]]], budget, $Aborted];
   Which[
     res === $Aborted, cls = "TIMEOUT",
     ListQ[res] && First[res] === "not elementary", cls = "FALSE_CERT"; det = res[[2]],
     ListQ[res], cls = "DECLINE"; det = res[[2]],
     !FreeQ[res, Integrate], cls = "UNEVAL",
     TrueQ[Quiet[TimeConstrained[PossibleZeroQ[D[res, x] - f], 60, False]]], cls = "OK",
     True, cls = "WRONG"; det = ToString[InputForm[res]]];
   Print["RESULT\t", id, "\t", fam, "\t", cls, "\t", Round[AbsoluteTime[] - t0, 0.01], "\t", det]],
  {ln, lines}];
