(* run.m -- integrate the cases of the file CASES (lines CASE<TAB>id<TAB>fam<TAB>f)
   with BronsteinIntegrate (METHOD=Bronstein) or plain Integrate (METHOD=Automatic),
   and print one RESULT line per case:
     RESULT<TAB>id<TAB>fam<TAB>class<TAB>seconds<TAB>detail
   class: OK | SLOWCHECK (verified numerically, PossibleZeroQ out of time) |
          WRONG | FALSE_CERT | DECLINE | TIMEOUT | UNEVAL.            *)
LoadModule["bronstein/Bronstein.m"];
method = Environment["METHOD"]; If[!StringQ[method], method = "Bronstein"];
budget = ToExpression[Environment["BUDGET"]]; If[!NumberQ[budget], budget = 30];
lines = Select[ReadList[Environment["CASES"], String], StringLength[#] > 5 && StringTake[#, 5] === "CASE\t" &];
Do[Module[{parts = StringSplit[ln, "\t"], id, fam, f, t0, res, cls, det = ""},
   {id, fam} = parts[[2 ;; 3]]; f = ToExpression[parts[[4]]];   (* parts[[5]], when present, is the source F *)
   t0 = AbsoluteTime[];
   res = TimeConstrained[Quiet[If[method === "Bronstein", BronsteinIntegrate[f, x], Integrate[f, x]]], budget, $Aborted];
   Which[
     res === $Aborted, cls = "TIMEOUT",
     ListQ[res] && First[res] === "not elementary", cls = "FALSE_CERT"; det = res[[2]],
     ListQ[res], cls = "DECLINE"; det = res[[2]],
     !FreeQ[res, Integrate], cls = "UNEVAL",
     (pz = Quiet[TimeConstrained[PossibleZeroQ[D[res, x] - f], 60, $Aborted]]) === True, cls = "OK",
     (* PossibleZeroQ out of time on a large answer: a 30-digit check at three
        points separates a slow verification (SLOWCHECK) from a wrong answer *)
     pz === $Aborted && AllTrue[{7/10, 13/10, 23/10},
        Abs[Quiet[N[(D[res, x] - f) /. x -> #, 30]]] < 10^-20 &], cls = "SLOWCHECK",
     True, cls = "WRONG"; det = ToString[InputForm[res]]];
   Print["RESULT\t", id, "\t", fam, "\t", cls, "\t", Round[AbsoluteTime[] - t0, 0.01], "\t", det]],
  {ln, lines}];
