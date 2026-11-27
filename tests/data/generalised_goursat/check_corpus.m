(* tests/data/generalised_goursat/check_corpus.m -- judge one GeneralisedGoursat corpus.

   ggCheckCorpus[file, timeout] runs every case of the corpus file (a list of
   <|"Label", "Kind" -> "FRp" | "expr", "F", "R", "p", "Integrand", "Truth"|>,
   written by tools/goursat/export_corpora.py) and returns
     {ok, {{label, verdict, status, seconds}, ...}}
   listing every case that is not OK.  Verdicts (as stress_goursat.py judges):
     OK       the verdict agrees with the truth (or there is none) and the exact
              verifier accepts the answer
     MISSED   "undecided" although the truth is known (incompleteness, not a bug)
     WRONG    the verdict contradicts the truth, or the answer does not verify
     UNEVAL   "elementary" with an unevaluated piece
     DECLINE  $Failed (out of shape, or an internal failure)
     TIMEOUT  over the per-case budget
   Used by tests/test_generalised_goursat_stress.c. *)

ggJudge[c_, timeout_] := Module[{r, secs, ok, st, truth = c["Truth"]},
  {secs, r} = AbsoluteTiming[MemoryConstrained[TimeConstrained[
    If[c["Kind"] === "FRp",
      GeneralisedGoursat`GoursatReduce[c["F"], c["R"], t, c["p"]],
      GeneralisedGoursat`GoursatIntegrateAlgebraic[c["Integrand"], t]],
    timeout, "TIMEOUT"], 4*10^9, "MEMORY"]];
  Which[
    r === "TIMEOUT" || r === "MEMORY", {c["Label"], r, r, secs},
    !AssociationQ[r], {c["Label"], "DECLINE", "$Failed", secs},
    True,
      st = r["Status"];
      ok = TimeConstrained[If[c["Kind"] === "FRp",
        GeneralisedGoursat`GoursatVerify[r, c["F"], t],
        GeneralisedGoursat`GoursatVerifyExpr[r, c["Integrand"], t]], timeout, "TIMEOUT"];
      {c["Label"],
       Which[
         st === "elementary" && !FreeQ[r["Antiderivative"], GeneralisedGoursat`Private`ggIntegral],
           "UNEVAL",
         ok =!= True, "WRONG",
         truth === None || st === truth, "OK",
         st === "undecided", "MISSED",
         True, "WRONG"],
       st, secs}]];

ggCheckCorpus[file_String, timeout_] := Module[{cases = Get[file], rows},
  rows = ggJudge[#, timeout] & /@ cases;
  {Count[rows, {_, "OK", _, _}], Select[rows, #[[2]] =!= "OK" &]}];
