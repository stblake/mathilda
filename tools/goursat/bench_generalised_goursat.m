(* tools/goursat/bench_generalised_goursat.m -- run one GeneralisedGoursat corpus.

     GG_CORPUS=stress_2026 GG_TIMEOUT=120 GG_METHOD=new ./Mathilda -file tools/goursat/bench_generalised_goursat.m

   Reads tests/data/generalised_goursat/$GG_CORPUS.m and writes
   tests/data/generalised_goursat/mathilda_$GG_CORPUS.tsv (GG_METHOD=new, the .m
   package) or c_$GG_CORPUS.tsv (GG_METHOD=c, the C GoursatAlgebraic, through
   Integrate[..., Method -> "GoursatAlgebraic"]), one row per case:

     label  status  seconds  verified  unevaluated  residual

   "seconds" times the reduction only (GoursatReduce / GoursatIntegrateAlgebraic),
   not the verification, matching tools/goursat/bench_sympy.py.  "verified" is the
   exact formal check (GoursatVerify / GoursatVerifyExpr).  Each case runs under
   TimeConstrained[GG_TIMEOUT] and MemoryConstrained[4 GB]. *)

LoadModule["goursat/GeneralisedGoursat.m"];

ggCorpus = Environment["GG_CORPUS"];
ggTimeout = ToExpression[Replace[Environment["GG_TIMEOUT"], $Failed -> "120"]];
ggMethod = Replace[Environment["GG_METHOD"], $Failed -> "new"];
ggDir = "tests/data/generalised_goursat/";
ggCases = Get[ggDir <> ggCorpus <> ".m"];

ggClean[s_String] := StringReplace[s, {"\t" -> " ", "\n" -> " "}];

ggRunNew[c_] := Module[{r, secs, ok, st, resid},
  {secs, r} = AbsoluteTiming[
    MemoryConstrained[TimeConstrained[
      If[c["Kind"] === "FRp",
        GeneralisedGoursat`GoursatReduce[c["F"], c["R"], t, c["p"]],
        GeneralisedGoursat`GoursatIntegrateAlgebraic[c["Integrand"], t]],
      ggTimeout, "TIMEOUT"], 4*10^9, "MEMORY"]];
  Which[
    r === "TIMEOUT" || r === "MEMORY", {r, secs, "", False, ""},
    !AssociationQ[r], {"DECLINE", secs, "", False, ""},
    True,
      ok = TimeConstrained[If[c["Kind"] === "FRp",
             GeneralisedGoursat`GoursatVerify[r, c["F"], t],
             GeneralisedGoursat`GoursatVerifyExpr[r, c["Integrand"], t]], ggTimeout, "TIMEOUT"];
      {r["Status"], secs, ToString[ok],
       !FreeQ[r["Antiderivative"], GeneralisedGoursat`Private`ggIntegral | Integrate],
       ToString[r["Residual"], InputForm]}]];

ggRunC[c_] := Module[{r, secs, ok},
  {secs, r} = AbsoluteTiming[TimeConstrained[
      Integrate[c["Integrand"], t, Method -> "GoursatAlgebraic"], ggTimeout, "TIMEOUT"]];
  Which[
    r === "TIMEOUT", {"TIMEOUT", secs, "", False, ""},
    !FreeQ[r, Integrate], {"DECLINE", secs, "", True, ""},
    True,
      ok = TimeConstrained[PossibleZeroQ[D[r, t] - c["Integrand"]], ggTimeout, "TIMEOUT"];
      {"elementary", secs, ToString[ok], False, ""}]];

ggOut = OpenWrite[ggDir <> If[ggMethod === "c", "c_", "mathilda_"] <> ggCorpus <> ".tsv"];
WriteString[ggOut, "label\tstatus\tseconds\tverified\tunevaluated\tresidual\n"];
Do[
  With[{row = If[ggMethod === "c", ggRunC[c], ggRunNew[c]]},
    WriteString[ggOut, ggClean[c["Label"]], "\t", row[[1]], "\t",
      ToString[NumberForm[row[[2]], {8, 4}]], "\t", row[[3]], "\t", ToString[row[[4]]], "\t",
      ggClean[StringTake[row[[5]], UpTo[200]]], "\n"]],
  {c, ggCases}];
Close[ggOut];
Print[ggCorpus, ": ", Length[ggCases], " cases"];
