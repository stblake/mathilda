(* ::Package:: *)

(* ========================================================================== *)
(* tools/trager_profile.m -- where does Integrate[f, x, Method -> "Trager"]    *)
(* spend its time?  A per-function timer over the Trager package's stages.     *)
(*                                                                            *)
(*   TRAGER_PROFILE_F='1/(x + 1) (x^3 - 2 x)^(1/4)' \                         *)
(*     ./Mathilda -file tools/trager_profile.m                                *)
(*                                                                            *)
(* Optional environment:                                                      *)
(*   TRAGER_PROFILE_TIME    TimeConstrained limit in seconds (default 300)     *)
(*   TRAGER_PROFILE_NATIVE  "0" runs with the FLINT` kernels switched off      *)
(*                          (tragerOption["NativeKernels"] = False): A/B      *)
(*                                                                            *)
(* Each listed function gets a DownValue wrapper that adds the wall time of   *)
(* its outermost active call (recursion counted once) to a total.  Totals are *)
(* INCLUSIVE: a stage's time contains the stages it calls, so the table reads *)
(* top-down.  The wrappers cost a few microseconds a call, which matters only *)
(* for the very frequent leaf functions (kc, gfMod).                          *)
(* ========================================================================== *)

trProfFunctions = {
  "integrateComponent", "hermiteReduce", "componentResidues", "residueField",
  "residueFieldTower", "factorOverField", "tragerLogPart", "tragerBasisLogPart",
  "tragerOrbitLogTerms", "tragerLogTerms", "torsionOrderSpec", "orderModPSpec",
  "torsionScreen", "traceCertificate", "traceDivisorTest", "factorIdeal", "generatorTerms", "principalGenerator", "logSumResidual",
  "jacAdd", "jacReduce", "jacMultiple", "reducedRows", "elementInverse",
  "idealTimes", "idealTimesGenerators", "nativeIdeal", "idealReduce",
  "divisorIdeal", "specIdealX", "hnf", "hnfModDet", "kGCD", "kExtGCD",
  "kDivRem", "afTimesPoly", "kc", "gfMod", "rc", "presentable", "spotCheck"};

(* load the package (lazily loaded by the first Trager call) *)
Quiet[Integrate[Sqrt[x^3 + 1]/(x - 2), x, Method -> "Trager"]];

trProfTotal[_] = 0.;
trProfCount[_] = 0;
Do[
  With[{s = ToExpression["Trager`Private`" <> nm], nm = nm},
    If[DownValues[s] =!= {},
      PrependTo[DownValues[s],
        HoldPattern[s[args___]] /; !TrueQ[trProfInside[nm]] :>
          Block[{trProfInside}, trProfInside[nm] = True;
            Module[{t0, r},
              {t0, r} = AbsoluteTiming[s[args]];
              trProfTotal[nm] = trProfTotal[nm] + t0;
              trProfCount[nm] = trProfCount[nm] + 1;
              r]]]]],
  {nm, trProfFunctions}];

trProfF = Environment["TRAGER_PROFILE_F"];
If[!StringQ[trProfF],
  Print["usage: TRAGER_PROFILE_F='<integrand in x>' Mathilda -file tools/trager_profile.m"],
  Module[{f = ToExpression[trProfF], lim, tt, rr, rows},
    lim = If[StringQ[Environment["TRAGER_PROFILE_TIME"]], ToExpression[Environment["TRAGER_PROFILE_TIME"]], 300];
    If[Environment["TRAGER_PROFILE_NATIVE"] === "0", Trager`Private`tragerOption["NativeKernels"] = False];
    {tt, rr} = AbsoluteTiming[TimeConstrained[Trager`TragerIntegrate[f, x], lim, $Aborted]];
    Print["integrand: ", InputForm[f]];
    Print["result:    ", If[rr === $Aborted, "TIMEOUT after " <> ToString[lim] <> " s",
      If[ListQ[rr], rr, "antiderivative (LeafCount " <> ToString[LeafCount[rr]] <> ")"]]];
    Print["total:     ", tt, " s   (native kernels: ", Trager`Private`tragerOption["NativeKernels"], ")"];
    rows = Reverse[SortBy[Select[{#, trProfTotal[#], trProfCount[#]} & /@ trProfFunctions, #[[3]] > 0 &], #[[2]] &]];
    Print["function                     seconds     calls"];
    Do[Print[StringPadRight[r[[1]], 28], StringPadLeft[ToString[NumberForm[r[[2]], {8, 3}]], 10],
             StringPadLeft[ToString[r[[3]]], 10]], {r, rows}]]];
