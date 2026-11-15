(* tests/bronstein_bench/ideal_ops.m -- the ideal and field operations on the
   frozen inputs of capture.m (Phase 4.5): each call timed through the package
   routines (tragerOption["NativeFn"] = False) and through the native kernels
   (True), with the two results compared exactly.
   ./Mathilda -file tests/bronstein_bench/ideal_ops.m [maxSeconds per call] *)
Integrate[Sqrt[x], x, Method -> "Trager"];
ioRecs = Get["tests/bronstein_bench/ideal_inputs.m"];
ioOps = {"hnfModDet", "hnf", "idealTimes", "idealPower", "idealInverse", "idealReduce", "idealTimesGenerators", "principalGenerator", "normalizeAtInfinity"};
ioRun[r_, native_] := Module[{res},
  Trager`Private`tragerOption["NativeFn"] = native;
  res = AbsoluteTiming[TimeConstrained[Apply[Symbol["Trager`Private`" <> r["op"]], r["args"]], 120, "TIMEOUT"]];
  Trager`Private`tragerOption["NativeFn"] = False;
  res];
(* equality by value: over Q(theta)(vars) a canonical form is unique only up to
   the factored / expanded shape Together leaves (Trager's rc), so two results
   agree when every entry difference is 0; an ideal is compared part by part *)
ioZero[e_] := TrueQ[Quiet[Together[Expand[e]] === 0]];
ioSame[a_, b_] := Which[
  a === b, True,
  a === "TIMEOUT" || b === "TIMEOUT", "timeout",
  AssociationQ[a] && AssociationQ[b], Keys[a] === Keys[b] && AllTrue[Keys[a], ioSame[a[#], b[#]] === True &],
  ListQ[a] && ListQ[b], Dimensions[a] === Dimensions[b] && AllTrue[Flatten[a - b], ioZero],
  True, ioZero[a - b]];
Do[
  If[MemberQ[ioOps, r["op"]],
    With[{m = ioRun[r, False], n = ioRun[r, True]},
      Print["IOP\t", r["op"], "\t", r["field"], "\t", NumberForm[m[[1]], {6, 3}], "\t", NumberForm[n[[1]], {6, 3}], "\t",
        If[n[[1]] > 0, NumberForm[m[[1]]/n[[1]], {6, 1}], "-"], "\t", ioSame[m[[2]], n[[2]]]]]],
  {r, ioRecs}];
