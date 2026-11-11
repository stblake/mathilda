fs = Get["tasks/polygcd_together_deficiencies/d3_stress_inputs.m"];
pts = {1/3, 7/5, -2/7};
ops = {{"Together", Together}, {"TogetherExt", Together[#, Extension -> Automatic] &}, {"Cancel", Cancel}, {"CancelExt", Cancel[#, Extension -> Automatic] &}, {"Simplify", Simplify}};
res = <||>;
Do[f = fs[[i]];
  Do[r = AbsoluteTiming[TimeConstrained[op[[2]][f], 10, $Aborted]];
    st = If[r[[2]] === $Aborted, "TIMEOUT",
      If[And @@ Table[With[{a = N[f /. x -> p, 30], b = N[r[[2]] /. x -> p, 30]}, !NumberQ[a] || Abs[a] > 10^6 || Abs[a - b] < 10^-15 (1 + Abs[a])], {p, pts}], "OK", "WRONG"]];
    res[{op[[1]], st}] = Lookup[res, Key[{op[[1]], st}], 0] + 1;
    If[st === "WRONG", Print["WRONG ", i, " ", op[[1]], " ", InputForm[f]]],
   {op, ops}], {i, Length[fs]}];
Print[KeySort[res]];
