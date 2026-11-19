(* RDE false "NoSolution" hunt: y known in the tower, g = D[y] + f y *)
fams = {
 {"Q", {Log[x], Log[1 + x], E^x, x, 1/(1 + x)}, {1, 2, x, 1/x, 1/(1 + x), 2 x}},
 {"QI", {Log[x], Log[x + I], E^(I x), x, 1/(x - I)}, {I, 2 I/(1 + x^2), x, 1/x, I x}},
 {"Exp", {E^x, E^(x^2), Log[x], x, 1/(1 + E^x)}, {1, 2 x, E^x, x, 1 + 1/x}},
 {"Qa", {Log[x], Log[x + a], E^(a x), x}, {a, 1/x, a/(x + a), x}}};
Do[{nm, ks, fs} = fam; SeedRandom[4040 + StringLength[nm]];
 res = Table[Module[{y, f, g, r},
    y = RandomChoice[{1, 2, x, a}] RandomChoice[ks]^RandomInteger[{1, 2}] RandomChoice[ks] + RandomChoice[{0, 1, x}] RandomChoice[ks];
    f = RandomChoice[fs] + RandomChoice[{0, 0, D[RandomChoice[ks], x]}];
    g = Together[D[y, x] + f y];
    (* the solution field is the tower of y's kernels too: without them a y in
       Q(x, Log[1 + x]) is rightly "NoSolution" over the field of f and g *)
    r = TimeConstrained[Risch`RischDECertified[f, g, x, Cases[{y}, Log[_] | Power[E, _], Infinity], If[FreeQ[{y, f}, a], {}, {a}]], 30, $TO];
    {y, f, r}], {80}];
 Print[nm, ": NoSolution(false)=", Count[res, {_, _, "NoSolution"}], " solved=", Count[res, {_, _, r_ /; ! StringQ[r] && r =!= $Failed && r =!= $TO}], " failed=", Count[res, {_, _, $Failed}], " TO=", Count[res, {_, _, $TO}]];
 Scan[Print["  FALSE: y=", InputForm[#[[1]]], " f=", InputForm[#[[2]]]] &, Select[res, #[[3]] === "NoSolution" &]], {fam, fams}];
