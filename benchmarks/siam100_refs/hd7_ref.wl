(* Holdout-A D7 family (experiment 212), method B: island ends and kinks by
   FindRoot from the scan in hd7_ref.py, then NIntegrate (global adaptive
   Gauss-Kronrod) of the same 1-D forms at WorkingPrecision 100; plus the exact
   closed forms (Integrate / Gamma) for hd7.1, 2, 7, 8, 9 where WL finds one. *)
wp = 100;
fr[eq_, y0_] := y /. FindRoot[eq, {y, SetPrecision[y0, wp]}, WorkingPrecision -> wp, AccuracyGoal -> 95, PrecisionGoal -> 95];
ni[f_, rng_] := NIntegrate[f, rng, WorkingPrecision -> wp, PrecisionGoal -> 80, MaxRecursion -> 60];
w[g_, inv_, sc_][y_] := With[{s = Sqrt[-sc Log[g[y]]]}, inv[-y + s] - inv[-y - s]];
g32[y_] := y^2/32 + Exp[Sin[y]]; g16[y_] := y^2/16 + Exp[Sin[y]]; g100[y_] := y^2/100 + Exp[Sin[y]];
c3[t_] := Surd[t, 3]; c5[t_] := Surd[t, 5];
(* {g, inv, sc, island ends (starts), kink starts per island} *)
area[g_, inv_, sc_, isl_] := Total[Function[i, Module[{a, b, ks},
    a = If[i[[1]] === 0, 0, fr[g[y] == 1, i[[1]]]];
    b = If[i[[2]] === 0, 0, fr[g[y] == 1, i[[2]]]];
    ks = fr[-sc Log[g[y]] == y^2, #] & /@ i[[3]];
    ni[w[g, inv, sc][y], Join[{y, a}, ks, {b}]]]] /@ isl];
r = <||>;
r[1] = ni[2 (1 - x^4)^(1/4), {x, -1, 0, 1}];
a2 = x /. FindRoot[Cos[x] == x^2, {x, 8/10}, WorkingPrecision -> wp, AccuracyGoal -> 95];
r[2] = ni[Cos[x] - x^2, {x, -a2, a2}];
r[3] = area[g16, c3, 1, {{-2.576, 0, {-0.8}}}];
r[4] = area[g32, ArcSinh, 1, {{-2.85, 0, {}}, {3.77, 4.45, {}}}];
r[5] = area[g32, c5, 1, {{-2.85, 0, {-0.83}}, {3.77, 4.45, {}}}];
r[6] = area[g32, c3, 2, {{-2.85, 0, {-1.29}}, {3.77, 4.45, {}}}];
r[7] = ni[Pi - 2 ArcSin[1 - Sin[x]], {x, 0, Pi/2, Pi}];
r[8] = ni[Log[4 - Exp[x]] + x, {x, -Log[2 + Sqrt[3]], Log[2 + Sqrt[3]]}];
a9 = x /. FindRoot[3 - x^2 - Exp[x] == 0, {x, -17/10}, WorkingPrecision -> wp, AccuracyGoal -> 95];
b9 = x /. FindRoot[3 - x^2 - Exp[x] == 0, {x, 8/10}, WorkingPrecision -> wp, AccuracyGoal -> 95];
r[9] = ni[3 - x^2 - Exp[x], {x, a9, b9}];
r[10] = area[g100, c3, 1, {{-7.94, -7.04, {}}, {-3.04, 0, {-0.86}}, {3.25, 5.85, {}}}];
Do[Print["hd7.", k, "  ", ToString[N[r[k], 75], InputForm]], {k, 10}];
Print["closed1 ", N[4 Gamma[5/4]^2/Gamma[3/2], 75]];
Print["closed2 ", N[2 Sin[a2] - 2 a2^3/3, 75]];
Print["closed9 ", N[3 (b9 - a9) - (b9^3 - a9^3)/3 - (Exp[b9] - Exp[a9]), 75]];
Do[With[{f = {Pi - 2 ArcSin[1 - Sin[x]], Log[4 - Exp[x]] + x}[[j]],
         rg = {{x, 0, Pi}, {x, -Log[2 + Sqrt[3]], Log[2 + Sqrt[3]]}}[[j]]},
   c = TimeConstrained[Integrate[f, rg], 120, $Failed];
   Print["closed", {7, 8}[[j]], " ", InputForm[c], " ", If[NumericQ[c], N[c, 75], c]]], {j, 2}];
