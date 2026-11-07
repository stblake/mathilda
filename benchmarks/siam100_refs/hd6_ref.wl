(* Holdout-A D6 family (experiment 211), method B: NIntegrate over the whole
   interval (no splitting), global adaptive Gauss-Kronrod, WorkingPrecision 90;
   plus the exact Integrate closed forms where WL has one (hd6.4, 7, 9, 10). *)
fs = {
  {Cos[Exp[x^2]], 0, 2},
  {Exp[x + Cos[Exp[Exp[x]]]], -1, 1},
  {Sin[Exp[Exp[2 x]]], 0, 1},
  {Sin[x^3], 0, 3},
  {Exp[Sin[Exp[x]]], 0, Pi},
  {Sin[Exp[4 x]]/(2 + x), -1, 1},
  {x Exp[x^2] Sin[Exp[x^2]], 0, 2},
  {Exp[x + Sin[Exp[Exp[x + 1/2]]]], -1, 1},
  {Cos[x^2 + x], 0, 8},
  {Cos[1/(x + 1/100)], 0, 1}};
Do[With[{f = fs[[k, 1]], a = fs[[k, 2]], b = fs[[k, 3]]},
  v = NIntegrate[f, {x, a, b}, WorkingPrecision -> 90, PrecisionGoal -> 70,
        MaxRecursion -> 60, Method -> "GaussKronrodRule"];
  Print["hd6.", k, "  ", NumberForm[v, 70]];
  If[MemberQ[{4, 7, 9, 10}, k],
    c = Integrate[f, {x, a, b}];
    Print["   closed ", InputForm[c], "  ", NumberForm[N[c, 70], 70]]]], {k, 10}];
