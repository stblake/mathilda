(* Reference for d7, second method: WL NIntegrate (Gauss-Kronrod/DoubleExponential
   chosen by NIntegrate) of the x-interval length at WorkingPrecision 100, band
   ends by FindRoot at WorkingPrecision 100. *)
g[y_] := y^2/32 + Exp[Sin[y]];
w[y_?NumericQ] := Surd[-y + Sqrt[-Log[g[y]]], 3] - Surd[-y - Sqrt[-Log[g[y]]], 3];
wp = 100;
y1 = y /. FindRoot[g[y] == 1, {y, -285/100}, WorkingPrecision -> wp];
y2 = y /. FindRoot[g[y] == 1, {y, 377/100}, WorkingPrecision -> wp];
y3 = y /. FindRoot[g[y] == 1, {y, 445/100}, WorkingPrecision -> wp];
yc = y /. FindRoot[-Log[g[y]] == y^2, {y, -83/100}, WorkingPrecision -> wp];
a = NIntegrate[w[y], {y, y1, yc, 0}, WorkingPrecision -> wp, PrecisionGoal -> 80, MaxRecursion -> 200] +
    NIntegrate[w[y], {y, y2, y3}, WorkingPrecision -> wp, PrecisionGoal -> 80, MaxRecursion -> 200];
Print[N[a, 80]];
