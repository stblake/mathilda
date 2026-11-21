dead[m_] := StringJoin[Table["q" <> ToString[k] <> "[u, v, {1, 2, 3}, w + " <> ToString[k] <> "]; ", {k, m}]];
ToExpression["fM[x_] := Module[{a = x, b}, If[x > 0, b = a + 1, " <> dead[1000] <> "Null]; b]"];
t0 = TimeUsed[]; Do[fM[1], {2000}]; Print["Module us/call ", 10^6 (TimeUsed[] - t0)/2000];
