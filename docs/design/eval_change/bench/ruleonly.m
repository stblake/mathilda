dead[m_] := StringJoin[Table["q" <> ToString[k] <> "[u, v, {1, 2, 3}, w + " <> ToString[k] <> "]; ", {k, m}]];
ToExpression["fP[x_] := If[x > 0, x + 1, " <> dead[1000] <> "Null]"];
Print[AbsoluteTiming[Do[fP[1], {30000}]][[1]]/30000*10^6, " us/call"];
