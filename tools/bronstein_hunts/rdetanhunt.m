(* The Risch DE over a split tangent level (campaign step 4): Tanh, and Tan
   over k(i), solved by the Moebius change to the exponential theta
   (rde_tower_mobius).  y random in the tower, g = D[y] + f y: a solution exists
   by construction, so "NoSolution" is a false certificate, and an answer must
   be y (the solution is unique when f is not a log-derivative). *)
tk = {Tanh[x], Tanh[x^2], Tanh[1 + x], Tanh[x]/x};
ti = {Tan[x], Tan[1 + x], Tan[2 x]};
mky[T_, c_] := RandomChoice[{1, x, c, x c}] RandomChoice[{1, T, T^2, 1/(1 + T), 1/(T - 1), 1/(T + c), T/(x + T), E^T, Log[x]}] / RandomChoice[{1, x, 1 + x, 1 + T^2, T}];
Do[SeedRandom[sd];
 res = Table[Module[{T, c, y, f, g, r},
    If[k == "tanh", T = RandomChoice[tk]; c = 1, T = RandomChoice[ti]; c = I];
    y = mky[T, c]; f = RandomChoice[{1, x, c, T, 1 - T^2, 1/x, x T}];
    g = Together[D[y, x] + f y];
    r = TimeConstrained[Quiet[Risch`RischDECertified[f, g, x]], 30, $TO];
    Which[r === "NoSolution", {"FALSE", y, f}, r === $TO, {"TO"}, r === $Failed, {"failed"},
      Head[r] === Risch`RischDECertified, {"undec"},
      PossibleZeroQ[Together[r - y]] || PossibleZeroQ[Together[D[r, x] + f r - g]], {"solved"}, True, {"WRONG", y, f, r}]],
   {60}];
 Print[k, " seed ", sd, ": ", Tally[First /@ res]];
 Scan[Print["  ", #[[1]], ": y = ", InputForm[#[[2]]], ", f = ", InputForm[#[[3]]]] &, Select[res, MatchQ[#[[1]], "FALSE" | "WRONG"] &]],
 {k, {"tanh", "tan-i"}}, {sd, {1801, 1802}}];
