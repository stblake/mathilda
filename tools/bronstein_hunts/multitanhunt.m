(* Several tangents (campaign step 6).  F random with TWO tangent kernels,
   independent (Tan[u], Tan[u^2]) or dependent (Tan[u], Tan[2 u], since Tan[2 u]
   is rational in Tan[u]); D[F] is elementary, so a False is a false
   certificate. *)
pairs = {{Tan[u], Tan[u^2]}, {Tan[u], Tan[2 u]}, {Tan[u], Tan[u + 1]}, {Tan[u], Tan[u/2]}, {Tan[u^2], Tan[u^3]}};
mk[{A_, B_}] := {A, B, A B, 1/(1 + A^2), 1/(1 + B), Log[1 + A^2], Log[A + B], E^A, u, 1};
Do[SeedRandom[sd];
 res = Table[Module[{P = RandomChoice[pairs], ks, F},
    ks = mk[P];
    F = Sum[RandomChoice[{-1, 1, 2, 1/2}] u^RandomInteger[{0, 1}] RandomChoice[ks] RandomChoice[ks]/RandomChoice[{1, u, 1 + u, 1 + P[[1]]}], {RandomInteger[{1, 2}]}];
    {F, P, TimeConstrained[Quiet[Risch`ElementaryIntegralQ[D[F, u], u]], 30, $TO]}], {80}];
 Print["seed ", sd, ": True=", Count[res, {_, _, True}], " False=", Count[res, {_, _, False}], " TO=", Count[res, {_, _, $TO}], " undec=", Count[res, {_, _, _Risch`ElementaryIntegralQ}]];
 Print["  dependent pairs (Tan[2 u], Tan[u + 1], Tan[u/2]) undec=", Count[res, {_, P_ /; !MatchQ[P, {Tan[u], Tan[u^2]} | {Tan[u^2], Tan[u^3]}], _Risch`ElementaryIntegralQ}]];
 Scan[Print["  FALSE CERT: ", InputForm[#[[1]]]] &, Select[res, #[[3]] === False &]], {sd, {1601, 1602}}];
