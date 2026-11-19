(* Tangents over k(i) (campaign step 4: t^2 + 1 = (t - I)(t + I) splits).  F
   random in {u, Tan[v], Log[u + I], E^(I u), ...} with Gaussian coefficients;
   D[F] is elementary, so a False is a false certificate. *)
tk = {Tan[u], Tan[u^2], Tan[2 u], Tan[u + 1]};
mk[T_] := {Log[u + I], Log[T - I], Log[T + I], 1/(T - I), 1/(T + I), E^(I u), E^T, Log[u], 1, T, T^2, 1/(1 + T^2)};
Do[SeedRandom[sd];
 res = Table[Module[{T = RandomChoice[tk], ks, F},
    ks = mk[T];
    F = Sum[RandomChoice[{-1, 1, I, -I, 1 + I, 1/2}] u^RandomInteger[{0, 1}] RandomChoice[ks] RandomChoice[ks]/RandomChoice[{1, u, u + I, 1 + T, T - I}], {RandomInteger[{1, 2}]}];
    {F, TimeConstrained[Quiet[Risch`ElementaryIntegralQ[D[F, u], u]], 30, $TO]}], {80}];
 Print["seed ", sd, ": True=", Count[res, {_, True}], " False=", Count[res, {_, False}], " TO=", Count[res, {_, $TO}], " undec=", Count[res, {_, _Risch`ElementaryIntegralQ}]];
 Scan[Print["  FALSE CERT: ", InputForm[#[[1]]]] &, Select[res, #[[2]] === False &]], {sd, {1501, 1502}}];
SeedRandom[1510];
res2 = Table[Module[{T = RandomChoice[tk], ks, f},
    ks = mk[T];
    f = RandomChoice[{1, u, I, T - I, 1/(u + I)}] RandomChoice[ks] RandomChoice[{1, T, I u, T + I}];
    {f, TimeConstrained[Quiet[Risch`ElementaryIntegralQ[f, u]], 30, $TO], TimeConstrained[Quiet[Integrate`RischTranscendental[f, u]], 30, $TO]}], {80}];
Print["random: True=", Count[res2, {_, True, _}], " False=", Count[res2, {_, False, _}], " TO=", Count[res2, {_, $TO, _}], " undec=", Count[res2, {_, _Risch`ElementaryIntegralQ, _}]];
Scan[If[#[[2]] === False && Head[#[[3]]] =!= Integrate`RischTranscendental && #[[3]] =!= $TO && PossibleZeroQ[D[#[[3]], u] - #[[1]]], Print["  CONTRADICTION: ", InputForm[#[[1]]]]] &, res2];
