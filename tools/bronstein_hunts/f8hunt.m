(* F8: false-certificate hunt over hypertangent towers.  F random in the tower
   {u, Tan[v], E^w, Log[..]}; D[F] is elementary, so a False is a false
   certificate.  Second half: random integrands (ground truth unknown) whose
   False verdicts are re-checked by a cheap numeric sanity: none required. *)
tk = {Tan[u], Tan[u^2], Tan[u]/u, Tan[2 u], Tan[u + 1]};
mk[T_] := {E^T, E^(T^2), E^(u T), E^(T/u), Log[1 + T^2], Log[T], Log[u + T], Log[u], E^u, 1, T, T^2, 1/(1 + T), 1/(1 + T^2), 1/T};
Do[SeedRandom[sd];
 res = Table[Module[{T = RandomChoice[tk], ks, F},
    ks = mk[T];
    F = Sum[RandomChoice[{-1, 1, 2, 1/2, a}] u^RandomInteger[{0, 1}] RandomChoice[ks] RandomChoice[ks]/RandomChoice[{1, u, 1 + u, 1 + T, 1 + T^2, T}], {RandomInteger[{1, 2}]}];
    {F, TimeConstrained[Quiet[Risch`ElementaryIntegralQ[D[F, u], u]], 30, $TO]}], {100}];
 Print["seed ", sd, ": True=", Count[res, {_, True}], " False=", Count[res, {_, False}], " TO=", Count[res, {_, $TO}], " undec=", Count[res, {_, _Risch`ElementaryIntegralQ}]];
 Scan[Print["  FALSE CERT: ", InputForm[#[[1]]]] &, Select[res, #[[2]] === False &]], {sd, {801, 802, 803}}];
(* non-elementary side: random integrands; count decisions *)
SeedRandom[810];
res2 = Table[Module[{T = RandomChoice[tk], ks, f},
    ks = mk[T];
    f = RandomChoice[{1, u, T, 1 + T^2, 1/u, u T}] RandomChoice[ks] RandomChoice[{1, T, u, 1 + T^2}];
    {f, TimeConstrained[Quiet[Risch`ElementaryIntegralQ[f, u]], 30, $TO], TimeConstrained[Quiet[Integrate`RischTranscendental[f, u]], 30, $TO]}], {100}];
Print["random: True=", Count[res2, {_, True, _}], " False=", Count[res2, {_, False, _}], " TO=", Count[res2, {_, $TO, _}]];
(* consistency: an exhibited verified antiderivative with a False verdict is a false certificate *)
Scan[If[#[[2]] === False && Head[#[[3]]] =!= Integrate`RischTranscendental && #[[3]] =!= $TO && PossibleZeroQ[D[#[[3]], u] - #[[1]]], Print["  CONTRADICTION: ", InputForm[#[[1]]]]] &, res2];
Scan[Print["  False: ", InputForm[#[[1]]]] &, Select[res2, #[[2]] === False &]];
