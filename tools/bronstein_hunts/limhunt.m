(* LimitedIntegrate below a logarithm (campaign step 4): a primitive polynomial
   p_1 t + p_0, t = Log[v], over an exponential or hyperbolic tower, whose
   theta^1 coefficient integrates to a new logarithm.  F random with products
   of Log[v] and functions of E^w / Tanh[w] / Tan[w]; D[F] is elementary, so a
   False is a false certificate.  Second half: random integrands Log[v] g;
   a False contradicted by a verified antiderivative is reported. *)
lg = {Log[u], Log[1 + u], Log[u^2 + 1], Log[1 + E^u], Log[u + E^u]};
ek = {E^u, E^(2 u), E^(u^2), Tanh[u], Tanh[2 u], Tan[u], Coth[u], 1/(1 + E^u), E^u/(1 + E^(2 u))};
Do[SeedRandom[sd];
 res = Table[Module[{G = RandomChoice[ek], P = RandomChoice[lg], F},
    F = Sum[RandomChoice[{-1, 1, 2, 1/2}] u^RandomInteger[{0, 1}] RandomChoice[{P, P^2, 1}] RandomChoice[{G, G^2, Log[1 + G^2], Log[1 + G], 1/(1 + G)}]/RandomChoice[{1, u, 1 + u, 1 + G^2}], {RandomInteger[{1, 2}]}];
    {F, TimeConstrained[Quiet[Risch`ElementaryIntegralQ[D[F, u], u]], 30, $TO]}], {80}];
 Print["seed ", sd, ": True=", Count[res, {_, True}], " False=", Count[res, {_, False}], " TO=", Count[res, {_, $TO}], " undec=", Count[res, {_, _Risch`ElementaryIntegralQ}]];
 Scan[Print["  FALSE CERT: ", InputForm[#[[1]]]] &, Select[res, #[[2]] === False &]], {sd, {1801, 1802}}];
SeedRandom[1810];
res2 = Table[Module[{G = RandomChoice[ek], P = RandomChoice[lg], f},
    f = RandomChoice[{1, u, 1/u, 1/(1 + u)}] RandomChoice[{P, P^2}] RandomChoice[{G, G^2, 1/(1 + G), D[G, u], 1 - G^2}];
    {f, TimeConstrained[Quiet[Risch`ElementaryIntegralQ[f, u]], 30, $TO], TimeConstrained[Quiet[Integrate`RischTranscendental[f, u]], 30, $TO]}], {80}];
Print["random: True=", Count[res2, {_, True, _}], " False=", Count[res2, {_, False, _}], " TO=", Count[res2, {_, $TO, _}], " undec=", Count[res2, {_, _Risch`ElementaryIntegralQ, _}]];
Scan[If[#[[2]] === False && Head[#[[3]]] =!= Integrate`RischTranscendental && #[[3]] =!= $TO && PossibleZeroQ[D[#[[3]], u] - #[[1]]], Print["  CONTRADICTION: ", InputForm[#[[1]]]]] &, res2];
Scan[Print["  False: ", InputForm[#[[1]]]] &, Select[res2, #[[2]] === False &]];
