(* RT_PRIM towers: D[F] is elementary over the Liouvillian tower of F, so False is a false certificate *)
ks = {Erf[u], Erfi[u], ExpIntegralEi[u], LogIntegral[u], Erf[u^2], ExpIntegralEi[2 u], E^(-u^2), E^u, Log[u], u, 1};
Do[SeedRandom[sd];
 res = Table[With[{F = Sum[RandomChoice[{-1, 1, 2, 1/2}] u^RandomInteger[{0, 2}] RandomChoice[ks] RandomChoice[ks]/RandomChoice[{1, u, 1 + u}], {RandomInteger[{1, 2}]}]},
   {F, TimeConstrained[Quiet[Risch`ElementaryIntegralQ[D[F, u], u]], 30, $TO]}], {80}];
 (* oracle: a False is a false certificate only when every Liouvillian kernel of F
    occurs in D[F] (2 ExpIntegralEi[u] has D = 2 E^u/u: False is right) *)
 lk = (Erf | Erfi | Erfc | ExpIntegralEi | LogIntegral)[_];
 res = Map[If[#[[2]] === False && !SubsetQ[Union[Cases[{D[#[[1]], u]}, lk, Infinity]], Union[Cases[{#[[1]]}, lk, Infinity]]],
    {#[[1]], "FalseOK"}, #] &, res];
 Print["seed ", sd, ": True=", Count[res, {_, True}], " False=", Count[res, {_, False}], " FalseOK=", Count[res, {_, "FalseOK"}], " TO=", Count[res, {_, $TO}], " undec=", Count[res, {_, _Risch`ElementaryIntegralQ}]];
 Scan[Print["  FALSE CERT: ", InputForm[#[[1]]]] &, Select[res, #[[2]] === False &]], {sd, {901, 902}}];
