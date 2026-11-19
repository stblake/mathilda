ks = {E^u, E^(u^2), E^(u + u^2), E^(2 u + u^2), E^(u - u^2), E^(a u), E^(a u + u^2), Log[u], Log[1 + E^u], E^(u^2 + u^3), E^(u^3)};
Do[SeedRandom[sd];
 res = Table[With[{F = Sum[RandomChoice[{-1, 1, 2, a, 1/2}] u^RandomInteger[{0, 2}] RandomChoice[ks]^RandomInteger[{1, 2}] RandomChoice[ks]/RandomChoice[{1, 1 + u, u, 1 + E^u, 1 + E^(u^2)}], {RandomInteger[{1, 3}]}]},
   {F, TimeConstrained[Quiet[Risch`ElementaryIntegralQ[Expand[D[F, u]], u]], 20, $TO]}], {120}];
 Print["seed ", sd, ": True=", Count[res, {_, True}], " False=", Count[res, {_, False}], " TO=", Count[res, {_, $TO}]];
 Scan[Print["  FALSE CERT: ", InputForm[#[[1]]]] &, Select[res, #[[2]] === False &]], {sd, {1201, 1202}}];
Scan[Print[InputForm[#], " -> ", Quiet[Risch`ElementaryIntegralQ[#, u]]] &, {E^(u + u^2)/u + E^u, E^(2 u + u^2) + E^(u^2), E^(u + u^2) (1 + 2 u) + E^u/u, E^(u + u^2)/(1 + u)}];
