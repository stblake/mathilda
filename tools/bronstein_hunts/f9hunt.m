SeedRandom[909];
pool = {Sqrt[2] u, 2 u, u, Sqrt[3] u^2, u^2, Pi u, Sqrt[2] + u, 1/(1 + u), Sqrt[2]/(1 + u), a u, Sqrt[2] a u, (1 + Sqrt[2]) u, Log[2] u};
bad = 0; n = 0;
Do[gs = RandomSample[pool, RandomInteger[{1, 4}]];
   rs = RandomChoice[{-2, -1, 0, 1, 1/2, 3}, Length[gs]];
   th = Expand[rs . gs]; r = Risch`RationalSpan[th, gs, {u}]; n++;
   If[r === False, bad++; Print["FALSE: ", InputForm[{th, gs}]]], {300}];
Print["span: ", bad, " false of ", n];
ks = {E^(Sqrt[2] u), E^(2 u), E^u, u^Sqrt[2], Log[u], E^((1 + Sqrt[2]) u), Log[1 + u], E^(Sqrt[3] u)};
SeedRandom[910];
res = Table[With[{F = Sum[RandomChoice[{-1, 1, 2, Sqrt[2], a}] u^RandomInteger[{0, 2}] RandomChoice[ks]^RandomInteger[{1, 2}] RandomChoice[ks]/RandomChoice[{1, 1 + u, u}], {RandomInteger[{1, 3}]}]},
   {F, TimeConstrained[Quiet[Risch`ElementaryIntegralQ[Expand[D[F, u]], u]], 20, $TO]}], {150}];
Print["EIQ: True=", Count[res, {_, True}], " False=", Count[res, {_, False}], " TO=", Count[res, {_, $TO}]];
Scan[Print["  FALSE CERT: ", InputForm[#[[1]]]] &, Select[res, #[[2]] === False &]];
