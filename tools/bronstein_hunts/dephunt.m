fams = {
 {"DI", {Log[u], Log[1 + I u], Log[u - I], Log[u + I], Log[1 - I u]}, {-1, 1, 2, I}},
 {"Dneg", {Log[u], Log[-u], Log[1 + u], Log[-1 - u], Log[2 u]}, {-1, 1, 2, 1/2}},
 {"Dprod", {Log[u - 1], Log[u + 1], Log[u^2 - 1], Log[u], Log[2 + 2 u]}, {-1, 1, 2}},
 {"Dexp", {Log[u], E^u, Log[u E^u], Log[1 + u], Log[(1 + u) E^(2 u)]}, {-1, 1, 2}}
};
Do[
  {nm, ks, cs} = fam; SeedRandom[4242 + StringLength[nm]];
  res = Table[With[{F = Sum[RandomChoice[cs] u^RandomInteger[{0, 2}] RandomChoice[ks]^RandomInteger[{1, 2}] RandomChoice[ks], {RandomInteger[{1, 3}]}]},
     {F, TimeConstrained[Quiet[Risch`ElementaryIntegralQ[Expand[D[F, u]], u]], 20, $TO]}], {80}];
  bad = Select[res, #[[2]] === False &];
  Print[nm, ": True=", Count[res, {_, True}], " False=", Length[bad], " TO=", Count[res, {_, $TO}],
     " undec=", Count[res, {_, _Risch`ElementaryIntegralQ}]];
  Scan[Print["   FALSE CERT: ", InputForm[#[[1]]]] &, bad],
  {fam, fams}]
