fams = {
 {"Q", {Log[u], Log[1 + u], Log[u - 2], ArcTan[u], Log[1 + u^2], ArcTan[2 u], Log[u^2 + u + 1], Log[3 + u]}, {-2, -1, 1, 2, 1/2}},
 {"QI", {Log[u], Log[u + I], Log[u - I], Log[1 + I u], Log[1 + u], ArcTan[u]}, {-1, 1, 2, I, -I}},
 {"Qa", {Log[u], Log[u + a], Log[1 + u], Log[u - b], Log[a u + 1], ArcTan[u], Log[u^2 + a], Log[a + b u]}, {-2, -1, 1, 2, a, 1/a, b}},
 {"Exp", {Log[u], Log[1 + u], E^u, Log[1 + E^u], Log[u + a], ArcTan[E^u]}, {-1, 1, 2, a}},
 {"Nest", {Log[u], Log[Log[u]], Log[1 + Log[u]], Log[1 + u], Log[u + a]}, {-1, 1, 2, a}}
};
Do[
  {nm, ks, cs} = fam; SeedRandom[2026 + StringLength[nm]];
  res = Table[With[{F = Sum[RandomChoice[cs] u^RandomInteger[{0, 2}] RandomChoice[ks]^RandomInteger[{1, 2}] RandomChoice[ks], {RandomInteger[{1, 3}]}]},
     {F, TimeConstrained[Quiet[Risch`ElementaryIntegralQ[Expand[D[F, u]], u]], 20, $TO]}], {120}];
  bad = Select[res, #[[2]] === False &];
  Print[nm, ": True=", Count[res, {_, True}], " False=", Length[bad], " TO=", Count[res, {_, $TO}],
     " undec=", Count[res, {_, _Risch`ElementaryIntegralQ}]];
  Scan[Print["   FALSE CERT: ", InputForm[#[[1]]]] &, bad],
  {fam, fams}]
