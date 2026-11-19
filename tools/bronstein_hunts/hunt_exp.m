ks = {Log[u], Log[1 + u], E^u, Log[1 + E^u], Log[u + a], ArcTan[E^u]}; cs = {-1, 1, 2, a};
SeedRandom[2026 + StringLength["Exp"]];
Do[With[{F = Sum[RandomChoice[cs] u^RandomInteger[{0, 2}] RandomChoice[ks]^RandomInteger[{1, 2}] RandomChoice[ks], {RandomInteger[{1, 3}]}]},
   Print[k, " ", Head[TimeConstrained[Quiet[Risch`ElementaryIntegralQ[Expand[D[F, u]], u]], 20, $TO]] /. Risch`ElementaryIntegralQ -> U, " ", InputForm[F]]], {k, 120}]
