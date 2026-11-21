# Ch.7 -- non-local exits: Return, Throw/Catch, Abort/CheckAbort/AbortProtect
sign[x_] := If[x < 0, Return["negative"], "non-negative"]
{sign[-3], sign[4]}
firstPrime[l_] := Do[If[PrimeQ[x], Return[x]], {x, l}]
firstPrime[{4, 6, 9, 11, 13}]
Catch[Do[If[i^2 > 50, Throw[i]], {i, 100}]]
CheckAbort[Abort[]; "never reached", "aborted"]
{1, CheckAbort[Abort[], 2], 3}
CheckAbort[AbortProtect[Abort[]; done = True; 7], "abort delivered at the end"]
done
