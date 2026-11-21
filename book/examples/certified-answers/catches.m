# A near miss: an identity that fails by 10^-25 x.
e = (1 + x)^20 - Expand[(1 + x)^20] + 10^-25 x;
PossibleZeroQ[e]
Quiet[Certify[(1 + x)^20 == Expand[(1 + x)^20] - 10^-25 x, "Engine" -> {"Mathilda", "Exact"}]]["Engines"]
