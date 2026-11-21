# Antiderivatives: the cleared derivative identity.
SetOptions[Certify, "Verify" -> False];
c = Certify[Integrate[1/(x^2 + 2), x]];
c["Result"]
c["LeanTheorem"]
Length[c["TrustedSteps"]]
Certify[Integrate[(x^4 - 3 x^2 + 6)/(x^6 - 5 x^4 + 5 x^2 + 4), x], "Engine" -> "Exact", "Verify" -> True]["Verified"]
Certify[Integrate[x, x] == x^2/2 + 7, "Engine" -> "Exact", "Verify" -> True]["Theorem"]
