# Sums, roots, matrices, ideals.
SetOptions[Certify, "Verify" -> False];
Certify[Sum[k 2^k, {k, 0, n}]]["LeanTheorem"]
Certify[Solve[x^3 - 6 x^2 + 11 x - 6 == 0, x]]["LeanTheorem"]
Certify[Inverse[{{2, 1}, {1, 1}}]]["LeanTheorem"]
c = Certify[Implies[x^2 + y^2 == 1 && x == y, 2 y^2 == 1]];
c["Certificate"]["Cofactors"]
c["LeanTheorem"]
