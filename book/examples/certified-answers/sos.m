# Sums of squares and an LP optimum.
c = Certify[ForAll[{x, y}, x^2 + y^2 >= 2 x y], "Verify" -> False];
StringSplit[c["ProofScript"], "\n"]
c = Certify[ForAll[{x, y}, x^4 + y^4 - x^2 y^2 + 1 >= 0], "Engine" -> "Exact"];
{c["Verified"], c["Certificate"]["Method"]}
c = Certify[Maximize[{3 x + 2 y, x + y <= 4, x + 3 y <= 6, x >= 0, y >= 0}, {x, y}], "Engine" -> "Exact"];
{c["Result"], c["Verified"]}
c["Certificate"]["Multipliers"]
c["Certificate"]["LinearProgram"]
