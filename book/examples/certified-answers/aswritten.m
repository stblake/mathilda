# The claim is certified as written.
Certify[Simplify[(x^2 - 1)/(x - 1)] == x + 1]
c = Certify[(x^2 - 1)/(x - 1) == x + 1, "Verify" -> False];
c["LeanTheorem"]
c["Hypotheses"]
c["TrustedSteps"]
x = 5;
Certify[x^2 == 25]
Clear[x]
Certify[Sin[x]^2 + Cos[x]^2 == 1]
c = Certify[(x + 1)^2 == x^2 + 2 x + 1, Assumptions -> x^2 == 2, "Verify" -> False];
c["LeanTheorem"]
c["Domain"]
