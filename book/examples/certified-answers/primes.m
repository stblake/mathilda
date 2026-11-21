# Primality: Pratt certificates.
c = Certify[FactorInteger[2^64 + 1], "Verify" -> False];
c["LeanTheorem"]
Certify[FactorInteger[2^64 + 1], "Engine" -> "Exact"]["Verified"]
Certify[PrimeQ[91], "Verify" -> False]["LeanTheorem"]
Lean`Internal`Oracle["Pratt", 274177]
