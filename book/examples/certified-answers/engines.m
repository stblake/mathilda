# Engines: the exact re-check, the heuristic, and consensus.
c = Certify[(x^2 - 1)/(x - 1) == x + 1, "Engine" -> "Exact"];
c["Verified"]
c["Checks"]
Quiet[Certify[(x + 1)^3 == x^3 + 3 x^2 + 2 x + 1, "Engine" -> "Exact"]["Verified"]]
Certify[(x + 1)^2 == x^2 + 2 x + 1, "Engine" -> "Mathilda"]["Verified"]
m = Quiet[Certify[Factor[x^10 - 1], "Engine" -> {"Lean", "Exact"}]];
m["Engines"]
m["Verified"]
Certify[Factor[x^10 - 1], "Engine" -> {"Exact", "Mathilda"}]["Verified"]
fake = ProofObject[<|"Kind" -> "Identity", "Verified" -> True|>];
fake["Verified"]
