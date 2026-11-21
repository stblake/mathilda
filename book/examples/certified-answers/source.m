# The Lean file Certify writes, generated without running Lean.
c = Certify[Factor[x^4 - 1], "Verify" -> False];
c["Verified"]
c["LeanTheorem"]
StringSplit[c["ProofScript"], "\n"]
StringSplit[Last[StringSplit[c["LeanSource"], "end Mathilda\n\n"]], "\n"]
StringCount[c["LeanSource"], "\n"]
