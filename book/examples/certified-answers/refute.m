# Refute: a counterexample, and a kernel-checkable negation.
r = Refute[(x + 1)^3 == x^3 + 3 x^2 + 2 x + 1, "Verify" -> False];
r["Counterexample"]
r["LeanTheorem"]
r["ProofScript"]
Refute[ForAll[{x, y}, x^2 + y^2 >= 2 x y + 1/10], "Engine" -> "Exact"]["Counterexample"]
Refute[1/x == 1, "Engine" -> "Exact"]["Counterexample"]
r = Refute[Implies[x^2 + 4 == 0, x == 7], "Engine" -> "Exact"];
{r["Counterexample"], r["Verified"]}
Refute[x^2 >= x, Assumptions -> x > 0, "Engine" -> "Exact"]["Counterexample"]
Refute[ForAll[{x, y}, x^4 y^2 + x^2 y^4 - 3 x^2 y^2 + 1 >= 0], "Engine" -> "Exact"]
