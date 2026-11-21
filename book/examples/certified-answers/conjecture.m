# Conjecture: experiment, guess, test, certify.
c = Conjecture[Sum[k^4, {k, 1, n}], n];
c["Guess"]
{c["Form"], c["FittedOn"], c["TestedOn"]}
c["Status"]
c["Proof"]["Engines"]
c = Conjecture[Sum[k 3^k, {k, 0, n}], n];
{c["Guess"], c["Form"], c["Status"]}
c = Conjecture[Sum[Fibonacci[k]^2, {k, 1, n}], n];
{c["Guess"], c["Form"], c["Status"]}
c = Conjecture[Sum[Binomial[n, k]^2, {k, 0, n}], n];
{c["Guess"], c["Status"], c["Verified"]}
