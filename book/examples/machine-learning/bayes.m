# Gaussian naive Bayes, checked against the densities written out by hand.
lohi = {1. -> "lo", 2. -> "lo", 3. -> "lo", 4. -> "hi", 5. -> "lo", 6. -> "hi", 7. -> "hi", 8. -> "hi"};
nb = Classify[lohi, Method -> "NaiveBayes"];
nb[[2]]
nb[3.9, "Probabilities"]
{p, q} = {PDF[NormalDistribution[2.75, Sqrt[2.1875]], 3.9], PDF[NormalDistribution[6.25, Sqrt[2.1875]], 3.9]}; {p/(p + q), q/(p + q)}
Map[nb, {4.4, 4.5, 4.6}]
skew = Classify[{1. -> "lo", 2. -> "lo", 3. -> "lo", 5. -> "hi", 7. -> "hi", 9. -> "hi", 11. -> "hi", 13. -> "hi"}, Method -> "NaiveBayes"];
Map[skew, {-10., 0., 4., 20.}]
