# The fitted logistic curve for the one-feature "lo"/"hi" data; it crosses 1/2 at x = 4.5.
c1 = Classify[{1. -> "lo", 2. -> "lo", 3. -> "lo", 4. -> "hi", 5. -> "lo", 6. -> "hi", 7. -> "hi", 8. -> "hi"}, Method -> "LogisticRegression"];
fig = Plot[Last[Last[c1[x, "Probabilities"]]], {x, 0, 9}]
