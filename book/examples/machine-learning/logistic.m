# Logistic regression by IRLS on one feature, then one-vs-rest on three classes.
lohi = {1. -> "lo", 2. -> "lo", 3. -> "lo", 4. -> "hi", 5. -> "lo", 6. -> "hi", 7. -> "hi", 8. -> "hi"};
lr = Classify[lohi, Method -> "LogisticRegression"];
beta = lr[[2, 2]]
x0 = -beta[[1]]/beta[[2]]
lr[x0, "Probabilities"]
1/(1 + Exp[-(beta[[1]] + beta[[2]] 6.)])
lr[6., "Probabilities"]
sep = Classify[{1. -> "lo", 2. -> "lo", 8. -> "hi", 9. -> "hi"}, Method -> "LogisticRegression"];
sep[[2, 2]]
sep[5., "Probabilities"]
three = Classify[{1. -> "a", 2. -> "a", 5. -> "b", 6. -> "b", 9. -> "c", 10. -> "c"}, Method -> "LogisticRegression"];
three[[2, 2]]
three[5.5, "Probabilities"]
