# A CART tree read from its payload, then a random forest and its reproducibility.
lohi = {1. -> "lo", 2. -> "lo", 3. -> "lo", 4. -> "hi", 5. -> "lo", 6. -> "hi", 7. -> "hi", 8. -> "hi"};
t = Classify[lohi, Method -> "DecisionTree"];
{t[[4]], t[[2, 2]]}
t[[2, 3]]
Map[t, {3.9, 4.2, 5.2}]
SeedRandom[7]; blob[m_, s_, n_] := Transpose[{RandomVariate[NormalDistribution[m[[1]], s], n], RandomVariate[NormalDistribution[m[[2]], s], n]}];
{a, b, c} = {blob[{0, 0}, 0.8, 30], blob[{3, 3}, 0.8, 30], blob[{0, 4}, 0.8, 30]};
data = Join[Thread[a -> "A"], Thread[b -> "B"], Thread[c -> "C"]];
tree = Classify[data, Method -> "DecisionTree"]; tree[[4]]
Count[Map[tree[First[#]] === Last[#] &, data], True]
SeedRandom[42]; f1 = Classify[data, Method -> "RandomForest"]; SeedRandom[42]; f2 = Classify[data, Method -> "RandomForest"];
{f1 === f2, f1[[4]], Length[f1[[2, 2]]]}
Map[Length[First[#]] &, Take[f1[[2, 2]], 8]]
f1[{1.5, 2.}, "Probabilities"]
SeedRandom[43]; Classify[data, Method -> "RandomForest"][{1.5, 2.}, "Probabilities"]
