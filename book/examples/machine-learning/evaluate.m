# Held-out accuracy, a confusion matrix, a method comparison, and cross-validation.
SeedRandom[7]; blob[m_, s_, n_] := Transpose[{RandomVariate[NormalDistribution[m[[1]], s], n], RandomVariate[NormalDistribution[m[[2]], s], n]}];
{a, b, c} = {blob[{0, 0}, 0.8, 30], blob[{3, 3}, 0.8, 30], blob[{0, 4}, 0.8, 30]};
data = Join[Thread[a -> "A"], Thread[b -> "B"], Thread[c -> "C"]];
SeedRandom[8]; test = Join[Thread[blob[{0, 0}, 0.8, 100] -> "A"], Thread[blob[{3, 3}, 0.8, 100] -> "B"], Thread[blob[{0, 4}, 0.8, 100] -> "C"]];
accuracy[f_, d_] := N[Mean[Boole[Map[f[First[#]] === Last[#] &, d]]]]
nn = Classify[data]; {accuracy[nn, data], accuracy[nn, test]}
pairs = Map[{Last[#], nn[First[#]]} &, test];
Table[Count[pairs, {i, j}], {i, {"A", "B", "C"}}, {j, {"A", "B", "C"}}]
methods = {"NearestNeighbors", "NaiveBayes", "LogisticRegression", "DecisionTree", "RandomForest"};
SeedRandom[1]; Table[With[{f = Classify[data, Method -> m]}, {m, accuracy[f, data], accuracy[f, test]}], {m, methods}]
SeedRandom[2]; folds = Partition[RandomSample[data], 18]; Length[folds]
cv[m_] := Mean[Table[accuracy[Classify[Flatten[Delete[folds, i]], Method -> m], folds[[i]]], {i, Length[folds]}]]
Map[cv, Most[methods]]
