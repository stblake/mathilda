# A boundary no straight line can draw, and a feature that makes it straight.
SeedRandom[5]; ring[r_, n_] := Table[With[{t = RandomReal[{0, 2 Pi}], s = r + RandomVariate[NormalDistribution[0, 0.25]]}, {s Cos[t], s Sin[t]}], {n}];
rd = Join[Thread[ring[1, 100] -> "in"], Thread[ring[2.5, 100] -> "out"]];
SeedRandom[6]; rt = Join[Thread[ring[1, 300] -> "in"], Thread[ring[2.5, 300] -> "out"]];
accuracy[f_, d_] := N[Mean[Boole[Map[f[First[#]] === Last[#] &, d]]]]
Map[{#, accuracy[Classify[rd, Method -> #], rt]} &, {"NearestNeighbors", "NaiveBayes", "LogisticRegression", "DecisionTree"}]
lift[p_ -> c_] := Append[p, p . p] -> c
lr3 = Classify[Map[lift, rd], Method -> "LogisticRegression"];
accuracy[lr3, Map[lift, rt]]
lr3[[2, 2]]
Sqrt[-lr3[[2, 2, 1]]/lr3[[2, 2, 4]]]
