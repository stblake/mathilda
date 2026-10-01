# k-nearest neighbours: votes, and the effect of feature scale.
SeedRandom[7]; blob[m_, s_, n_] := Transpose[{RandomVariate[NormalDistribution[m[[1]], s], n], RandomVariate[NormalDistribution[m[[2]], s], n]}];
{a, b, c} = {blob[{0, 0}, 0.8, 30], blob[{3, 3}, 0.8, 30], blob[{0, 4}, 0.8, 30]};
data = Join[Thread[a -> "A"], Thread[b -> "B"], Thread[c -> "C"]];
k5 = Classify[data, Method -> {"NearestNeighbors", "NeighborsNumber" -> 5}];
k5[{1.5, 2.}, "Probabilities"]
mm = Map[{1000 First[#][[1]], First[#][[2]]} -> Last[#] &, data];
Classify[mm][{0, 4.}]
X = Keys[mm]; mu = Mean[X]; sd = StandardDeviation[X];
std = Classify[Thread[Standardize[X] -> Values[mm]]];
std[({0, 4.} - mu)/sd]
