# Fisher's iris data (first fifteen flowers of each species): sepal length, sepal width, petal length, petal width in cm.
setosa = {{5.1, 3.5, 1.4, 0.2}, {4.9, 3.0, 1.4, 0.2}, {4.7, 3.2, 1.3, 0.2}, {4.6, 3.1, 1.5, 0.2}, {5.0, 3.6, 1.4, 0.2}, {5.4, 3.9, 1.7, 0.4}, {4.6, 3.4, 1.4, 0.3}, {5.0, 3.4, 1.5, 0.2}, {4.4, 2.9, 1.4, 0.2}, {4.9, 3.1, 1.5, 0.1}, {5.4, 3.7, 1.5, 0.2}, {4.8, 3.4, 1.6, 0.2}, {4.8, 3.0, 1.4, 0.1}, {4.3, 3.0, 1.1, 0.1}, {5.8, 4.0, 1.2, 0.2}};
versicolor = {{7.0, 3.2, 4.7, 1.4}, {6.4, 3.2, 4.5, 1.5}, {6.9, 3.1, 4.9, 1.5}, {5.5, 2.3, 4.0, 1.3}, {6.5, 2.8, 4.6, 1.5}, {5.7, 2.8, 4.5, 1.3}, {6.3, 3.3, 4.7, 1.6}, {4.9, 2.4, 3.3, 1.0}, {6.6, 2.9, 4.6, 1.3}, {5.2, 2.7, 3.9, 1.4}, {5.0, 2.0, 3.5, 1.0}, {5.9, 3.0, 4.2, 1.5}, {6.0, 2.2, 4.0, 1.0}, {6.1, 2.9, 4.7, 1.4}, {5.6, 2.9, 3.6, 1.3}};
virginica = {{6.3, 3.3, 6.0, 2.5}, {5.8, 2.7, 5.1, 1.9}, {7.1, 3.0, 5.9, 2.1}, {6.3, 2.9, 5.6, 1.8}, {6.5, 3.0, 5.8, 2.2}, {7.6, 3.0, 6.6, 2.1}, {4.9, 2.5, 4.5, 1.7}, {7.3, 2.9, 6.3, 1.8}, {6.7, 2.5, 5.8, 1.8}, {7.2, 3.6, 6.1, 2.5}, {6.5, 3.2, 5.1, 2.0}, {6.4, 2.7, 5.3, 1.9}, {6.8, 3.0, 5.5, 2.1}, {5.7, 2.5, 5.0, 2.0}, {5.8, 2.8, 5.1, 2.4}};
species = {"setosa", "versicolor", "virginica"}; groups = {setosa, versicolor, virginica};
train = Flatten[Table[Thread[Take[groups[[i]], 10] -> species[[i]]], {i, 3}]]; test = Flatten[Table[Thread[Take[groups[[i]], -5] -> species[[i]]], {i, 3}]];
flower = Classify[train, Method -> "LogisticRegression"]
Map[flower[First[#]] &, test]
Count[Map[flower[First[#]] === Last[#] &, test], True]
flower[{6.0, 2.9, 4.9, 1.7}, "Probabilities"]
Map[Function[m, Count[Map[Classify[train, Method -> m][First[#]] === Last[#] &, test], True]], {"NearestNeighbors", "NaiveBayes", "DecisionTree"}]
tree = Classify[train, Method -> "DecisionTree"]; tree[[2, 2]]
tree[[2, 3]]
all = Join[setosa, versicolor, virginica]; pc = PrincipalComponents[all];
v = Map[Variance, Transpose[pc]]; v/Total[v]
DimensionReduction[all, 2][{5.9, 3.0, 5.1, 1.8}]
km = FindClusters[all, 3, Method -> "KMeans"]; Table[Length[Intersection[km[[i]], groups[[j]]]], {i, 3}, {j, 3}]
